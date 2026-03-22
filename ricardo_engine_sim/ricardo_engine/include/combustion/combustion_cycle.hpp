#pragma once
#include "utils/types.hpp"
#include "chemistry/blend_chemistry.hpp"
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>

namespace Ricardo {
namespace Combustion {

class WiebeFunction {
public:
    struct Params {
        double a=5.0, m=2.0, theta_0=0.0, delta_theta=60.0;
        Params()=default;
    };
    explicit WiebeFunction() : p_(Params()) {}
    explicit WiebeFunction(const Params& p) : p_(p) {}
    double xb(double theta) const {
        double xi=(theta-p_.theta_0)/p_.delta_theta;
        if(xi<=0.0)return 0.0; if(xi>=1.0)return 1.0;
        return 1.0-std::exp(-p_.a*std::pow(xi,p_.m+1.0));
    }
    void set_start(double t){p_.theta_0=t;}
    void set_duration(double d){p_.delta_theta=d;}
    double start()const{return p_.theta_0;}
    double duration()const{return p_.delta_theta;}
    void calibrate(double Su,double St,double bore,const OperatingConditions& ops){
        double S_avg=0.5*(Su+St);
        double t_burn=0.75*bore/std::max(S_avg,0.3);
        double dps=ops.rpm*6.0;
        p_.delta_theta=std::clamp(t_burn*dps,20.0,80.0);
    }
private:
    Params p_;
};

class CombustionCycle {
public:
    struct Config {
        double dtheta=0.25;
        bool heat_transfer=true,track_knock=true,track_emissions=true;
        Config()=default;
    };

    struct CycleResult {
        double IMEP_bar=0,W_indicated_J=0,thermal_efficiency=0;
        double comb_efficiency=0,BSFC_g_kWh=0;
        double peak_pressure_bar=0,peak_pressure_CA=0,peak_temperature_K=0;
        double MFB10=0,MFB50=0,MFB90=0,CA10_90=0,spark_CA=0;
        double Su_lam_m_s=0,St_turb_m_s=0,delta_L_um=0;
        double knock_index=0; bool knock_detected=false;
        double NOx_g_kWh=0,CO_g_kWh=0,HC_g_kWh=0;
        std::vector<double> CA_deg,P_bar,T_u_K,T_b_K,xb_vec,ROHR_J_deg,NOx_ppm,knock_I;
    };

    CombustionCycle(const EngineGeometry& g,const FuelBlend& b,
                    const OperatingConditions& o,const Config& c)
        :g_(g),b_(b),o_(o),cfg_(c),chem_(b){}

    CycleResult run(){
        CycleResult res;
        res.spark_CA=o_.spark_timing_deg;
        const double Ru=287.5,Rb=290.0,Cv_u=740.0,Cv_b=850.0;
        const double gamma_u=1.35,gamma_b=1.28;

        // Initial state at BDC (-360 deg)
        double P=o_.intake_pressure;
        double V=g_.bdc_volume();
        double T_u=o_.intake_temp;
        double m=P*V/(Ru*T_u);
        double AFR=b_.stoich_afr()*o_.lambda;
        double m_fuel=m/(1.0+AFR);
        double Q_fuel=m_fuel*b_.lhv();

        double T_b=T_u,m_b=0.0,xb=0.0,Q_released=0.0;
        double knock_I_val=0.0;
        double NOx_total_kg=0.0,CO_total_kg=0.0;
        double phi=1.0/o_.lambda;

        // Flame speed at estimated spark conditions
        double ratio_c=g_.bdc_volume()/g_.tdc_volume();
        double T_spark=T_u*std::pow(ratio_c,gamma_u-1.0);
        double P_spark=P*std::pow(ratio_c,gamma_u);
        double Su=chem_.laminar_flame_speed(T_spark,P_spark,phi);
        double dL=chem_.laminar_flame_thickness(T_spark,P_spark,phi);
        double u_prime=0.5*2.0*g_.stroke*o_.rpm/60.0;
        double l_t=0.06*g_.bore;
        double St=chem_.turbulent_flame_speed(Su,u_prime,l_t,dL);

        WiebeFunction wiebe;
        wiebe.set_start(o_.spark_timing_deg);
        wiebe.calibrate(Su,St,g_.bore,o_);

        res.Su_lam_m_s=Su; res.St_turb_m_s=St; res.delta_L_um=dL*1e6;

        bool f10=false,f50=false,f90=false;
        double W_net=0.0;
        double V_prev=V,P_prev=P;

        for(double ca=-360.0+cfg_.dtheta; ca<=180.0+0.001; ca+=cfg_.dtheta){
            V=g_.cylinder_volume(ca);
            double dV=V-V_prev;

            // Heat release
            double xb_new=wiebe.xb(ca);
            double dxb=xb_new-xb; xb=xb_new;
            double dQ=dxb*Q_fuel*0.995; // 99.5% combustion efficiency // 97% combustion efficiency
            Q_released+=dQ;
            double dm_b=dxb*m;
            m_b=std::min(m_b+dm_b,m);
            double m_u=std::max(0.0,m-m_b);

            // Effective gamma for mixed zones
            double g_eff=(m_b>0.001*m)?
                (m_u*Cv_u*gamma_u+m_b*Cv_b*gamma_b)/(m_u*Cv_u+m_b*Cv_b+1e-30)
                :gamma_u;

            // First law: dP = (g-1)/V*dQ - g*P/V*dV
            double dP=(g_eff-1.0)/V*dQ - g_eff*P/V*dV;
            P=std::max(500.0,P+dP);

            // Woschni heat loss
            if(cfg_.heat_transfer){
                double bore=g_.bore;
                double h_p=V/(Constants::PI/4.0*bore*bore);
                double A_w=2.0*Constants::PI/4.0*bore*bore+Constants::PI*bore*h_p;
                double w=2.28*2.0*g_.stroke*o_.rpm/60.0;
                double htc=3.26*std::pow(bore,-0.2)*std::pow(P/1000.0,0.8)*
                           std::pow(std::max(T_u,300.0),-0.55)*std::pow(w,0.8);
                double dt=cfg_.dtheta/o_.crank_deg_per_second();
                double dQ_loss=htc*A_w*(T_u-o_.coolant_temp)*dt;
                P=std::max(500.0,P-(g_eff-1.0)/V*dQ_loss);
            }

            // Zone temperatures
            if(ca<=o_.spark_timing_deg||m_b<1e-12*m){
                // isentropic unburned
                if(std::abs(dV)>0&&V_prev>0)
                    T_u=std::max(200.0,T_u*std::pow(V_prev/V,gamma_u-1.0));
            } else {
                // Burned zone: adiabatic mixing
                if(dm_b>1e-15){
                    double h_comb=b_.lhv()/(1.0+AFR);
                    double T_ad=T_u+dm_b*h_comb/(m_b*Cv_b+1e-20);
                    T_b=(T_b*(m_b-dm_b)+dm_b*T_ad)/std::max(m_b,1e-20);
                    T_b=std::clamp(T_b,T_u,3000.0);
                }
                // Unburned: from ideal gas (V_u = V - V_b)
                double V_b=m_b*Rb*T_b/(P+1e-9);
                V_b=std::clamp(V_b,0.0,V*0.999);
                double V_u=V-V_b;
                double m_u2=std::max(m_u,1e-20);
                T_u=P*V_u/(m_u2*Ru);
                T_u=std::clamp(T_u,200.0,2500.0);
            }

            // Work (trapezoidal)
            W_net+=0.5*(P+P_prev)*dV;

            // Peak tracking
            if(P/1e5>res.peak_pressure_bar){
                res.peak_pressure_bar=P/1e5;
                res.peak_pressure_CA=ca;
            }
            double T_max=std::max(T_u,T_b);
            if(T_max>res.peak_temperature_K) res.peak_temperature_K=T_max;

            // MFB
            if(!f10&&xb>=0.10){res.MFB10=ca;f10=true;}
            if(!f50&&xb>=0.50){res.MFB50=ca;f50=true;}
            if(!f90&&xb>=0.90){res.MFB90=ca;f90=true;}

            // Knock (Livengood-Wu)
            if(cfg_.track_knock&&!res.knock_detected&&ca>o_.spark_timing_deg&&m_u>0.001*m){
                double dt=cfg_.dtheta/o_.crank_deg_per_second();
                knock_I_val+=chem_.knock_integral_increment(T_u,P,phi,dt);
                if(knock_I_val>=1.0) res.knock_detected=true;
            }

            // Emissions
            if(cfg_.track_emissions&&m_b>1e-10&&T_b>1600.0){
                double dt=cfg_.dtheta/o_.crank_deg_per_second();
                double V_b_m3=m_b*Rb*T_b/(P+1e-9);
                double rNO=chem_.NOx_formation_rate(T_b,P,phi);
                NOx_total_kg+=rNO*0.030*V_b_m3*dt;
                if(phi>1.0){
                    double rCO=std::max(0.0,chem_.CO_rate(T_b,P,phi,0.005));
                    CO_total_kg+=rCO*0.028*V_b_m3*dt;
                }
            }

            res.CA_deg.push_back(ca);
            res.P_bar.push_back(P/1e5);
            res.T_u_K.push_back(T_u);
            res.T_b_K.push_back(T_b);
            res.xb_vec.push_back(xb);
            res.ROHR_J_deg.push_back(dQ/cfg_.dtheta);
            res.NOx_ppm.push_back(NOx_total_kg/(m*0.030+1e-20)*1e6);
            res.knock_I.push_back(knock_I_val);

            V_prev=V; P_prev=P;
        }

        res.knock_index=knock_I_val;
        res.CA10_90=res.MFB90-res.MFB10;
        res.W_indicated_J=W_net;
        res.IMEP_bar=W_net/g_.displacement()/1e5;
        res.thermal_efficiency=W_net/(Q_fuel+1e-9);
        res.comb_efficiency=Q_released/(Q_fuel+1e-9);

        double W_kWh=std::max(W_net,10.0)/3.6e6;
        res.NOx_g_kWh=NOx_total_kg*1000.0/W_kWh;
        res.CO_g_kWh=CO_total_kg*1000.0/W_kWh;
        res.HC_g_kWh=std::max(0.0,(1.0-res.comb_efficiency))*m_fuel*1000.0/W_kWh;
        res.BSFC_g_kWh=m_fuel*1000.0/W_kWh;

        return res;
    }

private:
    EngineGeometry g_; FuelBlend b_; OperatingConditions o_; Config cfg_;
    Chemistry::BlendChemistry chem_;
};

} // namespace Combustion
} // namespace Ricardo
