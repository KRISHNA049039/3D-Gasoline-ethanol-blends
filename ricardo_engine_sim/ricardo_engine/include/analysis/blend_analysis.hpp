#pragma once
#include "utils/types.hpp"
#include "combustion/combustion_cycle.hpp"
#include "chemistry/blend_chemistry.hpp"
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <algorithm>

namespace Ricardo {
namespace Analysis {

struct BlendSweepPoint {
    FuelBlend blend;
    OperatingConditions ops;
    Combustion::CombustionCycle::CycleResult result;
    double knock_safety_margin;
};

class BlendSweep {
public:
    struct Config {
        std::vector<double> ethanol_fractions={0.0,0.05,0.10,0.15,0.20,0.30,0.50,0.85};
        OperatingConditions base_ops;
        EngineGeometry geom;
        bool adjust_spark=true;
    };
    explicit BlendSweep(const Config& cfg):cfg_(cfg){}

    std::vector<BlendSweepPoint> run(){
        std::vector<BlendSweepPoint> pts;
        for(double ef:cfg_.ethanol_fractions){
            FuelBlend b; b.ethanol_fraction=ef; b.gasoline_fraction=1.0-ef;
            OperatingConditions ops=cfg_.base_ops;
            if(cfg_.adjust_spark)
                ops.spark_timing_deg=cfg_.base_ops.spark_timing_deg+ef*4.0;

            Combustion::CombustionCycle::Config cyc_cfg;
            Combustion::CombustionCycle cycle(cfg_.geom,b,ops,cyc_cfg);
            auto res=cycle.run();

            BlendSweepPoint pt;
            pt.blend=b; pt.ops=ops; pt.result=res;
            pt.knock_safety_margin=1.0-res.knock_index;
            pts.push_back(std::move(pt));
        }
        return pts;
    }
private:
    Config cfg_;
};

struct EmissionsReport {
    double NOx_g_kWh,CO_g_kWh,HC_g_kWh,CO2_g_kWh,PM_mg_per_km;
    std::string euro_class;
    static EmissionsReport from_result(
            const Combustion::CombustionCycle::CycleResult& r,
            const FuelBlend& b){
        EmissionsReport e;
        e.NOx_g_kWh=r.NOx_g_kWh; e.CO_g_kWh=r.CO_g_kWh; e.HC_g_kWh=r.HC_g_kWh;
        double cfrac=b.carbon_atoms()*0.012/b.molar_mass();
        e.CO2_g_kWh=r.BSFC_g_kWh*cfrac*(44.0/12.0);
        e.PM_mg_per_km=(b.gasoline_fraction>0.5)?2.0+3.0*r.HC_g_kWh:0.5;
        bool eu6=r.NOx_g_kWh<0.06&&r.CO_g_kWh<1.0&&r.HC_g_kWh<0.10;
        e.euro_class=eu6?"Euro 6":"Pre-Euro 6";
        return e;
    }
};

class ReportGenerator {
public:
    static std::string text_report(
        const std::vector<BlendSweepPoint>& sw,
        const EngineGeometry& g,
        const OperatingConditions& ops)
    {
        std::ostringstream ss;
        ss<<std::fixed;
        auto hr=[&](int n){ss<<std::string(n,'-')<<"\n";};

        ss<<"=================================================================\n";
        ss<<"  RICARDO ENGINE COMBUSTION ANALYSIS - ETHANOL/GASOLINE BLENDS\n";
        ss<<"=================================================================\n\n";

        ss<<"ENGINE: Ricardo E6  Bore="<<g.bore*1000<<"mm  Stroke="<<g.stroke*1000
          <<"mm  CR="<<g.compression_ratio<<":1  Vd="
          <<std::setprecision(1)<<g.displacement()*1e6<<" cm3\n";
        ss<<"OPS: "<<ops.rpm<<" rpm  lambda="<<ops.lambda
          <<"  P_in="<<ops.intake_pressure/1000<<" kPa  T_in="<<ops.intake_temp<<" K\n\n";

        // ── Performance table ────────────────────────────────────────────
        ss<<"COMBUSTION PERFORMANCE\n"; hr(80);
        ss<<std::setw(6)<<"Blend"<<std::setw(9)<<"IMEP"<<std::setw(9)<<"eta_th"
          <<std::setw(10)<<"Ppeak"<<std::setw(9)<<"CA_Pk"<<std::setw(9)<<"MFB50"
          <<std::setw(9)<<"D10-90"<<std::setw(10)<<"Tpeak"<<std::setw(10)<<"BSFC\n";
        ss<<std::setw(6)<<""<<std::setw(9)<<"[bar]"<<std::setw(9)<<"[%]"
          <<std::setw(10)<<"[bar]"<<std::setw(9)<<"[°CA]"<<std::setw(9)<<"[°CA]"
          <<std::setw(9)<<"[°CA]"<<std::setw(10)<<"[K]"<<std::setw(10)<<"[g/kWh]\n";
        hr(80);
        for(auto& pt:sw){auto& r=pt.result;
            ss<<std::setw(6)<<pt.blend.name()
              <<std::setprecision(2)<<std::setw(9)<<r.IMEP_bar
              <<std::setw(9)<<r.thermal_efficiency*100.0
              <<std::setw(10)<<r.peak_pressure_bar
              <<std::setw(9)<<r.peak_pressure_CA
              <<std::setw(9)<<r.MFB50
              <<std::setw(9)<<r.CA10_90
              <<std::setprecision(0)<<std::setw(10)<<r.peak_temperature_K
              <<std::setprecision(1)<<std::setw(10)<<r.BSFC_g_kWh<<"\n";
        }

        // ── Flame speed ──────────────────────────────────────────────────
        ss<<"\nFLAME SPEED & COMBUSTION CHARACTERISTICS\n"; hr(65);
        ss<<std::setw(6)<<"Blend"<<std::setw(12)<<"Su [m/s]"<<std::setw(13)
          <<"St [m/s]"<<std::setw(12)<<"dL [um]"<<std::setw(10)<<"RON"
          <<std::setw(12)<<"Dur [deg]\n"; hr(65);
        for(auto& pt:sw){auto& r=pt.result;
            ss<<std::setw(6)<<pt.blend.name()
              <<std::setprecision(3)<<std::setw(12)<<r.Su_lam_m_s
              <<std::setw(13)<<r.St_turb_m_s
              <<std::setprecision(1)<<std::setw(12)<<r.delta_L_um
              <<std::setw(10)<<pt.blend.octane_number()
              <<std::setw(12)<<r.CA10_90<<"\n";
        }

        // ── Knock ────────────────────────────────────────────────────────
        ss<<"\nKNOCK ANALYSIS (Livengood-Wu integral)\n"; hr(55);
        ss<<std::setw(6)<<"Blend"<<std::setw(14)<<"Knock Index"
          <<std::setw(14)<<"Safety Margin"<<std::setw(20)<<"Status\n"; hr(55);
        for(auto& pt:sw){auto& r=pt.result;
            std::string stat=r.knock_detected?"** KNOCK **":
                             (r.knock_index>0.7?"Borderline":"Safe");
            ss<<std::setw(6)<<pt.blend.name()
              <<std::setprecision(3)<<std::setw(14)<<r.knock_index
              <<std::setw(14)<<(1.0-r.knock_index)
              <<std::setw(20)<<stat<<"\n";
        }

        // ── Emissions ────────────────────────────────────────────────────
        ss<<"\nEMISSIONS ANALYSIS\n"; hr(80);
        ss<<std::setw(6)<<"Blend"<<std::setw(13)<<"NOx[g/kWh]"<<std::setw(12)
          <<"CO[g/kWh]"<<std::setw(12)<<"HC[g/kWh]"<<std::setw(14)<<"CO2[g/kWh]"
          <<std::setw(12)<<"Euro\n"; hr(80);
        for(auto& pt:sw){
            auto er=EmissionsReport::from_result(pt.result,pt.blend);
            ss<<std::setw(6)<<pt.blend.name()
              <<std::setprecision(3)<<std::setw(13)<<er.NOx_g_kWh
              <<std::setw(12)<<er.CO_g_kWh
              <<std::setw(12)<<er.HC_g_kWh
              <<std::setprecision(1)<<std::setw(14)<<er.CO2_g_kWh
              <<std::setw(12)<<er.euro_class<<"\n";
        }

        // ── Insights ─────────────────────────────────────────────────────
        ss<<"\nKEY INSIGHTS\n"; hr(55);
        auto bIMEP=std::max_element(sw.begin(),sw.end(),
            [](auto& a,auto& b){return a.result.IMEP_bar<b.result.IMEP_bar;});
        auto bEff=std::max_element(sw.begin(),sw.end(),
            [](auto& a,auto& b){return a.result.thermal_efficiency<b.result.thermal_efficiency;});
        auto bKnk=std::max_element(sw.begin(),sw.end(),
            [](auto& a,auto& b){return a.knock_safety_margin<b.knock_safety_margin;});

        ss<<"  Best IMEP          : "<<bIMEP->blend.name()<<" ("
          <<std::setprecision(2)<<bIMEP->result.IMEP_bar<<" bar)\n";
        ss<<"  Best Thermal Eff.  : "<<bEff->blend.name()<<" ("
          <<std::setprecision(1)<<bEff->result.thermal_efficiency*100.0<<" %)\n";
        ss<<"  Best Knock Margin  : "<<bKnk->blend.name()<<" (margin="
          <<std::setprecision(2)<<bKnk->knock_safety_margin<<")\n";

        ss<<"\n  Laminar Flame Speed E0->"<<sw.back().blend.name()<<": "
          <<std::setprecision(2)<<sw.front().result.Su_lam_m_s<<" -> "
          <<sw.back().result.Su_lam_m_s<<" m/s (+";
        double pct=(sw.back().result.Su_lam_m_s/sw.front().result.Su_lam_m_s-1.0)*100.0;
        ss<<std::setprecision(1)<<pct<<"%)\n";

        auto er0=EmissionsReport::from_result(sw.front().result,sw.front().blend);
        auto er85=EmissionsReport::from_result(sw.back().result,sw.back().blend);
        double co2red=(er0.CO2_g_kWh-er85.CO2_g_kWh)/(er0.CO2_g_kWh+1e-9)*100.0;
        ss<<"  Fossil CO2 reduction E0->"<<sw.back().blend.name()<<": "
          <<std::setprecision(1)<<co2red<<" %\n";

        ss<<"=================================================================\n";
        return ss.str();
    }

    static std::string pressure_trace_csv(
        const Combustion::CombustionCycle::CycleResult& r)
    {
        std::ostringstream ss;
        ss<<"CA_deg,P_bar,T_u_K,T_b_K,xb,ROHR_J_deg,NOx_ppm,knock_integral\n";
        for(size_t i=0;i<r.CA_deg.size();i++){
            ss<<std::fixed<<std::setprecision(3)
              <<r.CA_deg[i]<<","<<r.P_bar[i]<<","
              <<(i<r.T_u_K.size()?r.T_u_K[i]:0.0)<<","
              <<(i<r.T_b_K.size()?r.T_b_K[i]:0.0)<<","
              <<(i<r.xb_vec.size()?r.xb_vec[i]:0.0)<<","
              <<(i<r.ROHR_J_deg.size()?r.ROHR_J_deg[i]:0.0)<<","
              <<(i<r.NOx_ppm.size()?r.NOx_ppm[i]:0.0)<<","
              <<(i<r.knock_I.size()?r.knock_I[i]:0.0)<<"\n";
        }
        return ss.str();
    }
};

} // namespace Analysis
} // namespace Ricardo
