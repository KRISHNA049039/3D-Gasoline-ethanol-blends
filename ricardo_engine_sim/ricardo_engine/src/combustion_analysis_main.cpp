// =========================================================================
//  Ricardo Engine 3D CFD & Combustion Analysis Suite
//  Ethanol-Gasoline Blend Analysis  |  C++ Simulation Framework
// =========================================================================
#include "utils/types.hpp"
#include "mesh/cylinder_mesh.hpp"
#include "physics/turbulence.hpp"
#include "chemistry/blend_chemistry.hpp"
#include "combustion/combustion_cycle.hpp"
#include "analysis/blend_analysis.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <chrono>
#include <cstdlib>

using namespace Ricardo;

void section(const std::string& t){
    std::cout<<"\n\033[1;36m"<<t<<"\033[0m\n"<<std::string(t.size(),'=')<<"\n";
}

// ── 3D Mesh demonstration ─────────────────────────────────────────────────
void run_mesh_demo(const EngineGeometry& g){
    section("3D CYLINDER MESH GENERATION (O-grid structured)");
    Mesh::CylinderMeshGenerator::Config cfg;
    cfg.nr=10; cfg.nth=20; cfg.nz=25;
    Mesh::CylinderMeshGenerator mg(g,cfg);

    std::cout<<"  CA [deg]  | Volume [cm3] | Nodes  | Cells  | V_mesh [m3]\n";
    std::cout<<"  "<<std::string(60,'-')<<"\n";
    for(double ca:{-360.0,-180.0,-90.0,0.0,90.0,180.0}){
        mg.generate(ca);
        auto s=Mesh::MeshStats::compute(mg);
        std::cout<<std::fixed<<std::setprecision(2)
                 <<"  "<<std::setw(8)<<ca
                 <<"  | "<<std::setw(10)<<g.cylinder_volume(ca)*1e6
                 <<"   | "<<std::setw(6)<<s.n_nodes
                 <<"  | "<<std::setw(6)<<s.n_cells
                 <<"  | "<<std::scientific<<std::setprecision(3)<<s.total_volume<<"\n";
    }
    mg.generate(0.0);
    std::cout<<"\n  TDC zones:\n";
    for(auto& z:mg.zones())
        std::cout<<"    ["<<z.id<<"] "<<z.name<<" | type="<<z.type
                 <<" | cells="<<z.cell_ids.size()
                 <<" | bnd_faces="<<z.face_ids.size()<<"\n";
}

// ── Turbulence k-eps ──────────────────────────────────────────────────────
void run_turbulence_demo(const EngineGeometry& g, const OperatingConditions& o){
    section("k-epsilon TURBULENCE MODEL");
    Physics::EngineTurbulence turb(g,o);
    turb.initialize();
    std::cout<<"  CA[deg]  |  k[m2/s2]  |  eps[m2/s3] | u'[m/s] | l_t[mm]\n";
    std::cout<<"  "<<std::string(60,'-')<<"\n";
    double prev=-360.0;
    for(double ca:{-360.0,-300.0,-240.0,-180.0,-120.0,-60.0,-30.0,0.0,30.0}){
        if(ca>prev)turb.advance(prev,ca,1.2);
        auto& ts=turb.state();
        std::cout<<std::fixed<<std::setprecision(4)
                 <<"  "<<std::setw(7)<<ca
                 <<"  |  "<<std::setw(8)<<ts.k
                 <<"  |  "<<std::setw(9)<<ts.eps
                 <<"  | "<<std::setw(7)<<ts.u_prime
                 <<"  | "<<std::setw(7)<<ts.l_int*1000<<"\n";
        prev=ca;
    }
}

// ── Chemistry properties ──────────────────────────────────────────────────
void run_chemistry_demo(){
    section("FUEL BLEND CHEMISTRY PROPERTIES");
    std::cout<<"  T=700K, P=15bar, phi=1.0\n";
    std::cout<<"  Blend | LHV[MJ/kg] | AFR_st | RON   | Su[m/s] | dL[um] | rho[kg/m3]\n";
    std::cout<<"  "<<std::string(72,'-')<<"\n";
    for(double ef:{0.0,0.10,0.20,0.30,0.50,0.85}){
        FuelBlend b; b.ethanol_fraction=ef; b.gasoline_fraction=1-ef;
        Chemistry::BlendChemistry ch(b);
        double Su=ch.laminar_flame_speed(700,15e5,1.0);
        double dL=ch.laminar_flame_thickness(700,15e5,1.0)*1e6;
        std::cout<<std::fixed<<std::setprecision(2)
                 <<"  "<<std::setw(5)<<b.name()
                 <<" | "<<std::setw(9)<<b.lhv()/1e6
                 <<" | "<<std::setw(6)<<b.stoich_afr()
                 <<"  | "<<std::setw(5)<<b.octane_number()
                 <<"  | "<<std::setw(7)<<Su
                 <<"   | "<<std::setw(6)<<dL
                 <<"  | "<<std::setw(8)<<b.density_liquid()<<"\n";
    }
}

// ── Flame speed map ───────────────────────────────────────────────────────
void run_flame_speed_study(){
    section("LAMINAR FLAME SPEED PARAMETRIC STUDY  [m/s] @ P=15bar,T=700K");
    std::vector<double> phis={0.7,0.8,0.9,1.0,1.1,1.2,1.3};
    std::vector<double> blends={0.0,0.10,0.30,0.50,0.85};
    std::cout<<"  Blend  |";
    for(double phi:phis) std::cout<<std::setw(8)<<std::fixed<<std::setprecision(1)<<phi;
    std::cout<<"\n  "<<std::string(8+8*phis.size(),'-')<<"\n";
    std::ofstream fs("results/flame_speed_map.csv");
    fs<<"blend,eth_pct,phi,T_K,P_Pa,Su_m_s,dL_um\n";
    for(double ef:blends){
        FuelBlend b; b.ethanol_fraction=ef; b.gasoline_fraction=1-ef;
        Chemistry::BlendChemistry ch(b);
        std::cout<<"  "<<std::setw(5)<<b.name()<<"  |";
        for(double phi:phis){
            double Su=ch.laminar_flame_speed(700,15e5,phi);
            std::cout<<std::setw(8)<<std::fixed<<std::setprecision(3)<<Su;
            double dL=ch.laminar_flame_thickness(700,15e5,phi)*1e6;
            for(double T:{400,500,600,700,800})
                fs<<b.name()<<","<<ef*100<<","<<phi<<","<<T<<",15e5,"
                  <<ch.laminar_flame_speed(T,15e5,phi)<<","
                  <<ch.laminar_flame_thickness(T,15e5,phi)*1e6<<"\n";
            (void)dL;
        }
        std::cout<<"\n";
    }
    fs.close();
    std::cout<<"  -> results/flame_speed_map.csv\n";
}

// ── Knock study ───────────────────────────────────────────────────────────
void run_knock_study(const EngineGeometry& g, const OperatingConditions& base){
    section("KNOCK INDEX SWEEP (Livengood-Wu integral, threshold=1.0)");
    std::vector<double> rpms={1000,1500,2000,3000,4000};
    std::vector<double> blends={0.0,0.10,0.30,0.85};
    std::cout<<"  Blend\\RPM |";
    for(double r:rpms) std::cout<<std::setw(10)<<(int)r;
    std::cout<<"\n  "<<std::string(11+10*rpms.size(),'-')<<"\n";
    for(double ef:blends){
        FuelBlend b; b.ethanol_fraction=ef; b.gasoline_fraction=1-ef;
        std::cout<<"  "<<std::setw(8)<<b.name()<<"  |";
        for(double rpm:rpms){
            OperatingConditions o=base; o.rpm=rpm;
            Combustion::CombustionCycle::Config cfg;
            cfg.heat_transfer=true; cfg.track_knock=true; cfg.track_emissions=false;
            auto r=Combustion::CombustionCycle(g,b,o,cfg).run();
            std::string flag=(r.knock_detected?"*":"");
            std::cout<<std::setw(9)<<std::fixed<<std::setprecision(3)<<r.knock_index<<flag;
        }
        std::cout<<"\n";
    }
    std::cout<<"  (* = knock detected)\n";
}

// ── Full combustion sweep ─────────────────────────────────────────────────
void run_combustion_analysis(const EngineGeometry& g, const OperatingConditions& ops){
    section("FULL COMBUSTION CYCLE ANALYSIS — ETHANOL/GASOLINE BLEND SWEEP");

    Analysis::BlendSweep::Config cfg;
    cfg.geom=g; cfg.base_ops=ops; cfg.adjust_spark=true;
    cfg.ethanol_fractions={0.0,0.05,0.10,0.15,0.20,0.30,0.50,0.85};

    std::cout<<"  Sweeping "<<cfg.ethanol_fractions.size()<<" blends...\n";
    auto t0=std::chrono::steady_clock::now();
    Analysis::BlendSweep sw(cfg);
    auto pts=sw.run();
    auto t1=std::chrono::steady_clock::now();
    std::cout<<"  Done in "<<std::fixed<<std::setprecision(3)
             <<std::chrono::duration<double>(t1-t0).count()<<" s\n\n";

    // Print report
    std::cout<<Analysis::ReportGenerator::text_report(pts,g,ops);

    // CSV summary
    std::ofstream sum("results/blend_summary.csv");
    sum<<"blend,eth_pct,IMEP_bar,eta_th_pct,Ppeak_bar,CA_Ppeak,"
       <<"MFB10,MFB50,MFB90,CA10_90,Tpeak_K,BSFC_g_kWh,"
       <<"Su_lam,St_turb,dL_um,"
       <<"knock_idx,knock_margin,knocked,"
       <<"NOx_g_kWh,CO_g_kWh,HC_g_kWh,CO2_g_kWh\n";
    for(auto& pt:pts){
        auto& r=pt.result;
        Analysis::EmissionsReport er=Analysis::EmissionsReport::from_result(r,pt.blend);
        sum<<std::fixed<<std::setprecision(4)
           <<pt.blend.name()<<","<<pt.blend.ethanol_fraction*100<<","
           <<r.IMEP_bar<<","<<r.thermal_efficiency*100<<","
           <<r.peak_pressure_bar<<","<<r.peak_pressure_CA<<","
           <<r.MFB10<<","<<r.MFB50<<","<<r.MFB90<<","<<r.CA10_90<<","
           <<r.peak_temperature_K<<","<<r.BSFC_g_kWh<<","
           <<r.Su_lam_m_s<<","<<r.St_turb_m_s<<","<<r.delta_L_um<<","
           <<r.knock_index<<","<<pt.knock_safety_margin<<","<<(r.knock_detected?1:0)<<","
           <<r.NOx_g_kWh<<","<<r.CO_g_kWh<<","<<r.HC_g_kWh<<","<<er.CO2_g_kWh<<"\n";
    }
    sum.close();

    // P-CA traces
    for(auto& pt:pts){
        std::ofstream tr("results/P_CA_"+pt.blend.name()+".csv");
        tr<<Analysis::ReportGenerator::pressure_trace_csv(pt.result);
        tr.close();
    }
    std::cout<<"\n  Exported:\n    results/blend_summary.csv\n    results/P_CA_E*.csv\n"
             <<"    results/flame_speed_map.csv\n";
}

int main(){
    std::cout<<"\033[1;32m"
             <<"+=============================================================+\n"
             <<"|   RICARDO ENGINE 3D CFD & COMBUSTION ANALYSIS SUITE        |\n"
             <<"|   Ethanol-Gasoline Blend Simulation  |  C++ Framework       |\n"
             <<"+=============================================================+\n"
             <<"\033[0m\n";

    if(std::system("mkdir -p results")!=0){}

    // Ricardo E6 research engine geometry
    EngineGeometry g;
    g.bore=0.0760; g.stroke=0.1110; g.con_rod_length=0.1900;
    g.compression_ratio=8.0; g.num_cylinders=1;

    OperatingConditions ops;
    ops.rpm=2000; ops.lambda=1.0; ops.spark_timing_deg=-20.0;
    ops.intake_pressure=101325; ops.intake_temp=320; ops.coolant_temp=353;

    std::cout<<"Ricardo E6 Research Engine\n";
    std::cout<<"  Bore x Stroke = "<<g.bore*1000<<" x "<<g.stroke*1000
             <<" mm  |  CR="<<g.compression_ratio
             <<":1  |  Vd="<<std::fixed<<std::setprecision(1)
             <<g.displacement()*1e6<<" cm3\n";
    std::cout<<"  BDC volume="<<g.bdc_volume()*1e6<<"  TDC volume="
             <<g.tdc_volume()*1e6<<" cm3\n";

    run_mesh_demo(g);
    run_turbulence_demo(g,ops);
    run_chemistry_demo();
    run_flame_speed_study();
    run_knock_study(g,ops);
    run_combustion_analysis(g,ops);

    std::cout<<"\n\033[1;32m✓ All analyses complete.\033[0m\n\n";
    return 0;
}
