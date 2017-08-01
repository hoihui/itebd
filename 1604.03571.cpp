#include <iostream>
#include <fstream>
#include <string>
#include "itensor/all.h"
#include "itebd.h"
using namespace itensor;


int main ( int argc, char *argv[] ) {
    std::string output_file("fig2a.dat");
    double hx,hz;
    SpinHalf sites(2);

    hx=0; hz=0.5;
    AutoMPO ampo(sites);
    ampo += -4,"Sz",1,"Sz",2;
    ampo += -hz*2/2,"Sx",1;
    ampo += -hz*2/2,"Sx",2;
    ampo += -hx*2/2,"Sz",1;
    ampo += -hx*2/2,"Sz",2;
    itebd<Index> Sys(ampo);// Sys.randomize();
    for (int e=0;e>-8;e-=1) std::cout<<Sys.step_imag(std::pow(10,e),1000)<<std::endl;


    hx=0.1;hz=0.25;
    double T=20, dt=.001, thres=1E-12, M=100;
    ampo = AutoMPO(sites);
    ampo += -4,"Sz",1,"Sz",2;
    ampo += -hz*2/2,"Sx",1;
    ampo += -hz*2/2,"Sx",2;
    ampo += -hx*2/2,"Sz",1;
    ampo += -hx*2/2,"Sz",2;
    auto op=sites.op("Sz",1)*sites.op("Id",2)+sites.op("Id",1)*sites.op("Sz",2);
    std::vector<double> opv(1,std::real(Sys.measure(op)));
    Sys.resetTime();
    Sys.setH(ampo);

    while (Sys.tv.back()<T) {
        Sys.step(dt,10,thres,M);
        opv.push_back(std::real(Sys.measure(op)));
        std::cout<<Sys.tv.back()<<":  "<<Sys.bonddim<<"  "<<opv.back()<<std::endl;
    }


    std::ofstream file;
    file.open(output_file);
    for(auto i=0; i<Sys.tv.size(); i++)
        file<<Sys.tv[i]<<" "<<opv[i]<<std::endl;
    file.close();
    return 0;
}
