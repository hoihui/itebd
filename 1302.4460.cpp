#include <iostream>
#include <fstream>
#include <complex>
#include <string>
#include "itensor/all.h"
#include "itebd.h"
using namespace itensor;
const std::complex<double> I(0.0,1.0);

int main ( int argc, char *argv[] ) {
    std::string output_file("fig2a.dat");
    int z(2);
    double D(0.25),A(1),Omega(1.4);
    double thres=1E-8, M=100;
    SpinOne sites(2);

    AutoMPO ampo0(sites);
    ampo0 += "Sz",1,"Sz",2;
    ampo0 += .5,"S+",1,"S-",2;
    ampo0 += .5,"S-",1,"S+",2;
    ampo0 += D/z,"Sz2",1;
    ampo0 += D/z,"Sz2",2;
    itebd<Index> Sys(ampo0,2,2);// Sys.randomize();
    for (double e=0;e>-6;e-=.5) std::cout<<Sys.step_imag(std::pow(10,e),10,thres,M)<<std::endl;

    double T=9.6, dt=.01;
    auto mz=0.5*sites.op("Sz",1)*sites.op("Id",2)+0.5*sites.op("Id",1)*sites.op("Sz",2);
    std::vector<double> mzv(1,std::real(Sys.measure(mz)));
    Sys.resetTime();

    while (Sys.tv.back()<T) {
        double t=Sys.tv.back();
        AutoMPO ampo(sites);
        ampo += "Sz",1,"Sz",2;
        ampo += .5,"S+",1,"S-",2;
        ampo += .5,"S-",1,"S+",2;
        ampo += D/z,"Sz2",1;
        ampo += D/z,"Sz2",2;
        ampo += -A*std::exp(-I*Omega*t)/z,"S+",1;
        ampo += -A*std::exp(-I*Omega*t)/z,"S+",2;
        ampo += -A*std::exp(I*Omega*t)/z,"S-",1;
        ampo += -A*std::exp(I*Omega*t)/z,"S-",2;
        Sys.setH(ampo);
        Sys.step(dt,1,thres,M);
        mzv.push_back(std::real(Sys.measure(mz)));
        std::cout<<Sys.tv.back()<<":  "<<Sys.bonddim<<"  "<<mzv.back()<<std::endl;
    }


    std::ofstream file;
    file.open(output_file);
    for(auto i=0; i<Sys.tv.size(); i++)
        file<<Sys.tv[i]<<" "<<mzv[i]<<std::endl;
    file.close();
    return 0;
}
