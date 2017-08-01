#include <iostream>
#include <fstream>
#include <complex>
#include <string>
#include "itensor/all.h"
#include "itebd.h"
using namespace itensor;
//const std::complex<double> I(0.0,1.0);


int main ( int argc, char *argv[] ) {
    std::string out("");
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "-out") out = argv[++i];
    }

    double J(-1),h(.5),dt(.001),thres(1E-5);  //Ising with transverse field
    unsigned z(4);
    std::string output_file="TFIM.dat";
    SpinHalf sites(2);
    auto mz=1.0*sites.op("Sz",1)*sites.op("Id",2)+1.0*sites.op("Id",1)*sites.op("Sz",2);
    auto mz1=sites.op("Sz",1)*sites.op("Sz",2);
    AutoMPO ampo(sites);
    ampo += 4*J,"Sz",1,"Sz",2;
    ampo += h*2/z,"Sx",1;
    ampo += h*2/z,"Sx",2;

    std::vector<double> hv, mzv;
    double dh=0.1,H=5;
    for (int i=0;i<=H/dh;i++){
        h=i*dh;
        AutoMPO ampo(sites);
        ampo += 4*J,"Sz",1,"Sz",2;
        ampo += h*2/z,"Sx",1;
        ampo += h*2/z,"Sx",2;
        itebd<Index> Sys(z,ampo,thres);
        for (double e=0;e>-5;e-=.9) Sys.step_imag(std::pow(10,e),100);
        hv.push_back(h);
        mzv.push_back(std::abs(Sys.measure(mz)));
        std::cout<<hv.back()<<": "<<mzv.back()<<std::endl;
    }

    std::ofstream file;
    file.open(output_file);
    for(auto i=0; i<hv.size(); i++){
        file<<hv[i]<<" "<<mzv[i]<<std::endl;
    }
    file.close();
    return 0;
}
