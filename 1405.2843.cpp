#include <iostream>
#include <fstream>
#include <complex>
#include <string>
#include "itensor/all.h"
#include "itebd.h"
using namespace itensor;


int main ( int argc, char *argv[] ) {
    std::string output_file("tableS1.dat");
    double thres(1E-8),M(400),dt(0.001),T(3);
    double Delta(3);
    SpinHalf sites(2);

    AutoMPO ampo(sites);
    auto op=4*sites.op("Sz",1)*sites.op("Sz",2);
    ampo += -4*Delta,"Sz",1,"Sz",2;
    ampo += -4/2,"Sp",1,"Sm",2;
    ampo += -4/2,"Sm",1,"Sp",2;
    itebd<IQIndex> Sys(ampo,1,2); //Neel

    std::vector<double> opv(1,std::real(Sys.measure(op)));
    std::cout<<0<<":  "<<opv.back()<<std::endl;
    while (Sys.tv.back()<T) {
        Sys.step(dt,10,thres,M);
        opv.push_back(Sys.measure(op));
        std::cout<<Sys.tv.back()<<":  "<<Sys.bonddim<<"  "<<opv.back()<<std::endl;
    }

//    std::ofstream file;
//    file.open(output_file);
//    for(auto i=0; i<Sys.tv.size(); i++)
//        file<<Sys.tv[i]<<" "<<opv[i]<<std::endl;
//    file.close();
    return 0;
}
