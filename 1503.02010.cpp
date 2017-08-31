#include <iostream>
#include <fstream>
#include "itensor/all.h"
#include "itebd.h"
using namespace itensor;


int main ( int argc, char *argv[] ) {
    std::string output_file("fig1_U=4.dat");
    double thres=1E-8,M=300, dt=0.002, T=3;
    double U(4);
    unsigned z(2); //coordination number

    //2-site Hamiltonian
    Hubbard sites(2);
    auto s1=sites(1), s2=sites(2);
    AutoMPO ampo(sites);
    ampo += -1,"Cdagup",1,"Cup",2;
    ampo += -1,"Cdagup",2,"Cup",1;
    ampo += -1,"Cdagdn",1,"Cdn",2;
    ampo += -1,"Cdagdn",2,"Cdn",1;
    ampo += U/z,"Nupdn",1;
    ampo += U/z,"Nupdn",2;
    itebd<IQIndex> Sys(ampo,2,3);

    //observable (doublon)
    auto D=0.5*(sites.op("Nupdn",1)*sites.op("Id",2)+sites.op("Id",1)*sites.op("Nupdn",2));
    auto mz=0.5*(sites.op("Sz",1)*sites.op("Id",2)-sites.op("Id",1)*sites.op("Sz",2));
    std::vector<double> Dv(1,std::real(Sys.measure(D))), mzv(1,std::real(Sys.measure(mz)));

    //evolution
    while (Sys.tv.back()<T)  {
        Sys.step(dt,5,thres,M);
        Dv.push_back(std::real(Sys.measure(D)));
        mzv.push_back(std::real(Sys.measure(mz)));
        std::cout<<Sys.tv.back()<<":  "<<Sys.bonddim<<"  "<<Dv.back()<<"  "<<mzv.back()<<std::endl;
    }

    //output to file
    std::ofstream file;
    file.open(output_file);
    for(auto i=0; i<Sys.tv.size(); i++)
        file<<Sys.tv[i]<<" "<<Dv[i]<<" "<<mzv[i]<<std::endl;
    file.close();

    return 0;
}
