#include <iostream>
#include <fstream>
#include <complex>
#include "itensor/all.h"
#include "itebd.h"
using namespace itensor;


int main ( int argc, char *argv[] ) {
    double gi,gf;
    unsigned z(2); //coordination number
    SpinHalf sites(2);

    {
        std::string output_file("fig3b.dat");
        gi = 0.4; gf=-0.4;
        int M=100;
        double thres=1E-9,dt=0.01,T=150;
        auto ampo0 = AutoMPO(sites);
        ampo0 +=  -4,"Sz", 1, "Sz", 2;
        ampo0 += (1+gi)*2/z, "Sx", 1;
        ampo0 += (1+gi)*2/z, "Sx", 2;
        itebd<Index> Sys(ampo0);
        Sys.randomize();
        for (double e = 0; e > -5.5; e -= 0.1) {
            double dt = std::pow(10, e);
            std::cout << dt << "  " << Sys.step_imag(dt, 20, 0, M) << "  " << Sys.bonddim << std::endl;
        }
        Sys.backup();


        std::vector<decltype(Sys.Sv)> Entropies;
        std::vector<double> Gammas={0.01,0.02,0.03};
        for (auto Gamma : Gammas) {
            Sys.restore();
            Sys.resetTime();
            while (Sys.tv.back() < T) {
                double g = gi - Gamma * Sys.tv.back();
                if (g < gf) g = gf;
                auto ampo = AutoMPO(sites);
                ampo += -4, "Sz", 1, "Sz", 2;
                ampo += (1 + g) * 2 / z, "Sx", 1;
                ampo += (1 + g) * 2 / z, "Sx", 2;
                Sys.setH(ampo);
                Sys.step(dt, 10, thres, M);
                std::cout << Sys.tv.back() << ": " << Sys.Sv.back() << ": " << Sys.bonddim << std::endl;
            }
            Entropies.push_back(Sys.Sv);
        }

        std::ofstream file;
        file.open(output_file);
        for (auto i = 0; i < Sys.tv.size(); i++) {
            file << Sys.tv[i] << " " << Entropies[0][i]<< " " << Entropies[1][i]<< " " << Entropies[2][i]<< std::endl;
        }
        file.close();
    }

    return 0;
}
//
// Created by Hoi Hui on 7/28/17.
//

