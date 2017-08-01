#include <iostream>
#include <fstream>
#include <complex>
#include "itensor/all.h"
#include "itebd.h"
using namespace itensor;


int main ( int argc, char *argv[] ) {
    double h;
    unsigned z(2); //coordination number
    SpinHalf sites(2);

    {
        int Chi=50;
        std::string output_file("fig5.dat");
        h = 1;
        auto ampo = AutoMPO(sites);
        ampo += 4, "Sx", 1, "Sx", 2;
        ampo += -h * 2 / z, "Sz", 1;
        ampo += -h * 2 / z, "Sz", 2;
        itebd<Index> Sys(ampo);
        Sys.randomize();
        for (double e = 0; e > -5.5; e -= 0.1) {
            double dt = std::pow(10, e);
            std::cout << dt << "  " << Sys.step_imag(dt, 20, 0, Chi) << "  " << Sys.bonddim << std::endl;
        }

        auto mz = 1 * sites.op("Sz", 1) * sites.op("Id", 2) + 1 * sites.op("Id", 1) * sites.op("Sz", 2);
        std::vector<int> rv(300);
        std::generate(rv.begin(), rv.end(), [] {static int i{1}; return i++;});
        auto Sz = std::real(Sys.measure(mz));
        std::cout << Sz << std::endl;
        auto Cv = Sys.measure(2 * sites.op("Sz", 1), 2 * sites.op("Sz", 2), rv);

        std::ofstream file;
        file.open(output_file);
        for (auto i = 0; i < rv.size(); i++) {
            std::cout << rv[i] << ": " << std::real(Cv[i]) - Sz * Sz << std::endl;
            file << rv[i] << " " << std::real(Cv[i]) - Sz * Sz << std::endl;
        }
        file.close();
    }


    {
        std::string output_file("fig6.dat");
        h = 10;
        auto ampo = AutoMPO(sites);
        ampo += 4, "Sx", 1, "Sx", 2;
        ampo += -h * 2 / z, "Sz", 1;
        ampo += -h * 2 / z, "Sz", 2;
        itebd<Index> Sys(ampo);// Sys.randomize();
        for (int e = 0; e > -5; e -= 1) std::cout << Sys.step_imag(std::pow(10, e), 500) << std::endl;
        h = 3;
        double T = 10, dt = .01, thres = 0, M = 100;
        ampo = AutoMPO(sites);
        ampo += 4, "Sx", 1, "Sx", 2;
        ampo += -h * 2 / z, "Sz", 1;
        ampo += -h * 2 / z, "Sz", 2;
        auto mz = 0.5 * sites.op("Sz", 1) * sites.op("Id", 2) + 0.5 * sites.op("Id", 1) * sites.op("Sz", 2);
        std::vector<double> mzv(1, std::real(Sys.measure(mz)));
        Sys.resetTime();
        Sys.setH(ampo);

        while (Sys.tv.back() < T) {
            Sys.step(dt, 1, thres, M);
            mzv.push_back(std::real(Sys.measure(mz)));
            std::cout << Sys.tv.back() << ":  " << Sys.bonddim << "  " << mzv.back() << std::endl;
        }
        std::ofstream file;
        file.open(output_file);
        for (auto i = 0; i < Sys.tv.size(); i++)
            file << Sys.tv[i] << " " << mzv[i] << std::endl;
        file.close();
    }
    return 0;
}
