#include <algorithm>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <complex>
#include <iostream>
#include "itensor/all.h"
#pragma once

namespace itensor {

    template<typename I, unsigned z=2>
    struct itebd {
        using T = ITensorT<I>;

        itebd(const AutoMPO &ampo, int QNA=1, int QNB=1, int Chi0=1);

        double step(std::complex<double> dt, size_t steps = 1, double thres = 1E-10, int maxm = 0);

        double step_imag(std::complex<double> dt, size_t steps = 1, double thres = 1E-10, int maxm = 0);

        ITensorT<I> measure();

        std::complex<double> measure(const T &op);

        std::vector<std::complex<double>> measure(const T &, const T &,  const std::vector<int> &rv);

        void randomize();

        void resetTime();

        void setH(const AutoMPO &ampo);

        void backup();
        void restore();

        T H,U, swapAtoB, swapBtoA, GA,GAback, GB,GBback, psi,psiback;
        std::complex<double> dt_used_for_U;
        std::vector<T> L,Lback;
        AutoMPO ampo;
        IQIndex sA, sB;
        std::vector<double> tv, Sv;
        double error;
        int bonddim;
    private:

        double add_time_entropy(std::complex<double> dt = 0);
    };

};

namespace itensor {
/*----------------------------Constructors----------------------------------*/
    template<typename I,unsigned z>
    itebd<I,z>::itebd(const AutoMPO &ampo, int QNA, int QNB, int Chi0):
            ampo(ampo), sA(ampo.sites()(1)), sB(ampo.sites()(2))
    {
        using IVtype = typename I::indexval_type;
        std::vector<std::vector<IndexQN>> iqA(z), iqB(z);
        for (auto e : sA) {
            for (unsigned i = 0; i < z; i++) {
                iqA[i].emplace_back(Index(e.index.name(), Chi0), e.qn);
                iqB[i].emplace_back(Index(e.index.name(), Chi0), e.qn);
            }
        }
        std::vector<IQIndex> A(z), B(z);
        auto GAinds = std::vector<IQIndex>(z + 1);
        auto GBinds = std::vector<IQIndex>(z + 1);
        L=std::vector<T>(z);
        for (unsigned i = 0; i < z; i++) {
            A[i] = IQIndex("A", std::move(iqA[i]));
            B[i] = IQIndex("B", std::move(iqB[i]));
            L[i] = IQTensor(A[i], B[i]);  //automatically converted to ITensor if I=Index
            GAinds[i] = dag(A[i]);
            GBinds[i] = dag(B[i]);
        }
        GAinds[z] = sA;
        GBinds[z] = sB;
        GA = IQTensor(GAinds);
        GB = IQTensor(GBinds);

        std::vector<IVtype> GAindvals = {IVtype(sA,QNA),IVtype(A[0],QNA),IVtype(A[1],QNB)};
        std::vector<IVtype> GBindvals = {IVtype(sB,QNB),IVtype(B[0],QNA),IVtype(B[1],QNB)};
        L[0].set(A[0](QNA), B[0](QNA), 1);
        L[1].set(A[1](QNB), B[1](QNB), 1);
        for (unsigned i = 2; i < z; i++) {
            GAindvals.push_back(IVtype(A[i],1));
            GBindvals.push_back(IVtype(B[i],1));
            L[i].set(A[i](1), B[i](1), 1);
        }
        GA.set(GAindvals,1);
        GB.set(GBindvals,1);

        // if (z == 2) {
        //     GA.set(sA(QNA), A[0](QNA), A[1](QNB), 1);
        //     GB.set(sB(QNB), B[0](QNA), B[1](QNB), 1);
        //     L[0].set(A[0](QNA), B[0](QNA), 1);
        //     L[1].set(A[1](QNB), B[1](QNB), 1);
        // } else if (z == 3) {
        //     GA.set(sA(QNA), A[0](QNA), A[1](1), A[2](QNB), 1);
        //     GB.set(sB(QNB), B[0](QNA), B[1](1), B[2](QNB), 1);
        //     L[0].set(A[0](QNA), B[0](QNA), 1);
        //     L[1].set(A[1](1), B[1](1), 1);
        //     L[2].set(A[2](QNB), B[2](QNB), 1);
        // }

        auto Hmpo = toMPO<T>(ampo);
        H = Hmpo.A(1) * Hmpo.A(2);
        dt_used_for_U = 0;
        swapAtoB = delta(dag(sA), sB);
        swapBtoA = delta(sA, dag(sB));
        add_time_entropy();
    }

/*----------------------------Stepper----------------------------------*/
    template<typename I,unsigned z>
    double itebd<I,z>::step(std::complex<double> dt, size_t steps, double thres, int maxm) {
        // auto Uev = toExpH<T>(ampo, dt * Cplx_i);
        // auto U = Uev.A(1) * Uev.A(2);
        if (dt_used_for_U!=dt){
            U = expHermitian(H,-dt * Cplx_i);
            dt_used_for_U = dt;
        }
        T Upsi;

        auto invL=std::vector<T>(z);
        for(unsigned i=1;i<z;i++){ //as L[0] is recomputed immediately below
            auto Ai = commonIndex(L[i % z], GA),
                 Bi = commonIndex(L[i % z], GB);
            invL[i]=T(Ai, Bi);
            for (int comp = 1; comp <= Ai.m(); comp++)
                if (L[i].real(Ai(comp), Bi(comp)) != 0)
                    invL[i].set(Ai(comp), Bi(comp), 1. / L[i].real(Ai(comp), Bi(comp)));
        }

        Args svdargs;
        if (thres) svdargs.add("Cutoff",thres);
        if (maxm) svdargs.add("Maxm",maxm);
        std::complex<double> E;
        for (size_t run = 0; run < steps; run++) {
            E = 0;
            for (unsigned sh = 0; sh < z; sh++) {
                if (sh % 2) {
                    GA *= swapAtoB;
                    GB *= swapBtoA;
                }
                psi = GA * L[sh] * GB;
                for (unsigned i = sh + 1; i < sh + z; i++) {
                    auto Bi = commonIndex(L[i % z], GB);
                    psi = L[i % z] * prime(psi, Bi) * prime(L[i % z], Bi);
                }
                Upsi = psi * U;
                Upsi.noprime(Site);

                auto GAinds = std::vector<I>(z);
                GAinds[0] = (sh % 2) ? sB : sA;
                for (unsigned i = 1; i < z; i++) GAinds[i] = commonIndex(L[(sh + i) % z], GB);
                GA = T(GAinds);
                GB = T();
                svd(Upsi, GA, L[sh], GB, svdargs);
                L[sh] /= norm(L[sh]);

                auto Ai = commonIndex(L[sh], GA),
                     Bi = commonIndex(L[sh], GB);
                bonddim = Ai.m();
                invL[sh]=T(Ai, Bi);
                for (int comp = 1; comp <= Ai.m(); comp++)
                    if (L[sh].real(Ai(comp), Bi(comp)) != 0)
                        invL[sh].set(Ai(comp), Bi(comp), 1. / L[sh].real(Ai(comp), Bi(comp)));

                for (unsigned i = sh + 1; i < sh + z; i++) {
                    GA *= dag(invL[i%z]);
                    GB *= dag(invL[i%z]);
                }
                if (sh % 2) {
                    GA *= swapBtoA;
                    GB *= swapAtoB;
                }
                auto psiUpsi = Upsi * dag(psi);
                auto psipsi = psi * dag(psi);
                E += (std::log(psiUpsi.cplx()) / (-dt * Cplx_i)) / 2.0 / psipsi.cplx();
            }
        }
        add_time_entropy(static_cast<double>(steps) * dt);
        return std::real(E);
    }

    template<typename I,unsigned z>
    double itebd<I,z>::step_imag(std::complex<double> dt, size_t steps, double thres, int maxm) {
        return itebd<I,z>::step(-Cplx_i * dt, steps, thres, maxm);
    }


/*----------------------------Helpers----------------------------------*/
    template<typename I,unsigned z>
    double itebd<I,z>::add_time_entropy(std::complex<double> dt) {
        auto A0 = commonIndex(L[0], GA);
        auto B0 = commonIndex(L[0], GB);
        if (tv.empty()) tv.push_back(0);
        else tv.push_back(tv.back() + std::abs(dt));
        double en = 0;
        for (int i = 1; i <= A0.m(); i++) {
            double eig = L[0].real(A0(i), B0(i));
            eig = eig * eig;
            if (eig > 0 && eig < 1) en += -eig * std::log(eig);
        }
        Sv.push_back(en);
        return en;
    }


    template<typename I,unsigned z>
    void itebd<I,z>::randomize() {
        itensor::randomize(GA);
        itensor::randomize(GB);
    }
    template<typename I,unsigned z>
    void itebd<I,z>::resetTime() {
        tv=std::vector<double>();
        Sv=std::vector<double>();
        add_time_entropy();
    }
    template<typename I,unsigned z>
    void itebd<I,z>::backup() {
        Lback = std::vector<T>(L);
        GAback=GA;
        GBback=GB;
        psiback=psi;
    }
    template<typename I,unsigned z>
    void itebd<I,z>::restore() {
        L = std::vector<T>(Lback);
        GA=GAback;
        GB=GBback;
        psi=psiback;
    }
    template<typename I,unsigned z>
    void itebd<I,z>::setH(const AutoMPO &in) {
        ampo=in; 
        auto Hmpo = toMPO<T>(ampo);
        H = Hmpo.A(1) * Hmpo.A(2);
        dt_used_for_U = 0;
    }

    template<typename I,unsigned z>
    ITensorT<I> itebd<I,z>::measure() {
        psi = GA * L[0] * GB;
        for (unsigned i = 1; i < z; i++) {
            auto Bi = commonIndex(L[i], GB);
            psi = L[i] * prime(psi, Bi) * prime(L[i], Bi);
        }
        return psi;
    }

    template<typename I,unsigned z>
    std::complex<double> itebd<I,z>::measure(const T &op) {
        measure();
        auto num = noprime(psi * op) * dag(psi), den = noprime(psi) * dag(psi);
//        std::cout<<den.cplx()<<std::endl;
        return num.cplx() / den.cplx();
    }

    template<typename I,unsigned z>
    std::vector<std::complex<double>> itebd<I,z>::measure(const T &op1A, const T &op2, const std::vector<int> &rv){
        std::vector<std::complex<double>> Cv(rv.size());
        T op2A,op2B;
        if (hasindex(op2,sA)){
            op2A = op2;
            op2B = op2*swapBtoA*prime(dag(swapBtoA));
        } else{
            op2B = op2;
            op2A = op2*swapAtoB*prime(dag(swapAtoB));
        }
        auto B0=commonIndex(L[0],GB);
        auto A1=commonIndex(L[1],GA);
        auto GAL0 = GA*L[0],
             GBL1 = GB*L[1];
//        auto L0block = GAL0*L[1]*prime(dag(L[1]),A1)*prime(dag(GAL0),A1,B0);
        auto L1block = GAL0*L[1]*prime(dag(L[1]),A1)*op1A*prime(dag(GAL0));
//        auto R0block = GAL0*prime(dag(GAL0),A1);
//             R0block = R0block*GBL1*prime(dag(GBL1),A1,B0);
        auto R1block = GAL0*prime(dag(GAL0),A1);
             R1block = R1block*GBL1*op2B*prime(dag(GBL1));
        auto R2block = GAL0*op2A*prime(dag(GAL0),A1,sA);
             R2block = R2block*GBL1*prime(dag(GBL1),A1,B0);
        int Lblocksize=0;
        auto B1=commonIndex(L[1],GB);
        auto A0=commonIndex(L[0],GA);
        for (int ri=0;ri<rv.size();ri++){
            while (Lblocksize<(rv[ri]+1)/2-1){
//                L0block=L0block*GBL1*prime(dag(GBL1),B0,A1);
//                L0block=L0block*GAL0*prime(dag(GAL0),A1,B0);
                L1block=L1block*GBL1*prime(dag(GBL1),B0,A1);
                L1block=L1block*GAL0*prime(dag(GAL0),A1,B0);
                ++Lblocksize;
            }
            T num;
            if (rv[ri]%2) num=L1block*R1block;
            else num = L1block*R2block;
            Cv[ri]=num.cplx();
        }
        return Cv;
    }
}
