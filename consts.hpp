//
// Created by tusa on 19/12/20.
//

#ifndef KERR2KITAEV2_CONSTS_HPP
#define KERR2KITAEV2_CONSTS_HPP
#include <boost/math/constants/constants.hpp>
#include <vector>
#include <complex>
#include <thread>
#include <string>
class CONSTS {


public:
    CONSTS() {




        int tmp1 = (int) std::thread::hardware_concurrency();
        int tmp2 = 12;
        this->threadNum = tmp1 > tmp2 ? tmp1 : tmp2;
    }

public:
    void setN(const int&nVal){
        this->N=nVal;
    }
    void setDk(){
        this->dk=  2 * M_PI / ((double) N);
    }
    void initKInd() {
        for (int i = 0; i < this->N + 1; i++) {
            kIndAll.push_back(i);
        }
    }
    void initDir(){
        this->dir+=boost::lexical_cast<std::string>(N);

    }

public:
    int N ;//lattice number
    double dk;

    double lmd = 0;


    /*
     * Indices of momentum space, 0,1,...,N-1,N
     * */
    std::vector<int> kIndAll;

    /*
     * Parameters before the quench
     * */
    double mu0 = 0;
    double t0 = 1.0;
    double d0 = -1.0;

    /*
     * Parameters after the quench
     * */
    double mu1 = -3;
    double t1 = t0;
    double d1 = d0;
    /*
     * consts for time*/

    int R = 160;//small time step number

    int Q = 1000;//large time step number
    double dt = 0.000125;//small time step
    double ds = (double) R * dt;//large time step


    double tol = 1e-16;
    double cutOff = 1.4;

    int threadNum;//= std::thread::hardware_concurrency();

    std::string dir = "/home/disk2/Documents/cppCode/kerr2Kitaev2/testkNum25/kNum" ;


};


#endif //KERR2KITAEV2_CONSTS_HPP
