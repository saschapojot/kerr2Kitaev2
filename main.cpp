#include "solver.hpp"
#include <chrono>
#include <iostream>
int main() {
//    CONSTS con1;
//    con1.lmd=6;
//    con1.mu0=-6;
//    con1.mu1=-6;
//    con1.setN(1000);
//    con1.setDk();
//    con1.initKInd();
//    con1.initDir();
//    solver sol0=solver(con1);
//    auto start=std::chrono::high_resolution_clock::now();
//    sol0.runCalc();
//    auto stop=std::chrono::high_resolution_clock::now();
//    auto duration=std::chrono::duration_cast<std::chrono::microseconds>(stop-start);
//    std::cout<<duration.count()/1e6<<" seconds"<<std::endl;

std::vector<int>NVec={250,500,1000,1500,3000,6000,12000,24000,48000,96000};

for(const int &n:NVec){
    CONSTS conTmp;
    conTmp.setN(n);
    conTmp.setDk();
    conTmp.initKInd();
    conTmp.initDir();
    solver sol0=solver(conTmp);
        auto start=std::chrono::high_resolution_clock::now();
    sol0.runCalc();
    auto stop=std::chrono::high_resolution_clock::now();
    auto duration=std::chrono::duration_cast<std::chrono::microseconds>(stop-start);
    std::cout<<duration.count()/1e6<<" seconds"<<std::endl;
}




    return 0;
}    
