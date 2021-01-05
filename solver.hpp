//
// Created by tusa on 19/12/20.
//

#ifndef KERR2KITAEV2_SOLVER_HPP
#define KERR2KITAEV2_SOLVER_HPP
#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/MatrixFunctions>
#include <map>
#include "consts.hpp"
#include <complex>
#include <fstream>
#include <cmath>
#include <utility>
//#include "matplotlib-cpp/matplotlibcpp.h"
#include <memory>
#include <cmath>
//namespace plt=matplotlibcpp;
#include <iostream>

class solver{

public:
    solver(const CONSTS&);

    /*
     * Functions to initialize
     * */
    double b0(const int &k);
    double c0(const int &k);

    Eigen::Vector2cd initVec(const int&k);

    Eigen::Matrix2cd H0(const int &k);//linear part of the Hamiltonian

    Eigen::Matrix2cd expH0(const int&k, const double & stept);//exponentiates the H0 matrix in S2
    Eigen::Vector2cd S2(const int&k, const Eigen::Vector2cd& vecStart, const double&stept);//one step S2

    //0th
    std::vector<Eigen::Vector2cd> calculateOneRow(const int &k);//calculates state vectors for one k

    //1st
    std::complex<double> Jkab(const std::vector<Eigen::Vector2cd> & veck, const int &k, const int &a, const int &b);
    std::complex<double> simpsonD(const std::vector<Eigen::Vector2cd> & veck, const int &k, const int &a);
    std::vector<std::complex<double>> simpTabOneRow(const std::vector<Eigen::Vector2cd> &veck, const int &k);
    //2nd

    std::vector<std::complex<double>> simpTabOneRowSum(const  std::vector<Eigen::Vector2cd> &veck, const int &k);
    void thetaDTabOneRow(const  std::vector<Eigen::Vector2cd> &veck, const int &k);

    //3rd
    void thetaTotTabOneRow(const std::vector<Eigen::Vector2cd> &veck, const int &k);

    //4th
    void thetaGOneRow(const std::vector<Eigen::Vector2cd> &veck, const int &k);

    //5th
//
  std::vector<std::vector<int>> partitionOfKIndices(const std::vector<int> & kinds);

    //6th
    void worker(const std::vector<int> & onePart);

    //7th
    void assignWork();


    //8th
    double jumpDecision(const double&incr);
    void calcIncrOfThetaG();
    //9th
    void writeBeta();
    void writeThetaGTab();


    //10th
    void runCalc();



public:
    CONSTS CON;
    std::vector<std::vector<double>>thetaDTab;//K=0,1,...,N-1,N; q=0,1,...,Q
    std::vector<std::vector<double>>thetaTotTab;//k=0,1,...,N; q=0,1,..., Q;
    std::vector<std::vector<double>> thetaGTab;//k=0,1,...,N; q=0,1,...,Q;
    std::vector<std::vector<double>> beta;//q=0,1,...,Q; k=0,1,...,N-1(increment at each k);
    std::vector<double>W;//q=0,1,...,Q;

    std::vector<double>rateFunction;//q=0,1,...,Q;
};










#endif //KERR2KITAEV2_SOLVER_HPP
