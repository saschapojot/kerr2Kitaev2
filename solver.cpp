//
// Created by tusa on 19/12/20.
//

#include "solver.hpp"
solver::solver(const CONSTS &conTmp) {
    //params
    this->CON = conTmp;
    std::vector<double> tmp(this->CON.Q + 1, 0.0);//q=0,1,...,Q
    //init thetaDTab, thetaTotTab, thetaGTab
    for (int i = 0; i < this->CON.N + 1; i++) {

        this->thetaDTab.push_back(tmp);
        this->thetaTotTab.push_back(tmp);
        this->thetaGTab.push_back(tmp);
    }
    std::vector<double> tmp1(this->CON.N,0.0);//k=0,1,...,N-1;
    //init beta
    for(int j=0;j<this->CON.Q+1;j++){
        this->beta.push_back(tmp1);
    }


}


double solver::b0(const int &k) {
    return this->CON.d0 * std::sin((double) k * this->CON.dk);
}

double solver::c0(const int &k) {

    return this->CON.t0 * std::cos((double) k * this->CON.dk) + this->CON.mu0 / 2;

}

///
/// \param k
/// \return init vector for each k

Eigen::Vector2cd solver::initVec(const int &k) {

    Eigen::Vector2cd rst;
    double b0Val = this->b0(k);
    double c0Val = this->c0(k);

    double denom2 = 2 * b0Val * b0Val + 2 * c0Val * c0Val - 2 * c0Val * std::sqrt(b0Val * b0Val + c0Val * c0Val);
    if (std::abs(denom2) >= this->CON.tol) {
        double denom = std::sqrt(denom2);
        std::complex<double> denomZ{denom, 0};

        std::complex<double> numer1{0, b0Val};
        std::complex<double> numer2{c0Val - std::sqrt(c0Val * c0Val + b0Val * b0Val), 0};
        rst[0] = numer1 / denomZ;
        rst[1] = numer2 / denomZ;
        return rst;


    } else {

        if (c0Val > 0) {
            rst[0] = std::complex<double>(0, 1);
            rst[1] = std::complex<double>(0, 0);
            return rst;
        } else {
            rst[0] = std::complex<double>(0, 0);
            rst[1] = std::complex<double>(1, 0);
            return rst;

        }

    }

}


///
/// \param k
/// \return linear part of the Hamiltonian after the quench
Eigen::Matrix2cd solver::H0(const int &k) {
    Eigen::Matrix2cd rst;


    rst(0, 0) = std::complex<double>(-this->CON.mu1 - this->CON.t1 * std::cos((double) k * this->CON.dk), 0);

    double tmp01 = this->CON.d1 * std::sin(this->CON.dk * (double) k);

    rst(0, 1) = std::complex<double>(0, tmp01);

    rst(1, 0) = std::complex<double>(0, -tmp01);
    rst(1, 1) = std::complex<double>(this->CON.t1 * std::cos((double) k * this->CON.dk), 0);

    return rst;


}

Eigen::Matrix2cd solver::expH0(const int &k, const double &stept) {
    Eigen::Matrix2cd h0Val = this->H0(k);
    Eigen::Matrix2cd z0 = -1.0j * stept / 2.0 * h0Val;
    return z0.exp();

}

///
/// \param k
/// \param vecStart
/// \param stept
/// \return state vector after 1 step of 2nd Strang splitting S2
Eigen::Vector2cd solver::S2(const int &k, const Eigen::Vector2cd &vecStart, const double &stept) {
//step 1
    Eigen::Matrix2cd exph0Val = this->expH0(k, stept);
    Eigen::Vector2cd vec1 = exph0Val * vecStart;

//step 2
    std::complex<double> vkm = vec1[0];
    std::complex<double> wkm = vec1[1];
    //calculate etakm
    std::complex<double> tmp1 = -1.0j * this->CON.lmd * std::pow(std::abs(vkm), 2) * stept;
    std::complex<double> etakm = vkm * std::exp(tmp1);
    // calculate zetakm
    std::complex<double> tmp2 = -1.0j * this->CON.lmd * std::pow(std::abs(wkm), 2) * stept;
    std::complex<double> zetakm = wkm * std::exp(tmp2);

    //step 3
    Eigen::Vector2cd vec2;
    vec2 << etakm, zetakm;
    return exph0Val * vec2;


}


///
/// \param k
/// \return state vector [z1, z2] at momentum k, at each time t
std::vector<Eigen::Vector2cd> solver::calculateOneRow(const int &k) {
    std::vector<Eigen::Vector2cd> rst;

    //initialization
    Eigen::Vector2cd veck0 = this->initVec(k);
    rst.push_back(veck0);

    //time step number: 0,1,...,Q*R-1
    for (int m = 0; m < this->CON.Q * this->CON.R; m++) {
        auto vecCurr = rst.back();
        //Use 2nd order
        auto vecNext = this->S2(k, vecCurr, this->CON.dt);
        rst.push_back(vecNext);
    }

    return rst;


}


std::complex<double> solver::Jkab(const std::vector<Eigen::Vector2cd> &veck, const int &k, const int &a,
                                  const int &b) {

    //subinterval  starting num a=0,1,...,Q-1, within each subinterval b=0,1,...,R


    Eigen::Vector2cd vec=veck[a*this->CON.R+b];
    std::complex<double> yk = vec[0];
    std::complex<double> zk = vec[1];

    Eigen::Matrix2cd H1;
    H1(0, 0) = std::complex<double>(-this->CON.mu1 - this->CON.t1 * std::cos(this->CON.dk * (double) k) +
                                    this->CON.lmd * std::pow(std::abs(yk), 2), 0);
    double tmp01 = this->CON.d1 * std::sin(this->CON.dk * (double) k);
    H1(0, 1) = std::complex<double>(0, tmp01);
    H1(1, 0) = std::complex<double>(0, -tmp01);
    H1(1, 1) = std::complex<double>(
            this->CON.t1 * std::cos(this->CON.dk * (double) k) + this->CON.lmd * std::pow(std::abs(zk), 2), 0);

    std::complex<double> rst;
    rst = vec.adjoint() * H1 * vec;
    rst /= (std::pow(std::abs(yk), 2) + std::pow(std::abs(zk), 2));
    return rst;



}


///
/// \param veck
/// \param k
/// \param a
/// \return Simpson integration
std::complex<double> solver::simpsonD(const std::vector<Eigen::Vector2cd> &veck, const int &k, const int &a) {
    ///a=0,1,...,Q-1
    //integration over [ads, (a+1)ds]
    std::complex<double> evenSum, oddSum;
    oddSum=std::complex<double>(0,0);
    evenSum=std::complex<double>(0,0);

    //compute odd sums
    for(int b=1;b<this->CON.R;b+=2){
        oddSum+=this->Jkab(veck,k,a,b);
    }

    //compute even sums
    for(int b=2;b<this->CON.R;b+=2){
        evenSum+=this->Jkab(veck,k,a,b);

    }


    std::complex<double> rst=this->Jkab(veck,k,a,0)+4.0*oddSum+2.0*evenSum+this->Jkab(veck,k,a,this->CON.R);
    rst*=this->CON.dt/3.0;
    return  rst;
}


std::vector<std::complex<double>> solver::simpTabOneRow(const std::vector<Eigen::Vector2cd> &veck, const int &k) {
    //a=0,1,...,Q-1 th subinterval
    std::vector<std::complex<double>> rst;
    for (int a = 0; a < this->CON.Q; a++) {
        rst.push_back(this->simpsonD(veck, k, a));
    }

    return rst;


}

std::vector<std::complex<double>> solver::simpTabOneRowSum(const std::vector<Eigen::Vector2cd> &veck, const int &k) {
    std::vector<std::complex<double>> simpOneRowVals = this->simpTabOneRow(veck,
                                                                           k);//elem number: 0, 1, ..., Q-1 th subinterval

    std::vector<std::complex<double>> rst;//elem number: 0, 1, ..., Q
    std::complex<double> tmp = std::complex<double>(0, 0);
    rst.push_back(tmp);
    for (const auto &elem:simpOneRowVals) {
        tmp += elem;
        rst.push_back(tmp);
    }

    return rst;
}

void solver::thetaDTabOneRow(const std::vector<Eigen::Vector2cd> &veck, const int &k) {


    std::vector<std::complex<double>> oneRowSum = this->simpTabOneRowSum(veck, k);
//q=0,1,...,Q
    Eigen::Vector2cd vec0 = veck[0];
    std::complex<double> tmp2D = vec0.adjoint() * vec0;
    for (int q = 0; q < this->CON.Q + 1; q++) {
        std::complex<double> tmpSum = oneRowSum[q];
        tmpSum *= -1;
        Eigen::Vector2cd vecqR = veck[q * this->CON.R];
        std::complex<double> tmp2N = (vecqR.adjoint() * vecqR);
        std::complex<double> term2 = std::complex<double>(0, 0.5) * std::log(tmp2N / tmp2D);
        double tmpThetaD = (tmpSum + term2).real();
        this->thetaDTab[k][q] = tmpThetaD;

    }

}

void solver::thetaTotTabOneRow(const std::vector<Eigen::Vector2cd> &veck, const int &k) {

    //k=0,1,...,N
    //q=0,1,...,Q
    Eigen::Vector2cd vec0 = veck[0];

    for (int q = 0; q < this->CON.Q + 1; q++) {
        Eigen::Vector2cd vecqR = veck[q * this->CON.R];
        std::complex<double> tmp = vec0.adjoint() * vecqR;
        std::complex<double> rst = std::complex<double>(0, -1) * std::log(tmp / std::abs(tmp));
        this->thetaTotTab[k][q] = rst.real();


    }


}

void solver::thetaGOneRow(const std::vector<Eigen::Vector2cd> &veck, const int &k) {
    this->thetaTotTabOneRow(veck, k);
    this->thetaDTabOneRow(veck, k);

    for (int q = 0; q < this->CON.Q + 1; q++) {
        this->thetaGTab[k][q] = this->thetaTotTab[k][q] - this->thetaDTab[k][q];
    }

}


std::vector<std::vector<int>> solver::partitionOfKIndices(const std::vector<int> &kinds) {

    int latticeSize = (int) kinds.size();
    int firstPartSize = std::floor(latticeSize / this->CON.threadNum);
    int leftNum = latticeSize - (this->CON.threadNum - 1) * firstPartSize;
    std::vector<std::vector<int>> rst;
    for (int i = 0; i < this->CON.threadNum - 1; i++) {
        std::vector<int> tmp0;
        int startTmp = i * firstPartSize;
        for (int j = 0; j < firstPartSize; j++) {
            tmp0.push_back(startTmp + j);
        }
        rst.push_back(tmp0);
    }
    int startLastPart = (this->CON.threadNum - 1) * firstPartSize;
    std::vector<int> tmp1;
    for (int j = 0; j < leftNum; j++) {
        tmp1.push_back(startLastPart + j);
    }

    if (tmp1.size() > 0) { rst.push_back(tmp1); }
    return rst;


}

void solver::worker(const std::vector<int> &onePart) {
    for (const auto &kVal:onePart) {
        std::vector<Eigen::Vector2cd> solutionOneRowTmp = this->calculateOneRow(kVal);
        this->thetaGOneRow(solutionOneRowTmp, kVal);
    }


}

void solver::assignWork() {
    std::vector<std::vector<int>> kindPartsAll = this->partitionOfKIndices(this->CON.kIndAll);
    std::vector<std::thread> thrds;
    for (const auto &kVecTmp:kindPartsAll) {
        thrds.emplace_back(&solver::worker, this, kVecTmp);
    }

    for (auto &th:thrds) {
        th.join();
    }

}

double solver::jumpDecision(const double &incr) {
    double tmp = incr / M_PI;
    if (tmp >= this->CON.cutOff) {
        return incr - 2 * M_PI;
    } else if (tmp <= -this->CON.cutOff) {
        return incr + 2 * M_PI;
    } else {
        return incr;
    }


}


void solver::calcIncrOfThetaG() {
    for (int q = 0; q < this->CON.Q + 1; q++) {
        for (int k = 0; k < this->CON.N; k++) {
            this->beta[q][k] = this->jumpDecision(this->thetaGTab[k+1][q] - this->thetaGTab[k][q]);
        }
    }

}

void solver::writeBeta() {
    std::string outFileName = this->CON.dir + "beta";
    outFileName += "mu0" + boost::lexical_cast<std::string>(this->CON.mu0)
                   + "t0" + boost::lexical_cast<std::string>(this->CON.t0)
                   + "d0" + boost::lexical_cast<std::string>(this->CON.d0)
                   + "mu1" + boost::lexical_cast<std::string>(this->CON.mu1)
                   + "t1" + boost::lexical_cast<std::string>(this->CON.t1)
                   + "d1" + boost::lexical_cast<std::string>(this->CON.d1)
                   + "lmd" + boost::lexical_cast<std::string>(this->CON.lmd) + ".csv";

    std::ofstream ofPtr;
    ofPtr.open(outFileName);
    for (const auto &vec:this->beta) {
        size_t colN = vec.size();
        for (auto i = 0; i < colN - 1; i++) {
            ofPtr << vec[i] << ",";
        }
        ofPtr << vec[colN - 1] << "\n";
    }
    ofPtr.close();


}
void solver::writeThetaGTab() {

    std::string outFileName = this->CON.dir + "G" + "mu0" + boost::lexical_cast<std::string>(this->CON.mu0)
                              + "t0" + boost::lexical_cast<std::string>(this->CON.t0)
                              + "d0" + boost::lexical_cast<std::string>(this->CON.d0)
                              + "mu1" + boost::lexical_cast<std::string>(this->CON.mu1)
                              + "t1" + boost::lexical_cast<std::string>(this->CON.t1)
                              + "d1" + boost::lexical_cast<std::string>(this->CON.d1)
                              + "lmd" + boost::lexical_cast<std::string>(this->CON.lmd) + ".csv";
    std::ofstream ofPtr;
    ofPtr.open(outFileName);
    for (const auto &vec:this->thetaGTab) {
        size_t colN = vec.size();
        for (auto i = 0; i < colN - 1; i++) {
            ofPtr << vec[i] << ",";
        }
        ofPtr << vec[colN - 1] << "\n";
    }
    ofPtr.close();
}
void solver::runCalc() {

    this->assignWork();
    this->calcIncrOfThetaG();

    this->writeBeta();
    this->writeThetaGTab();
}
