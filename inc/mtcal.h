//
// Created by wang- on 25-1-7.
//

#ifndef MTCAL_H
#define MTCAL_H
#include <iostream>
#include <fstream>
#include <array>
#include <cmath>
#include <iomanip>
#include "moe.h"
#include "matcal.cpp"
#include "getbasis.h"
void normalizeBySecondColumn(std::vector<std::vector<double>>& matrix);
void removeRange(std::vector<int> &vec, int i, int k);
std::vector<std::vector<double> > ychangeb(
    const std::vector<std::vector<std::vector<double> > > &ystrall1,
    std::vector<std::vector<std::vector<double> > > &qstrall,
    int k, int rkp, int tnum, std::vector<int> torder, basis rjl);

std::vector<std::vector<double> > ychangep(
    std::vector<std::vector<std::vector<double> > > ystrall1,
    std::vector<std::vector<std::vector<double> > > qstrall,
    int k, int rkp, int tnum, std::vector<int> torder, basis rjl);

std::vector<std::vector<double> > ystrchaa
(std::vector<std::vector<std::vector<double> > > & ystrall1,
std::vector<std::vector<std::vector<double> > > & ystrall2,
int i, int k, int rip,  basis rjlinv, int snum, std::vector<int> sorder,int t);

double qt2term(basis rjl,basis rjlinv,
std::vector<std::vector<std::vector<double> > > ystrall,
std::vector<std::vector<std::vector<double> > > ystrallinv,
int i,int ii,int k, std::vector<int>rkp,std::vector<int>lall,int rip,
std::vector<std::vector<std::vector<double>> > qstrall,std::vector<int>torder,int tnum);

double as2term(basis rjl,basis rjlinv,
std::vector<std::vector<std::vector<double> > > ystrall,
std::vector<std::vector<std::vector<double> > > ystrallinv,
int i,int ii,int k,std::vector<int>lall,int rip,
std::vector<std::vector<std::vector<double>> > qstrall,std::vector<int>sorder,int snum,int t,int contj);

double h0cal(basis rjl, basis rjlinv,
             std::vector<std::vector<std::vector<double> > > ystrall,
             std::vector<std::vector<std::vector<double> > > ystrallinv,
             std::vector<double> energy);

double qtqt(basis rjl, basis rjlinv,
            std::vector<std::vector<std::vector<double> > > ystrall,
            std::vector<std::vector<std::vector<double> > > ystrallinv,
            std::vector<std::vector<std::vector<double> > > qstrall, int tnum,
            std::vector<int> qorder);
double qtqttest(basis rjl, basis rjlinvcha,
            std::vector<std::vector<std::vector<double> > > ystrall,
            std::vector<std::vector<std::vector<double> > > ystrallinv,
            std::vector<std::vector<std::vector<double> > > qstrall, int tnum,
            std::vector<int> qorder);

double atat(basis rjl, basis rjlinv,
            std::vector<std::vector<std::vector<double> > > ystrall,
            std::vector<std::vector<std::vector<double> > > ystrallinv,
            std::vector<std::vector<std::vector<double> > > qstrall, int snum,
            std::vector<int> sorder);
double sigqtterm(basis rjl,basis rjlinv,
std::vector<std::vector<std::vector<double> > > ystrall,
std::vector<std::vector<std::vector<double> > > ystrallinv,
int k,int knum, int rkp,std::vector<int>lall,
std::vector<std::vector<std::vector<double>> > qstrall,std::vector<int>torder,int tnum);
double sigqt(basis rjl,basis rjlinv,
        std::vector<std::vector<std::vector<double> > > ystrall,
        std::vector<std::vector<std::vector<double> > > ystrallinv,
        std::vector<std::vector<std::vector<double> > > qstrall, int tnum,
        std::vector<int> torder);
// std::vector<std::vector<std::vector<double>>>matspl(std::vector<basis>rjl,
//     std::vector<std::vector<std::vector<std::vector<double>>>>ystrall,std::vector<int> torder,
//     std::vector<std::vector<std::vector<double>>>qstrall);

std::vector<std::vector<std::vector<double>>> matspl(const std::vector<basis>& rjl,
    const std::vector<std::vector<std::vector<std::vector<double>>>>& ystrall,
    const std::vector<int>& torder,
    const std::vector<std::vector<std::vector<double>>>& qstrall);
double pvmatcal(basis rjlp,basis rjlpinv,
    int tpnum,std::vector<int> tporder,
    basis rjlv,basis rjlvinv,
    int tvnum,std::vector<int> tvorder,
    int jall,std::vector<std::vector<std::vector<double>>> matsliqn,
    std::vector<std::vector<std::vector<double>>> matsliqp,
    int p1,int p2,int n1,int n2
    );
std::vector<std::vector<double>> matcalsim(std::vector<basis>rjl,
    std::vector<std::vector<std::vector<std::vector<double>>>>ystrall,
    std::vector<int> torder,std::vector<int>sorder,std::vector<double> energy,
    std::vector<std::vector<std::vector<double>>>qstrall,std::vector<double> kt,
    std::vector<std::vector<std::vector<double>>>sstrall,std::vector<double> gs);
std::vector<std::vector<double>> overlapmatcal(std::vector<basis>&rjl,
    std::vector<std::vector<std::vector<std::vector<double>>>>&ystrall);
std::map<int,std::vector<std::vector<double>>> cpmatcal(std::vector<std::vector<double>>slimmatp,
    std::vector<std::vector<double>>slimmatn,std::map<int, std::vector<CoupledBasis>>cpbasisall,
    std::vector<basis>rjln,std::vector<basis>rjlp,
    std::vector<int> tporder,std::vector<int> tvorder,
    std::vector<std::vector<std::vector<double>>> matsliqn,
    std::vector<std::vector<std::vector<double>>> matsliqp,std::vector<double> ktn
    );

#endif //MTCAL_H
