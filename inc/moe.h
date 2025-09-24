//
// Created by wang- on 2024/12/8.
//

#ifndef MOE_H
#define MOE_H
#include <iostream>
#include <vector>
#include "WignerSymbol.hpp"
#include <cmath>
#include <algorithm>
#include "global.h"
#include <unordered_set>
#include <unordered_map>



struct Nucleus
{
    int n; // 主量子数
    int l; // 轨道量子数
    int j; // 总角动量量子数
    int parity; // 宇称

    // 构造函数，根据 l 的奇偶性计算宇称
    Nucleus(int n_val, int l_val, int j_val)
        : n(n_val), l(l_val), j(j_val), parity(l_val % 4 == 0 ? 1 : -1)
    {
    }
};
struct Nucleusm
{
    int n; // 主量子数
    int l; // 轨道量子数
    int j; // 总角动量量子数
    int parity; // 宇称
    int m;
    int num;

    // 构造函数，根据 l 的奇偶性计算宇称
    Nucleusm(int n_val, int l_val, int j_val,int m_val,int num_val)
        : n(n_val), l(l_val), j(j_val), parity(l_val % 4 == 0 ? 1 : -1),m(m_val),num(num_val)
    {
    }
};
std::vector<Nucleus> nucleus{
                            {0, 4, 5},
                            {0, 0, 1},
                            {0, 4, 3}
};
std::vector<Nucleus> nucleus2{
                                {0, 4, 5},
                                {0, 0, 1},
                                {0, 4, 3}
};
std::vector<Nucleusm> nucleusm;
std::vector<Nucleusm> nucleusm2;
int jmax = 5;
std::vector<std::vector<int> > allnucleus = {
    {0, 4, 5},
    {0, 0, 1},
    {0, 4, 3}
};
std::vector<std::vector<int> > allnucleus2 = {
    {0, 4, 5},
    {0, 0, 1},
    {0, 4, 3}
};
std::vector<int> rvecall={0,4};
std::vector<int> rparityvec={1,1};
std::vector<int> rvecall2={0,4};
std::vector<int> rparityvec2={1,1};

extern std::vector<Nucleus> nucleus;
extern std::vector<Nucleus> nucleus2;
extern int jmax;
extern std::vector<std::vector<int> > allnucleus;
extern std::vector<std::vector<int> > allnucleus2;
extern std::vector<int> rvecall;
extern std::vector<int> rparityvec;
extern std::vector<int> rvecall2;
extern std::vector<int> rparityvec2;
extern std::vector<Nucleusm> nucleusm;
extern std::vector<Nucleusm> nucleusm2;
extern std::ofstream outfile;


struct basis {
    int j; // 单粒子
    int jn;
    std::vector<int> r; // 粒子对角动量
    std::vector<int> rn;
    std::vector<int> sj; // 耦合角动量
    std::vector<int> rparity;
    int parity;

    bool operator<(const basis& other) const {
        // 先比较 sj（如果不为空）

        if (parity != other.parity) {
            return parity > other.parity;
        }
        if (sj.empty() != other.sj.empty()) {
            return sj.empty(); // 如果 sj 为空，则排在前面
        }
        if (!sj.empty() && !other.sj.empty() && sj.back() != other.sj.back()) {
            return sj.back() < other.sj.back();
        }
        if (sj != other.sj) {
            return sj < other.sj;
        }


        // 按 jn 排序
        if (jn != other.jn) {
            return jn < other.jn;
        }

        // 按 rn 字典序排序
        if (rn != other.rn) {
            return rn < other.rn;
        }

        // 按 r 字典序排序
        if (r != other.r) {
            return r < other.r;
        }

        // 按 j 排序
        if (j != other.j) {
            return j < other.j;
        }

        // 按 rparity 排序
        if (rparity != other.rparity) {
            return rparity < other.rparity;
        }

        return false; // 如果所有字段都相等，则它们不小于彼此
    }
};

struct basism {
    std::vector<int> r; // 粒子对角动量
    std::vector<int> rn;
    std::vector<int> rparity;
    int parity;
    std::vector<std::vector<int>> mvec;
    std::vector<int> rvecnum;
    std::vector<int> m_rvec;
    std::vector<int> bigmvec;
    bool operator<(const basism& other) const {



        if (parity != other.parity) {
            return parity > other.parity;
        }

        // 按 rn 字典序排序
        if (rn != other.rn) {
            return rn < other.rn;
        }

        // 按 r 字典序排序
        if (r != other.r) {
            return r < other.r;
        }
        // 按 rparity 排序
        if (rparity != other.rparity) {
            return rparity < other.rparity;
        }
        if (mvec != other.mvec) {
            return mvec > other.mvec;
        }
        if (rvecnum != other.rvecnum) {
            return rvecnum > other.rvecnum;
        }

        return false; // 如果所有字段都相等，则它们不小于彼此
    }
};

#include "moeson.cpp"
double ceshi(double x1);

bool isOdd(int x);

void removestr(std::vector<std::vector<std::vector<double> > > &arr, int i);

int sign_func(std::vector<int> x);

double mysqrt(int x);

int deltatwo(int a, int b);

double calculate_sum(const std::vector<std::vector<std::vector<double> > > &v1,
                           const std::vector<std::vector<std::vector<double> > > &v2,
                           int rk, int sn, int rip, int j, int t);

bool is_in_array(const std::vector<int> &arr, int value);


std::vector<int> is_get_in(const std::vector<std::vector<int> > &arr, int num, int jip);



std::vector<int> genmvec(int a, int b);

std::vector<int> genmrivec(std::vector<int> x, int y);



struct getlks
{
    double H;
    std::vector<int> l;
};

struct getlis
{
};

std::vector<std::pair<int, int> > findcpy(const std::vector<Nucleus> &nuclei, int r);



double getphi(const std::vector<std::vector<std::vector<double> > > &A,
                     const std::vector<std::vector<std::vector<double> > > &B,
                     basis rjl, int k);


double ufrom6_j(int a, int b, int d, int c, int e, int f);

double hfrom6_j(struct basis rjl, int k, int s, int L1, int L2);

double gfrom6_j(struct basis rjl, int k, int s, int t, int lk_1);

double mfrom6_j(struct basis rjl, int k, int t, int rkp, int lii);

double qfrom6_j(struct basis rjl, int i, int t, int li, int li_1);

std::vector<int> getsamerange(const std::vector<int> &vec1, const std::vector<int> &vec2);

std::vector<std::vector<double> > getz(std::vector<std::vector<std::vector<double> > > & ystrall1,
                                              std::vector<std::vector<std::vector<double> > > & ystrall2,
                                              int i, int k, int rip, int rk, int ri, basis rjl, int n, int sn, int t);


double getLk(basis rjl, basis rjlinv, int sn, int k, int m, std::vector<int> lall, double h,
             std::vector<std::vector<std::vector<double>>> ystrall,
             std::vector<std::vector<std::vector<double>>> ystrallystrallinv);
#include"moe.cpp"
double overlap(basis x, basis y, std::vector<std::vector<std::vector<double> > > ystrall,
                      std::vector<std::vector<std::vector<double> > > ystrallystrallinv);



double getlis(basis rjl, basis rjlinv, int k, int i, int ii, int t, std::vector<int> lii, int lk_1, int rip,
              std::vector<int> lkall,
              std::vector<std::vector<std::vector<double>>> ystrall,
              std::vector<std::vector<std::vector<double>>> ystrallinv, int contj);
double pair2cal(
    const std::vector<std::vector<std::vector<double>>>& ystra1_1,
    const std::vector<std::vector<std::vector<double>>>& ystra2_1, basis rjl,basis rjlinv);



// 声明 getlis 函数


// 声明 getLk 函数






#endif //MOE_H
