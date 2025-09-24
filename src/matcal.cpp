//
// Created by wang- on 25-1-6.
//
#include "../inc/getbasis.h"
#include <map>
#include "../inc/writeout.h"

#include"../inc/mtcal.h"

#include <fstream>
#include <iostream>
#include <omp.h>
#include <vector>
#include <execution>  // For parallel execution
#include <algorithm>  // For std::for_each with parallel execution
#include <thread>
#include"../inc/moe.h"
std::ofstream outfile("basis_output.txt");

double computeSquareSum(const std::vector<double>& arr) {
    double squareSum = 0.0;
    for (double value : arr) {
        squareSum += value * value;
    }
    return squareSum;
}

void normalizeBySecondColumn(std::vector<std::vector<double>>& matrix) {
    int rows = matrix.size();

    // 找到所有不同的第二列的值
    std::vector<int> secondColValues;
    for (int i = 0; i < rows; ++i) {
        secondColValues.push_back(matrix[i][2]);  // 第二列的索引是2
    }

    // 去重
    std::sort(secondColValues.begin(), secondColValues.end());
    secondColValues.erase(std::unique(secondColValues.begin(), secondColValues.end()), secondColValues.end());

    // 对每个分组（第二列相同的行）进行归一化
    for (int value : secondColValues) {
        // 收集该组的所有行
        std::vector<int> indices;
        for (int i = 0; i < rows; ++i) {
            if (matrix[i][2] == value) {  // 如果第二列的值与当前组值相同
                indices.push_back(i);
            }
        }

        // 对该组的第三列进行归一化
        std::vector<double> groupThirdCol;
        for (int idx : indices) {
            groupThirdCol.push_back(matrix[idx][3]);  // 第三列的索引是3
        }
        double ssum=computeSquareSum(groupThirdCol);
        ssum=sqrt(ssum);
        // 找到该组的第三列的最大值和最小值

        // 进行归一化操作
        for (int idx : indices) {
            matrix[idx][3] = matrix[idx][3]/ssum ;  // 更新第三列的值
        }
    }
}

void removeRange(std::vector<int> &vec, int i, int k) {
    // 将索引转换为从 0 开始数
    int start = i - 1; // 第 i 个元素（从 1 开始数）对应 vec[start]
    int end = k - 1;   // 第 k 个元素（从 1 开始数）对应 vec[end]

    // 检查索引范围是否有效
    if (start >= 0 && end < vec.size() && start <= end) {
        vec.erase(vec.begin() + start, vec.begin() + end+1 ); // 删除索引范围
    } else {
        //std::cout << "error [" << i << ", " << k << "]" << std::endl;
    }
}
std::vector<std::vector<double> > ychangeb(
    const std::vector<std::vector<std::vector<double> > > &ystrall1,
    std::vector<std::vector<std::vector<double> > > &qstrall,
    int k, int rkp, int tnum, std::vector<int> torder, basis rjl)
{
    int nucnum = nucleus.size();
    int N2 = rjl.r.size();
    int rk = rjl.r[k - 1];
    int t = torder[tnum];
    double c4 =  mysqrt(t + 1) * mysqrt(rkp + 1);
    std::vector<std::vector<double> > ychange;
    ychange.resize(nucnum);
    for (int i_m = 0; i_m < nucnum; i_m++)
    {
        ychange[i_m].resize(nucnum, 0.0); // 初始化二维数组的每一行
    }
    for (int a = 0; a < nucnum; a++)
    {
        for (int b = 0; b < nucnum; b++)
        {
            double sum = 0;
            int l_a = nucleus[a].l;
            int l_b = nucleus[b].l;
            int j_a = nucleus[a].j;
            int j_b = nucleus[b].j;


                for (int d = 0; d < nucnum; d++)
                {
                    for (int d0 = 0; d0 < nucnum; d0++)
                    {
                        int jd = nucleus[d].j;
                        int jd0 = nucleus[d0].j;
                        double y1 = ystrall1[k - 1][d][d0] * qstrall[tnum][d0][b] * qstrall[tnum][d][a];
                        y1 = y1 * util::wigner_6j(rk, t, rkp, j_b, jd, jd0) * util::wigner_6j(rk, t, rkp, jd, j_b, j_a);
                        double y2 = ystrall1[k - 1][d][d0] * qstrall[tnum][d0][a] * qstrall[tnum][d][b];
                        y2 = y2 * util::wigner_6j(rk, t, rkp, j_a, jd, jd0) * util::wigner_6j(rk, t, rkp, jd, j_a, j_b);
                        y2 = (std::pow(-1, (j_a + j_b + rk) / 2.0)) * y2;
                        sum = sum + y1 - y2;
                    }
                }

            ychange[a][b] = c4 * sum;
        }
    }
    return ychange;
}

std::vector<std::vector<double> > ychangep(
    std::vector<std::vector<std::vector<double> > > ystrall1,
    std::vector<std::vector<std::vector<double> > > qstrall,
    int k, int rkp, int tnum, std::vector<int> torder, basis rjl)
{
    int nucnum = nucleus.size();
    int N2 = rjl.r.size();
    int rk = rjl.r[k - 1];
    int t = torder[tnum];
    double c4 = mysqrt(rk + 1) * mysqrt(t + 1) ;
    std::vector<std::vector<double> > ychange;
    ychange.resize(nucnum);
    for (int i_m = 0; i_m < nucnum; i_m++)
    {
        ychange[i_m].resize(nucnum, 0.0); // 初始化二维数组的每一行
    }
    for (int d = 0; d < nucnum; d++)
    {
        for (int a = 0; a < nucnum; a++)
        {
            double sum = 0;
            int l_d = nucleus[d].l;
            int l_a = nucleus[a].l;
            int j_d = nucleus[d].j;
            int j_a = nucleus[a].j;
            if ( (j_d + j_a) >= rkp
                && std::abs(j_d - j_a) <= rkp)
            {
                for (int bk = 0; bk < nucnum; bk++)
                {

                        int jb = nucleus[bk].j;

                        double y1 = ystrall1[k - 1][a][bk] * qstrall[tnum][bk][d] ;

                        y1 = y1 * util::wigner_6j(rk, t, rkp, j_d, j_a, jb) ;
                        double y2 = ystrall1[k - 1][d][bk] * qstrall[tnum][bk][a] ;

                        y2 = y2 * util::wigner_6j(rk, t, rkp, j_a, j_d, jb) ;
                        y2 = (std::pow(-1, ((j_a + j_d + rkp) / 2.0)) )* y2;
                        sum = sum + y1 - y2;

                }
            }
            ychange[d][a] = c4 * sum;
        }
    }
    return ychange;
}

std::vector<std::vector<double> > ystrchaa
(std::vector<std::vector<std::vector<double> > > & ystrall1,
std::vector<std::vector<std::vector<double> > > & ystrall2,
int i, int k, int rip,  basis rjlinv, int snum, std::vector<int> sorder,int t)
{
    int s=sorder[snum];
    int ri;
    if (i==0)
    {
        ri=rjlinv.j;
    }
    else
    {
        ri=rjlinv.r[i-1];
    }
    int rk=rjlinv.r[k-1];

    int nucnum = nucleus.size();
    int N2 = rjlinv.r.size();
    double c4 = 4 * mysqrt(rk + 1) * mysqrt(ri + 1) * mysqrt(s + 1) * mysqrt(t + 1);
    std::vector<std::vector<double> > ychange;
    ychange.resize(nucnum);
    for (int i_m = 0; i_m < nucnum; i_m++) {
        ychange[i_m].resize(nucnum, 0.0);  // 初始化二维数组的每一行
    }
    for (int i_m = 0; i_m < nucnum; i_m++)
    {
        for (int i_n = 0; i_n < nucnum; i_n++)
        {
            double sum = 0;
            int l_n = nucleus[i_n].l;
            int l_m = nucleus[i_m].l;
            int j_n = nucleus[i_n].j;
            int j_m = nucleus[i_m].j;
            if ( (j_n + j_m) >= rip && std::abs(j_n - j_m) <= rip)
            {
                for (int b_m = 0; b_m < nucnum; b_m++)
                {
                    for (int b_n = 0; b_n < nucnum; b_n++)
                    {
                        int jb = nucleus[b_m].j;
                        int jbp = nucleus[b_n].j;
                        double y1 = ystrall2[snum][b_m][b_n] * ystrall1[k - 1][i_m][b_m] * ystrall1[i -
                        1][i_n][b_n];
                        y1 = y1 * util::wigner_6j(rk, s, t, jbp, j_m, jb) * util::wigner_6j(ri, t, rip, j_m, j_n, jbp);
                        double y2=ystrall2[snum ][b_m][b_n] * ystrall1[k - 1][i_n][b_m] * ystrall1[i - 1][i_m][b_n];
                        y2=y2*util::wigner_6j(rk, s, t, jbp, j_n, jb) * util::wigner_6j(ri, t, rip, j_n, j_m, jbp);
                        y2=(std::pow(-1,(j_n+j_m+rip)/2.0))*y2;
                        sum = sum + y1-y2;
                    }
                }

            }
            ychange[i_m][i_n] = c4*sum;
        }
    }
    return ychange;
}



double h0cal(basis rjl, basis rjlinv,
             std::vector<std::vector<std::vector<double> > > ystrall,
             std::vector<std::vector<std::vector<double> > > ystrallinv,
             std::vector<double> energy)
{
    // if (!outfile.is_open())
    // {
    //     // 如果文件未打开，则尝试打开文件
    //     outfile.open("basis_output.txt",std::ios::app);
    // }
    // outfile<<"____________________h0cal_____________________"<<std::endl;
    // outfile<<"rjl1"<<std::endl;
    // printBasisOneLinefile(outfile,rjl);
    // outfile<<"rjl2"<<std::endl;
    // printBasisOneLinefile(outfile,rjlinv);
    // if (!outfile.is_open())
    // {
    //     // 如果文件未打开，则尝试打开文件
    //     outfile.open("basis_output.txt",std::ios::app);
    // }
    // outfile<<"ystrall1_1"<<std::endl;
    // writeMatrixToFile(ystrall[0],outfile);
    // if (!outfile.is_open())
    // {
    //     // 如果文件未打开，则尝试打开文件
    //     outfile.open("basis_output.txt",std::ios::app);
    // }
    // outfile<<"ystrall1_2"<<std::endl;
    // writeMatrixToFile(ystrall[1],outfile);
    // if (!outfile.is_open())
    // {
    //     // 如果文件未打开，则尝试打开文件
    //     outfile.open("basis_output.txt",std::ios::app);
    // }
    // outfile<<"ystrall1_3"<<std::endl;
    // writeMatrixToFile(ystrall[2],outfile);
    // if (!outfile.is_open())
    // {
    //     // 如果文件未打开，则尝试打开文件
    //     outfile.open("basis_output.txt",std::ios::app);
    // }
    // outfile<<"ystrall2_1"<<std::endl;
    // writeMatrixToFile(ystrallinv[0],outfile);
    // if (!outfile.is_open())
    // {
    //     // 如果文件未打开，则尝试打开文件
    //     outfile.open("basis_output.txt",std::ios::app);
    // }
    // outfile<<"ystrall2_2"<<std::endl;
    // writeMatrixToFile(ystrallinv[1],outfile);
    // if (!outfile.is_open())
    // {
    //     // 如果文件未打开，则尝试打开文件
    //     outfile.open("basis_output.txt",std::ios::app);
    // }
    // outfile<<"ystrall2_3"<<std::endl;
    // writeMatrixToFile(ystrallinv[2],outfile);
    int N = rjl.r.size();
    int nucnum = nucleus.size();
    double sum = 0.0;
    for (int k = 0; k <= N; k++)
    {
        std::vector<std::vector<std::vector<double> > > ystrall1 = ystrallinv;
        if (k == 0)
        {
            if (rjl.j==0)
            {
                continue;
            }
            double e=energy[rjlinv.jn];
            sum = sum + (e*overlap(rjl, rjlinv, ystrall, ystrallinv));
        } else
        {
            std::vector<std::vector<double> > y = ystrallinv[k - 1];
            for (int i_m = 0; i_m < nucnum; i_m++)
            {
                for (int j_m = 0; j_m < nucnum; j_m++)
                {
                    double change=0;

                    change=y[i_m][j_m]* (energy[i_m]+ energy[j_m]);

                    ystrall1[k - 1][i_m][j_m] = change ;

                }
            }
            sum = sum + overlap(rjl, rjlinv, ystrall, ystrall1);




            // if (!outfile.is_open())
            // {
            //     // 如果文件未打开，则尝试打开文件
            //     outfile.open("basis_output.txt",std::ios::app);
            // }
            // outfile<<"_______________k______________="<<k<<std::endl;
            // outfile<<"rjl1"<<std::endl;
            // printBasisOneLinefile(outfile,rjl);
            // outfile<<"rjl2"<<std::endl;
            // printBasisOneLinefile(outfile,rjlinv);
            // if (!outfile.is_open())
            // {
            //     // 如果文件未打开，则尝试打开文件
            //     outfile.open("basis_output.txt",std::ios::app);
            // }
            // outfile<<"ystrall1_1"<<std::endl;
            // writeMatrixToFile(ystrall[0],outfile);
            // if (!outfile.is_open())
            // {
            //     // 如果文件未打开，则尝试打开文件
            //     outfile.open("basis_output.txt",std::ios::app);
            // }
            // outfile<<"ystrall1_2"<<std::endl;
            // writeMatrixToFile(ystrall[1],outfile);
            // if (!outfile.is_open())
            // {
            //     // 如果文件未打开，则尝试打开文件
            //     outfile.open("basis_output.txt",std::ios::app);
            // }
            // outfile<<"ystrall1_3"<<std::endl;
            // writeMatrixToFile(ystrall[2],outfile);
            // if (!outfile.is_open())
            // {
            //     // 如果文件未打开，则尝试打开文件
            //     outfile.open("basis_output.txt",std::ios::app);
            // }
            // outfile<<"ystrall2_1"<<std::endl;
            // writeMatrixToFile(ystrall1[0],outfile);
            // if (!outfile.is_open())
            // {
            //     // 如果文件未打开，则尝试打开文件
            //     outfile.open("basis_output.txt",std::ios::app);
            // }
            // outfile<<"ystrall2_2"<<std::endl;
            // writeMatrixToFile(ystrall1[1],outfile);
            // if (!outfile.is_open())
            // {
            //     // 如果文件未打开，则尝试打开文件
            //     outfile.open("basis_output.txt",std::ios::app);
            // }
            // outfile<<"ystrall2_3"<<std::endl;
            // writeMatrixToFile(ystrall1[2],outfile);
            //
            // if (!outfile.is_open())
            // {
            //     // 如果文件未打开，则尝试打开文件
            //     outfile.open("basis_output.txt",std::ios::app);
            // }
            // outfile<<"sum1"<< "hocalk= " << k<<" hocalsum= " <<  sum<<std::endl;

        }
        // std::cout  << "hocalk= " << k<<" hocalsum= " <<  sum << std::endl;
    }
    return sum;
}

double qt2term(basis rjl,basis rjlinv,
std::vector<std::vector<std::vector<double> > > ystrall,
std::vector<std::vector<std::vector<double> > > ystrallinv,
int i,int ii,int k, std::vector<int>rkp,std::vector<int>lall,int rip,
std::vector<std::vector<std::vector<double>> > qstrall,std::vector<int>torder,int tnum)
{
    int t=torder[tnum];
    int rk=rjlinv.r[k-1];
    if (ii==k-1)
    {
        int lk_1=0;
        double sum=0;
        for (int rkpp:rkp)
        {
            if (!lall.empty())
            {
                 lk_1 = lall.back();
            }else
            {
                return 0;
            }
            std::vector<int>lk_1line=genmvec(rjlinv.sj[k-1],rkpp);
            if (is_in_array(lk_1line,lk_1))
            {
                std::vector<std::vector<std::vector<double> > >ystrallinvkp=ystrallinv;
                std::vector<std::vector<double>>ystrallchp;
                ystrallchp=ychangep(ystrallinv, qstrall, k, rkpp, tnum, torder, rjlinv);
                ystrallinvkp[k-1]=ystrallchp;
                int liinm=lall[0];
                double m ;
                if (i==0)
                {
                    m=1;
                }else
                {
                    m=mfrom6_j(rjlinv, i, t, rip, liinm);
                }
                std::vector<int>lk1;
                std::vector<int>lk2;
                lk1=genmvec(lk_1,rk);
                lk2=genmvec(t,rjlinv.sj[k-1]);
                std::vector<int>lk=getsamerange(lk1,lk2);
                double mr=0;
                for (int lkk:lk)
                {
                    mr=mr+ufrom6_j(rk,t,lk_1,rjlinv.sj[k-1],rkpp,lkk);
                }
                double q = 1;
                if (i == k - 1)
                {
                    q = 1;
                } else if (i == k - 2)
                {
                    q = qfrom6_j(rjlinv, k - 1, t, lk_1, lall[lall.size()-2]);
                } else
                {
                    int i_cc = 0;
                    for (int i_c = i + 1; i_c <= k - 1; i_c++)
                    {
                        int li_2 = lall[i_cc];
                        int li_1 = lall[i_cc + 1];
                        i_cc++;
                        q = q * qfrom6_j(rjlinv, i_c, t, li_1, li_2);
                    }
                }
                int jk_1;
                if (k==1)
                {
                    jk_1=rjlinv.j;
                }else
                {
                    jk_1=rjlinv.sj[k-2];
                }
                int jk=rjlinv.sj[k-1];
                double g=2*mysqrt(rkpp+1)*ufrom6_j(jk_1,t,jk,rkpp,lk_1,rk)/mysqrt(rk+1);
                std::vector<int>liall=lall;
                basis rjlinvcha=rjlinv;
                if (i==0)
                {
                    rjlinvcha.j=rip;
                    if (!liall.empty())
                    {
                        // 检查是否为空
                        liall.erase(liall.begin()); // 删除第一个元素
                    }else
                    {
                        //std::cout<<"error"<<std::endl;
                    }
                    removeRange(rjlinvcha.sj,1,k-1);
                    rjlinvcha.r[k-1]=rkpp;
                    rjlinvcha.sj.insert(rjlinvcha.sj.begin()+0,liall.begin(),liall.end());

                }else
                    {
                    rjlinvcha.r[i-1]=rip;
                    removeRange(rjlinvcha.sj,i,k-1);
                    rjlinvcha.sj.insert(rjlinvcha.sj.begin()+i-1,liall.begin(),liall.end());
                    }
                rjlinvcha.r[k-1]=rkpp;
                mr=1;
                sum=sum+(q*g*m*mr*overlap(rjl,rjlinvcha,ystrall,ystrallinvkp));
            }
        }
        return sum;
    }else
    {
        std::vector<int>li1all1=genmvec(lall.back(),rjlinv.r[ii]);
        std::vector<int>li1all2=genmvec(rjlinv.sj[ii],t);
        std::vector<int>li1all=getsamerange(li1all1,li1all2);
        std::vector<int>lallchange=lall;
        double sum=0;
        for (int li1:li1all)
        {
            lallchange.push_back(li1);
            sum=sum+qt2term( rjl, rjlinv,ystrall,ystrallinv,
                i, ii+1, k, rkp,lallchange, rip,qstrall,torder, tnum);
            lallchange.pop_back();
        }
        return sum;
    }

}

double qtqt(basis rjl, basis rjlinvcha,
            std::vector<std::vector<std::vector<double> > > ystrall,
            std::vector<std::vector<std::vector<double> > > ystrallinv,
            std::vector<std::vector<std::vector<double> > > qstrall, int tnum,
            std::vector<int> qorder)
{
    if (rjlinvcha.sj.back()!=rjl.sj.back())
    {
        return 0;
    }
    basis rjlinv=rjlinvcha;
    int N = rjl.r.size();
    int jinv = rjlinv.j;
    int jinvn = rjlinv.jn;
    int t = qorder[tnum];
    double sum = 0.0;
    for (int k = 0; k <= N; k++)
    {
        std::vector<int> rkp;
        if (k == 0)
        {
            if (rjl.j==0)
            {
                continue;
            }
            rkp = genmvec(t, jinv);
            for (int rkpp: rkp)
            {
                int contj = -1;
                std::vector<int> rkp2 = is_get_in(allnucleus, 2, rkpp);
                int par = nucleus[rjlinv.jn].parity;
                if (isOdd(t / 2)) { par = par * (-1); } else { par = par * (1); }
                if (rkp2.empty())
                {
                    continue;
                } else
                {
                    for (const auto &elementnum: rkp2)
                    {
                        if (nucleus[elementnum].parity == par)
                        {
                            contj = elementnum;
                        }
                    }
                    if (contj == -1)
                    {
                        continue;
                    }
                }
                double c = 0;
                double c1 = 0;
                std::vector<int> sig;
                sig = {rkpp, -jinv, -t};
                c = sign_func(sig);
                c = c * mysqrt(rkpp + 1) / mysqrt(jinv + 1);
                c1 = (t + 1) / (mysqrt(jinv + 1) * mysqrt(rkpp + 1));
                c1 = c1 * qstrall[tnum][contj][jinvn] * qstrall[tnum][jinvn][contj];
                sum = sum + c1 * c * overlap(rjl, rjlinv, ystrall, ystrallinv);
            }
        } else
        {
            int rk = rjlinv.r[k - 1];
            rkp = genmvec(t, rjlinv.r[k - 1]);
            for (int rkpp: rkp)
            {
                std::vector<std::vector<double> > ystrchange;
                ystrchange = ychangeb(ystrallinv, qstrall, k, rkpp, tnum, qorder, rjlinv);
                std::vector<std::vector<std::vector<double> > > ystrallch = ystrallinv;
                ystrallch[k - 1] = ystrchange;
                double c = 0;
                std::vector<int> sig;
                sig = {rkpp, -rk, -t};
                c = sign_func(sig);
                c = c * mysqrt(rkpp + 1) / mysqrt(rk + 1);
                sum = sum + c * overlap(rjl, rjlinv, ystrall, ystrallch);
            }
        }
    }

    for (int k = 1; k <= N; k++)
    {
        for (int i = 0; i <= k-1; i++)
        {   if (i==0)
        {
            if (rjlinv.j==0)
            {
                continue;
            }
        }
            basis rjlinvp=rjlinvcha;
            std::vector<int> r=rjlinvcha.r;
            std::vector<int> rip;
            std::vector<int> rkp;
            double c1=1;
            double c2=1;
            rkp=genmvec(r[k-1], t);
            std::vector<std::vector<std::vector<double> > > ystrallinvchan=ystrallinv;
            if (i==0)
            {
                rip=genmvec(rjlinvcha.j,t);

                for (int ripp: rip)
                {
                    int contj = -1;
                    std::vector<int> rip2= is_get_in(allnucleus, 2, ripp);
                    int par = nucleus[rjlinvcha.jn].parity;
                    if (isOdd(t / 2))
                    {
                        par = par * (-1);
                    } else
                    {
                        par = par * (1);
                    }
                    if (rip2.empty())
                    {
                        continue;
                    } else
                    {
                        for (const auto &elementnum: rip2)
                        {
                            if (nucleus[elementnum].parity == par)
                            {
                                contj = elementnum;
                            }
                        }
                        if (contj == -1)
                        {
                            continue;
                        }
                    }

                    rjlinv.jn=contj;
                    c1=-qstrall[tnum][contj][rjlinv.jn];
                    c2=mysqrt(t + 1) / mysqrt(ripp + 1);
                    sum=sum+c1*c2*qt2term(rjl,rjlinv,ystrall,ystrallinv,
                            i,i,k,rkp,{ripp},ripp,qstrall,qorder,tnum);
                }
            }else
            {
                rip=genmvec(r[i-1],t);
                for (int ripp: rip)
                {
                    std::vector<std::vector<double>> ystrchangpp=
                        ychangep(ystrallinv, qstrall, i, ripp, tnum, qorder, rjlinvcha);
                    ystrallinvchan[i-1]=ystrchangpp;
                    std::vector<int>lii1;
                    lii1=genmvec(rjlinv.sj[i-1],t);
                    std::vector<int>lii2;
                    lii2=lii1;
                    if (i!=1)
                    {
                        lii2=genmvec(rjlinv.sj[i-2],ripp);
                    }
                    else
                    {
                        lii2=genmvec(rjlinv.j,ripp);
                    }
                    std::vector<int>lii=getsamerange(lii1, lii2);
                    for (int lii1_1: lii)
                    {
                        sum=sum+qt2term(rjl,rjlinv,ystrall,ystrallinvchan,
                            i,i,k,rkp,{lii1_1},ripp,qstrall,qorder,tnum);
                    }
                }
            }
        }
    }
    return sum;
}


double qtqteventerm2(basis rjl,basis rjlinv,
std::vector<std::vector<std::vector<double> > > ystrall,
std::vector<std::vector<std::vector<double> > > ystrallinv,
int i,int ii,std::vector<int>lall,int ripp,
std::vector<std::vector<std::vector<double>> > qstrall,std::vector<int>torder,int tnum,int l1)
{
    int n=rjl.r.size();
    int t=torder[tnum];
    if (ii==rjlinv.r.size())
    {
        double sum=0;
        int ln_1;
        int lir;
        if (i==n)
        {
            ln_1=rjlinv.sj[i-2];
            lir=ripp;
        }else
        {
            ln_1=lall.back();
            lir=rjlinv.r[n-1];
        }
        std::vector<int>lkm1=genmvec(ln_1,lir);
        std::vector<int>lkm2=genmvec(rjlinv.sj[n-1],t);
        std::vector<int>lkm=getsamerange(lkm1,lkm2);
        if (is_in_array(lkm,l1))
        {
            lall.push_back(l1);
            std::vector<int> sjcha=rjlinv.sj;
            removeRange(sjcha,i,n);
            sjcha.insert(sjcha.end(),lall.begin(),lall.end());
            double q = 1;
            if (i == n)
            {
                q = 1;
            } else if (i == n - 1)
            {
                q = qfrom6_j(rjlinv, n, t, lall.back(), lall[lall.size()-2]);
            } else
            {
                int i_cc = 0;
                for (int i_c = i + 1; i_c <= n; i_c++)
                {
                    int li_2 = sjcha[i_c-2];
                    int li_1 = sjcha[i_c -1];
                    i_cc++;
                    q = q * qfrom6_j(rjlinv, i_c, t, li_1, li_2);
                }
            }
            double m=mfrom6_j(rjlinv,i,t,ripp,lall[0]);
            ystrallinv[i-1]=ychangep(ystrallinv,qstrall,i,ripp,tnum,torder,rjlinv);
            basis rjl2=rjlinv;
            rjl2.r[i-1]=ripp;
            removeRange(rjl2.sj,i,n);
            rjl2.sj.insert(rjl2.sj.end(),lall.begin(),lall.end());


            sum=q*m*overlap(rjl,rjl2,ystrall,ystrallinv);
            if (!outfile.is_open())
            {
                // 如果文件未打开，则尝试打开文件
                outfile.open("basis_output.txt",std::ios::app);
            }
            outfile<<"i="<<i<<std::endl;

            outfile<< "q2="<<q<<std::endl;
            outfile<< "m2="<<m<<std::endl;
            outfile<<"basis1"<<std::endl;
            printBasisOneLinefile(outfile,rjl);
            outfile<<"basis2"<<std::endl;
            printBasisOneLinefile(outfile,rjl2);
            outfile<< "product="<<sum<<std::endl;
            outfile<<"ystrall"<<std::endl;
            writeMatrixToFile(ystrall[0],outfile);
            if (!outfile.is_open())
            {
                // 如果文件未打开，则尝试打开文件
                outfile.open("basis_output.txt",std::ios::app);
            }
            outfile<<"ystrall1_2"<<std::endl;
            writeMatrixToFile(ystrall[1],outfile);
            if (!outfile.is_open())
            {
                // 如果文件未打开，则尝试打开文件
                outfile.open("basis_output.txt",std::ios::app);
            }
            outfile<<"ystrall2_1"<<std::endl;
            if (!outfile.is_open())
            {
                // 如果文件未打开，则尝试打开文件
                outfile.open("basis_output.txt",std::ios::app);
            }
            writeMatrixToFile(ystrallinv[0],outfile);
            if (!outfile.is_open())
            {
                // 如果文件未打开，则尝试打开文件
                outfile.open("basis_output.txt",std::ios::app);
            }
            outfile<<"ystrall2_2"<<std::endl;
            writeMatrixToFile(ystrallinv[1],outfile);
            outfile.close();
            return sum;

        }else{return 0;}

    }
    else
    {
        double sum=0.0;
        std::vector<int> li1;
        std::vector<int> li2;
        std::vector<int> li;
        int ln_1;
        if (ii==i)
        {
            if (i==1)
            {
                li1={ripp};
            }else
            {
                ln_1=rjlinv.sj[i-2];
                li1=genmvec(ln_1,ripp);
            }

        }else
        {
            ln_1=lall.back();
            li1=genmvec(ln_1,rjlinv.r[ii-1]);
        }
        li2=genmvec(rjlinv.sj[ii-1],t);

        li=getsamerange(li1,li2);
        for (int li1_1: li)
        {
            lall.push_back(li1_1);

            sum=sum+qtqteventerm2(rjl,rjlinv,ystrall,ystrallinv,i,ii+1,lall,ripp,qstrall,torder,tnum,l1);
            lall.pop_back();
            double m;
        }

        return sum;
    }

}




double qtqteventerm1(basis rjl,basis rjlinv,
std::vector<std::vector<std::vector<double> > > ystrall,
std::vector<std::vector<std::vector<double> > > ystrallinv,
int k,int kk,std::vector<int>lall,int rkpp,
std::vector<std::vector<std::vector<double>> > qstrall,std::vector<int>torder,int tnum,int l1)
{
   int n=rjl.r.size();
    int t=torder[tnum];
    if (kk==rjl.r.size())
    {
        double sum=0;
        int ln_1;
        int lkr;
        if (k==n)
        {
            ln_1=rjl.sj[k-2];
             lkr=rkpp;
        }else
        {
            ln_1=lall.back();
             lkr=rjl.r[n-1];
        }
        std::vector<int>lkm1=genmvec(ln_1,lkr);
        std::vector<int>lkm2=genmvec(rjl.sj[n-1],t);
        std::vector<int>lkm=getsamerange(lkm1,lkm2);
        if (is_in_array(lkm,l1))
        {
            lall.push_back(l1);
            std::vector<int> sjcha=rjl.sj;
            removeRange(sjcha,k,n);
            sjcha.insert(sjcha.end(),lall.begin(),lall.end());
            double q = 1;
            if (k == n)
            {
                q = 1;
            } else if (k == n - 1)
            {
                q = qfrom6_j(rjl, n , t, lall.back(), lall[lall.size()-2]);
            } else
            {
                int i_cc = 0;
                for (int i_c = k + 1; i_c <= n; i_c++)
                {
                    int li_2 = sjcha[i_c-2];
                    int li_1 = sjcha[i_c -1];
                    i_cc++;
                    q = q * qfrom6_j(rjl, i_c, t, li_1, li_2);
                }
            }
            double m=mfrom6_j(rjl,k,t,rkpp,lall[0]);
            if (!outfile.is_open())
            {
                // 如果文件未打开，则尝试打开文件
                outfile.open("basis_output.txt",std::ios::app);
            }
            outfile<<"k="<<k<<std::endl;
            outfile<< "q1="<<q<<std::endl;
            outfile<< "m1="<<m<<std::endl;
            ystrall[k-1]=ychangep(ystrall,qstrall,k,rkpp,tnum,torder,rjl);
            basis rjl1=rjl;
            rjl1.r[k-1]=rkpp;
            removeRange(rjl1.sj,k,n);

            rjl1.sj.insert(rjl1.sj.end(),lall.begin(),lall.end());
            double ovrjl1;
            ovrjl1=overlap(rjl1,rjl1,ystrall,ystrall);
            // if (abs(ovrjl1)<=1e- )
            // {
            //     return 0;
            // }

            for (int i=1; i<=n; i++)
            {
                std::vector<int> rip;
                rip=genmvec(rjlinv.r[i-1],t);
                for (int ripp: rip)
                {
                    std::vector<int>lall1={};
                    sum=sum+(m*q*qtqteventerm2(rjl1,rjlinv,ystrall,ystrallinv,i,i,lall1,ripp,qstrall,torder,tnum,l1));
                }
            }
            return sum;

        }else{return 0;}

    }
    else
    {
        double sum=0;
        std::vector<int> li1;
        std::vector<int> li2;
        std::vector<int> li;
        int ln_1;
        if (kk==k)
        {
            if (k==1)
            {
                li1={rkpp};
            }else
            {
                ln_1=rjl.sj[k-2];
                li1=genmvec(ln_1,rkpp);
            }

        }else
        {
            ln_1=lall.back();
            li1=genmvec(ln_1,rjl.r[kk-1]);
        }
        li2=genmvec(rjl.sj[kk-1],t);
        li=getsamerange(li1,li2);
        for (int li1_1: li)
        {
            lall.push_back(li1_1);
            sum=sum+qtqteventerm1(rjl,rjlinv,ystrall,ystrallinv,k,kk+1,lall,rkpp,qstrall,torder,tnum,l1);
            lall.pop_back();
        }

    return sum;
    }
}

double qtqteven(basis rjl, basis rjlinv,
            std::vector<std::vector<std::vector<double> > > ystrall,
            std::vector<std::vector<std::vector<double> > > ystrallinv,
            std::vector<std::vector<std::vector<double> > > qstrall, int tnum,
            std::vector<int> qorder)
{
    if (!outfile.is_open())
    {
        // 如果文件未打开，则尝试打开文件
        outfile.open("basis_output.txt",std::ios::app);
    }
    outfile<<"basisrjl"<<std::endl;
    printBasisOneLinefile(outfile,rjl);
    outfile<<"basisrjlinv"<<std::endl;
    printBasisOneLinefile(outfile,rjlinv);

    double sum = 0;
    if (rjl.sj.back()!=rjlinv.sj.back())
    {
        return 0;
    }
    int t=qorder[tnum];
    std::vector<int> L1;
    L1=genmvec(rjl.sj.back(),t);
    for (int l1: L1)
    {
        outfile<<"LN="<<l1<<std::endl;
        double c1=(l1+1.0)/(rjl.sj.back()+1.0);
        for (int k=1;k<=rjl.r.size();k++)
        {
            std::vector<int> rkp;
            rkp=genmvec(rjl.r[k-1],t);
            for (int rkpp: rkp)
            {

                sum=sum+c1*qtqteventerm1(rjl,rjlinv,ystrall,ystrallinv,k,k,{},rkpp,qstrall,qorder,tnum,l1);
            }
        }

    }
    if (!outfile.is_open())
    {
        // 如果文件未打开，则尝试打开文件
        outfile.open("basis_output.txt",std::ios::app);
    }
    outfile<<"sum="<<sum<<std::endl;
    outfile<<"_________________________________________________________________________________________________________"<<std
    ::endl;
    return sum;

}


double getzinq (const std::vector<std::vector<std::vector<double> > > &ystrall1,
    std::vector<std::vector<std::vector<double> > > &qstrall,
    int k, int rkp, int tnum, std::vector<int> torder, basis rjl,int a,int b,int d)
{
    int nucnum = nucleus.size();
    double t=torder[tnum];
    int rk=rjl.r[k-1];
    double c1=mysqrt(rk+1)*mysqrt(t+1);
    double sum=0.0;

    for(int d0=0;d0<nucnum;d0++)
    {
        int j_a = nucleus[a].j;
        int j_b = nucleus[b].j;
        int j_d = nucleus[d].j;
        int j_d0=nucleus[d0].j;
        double re=ystrall1[k-1][d][d0];
        re=re*qstrall[tnum][d0][b];
        re=re*util::wigner_6j(rk,t,rkp,j_b,j_d,j_d0);
        sum=sum+re;
    }
    sum=sum*c1;
    return sum;

}

double getyinq (const std::vector<std::vector<std::vector<double> > > &ystrall1,
    std::vector<std::vector<std::vector<double> > > &qstrall,
    int k, int rkp, int tnum, std::vector<int> torder, basis rjl,int a,int b,int d)
{
    int nucnum = nucleus.size();
    double t=torder[tnum];
    double y1=getzinq(ystrall1,qstrall,k,rkp,tnum,torder,rjl,a,b,d);
    double y2=getzinq(ystrall1,qstrall,k,rkp,tnum,torder,rjl,a,d,b);
    int j_b = nucleus[b].j;
    int j_d = nucleus[d].j;
    std::vector<int> signorder={j_b,j_d,rkp};
    int sign=sign_func(signorder);
    double re=y1-(sign*y2);
    return re;

}



double getzinq (const std::vector<std::vector<std::vector<double> > > &ystrall1,
    std::vector<std::vector<std::vector<double> > > &qstrall,
    int k, int rkp, int tnum, std::vector<int> torder, basis rjl,int a,int b)
{
    int nucnum = nucleus.size();
    double t=torder[tnum];
    double c1=mysqrt(rkp+1)*mysqrt(t+1);
    double sum=0.0;
    int rk=rjl.r[k-1];

    for (int d=0; d<nucnum; d++)
    {
        int j_a = nucleus[a].j;
        int j_b = nucleus[b].j;
        int j_d = nucleus[d].j;
        double y1=getyinq(ystrall1,qstrall,k,rkp,tnum,torder,rjl,a,b,d);
        double re=y1*qstrall[tnum][d][a];
        re=re*util::wigner_6j(rkp,t,rk,j_a,j_b,j_d);
        sum=sum+re;
    }
    sum=sum*c1;
    return sum;
}


std::vector<std::vector<double> > ychangeb1(
    const std::vector<std::vector<std::vector<double> > > &ystrall1,
    std::vector<std::vector<std::vector<double> > > &qstrall,
    int k, int rkp, int tnum, std::vector<int> torder, basis rjl)
{
    int nucnum = nucleus.size();
    int N2 = rjl.r.size();
    int rk = rjl.r[k - 1];
    int t = torder[tnum];
    double c4 =  mysqrt(t + 1) * mysqrt(rkp + 1);
    std::vector<std::vector<double> > ychange;
    ychange.resize(nucnum);
    for (int i_m = 0; i_m < nucnum; i_m++)
    {
        ychange[i_m].resize(nucnum, 0.0); // 初始化二维数组的每一行
    }
    for (int a = 0; a < nucnum; a++)
    {
        for (int b = 0; b < nucnum; b++)
        {
            int j_a = nucleus[a].j;
            int j_b = nucleus[b].j;
            std::vector<int> signorder={j_a,j_b,rk};
            int sign=sign_func(signorder);
            double z1=getzinq(ystrall1,qstrall,k,rkp,tnum,torder,rjl,a,b);
            double z2=getzinq(ystrall1,qstrall,k,rkp,tnum,torder,rjl,b,a);
            double re=z1-(sign*z2);
            ychange[a][b]=re;
        }
    }
    return ychange;
}

std::vector<std::vector<double> > ychangep1(
    std::vector<std::vector<std::vector<double> > > ystrall1,
    std::vector<std::vector<std::vector<double> > > qstrall,
    int k, int rkp, int tnum, std::vector<int> torder, basis rjl)
{
    int nucnum = nucleus.size();
    int N2 = rjl.r.size();
    int rk = rjl.r[k - 1];
    int t = torder[tnum];
    double c4 = mysqrt(rk + 1) * mysqrt(t + 1) ;
    std::vector<std::vector<double>> ychange(nucnum, std::vector<double>(nucnum, 0.0));

    for (int a = 0; a < nucnum; a++)
    {
        for (int b=0;b<nucnum;b++)
        {
            double re=getyinq(ystrall1,qstrall,k,rkp,tnum,torder,rjl,1,a,b);
            ychange[a][b]=re;
        }
    }
    return ychange;
}



double qt2term1(basis rjl,basis rjlinv,
std::vector<std::vector<std::vector<double> > > ystrall,
std::vector<std::vector<std::vector<double> > > ystrallinv,
int i,int ii,int k, std::vector<int>rkp,std::vector<int>lall,int rip,
std::vector<std::vector<std::vector<double>> > qstrall,std::vector<int>torder,int tnum)
{
    int t=torder[tnum];
    int rk=rjlinv.r[k-1];
    if (ii==k-1)
    {
        int lk_1=0;
        double sum=0;
        for (int rkpp:rkp)
        {
            lk_1 = lall.back();
            std::vector<int>lk_1line=genmvec(rjlinv.sj[k-1],rkpp);
            if (is_in_array(lk_1line,lk_1))
            {
                std::vector<std::vector<std::vector<double> > >ystrallinvkp=ystrallinv;
                std::vector<std::vector<double>>ystrallchp;
                ystrallchp=ychangep1(ystrallinv, qstrall, k, rkpp, tnum, torder, rjlinv);
                ystrallinvkp[k-1]=ystrallchp;
                int liinm=lall[0];
                double m ;
                if (i==0)
                {
                    m=1;
                }else
                {
                    m=mfrom6_j(rjlinv, i, t, rip, liinm);
                }
                std::vector<int>lk1;
                std::vector<int>lk2;
                lk1=genmvec(lk_1,rk);
                lk2=genmvec(t,rjlinv.sj[k-1]);
                std::vector<int>lk=getsamerange(lk1,lk2);
                // for (int lkk:lk)
                // {
                //     mr=mr+ufrom6_j(rk,t,lk_1,rjlinv.sj[k-1],rkpp,lkk);
                // }
                double q = 1;
                if (i == k - 1)
                {
                    q = 1;
                } else if (i == k - 2)
                {
                    q = qfrom6_j(rjlinv, k - 1, t, lk_1, lall[lall.size()-2]);

                } else
                {
                    int i_cc = 0;
                    for (int i_c = i + 1; i_c <= k - 1; i_c++)
                    {
                        int li_2 = lall[i_cc];
                        int li_1 = lall[i_cc + 1];
                        q = q * qfrom6_j(rjlinv, i_c, t, li_1, li_2);
                        i_cc++;
                    }
                }

                int jk_1;
                if (k==1)
                {
                    jk_1=rjlinv.j;
                }else
                {
                    jk_1=rjlinv.sj[k-2];
                }
                int jk=rjlinv.sj[k-1];
                double g=2*mysqrt(rkpp+1)*ufrom6_j(jk_1,t,jk,rkpp,lk_1,rk)/mysqrt(rk+1);
                std::vector<int>liall=lall;
                basis rjlinvcha=rjlinv;
                if (i==0)
                {
                    rjlinvcha.j=rip;
                    if (!liall.empty())
                    {
                        // 检查是否为空
                        liall.erase(liall.begin()); // 删除第一个元素
                    }else
                    {
                        //std::cout<<"error"<<std::endl;
                    }
                    removeRange(rjlinvcha.sj,1,k-1);
                    rjlinvcha.r[k-1]=rkpp;
                    rjlinvcha.sj.insert(rjlinvcha.sj.begin()+0,liall.begin(),liall.end());

                }else
                    {
                    rjlinvcha.r[i-1]=rip;
                    removeRange(rjlinvcha.sj,i,k-1);
                    rjlinvcha.sj.insert(rjlinvcha.sj.begin()+i-1,liall.begin(),liall.end());
                    }
                rjlinvcha.r[k-1]=rkpp;
                sum=sum+(q*g*m*overlap(rjl,rjlinvcha,ystrall,ystrallinvkp));



                // if (!outfile.is_open())
                // {
                //     // 如果文件未打开，则尝试打开文件
                //     outfile.open("basis_output.txt",std::ios::app);
                // }
                // outfile<<"rkp="<<rkpp<<std::endl;
                // outfile<<"basis1"<<std::endl;
                // printBasisOneLinefile(outfile,rjl);
                // outfile<<"basis2"<<std::endl;
                // printBasisOneLinefile(outfile,rjlinvcha);
                double ov=overlap(rjl, rjlinvcha, ystrall, ystrallinvkp);
                int sizey=ystrall.size();
                // for (int iy=0; iy<sizey; iy++)
                // {
                //     if (!outfile.is_open())
                //     {
                //         // 如果文件未打开，则尝试打开文件
                //         outfile.open("basis_output.txt",std::ios::app);
                //     }
                //     outfile<<"ystrall1"<<"_"<<iy<<std::endl;
                //     writeMatrixToFile(ystrall[iy],outfile);
                //     if (!outfile.is_open())
                //     {
                //         // 如果文件未打开，则尝试打开文件
                //         outfile.open("basis_output.txt",std::ios::app);
                //     }
                //     outfile<<"ystrall2"<<"_"<<iy<<std::endl;
                //     writeMatrixToFile(ystrallinvkp[iy],outfile);
                //
                // }
                // if (!outfile.is_open())
                // {
                //     // 如果文件未打开，则尝试打开文件
                //     outfile.open("basis_output.txt",std::ios::app);
                // }
                // outfile<< "overlap="<<ov<<std::endl;
                // outfile<< "g="<<g<<std::endl;
                // outfile<< "m="<<m<<std::endl;
                // outfile<< "q="<<q<<std::endl;
                // if (!outfile.is_open())
                // {
                //     // 如果文件未打开，则尝试打开文件
                //     outfile.open("basis_output.txt",std::ios::app);
                // }
                // outfile<<"littlesum="<<sum<<std::endl;
            }
        }
        return sum;
    }else
    {
        std::vector<int>li1all1=genmvec(lall.back(),rjlinv.r[ii]);
        std::vector<int>li1all2=genmvec(rjlinv.sj[ii],t);
        std::vector<int>li1all=getsamerange(li1all1,li1all2);
        std::vector<int>lallchange=lall;
        double sum=0;
        for (int li1:li1all)
        {
            lallchange.push_back(li1);
            sum=sum+qt2term1( rjl, rjlinv,ystrall,ystrallinv,
                i, ii+1, k, rkp,lallchange, rip,qstrall,torder, tnum);
            lallchange.pop_back();
        }
        // if (!outfile.is_open())
        // {
        //     // 如果文件未打开，则尝试打开文件
        //     outfile.open("basis_output.txt",std::ios::app);
        // }
        // outfile<<"middlesum="<<sum<<std::endl;
        return sum;
    }

}

double qtqttest(basis rjl, basis rjlinvcha,
            std::vector<std::vector<std::vector<double> > > ystrall,
            std::vector<std::vector<std::vector<double> > > ystrallinv,
            std::vector<std::vector<std::vector<double> > > qstrall, int tnum,
            std::vector<int> qorder)
{
    // if (!outfile.is_open())
    // {
    //     // 如果文件未打开，则尝试打开文件
    //     outfile.open("basis_output.txt",std::ios::app);
    // }
    // outfile<<"--------------begincal_______________________________"<<std::endl;
    // if (!outfile.is_open())
    // {
    //     // 如果文件未打开，则尝试打开文件
    //     outfile.open("basis_output.txt",std::ios::app);
    // }
    // outfile<<"basis1"<<std::endl;
    // printBasisOneLinefile(outfile,rjl);
    // outfile<<"basis2"<<std::endl;
    // printBasisOneLinefile(outfile,rjlinvcha);
    int sizey1=ystrall.size();
    // for (int iy1=0; iy1<sizey1; iy1++)
    // {
    //     if (!outfile.is_open())
    //     {
    //         // 如果文件未打开，则尝试打开文件
    //         outfile.open("basis_output.txt",std::ios::app);
    //     }
    //     outfile<<"ystrall1"<<"_"<<iy1<<std::endl;
    //     writeMatrixToFile(ystrall[iy1],outfile);
    //     if (!outfile.is_open())
    //     {
    //         // 如果文件未打开，则尝试打开文件
    //         outfile.open("basis_output.txt",std::ios::app);
    //     }
    //     outfile<<"ystrall2"<<"_"<<iy1<<std::endl;
    //     writeMatrixToFile(ystrallinv[iy1],outfile);
    //
    // }
    // if (!outfile.is_open())
    // {
    //     // 如果文件未打开，则尝试打开文件
    //     outfile.open("basis_output.txt",std::ios::app);
    // }
    // outfile<<"qstrall_"<<std::endl;
    // writeMatrixToFile(qstrall[0],outfile);

    if (rjlinvcha.sj.back()!=rjl.sj.back())
    {
        return 0;
    }
    basis rjlinv=rjlinvcha;
    int N = rjl.r.size();
    int jinv = rjlinv.j;
    int jinvn = rjlinv.jn;
    int t = qorder[tnum];
    double sum = 0.0;
    int sigt=sign_func({t});
    // if (!outfile.is_open())
    // {
    //     // 如果文件未打开，则尝试打开文件
    //     outfile.open("basis_output.txt",std::ios::app);
    // }
    // outfile<<"term1_______________________________"<<std::endl;
    for (int k = 0; k <= N; k++)
    {
        // if (!outfile.is_open())
        // {
        //     // 如果文件未打开，则尝试打开文件
        //     outfile.open("basis_output.txt",std::ios::app);
        // }
        // outfile<<"k="<<k<<std::endl;
        std::vector<int> rkp;
        if (k == 0)
        {
            if (rjl.j==0)
            {
                continue;
            }
            rkp = genmvec(t, jinv);
            for (int rkpp: rkp)
            {
                int contj = -1;
                std::vector<int> rkp2 = is_get_in(allnucleus, 2, rkpp);
                int par = rjlinvcha.parity;
                if (isOdd(t / 2)) { par = par * (-1); } else { par = par * (1); }
                if (rkp2.empty())
                {
                    continue;
                } else
                {
                    for (const auto &elementnum: rkp2)
                    {
                        if (nucleus[elementnum].parity == par)
                        {
                            contj = elementnum;
                        }
                    }
                    if (contj == -1)
                    {
                        continue;
                    }
                }
                double c = 0;
                double c1 = 0;
                std::vector<int> sig;
                sig = {rkpp, -jinv};
                c = sign_func(sig);
                c = c * mysqrt(rkpp + 1) / mysqrt(jinv + 1);
                c1 = (t + 1) / (mysqrt(jinv + 1) * mysqrt(rkpp + 1));
                c1 = c1 * qstrall[tnum][contj][jinvn] * qstrall[tnum][jinvn][contj];
                sum = sum + c1 * c * overlap(rjl, rjlinv, ystrall, ystrallinv);
                // if (!outfile.is_open())
                // {
                //     // 如果文件未打开，则尝试打开文件
                //     outfile.open("basis_output.txt",std::ios::app);
                // }
                // outfile<<"rkp="<<rkpp<<std::endl;
                // outfile<<"c1="<<c1<<std::endl;
                // outfile<<"c="<<c<<std::endl;
                // outfile<<"basis1"<<std::endl;
                // printBasisOneLinefile(outfile,rjl);
                // outfile<<"basis2"<<std::endl;
                // printBasisOneLinefile(outfile,rjlinv);
                double ov=overlap(rjl, rjlinv, ystrall, ystrallinv);
                int sizey=ystrall.size();
                // for (int iy=0; iy<sizey; iy++)
                // {
                //     if (!outfile.is_open())
                //     {
                //         // 如果文件未打开，则尝试打开文件
                //         outfile.open("basis_output.txt",std::ios::app);
                //     }
                //     outfile<<"ystrall1"<<"_"<<iy<<std::endl;
                //     writeMatrixToFile(ystrall[iy],outfile);
                //     if (!outfile.is_open())
                //     {
                //         // 如果文件未打开，则尝试打开文件
                //         outfile.open("basis_output.txt",std::ios::app);
                //     }
                //     outfile<<"ystrall2"<<"_"<<iy<<std::endl;
                //     writeMatrixToFile(ystrallinv[iy],outfile);
                //
                // }
                // if (!outfile.is_open())
                // {
                //     // 如果文件未打开，则尝试打开文件
                //     outfile.open("basis_output.txt",std::ios::app);
                // }
                // outfile<< "overlap="<<ov<<std::endl;
                // if (!outfile.is_open())
                // {
                //     // 如果文件未打开，则尝试打开文件
                //     outfile.open("basis_output.txt",std::ios::app);
                // }
                // outfile<<"largesum="<<sum<<std::endl;

            }
        } else
        {
            int rk = rjlinv.r[k - 1];
            rkp = genmvec(t, rjlinv.r[k - 1]);
            for (int rkpp: rkp)
            {
                std::vector<std::vector<double> > ystrchange;

                ystrchange = ychangeb1(ystrallinv, qstrall, k, rkpp, tnum, qorder, rjlinv);
                std::vector<std::vector<std::vector<double> > > ystrallch = ystrallinv;
                ystrallch[k - 1] = ystrchange;
                double c = 0;
                std::vector<int> sig;
                sig = {rkpp, -rk};
                c = sign_func(sig);
                c = c * mysqrt(rkpp + 1) / mysqrt(rk + 1);
                sum = sum + c * overlap(rjl, rjlinv, ystrall, ystrallch);
                // if (!outfile.is_open())
                // {
                //     // 如果文件未打开，则尝试打开文件
                //     outfile.open("basis_output.txt",std::ios::app);
                // }
                // outfile<<"rkp="<<rkpp<<std::endl;
                // outfile<<"basis1"<<std::endl;
                // printBasisOneLinefile(outfile,rjl);
                // outfile<<"basis2"<<std::endl;
                // printBasisOneLinefile(outfile,rjlinv);
                double ov=overlap(rjl, rjlinv, ystrall, ystrallch);
                int sizey=ystrall.size();
                // for (int iy=0; iy<sizey; iy++)
                // {
                //     if (!outfile.is_open())
                //     {
                //         // 如果文件未打开，则尝试打开文件
                //         outfile.open("basis_output.txt",std::ios::app);
                //     }
                //     outfile<<"ystrall1"<<"_"<<iy<<std::endl;
                //     writeMatrixToFile(ystrall[iy],outfile);
                //     if (!outfile.is_open())
                //     {
                //         // 如果文件未打开，则尝试打开文件
                //         outfile.open("basis_output.txt",std::ios::app);
                //     }
                //     outfile<<"ystrall2"<<"_"<<iy<<std::endl;
                //     writeMatrixToFile(ystrallch[iy],outfile);
                //
                // }
                // if (!outfile.is_open())
                // {
                //     // 如果文件未打开，则尝试打开文件
                //     outfile.open("basis_output.txt",std::ios::app);
                // }
                // outfile<< "overlap="<<ov<<std::endl;
                // outfile<< "c="<<c<<std::endl;
                // if (!outfile.is_open())
                // {
                //     // 如果文件未打开，则尝试打开文件
                //     outfile.open("basis_output.txt",std::ios::app);
                // }
                // outfile<<"largesum="<<sum<<std::endl;
            }
        }
    }
    // if (!outfile.is_open())
    // {
    //     // 如果文件未打开，则尝试打开文件
    //     outfile.open("basis_output.txt",std::ios::app);
    // }
    // outfile<<"term2"<<"_______________________"<<std::endl;
    for (int k = 1; k <= N; k++)
    {
        // if (!outfile.is_open())
        // {
        //     // 如果文件未打开，则尝试打开文件
        //     outfile.open("basis_output.txt",std::ios::app);
        // }
        // outfile<<"k"<<"="<<k<<std::endl;
        for (int i = 0; i <= k-1; i++)
        {
            rjlinv=rjlinvcha;
            // if (!outfile.is_open())
            // {
            //     // 如果文件未打开，则尝试打开文件
            //     outfile.open("basis_output.txt",std::ios::app);
            // }
            // outfile<<"i"<<"="<<i<<std::endl;
            if (i==0)
        {
            if (rjlinv.j==0)
            {
                continue;
            }
        }
            basis rjlinvp=rjlinvcha;
            std::vector<int> r=rjlinvcha.r;
            std::vector<int> rip;
            std::vector<int> rkp;
            double c1=1;
            double c2=1;
            rkp=genmvec(r[k-1], t);
            std::vector<std::vector<std::vector<double> > > ystrallinvchan=ystrallinv;
            if (i==0)
            {

                rip=genmvec(rjlinvcha.j,t);

                for (int ripp: rip)
                {
                    // if (!outfile.is_open())
                    // {
                    //     // 如果文件未打开，则尝试打开文件
                    //     outfile.open("basis_output.txt",std::ios::app);
                    // }
                    // outfile<<"rip="<<ripp<<std::endl;

                    int contj = -1;
                    std::vector<int> rip2= is_get_in(allnucleus, 2, ripp);
                    int par = rjlinvcha.parity;
                    if (isOdd(t / 2))
                    {
                        par = par * (-1);
                    } else
                    {
                        par = par * (1);
                    }
                    if (rip2.empty())
                    {
                        continue;
                    }
                    else
                    {
                        for (const auto &elementnum: rip2)
                        {
                            if (nucleus[elementnum].parity == par)
                            {
                                contj = elementnum;
                            }
                        }
                        if (contj == -1)
                        {
                            continue;
                        }
                    }
                    double c=0;
                    std::vector<int> sig;
                    sig = {t,-ripp, -rjlinvcha.j};
                    c = sign_func(sig);
                    rjlinv.jn=contj;
                    c1=qstrall[tnum][rjlinvcha.jn][contj];
                    c2=mysqrt(t + 1) / mysqrt(ripp + 1);
                    // if (!outfile.is_open())
                    // {
                    //     // 如果文件未打开，则尝试打开文件
                    //     outfile.open("basis_output.txt",std::ios::app);
                    // }
                    // outfile<<"ripp="<<ripp<<std::endl;
                    // outfile<<"c="<<c<<std::endl;
                    // outfile<<"c1="<<c1<<std::endl;
                    // outfile<<"c2="<<c2<<std::endl;
                    // outfile<<"sigt="<<sigt<<std::endl;
                    sum=sum+sigt*c*c1*c2*qt2term1(rjl,rjlinv,ystrall,ystrallinv,
                            i,i,k,rkp,{ripp},ripp,qstrall,qorder,tnum);
                    // if (!outfile.is_open())
                    // {
                    //     // 如果文件未打开，则尝试打开文件
                    //     outfile.open("basis_output.txt",std::ios::app);
                    // }
                    // outfile<<"largesum="<<sum<<std::endl;
                }
            }else
            {

                rip=genmvec(r[i-1],t);
                for (int ripp: rip)
                {
                    // if (!outfile.is_open())
                    // {
                    //     // 如果文件未打开，则尝试打开文件
                    //     outfile.open("basis_output.txt",std::ios::app);
                    // }
                    // outfile<<"rip="<<ripp<<std::endl;

                    std::vector<int>lii1;
                    lii1=genmvec(rjlinv.sj[i-1],t);
                    std::vector<int>lii2;
                    lii2=lii1;
                    if (i!=1)
                    {
                        lii2=genmvec(rjlinv.sj[i-2],ripp);
                    }
                    else
                    {
                        lii2=genmvec(rjlinv.j,ripp);
                    }
                    std::vector<int>lii=getsamerange(lii1, lii2);
                    if (lii.empty())
                    {
                        continue;
                    }
                    std::vector<std::vector<double>> ystrchangpp=
                        ychangep1(ystrallinv, qstrall, i, ripp, tnum, qorder, rjlinvcha);
                    ystrallinvchan[i-1]=ystrchangpp;
                    for (int lii1_1: lii)
                    {
                        sum=sum+(qt2term1(rjl,rjlinv,ystrall,ystrallinvchan,
                            i,i,k,rkp,{lii1_1},ripp,qstrall,qorder,tnum));
                        // if (!outfile.is_open())
                        // {
                        //     // 如果文件未打开，则尝试打开文件
                        //     outfile.open("basis_output.txt",std::ios::app);
                        // }
                        // outfile<<"largesum="<<sum<<std::endl;
                    }
                }
            }
        }
    }
    sum=sigt*sum;
    // if (!outfile.is_open())
    // {
    //     // 如果文件未打开，则尝试打开文件
    //     outfile.open("basis_output.txt",std::ios::app);
    // }
    // outfile<<"sum="<<sum<<std::endl;
    return sum;
}

double atatterm2(basis rjl, basis rjlinv,std::vector<std::vector<std::vector<double> > > ystrall,
            std::vector<std::vector<std::vector<double> > > ystrallinv,
            std::vector<std::vector<std::vector<double> > > qstrall, int snum,
            std::vector<int> sorder,int s,int k,int t,std::vector<int> lii,std::vector<int> lall,int i,int ii,int
            lk_1,int rip)
{
    int rk;
    rk = rjl.r[k - 1];
    int ri;
    ri=rjl.r[i-1];
    std::vector<int> r = rjl.r;
    std::vector<int> sj = rjl.sj;
    std::vector<int> sjinv = rjlinv.sj;
    std::vector<int> rinv = rjlinv.r;
    int ji=rjl.sj[ii - 1];
    int N1 = rinv.size();
    int N2 = r.size();
    std::vector<int> li1;
    std::vector<int> li2;
    std::vector<int> li;
    double sum = 0;
    std::vector<int> liitem = lii;
    li1 = genmvec(t, ji);
    li2 = li1;
    double rii;
    double jii;
    if (ii==i)
    {
        rii=rip;
        if (i==0)
        {
            jii=0;
        }else if (i==1)
        {
            jii=rjl.j;
        }else
        {
            jii=rjl.sj[i - 2];
        }
    }else
    {
        rii=r[ii-1];
        jii=lii.back();
    }
    li2 = genmvec(rii,jii);

    li = getsamerange(li1, li2);
    if (ii == k - 1)
    {
        if (!outfile.is_open())
        {
            // 如果文件未打开，则尝试打开文件
            outfile.open("basis_output.txt",std::ios::app);
        }
        outfile<<"i="<<i<<std::endl;
        if (is_in_array(li, lk_1))
        {
            double g = gfrom6_j(rjl, k, s, t, lk_1);
            int liinm;
            if (i==k-1)
            {
                liinm=lk_1;
            }else{liinm=lii[0];}
            double m;
            if (i==0)
            {
                m=1;
            }else
            {
                m= mfrom6_j(rjl, i, t, rip, liinm);
            }
            double q = 1;
            if (i == k - 1)
            {
                q = 1;
                liitem.push_back(lk_1);
            } else if (i == k - 2)
            {
                q = qfrom6_j(rjl, k - 1, t, lk_1, lii.back());
                liitem.push_back(lk_1);
            } else
            {
                liitem.push_back(lk_1);
                int i_cc = 0;
                for (int i_c = i + 1; i_c <= k - 1; i_c++)
                {
                    int li_2 = liitem[i_cc];
                    int li_1 = liitem[i_cc + 1];

                    q = q * qfrom6_j(rjl, i_c, t, li_1, li_2);
                    i_cc++;
                }
            }
            double c = g * m * q;
            std::cout << "gmq" << std::endl;
            if (g==0)
            {
                double ggg=1;
            }
            std::vector<std::vector<double> > ystrallchange =
                getz(ystrall, qstrall,
            i, k, rip, rk, ri, rjl, snum+1, s, t);
            std::vector<std::vector<std::vector<double> > > ystrall1 = ystrall;
            std::vector<std::vector<std::vector<double> > > ystrall2 = ystrallinv;
            ystrall1[i - 1] = ystrallchange;
            ystrall1.erase(ystrall1.begin() + k - 1);
            ystrall1.push_back(qstrall[snum]);
            basis rjl1=rjl;
            basis rjl2=rjlinv;
            rjl1.r.erase(rjl1.r.begin()+k-1);
            rjl1.r.push_back(s);
            rjl1.r[i-1]=rip;
            int n=rjl1.sj.size();
            removeRange(rjl1.sj,i,n-1);
            liitem=lii;
            liitem.push_back(lk_1);
            rjl1.sj.insert(rjl1.sj.end()-1, liitem.begin(), liitem.end());
            std::vector<int>lkall1=lall;
            reverse(lkall1.begin(), lkall1.end());
            if (k!=n)
            {
                rjl1.sj.insert(rjl1.sj.end()-1, lkall1.begin(), lkall1.end());
            }

            double ov=overlap(rjl1,rjl2,ystrall1,ystrall2);


            if (!outfile.is_open())
            {
                // 如果文件未打开，则尝试打开文件
                outfile.open("basis_output.txt",std::ios::app);
            }
            outfile<<"basis1"<<std::endl;
            printBasisOneLinefile(outfile,rjl1);
            outfile<<"basis2"<<std::endl;
            printBasisOneLinefile(outfile,rjl2);
            outfile<<"g="<<g<<std::endl;
            outfile<<"m="<<m<<std::endl;
            outfile<<"q="<<q<<std::endl;
            outfile<<"overlap="<<ov<<std::endl;
            outfile<<"ystrall"<<std::endl;
            writeMatrixToFile(ystrall1[0],outfile);
            if (!outfile.is_open())
            {
                // 如果文件未打开，则尝试打开文件
                outfile.open("basis_output.txt",std::ios::app);
            }
            outfile<<"ystrall1_2"<<std::endl;
            writeMatrixToFile(ystrall1[1],outfile);
            if (!outfile.is_open())
            {
                // 如果文件未打开，则尝试打开文件
                outfile.open("basis_output.txt",std::ios::app);
            }
            outfile<<"ystrall1_3"<<std::endl;
            writeMatrixToFile(ystrall1[2],outfile);
            if (!outfile.is_open())
            {
                // 如果文件未打开，则尝试打开文件
                outfile.open("basis_output.txt",std::ios::app);
            }
            outfile<<"ystrall2_1"<<std::endl;
            if (!outfile.is_open())
            {
                // 如果文件未打开，则尝试打开文件
                outfile.open("basis_output.txt",std::ios::app);
            }
            writeMatrixToFile(ystrall2[0],outfile);
            if (!outfile.is_open())
            {
                // 如果文件未打开，则尝试打开文件
                outfile.open("basis_output.txt",std::ios::app);
            }
            outfile<<"ystrall2_2"<<std::endl;
            writeMatrixToFile(ystrall2[1],outfile);
            if (!outfile.is_open())
            {
                // 如果文件未打开，则尝试打开文件
                outfile.open("basis_output.txt",std::ios::app);
            }
            outfile<<"ystrall2_3"<<std::endl;
            writeMatrixToFile(ystrall2[2],outfile);
            outfile.close();


            return c*ov;
        }else
        {
            return 0;
        }

    }else
    {
        for(int li_1: li)
        {

            std::vector<int> lis;
            lis = lii;
            lis.push_back(li_1);
            sum = sum + atatterm2(rjl, rjlinv,ystrall,ystrallinv ,qstrall, snum, sorder,
                                s, k,t, lis,lall,i,ii+1,lk_1,rip);
            lis.pop_back();
        }

        return sum;
    }
}

double atatterm1( basis rjl,basis rjlinv,std::vector<std::vector<std::vector<double> > > ystrall,
            std::vector<std::vector<std::vector<double> > > ystrallinv,
            std::vector<std::vector<std::vector<double> > > qstrall, int snum,
            std::vector<int> sorder,int s,int k, int kk,std::vector<int> lall)
{
    int n=rjl.r.size();
    std::vector<int> lk_1_1;
    std::vector<int> lk_12;
    std::vector<int> lk_1;
    if (kk==k-1)
    {
        if (!outfile.is_open())
        {
            // 如果文件未打开，则尝试打开文件
            outfile.open("basis_output.txt",std::ios::app);
        }
        outfile<<"k="<<k<<std::endl;
        int lk;
        if (k==n)
        {
            lk=rjl.sj.back();
            lk_1_1=genmvec(lk,s);
            lk_12=lk_1_1;
        }else
        {
            lk=lall.back();
            lk_1_1=genmvec(lk,rjl.r[k]);
            lk_12=genmvec(rjl.sj[k-1],s);

        }
        lk_1=getsamerange(lk_1_1,lk_12);
        if (lk_1.empty())
        {
            return 0;
        }
        double sum1=0;
        for (int lk_11: lk_1)
        {

            int ln_1;
            if (k==n)
            {
                ln_1=lk_11;
            }else
            {
                ln_1=lall[0];
            };
            std::vector<int> sigv;
            sigv={rjl.sj.back(),s,-ln_1};
            int sig=sign_func(sigv);
            double c1=mysqrt(ln_1+1);
            int h=1;
            if (k==n)
            {
                h=1;
            }else if (k==n-1)
            {
                h=hfrom6_j(rjl,n,s,ln_1,lk_11);
            }else
            {
                int ic=1;
                for (int i=n;i>=k+1;i--)
                {
                    if (i>k+1)
                    {
                        h=h*hfrom6_j(rjl,i,s,lall[ic],lall[ic+1]);
                    }else
                    {
                        h=h*hfrom6_j(rjl,i,s,lall[ic],lk_11);
                    }


                }
            }
            if (!outfile.is_open())
            {
                // 如果文件未打开，则尝试打开文件
                outfile.open("basis_output.txt",std::ios::app);
            }
            outfile<<"c="<<c1<<std::endl;
            outfile<<"sig="<<sig<<std::endl;
            outfile<<"h="<<h<<std::endl;

            double term1;
            double term2;
            int jk_1=rjl.j;
            if (k!=1)
            {
                jk_1=rjl.sj[k-2];
            }
            int pd1=deltatwo(s,rjl.r[k-1])*deltatwo(lk_11,jk_1);
            if (pd1==0)
            {
                term1=0;
            }else
            {
                double phik=0;
                double sum = 0.0;

                // 遍历 A 和 B 的第一维
                for (size_t a = 0; a < ystrall[k-1].size(); ++a)
                {
                    // 第一维键 r
                    for (size_t b = 0; b < qstrall[snum][a].size(); ++b)
                    {
                        sum=sum+ystrall[k-1][a][b]*qstrall[snum][a][b];
                    }
                }
                int jk = rjl.sj[k - 1];
                int rk = rjl.r[k - 1];
                double cp = 2 * std::pow(-1, (jk - jk_1 - rk)/2.0) * mysqrt(jk + 1) / mysqrt(jk_1 + 1);

                phik=cp * sum;
                basis rjl1=rjl;
                basis rjl2=rjlinv;
                rjl1.r.erase(rjl1.r.begin()+k-1);
                rjl1.r.push_back(s);
                removeRange(rjl1.sj,k,n-1);
                std::vector<int>lkall1=lall;
                reverse(lkall1.begin(), lkall1.end());
                if (k!=n)
                {
                    rjl1.sj.insert(rjl1.sj.end()-1, lkall1.begin(), lkall1.end());
                }
                std::vector<std::vector<std::vector<double> > > ystrall1=ystrall;
                std::vector<std::vector<std::vector<double> > > ystrall2=ystrallinv;
                ystrall1.erase(ystrall1.begin() + k - 1);
                ystrall1.push_back(qstrall[snum]);
                double ov=overlap(rjl1,rjl2,ystrall1,ystrall2);
                term1=ov*phik;

                if (!outfile.is_open())
                {
                    // 如果文件未打开，则尝试打开文件
                    outfile.open("basis_output.txt",std::ios::app);
                }
                outfile<<"basis1"<<std::endl;
                printBasisOneLinefile(outfile,rjl1);
                outfile<<"basis2"<<std::endl;
                printBasisOneLinefile(outfile,rjl2);
                outfile<<"term1="<<term1<<std::endl;
                outfile<<"phik="<<phik<<std::endl;
                outfile<<"overlap="<<ov<<std::endl;
                outfile<<"ystrall"<<std::endl;
                writeMatrixToFile(ystrall1[0],outfile);
                if (!outfile.is_open())
                {
                    // 如果文件未打开，则尝试打开文件
                    outfile.open("basis_output.txt",std::ios::app);
                }
                outfile<<"ystrall1_2"<<std::endl;
                writeMatrixToFile(ystrall1[1],outfile);
                if (!outfile.is_open())
                {
                    // 如果文件未打开，则尝试打开文件
                    outfile.open("basis_output.txt",std::ios::app);
                }
                outfile<<"ystrall1_3"<<std::endl;
                writeMatrixToFile(ystrall1[2],outfile);
                if (!outfile.is_open())
                {
                    // 如果文件未打开，则尝试打开文件
                    outfile.open("basis_output.txt",std::ios::app);
                }
                outfile<<"ystrall2_1"<<std::endl;
                if (!outfile.is_open())
                {
                    // 如果文件未打开，则尝试打开文件
                    outfile.open("basis_output.txt",std::ios::app);
                }
                writeMatrixToFile(ystrall2[0],outfile);
                if (!outfile.is_open())
                {
                    // 如果文件未打开，则尝试打开文件
                    outfile.open("basis_output.txt",std::ios::app);
                }
                outfile<<"ystrall2_2"<<std::endl;
                writeMatrixToFile(ystrall2[1],outfile);
                if (!outfile.is_open())
                {
                    // 如果文件未打开，则尝试打开文件
                    outfile.open("basis_output.txt",std::ios::app);
                }
                outfile<<"ystrall2_3"<<std::endl;
                writeMatrixToFile(ystrall2[2],outfile);
                outfile.close();


            }
            std::vector<int> t;
            t = genmvec(rjl.r[k - 1], s);
            for (int i = k - 1; i >= 0; i = i - 1)
            {
                double sum2=0.0;
                std::vector<int> rip;
                if (i == 0)
                {
                    rip = genmrivec(t, rjl.j);
                    if (rjl.j==0)
                    {
                        continue;
                    }
                } else
                {
                    rip = genmrivec(t, rjl.r[i - 1]);
                }
                for (int ripp: rip)
                {
                    if (!outfile.is_open())
                    {
                        // 如果文件未打开，则尝试打开文件
                        outfile.open("basis_output.txt",std::ios::app);
                    }
                    outfile<<"ripp="<<ripp<<std::endl;
                    int contj = -1;
                    if (i == 0)
                    {
                        std::vector<int> rip2;
                        rip2 = is_get_in(allnucleus, 2, ripp);
                        if (rip2.empty())
                        {
                            continue;
                        } else
                        {
                            for (const auto &elementnum: rip2)
                            {
                                contj = elementnum;
                                /*
                                if (nucleus[elementnum].parity == nucleus[rjl.jn].parity)
                                {
                                    contj = elementnum;
                                }
                                */
                            }
                            if (contj == -1)
                            {
                                continue;
                            }
                        }
                    }
                    std::vector<int> ript;
                    double ri;
                    if (i==0)
                    {
                        ri=rjl.j;
                    }else
                    {
                        ri=rjl.r[i-1];
                    }
                    ript = genmvec(ripp, ri);
                    ript=getsamerange(ript,t);
                    for (int tri: ript)
                    {
                        if (!outfile.is_open())
                        {
                            // 如果文件未打开，则尝试打开文件
                            outfile.open("basis_output.txt",std::ios::app);
                        }
                        outfile<<"t="<<tri<<std::endl;


                        double ji_1=rjl.j;
                        if (i>1)
                        {
                            ji_1 = rjl.sj[i-2];
                        }else if (i==0)
                        {
                            ji_1 = 0;
                        }
                        double sj0;
                        if (i==0)
                        {
                            sj0=rjl.j;
                        }else
                        {
                            sj0=rjl.sj[i-1];
                        }
                        std::vector<int> li1 = genmvec(tri, sj0);

                        std::vector<int> li2 = genmvec(ripp, ji_1);
                        li1 = getsamerange(li1, li2);
                        int ii = i;



                        std::vector<int> lii={};
                        double atatterm22=atatterm2(rjl, rjlinv,ystrall,ystrallinv ,qstrall, snum, sorder,
                            s, k,tri, lii,lall,i,ii,lk_11,ripp);
                             term2= term2 + atatterm22;
                    }
                }
            }
            sum1=sum1+(c1*sig*h*(term1+term2));

        }
        return sum1;
    }else
    {
        std::vector<int> lkk1;
        std::vector<int> lkk2;
        std::vector<int> lkk;
        int l1;
        int r1;
        int jk;
        jk=rjl.sj[kk];
        lkk1=genmvec(jk,s);
        if (lall.empty())
        {
            lkk2=lkk1;
        }else
        {
            l1=lall.back();
            lkk2=genmvec(l1,rjl.r[kk+1]);
        }
        lkk=getsamerange(lkk1,lkk2);
        double sum1=0.0;
        for (int lkk_1: lkk)
        {
            std::vector<int> lallc;
            lallc=lall;
            lallc.push_back(lkk_1);
            sum1=sum1+atatterm1(rjl,rjlinv,ystrall,ystrallinv,qstrall,snum,sorder,s,k,kk-1,lallc);
        }
        return sum1;

    }
}


double atateven(basis rjl,basis rjlinv,std::vector<std::vector<std::vector<double> > > ystrall,
            std::vector<std::vector<std::vector<double> > > ystrallinv,
            std::vector<std::vector<std::vector<double> > > qstrall, int snum,
            std::vector<int> sorder)
{
    if (!outfile.is_open())
    {
        // 如果文件未打开，则尝试打开文件
        outfile.open("basis_output.txt",std::ios::app);
    }
    outfile<<"basisrjl"<<std::endl;
    printBasisOneLinefile(outfile,rjl);
    outfile<<"basisrjlinv"<<std::endl;
    printBasisOneLinefile(outfile,rjlinv);
    double c1=mysqrt(rjl.sj.back()+1);
    int s=sorder[snum];
    double sum=0.0;
    if (rjl.sj.back()!=rjlinv.sj.back())
    {
        return 0;
    }
    int n=rjl.r.size();
    for (int k=1;k<=rjl.r.size();k++)
    {
        sum=sum+atatterm1(rjl,rjlinv,ystrall,ystrallinv,
            qstrall,snum,sorder,s,k,n-1,{});
    }
    sum=sum/c1;

    if (!outfile.is_open())
    {
        // 如果文件未打开，则尝试打开文件
        outfile.open("basis_output.txt",std::ios::app);
    }
    outfile<<"sum="<<sum<<std::endl;
    outfile<<"______________________"
             "_________________________________________________________________________"<<std::endl;

    return sum;
}



/*
double atat1term(basis rjl2,int k,int s,basis rjl1,std::vector<std::vector<std::vector<double> > > ystr1,
    std::vector<std::vector<std::vector<double> > > ystr2,std::vector<int> lii,int kk)
{
    std::vector<int> lii1_1;

    int li_1;
    int N=rjl2.r.size();
    if (k==1)
    {
        li_1=rjl2.j;
    }else if (kk<=N)
    {
        li_1=lii.back();
    }

}
*/
double as2term(basis rjl,basis rjlinv,
std::vector<std::vector<std::vector<double> > > ystrall,
std::vector<std::vector<std::vector<double> > > ystrallinv,
int i,int ii,int k,std::vector<int>lall,int rip,
std::vector<std::vector<std::vector<double>> > qstrall,std::vector<int>sorder,int snum,int t,int contj)
{

    int s=sorder[snum];
    int rk=rjlinv.r[k-1];
    if (ii==k-1)
    {
        int lk_1=0;
        double sum=0;
        if (!lall.empty())
        {
            lk_1 = lall.back();
        }else
        {
            return 0;
        }
        std::vector<int>lk_1line=genmvec(rjlinv.sj[k-1],s);
        if (is_in_array(lk_1line,lk_1))
        {
            int liinm=lall[0];
            double m ;
            if (i==0)
            {
                m=1;
            }else
            {
                m=mfrom6_j(rjlinv, i, t, rip, liinm);
            }
            double q = 1;
            if (i == k - 1)
            {
                q = 1;
            } else if (i == k - 2)
            {
                q = qfrom6_j(rjlinv, k - 1, t, lk_1, lall[lall.size()-2]);
            } else
            {
                int i_cc = 0;
                for (int i_c = i + 1; i_c <= k - 1; i_c++)
                {
                    int li_2 = lall[i_cc];
                    int li_1 = lall[i_cc + 1];
                    i_cc++;
                    q = q * qfrom6_j(rjlinv, i_c, t, li_1, li_2);
                }
            }
            int jk_1;
            if (k==1)
            {
                jk_1=rjlinv.j;
            }else
            {
                jk_1=rjlinv.sj[k-2];
            }
            int jk=rjlinv.sj[k-1];
            int signx;
            std::vector<int>signxline={t,-s,-rk};
            signx=sign_func(signxline);
            double g=mysqrt(t+1)*ufrom6_j(jk_1,t,jk,s,lk_1,rk)/mysqrt(rk+1);
            std::vector<int>liall=lall;
            basis rjlinvcha=rjlinv;
            if (i==0)
            {
                rjlinvcha.j=rip;
                if (!liall.empty())
                {
                    // 检查是否为空
                    liall.erase(liall.begin()); // 删除第一个元素
                }else
                {
                    //std::cout<<"error"<<std::endl;
                }
                removeRange(rjlinvcha.sj,1,k-1);
                rjlinvcha.sj.insert(rjlinvcha.sj.begin()+0,liall.begin(),liall.end());
                rjlinvcha.jn=contj;

            }else
            {
                rjlinvcha.r[i-1]=rip;
                removeRange(rjlinvcha.sj,i,k-1);
                rjlinvcha.sj.insert(rjlinvcha.sj.begin()+i-1,liall.begin(),liall.end());
            }
            rjlinvcha.r[k-1]=s;
            double ov1=overlap(rjl,rjlinvcha,ystrall,ystrallinv);
            // if (!outfile.is_open())
            // {
            //     // 如果文件未打开，则尝试打开文件
            //     outfile.open("basis_output.txt",std::ios::app);
            // }
            //
            // outfile<<"g="<<g<<std::endl;
            // outfile<<"m="<<m<<std::endl;
            // outfile<<"q="<<q<<std::endl;
            // outfile<<"sigx="<<signx<<std::endl;
            // outfile<<"ov1="<<ov1<<std::endl;
            // if (!outfile.is_open())
            // {
            //     // 如果文件未打开，则尝试打开文件
            //     outfile.open("basis_output.txt",std::ios::app);
            // }
            // outfile<<"basis1"<<std::endl;
            // printBasisOneLinefile(outfile,rjl);
            // if (!outfile.is_open())
            // {
            //     // 如果文件未打开，则尝试打开文件
            //     outfile.open("basis_output.txt",std::ios::app);
            // }
            // outfile<<"basis2"<<std::endl;
            // printBasisOneLinefile(outfile,rjlinvcha);
            // int sizey=ystrall.size();
            // for (int iy=0; iy<sizey; iy++)
            // {
            //     if (!outfile.is_open())
            //     {
            //         // 如果文件未打开，则尝试打开文件
            //         outfile.open("basis_output.txt",std::ios::app);
            //     }
            //     outfile<<"ystrall1"<<"_"<<iy<<std::endl;
            //     writeMatrixToFile(ystrall[iy],outfile);
            //     if (!outfile.is_open())
            //     {
            //         // 如果文件未打开，则尝试打开文件
            //         outfile.open("basis_output.txt",std::ios::app);
            //     }
            //     outfile<<"ystrall2"<<"_"<<iy<<std::endl;
            //     writeMatrixToFile(ystrallinv[iy],outfile);
            //
            // }

            sum=sum+(signx*q*g*m*ov1);
            // if (!outfile.is_open())
            // {
            //     // 如果文件未打开，则尝试打开文件
            //     outfile.open("basis_output.txt",std::ios::app);
            // }
            // outfile<<"littlesum="<<sum<<std::endl;
        }
        return sum;
    }else
    {
        std::vector<int>li1all1=genmvec(lall.back(),rjlinv.r[ii]);
        std::vector<int>li1all2=genmvec(rjlinv.sj[ii],t);
        std::vector<int>li1all=getsamerange(li1all1,li1all2);
        std::vector<int>lallchange=lall;
        double sum=0;
        for (int li1:li1all)
        {
            lallchange.push_back(li1);
            sum=sum+as2term( rjl, rjlinv,ystrall,ystrallinv,
            i, ii+1, k,lallchange, rip,qstrall,sorder, snum, t,contj);
            lallchange.pop_back();
            // if (!outfile.is_open())
            // {
            //     // 如果文件未打开，则尝试打开文件
            //     outfile.open("basis_output.txt",std::ios::app);
            // }
            // outfile<<"middlesum="<<sum<<std::endl;
        }
        return sum;
    }

}


double atat(basis rjl, basis rjlinv,
            std::vector<std::vector<std::vector<double> > > ystrall,
            std::vector<std::vector<std::vector<double> > > ystrallinv,
            std::vector<std::vector<std::vector<double> > > qstrall, int snum,
            std::vector<int> sorder)
{
    if (rjlinv.sj.back()!=rjl.sj.back())
    {
        return 0;
    }
    // if (!outfile.is_open())
    // {
    //     // 如果文件未打开，则尝试打开文件
    //     outfile.open("basis_output.txt",std::ios::app);
    // }
    // outfile<<"---------------------begincal-------------------"<<std::endl;
    // if (!outfile.is_open())
    // {
    //     // 如果文件未打开，则尝试打开文件
    //     outfile.open("basis_output.txt",std::ios::app);
    // }
    // outfile<<"basis"<<std::endl;
    // printBasisOneLinefile(outfile,rjl);
    // outfile<<"basisinv"<<std::endl;
    // printBasisOneLinefile(outfile,rjlinv);
    int N = rjl.r.size();
    int jinv = rjlinv.j;
    int jinvn = rjlinv.jn;
    int s = sorder[snum];
    double sum = 0.0;
    int nuclnum=nucleus.size();
    for (int k = 1; k <= N; k++)
    {
        // if (!outfile.is_open())
        // {
        //     // 如果文件未打开，则尝试打开文件
        //     outfile.open("basis_output.txt",std::ios::app);
        // }
        // outfile<<"k="<<k<<std::endl;
        // outfile<<"term1_______________________________"<<std::endl;

        basis rjl1=rjl;
        basis rjl2=rjlinv;
        int rk=rjl2.r[k-1];
        int boolrks=deltatwo(rk,s);
        if (boolrks!=0)
        {
            std::vector<std::vector<std::vector<double> > > ystr1=ystrall;
            std::vector<std::vector<std::vector<double> > > ystr2=ystrallinv;
            double c=0;
            c=mysqrt(s+1);
            double cphi=2/mysqrt(s+1);
            double phi=0;
            for (int i_m=0; i_m < nuclnum; i_m++)
            {
                for (int i_k=0; i_k < nuclnum; i_k++)
                {
                    phi=phi+(ystr2[k-1][i_m][i_k]*qstrall[snum][i_m][i_k]);
                }
            }
            phi=phi*cphi;
            ystr2[k-1]=qstrall[snum];
            rjl2.r[k-1]=s;


            double ov1=overlap(rjl1,rjl2,ystr1,ystr2);
            // if (!outfile.is_open())
            // {
            //     // 如果文件未打开，则尝试打开文件
            //     outfile.open("basis_output.txt",std::ios::app);
            // }
            // outfile<<"c=sqrt(s+1)="<<c<<std::endl;
            // outfile<<"phi="<<phi<<std::endl;
            // outfile<<"ov1="<<ov1<<std::endl;
            // if (!outfile.is_open())
            // {
            //     // 如果文件未打开，则尝试打开文件
            //     outfile.open("basis_output.txt",std::ios::app);
            // }
            // outfile<<"basis1"<<std::endl;
            // printBasisOneLinefile(outfile,rjl1);
            // if (!outfile.is_open())
            // {
            //     // 如果文件未打开，则尝试打开文件
            //     outfile.open("basis_output.txt",std::ios::app);
            // }
            // outfile<<"basis2"<<std::endl;
            // printBasisOneLinefile(outfile,rjl2);
            int sizey=ystr1.size();
            // for (int iy=0; iy<sizey; iy++)
            // {
            //     if (!outfile.is_open())
            //     {
            //         // 如果文件未打开，则尝试打开文件
            //         outfile.open("basis_output.txt",std::ios::app);
            //     }
            //     outfile<<"ystrall1"<<"_"<<iy<<std::endl;
            //     writeMatrixToFile(ystr1[iy],outfile);
            //     if (!outfile.is_open())
            //     {
            //         // 如果文件未打开，则尝试打开文件
            //         outfile.open("basis_output.txt",std::ios::app);
            //     }
            //     outfile<<"ystrall2"<<"_"<<iy<<std::endl;
            //     writeMatrixToFile(ystr2[iy],outfile);
            //
            // }

            sum=sum+c*phi*ov1;
            // if (!outfile.is_open())
            // {
            //     // 如果文件未打开，则尝试打开文件
            //     outfile.open("basis_output.txt",std::ios::app);
            // }
            // outfile<<"largesum="<<sum<<std::endl;
        }

        // if (!outfile.is_open())
        // {
        //     // 如果文件未打开，则尝试打开文件
        //     outfile.open("basis_output.txt",std::ios::app);
        // }
        // outfile<<"term2_______________________________"<<std::endl;

        std::vector<int> t;
        t = genmvec(rjlinv.r[k-1], s);
        std::vector<int> r=rjlinv.r;
        for (int i=k-1;i>=0;i=i-1)
        {
            // if (!outfile.is_open())
            // {
            //     // 如果文件未打开，则尝试打开文件
            //     outfile.open("basis_output.txt",std::ios::app);
            // }
            // outfile<<"i="<<i<<std::endl;
            std::vector<int> rip;
            if (i == 0)
            {
                rip = genmrivec(t, rjlinv.j);
                if (rjlinv.j==0)
                {
                    continue;
                }
            }
            else
            {
                rip=genmrivec(t, r[i - 1]);
            }
            for (int ripp: rip)
            {
                // if (!outfile.is_open())
                // {
                //     // 如果文件未打开，则尝试打开文件
                //     outfile.open("basis_output.txt",std::ios::app);
                // }
                // outfile<<"ripp="<<ripp<<std::endl;
                int contj = -1;
                if (i == 0)
                {
                    std::vector<int> rip2;
                    rip2 = is_get_in(allnucleus, 2, ripp);

                    if (rip2.empty())
                    {
                        continue;
                    }
                }
                std::vector<int> ript;
                double ri;
                if (i==0)
                {
                    ri=rjlinv.j;
                }else
                {
                    ri=r[i-1];
                }
                ript = genmvec(ripp, ri);
                ript=getsamerange(ript,t);
                for (int tri: ript)
                {
                    // if (!outfile.is_open())
                    // {
                    //     // 如果文件未打开，则尝试打开文件
                    //     outfile.open("basis_output.txt",std::ios::app);
                    // }
                    // outfile<<"tri="<<tri<<std::endl;
                    if (i == 0)
                    {
                        std::vector<int> rip2;
                        rip2 = is_get_in(allnucleus, 2, ripp);
                        int par = nucleus[rjlinv.jn].parity;
                        if (isOdd((s+rk) / 2))
                        {
                            par = par * (-1);
                        } else
                        {
                            par = par * (1);
                        }

                        if (rip2.empty())
                        {
                            continue;
                        } else
                        {
                            for (const auto &elementnum: rip2)
                            {

                                contj = elementnum;

                            }
                            if (contj == -1)
                            {
                                continue;
                            }
                        }
                        int c1;
                        std::vector<int> signx;
                        signx={tri,-rjlinv.j,-ripp};
                        c1=sign_func(signx);
                        double c33;
                        c33=c1*4*mysqrt(s+1)*mysqrt(rk+1)*mysqrt(tri+1)/mysqrt(ripp+1);
                        double c2=0;
                        for (int i_m=0;i_m<nuclnum;i_m++)
                        {
                            int j_m=nucleus[i_m].j;
                            double c3=0;
                            c3=ystrallinv[k-1][contj][i_m]*qstrall[snum][i_m][rjlinv.jn];
                            c3=c3*util::wigner_6j(rk,s,tri,rjlinv.j,ripp,j_m);
                            c2=c2+c3;
                        }
                        // if (!outfile.is_open())
                        // {
                        //     // 如果文件未打开，则尝试打开文件
                        //     outfile.open("basis_output.txt",std::ios::app);
                        // }
                        // outfile<<"c2="<<c2<<std::endl;
                        // outfile<<"c33="<<c33<<std::endl;
                        std::vector<std::vector<std::vector<double> > > ystr2_1=ystrallinv;
                        ystr2_1[k-1]=qstrall[snum];
                        sum=sum-c2*c33*as2term( rjl, rjlinv,ystrall,ystr2_1,
                            i, i, k,{ripp}, ripp,qstrall,sorder, snum, tri,contj);
                        // if (!outfile.is_open())
                        // {
                        //     // 如果文件未打开，则尝试打开文件
                        //     outfile.open("basis_output.txt",std::ios::app);
                        // }
                        // outfile<<"largesum="<<sum<<std::endl;

                    }
                    else
                    {
                        basis rjl1_1=rjl;
                        basis rjl2_1=rjlinv;
                        std::vector<std::vector<std::vector<double> > > ystr2_1=ystrallinv;
                        std::vector<std::vector<std::vector<double> > > ystr1_1=ystrall;
                        ystr2_1[k-1]=qstrall[snum];
                        std::vector<std::vector<double> > ystr2_ch=
                            ystrchaa(ystrallinv,qstrall,i,k,ripp,rjlinv,snum,sorder,tri);
                        ystr2_1[i-1]=ystr2_ch;
                        std::vector<int>lii1;
                        lii1=genmvec(rjlinv.sj[i-1],tri);
                        std::vector<int>lii2;
                        lii2=lii1;
                        if (i!=1)
                        {
                            lii2=genmvec(rjlinv.sj[i-2],ripp);
                        }
                        else
                        {
                            lii2=genmvec(rjlinv.j,ripp);
                        }
                        std::vector<int>lii=getsamerange(lii1, lii2);
                        for (int lii1_1: lii)
                        {
                            sum=sum-as2term( rjl, rjlinv,ystrall,ystr2_1,
                            i, i, k,{lii1_1}, ripp,qstrall,sorder, snum, tri,0);
                            // if (!outfile.is_open())
                            // {
                            //     // 如果文件未打开，则尝试打开文件
                            //     outfile.open("basis_output.txt",std::ios::app);
                            // }
                            // outfile<<"largesum="<<sum<<std::endl;
                        }
                    }
                }
            }
        }
    }
    return sum;
}
double sigqtterm(basis rjlinv,basis rjl,
std::vector<std::vector<std::vector<double> > > ystrallinv,
std::vector<std::vector<std::vector<double> > > ystrall,
int k,int knum, int rkp,std::vector<int>lall,
std::vector<std::vector<std::vector<double>> > qstrall,std::vector<int>torder,int tnum)
{
    int t=torder[tnum];
    int rk=rjlinv.j;
    basis rjl2_1=rjlinv;
    int N=rjlinv.r.size();
    int liinm=lall[0];
    if (k!=0)
    {
        rk=rjlinv.r[k-1];
    }
    if (knum==N)
    {
        double sum=0.0;
        int lk=lall.back();
        if (lk!=rjl.sj.back())
        {
            return 0.0;
        }
        double m=1;
        double q;
        std::vector<std::vector<std::vector<double> > >ystrallinvkp=ystrallinv;
        if (k==0)
        {
            m=1;
        }else
        {
            m=mfrom6_j(rjlinv,k,t,rkp,liinm);
        }
        q = 1;
        if (k == N)
        {
            q = 1;
        } else if (k == N - 1)
        {
            q = qfrom6_j(rjlinv, N, t, lk, lall[lall.size()-2]);
        } else
        {
            int i_cc = 0;
            for (int i_c = k + 1; i_c <= N; i_c++)
            {
                int li_2 = lall[i_cc];
                int li_1 = lall[i_cc + 1];
                i_cc++;
                q = q * qfrom6_j(rjlinv, i_c, t, li_1, li_2);
            }
        }
        std::vector<int>liall=lall;
        if (k==0)
        {
            rjl2_1.j=rkp;
            if (!liall.empty())
            {
                // 检查是否为空
                liall.erase(liall.begin()); // 删除第一个元素
            }else
            {
                //std::cout<<"error"<<std::endl;
            }
            removeRange(rjl2_1.sj,1,N);
            rjl2_1.sj.insert(rjl2_1.sj.begin()+0,liall.begin(),liall.end());

        }else
        {

            std::vector<std::vector<double>>ystrallchp;
            ystrallchp=ychangep(ystrallinv, qstrall, k, rkp, tnum, torder, rjlinv);
            double siginych=std::pow(-1, (rk + rkp)/2.0);
            for (auto& row : ystrallchp) {
                for (auto& element : row) {
                    element *= siginych;
                }
            }
            ystrallinvkp[k-1]=ystrallchp;

            rjl2_1.r[k-1]=rkp;
            removeRange(rjl2_1.sj,k,N);
            rjl2_1.sj.insert(rjl2_1.sj.begin()+k-1,liall.begin(),liall.end());
        }
        double re=(q*m*overlap(rjl2_1,rjl,ystrallinvkp,ystrall));

        // std::cout<< "re="<<re<<std::endl;
        // std::cout<< "q="<<q<<std::endl;
        // std::cout<< "m="<<m<<std::endl;

                // if (!outfile.is_open())
                // {
                //     // 如果文件未打开，则尝试打开文件
                //     outfile.open("basis_output.txt",std::ios::app);
                // }
                // outfile<<"k="<<k<<std::endl;
                // outfile<< "re="<<re<<std::endl;
                // outfile<< "q="<<q<<std::endl;
                // outfile<< "m="<<m<<std::endl;
                // outfile<<"basis1"<<std::endl;
                // printBasisOneLinefile(outfile,rjl2_1);
                // outfile<<"basis2"<<std::endl;
                // printBasisOneLinefile(outfile,rjl);
                // outfile<<"ystrall1_1"<<std::endl;
                // writeMatrixToFile(ystrallinvkp[0],outfile);
                // if (!outfile.is_open())
                // {
                //     // 如果文件未打开，则尝试打开文件
                //     outfile.open("basis_output.txt",std::ios::app);
                // }
                // outfile<<"ystrall1_2"<<std::endl;
                // writeMatrixToFile(ystrallinvkp[1],outfile);
                // if (!outfile.is_open())
                // {
                //     // 如果文件未打开，则尝试打开文件
                //     outfile.open("basis_output.txt",std::ios::app);
                // }
                // outfile<<"ystrall2_1"<<std::endl;
                // if (!outfile.is_open())
                // {
                //     // 如果文件未打开，则尝试打开文件
                //     outfile.open("basis_output.txt",std::ios::app);
                // }
                // writeMatrixToFile(ystrall[0],outfile);
                // if (!outfile.is_open())
                // {
                //     // 如果文件未打开，则尝试打开文件
                //     outfile.open("basis_output.txt",std::ios::app);
                // }
                // outfile<<"ystrall2_2"<<std::endl;
                // writeMatrixToFile(ystrall[1],outfile);
                // outfile.close();
        sum=sum+re;
        return sum;
    }else
    {
        std::vector<int>li1all1=genmvec(lall.back(),rjlinv.r[knum]);
        std::vector<int>li1all2=genmvec(rjlinv.sj[knum],t);
        std::vector<int>li1all=getsamerange(li1all1,li1all2);
        std::vector<int>lallchange=lall;
        double sum=0.0;
        for (int li1: li1all)
        {
            lallchange.push_back(li1);
            sum=sum+sigqtterm( rjlinv, rjl,ystrallinv,ystrall,
                k, knum+1, rkp,lallchange,qstrall,torder, tnum);
            lallchange.pop_back();
            // if (!outfile.is_open())
            // {
            //     // 如果文件未打开，则尝试打开文件
            //     outfile.open("basis_output.txt",std::ios::app);
            // }
            // outfile<<"little sum="<<sum<<std::endl;
        }
        return sum;
    }
}
double sigqt(basis rjlinv,basis rjl,
        std::vector<std::vector<std::vector<double> > > ystrallinv,
        std::vector<std::vector<std::vector<double> > > ystrall,
        std::vector<std::vector<std::vector<double> > > qstrall, int tnum,
        std::vector<int> torder)
{
    // if (!outfile.is_open())
    // {
    //     // 如果文件未打开，则尝试打开文件
    //     outfile.open("basis_output.txt",std::ios::app);
    // }
    // outfile<<"begincal"<<std::endl;
    // outfile<<"_________________________________________________________________________________________"<<std::endl;

    // if (!outfile.is_open())
    // {
    //     // 如果文件未打开，则尝试打开文件
    //     outfile.open("basis_output.txt",std::ios::app);
    // }
    // outfile<<"basisrjl"<<std::endl;
    // printBasisOneLinefile(outfile,rjlinv);
    // printBasisOneLinefile(outfile,rjl);



    int t=torder[tnum];

    int num=rjlinv.sj.size();
    double sum=0.0;
    for (int k=0;k<=num;k++)
    {
        if (k==0&& rjl.j==0)
        {
            continue;
        }
        basis rjlinvch=rjlinv;
        std::vector<int> rkp;
        double c1=1;
        double c2=1;
        int rk=rjlinv.j;
        if (k==0)
        {
            if (rjlinv.j==0)
            {
                continue;
            }
        }else
        {
            rk=rjlinv.r[k-1];
        }
        rkp=genmvec(rk, t);
        std::vector<std::vector<std::vector<double> > > ystrallinvchan=ystrallinv;
        if (k==0)
        {

            for (int rkpp:rkp)
            {
                std::vector<int> sigqv={t,rkpp,rk} ;
                int sigq=sign_func(sigqv);
                int contj = -1;
                std::vector<int> rip2= is_get_in(allnucleus, 2, rkpp);
                int par = nucleus[rjlinv.jn].parity;
                if (isOdd(t / 2))
                {
                    par = par * (-1);
                } else
                {
                    par = par * (1);
                }
                if (rip2.empty())
                {
                    continue;
                } else
                {
                    for (const auto &elementnum: rip2)
                    {
                        if (nucleus[elementnum].parity == par)
                        {
                            contj = elementnum;
                        }
                    }
                    if (contj == -1)
                    {
                        continue;
                    }
                }
                rjlinvch.jn=contj;
                std::vector<int> sigaav={t,-rk,-rkpp};
                double sigaa=sign_func(sigaav);
                std::vector<int> siga2v={rk,rkpp};
                double siga2=sign_func(siga2v);
                c1=qstrall[tnum][rjlinv.jn][contj];
                c2=mysqrt(t + 1) / mysqrt(rkpp + 1)*sigaa*siga2;
                // outfile<<"c1="<<c1<<std::endl;
                // outfile<<"c2="<<c2<<std::endl;
                // outfile<<"sigq="<<sigq<<std::endl;
                sum=sum+sigq*c1*c2*sigqtterm( rjlinvch, rjl,ystrallinv,ystrall,
                k, k, rkpp,{rkpp},qstrall,torder, tnum);
                // if (!outfile.is_open())
                // {
                //     // 如果文件未打开，则尝试打开文件
                //     outfile.open("basis_output.txt",std::ios::app);
                // }
                // outfile<<"middle sum2="<<sum<<std::endl;

            }
        }
        else
        {
            for (int rkpp:rkp)
            {
                std::vector<int> sigqv={t,rkpp,rk} ;
                int sigq=sign_func(sigqv);
                std::vector<int>lii1;
                lii1=genmvec(rjlinv.sj[k-1],t);
                std::vector<int>lii2;
                lii2=lii1;
                if (k!=1)
                {
                    lii2=genmvec(rjlinv.sj[k-2],rkpp);
                }
                else
                {
                    lii2=genmvec(rjlinv.j,rkpp);
                }
                std::vector<int>lii=getsamerange(lii1, lii2);
                for (int lii1_1: lii)
                {
                    sum=sum+sigq*sigqtterm( rjlinv, rjl,ystrallinv,ystrall,k,k,rkpp,{lii1_1},qstrall,torder, tnum);
                    // if (!outfile.is_open())
                    // {
                    //     // 如果文件未打开，则尝试打开文件
                    //     outfile.open("basis_output.txt",std::ios::app);
                    // }
                    // outfile<<"middle sum="<<sum<<std::endl;


                }
            }
        }
    }
    // if (!outfile.is_open())
    // {
    //     // 如果文件未打开，则尝试打开文件
    //     outfile.open("basis_output.txt",std::ios::app);
    // }
    // outfile<<"finalsum="<<sum<<std::endl;
    // outfile<<"_________________________________________________________________________________________"<<sum<<std::endl;

    return sum;
}


// std::vector<std::vector<std::vector<double>>>matspl(std::vector<basis>rjl,
//     std::vector<std::vector<std::vector<std::vector<double>>>>ystrall,std::vector<int> torder,
//     std::vector<std::vector<std::vector<double>>>qstrall)
// {
//     outfile.close();
//     std::ofstream outfile("basis_output.txt", std::ios::trunc);
//     int n=rjl.size();
//     std::vector<std::vector<std::vector<double>>>matresult;
//     for (int tnum=0;tnum<torder.size();tnum++)
//     {
//         std::vector<std::vector<double>> mat;
//         mat.assign(n, std::vector<double>(n, 0.0));
//         for (int i=0;i<n;i++)
//         {
//             basis rjl1=rjl[i];
//             std::vector<std::vector<std::vector<double>>> ystra1;
//             ystra1=ystrall[i];
//                 for (int j=0;j<n;j++)
//                 {
//                     basis rjl2=rjl[j];
//                     std::vector<std::vector<std::vector<double>>> ystra2;
//                     ystra2=ystrall[j];
//                     double c2=0;
//                     double c3=0;
//                     c2=-sigqt(rjl1,rjl2,ystra1,ystra2,qstrall,tnum,torder);
//                     c2=c2*std::sqrt(rjl2.sj.back()+1);
//                     mat[i][j]=c2;
//                 }
//
//         }
//         matresult.push_back(mat);
//     }
//     return matresult;
// }
std::mutex mat_mutex; // 用于保护对 matresult 的写入

// 处理单个矩阵元素的计算任务
void process_element(int i, int j, const std::vector<basis>& rjl,
                     const std::vector<std::vector<std::vector<std::vector<double>>>>& ystrall,
                     const std::vector<std::vector<std::vector<double>>>& qstrall,
                     int tnum, const std::vector<int>& torder,
                     std::vector<std::vector<double>>& mat) {
    const basis& rjl1 = rjl[i];
    const std::vector<std::vector<std::vector<double>>>& ystra1 = ystrall[i];
    const basis& rjl2 = rjl[j];
    const std::vector<std::vector<std::vector<double>>>& ystra2 = ystrall[j];

    double c2 = -sigqt(rjl1, rjl2, ystra1, ystra2, qstrall, tnum, torder);
    // c2 *= std::sqrt(rjl2.sj.back() + 1);

    // 使用互斥锁保护对 mat 的写入
    std::lock_guard<std::mutex> lock(mat_mutex);
    mat[i][j] = c2;
}

std::vector<std::vector<std::vector<double>>> matspl(const std::vector<basis>& rjl,
    const std::vector<std::vector<std::vector<std::vector<double>>>>& ystrall,
    const std::vector<int>& torder,
    const std::vector<std::vector<std::vector<double>>>& qstrall) {
    std::ofstream outfile("basis_output.txt", std::ios::trunc);
    int n = rjl.size();
    std::vector<std::vector<std::vector<double>>> matresult;

    for (int tnum = 0; tnum < torder.size(); tnum++) {
        std::vector<std::vector<double>> mat(n, std::vector<double>(n, 0.0));
        std::vector<std::thread> threads;

        // 并行化 i 和 j 的循环
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                threads.push_back(std::thread(process_element, i, j, std::ref(rjl), std::ref(ystrall),
                                              std::ref(qstrall), tnum, std::ref(torder), std::ref(mat)));
            }
        }

        // 等待所有线程完成
        for (auto& t : threads) {
            t.join();
        }

        matresult.push_back(mat);
    }

    outfile.close();
    return matresult;
}

double pvmatcal(basis rjlp,basis rjlpinv,
    int tpnum,std::vector<int> tporder,
    basis rjlv,basis rjlvinv,
    int tvnum,std::vector<int> tvorder,
    int jall,std::vector<std::vector<std::vector<double>>> matsliqn,
    std::vector<std::vector<std::vector<double>>> matsliqp,
    int p1,int p2,int n1,int n2
    )
{
    int jpi=rjlp.sj.back();
    int jpip=rjlpinv.sj.back();
    int jv=rjlv.sj.back();
    int jvp=rjlvinv.sj.back();
    int signx;
    int t=tporder[tpnum];
    std::vector<int> signxline={jpip,jv,jall,t};
    signx=sign_func(signxline);
    double c1=util::wigner_6j(jvp,jpip,jall,jpi,jv,t)*mysqrt(jvp+1)*mysqrt(jpip+1);
    double c2=matsliqn[tvnum][n1][n2];
    double c3=matsliqp[tvnum][p1][p2];
    return signx*c1*c2*c3;
}

// std::vector<std::vector<double>> matcalsim(std::vector<basis>rjl,
//     std::vector<std::vector<std::vector<std::vector<double>>>>ystrall,
//     std::vector<int> torder,std::vector<int>sorder,std::vector<double> energy,
//     std::vector<std::vector<std::vector<double>>>qstrall,std::vector<double> kt,
//     std::vector<std::vector<std::vector<double>>>sstrall,std::vector<double> gs)
// {
//     int n=rjl.size();
//     std::vector<std::vector<double>> mat;
//     mat.assign(n, std::vector<double>(n, 0.0));
//
//     for (int i=0;i<n;i++)
//     {
//         basis rjl1=rjl[i];
//         std::vector<std::vector<std::vector<double>>> ystra1;
//         ystra1=ystrall[i];
//         for (int j=0;j<n;j++)
//         {
//             // if (!outfile.is_open())
//             // {
//             //     // 如果文件未打开，则尝试打开文件
//             //     outfile.open("basis_output.txt",std::ios::app);
//             // }
//             // outfile<<"________________matcal_______________________"<<std::endl;
//             // outfile<<"num1="<<i<<"num2="<<j<<std::endl;
//
//             basis rjl2=rjl[j];
//             std::vector<std::vector<std::vector<double>>> ystra2;
//             ystra2=ystrall[j];
//             if (rjl1.sj.back()!=rjl2.sj.back())
//             {
//                 continue;
//             }
//
//             double c1=0;
//                 c1=h0cal(rjl1,rjl2,ystra1,ystra2,energy);
//
//             double c2=0;
//             double c3=0;
//             for (int tnum=0;tnum<torder.size();tnum++)
//             {
//                 c2=c2+(kt[tnum]*qtqttest(rjl1,rjl2,ystra1,ystra2,qstrall,tnum,torder));
//             }
//             for (int snum=0;snum<sorder.size();snum++)
//             {
//                 c3=c3+(gs[snum]*atat(rjl1,rjl2,ystra1,ystra2,sstrall,snum,sorder));
//             }
//             mat[i][j]=c1+c2+c3;
//         }
//     }
//     return mat;
// }

// std::vector<std::vector<double>> matcalsim(std::vector<basis> rjl,
//     std::vector<std::vector<std::vector<std::vector<double>>>> ystrall,
//     std::vector<int> torder, std::vector<int> sorder, std::vector<double> energy,
//     std::vector<std::vector<std::vector<double>>> qstrall, std::vector<double> kt,
//     std::vector<std::vector<std::vector<double>>> sstrall, std::vector<double> gs)
// {
//     int n = rjl.size();
//     std::vector<std::vector<double>> mat;
//     mat.assign(n, std::vector<double>(n, 0.0));
//
//     // 并行化外层循环
// #pragma omp parallel for collapse(2) // 使用collapse(2)来并行化嵌套循环
//     for (int i = 0; i < n; i++)
//     {
//         for (int j = 0; j < n; j++)
//         {
//             basis rjl1 = rjl[i];
//             std::vector<std::vector<std::vector<double>>> ystra1 = ystrall[i];
//
//             basis rjl2 = rjl[j];
//             std::vector<std::vector<std::vector<double>>> ystra2 = ystrall[j];
//
//             if (rjl1.sj.back() != rjl2.sj.back())
//             {
//                 continue;
//             }
//
//             // double c1 = h0cal(rjl1, rjl2, ystra1, ystra2, energy);
//             double c1=0;
//             double c2 = 0;
//             double c3 = 0;
//
//             for (int tnum = 0; tnum < torder.size(); tnum++)
//             {
//                 c2 += (kt[tnum] * qtqttest(rjl1, rjl2, ystra1, ystra2, qstrall, tnum, torder));
//             }
//
//             // for (int snum = 0; snum < sorder.size(); snum++)
//             // {
//             //     c3 += (gs[snum] * atat(rjl1, rjl2, ystra1, ystra2, sstrall, snum, sorder));
//             // }
//
//             mat[i][j] = c1 + c2 + c3;
//         }
//     }
//
//     return mat;
// }

std::mutex output_mutex; // 用于同步输出

// 处理计算的任务
void process_task(int i, int j, std::vector<basis>& rjl,
                  std::vector<std::vector<std::vector<std::vector<double>>>>& ystrall,
                  std::vector<int>& torder, std::vector<int>& sorder, std::vector<double>& energy,
                  std::vector<std::vector<std::vector<double>>>& qstrall, std::vector<double>& kt,
                  std::vector<std::vector<std::vector<double>>>& sstrall, std::vector<double>& gs,
                  std::vector<std::vector<double>>& mat) {

    basis rjl1 = rjl[i];
    std::vector<std::vector<std::vector<double>>> ystra1 = ystrall[i];
    basis rjl2 = rjl[j];
    std::vector<std::vector<std::vector<double>>> ystra2 = ystrall[j];

    // 如果条件不满足，跳过
    if (rjl1.sj.back() != rjl2.sj.back()) {
        return;
    }

    double c1 = 0, c2 = 0, c3 = 0;
    c1 = h0cal(rjl1, rjl2, ystra1, ystra2, energy); // Uncomment and implement if needed
    for (int tnum = 0; tnum < torder.size(); tnum++) {
        c2 += (kt[tnum] * qtqttest(rjl1, rjl2, ystra1, ystra2, qstrall, tnum, torder));
    }
    // 并行化计算部分（可以添加更多并行计算任务）
    for (int snum = 0; snum < sorder.size(); snum++) {
        c3 += (gs[snum] * atat(rjl1, rjl2, ystra1, ystra2, sstrall, snum, sorder));
    }

    mat[i][j] = c1 + c2 + c3; // 计算并存储到矩阵

    // 确保线程安全的输出
    {
        std::lock_guard<std::mutex> lock(output_mutex);

    }
}

std::vector<std::vector<double>> matcalsim(std::vector<basis> rjl,
    std::vector<std::vector<std::vector<std::vector<double>>>> ystrall,
    std::vector<int> torder, std::vector<int> sorder, std::vector<double> energy,
    std::vector<std::vector<std::vector<double>>> qstrall, std::vector<double> kt,
    std::vector<std::vector<std::vector<double>>> sstrall, std::vector<double> gs) {
    int n = rjl.size();
    std::vector<std::vector<double>> mat(n, std::vector<double>(n, 0.0));

    std::vector<std::thread> threads;

    // 并行化外层循环
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            threads.push_back(std::thread(process_task, i, j, std::ref(rjl), std::ref(ystrall),
                                          std::ref(torder), std::ref(sorder), std::ref(energy),
                                          std::ref(qstrall), std::ref(kt), std::ref(sstrall),
                                          std::ref(gs), std::ref(mat)));
        }
    }

    // 等待所有线程完成
    for (auto& t : threads) {
        t.join();
    }

    return mat;
}
// std::vector<std::vector<double>> matcalsim(std::vector<basis> rjl,
//     std::vector<std::vector<std::vector<std::vector<double>>>> ystrall,
//     std::vector<int> torder, std::vector<int> sorder, std::vector<double> energy,
//     std::vector<std::vector<std::vector<double>>> qstrall, std::vector<double> kt,
//     std::vector<std::vector<std::vector<double>>> sstrall, std::vector<double> gs) {
//
//     int n = rjl.size();
//     std::vector<std::vector<double>> mat(n, std::vector<double>(n, 0.0));
//
//     // 并行化 i 和 j 双层循环
// #pragma omp parallel for collapse(2)
//     for (int i = 0; i < n; ++i) {
//         for (int j = 0; j < n; ++j) {
//             basis rjl1 = rjl[i];
//             std::vector<std::vector<std::vector<double>>> ystra1 = ystrall[i];
//             basis rjl2 = rjl[j];
//             std::vector<std::vector<std::vector<double>>> ystra2 = ystrall[j];
//
//             if (rjl1.sj.back() != rjl2.sj.back()) {
//                 continue;
//             }
//
//             double c1 = 0, c2 = 0, c3 = 0;
//             c1 = h0cal(rjl1, rjl2, ystra1, ystra2, energy);
//
//             for (int tnum = 0; tnum < torder.size(); tnum++) {
//                 c2 += (kt[tnum] * qtqttest(rjl1, rjl2, ystra1, ystra2, qstrall, tnum, torder));
//             }
//
//             for (int snum = 0; snum < sorder.size(); snum++) {
//                 c3 += (gs[snum] * atat(rjl1, rjl2, ystra1, ystra2, sstrall, snum, sorder));
//             }
//
//             mat[i][j] = c1 + c2 + c3;
//         }
//     }
//
//     return mat;
// }


// std::vector<std::vector<double>> overlapmatcal(std::vector<basis>&rjl,
//     std::vector<std::vector<std::vector<std::vector<double>>>>&ystrall)
// {
//
//
//     const int n=rjl.size();
//     std::vector<std::vector<double>> overlapmat(n, std::vector<double>(n, 0.0));
//     for (int i=0;i<n;i++)
//     {
//         std::cout<<'i'<<i<<std::endl;
//         for (int j=0;j<n;j++)
//         {
//             // if (!outfile.is_open())
//             // {
//             //     // 如果文件未打开，则尝试打开文件
//             //     outfile.open("basis_output.txt",std::ios::app);
//             // }
//             // outfile<<"____________________begin_____________________"<<std::endl;
//             // outfile<<"rjl1"<<std::endl;
//             // printBasisOneLinefile(outfile,rjl[i]);
//             // outfile<<"rjl2"<<std::endl;
//             // printBasisOneLinefile(outfile,rjl[j]);
//             // if (!outfile.is_open())
//             // {
//             //     // 如果文件未打开，则尝试打开文件
//             //     outfile.open("basis_output.txt",std::ios::app);
//             // }
//             // outfile<<"ystrall1_1"<<std::endl;
//             // writeMatrixToFile(ystrall[i][0],outfile);
//             // if (!outfile.is_open())
//             // {
//             //     // 如果文件未打开，则尝试打开文件
//             //     outfile.open("basis_output.txt",std::ios::app);
//             // }
//             // outfile<<"ystrall1_2"<<std::endl;
//             // writeMatrixToFile(ystrall[i][1],outfile);
//             // if (!outfile.is_open())
//             // {
//             //     // 如果文件未打开，则尝试打开文件
//             //     outfile.open("basis_output.txt",std::ios::app);
//             // }
//             // outfile<<"ystrall1_3"<<std::endl;
//             // writeMatrixToFile(ystrall[i][2],outfile);
//             // if (!outfile.is_open())
//             // {
//             //     // 如果文件未打开，则尝试打开文件
//             //     outfile.open("basis_output.txt",std::ios::app);
//             // }
//             // outfile<<"ystrall2_1"<<std::endl;
//             // writeMatrixToFile(ystrall[j][0],outfile);
//             // if (!outfile.is_open())
//             // {
//             //     // 如果文件未打开，则尝试打开文件
//             //     outfile.open("basis_output.txt",std::ios::app);
//             // }
//             // outfile<<"ystrall2_2"<<std::endl;
//             // writeMatrixToFile(ystrall[j][1],outfile);
//             // if (!outfile.is_open())
//             // {
//             //     // 如果文件未打开，则尝试打开文件
//             //     outfile.open("basis_output.txt",std::ios::app);
//             // }
//             // outfile<<"ystrall2_3"<<std::endl;
//             // writeMatrixToFile(ystrall[j][2],outfile);
//
//             double c=0;
//             if (i==2 && j==2)
//             {
//                 double mmmm=1;
//             }
//             c=overlap(rjl[i],rjl[j],ystrall[i],ystrall[j]);
//             // if (!outfile.is_open())
//             // {
//             //     // 如果文件未打开，则尝试打开文件
//             //     outfile.open("basis_output.txt",std::ios::app);
//             // }
//             // outfile<<"all result"<< c <<std::endl;
//             // outfile<<"____________________end_____________________"<<std::endl;
//             std::cout<<'j'<<j<<std::endl;
//
//             overlapmat[i][j]=c;
//         }
//     }
//
//     return overlapmat;
// }

std::mutex output_mutex2; // 用于同步输出


// 处理每一对 (i, j) 的计算任务
void process_overlap(int i, int j, std::vector<basis>& rjl,
                     std::vector<std::vector<std::vector<std::vector<double>>>>& ystrall,
                     std::vector<std::vector<double>>& overlapmat) {
    basis rjl1 = rjl[i];
    std::vector<std::vector<std::vector<double>>> ystra1 = ystrall[i];
    basis rjl2 = rjl[j];
    std::vector<std::vector<std::vector<double>>> ystra2 = ystrall[j];
    if (rjl1.sj.back() != rjl2.sj.back()) {
        return;
    }

    double c = overlap(rjl1 , rjl2 , ystra1, ystra2);
    overlapmat[i][j] = c;
    // 确保线程安全的输出
    {
        std::lock_guard<std::mutex> lock(output_mutex2);
    }
}

// 并行计算重叠矩阵的函数
std::vector<std::vector<double>> overlapmatcal(std::vector<basis>& rjl,
                                               std::vector<std::vector<std::vector<std::vector<double>>>>& ystrall) {
    const int n = rjl.size();
    std::vector<std::vector<double>> overlapmat(n, std::vector<double>(n, 0.0));
    std::vector<std::thread> threads;

    // 并行化 (i, j) 双重循环
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            threads.push_back(std::thread(process_overlap, i, j, std::ref(rjl), std::ref(ystrall), std::ref(overlapmat)));
        }
    }
    // 等待所有线程完成
    for (auto& t : threads) {
        t.join();
    }
    return overlapmat;
}
// std::vector<std::vector<double>> overlapmatcal(std::vector<basis>& rjl,
//     std::vector<std::vector<std::vector<std::vector<double>>>>& ystrall) {
//
//     const int n = rjl.size();
//     std::vector<std::vector<double>> overlapmat(n, std::vector<double>(n, 0.0));
//
//     // 并行化外层和内层循环
// #pragma omp parallel for collapse(2)
//     for (int i = 0; i < n; i++) {
//         for (int j = 0; j < n; j++) {
//             double c = 0;
//
//             // 可以在此处保留调试信息，但最好避免过多的 I/O 操作
//             if (i == 2 && j == 2) {
//                 double mmmm = 1;
//             }
//
//             // 计算重叠矩阵的元素
//             c = overlap(rjl[i], rjl[j], ystrall[i], ystrall[j]);
//
//             // 保存结果到矩阵
//             overlapmat[i][j] = c;
//
//             // 在并行代码中避免过多的输出操作，如果需要调试，可以使用锁或延迟输出
//             // 如果需要输出调试信息，可以选择将信息收集起来，在循环外统一输出
//             // std::cout << 'i' << i << std::endl;
//             // std::cout << 'j' << j << std::endl;
//         }
//     }
//
//     return overlapmat;
// }

std::map<int,std::vector<std::vector<double>>> cpmatcal(std::vector<std::vector<double>>slimmatp,
    std::vector<std::vector<double>>slimmatn,std::map<int, std::vector<CoupledBasis>>cpbasisall,
    std::vector<basis>rjln,std::vector<basis>rjlp,
    std::vector<int> tporder,std::vector<int> tvorder,
    std::vector<std::vector<std::vector<double>>> matsliqn,
    std::vector<std::vector<std::vector<double>>> matsliqp,std::vector<double> ktn
    )
{
    std::map<int,std::vector<std::vector<double>>> cpmatcalresult;
    for (const auto& [key, coupledBasesVector] : cpbasisall)
    {
        int n=cpbasisall[key].size();
        std::vector<std::vector<double>> zeroMatrix(n, std::vector<double>(n, 0.0));
        for (int i=0;i<n;i++)
        {
            int num1p=coupledBasesVector[i].pi;
            int num1n=coupledBasesVector[i].ni;
            basis rjlp1=rjlp[num1p];
            basis rjln1=rjln[num1n];



            for (int j=0;j<n;j++)
            {

                int num2p=coupledBasesVector[j].pi;
                int num2n=coupledBasesVector[j].ni;
                double c1=slimmatp[num1p][num2p]+slimmatn[num1n][num2n];
                basis rjlp2=rjlp[num2p];
                basis rjln2=rjln[num2n];
                double ctow=0;
                for (int tpnum=0;tpnum<tporder.size();tpnum++)
                {
                    ctow=ctow+(ktn[tpnum]*pvmatcal( rjlp1,rjlp2,tpnum,
                    tporder,rjln1, rjln2,tpnum,tvorder,key,matsliqn,
                    matsliqp,num1p,num2p,num1n,num2n));
                }

                zeroMatrix[i][j]=ctow+c1;
            }
        }
        cpmatcalresult[key]=zeroMatrix;
    }
    return cpmatcalresult;
}
