
//
// Created by wang- on 2025/2/17.
//
#include <iostream>
#include <fstream>
#include <array>
#include <cmath>
#include <iomanip>
#include "mtcal.h"
#include "qcal.h"
#include <lambda_lanczos/lambda_lanczos.hpp>
#include <thread>
#include <vector>
#include <functional>
#include <chrono>  // 引入 chrono 库

// double getzinq (const std::vector<std::vector<std::vector<double> > > &ystrall1,
//     std::vector<std::vector<std::vector<double> > > &qstrall,
//     int k, int rkp, int tnum, std::vector<int> torder, basis rjl,int a,int b,int d)
// {
//     int nucnum = nucleus.size();
//     double t=torder[tnum];
//     int rk=rjl.r[k-1];
//     double c1=mysqrt(rk+1)*mysqrt(t+1);
//     double sum=0.0;
//
//     for(int d0=0;d0<nucnum;d0++)
//     {
//         int j_a = nucleus[a].j;
//         int j_b = nucleus[b].j;
//         int j_d = nucleus[d].j;
//         int j_d0=nucleus[d0].j;
//         double re=ystrall1[k-1][d][d0];
//         re=re*qstrall[tnum][d0][b];
//         re=re*util::wigner_6j(rk,t,rkp,j_b,j_d,j_d0);
//         sum=sum+re;
//     }
//     sum=sum*c1;
//     return sum;
//
// }
//
// double getyinq (const std::vector<std::vector<std::vector<double> > > &ystrall1,
//     std::vector<std::vector<std::vector<double> > > &qstrall,
//     int k, int rkp, int tnum, std::vector<int> torder, basis rjl,int a,int b,int d)
// {
//     int nucnum = nucleus.size();
//     double t=torder[tnum];
//     double y1=getzinq(ystrall1,qstrall,k,rkp,tnum,torder,rjl,a,b,d);
//     double y2=getzinq(ystrall1,qstrall,k,rkp,tnum,torder,rjl,a,d,b);
//     int j_b = nucleus[b].j;
//     int j_d = nucleus[d].j;
//     std::vector<int> signorder={j_b,j_d,rkp};
//     int sign=sign_func(signorder);
//     double re=y1-(sign*y2);
//     return re;
//
// }
//
//
//
// double getzinq (const std::vector<std::vector<std::vector<double> > > &ystrall1,
//     std::vector<std::vector<std::vector<double> > > &qstrall,
//     int k, int rkp, int tnum, std::vector<int> torder, basis rjl,int a,int b)
// {
//     int nucnum = nucleus.size();
//     double t=torder[tnum];
//     double c1=mysqrt(rkp+1)*mysqrt(t+1);
//     double sum=0.0;
//     int rk=rjl.r[k-1];
//
//     for (int d=0; d<nucnum; d++)
//     {
//         int j_a = nucleus[a].j;
//         int j_b = nucleus[b].j;
//         int j_d = nucleus[d].j;
//         double y1=getyinq(ystrall1,qstrall,k,rkp,tnum,torder,rjl,a,b,d);
//         double re=y1*qstrall[tnum][d][a];
//         re=re*util::wigner_6j(rkp,t,rk,j_a,j_b,j_d);
//         sum=sum+re;
//     }
//     sum=sum*c1;
//     return sum;
// }
//
//
// std::vector<std::vector<double> > ychangeb1(
//     const std::vector<std::vector<std::vector<double> > > &ystrall1,
//     std::vector<std::vector<std::vector<double> > > &qstrall,
//     int k, int rkp, int tnum, std::vector<int> torder, basis rjl)
// {
//     int nucnum = nucleus.size();
//     int N2 = rjl.r.size();
//     int rk = rjl.r[k - 1];
//     int t = torder[tnum];
//     double c4 =  mysqrt(t + 1) * mysqrt(rkp + 1);
//     std::vector<std::vector<double> > ychange;
//     ychange.resize(nucnum);
//     for (int i_m = 0; i_m < nucnum; i_m++)
//     {
//         ychange[i_m].resize(nucnum, 0.0); // 初始化二维数组的每一行
//     }
//     for (int a = 0; a < nucnum; a++)
//     {
//         for (int b = 0; b < nucnum; b++)
//         {
//             int j_a = nucleus[a].j;
//             int j_b = nucleus[b].j;
//             std::vector<int> signorder={j_a,j_b,rk};
//             int sign=sign_func(signorder);
//             double z1=getzinq(ystrall1,qstrall,k,rkp,tnum,torder,rjl,a,b);
//             double z2=getzinq(ystrall1,qstrall,k,rkp,tnum,torder,rjl,b,a);
//             double re=z1-(sign*z2);
//             ychange[a][b]=re;
//         }
//     }
//     return ychange;
// }
//
// std::vector<std::vector<double> > ychangep1(
//     std::vector<std::vector<std::vector<double> > > ystrall1,
//     std::vector<std::vector<std::vector<double> > > qstrall,
//     int k, int rkp, int tnum, std::vector<int> torder, basis rjl)
// {
//     int nucnum = nucleus.size();
//     int N2 = rjl.r.size();
//     int rk = rjl.r[k - 1];
//     int t = torder[tnum];
//     double c4 = mysqrt(rk + 1) * mysqrt(t + 1) ;
//     std::vector<std::vector<double>> ychange(nucnum, std::vector<double>(nucnum, 0.0));
//
//     for (int a = 0; a < nucnum; a++)
//     {
//         for (int b=0;b<nucnum;b++)
//         {
//             double re=getyinq(ystrall1,qstrall,k,rkp,tnum,torder,rjl,1,a,b);
//             ychange[a][b]=re;
//         }
//     }
//     return ychange;
// }
//
//
//
// double qt2term1(basis rjl,basis rjlinv,
// std::vector<std::vector<std::vector<double> > > ystrall,
// std::vector<std::vector<std::vector<double> > > ystrallinv,
// int i,int ii,int k, std::vector<int>rkp,std::vector<int>lall,int rip,
// std::vector<std::vector<std::vector<double>> > qstrall,std::vector<int>torder,int tnum)
// {
//     int t=torder[tnum];
//     int rk=rjlinv.r[k-1];
//     if (ii==k-1)
//     {
//         int lk_1=0;
//         double sum=0;
//         for (int rkpp:rkp)
//         {
//             if (!lall.empty())
//             {
//                  lk_1 = lall.back();
//             }else
//             {
//                 return 0;
//             }
//             std::vector<int>lk_1line=genmvec(rjlinv.sj[k-1],rkpp);
//             if (is_in_array(lk_1line,lk_1))
//             {
//                 std::vector<std::vector<std::vector<double> > >ystrallinvkp=ystrallinv;
//                 std::vector<std::vector<double>>ystrallchp;
//                 ystrallchp=ychangep1(ystrallinv, qstrall, k, rkpp, tnum, torder, rjlinv);
//                 ystrallinvkp[k-1]=ystrallchp;
//                 int liinm=lall[0];
//                 double m ;
//                 if (i==0)
//                 {
//                     m=1;
//                 }else
//                 {
//                     m=mfrom6_j(rjlinv, i, t, rip, liinm);
//                 }
//                 std::vector<int>lk1;
//                 std::vector<int>lk2;
//                 lk1=genmvec(lk_1,rk);
//                 lk2=genmvec(t,rjlinv.sj[k-1]);
//                 std::vector<int>lk=getsamerange(lk1,lk2);
//                 double mr=1;
//                 // for (int lkk:lk)
//                 // {
//                 //     mr=mr+ufrom6_j(rk,t,lk_1,rjlinv.sj[k-1],rkpp,lkk);
//                 // }
//                 double q = 1;
//                 if (i == k - 1)
//                 {
//                     q = 1;
//                 } else if (i == k - 2)
//                 {
//                     q = qfrom6_j(rjl, k - 1, t, lk_1, lall.back());
//                 } else
//                 {
//                     int i_cc = 0;
//                     for (int i_c = i + 1; i_c <= k - 1; i_c++)
//                     {
//                         int li_2 = lall[i_cc];
//                         int li_1 = lall[i_cc + 1];
//                         q = q * qfrom6_j(rjl, i_c, t, li_1, li_2);
//                         i_cc++;
//                     }
//                 }
//
//                 int jk_1;
//                 if (k==1)
//                 {
//                     jk_1=rjlinv.j;
//                 }else
//                 {
//                     jk_1=rjlinv.sj[k-2];
//                 }
//                 int jk=rjlinv.sj[k-1];
//                 double g=2*mysqrt(rkpp+1)*ufrom6_j(jk_1,t,jk,rkpp,lk_1,rk)/mysqrt(rk+1);
//                 std::vector<int>liall=lall;
//                 basis rjlinvcha=rjlinv;
//                 if (i==0)
//                 {
//                     rjlinvcha.j=rip;
//                     if (!liall.empty())
//                     {
//                         // 检查是否为空
//                         liall.erase(liall.begin()); // 删除第一个元素
//                     }else
//                     {
//                         //std::cout<<"error"<<std::endl;
//                     }
//                     removeRange(rjlinvcha.sj,1,k-1);
//                     rjlinvcha.r[k-1]=rkpp;
//                     rjlinvcha.sj.insert(rjlinvcha.sj.begin()+0,liall.begin(),liall.end());
//
//                 }else
//                     {
//                     rjlinvcha.r[i-1]=rip;
//                     removeRange(rjlinvcha.sj,i,k-1);
//                     rjlinvcha.sj.insert(rjlinvcha.sj.begin()+i-1,liall.begin(),liall.end());
//                     }
//                 rjlinvcha.r[k-1]=rkpp;
//                 mr=1;
//                 sum=sum+(q*g*m*mr*overlap(rjl,rjlinvcha,ystrall,ystrallinvkp));
//
//
//                 if (!outfile.is_open())
//                 {
//                     // 如果文件未打开，则尝试打开文件
//                     outfile.open("basis_output.txt",std::ios::app);
//                 }
//                 outfile<<"rkp="<<rkpp<<std::endl;
//                 outfile<<"basis1"<<std::endl;
//                 printBasisOneLinefile(outfile,rjl);
//                 outfile<<"basis2"<<std::endl;
//                 printBasisOneLinefile(outfile,rjlinvcha);
//                 double ov=overlap(rjl, rjlinvcha, ystrall, ystrallinvkp);
//                 int sizey=ystrall.size();
//                 for (int iy=0; iy<sizey; iy++)
//                 {
//                     if (!outfile.is_open())
//                     {
//                         // 如果文件未打开，则尝试打开文件
//                         outfile.open("basis_output.txt",std::ios::app);
//                     }
//                     outfile<<"ystrall1"<<"_"<<iy<<std::endl;
//                     writeMatrixToFile(ystrall[iy],outfile);
//                     if (!outfile.is_open())
//                     {
//                         // 如果文件未打开，则尝试打开文件
//                         outfile.open("basis_output.txt",std::ios::app);
//                     }
//                     outfile<<"ystrall2"<<"_"<<iy<<std::endl;
//                     writeMatrixToFile(ystrallinvkp[iy],outfile);
//
//                 }
//                 if (!outfile.is_open())
//                 {
//                     // 如果文件未打开，则尝试打开文件
//                     outfile.open("basis_output.txt",std::ios::app);
//                 }
//                 outfile<< "overlap="<<ov<<std::endl;
//                 outfile<< "g="<<g<<std::endl;
//                 outfile<< "m="<<m<<std::endl;
//                 outfile<< "q="<<q<<std::endl;
//             }
//         }
//         return sum;
//     }else
//     {
//         std::vector<int>li1all1=genmvec(lall.back(),rjlinv.r[ii]);
//         std::vector<int>li1all2=genmvec(rjlinv.sj[ii],t);
//         std::vector<int>li1all=getsamerange(li1all1,li1all2);
//         std::vector<int>lallchange=lall;
//         double sum=0;
//         for (int li1:li1all)
//         {
//             lallchange.push_back(li1);
//             sum=sum+qt2term( rjl, rjlinv,ystrall,ystrallinv,
//                 i, ii+1, k, rkp,lallchange, rip,qstrall,torder, tnum);
//             lallchange.pop_back();
//         }
//         return sum;
//     }
//
// }
//
// double qtqttest(basis rjl, basis rjlinvcha,
//             std::vector<std::vector<std::vector<double> > > ystrall,
//             std::vector<std::vector<std::vector<double> > > ystrallinv,
//             std::vector<std::vector<std::vector<double> > > qstrall, int tnum,
//             std::vector<int> qorder)
// {
//     if (rjlinvcha.sj.back()!=rjl.sj.back())
//     {
//         return 0;
//     }
//     basis rjlinv=rjlinvcha;
//     int N = rjl.r.size();
//     int jinv = rjlinv.j;
//     int jinvn = rjlinv.jn;
//     int t = qorder[tnum];
//     double sum = 0.0;
//     int sigt=sign_func({t});
//     if (!outfile.is_open())
//     {
//         // 如果文件未打开，则尝试打开文件
//         outfile.open("basis_output.txt",std::ios::app);
//     }
//     outfile<<"term1_______________________________"<<std::endl;
//     for (int k = 0; k <= N; k++)
//     {
//         if (!outfile.is_open())
//         {
//             // 如果文件未打开，则尝试打开文件
//             outfile.open("basis_output.txt",std::ios::app);
//         }
//         outfile<<"k="<<k<<std::endl;
//         std::vector<int> rkp;
//         if (k == 0)
//         {
//             if (rjl.j==0)
//             {
//                 continue;
//             }
//             rkp = genmvec(t, jinv);
//             for (int rkpp: rkp)
//             {
//                 int contj = -1;
//                 std::vector<int> rkp2 = is_get_in(allnucleus, 2, rkpp);
//                 int par = nucleus[rjlinv.jn].parity;
//                 if (isOdd(t / 2)) { par = par * (-1); } else { par = par * (1); }
//                 if (rkp2.empty())
//                 {
//                     continue;
//                 } else
//                 {
//                     for (const auto &elementnum: rkp2)
//                     {
//                         if (nucleus[elementnum].parity == par)
//                         {
//                             contj = elementnum;
//                         }
//                     }
//                     if (contj == -1)
//                     {
//                         continue;
//                     }
//                 }
//                 double c = 0;
//                 double c1 = 0;
//                 std::vector<int> sig;
//                 sig = {rkpp, -jinv, -t};
//                 c = sign_func(sig);
//                 c = c * mysqrt(rkpp + 1) / mysqrt(jinv + 1);
//                 c1 = (t + 1) / (mysqrt(jinv + 1) * mysqrt(rkpp + 1));
//                 c1 = c1 * qstrall[tnum][contj][jinvn] * qstrall[tnum][jinvn][contj];
//                 sum = sum + c1 * c * overlap(rjl, rjlinv, ystrall, ystrallinv);
//
//                 if (!outfile.is_open())
//                 {
//                     // 如果文件未打开，则尝试打开文件
//                     outfile.open("basis_output.txt",std::ios::app);
//                 }
//                 outfile<<"rkp="<<rkpp<<std::endl;
//                 outfile<<"basis1"<<std::endl;
//                 printBasisOneLinefile(outfile,rjl);
//                 outfile<<"basis2"<<std::endl;
//                 printBasisOneLinefile(outfile,rjlinv);
//                 double ov=overlap(rjl, rjlinv, ystrall, ystrallinv);
//                 int sizey=ystrall.size();
//                 for (int iy=0; iy<sizey; iy++)
//                 {
//                     if (!outfile.is_open())
//                     {
//                         // 如果文件未打开，则尝试打开文件
//                         outfile.open("basis_output.txt",std::ios::app);
//                     }
//                     outfile<<"ystrall1"<<"_"<<iy<<std::endl;
//                     writeMatrixToFile(ystrall[iy],outfile);
//                     if (!outfile.is_open())
//                     {
//                         // 如果文件未打开，则尝试打开文件
//                         outfile.open("basis_output.txt",std::ios::app);
//                     }
//                     outfile<<"ystrall2"<<"_"<<iy<<std::endl;
//                     writeMatrixToFile(ystrallinv[iy],outfile);
//
//                 }
//                 outfile<< "overlap="<<ov<<std::endl;
//
//             }
//         } else
//         {
//             int rk = rjlinv.r[k - 1];
//             rkp = genmvec(t, rjlinv.r[k - 1]);
//             for (int rkpp: rkp)
//             {
//                 std::vector<std::vector<double> > ystrchange;
//                 ystrchange = ychangeb1(ystrallinv, qstrall, k, rkpp, tnum, qorder, rjlinv);
//                 std::vector<std::vector<std::vector<double> > > ystrallch = ystrallinv;
//                 ystrallch[k - 1] = ystrchange;
//                 double c = 0;
//                 std::vector<int> sig;
//                 sig = {rkpp, -rk};
//                 c = sign_func(sig);
//                 c = c * mysqrt(rkpp + 1) / mysqrt(rk + 1);
//                 sum = sum + c * overlap(rjl, rjlinv, ystrall, ystrallch);
//                 if (!outfile.is_open())
//                 {
//                     // 如果文件未打开，则尝试打开文件
//                     outfile.open("basis_output.txt",std::ios::app);
//                 }
//                 outfile<<"rkp="<<rkpp<<std::endl;
//                 outfile<<"basis1"<<std::endl;
//                 printBasisOneLinefile(outfile,rjl);
//                 outfile<<"basis2"<<std::endl;
//                 printBasisOneLinefile(outfile,rjlinv);
//                 double ov=overlap(rjl, rjlinv, ystrall, ystrallch);
//                 int sizey=ystrall.size();
//                 for (int iy=0; iy<sizey; iy++)
//                 {
//                     if (!outfile.is_open())
//                     {
//                         // 如果文件未打开，则尝试打开文件
//                         outfile.open("basis_output.txt",std::ios::app);
//                     }
//                     outfile<<"ystrall1"<<"_"<<iy<<std::endl;
//                     writeMatrixToFile(ystrall[iy],outfile);
//                     if (!outfile.is_open())
//                     {
//                         // 如果文件未打开，则尝试打开文件
//                         outfile.open("basis_output.txt",std::ios::app);
//                     }
//                     outfile<<"ystrall2"<<"_"<<iy<<std::endl;
//                     writeMatrixToFile(ystrallch[iy],outfile);
//
//                 }
//                 if (!outfile.is_open())
//                 {
//                     // 如果文件未打开，则尝试打开文件
//                     outfile.open("basis_output.txt",std::ios::app);
//                 }
//                 outfile<< "overlap="<<ov<<std::endl;
//                 outfile<< "c="<<c<<std::endl;
//             }
//         }
//     }
//     if (!outfile.is_open())
//     {
//         // 如果文件未打开，则尝试打开文件
//         outfile.open("basis_output.txt",std::ios::app);
//     }
//     outfile<<"term2"<<"_______________________"<<std::endl;
//     for (int k = 1; k <= N; k++)
//     {
//         if (!outfile.is_open())
//         {
//             // 如果文件未打开，则尝试打开文件
//             outfile.open("basis_output.txt",std::ios::app);
//         }
//         outfile<<"k"<<"="<<k<<std::endl;
//         for (int i = 0; i <= k-1; i++)
//         {
//             if (!outfile.is_open())
//             {
//                 // 如果文件未打开，则尝试打开文件
//                 outfile.open("basis_output.txt",std::ios::app);
//             }
//             outfile<<"i"<<"="<<i<<std::endl;
//             if (i==0)
//         {
//             if (rjlinv.j==0)
//             {
//                 continue;
//             }
//         }
//             basis rjlinvp=rjlinvcha;
//             std::vector<int> r=rjlinvcha.r;
//             std::vector<int> rip;
//             std::vector<int> rkp;
//             double c1=1;
//             double c2=1;
//             rkp=genmvec(r[k-1], t);
//             std::vector<std::vector<std::vector<double> > > ystrallinvchan=ystrallinv;
//             if (i==0)
//             {
//                 rip=genmvec(rjlinvcha.j,t);
//
//                 for (int ripp: rip)
//                 {
//                     int contj = -1;
//                     std::vector<int> rip2= is_get_in(allnucleus, 2, ripp);
//                     int par = nucleus[rjlinvcha.jn].parity;
//                     if (isOdd(t / 2))
//                     {
//                         par = par * (-1);
//                     } else
//                     {
//                         par = par * (1);
//                     }
//                     if (rip2.empty())
//                     {
//                         continue;
//                     } else
//                     {
//                         for (const auto &elementnum: rip2)
//                         {
//                             if (nucleus[elementnum].parity == par)
//                             {
//                                 contj = elementnum;
//                             }
//                         }
//                         if (contj == -1)
//                         {
//                             continue;
//                         }
//                     }
//
//                     rjlinv.jn=contj;
//                     c1=-qstrall[tnum][contj][rjlinv.jn];
//                     c2=mysqrt(t + 1) / mysqrt(ripp + 1);
//                     sum=sum+c1*c2*qt2term(rjl,rjlinv,ystrall,ystrallinv,
//                             i,i,k,rkp,{ripp},ripp,qstrall,qorder,tnum);
//                 }
//             }else
//             {
//
//                 rip=genmvec(r[i-1],t);
//                 for (int ripp: rip)
//                 {
//                     if (!outfile.is_open())
//                     {
//                         // 如果文件未打开，则尝试打开文件
//                         outfile.open("basis_output.txt",std::ios::app);
//                     }
//                     outfile<<"rip="<<ripp<<std::endl;
//                     std::vector<std::vector<double>> ystrchangpp=
//                         ychangep1(ystrallinv, qstrall, i, ripp, tnum, qorder, rjlinvcha);
//                     ystrallinvchan[i-1]=ystrchangpp;
//                     std::vector<int>lii1;
//                     lii1=genmvec(rjlinv.sj[i-1],t);
//                     std::vector<int>lii2;
//                     lii2=lii1;
//                     if (i!=1)
//                     {
//                         lii2=genmvec(rjlinv.sj[i-2],ripp);
//                     }
//                     else
//                     {
//                         lii2=genmvec(rjlinv.j,ripp);
//                     }
//                     std::vector<int>lii=getsamerange(lii1, lii2);
//                     for (int lii1_1: lii)
//                     {
//                         sum=sum+(sigt*qt2term1(rjl,rjlinv,ystrall,ystrallinvchan,
//                             i,i,k,rkp,{lii1_1},ripp,qstrall,qorder,tnum));
//                     }
//                 }
//             }
//         }
//     }
//     return sum;
// }
void test_calculation(int num_calculations) {
    // 记录开始时间
    auto start_time = std::chrono::high_resolution_clock::now();

    // 创建线程容器
    std::vector<std::thread> threads;

    // 启动10个线程，每个线程计算一个6j系数，并执行num_calculations次
    for (int i = 2; i < 10; ++i) {  // 总共十个不同的6j系数
        threads.emplace_back([i, num_calculations] {
            for (int j = 0; j < num_calculations; ++j) {  // 每个系数计算num_calculations次
                double result = util::wigner_6j(i, i, i, i, i, i);
                // 可选：可以输出结果，但输出太多会影响时间
                // std::cout << "Result for 6j(" << i << ", " << i+1 << ", " << i+2 << ", "
                //           << i+3 << ", " << i+4 << ", " << i+5 << ") = " << result << std::endl;
            }
        });
    }

    // 等待所有线程完成
    for (auto& t : threads) {
        t.join();
    }

    // 记录结束时间
    auto end_time = std::chrono::high_resolution_clock::now();

    // 计算并输出计算时间
    std::chrono::duration<double> duration = end_time - start_time;
    std::cout << "Total computation time for " << num_calculations << " calculations: "
              << duration.count() << " seconds" << std::endl;
}


int main()
{
    // std::vector<std::vector<double> > ystrget;
    // nucleus={
    //     {0, 8, 7},
    //     {2, 4, 5},
    //     {2, 4, 3},
    //     {4, 0, 1},
    //     {0,10,11}
    // };
    // allnucleus={
    //     {0, 8, 7},
    //     {2, 4, 5},
    //     {2, 4, 3},
    //     {4, 0, 1},
    //     {0,10,11}
    // };
    // int size1=nucleus.size();
    // ystrget = {
    //     {4, 4, 0, -0.05916894},
    //     {0, 0, 0,  0.53584976},
    //     {1, 1, 0,  0.83771066},
    //     {2, 2, 0,  0.06838028},
    //     {3, 3, 0,  0.05412055},
    //     {4, 4, 1, -0.06582480},
    //     {0, 0, 1, -0.78392267},
    //     {0, 1, 1,  0.16358637},
    //     {0, 2, 1, -0.13229444},
    //     {1, 0, 1, -0.16358637},
    //     {1, 1, 1, -0.49316989},
    //     {1, 2, 1, -0.06876684},
    //     {1, 3, 1, -0.12267698},
    //     {2, 0, 1, -0.13229444},
    //     {2, 1, 1,  0.06876684},
    //     {2, 2, 1, -0.05965486},
    //     {2, 3, 1,  0.05327818},
    //     {3, 1, 1, -0.12267698},
    //     {3, 2, 1, -0.05327818}
    // };
    // basis rjl;
    // //rjl={0,0,{4,4,4},{1,1,1},{2,6,10}};
    // rjl={0,0,{0,0,0},{0,0,0},{0,0,0}};
    // //rjl={5,0,{4,4},{1,1},{9,13}};
    // std::vector<int> ro1=rjl.rn;
    // basis rjlinv;
    // rjlinv={0,0,{0,0,0},{0,0,0},{0,0,0}};
    //
    // //rjlinv={0,0,{4,4,4},{1,1,1},{2,6,10}};
    // //rjlinv={5,0,{4,4},{1,1},{9,13}};
    // std::vector<int> ro2=rjlinv.rn;
    // std::vector<int> rorder={0,4};
    // int size2=3;
    //
    // std::vector<std::vector<std::vector<double>>> ystr1(size2, std::vector<std::vector<double>>(size1, std::vector<double>(size1, 0)));
    // std::vector<std::vector<std::vector<double>>> ystr2(size2, std::vector<std::vector<double>>(size1, std::vector<double>(size1, 0)));
    //
    // getystrbasis(ystr1, ystrget, ro1);
    // getystrbasis(ystr2, ystrget, ro2);
    // std::vector<std::vector<double> > qstrgetp;
    // std::vector<int> qorderp={4};
    // qstrgetp=qstrall(qorderp);
    // int sizeq=qorderp.size();
    // std::vector<std::vector<std::vector<double>>> qstrallp(sizeq, std::vector<std::vector<double>>(size1,
    // std::vector<double>(size1, 0)));
    // getystrbasis(qstrallp, qstrgetp, {0});
    // int tnum=0;
    // double n=qtqttest(rjl,rjlinv,ystr1,ystr2,qstrallp,tnum,qorderp);
    // double ninv=qtqttest(rjlinv,rjl,ystr2,ystr1,qstrallp,tnum,qorderp);
    // double mm=1;
    // for (int num_calculations : {10, 100, 1000, 10000,100000,1000000,10000000,100000000,1000000000}) {
    //     test_calculation(num_calculations);
    // }




    return 0;



}
