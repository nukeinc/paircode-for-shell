//
// Created by wang- on 24-12-30.
//

#include "../inc/moe.h"

#include "../inc/moe.h"
#include "../inc/writeout.h"

#include"../inc/mtcal.h"
#include <omp.h> // OpenMP 并行加速
double overlap(basis x, basis y, std::vector<std::vector<std::vector<double> > > ystrall,
               std::vector<std::vector<std::vector<double> > > ystrallystrallinv);

double getlis(basis rjl, basis rjlinv, int k, int i, int ii, int t, std::vector<int> lii, int lk_1, int rip,
              std::vector<int> lkall,
              std::vector<std::vector<std::vector<double> > > ystrall,
              std::vector<std::vector<std::vector<double> > > ystrallinv, int contj)
{
    ////std::cout << "ii= " << ii << std::endl;
    if (i == 0 && rjl.j == 0) { return 0; }
    int rk;
    rk = rjl.r[k - 1];
    int ri;
    if (i == 0) { ri = rjl.j; } else { ri = rjl.r[i - 1]; }
    std::vector<int> r = rjl.r;
    std::vector<int> sj = rjl.sj;
    std::vector<int> sjinv = rjlinv.sj;
    std::vector<int> rinv = rjlinv.r;
    int ji;
    if (ii==0)
    {
        ji=rjl.j;
    }else
    {
        ji=rjl.sj[ii - 1];
    }
    int N1 = rinv.size();
    int N2 = r.size();
    std::vector<int> li1;
    std::vector<int> li2;
    std::vector<int> li;
    double sum = 0;
    std::vector<int> liitem = lii;
    li1 = genmvec(t, ji);
    li2 = li1;
    /*
    if (ii != i)
    {
        li2 = genmvec(lii.back(), rip);
    }else if (ii==0)
    {
        li2={rip};
    }else
    {
        int ji_1;
        if (i==1)
        {
            ji_1=rjl.j;
        }else
        {
            ji_1=rjl.sj[i-2];
        }
        li2=genmvec(rip,ji_1);
    }
    */
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
        if (is_in_array(li, lk_1))
        {

            double g = gfrom6_j(rjl, k, rinv[N1 - 1], t, lk_1);
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
            } else if (i == k - 2)
            {
                q = qfrom6_j(rjl, k - 1, t, lk_1, lii.back());
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
            // if (!outfile.is_open())
            // {
            //     // 如果文件未打开，则尝试打开文件
            //     outfile.open("basis_output.txt",std::ios::app);
            // }
            // outfile<< "g="<<g<<std::endl;
            // outfile<< "m="<<m<<std::endl;
            // outfile<<"q="<<q<<std::endl;
            // std::cout << "gmq" << std::endl;
            if (g==0)
            {
                double ggg=1;
            }
            // std::cout << g << std::endl;
            // std::cout << m << std::endl;
            // std::cout << q << std::endl;
            // std::cout << "Lk_1="<<lk_1 << std::endl;
            if (i == 0)
            {
                double c2;
                c2 = std::pow(-1, (t - rjl.j - rip)/2.0) * 4 * mysqrt(rk+1 ) * mysqrt(rinv[N2 - 1]+1) / mysqrt(rip+1);
                int sn = rinv[N2 - 1];
                int j = rjl.j;

                basis rjl1;
                basis rjl2;
                rjl1.j = rip;
                rjl2.j = rjlinv.j;
                rjl1.r = rjl.r;
                rjl2.r = rjlinv.r;
                rjl2.r.erase(rjl2.r.end()-1);
                rjl2.sj = rjlinv.sj;
                rjl2.sj.erase(rjl2.sj.end()-1);
                lii.push_back(lk_1);
                rjl1.sj = lii;
                rjl1.jn = contj;
                rjl2.jn = rjlinv.jn;

                std::vector lkall1 = lkall;

                reverse(lkall1.begin(), lkall1.end());
                rjl1.sj.insert(rjl1.sj.end(), lkall1.begin(), lkall1.end());
                rjl1.r.erase(rjl1.r.begin() + k - 1);
                if (i==0)
                {
                    rjl1.sj.erase(rjl1.sj.begin());
                }
                std::vector<std::vector<std::vector<double> > > ystrall1 = ystrall;

                std::vector<std::vector<std::vector<double> > > ystrall2 = ystrallinv;
                removestr(ystrall1,k-1);
                ystrall2.erase(ystrall2.end());
                double ov=overlap(rjl1, rjl2, ystrall1, ystrall2);
                //std::cout << "overlapmid"<<ov << std::endl;
                return c * overlap(rjl1, rjl2, ystrall1, ystrall2);

            } else
            {
                std::vector<std::vector<double> > ystrallchange = getz(ystrall, ystrallinv,
                                                                       i, k, rip, rk, ri, rjl, N2, rjlinv.r[N2 - 1], t);
                std::vector<std::vector<std::vector<double> > > ystrall1 = ystrall;
                std::vector<std::vector<std::vector<double> > > ystrall2 = ystrallinv;
                ystrall1[i - 1] = ystrallchange;
                basis rjl1;
                basis rjl2;
                rjl1.j = rjl.j;
                rjl2.j = rjlinv.j;
                rjl1.r = rjl.r;
                rjl2.r = rjlinv.r;
                rjl2.r.erase(rjl2.r.end());
                rjl2.sj = rjlinv.sj;
                rjl2.sj.erase(rjl2.sj.end());
                rjl1.sj = rjl.sj;
                rjl1.sj.erase(rjl1.sj.begin() + i-1, rjl1.sj.end());
                rjl1.jn = rjl.jn;
                rjl2.jn = rjlinv.jn;
                std::vector lkall1 = lkall;

                reverse(lkall1.begin(), lkall1.end());
                if(k==rjl.sj.size())
                {
                    lkall1.pop_back();
                }
                lii.push_back(lk_1);
                rjl1.sj.insert(rjl1.sj.end(), lii.begin(), lii.end());
                rjl1.sj.insert(rjl1.sj.end(), lkall1.begin(), lkall1.end());
                // rjl1.sj.erase(rjl1.sj.begin() + k - 1);
                rjl1.r.erase(rjl1.r.begin() + k - 1);
                rjl1.r[i-1]=rip;



                std::swap(ystrall1[k - 1], ystrall1.back());
                std::vector<std::vector<std::vector<double>>> ystrall1_1(ystrall1.begin(), ystrall1.end() - 1);
                ystrall2.erase(ystrall2.end()-1);
                //double ov=overlap(rjl1, rjl2, ystrall1_1, ystrall2);
                //std::cout << "overlapmid"<<ov << std::endl;
                return c * overlap(rjl1, rjl2, ystrall1_1, ystrall2);

            }
        } else
        {
            return 0;
        }
    }
    for (int li_1: li)
    {
        //std::cout << "li_1= " << li_1 << std::endl;
        std::vector<int> lis;
        lis = lii;
        lis.push_back(li_1);
        sum = sum + getlis(rjl, rjlinv, k, i, ii+1, t, lis, lk_1, rip, lkall, ystrall, ystrallinv, contj);
        lis.pop_back();
    }
    //std::cout << "getlis"<<sum << std::endl;
    return sum;
};


double getLk(basis rjl, basis rjlinv, int sn, int k, int m, std::vector<int> lall, double h,
             std::vector<std::vector<std::vector<double> > > ystrall,
             std::vector<std::vector<std::vector<double> > > ystrallystrallinv)
{
    //std::cout << "m= " << m<< "lall" << lall[0] << std::endl;
    int l = lall.back();
    std::vector<std::vector<int> > result;
    std::vector<int> r = rjl.r;
    std::vector<int> sj = rjl.sj;
    std::vector<int> sjinv = rjlinv.sj;
    std::vector<int> rinv = rjlinv.r;
    int N2 = r.size();
    std::vector<int> Leg1;
    std::vector<int> Leg2;
    std::vector<int> Leg;
    double sum = 0;
    std::vector<std::vector<std::vector<double> > > ystrall1 = ystrall;
    std::vector<std::vector<std::vector<double> > > ystrall2 = ystrallystrallinv;
    //std::cout << m << std::endl;
    if (m == k-1 || k==N2)
    {

        Leg1 = genmvec(sj[m + 1 - 1], sn);
        Leg2 = Leg1;
        if (m <=N2 - 2)
        {
            Leg2 = genmvec(l, r[m + 2 - 1]);
        }
        Leg = getsamerange(Leg1, Leg2);
        if (Leg.empty())
        {
            return 0;
        }
        if (k==1 && rjl.j==0)
        {
            Leg={0};
        }
        if (k==N2)
        {
            Leg=lall;
        }
        double sum1=0;
        for (int element: Leg)
        {   //std::cout << "leg" <<element<<  std::endl;
            ystrall1 = ystrall;
            ystrall2 = ystrallystrallinv;
            std::vector<int> x = lall;
            if (k!=N2)
            {
                x.push_back(element);
            }
            double h2 = h;
            auto it = x.rbegin();
            double phi = 1;
            if (k == N2 - 1)
            {
                h2 = hfrom6_j(rjl, N2, sn, x[0], x[1]);

            } else if (k == N2)
            {
                h2 = 1;
            } else
            {
                h2 = h * hfrom6_j(rjl, m + 2, sn, *(++it), *it);
            }
            double phik = getphi(ystrall, ystrallystrallinv, rjl, k);
            // std::cout <<"phik="<< phik << std::endl;



            double term1;
            basis rjl1;
            basis rjl2;
            std::vector<int> r1_1;
            std::vector<int> sj1_1;
            std::vector<int> r2_1;
            std::vector<int> sj2_1;
            r2_1 = rjlinv.r;
            sj2_1 = rjlinv.sj;
            r2_1.erase(r2_1.end() - 1);
            sj2_1.erase(sj2_1.end() - 1);
            r1_1 = r;
            r1_1.erase(r1_1.begin() + (k - 1));
            sj1_1 = sj;
            sj1_1.erase(sj1_1.begin() + (k - 1), sj1_1.end());
            std::vector<int>lkall1=lall;
            reverse(lkall1.begin(), lkall1.end());
            if (k!=N2)
            {
                sj1_1.insert(sj1_1.end(), lkall1.begin(), lkall1.end());
            }
            removestr(ystrall1, k - 1);
            removestr(ystrall2, N2 - 1);
            rjl1.j = rjl.j;
            rjl1.r = r1_1;
            rjl1.sj = sj1_1;
            rjl2.j = rjlinv.j;
            rjl2.r = r2_1;
            rjl2.sj = sj2_1;
            rjl1.jn = rjl.jn;
            rjl2.jn = rjlinv.jn;

            int jk_1 = rjl1.j;
            if (k != 1)
            {jk_1 = sj[k - 2];
            }

            term1 = phik * deltatwo(rinv[N2 - 1], r[k - 1]) * deltatwo(element, jk_1) * overlap(
                        rjl1, rjl2, ystrall1, ystrall2);

            // if (!outfile.is_open())
            // {
            //     // 如果文件未打开，则尝试打开文件
            //     outfile.open("basis_output.txt",std::ios::app);
            // }
            // outfile<<"phik=" <<phik <<std::endl;
            // outfile<<"term1=" <<term1 <<std::endl;

            double term2 = 0;

            std::vector<int> t;
            t = genmvec(r[k - 1], rinv[N2 - 1]);


            for (int i = k - 1; i >= 0; i = i - 1)
            {
                // if (!outfile.is_open())
                // {
                //     // 如果文件未打开，则尝试打开文件
                //     outfile.open("basis_output.txt",std::ios::app);
                // }
                // outfile<<"i=" <<i <<std::endl;

                // std::cout << "i="<< i << std::endl;
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
                    rip = genmrivec(t, r[i - 1]);
                }
                for (int ripp: rip)
                {
                    // std::cout << "ripp=" << ripp << std::endl;
                    // if (!outfile.is_open())
                    // {
                    //     // 如果文件未打开，则尝试打开文件
                    //     outfile.open("basis_output.txt",std::ios::app);
                    // }
                    // outfile<<"ripp=" << ripp<<std::endl;


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
                        // outfile<<"t=" << tri<<std::endl;



                        double ji_1=rjl.j;
                        if (i>1)
                        {
                            ji_1 = sj[i-2];
                        }else if (i==0)
                        {
                            ji_1 = 0;
                        }

                        // std::cout << "t= " << tri << std::endl;
                        double sj0;
                        if (i==0)
                        {
                            sj0=rjl.j;
                        }else
                        {
                            sj0=sj[i-1];
                        }
                        std::vector<int> li1 = genmvec(tri, sj0);

                        std::vector<int> li2 = genmvec(ripp, ji_1);
                        li1 = getsamerange(li1, li2);
                        /*
                        for (int lii2: li1)
                        {
                            //std::cout << "lii= " << lii2 << std::endl;
                            std::vector<int> lii;
                            lii = {lii2};
                            int ii = i;
                            double c1 = 1;
                            double c2sum = 1;
                            if (i == 0)
                            {
                                if (rjl.j == 0)
                                {
                                    continue;
                                }
                                c1 = std::pow(-1, (tri - rjl.j - ripp)/2.0) * 4 * mysqrt(rjl.r[k - 1] + 1) * mysqrt(
                                         rinv[N2 - 1]+1) / mysqrt(ripp + 1);
                                c2sum = 0;
                                for (int bcont = 0; bcont < ystrall[k - 1][contj].size(); ++bcont)
                                {   //std::cout << "bcont= " << bcont << std::endl;
                                    double c2 = 1;
                                    c2 = ystrall[k - 1][contj][bcont];
                                    c2 = c2 * ystrallystrallinv[N2 - 1][bcont][rjlinv.jn];
                                    c2 = c2 * util::wigner_6j(rjl.r[k - 1], rjlinv.r[N2 - 1], tri, rjl.j, ripp,
                                                              nucleus[bcont].j);
                                    c2sum = c2sum + c2;
                                }
                            }
                            if (c2sum == 0)
                            {
                                continue;
                            }
*/
                            int ii = i;
                            double c1 = 1;
                            double c2sum = 1;
                            if (i == 0)
                            {
                                if (rjl.j == 0)
                                {
                                    continue;
                                }
                                double sig=std::pow(-1, (tri - rjl.j - ripp)/2.0);
                                c1 = std::pow(-1, (tri - rjl.j - ripp)/2.0) * 4 * mysqrt(rjl.r[k - 1] + 1) * mysqrt(
                                         rinv[N2 - 1]+1)*mysqrt(tri+1) / mysqrt(ripp + 1);
                                c2sum = 0;
                                for (int bcont = 0; bcont < ystrall[k - 1][contj].size(); ++bcont)
                                {
                                    double c2 = 1;
                                    c2 = ystrall[k - 1][contj][bcont];
                                    c2 = c2 * ystrallystrallinv[N2 - 1][bcont][rjl.jn];
                                    c2 = c2 * util::wigner_6j(rjl.r[k - 1], rjlinv.r[N2 - 1], tri, rjl.j, ripp,
                                                              nucleus[bcont].j);
                                    c2sum = c2sum + c2;
                                }
                            }
                            if (c2sum == 0)
                            {
                                continue;
                            }

                            std::vector<int> lii;
                        double glisterm=getlis(rjl, rjlinv, k, i, ii, tri, lii, element, ripp, lall,
                                                                    ystrall,
                                                                    ystrallystrallinv, contj);

                        term2 = term2 + c1 * c2sum *glisterm;
                        // if (!outfile.is_open())
                        // {
                        //     // 如果文件未打开，则尝试打开文件
                        //     outfile.open("basis_output.txt",std::ios::app);
                        // }
                        // outfile<<"glisterm=" << glisterm<<std::endl;
                        // outfile<<"term2=" << term2<<std::endl;

                        //std::cout << "getlisterm"<< glisterm << std::endl;
                        //std::cout << "c1"<< c1 << std::endl;
                        //std::cout << "c2"<< c2sum << std::endl;




                    }
                }
            }
            // std::cout <<"term2"<< term2 << std::endl;
            // std::cout <<"term1"<< term1 << std::endl;
            // std::cout <<"h2"<< h2 << std::endl;
            // if (!outfile.is_open())
            // {
            //     // 如果文件未打开，则尝试打开文件
            //     outfile.open("basis_output.txt",std::ios::app);
            // }
            sum1=sum1+((term2+term1)*h2);
            // outfile<<"h2=" << h2<<std::endl;

            x.pop_back();
        }
        return sum1;
    }
    Leg1 = genmvec(sj[m + 1 - 1], sn);
    Leg2 = Leg1;
    if (m <=N2 - 2)
    {
        Leg2 = genmvec(l, r[m + 2 - 1]);
    }
    Leg = getsamerange(Leg1, Leg2);
    if (Leg.empty())
    {
        return 0;
    }
    int m1 = m - 1;
    for (int element: Leg)
    {
        //std::cout << "Leg " << element << std::endl;
        std::vector<int> x = lall;
        x.push_back(element);
        double h1 = h;
        if (m <= N2 - 2)
        {
            auto it = x.rbegin();
            h1 = h * hfrom6_j(rjl, m + 2, sn, *(++it), *it);
        }
        sum = sum + getLk(rjl, rjlinv, sn, k, m1, x, h1, ystrall, ystrallystrallinv);


    }
    return sum;
};

// double pair2cal(
//     const std::vector<std::vector<std::vector<double>>>& ystra1_1,
//     const std::vector<std::vector<std::vector<double>>>& ystra2_1, basis rjl,basis rjlinv)
// {
//
//     // outfile<<"basis1"<<std::endl;
//     // printBasisOneLinefile(outfile,rjl);
//     // outfile<<"basis2"<<std::endl;
//     // printBasisOneLinefile(outfile,rjlinv);
//     // outfile<<"ystrall"<<std::endl;
//     // writeMatrixToFile(ystra1_1[0],outfile);
//     // if (!outfile.is_open())
//     // {
//     //     // 如果文件未打开，则尝试打开文件
//     //     outfile.open("basis_output.txt",std::ios::app);
//     // }
//     // outfile<<"ystrall1_2"<<std::endl;
//     // writeMatrixToFile(ystra1_1[1],outfile);
//     // if (!outfile.is_open())
//     // {
//     //     // 如果文件未打开，则尝试打开文件
//     //     outfile.open("basis_output.txt",std::ios::app);
//     // }
//     // outfile<<"ystrall2_1"<<std::endl;
//     // if (!outfile.is_open())
//     // {
//     //     // 如果文件未打开，则尝试打开文件
//     //     outfile.open("basis_output.txt",std::ios::app);
//     // }
//     // writeMatrixToFile(ystra2_1[0],outfile);
//     // if (!outfile.is_open())
//     // {
//     //     // 如果文件未打开，则尝试打开文件
//     //     outfile.open("basis_output.txt",std::ios::app);
//     // }
//     // outfile<<"ystrall2_2"<<std::endl;
//     // writeMatrixToFile(ystra2_1[1],outfile);
//     double part1 = 0;
//     double part2 = 0;
//     double part3 = 0;
//     int S1_DIM=ystra1_1.size();
//     int A_DIM=ystra1_1[0].size();
//     int B_DIM=ystra1_1[0][0].size();
//
//
//     // 计算公式的第一个部分
//     if (rjl.r[0]==rjlinv.r[0]&& rjl.r[1]==rjlinv.r[1])
//     {
//
//         for (int a = 0; a < A_DIM; ++a) {
//             for (int b = 0; b < B_DIM; ++b) {
//                 for (int c = 0; c < A_DIM; ++c) {
//                     for (int d = 0; d < B_DIM; ++d) {
//                         part1 += 4 * ystra1_1[0][a][b] * ystra2_1[0][a][b] * ystra1_1[1][c][d]* ystra2_1[1][c][d];  // 使用 ystra1_1, ystra2_1
//                         // 和 ystra2_2
//                     }
//                 }
//             }
//         }
//     }
//     // if (!outfile.is_open())
//     // {
//     //     // 如果文件未打开，则尝试打开文件
//     //     outfile.open("basis_output.txt",std::ios::app);
//     // }
//
//
//     // 计算公式的第二个部分
//     if (rjl.r[0]==rjlinv.r[1]&& rjl.r[1]==rjlinv.r[0])
//     {
//         int sig2=std::pow(-1,((rjl.r[0]+rjl.r[1]+rjl.sj[1])/2));
//         for (int a = 0; a < A_DIM; ++a) {
//             for (int b = 0; b < B_DIM; ++b) {
//                 for (int c = 0; c < A_DIM; ++c) {
//                     for (int d = 0; d < B_DIM; ++d) {
//                         part2 += 4*sig2 * ystra1_1[0][a][b] * ystra2_1[1][a][b] * ystra1_1[1][c][d]* ystra2_1[0][c][d];
//                          // 使用 ystra1_1, ystra2_1
//                         // 和 ystra2_2
//                     }
//                 }
//             }
//         }
//     }
//
//
//     // 计算公式的第三个部分，求和操作
//     double cc=sqrt(rjlinv.r[0]+1)*sqrt(rjlinv.r[1]+1)*sqrt(rjl.r[0]+1)*sqrt(rjl.r[1]+1);
//     for (int a = 0; a < A_DIM; ++a) {
//         for (int b = 0; b < B_DIM; ++b) {
//             for (int c = 0; c < A_DIM; ++c) {
//                 for (int d = 0; d < B_DIM; ++d) {
//                     int ja=nucleus[a].j;
//                     int jb=nucleus[b].j;
//                     int jc=nucleus[c].j;
//                     int jd=nucleus[d].j;
//                     part3 += -16 * ystra2_1[0][a][b] * ystra2_1[1][c][d] * ystra1_1[0][a][c]*
//                     ystra1_1[1][b][d]*cc*util::wigner_9j(ja,jb,rjlinv.r[0],jc,jd,rjlinv.r[1],rjl.r[0],rjl.r[1],rjl.sj[1]);
//                     // 使用 ystra1_1, ystra2_1 和 ystra2_2
//                 }
//             }
//         }
//     }
//
//     // 最终结果是三个部分的和
//     return part1 + part2 + part3;
// }
// double pair2cal(
//     const std::vector<std::vector<std::vector<double>>>& ystra1_1,
//     const std::vector<std::vector<std::vector<double>>>& ystra2_1,
//     basis rjl, basis rjlinv)
// {
//     double part1 = 0, part2 = 0, part3 = 0;
//     int A_DIM = ystra1_1[0].size();
//     int B_DIM = ystra1_1[0][0].size();
//
//     // 直接计算 sqrt()，让编译器优化
//     double cc = mysqrt((rjlinv.r[0] + 1) * (rjlinv.r[1] + 1) * (rjl.r[0] + 1) * (rjl.r[1] + 1));
//
//     // 直接引用 vector，减少不必要的索引访问
//     auto& y0 = ystra1_1[0];
//     auto& y1 = ystra1_1[1];
//     auto& y2_0 = ystra2_1[0];
//     auto& y2_1 = ystra2_1[1];
//
//     bool cond1 = (rjl.r[0] == rjlinv.r[0] && rjl.r[1] == rjlinv.r[1]);
//     bool cond2 = (rjl.r[0] == rjlinv.r[1] && rjl.r[1] == rjlinv.r[0]);
//
//     if (cond1) {
//         for (int a = 0; a < A_DIM; ++a) {
//             for (int b = 0; b < B_DIM; ++b) {
//                 double y0ab = y0[a][b], y2_0ab = y2_0[a][b];
//                 double sum1=0;
//                 for (int c = 0; c < A_DIM; ++c) {
//                     for (int d = 0; d < B_DIM; ++d) {
//                         sum1 +=   y1[c][d] * y2_1[c][d];
//                     }
//                 }
//                 sum1=y0ab * y2_0ab * sum1;
//                 part1+=sum1;
//             }
//         }
//         part1=4*part1;
//     }
//
//
//     if (cond2) {
//         int sig2 = ((rjl.r[0] + rjl.r[1] + rjl.sj[1]) / 2) % 2 ? -1 : 1;
//         double coeff2 = 4 * sig2;
//
//         for (int a = 0; a < A_DIM; ++a) {
//             for (int b = 0; b < B_DIM; ++b) {
//                 double y0ab = y0[a][b], y2_1ab = y2_1[a][b];
//                 double sum2=0;
//                 for (int c = 0; c < A_DIM; ++c) {
//                     for (int d = 0; d < B_DIM; ++d) {
//                         sum2 +=y1[c][d] * y2_0[c][d];
//                     }
//                 }
//                 sum2=  y0ab * y2_1ab * sum2;
//                 part2+=sum2;
//             }
//         }
//         part2=coeff2 *part2;
//     }
//
//
//     for (int a = 0; a < A_DIM; ++a) {
//         for (int b = 0; b < B_DIM; ++b) {
//             double y2_0ab = y2_0[a][b];
//             for (int c = 0; c < A_DIM; ++c) {
//                 for (int d = 0; d < B_DIM; ++d) {
//                     int ja = nucleus[a].j, jb = nucleus[b].j, jc = nucleus[c].j, jd = nucleus[d].j;
//
//                     // 直接调用 `wigner_9j()`，让编译器优化
//                     double wigner_value = util::wigner_9j(ja, jb, rjlinv.r[0], jc, jd, rjlinv.r[1], rjl.r[0], rjl.r[1], rjl.sj[1]);
//
//                     part3 += -16 * y2_0ab * y2_1[c][d] * y0[a][c] * y1[b][d] * cc * wigner_value;
//                 }
//             }
//         }
//     }
//
//     return part1 + part2 + part3;
// }
double pair2cal(
    const std::vector<std::vector<std::vector<double>>>& ystra1_1,
    const std::vector<std::vector<std::vector<double>>>& ystra2_1,
    basis rjl, basis rjlinv)
{
    double part1 = 0;
    double part2 = 0;
    double part3 = 0;
    int A_DIM = ystra1_1[0].size();
    int B_DIM = ystra1_1[0][0].size();

    // 直接计算 sqrt()，让编译器优化
    double cc = mysqrt((rjlinv.r[0] + 1) * (rjlinv.r[1] + 1) * (rjl.r[0] + 1) * (rjl.r[1] + 1));

    // 直接引用 vector，减少不必要的索引访问
    auto& y0 = ystra1_1[0];
    auto& y1 = ystra1_1[1];
    auto& y2_0 = ystra2_1[0];
    auto& y2_1 = ystra2_1[1];

    bool cond1 = (rjl.r[0] == rjlinv.r[0] && rjl.r[1] == rjlinv.r[1]);
    bool cond2 = (rjl.r[0] == rjlinv.r[1] && rjl.r[1] == rjlinv.r[0]);

    if (cond1) {
        for (int a = 0; a < A_DIM; ++a) {
            for (int b = 0; b < B_DIM; ++b) {
                double y0ab = y0[a][b], y2_0ab = y2_0[a][b];
                double sum1=0;
                for (int c = 0; c < A_DIM; ++c) {
                    for (int d = 0; d < B_DIM; ++d) {
                        sum1 +=   y1[c][d] * y2_1[c][d];
                    }
                }
                sum1=y0ab * y2_0ab * sum1;
                part1+=sum1;
            }
        }
        part1=4*part1;
    }


    if (cond2) {
        int sig2 = ((rjl.r[0] + rjl.r[1] + rjl.sj[1]) / 2) % 2 ? -1 : 1;
        double coeff2 = 4 * sig2;

        for (int a = 0; a < A_DIM; ++a) {
            for (int b = 0; b < B_DIM; ++b) {
                double y0ab = y0[a][b], y2_1ab = y2_1[a][b];
                double sum2=0;
                for (int c = 0; c < A_DIM; ++c) {
                    for (int d = 0; d < B_DIM; ++d) {
                        sum2 +=y1[c][d] * y2_0[c][d];
                    }
                }
                sum2=  y0ab * y2_1ab * sum2;
                part2+=sum2;
            }
        }
        part2=coeff2 *part2;
    }


    for (int a = 0; a < A_DIM; ++a) {
        for (int b = 0; b < B_DIM; ++b) {
            double y2_0ab = y2_0[a][b];
            for (int c = 0; c < A_DIM; ++c) {
                for (int d = 0; d < B_DIM; ++d) {
                    int ja = nucleus[a].j, jb = nucleus[b].j, jc = nucleus[c].j, jd = nucleus[d].j;

                    // 直接调用 `wigner_9j()`，让编译器优化
                    double wigner_value = util::wigner_9j(ja, jb, rjlinv.r[0], jc, jd, rjlinv.r[1], rjl.r[0], rjl.r[1], rjl.sj[1]);

                    part3 += -16 * y2_0ab * y2_1[c][d] * y0[a][c] * y1[b][d] * cc * wigner_value;
                }
            }
        }
    }

    return part1 + part2 + part3;
}

double threecal(
    const std::vector<std::vector<std::vector<double>>>& ystra1_1,
    const std::vector<std::vector<std::vector<double>>>& ystra2_1,
    basis rjl, basis rjlinv)
{
    double result=0;
    double result1=0;
    double result2=0;
    int jp=rjlinv.j;
    int jpn=rjlinv.jn;
    int j=rjl.j;
    int jn=rjl.jn;
    int aj=rjl.sj[0];
    int s=rjlinv.r[0];
    int r=rjl.r[0];
    int c1=2*deltatwo(s,r)*deltatwo(j,jp);
    int nunum=nucleus.size();
    if (c1!=0)
    {
        double sum1=0;
        for (int a = 0; a < nunum; ++a)
        {
            for (int b = 0; b < nunum; ++b)
            {
                double c2=ystra1_1[0][a][b]*ystra2_1[0][a][b];
                sum1=sum1+c2;
            }
        }
        result1=sum1*c1;
    }
    double cc1=4*mysqrt(r+1)*mysqrt(s+1);
    std::vector<int> leg1=genmvec(s, r);
    std::vector<int> leg2=genmvec(j, jp);
    std::vector<int> leg3=getsamerange(leg1,leg2);
    double sum2=0;
    for (int t: leg3)
    {
        double g=gfrom6_j(rjl,1,s,t,jp);
        std::vector<int>sigpv={t,-j,-jp};
        int sigp=sign_func(sigpv);
        double cc2=cc1*mysqrt(t+1);
        double sum2_1=0;
        for (int b = 0; b < nunum; ++b)
        {
            int jb=nucleus[b].j;
            double cc3=ystra1_1[0][jpn][b]*ystra2_1[0][b][jn];
            double cc4=util::wigner_6j(r,s,t,j,jp,jb)/mysqrt(jp+1);
            sum2_1+=cc3*cc4*sigp*g;
        }
        sum2+=sum2_1*cc2;
    }
    result2=sum2;
    std::vector<int> sigv={aj,-jp,s};
    int sigp=sign_func(sigv);
    double c=mysqrt(jp+1)/mysqrt(aj+1);
    result2=sigp*c*result2;
    result=result1+result2;
    return result;

}

double overlap(basis x, basis y, std::vector<std::vector<std::vector<double> > > ystrall,
               std::vector<std::vector<std::vector<double> > > ystrallystrallinv)
{
    std::vector<int> r1 = x.r;
    std::vector<int> sj1 = x.sj;
    std::vector<int> r2 = y.r;
    std::vector<int> sj2 = y.sj;
    int j1 = x.j;
    int N1 = r1.size();
    int N2 = r2.size();
    int j2=y.j;
    if ( (sj1.back()) != (sj2.back()))
    {
        return 0;
    }
    if (N1==2 && x.j==0)
    {
        double sum=pair2cal(ystrall, ystrallystrallinv,x,y);
        return sum;
    }
    if (N1 == 1)
    {
        double result1;
        double result2;
        double y1 = 1;
        double y2 = 1;
        result1 = 2 * deltatwo(r1[0], r2[0]) * deltatwo(j1, j2);
        int nucnum = nucleus.size();
        double sum = 0;
        //std::cout << "calculate N=1j"  << sj1[0] << std::endl;
        for (int i_m = 0; i_m < nucnum; i_m++)
        {
            for (int i_n = 0; i_n < nucnum; i_n++)
            {
                sum = sum + ystrall[0][i_m][i_n] * ystrallystrallinv[0][i_m][i_n];
            }
        }
        result1 = sum * result1;
        double sum1 = 0;
        if (x.j==0)
        {
            result2=0;
        }else
        {
            result2 = 4 * mysqrt(x.r[0] + 1) * mysqrt(y.r[0] + 1);
            for (int i_m = 0; i_m < nucnum; i_m++)
            {
                int j11 = nucleus[i_m].j;
                int j1n = x.jn;
                int j2n = y.jn;
                double cc=util::wigner_6j(
                           x.r[0], y.j, j11, y.r[0], x.j, x.sj[0]);
                sum1 = sum1 + ystrall[0][i_m][j2n] * ystrallystrallinv[0][i_m][j1n] *cc;
            }
            if (std::abs(sum1)<0.01)
            {
                double mmmm=1;
            }
            result2=result2 * sum1;
        }
        double result;
        result = result1 + result2 ;
        if (std::abs(result-2)<0.01)
        {
            double mmmm=1;
        }
        // std::cout << "N1re " << result <<" r"<<x.r[0]<<" r"<<y.r[0]<<"j"<<x.j<<"j "<<y.j<<"sj "<<x.sj[0]<<"sj "<<y.sj[0]<<std::endl;
        return result;
        // double result=threecal(ystrall, ystrallystrallinv,x,y);
        // return result;
    } else
    {
        double jn1 = mysqrt(sj2[N1 - 2]+1);
        double jn2 = mysqrt(sj1[N2 - 1]+1);
        std::vector<int> signvec = {sj1[N2 - 1], -sj2[N2 - 2], r2[N1 - 1]};
        double c1 = (static_cast<double>(jn1) / jn2) * sign_func(signvec);
        double sum = 0;
        std::vector<int> ln_1;

        ln_1 = genmvec(sj1[N1 - 1], r2[N2 - 1]);
        if (!is_in_array(ln_1, sj2[N2 - 2]))
        {
            return 0;
        }
        for (int k = 1; k <= N1; ++k)
        {
            // if (!outfile.is_open())
            // {
            //     // 如果文件未打开，则尝试打开文件
            //     outfile.open("basis_output.txt",std::ios::app);
            // }
            // outfile<<"k="<<k<<std::endl;
            double c2;
            double hall = 1;
            std::vector<int> lall={sj2[N2 - 2]};

            // std::cout << "k= " << k << std::endl;
            sum = sum + getLk(x, y, y.r[N2 - 1], k, N1 - 2, lall, hall,
                              ystrall,
                              ystrallystrallinv);
            //std::cout << "sum= " << sum << std::endl;
            // outfile<<"sum="<<sum<<std::endl;

        }
        // outfile<<"c1="<<c1<<std::endl;
        // outfile<<"result="<<c1*sum<<std::endl;
        return c1 * sum;

    }
    return 0.0;
}

double ceshi(double x1)
{
    return x1;
}


