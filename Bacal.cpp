//
// Created by wang- on 25-7-11.
//
#define EIGEN_NO_DEBUG
#define EIGEN_UNROLLING_LIMIT 0
// #define EIGEN_USE_MKL_ALL
// #define EIGEN_USE_BLAS
// #define EIGEN_USE_LAPACKE
#include <iostream>
#include <fstream>
#include <array>
#include <chrono>
#include <cmath>
#include <iomanip>
#include "mtcal.h"
#include "qcal.h"
#include <cstdlib>
#include <lambda_lanczos/lambda_lanczos.hpp>
#include <Eigen/Dense>
#include <vector>
#include"writeout.h"
#include"bcscal.h"
#include <sstream>
#include <string>
#include <map>
#include <utility>
#include <maincal.h>


int main()
{


    nucleus={
        {0, 8, 7},
        {2, 4, 5},
        {2, 4, 3},
        {4, 0, 1},
        {0,10,11}
    };
    int size1=nucleus.size();
    allnucleus={
        {0, 8, 7},
        {2, 4, 5},
        {2, 4, 3},
        {4, 0, 1},
        {0,10,11}
    };
    nucleus2=
    {
        {0, 10,  9},
     {2,  6,  7},
     {2,  6,  5},
     {4,  2,  3},
     {4,  2,  1},
     {0, 12, 13}
    };
    int size2=nucleus.size();
    allnucleus2={
        {0, 10,  9},
     {2,  6,  7},
     {2,  6,  5},
     {4,  2,  3},
     {4,  2,  1},
     {0, 12, 13}
    };
    std::vector<double> energyp={-0.19392272 ,    -0.36722270  ,    2.45624118 ,     2.71147885 ,2.01885944     };
    std::vector<double> energyn={1.35292300 ,    -1.08623249 ,     1.69663730   ,   0.47103342    ,  1.27303342,1.13244657      };
    std::vector<std::vector<double> > ystrgetp;
    ystrgetp=readstrfile("D:/paircalpro/mschemecode/strfile1.txt");
    // ystrgetp = {
    //     {4, 4, 0, -0.05916894},
    //     {0, 0, 0,  0.53584976},
    //     {1, 1, 0,  0.83771066},
    //     {2, 2, 0,  0.06838028},
    //     {3, 3, 0,  0.05412055},
    //     {4, 4, 1, -0.06584450},
    //     {0, 0, 1, -0.78415728},
    //     {0, 1, 1,  0.16363533},
    //     {0, 2, 1, -0.13233403},
    //     {1, 0, 1, -0.16363533},
    //     {1, 1, 1, -0.49331749},
    //     {1, 2, 1, -0.06878742},
    //     {1, 3, 1, -0.12271369},
    //     {2, 0, 1, -0.13233403},
    //     {2, 1, 1,  0.06878742},
    //     {2, 2, 1, -0.05967271},
    //     {2, 3, 1,  0.05329413},
    //     {3, 1, 1, -0.12271369},
    //     {3, 2, 1, -0.05329413}
    // };
    std::vector<std::vector<double> > ystrgetp2;
    ystrgetp2=readstrfile("D:/paircalpro/mschemecode/strfile2.txt");
    normalizeBySecondColumn(ystrgetp);
    normalizeBySecondColumn(ystrgetp2);
    // nucleus=nucleus2;
    // allnucleus=allnucleus2;
    // size1=nucleus.size();


    std::vector<int> rorderp={0,4,8,12,6};
    std::vector<std::vector<int>>rorderrp={rorderp,{0,99,0,1},{0,1,2,3,4}};
    rvecall={0,4,8,12,6};
    rparityvec={1,1,1,1,-1};
    std::vector<std::pair<std::vector<int>, std::vector<int>>> catchpair1;
    catchpair1=generateValidPairs(4,rorderrp,51);
    std::vector<basis> allbasisp;
    allbasisp=calculateBasisWithRange(4, catchpair1);
    printBasisVectorOneLine(allbasisp);
    int sizep=allbasisp.size();
    int sizepr=allbasisp[0].r.size();
    std::vector<std::vector<std::vector<std::vector<double>>>> ystrallp(
        sizep,  // 第一维：大小为 sizep
        std::vector<std::vector<std::vector<double>>>(
            sizepr,  // 第二维：大小为 sizepr
            std::vector<std::vector<double>>(
                size1,  // 第三维：大小为 size1
                std::vector<double>(size1, 0.0)  // 第四维：大小为 size1，元素初始化为 0.0
            )
        )
    );
    calculateystr(ystrallp,ystrgetp,allbasisp);
    std::vector<basism> mtest;
    std::vector<basism> mtest_1;
    mtest=calculateBasisform(4,rorderrp);
    mtest_1=calculateBasisform_1(4,rorderrp);

    writeBasismvecToFile(mtest);
    Eigen::MatrixXd m1;
    Eigen::MatrixXd m1_1;
    m1=ComputeTransformationMatrix( mtest,allbasisp);
    m1_1=ComputeTransformationMatrix( mtest_1,allbasisp);

    std::vector<std::vector<double> >m1qget=m1qstrall(1.0,3.9102);
    std::vector<std::vector<std::vector<double>>> m1strallp1(1, std::vector<std::vector<double>>(size1,
    std::vector<double>(size1, 0)));
    getystrbasis(m1strallp1, m1qget, {0});

    std::cout <<"m1" << std::endl;
    std::cout << m1 << std::endl;
    std::cout <<"m1_1" << std::endl;
    std::cout << m1_1 << std::endl;
    std::vector<std::vector<Eigen::MatrixXd>>ystrm1= computeystrm( mtest, ystrgetp);
    std::vector<std::vector<Eigen::MatrixXd>>ystrm1_1= computeystrm( mtest_1, ystrgetp);


    Eigen::MatrixXd overlapm=caloverlapmmat(mtest,mtest,ystrm1,ystrm1);
    std::cout <<"overlapm" << std::endl;
    std::cout << overlapm << std::endl;
    std::vector<std::vector<double> > overlapmch=overlapchange(overlapm,m1);
    printMatrix(overlapmch);
    updateMatrices(allbasisp, overlapmch,ystrallp);
    std::vector<std::vector<double>>schmitmatp;

    schmitmatp=gramSchmidtInOverlap(overlapmch,allbasisp,ystrallp);
    printBasisVectorOneLine(allbasisp);




    std::vector<std::vector<double> > overlapp;
    overlapp=overlapmatcal(allbasisp,ystrallp);
    printMatrix(overlapp);

    m1=ComputeTransformationMatrix( mtest,allbasisp);
    m1_1=ComputeTransformationMatrix( mtest_1,allbasisp);
    std::vector<std::vector<double> > overlapmch1=overlapchange(overlapm,m1);
    std::cout <<"overlapmch1" << std::endl;
    printMatrix(overlapmch1);


    Eigen::MatrixXd schmitmat1=convertToEigenMatrix(schmitmatp);
    Eigen::MatrixXd overlap1=convertToEigenMatrix(overlapp);
    Eigen::MatrixXd overlap1cha=schmitmat1 * overlap1 * schmitmat1.transpose();
    // std::cout << overlap1cha << std::endl;


    std::vector<std::vector<double> > V1val;
    std::vector<DataRow> V1val1;
    std::vector<std::map<int, Matrix4D>>buildVValue2= processFile("D:/paircalpro/mschemecode/efct1.txt");
    std::ifstream input_file1("D:/paircalpro/mschemecode/efct3.txt"); // 替换为你的文件名
    if (!input_file1.is_open()) {
        std::cerr << "无法打开文件！" << std::endl;
        return 1;
    }
    std::vector<std::map<int, Matrix4D>>buildVValue3= processFilepn("D:/paircalpro/mschemecode/efct3.txt");
    std::map<int,Matrix4D>vpnmat= getvpnval(buildVValue3,{1.0});
    std::vector<std::vector<Eigen::MatrixXd>> q_pi;
    std::vector<Eigen::MatrixXd> V_it;
    std::vector<std::vector<Eigen::MatrixXd>> q_nu;
    std::vector<double> strengthvec={1.0};
    getsvdresult( vpnmat,q_pi,  // 输出: q_pi[i](j_alpha, j_beta)
        V_it,               // 输出: V_it(i, t)
        q_nu,strengthvec);
    printDiagonals(V_it);
    input_file1.close();
    // std::cout<<"transfermat"<<std::endl;
    // std::cout<<m1<<std::endl;
    // std::cout<<"q_pi"<<std::endl;
    // std::cout<<q_pi[2][0]<<std::endl;
    // std::cout<<q_nu[2][0]<<std::endl;


    auto nonZeropos= findNonZeroDiagonal(V_it);
    std::map<int,std::map<int,Eigen::MatrixXd>>qmatpi= qmatjcal( q_pi, m1,V_it,allbasisp,mtest,ystrm1);
    std::map<int,std::map<int,Eigen::MatrixXd>>qmatpi_1= qmatjcal_1( q_pi, m1,m1_1,
        V_it,allbasisp,mtest,mtest_1,ystrm1,ystrm1_1);
    std::cout<<"matre = "<<qmatpi[3][0]<<std::endl;
    std::cout<<"matre_1 = "<<qmatpi_1[3][0]<<std::endl;

    for (auto qmatcal : qmatpi)
    {
        int t=qmatcal.first;
        for (auto qmatcal2 : qmatcal.second)
        {
            int inum=qmatcal2.first;
            for (int i=0;i<allbasisp.size();++i)
            {
                int j1=allbasisp[i].sj.back();
                for (int j=0;j<allbasisp.size();++j)
                {
                    int j2=allbasisp[j].sj.back();
                    int num=(j1+j2+(t*2))/2;
                    if ((num & 1) != 0)
                    {
                        qmatpi[t][inum](i,j)=qmatpi_1[t][inum](i,j);
                    }
                }
            }
        }
    }
    std::cout<<"matreget = "<<qmatpi[2][0]<<std::endl;
    std::cout<<"matreget3 = "<<qmatpi[3][0]<<std::endl;



    std::map<int,std::map<int,Eigen::MatrixXd>> qmatpicha;
    for (auto qmatcal : qmatpi)
    {
        int t=qmatcal.first;
        std::cout<<"t = "<<t<<std::endl;
        for (auto qmatcal2 : qmatcal.second)
        {
            int inum=qmatcal2.first;
            qmatpicha[t][inum]=schmitmat1*qmatcal2.second* schmitmat1.transpose();
        }
    }




    double alpha1=4.33195830385423;
    std::vector<std::vector<double> > qstrgetp1;
    std::vector<int> qorderp1={4,6};
    qstrgetp1=qstrall(qorderp1,alpha1);

    int sizeq1=qorderp1.size();
    std::vector<std::vector<std::vector<double>>> qstrallp1(sizeq1, std::vector<std::vector<double>>(size1,
    std::vector<double>(size1, 0)));

    getystrbasis(qstrallp1, qstrgetp1, {0,1});


    qstrallp1.push_back(m1strallp1[0]);
    qorderp1.push_back(2);
    std::vector<Eigen::MatrixXd> bemematcal={};
    for (int tnum=0;tnum<qstrallp1.size();++tnum)
    {
        Eigen::MatrixXd qti=convertToEigenMatrix(qstrallp1[tnum]);
        Eigen::MatrixXd resultmat=bemecal(qti,m1,allbasisp,mtest,ystrm1,qorderp1[tnum]);
        Eigen::MatrixXd resultmat_1=bemecal_1(qti,m1,m1_1,allbasisp,
            mtest,mtest_1,ystrm1,ystrm1_1,qorderp1[tnum]);
        int t=qorderp1[tnum];
        for (int i=0;i<allbasisp.size();++i)
        {
            int j1=allbasisp[i].sj.back();
            for (int j=0;j<allbasisp.size();++j)
            {
                int j2=allbasisp[j].sj.back();
                int num=(j1+j2+t)/2;
                if ((num & 1) != 0)
                {
                    resultmat(i,j)=resultmat_1(i,j);
                }
            }
        }
        resultmat=schmitmat1*resultmat* schmitmat1.transpose();
        bemematcal.push_back(resultmat);
    }
    std::cout<<"bemematcal = "<<std::endl;
    for (const auto& bememat : bemematcal)
    {
        std::cout << bememat << std::endl;
    }














    std::vector<double>strength={1.0,0.66};
    Eigen::MatrixXd ham1cal=ham1( mtest,ystrm1,buildVValue2,strength);
    std::cout << "ham1cal" << std::endl;
    std::cout << ham1cal << std::endl;
    Eigen::MatrixXd changeham1=hamchange(mtest,allbasisp,m1, 0, 0, ham1cal);
    std::cout << "hamchangecal" << std::endl;
    std::cout << changeham1 << std::endl;


    Eigen::MatrixXd hamsige1= hamsige(mtest,ystrm1,energyp);
    std::cout << "hamsige" << std::endl;
    std::cout << hamsige1 << std::endl;
    Eigen::MatrixXd hamchangesig1=hamchange(mtest,allbasisp,m1, 0, 0, hamsige1);
    std::cout << "hamsigechange" << std::endl;
    std::cout << hamchangesig1 << std::endl;
    Eigen::MatrixXd hamp=changeham1+ hamchangesig1;
    std::cout << "hamp" << std::endl;
    std::cout << hamp << std::endl;
    Eigen::MatrixXd hampchange=schmitmat1*hamp* schmitmat1.transpose();
    std::cout << "hampchange" << std::endl;
    std::cout << hampchange << std::endl;
    std::vector<std::vector<double>>hampcc=eigenToNestedVector(hampchange);
    std::vector<int> jvecp;
    std::vector<int> parityvecp;
    for (int i=0;i<allbasisp.size();++i)
    {
        jvecp.push_back(allbasisp[i].sj.back());
        parityvecp.push_back(allbasisp[i].parity);
    }






    nucleus=nucleus2;
    allnucleus=allnucleus2;
    size1=nucleus.size();



    //2hezi
    std::vector<int> rordern={0,4,8,12,16,6};

    std::vector<std::vector<int>>rorderrn={rordern,};
    rvecall={0,4,8,12,16,6};
    rparityvec={1,1,1,1,1,-1};

    std::vector<std::pair<std::vector<int>, std::vector<int>>> catchpair2;
    catchpair2=generateValidPairs(4,rorderrn,51);
    std::vector<basis> allbasisp2;
    allbasisp2=calculateBasisWithRange(4, catchpair2);



    sizep=allbasisp2.size();
    sizepr=allbasisp2[0].r.size();
    std::vector<std::vector<std::vector<std::vector<double>>>> ystrallp2(
        sizep,  // 第一维：大小为 sizep
        std::vector<std::vector<std::vector<double>>>(
            sizepr,  // 第二维：大小为 sizepr
            std::vector<std::vector<double>>(
                size1,  // 第三维：大小为 size1
                std::vector<double>(size1, 0.0)  // 第四维：大小为 size1，元素初始化为 0.0
            )
        )
    );
    calculateystr(ystrallp2,ystrgetp2,allbasisp2);

    std::vector<basism> mtest2;
    std::vector<basism> mtest2_1;
    mtest2=calculateBasisform(4,rorderrn);
    mtest2_1=calculateBasisform_1(4,rorderrn);
    Eigen::MatrixXd m2;
    Eigen::MatrixXd m2_1;
    m2=ComputeTransformationMatrix( mtest2,allbasisp2);
    m2_1=ComputeTransformationMatrix( mtest2_1,allbasisp2);
    std::vector<std::vector<Eigen::MatrixXd>>ystrm2= computeystrm( mtest2, ystrgetp2);
    std::vector<std::vector<Eigen::MatrixXd>>ystrm2_1= computeystrm( mtest2_1, ystrgetp2);

    Eigen::MatrixXd overlapmn=caloverlapmmat(mtest2,mtest2,ystrm2,ystrm2);

    std::vector<std::vector<double> > overlapmch2=overlapchange(overlapmn,m2);
    std::cout << "overlapmch2" << std::endl;
    printMatrix(overlapmch2);


    std::vector<std::vector<double> > overlapp2;
    overlapp2=overlapmatcal(allbasisp2,ystrallp2);
    updateMatrices(allbasisp2, overlapmch2,ystrallp2);
    std::vector<std::vector<double>>schmitmatp2;
    schmitmatp2=gramSchmidtInOverlap(overlapmch2,allbasisp2,ystrallp2);
    Eigen::MatrixXd schmitmat2=convertToEigenMatrix(schmitmatp2);
    Eigen::MatrixXd overlap2=convertToEigenMatrix(overlapmch2);
    Eigen::MatrixXd overlap2cha=schmitmat2 * overlap2 * schmitmat2.transpose();
    std::cout << overlap2cha << std::endl;
    m2=ComputeTransformationMatrix( mtest2,allbasisp2);
    m2_1=ComputeTransformationMatrix( mtest2_1,allbasisp2);
    std::cout<< "m2" << std::endl;
    std::cout << m2 << std::endl;
    printBasisVectorOneLine(allbasisp2);



    std::map<int,std::map<int,Eigen::MatrixXd>>qmatnu= qmatjcal( q_nu, m2,V_it,allbasisp2,mtest2,ystrm2);
    std::map<int,std::map<int,Eigen::MatrixXd>>qmatnu_1= qmatjcal_1( q_nu, m2,m2_1,
        V_it,allbasisp2,mtest2,mtest2_1,ystrm2,ystrm2_1);
    std::cout<<"qmatnu = "<<qmatnu[2][0]<<std::endl;

    for (auto qmatcal : qmatnu)
    {
        int t=qmatcal.first;
        for (auto qmatcal2 : qmatcal.second)
        {
            int inum=qmatcal2.first;
            for (int i=0;i<allbasisp2.size();++i)
            {
                int j1=allbasisp2[i].sj.back();
                for (int j=0;j<allbasisp2.size();++j)
                {
                    int j2=allbasisp2[j].sj.back();
                    int num=(j1+j2+(t*2))/2;
                    if ((num & 1) != 0)
                    {
                        qmatnu[t][inum](i,j)=qmatnu_1[t][inum](i,j);
                    }
                }
            }
        }
    }

    std::map<int,std::map<int,Eigen::MatrixXd>> qmatnucha;
    for (auto qmatcal : qmatnu)
    {
        int t=qmatcal.first;
        std::cout<<"t = "<<t<<std::endl;
        for (auto qmatcal2 : qmatcal.second)
        {
            int inum=qmatcal2.first;
            qmatnucha[t][inum]=schmitmat2*qmatcal2.second* schmitmat2.transpose();
        }
    }


    double alpha2=4.33195830385423;
    std::vector<std::vector<double> > qstrgetp2;
    std::vector<int> qorderp2={4,6};
    qstrgetp2=qstrall(qorderp2,alpha2);

    int sizeq2=qorderp2.size();
    std::vector<std::vector<std::vector<double>>> qstrallp2(sizeq2, std::vector<std::vector<double>>(size1,
    std::vector<double>(size1, 0)));

    getystrbasis(qstrallp2, qstrgetp2, {0,1});

    std::vector<std::vector<double> >m1qget2=m1qstrall(0,-2.6782);
    std::vector<std::vector<std::vector<double>>> m1strallp2(1, std::vector<std::vector<double>>(size1,
    std::vector<double>(size1, 0)));
    getystrbasis(m1strallp2, m1qget2, {0});
    qstrallp2.push_back(m1strallp2[0]);
    qorderp2.push_back(2);
    sizeq2=qorderp2.size();
    std::vector<Eigen::MatrixXd> bemematcal2={};
    for (int tnum=0;tnum<qstrallp2.size();++tnum)
    {
        Eigen::MatrixXd qti=convertToEigenMatrix(qstrallp2[tnum]);
        Eigen::MatrixXd resultmat=bemecal(qti,m2,allbasisp2,mtest2,ystrm2,qorderp2[tnum]);
        Eigen::MatrixXd resultmat_1=bemecal_1(qti,m2,m2_1,allbasisp2,
            mtest2,mtest2_1,ystrm2,ystrm2_1,qorderp2[tnum]);
        int t=qorderp2[tnum];
        for (int i=0;i<allbasisp2.size();++i)
        {
            int j1=allbasisp2[i].sj.back();
            for (int j=0;j<allbasisp2.size();++j)
            {
                int j2=allbasisp2[j].sj.back();
                int num=(j1+j2+t)/2;
                if ((num & 1) != 0)
                {
                    resultmat(i,j)=resultmat_1(i,j);
                }
            }
        }
        resultmat=schmitmat2*resultmat* schmitmat2.transpose();
        bemematcal2.push_back(resultmat);
    }
    std::cout<<"bemematcal2 = "<<std::endl;
    for (const auto& bememat : bemematcal2)
    {
        std::cout << bememat << std::endl;
    }



    strength={1.0,0.57};
    std::vector<std::map<int, Matrix4D>>buildVValue2_2= processFile("D:/paircalpro/mschemecode/efct1copy.txt");
    Eigen::MatrixXd ham2cal=ham1( mtest2,ystrm2,buildVValue2_2,strength);
    std::cout << "ham2cal" << std::endl;
    std::cout << ham2cal << std::endl;
    Eigen::MatrixXd changeham2=hamchange(mtest2,allbasisp2,m2, 0, 0, ham2cal);
    std::cout << "hamchangecal2" << std::endl;
    std::cout << changeham2 << std::endl;
    Eigen::MatrixXd hamsige2= hamsige(mtest2,ystrm2,energyn);
    std::cout << "hamsige2" << std::endl;
    std::cout << hamsige2 << std::endl;
    Eigen::MatrixXd hamchangesig2=hamchange(mtest2,allbasisp2,m2, 0, 0, hamsige2);
    std::cout << "hamsigechange2" << std::endl;
    std::cout << hamchangesig2 << std::endl;
    Eigen::MatrixXd hamp2=changeham2+ hamchangesig2;
    std::cout << "hamp2" << std::endl;
    std::cout << hamp2 << std::endl;
    Eigen::MatrixXd hampchange2=schmitmat2*hamp2* schmitmat2.transpose();
    std::cout << "hampchange2" << std::endl;
    std::cout << hampchange2<< std::endl;
    std::vector<std::vector<double>>hampcc2=eigenToNestedVector(hampchange2);
    std::vector<int> jvecn;
    std::vector<int> parityvecn;
    for (int i=0;i<allbasisp2.size();++i)
    {
        jvecn.push_back(allbasisp2[i].sj.back());
        parityvecn.push_back(allbasisp2[i].parity);
    }
    std::map<int, std::vector<CoupledBasis>> coupleBasesall2;
    coupleBasesall2= coupleBases2(jvecn,jvecp,parityvecp,parityvecn);
    if (!outfile.is_open())
    {
        // 如果文件未打开，则尝试打开文件
        outfile.open("basis_output.txt",std::ios::app);
    }
    outfile<< "Coupled Bases for all2:" << std::endl;
    outfile<<coupleBasesall2<<std::endl;
    int sizecpbsize=coupleBasesall2.size();
    std::map<int, std::vector<std::vector<double>>> eigenre;
    for (const auto& [key, value] : coupleBasesall2)
    {
        // if (std::abs(key)>6)
        // {
        //     continue; // 只处理绝对值小于等于10的键
        // }

        int siz1=coupleBasesall2[key].size();
        int vitnum=0;

        std::vector<std::vector<double>> ham(siz1, std::vector<double>(siz1, 0));
        std::vector<std::vector<double>> ham2(siz1, std::vector<double>(siz1, 0));
        std::vector<std::vector<double>> ham3(siz1, std::vector<double>(siz1, 0));
        std::vector<std::vector<double>> bemcal1(siz1, std::vector<double>(siz1, 0));
        std::vector<std::vector<double>> bemcal2(siz1, std::vector<double>(siz1, 0));
        for (int i=0;i< siz1;i++)
        {
            for (int j=0;j< siz1;j++)
            {
                int n1=coupleBasesall2[key][i].ni;
                int n2=coupleBasesall2[key][j].ni;
                int p1=coupleBasesall2[key][i].pi;
                int p2=coupleBasesall2[key][j].pi;
                int jn=jvecn[n1];
                int jnp=jvecn[n2];
                int jp=jvecp[p1];
                int jpp=jvecp[p2];
                if (n1==n2)
                {
                    ham[i][j]= ham[i][j]+hampcc[p1][p2];
                }
                if (p1==p2)
                {
                    ham[i][j]=ham[i][j]+hampcc2[n1][n2];
                }
                int deltann=deltatwo(jn,jnp);
                int deltapp=deltatwo(jp,jpp);
                ham[i][j]=ham[i][j]*deltann*deltapp;
                for (auto qmatnuchat:qmatnucha)
                {
                    int tnum=qmatnuchat.first;
                    for (auto qmatcalsec:qmatnuchat.second)
                    {
                        int rownum=qmatcalsec.first;
                        double Vit=V_it[tnum](rownum,rownum);
                        int tt=2*tnum;
                        std::vector<int> sigl={jpp,jn,std::abs(key),tt};
                        int sig=sign_func(sigl);
                        double c2=util::wigner_6j(jnp,jpp,std::abs(key),jp,jn,tt)*mysqrt(jnp+1)*mysqrt(jpp+1)*Vit;
                        Eigen::MatrixXd matsliqnch= qmatnucha[tnum][rownum];
                        Eigen::MatrixXd matsliqpch= qmatpicha[tnum][rownum];
                        double c3=matsliqnch(n1,n2)*matsliqpch(p1,p2);
                        ham2[i][j]+=c2*c3*sig;
                    }
                }

                ham3[i][j]=ham[i][j]+ham2[i][j];



            }
        }

        // if (std::abs(key)<100)
        // {
        //     std::cout<<"ham"<<std::endl;
        //     printMatrix(ham);
        //     std::cout<<"ham2"<<std::endl;
        //     printMatrix(ham2);
        //     std::cout<<"ham3"<<std::endl;
        //     printMatrix(ham3);
        // }

        int nmat=ham3.size();

        auto mv_mul1 = [&](const std::vector<double>& in, std::vector<double>& out) {
            for(int i = 0; i < nmat; ++i) {
                for(int j = 0; j < nmat; ++j) {
                    out[i] +=ham3[i][j]*in[j];
                }
            }
        };

        LambdaLanczos<double> engine1(mv_mul1, nmat, false, 2); // true means to calculate the largest eigenvalue.
        std::vector<double> eigenvalues1;
        std::vector<std::vector<double>> eigenvectors1;
        engine1.run(eigenvalues1, eigenvectors1);
        eigenre[key]=eigenvectors1;
        std::cout << "J= " <<key<< std::endl;
        std::cout << "Eigenvalues: " << std::endl;
        for (size_t i = 0; i < eigenvalues1.size(); ++i) {
            std::cout << "Eigenvalue " << i + 1 << ": " << eigenvalues1[i] << std::endl;
        }

        // 输出特征向量
        std::cout << "Eigenvectors: " << std::endl;
        for (size_t i = 0; i < eigenvectors1.size(); ++i) {
            std::cout << "Eigenvector " << i + 1 << ": [ ";
            for (size_t j = 0; j < eigenvectors1[i].size(); ++j) {
                std::cout << eigenvectors1[i][j] << " ";
            }
            std::cout << "]" << std::endl;
        }
        //computeEigenvaluesAndEigenvectors(matchange1);

        double mmm=1;
    }


    int sizcpbsize=coupleBasesall2.size();
    int sizkey=eigenre.size();
    std::vector<Eigen::MatrixXd> bematrices(sizeq2, Eigen::MatrixXd::Zero(sizkey, sizkey));
    double ep=2.4;
    double en=1.3;
    int keynum=0;
    for (const auto& [key, value] : eigenre)
    {
        int jisum=std::abs(key);
        int siz1=coupleBasesall2[key].size();
        int vitnum=0;

        int keynum2=0;
        for (const auto& [key2, value2] : eigenre)
        {
            int jfsum=std::abs(key2);
            int siz2=coupleBasesall2[key2].size();
            std::vector<std::vector<std::vector<double>>> bmat1(
                sizeq2, std::vector<std::vector<double>>(siz1, std::vector<double>(siz2, 0.0)));
            std::vector<std::vector<std::vector<double>>> bmat2(
                sizeq2, std::vector<std::vector<double>>(siz1, std::vector<double>(siz2, 0.0)));
            for (int i=0;i< siz1;i++)
            {
                for (int j=0;j< siz2;j++)
                {
                    int n1=coupleBasesall2[key][i].ni;
                    int n2=coupleBasesall2[key2][j].ni;
                    int p1=coupleBasesall2[key][i].pi;
                    int p2=coupleBasesall2[key2][j].pi;
                    int jn=jvecn[n1];
                    int jnp=jvecn[n2];
                    int jp=jvecp[p1];
                    int jpp=jvecp[p2];


                    int deltann=deltatwo(jn,jnp);
                    int deltapp=deltatwo(jp,jpp);
                    for (int qmatnum=0; qmatnum<sizeq2;++qmatnum)
                    {
                        int t=qorderp2[qmatnum];
                        Eigen::MatrixXd bmatcal1=bemematcal[qmatnum];
                        Eigen::MatrixXd bmatcal2=bemematcal2[qmatnum];
                        double calpi=deltann*overlap2cha(n1,n2)
                        *ufrom6_j(jnp,jp,jfsum,t,jisum,jpp)*bmatcal1(p1,p2);
                        double calnu=deltapp*overlap1cha(p1,p2)*std::pow(-1,(jn-jnp+jfsum-jisum)/2)
                        *ufrom6_j(jpp,jn,jfsum,t,jisum,jnp)*bmatcal2(n1,n2);
                        bmat1[qmatnum][i][j]=calpi;
                        bmat2[qmatnum][i][j]=calnu;
                    }
                }
            }
            if (key==0 && key2==-6)
            {
                std::cout << "bmat1" << std::endl;
                for (int qmatnum=0; qmatnum<sizeq2;++qmatnum)
                {
                    std::cout << "qmatnum = " << qmatnum << std::endl;
                    printMatrix(bmat1[qmatnum]);
                }
                std::cout << "bmat2" << std::endl;
                for (int qmatnum=0; qmatnum<sizeq2;++qmatnum)
                {
                    std::cout << "qmatnum = " << qmatnum << std::endl;
                    printMatrix(bmat2[qmatnum]);
                }

            }
            std::vector<double> eigencal1=eigenre[key][0];
            std::vector<double> eigencal2=eigenre[key2][0];

            for(int qmatnum=0; qmatnum<sizeq2;++qmatnum)
            {
                double epp=ep;
                double enn=en;
                int t=qorderp2[qmatnum]/2;
                int fact=1;
                if ((t%2)==1 )
                {
                    fact=-1;
                }
                if (qmatnum==2 )
                {
                    epp=1;
                    enn=1;
                }
                double sum=0.0;
                Eigen::MatrixXd matcal111(eigencal1.size(),eigencal2.size());
                Eigen::MatrixXd matcal222(eigencal1.size(),eigencal2.size());

                for (int i=0;i<eigencal1.size();++i)
                {
                    for (int j=0;j<eigencal2.size();++j)
                    {
                        double cal1=eigencal1[i];
                        double cal2=eigencal2[j];
                        sum+=cal1*bmat1[qmatnum][i][j]*cal2*epp;
                        sum+=cal2*bmat2[qmatnum][i][j]*cal1*enn*fact;

                        if (key==0 && key2==-6)
                        {
                            matcal111(i,j)=cal1*bmat1[qmatnum][i][j]*cal2*epp;
                            matcal222(i,j)=cal2*bmat2[qmatnum][i][j]*cal1*enn*fact;
                            if (i==0 && j==0)
                            {
                                std::cout<< "cal1 = " << cal1 << std::endl;
                                std::cout<< "cal2 = " << cal2 << std::endl;
                                std::cout<<"enn = " << enn << std::endl;
                                std::cout<<"epp = " << epp << std::endl;
                                std::cout<<"bmat1"<<bmat1[qmatnum][i][j]<<std::endl;
                            }
                        }

                    }
                }
                if (key==0 && key2==-6)
                {
                    if (qmatnum==1)
                    {
                        std::cout << "matcal111" << std::endl;
                        std::cout << matcal111 << std::endl;
                        std::cout << "matcal222" << std::endl;
                        std::cout << matcal222 << std::endl;
                        std::cout<<"fact"<<std::endl;
                        std::cout << fact << std::endl;
                    }
                }
                sum=std::pow(sum,2.0)*(jfsum+1.0)/(jisum+1.0);
                bematrices[qmatnum](keynum,keynum2)=sum;

            }
            keynum2++;

        }
        keynum++;
    }
    std::cout << "bematrices" << std::endl;
    std::cout<< bematrices[0] << std::endl;
    std::cout << "bematrices" << std::endl;
    std::cout<< bematrices[1] << std::endl;
    std::cout << "bematrices" << std::endl;
    std::cout<< bematrices[2] << std::endl;




    return 0;



}