//
// Created by wang- on 24-12-31.
//
#include "../inc/moe.h"
bool isOdd(int x)
{
    return x % 2 != 0; // 如果 x 除以 2 的余数不为 0，则是奇数，返回 true
}

void removestr(std::vector<std::vector<std::vector<double> > > &arr, int i)
{
    // 检查索引是否有效
    if (i >= 0 && i < arr.size())
    {
        arr.erase(arr.begin() + i); // 删除第一维索引为 i 的元素
    } else
    {
        //std::cout << "Index out of bounds: " << i << std::endl;
    }
}
 int sign_func(std::vector<int> x)
{
    int sum = 0;
    for (int num: x)
    {
        // 对 vec 中的每个元素 num
        sum += num; // 累加元素
    }
    sum = sum / 2;
    if (isOdd(sum))
    {
        return -1;
    }  else { return 1; }
}
 double mysqrt(int x)
{

    if (x < 17)
    {
        return sqrttab[x];
    } else { return sqrt(x); };
}



int deltatwo(int a, int b)
{
    if (a == b) { return 1; } else { return 0; };
}
 double calculate_sum(const std::vector<std::vector<std::vector<double> > > &v1,
                            const std::vector<std::vector<std::vector<double> > > &v2,
                            int rk, int sn, int rip, int j, int t)
{
    double result = 0.0;

    // 遍历所有可能的 b
    for (int b = 0; b < jmax; ++b)
    {
        // 根据给定条件，k1 = j2 = b
        int k1 = b;
        int j2 = b;

        // 计算 v1[i1][j1][k1] 和 v2[i2][j2][k2]，并乘以系数 j_factor
        result += v1[rk][rip][k1] * v2[sn][b][j] * util::wigner_6j(rk, sn, t, j, rip, b);
    }

    return result;
}
bool is_in_array(const std::vector<int> &arr, int value)
{
    // 使用 std::find 查找元素
    return std::find(arr.begin(), arr.end(), value) != arr.end();
}
 std::vector<int> is_get_in(const std::vector<std::vector<int> > &arr, int num, int jip)
{
    std::vector<int> result = {};
    for (int i = 0; i < arr.size(); i++)
    {
        if (jip == arr[i][num])
        {
            result.push_back(i);
        }
    }
    return result;
}
// std::vector<int> genmvec(int a, int b) {
//     static std::array<int, 21> temp; // 预分配最多 21 个元素
//     int min_c = std::abs(a - b);
//     int max_c = a + b;
//
//     int count = 0;
//     for (int c = min_c; c <= max_c; c += 2) {
//         temp[count++] = c; // 直接赋值，避免 push_back()
//     }
//
//     return std::vector<int>(temp.begin(), temp.begin() + count); // 直接构造 vector
// }
std::vector<int> genmvec(int a, int b)
{
    std::vector<int> result;

    int min_c = std::abs(a - b);
    int max_c = a + b;

    int count = (max_c - min_c) / 2 + 1; // 计算元素个数
    result.reserve(count);  // 预分配内存，提高 push_back 性能

    for (int c = min_c; c <= max_c; c += 2)
    {
        result.emplace_back(c); // 直接构造，减少拷贝
    }

    return result;
}
//  std::vector<int> genmvec(int a, int b)
// {
//     std::vector<int> result;
//
//     // 计算角动量的最小值和最大值
//     int min_c = std::abs(a - b); // 最小值 |a - b|
//     int max_c = a + b; // 最大值 a + b
//
//     // 根据规则生成范围，步长为 1
//     for (int c = min_c; c <= max_c; c = c + 2)
//     {
//         result.push_back(c); // 将 c 添加到数组中
//     }
//
//     return result;
// }
 std::vector<int> genmrivec(std::vector<int> x, int y)
{
    std::vector<int> result;
    auto maxx = std::max_element(x.begin(), x.end());
    auto minx = std::min_element(x.begin(), x.end());
    int dmin = *minx - y;
    int dmax = *maxx - y;
    int minre;
    int maxre = *maxx + y;
    int admin = std::abs(dmin);
    int admax = std::abs(dmax);
    if (dmin <= 0 && dmax >=0)
    {
        if (isOdd(dmax))
        {
            minre =1;
        }else
        {
            minre = 0;
        }
    } else if (admin > admax)
    {
        minre = admax;
    } else
    {
        minre = admin;
    }
    for (int c = minre; c <= maxre; c = c + 2)
    {
        result.push_back(c);
    }
    return result;
}
std::vector<std::pair<int, int> > findcpy(const std::vector<Nucleus> &nuclei, int r)
{
    std::vector<std::pair<int, int> > results;

    // 遍历所有可能的核子对 (i, j)
    for (size_t i = 0; i < nuclei.size(); ++i)
    {
        for (size_t j = i + 1; j < nuclei.size(); ++j)
        {
            int j1 = nuclei[i].j;
            int j2 = nuclei[j].j;

            // 检查耦合规则：|j1 - j2| <= r <= j1 + j2
            int minR = std::abs(j1 - j2);
            int maxR = j1 + j2;
            if (r >= minR && r <= maxR)
            {
                // 检查宇称规则：只有两个核子的宇称和为 +1 时才有效
                int paritySum = nuclei[i].parity + nuclei[j].parity;
                if (paritySum == 2)
                {
                    // 宇称和为正
                    results.emplace_back(i, j); // 存储核子对的索引
                }
            }
        }
    }

    return results;
}
double getphi(const std::vector<std::vector<std::vector<double> > > &A,
                     const std::vector<std::vector<std::vector<double> > > &B,
                     basis rjl, int k)
{
    double sum = 0.0;
    int N=rjl.r.size();

    // 遍历 A 和 B 的第一维
    for (size_t a = 0; a < A[k-1].size(); ++a)
    {
        // 第一维键 r
        for (size_t b = 0; b < A[k-1][a].size(); ++b)
        {
            sum=sum+A[k-1][a][b]*B[N-1][a][b];
        }
    }
    int jk = rjl.sj[k - 1];
    int jk_1;
    if (k == 1) { jk_1 = rjl.j; } else { jk_1 = rjl.sj[k - 2]; }
    int rk = rjl.r[k - 1];
    double cp = 2 * std::pow(-1, (jk - jk_1 - rk)/2.0) * mysqrt(jk + 1) / mysqrt(jk_1 + 1);

    return cp * sum; // 返回求和的结果
}
double ufrom6_j(int a, int b, int d, int c, int e, int f)
{
    std::vector<int> x = {a, b, d, c};
    int signx = sign_func(x);
    double u = signx *mysqrt(e+1)*mysqrt(f+1) * util::wigner_6j(a, b, e, c, d, f);
    return u;
};
double hfrom6_j(struct basis rjl, int k, int s, int L1, int L2)
{
    if (s == 0) { return 1; }
    int N = rjl.sj.size();
    std::vector<int> r = rjl.r;
    std::vector<int> sj = rjl.sj;
    int jk_1 = rjl.j;
    if (k != 1)
    {
        jk_1 = sj[k - 2];
    }
    std::vector<int> x = {jk_1, L1, -sj[k - 1], -L2};
    int signx = sign_func(x);

    double u = signx * ufrom6_j(r[k - 1], L1, jk_1, s, L2, sj[k - 1]);
    return u;

};
 double gfrom6_j(struct basis rjl, int k, int s, int t, int lk_1)
{
    int jk_1;
    if (k == 0) { return 0; } else if (k == 1)
    {
        jk_1 = rjl.j;
    } else { jk_1 = rjl.sj[k - 2]; }
    int jk = rjl.sj[k - 1];
    int rk = rjl.r[k - 1];
     if (-ufrom6_j(rk, s, jk_1, lk_1, t, jk)==0)
     {
         double ggg=-ufrom6_j(rk, s, jk_1, lk_1, t, jk);
     }
    return -ufrom6_j(rk, s, jk_1, lk_1, t, jk);
};
double mfrom6_j(struct basis rjl, int k, int t, int rkp, int lii)
{
    int rk ;
    int jk_1;
    int jk ;
    if (k == 0) { return 1; }
    else
    {
        rk=rjl.r[k - 1];
        jk=rjl.sj[k - 1];
        if (k==1)
        {
            jk_1 = rjl.j;
        }else
        {
            jk_1 = rjl.sj[k - 2];
        }
    }
    double g=ufrom6_j(rk, t, jk_1, lii, rkp, jk);
    return ufrom6_j(rk, t, jk_1, lii, rkp, jk);
}
double qfrom6_j(struct basis rjl, int i, int t, int li, int li_1)
{
    int ji_1 = rjl.j;
    if (i != 1)
    {
        ji_1 = rjl.sj[i - 2];
    }
    int ji = rjl.sj[i - 1];
    int ri = rjl.r[i - 1];
    int sign = std::pow(-1, (ji_1 + li - ji - li_1)/2.0);
    return sign * ufrom6_j(ri, li, ji_1, t, li_1, ji);
}
std::vector<int> getsamerange(const std::vector<int> &vec1, const std::vector<int> &vec2)
{
    std::unordered_set<int> set1(vec1.begin(), vec1.end()); // 使用set来去重并快速查找
    std::vector<int> intersection;

    for (int num: vec2)
    {
        if (set1.find(num) != set1.end())
        {
            // 如果vec2中的元素在set1中存在
            intersection.push_back(num);
        }
    }

    return intersection;
};
std::vector<std::vector<double> > getz(std::vector<std::vector<std::vector<double> > > & ystrall1,
                                              std::vector<std::vector<std::vector<double> > > & ystrall2,
                                              int i, int k, int rip, int rk, int ri, basis rjl, int n, int sn, int t)
{
    int nucnum = nucleus.size();
    int N2 = rjl.r.size();
    double c4 = 4 * mysqrt(rk + 1) * mysqrt(ri + 1) * mysqrt(sn + 1) * mysqrt(t + 1);
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
                        double y1 = ystrall2[n - 1][b_m][b_n] * ystrall1[k - 1][i_m][b_m] * ystrall1[i - 1][i_n][b_n];
                        y1 = y1 * util::wigner_6j(rk, sn, t, jbp, j_m, jb) * util::wigner_6j(ri, t, rip, j_m, j_n, jbp);
                        double y2=ystrall2[n - 1][b_m][b_n] * ystrall1[k - 1][i_n][b_m] * ystrall1[i - 1][i_m][b_n];
                        y2=y2*util::wigner_6j(rk, sn, t, jbp, j_n, jb) * util::wigner_6j(ri, t, rip, j_n, j_m, jbp);
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

