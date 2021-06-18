#include <stdio.h>
#include <math.h>
#include "algorithm.h"

static int LUP_Descomposition(vector<vector<double>> &A, vector<vector<double>> &L, vector<vector<double>> &U, vector<int> &P);  //LUP分解
static vector<double> LUP_Solve(const vector<vector<double>> L, const vector<vector<double>> U, const vector<int> P, const vector<double> b);  //LUP求解方程

static int LUP_Descomposition(vector<vector<double>> &A, vector<vector<double>> &L, vector<vector<double>> &U, vector<int> &P)  //LUP分解
{
    int N = A.size();
    P.resize(N, 0);
    int row=0;
    for(int i=0;i<N;i++)
        P[i]=i;  //初等变换基准向量
    for(int i=0;i<N-1;i++)
    {
        {//1.通过初等变换交换行序, 使Gauss消去法有较大的主元
            double p=0.;
            for(int j=i;j<N;j++)
            {
                if(abs(A[j][i])>p)
                {
                    p=abs(A[j][i]);
                    row=j;
                }
            }
            if(p <= 0)  //全零行,矩阵奇异,无法计算逆
            {
                printf("Inverse matrix not exist.\n");
                return -1;
            }
            {//交换P[i]和P[row]
                int tmp=P[i];
                P[i]=P[row];
                P[row]=tmp;
            }
            {//交换A[i]和 A[row]
                vector<double> tmp;
                tmp=A[i];
                A[i]=A[row];
                A[row]=tmp;
            }
        }

        {//2.通过初等变换构造上三角矩阵
            double u=A[i][i],l=0.;
            for(int j=i+1;j<N;j++)
            {
                l=A[j][i]/u;
                A[j][i]=l;  //消去后结果应为0, 此处为暂存初等变换的比例值
                for(int k=i+1;k<N;k++)
                {
                    A[j][k]=A[j][k]-A[i][k]*l;
                }
            }
        }
    }

    // 构造L和U
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<=i;j++)
        {
            if(i!=j)
            {
                L[i][j]=A[i][j];
            }
            else
            {
                L[i][j]=1;
            }
        }
        for(int k=i;k<N;k++)
        {
            U[i][k]=A[i][k];
        }
    }
    return 0;
}

static vector<double> LUP_Solve(const vector<vector<double>> L, const vector<vector<double>> U, const vector<int> P, const vector<double> b)  //LUP求解方程
{
    int N = L.size();
    vector<double> x(N, 0);
    vector<double> y(N, 0);

    //正向替换-求解下三角矩阵Ly=b, L对角线元素为1
    for(int i = 0;i < N; i++)
    {
        y[i] = b[P[i]];
        for(int j = 0; j < i; j++)
        {
            y[i] = y[i] - L[i][j]*y[j];
        }
    }
    //反向替换-求解上三角矩阵Ux=y
    for(int i = N-1;i >= 0; i--)
    {
        x[i]=y[i];
        for(int j = N-1; j > i; j--)
        {
            x[i] = x[i] - U[i][j]*x[j];
        }
        x[i] /= U[i][i];
    }
    return x;
}

vector<vector<double>> LUP_solve_inverse(const vector<vector<double>> A)  //LUP求逆(将每列b求出的各列x进行组装)
{
    int N = A.size();
    vector<vector<double>> A_mirror = A;  //创建矩阵A的副本，注意不能直接用A计算，因为LUP分解算法已将其改变
    vector<vector<double>> inv_A(N, vector<double>(N, 0));  //最终的逆矩阵（还需要转置）
    vector<double> inv_A_each(N, 0);//矩阵逆的各列
    vector<double> b(N, 0);//b阵为单位矩阵B阵的列矩阵分量
    for(int i=0;i<N;i++)
    {
        vector<vector<double>> L(N, vector<double>(N, 0));
        vector<vector<double>> U(N, vector<double>(N, 0));
        vector<int> P(N, 0);

        //构造单位矩阵的每一列
        b = vector<double>(N, 0);
        b[i]=1;

        //每次都需要重新将A复制一份
        A_mirror = A;

        if(LUP_Descomposition(A_mirror,L,U,P) < 0)
            return vector<vector<double>>(0);

        inv_A_each=LUP_Solve(L,U,P,b);
        inv_A[i] = inv_A_each;//将各列拼接起来
    }
    inv_A = transpose(inv_A);//由于现在根据每列b算出的x按行存储，因此需转置
    return inv_A;
}

double absdet(const vector<vector<double>> A)
{
    double d = 1.;
    int N = A.size();
    vector<vector<double>> A_mirror = A;  //创建矩阵A的副本，注意不能直接用A计算，因为LUP分解算法已将其改变
    vector<vector<double>> L(N, vector<double>(N, 0));
    vector<vector<double>> U(N, vector<double>(N, 0));
    vector<int> P(N, 0);

    if(LUP_Descomposition(A_mirror,L,U,P) < 0)
        return 0;

    for(int i=0; i<N; i++)
        d *= U[i][i];
    return fabs(d);
}
/****************************************************************************************************/
vector<vector<double>> transpose(const vector<vector<double>> A)  //矩阵转置
{
    int m = A.size();
    if(m > 0)
    {
        int n = A[0].size();
        vector<vector<double>> B(n, vector<double>(m, 0));
        for(int i = 0; i < m; i++)
            for(int j = 0; j < n; j++)
                B[j][i] = A[i][j];
        return B;
    }
    return vector<vector<double>>(0);
}
/****************************************************************************************************/
vector<vector<double>> mul(const vector<vector<double>> A, const vector<vector<double>> B)  //矩阵乘法
{
    int M = A.size(), N = A[0].size(), P = B.size(), Q = B[0].size();
    if(N != P)
    {
        printf("Matrix dim error!");
        return vector<vector<double>>(0, vector<double>(0));
    }
    vector<vector<double>> C(M, vector<double>(Q, 0));
    for(int i=0;i<M;i++)
        for(int j=0;j<Q;j++)
            for(int k=0;k<N;k++)
                C[i][j] += A[i][k]*B[k][j];
    return C;
}

vector<double> mul(const vector<vector<double>> A, const vector<double> B)  //矩阵-列向量乘法
{
    int M = A.size(), N = A[0].size(), P = B.size();
    if(N != P)
    {
        printf("Vector dim error!");
        return vector<double>(0);
    }
    vector<double> C(M, 0);
    for(int i=0;i<M;i++)
        for(int j=0;j<N;j++)
            C[i] += A[i][j]*B[j];
    return C;
}
vector<double> mul(const vector<double> B, const vector<vector<double>> A)  {return mul(A, B);}  //向量-矩阵乘法

vector<vector<double>> mul(const double a, const vector<vector<double>> B)  //常数-矩阵乘法
{
    int M = B.size(), N = B[0].size();
    vector<vector<double>> C(M, vector<double>(N, 0));
    for(int i=0;i<M;i++)
        for(int j=0;j<N;j++)
            C[i][j] = a*B[i][j];
    return C;
}
vector<vector<double>> mul(const vector<vector<double>> B, const double a)  {return mul(a, B);}  //常数-矩阵乘法

vector<double> mul(const double a, const vector<double> B)  //常数-向量乘法
{
    int M = B.size();
    vector<double> C(M, 0);
    for(int i=0;i<M;i++)
        C[i] = a*B[i];
    return C;
}
vector<double> mul(const vector<double> B, const double a)  {return mul(a, B);}  //常数-向量乘法

double mul(const vector<double> A, const vector<double> B)  //向量-向量内积
{
    int M = A.size(), N = B.size();
    if(M != N)
    {
        printf("Vector dim error!");
        return 0;
    }
    double C = 0.;
    for(int i=0;i<M;i++)
        C += A[i]*B[i];
    return C;  
}  
/****************************************************************************************************/
vector<vector<double>> eye(double dim, double k)  //创建对角矩阵
{
    vector<vector<double>> A(dim, vector<double>(dim, 0));
    for (int i=0; i<dim; i++)
        A[i][i] = k;
    return A;
}
/****************************************************************************************************/
vector<vector<double>> add(const vector<vector<double>> A, const vector<vector<double>> B)  //矩阵加法
{
    int M = A.size(), N = A[0].size(), P = B.size(), Q = B[0].size();
    if(M != P || N != Q)
    {
        printf("Matrix dim error!");
        return vector<vector<double>>(0, vector<double>(0));
    }
    vector<vector<double>> C(M, vector<double>(N, 0));
    for(int i=0;i<M;i++)
        for(int j=0;j<N;j++)
                C[i][j] = A[i][j] + B[i][j];
    return C;
}
vector<double> add(const vector<double> A, const vector<double> B)  //向量加法
{
    int M = A.size(), N = B.size();
    if(M != N)
    {
        printf("Vector dim error!");
        return A;
    }
    vector<double> C(M, 0);
    for(int i=0;i<M;i++)
        C[i] = A[i]+B[i];
    return C;
}
/****************************************************************************************************/
/*
int main()
{
    const int N = 5;  //矩阵行数/列数
    vector<vector<double>> A(N, vector<double>(N,0));

    for(int i=0; i<N ;i++)
    {
        for(int j=0; j<N;j++)
        {
            A[i][j] = rand() * 0.001;
        }
    }

    vector<vector<double>> E_test(N, vector<double>(N,0));
    vector<vector<double>> invOfA(N, vector<double>(N,0));
    invOfA=LUP_solve_inverse(A);

    E_test=mul(A,invOfA);  //验证精确度

    printf("Matrix A:\n");
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            printf("%lf\t",A[i][j]);
        }
        putchar('\n');
    }

    printf("inv_A:\n");
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            printf("%lf\t",invOfA[i][j]);
        }
        putchar('\n');
    }

    printf("E_test:\n");
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            printf("%.1f\t",E_test[i][j]);
        }
        putchar('\n');
    }
    return 0;
}
*/