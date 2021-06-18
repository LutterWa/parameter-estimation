#ifndef __PARAMETER_IDENTIFICATION__
#define __PARAMETER_IDENTIFICATION__
#include <vector>
#include <algorithm>
#include <functional>

using std::vector;

enum ESTIMATE_METHOD {CK, MLE};

class PARAMETER_IDENTIFICATION
{
    private:
    enum ESTIMATE_METHOD method = CK;  //参数辨识方法，CK=C-k法，MLE=最大似然估计法
    double h = 0.002;  //采样时间间隔
    double Area = 0.018869;  //参考面积
    double m0 = 47.5;  //弹重
    vector<double> dery(const vector<double> y);  // 导弹运动方程组
    vector<double> rk4(const vector<double> y);  //四阶龙格库塔法
    vector<double> chapman_kirk(const vector<vector<double>> exp_batch, const int ck);  //C-K法
    double max_likelihood_estimate(const vector<vector<double>> exp_batch, const int ck);  //最大似然估计法
    
    public:
    PARAMETER_IDENTIFICATION();
    ~PARAMETER_IDENTIFICATION();
    PARAMETER_IDENTIFICATION(enum ESTIMATE_METHOD method, double h=0.002, double Area=0.018869, double m0=47.5);
    vector<vector<double>> ak;  //待辨识参数
    vector<double> ma;  //马赫数向量
    vector<double> x;
    vector<double> y;
    vector<double> z;
    vector<double> v;
    int load_data();  //加载辨识数据
    int save_data();  //保存辨识结果
    int parameter_estimation(void);
};

typedef struct{  //气象参数结构体
    double density;  //空气密度
    double soundsp;  //音速
    double density_dot;  //空气密度关于高度的导数
}ATMOSPHERE;

vector<vector<double>> LUP_solve_inverse(const vector<vector<double>> A);  //LUP求逆(将每列b求出的各列x进行组装)
vector<vector<double>> transpose(const vector<vector<double>> A);  //矩阵转置
double absdet(const vector<vector<double>> A); //矩阵行列式绝对值

vector<vector<double>> mul(const vector<vector<double>> A, const vector<vector<double>> B);  //矩阵乘法
vector<double> mul(const vector<vector<double>> A, const vector<double> B);  //矩阵-列向量乘法
vector<double> mul(const vector<double> B, const vector<vector<double>> A);  //行向量-矩阵乘法
vector<vector<double>> mul(const double a, const vector<vector<double>> B);  //常数-矩阵乘法
vector<vector<double>> mul(const vector<vector<double>> B, const double a);  //矩阵-常数乘法
vector<double> mul(const double a, const vector<double> B);  //常数-向量乘法
vector<double> mul(const vector<double> B, const double a);  //向量-常数乘法
double mul(const vector<double> A, const vector<double> B);  //向量-向量内积

vector<vector<double>> add(const vector<vector<double>> A, const vector<vector<double>> B);  //矩阵加法
vector<double> add(const vector<double> A, const vector<double> B);  //向量加法

vector<vector<double>> eye(const double dim, const double k=1.0);  //创建对角矩阵

#endif /*__PARAMETER_IDENTIFICATION__*/