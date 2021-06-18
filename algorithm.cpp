#include <stdio.h>
#include <math.h>
#include "algorithm.h"

#define g 9.81  //重力加速度

#define k_num 2//待辨识参数项数
#define y_num 6//系统状态参数项数
#define ma_num 0.1 //马赫数间隔


static ATMOSPHERE atmosStd(double height);  //获取气象参数
static double kalman_filter(const double z_measure, double &z_last, double &pz_last);  //卡尔曼滤波

PARAMETER_IDENTIFICATION::PARAMETER_IDENTIFICATION(){}
PARAMETER_IDENTIFICATION::~PARAMETER_IDENTIFICATION(){}
PARAMETER_IDENTIFICATION::PARAMETER_IDENTIFICATION(enum ESTIMATE_METHOD method, double h, double Area, double m0)
{
    this->method = method;
    this->h = h;
    this->Area = Area;
    this->m0 = m0;
}

vector<double> PARAMETER_IDENTIFICATION::dery(const vector<double> Y)  //待辨识参数的积分方程组
{
    vector<double> dy(y_num, 0);
    double v = Y[0];  //速度
    double theta = Y[1];  //弹道倾角
    double psi = Y[2];  //弹道偏角
    double y = Y[4];  //飞行高度
    const ATMOSPHERE atmosphere = atmosStd(y);
    double Ma = v / atmosphere.soundsp;

    vector<double>::iterator it;
    //提取待辨识系数初值
    if(Ma >= ma.back())
        it = ma.end() - 1;
    else if (Ma <= ma[0])
        it = ma.begin() + 1;
    else
        it = find_if(ma.begin(), ma.end(), std::bind2nd(std::greater<double>(),Ma));
    vector<double> k = ak[std::distance(ma.begin(), it)-1];

    double CD = k[0];  //阻力系数
    double CL = k[1];  //升力系数
    double CZ = k[1];  //侧向力系数
    double qS = 0.5*atmosphere.density*v*v*Area;  //动压×参考面积

	dy[0] = -(qS*CD + m0*g*sin(theta))/m0;  //v
	dy[1] = (qS*CL - m0*g*cos(theta))/m0/v;  //theta
	dy[2] = -qS*CZ / (m0*v*cos(theta));  //psi
	dy[3] = v*cos(theta)*cos(psi);  //x
	dy[4] = v*sin(theta);  //y
	dy[5] = -v*cos(theta)*sin(psi);  //z

    return dy;
}

vector<double> PARAMETER_IDENTIFICATION::rk4(const vector<double> Y)  //四阶龙格库塔
{
    vector<double>k1 = mul(h, dery(Y));
    vector<double>k2 = mul(h, dery(add(Y, mul(0.5, k1))));
    vector<double>k3 = mul(h, dery(add(Y, mul(0.5, k2))));
    vector<double>k4 = mul(h, dery(add(Y, k3)));

    return mul(1./6., add(add(k1, k4), mul(2., add(k2, k3))));
}

int PARAMETER_IDENTIFICATION::load_data()  //加载辨识数据
{
    FILE *fp = fopen("./output.dat", "r");
    if(fp == NULL)
    {
        perror("Open file error");
        return -1;
    }
    int ch=fgetc(fp);
    if(ch==EOF)
    {
       printf("The file is empty\n");
       fclose(fp);
       return -1;
    }
    rewind(fp);  //让文件内部的指针重新指向文件开头

    //卡尔曼滤波临时变量
    double x_last=0, px_last=1;
    double y_last=0, py_last=1;
    double z_last=0, pz_last=1;
    double v_last=0, pv_last=1;
    double ma_min = 3, ma_max = 0;
    int filter_flag = 0;

    while(!feof(fp))
    {
        double Y[8], alpha, CD, MA, q, aa, ad, aw, ba, bd, wm; //文件中每行共18个参数
        fscanf(fp, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf", \
        &Y[0],&Y[1],&Y[2],&Y[3],&Y[4],&Y[5],&Y[6],&Y[7],&alpha,&CD,&MA,&q,&aa,&ad,&aw,&ba,&bd,&wm);
        
        if(filter_flag == 0)
        {
            filter_flag = 1;
            x_last = Y[4];
            y_last = Y[5];
            z_last = Y[6];
            v_last = Y[1];
        }

        Y[4] = kalman_filter(Y[4], x_last, px_last);
        x.push_back(Y[4]);
        Y[5] = kalman_filter(Y[5], y_last, py_last);
        y.push_back(Y[5]);
        Y[6] = kalman_filter(Y[6], z_last, pz_last);
        z.push_back(Y[6]);
        Y[1] = kalman_filter(Y[1], v_last, pv_last);
        v.push_back(Y[1]);

        double Ma = v.back() / atmosStd(y.back()).soundsp;
        if(Ma < ma_min)  {ma_min = Ma;}
        if(Ma > ma_max)  {ma_max = Ma;}
    }
    fclose(fp);
    //意外的文件尾部处理:丢弃最后一行数据
    x.pop_back();
    y.pop_back();
    z.pop_back();
    v.pop_back();

    printf("data load_done! data size = %d\n", x.size());

    ma_max = ceil(ma_max / ma_num);  //向上取整
    ma_min = floor(ma_min / ma_num);  //向下取整

    for(int i = ma_min; i<= ma_max; i++)
    {
        ma.push_back(i * ma_num);
    }
    ak.resize(ma.size() - 1, vector<double>(k_num, 0.));  //为ak赋初值，结构：每一个Ma对应一组待辨识参数ak

    return 0;
}

int PARAMETER_IDENTIFICATION::save_data()  //保存辨识结果
{
    FILE *fp = fopen("./result.dat", "w");
    if(fp == NULL)
    {
        perror("Open file error");
        return -1;
    }
    fprintf(fp, "Current method is %s\nMa\t\tCD\t\tCL\n", (method == CK ? "CK" : "MLE"));
    for(int i=0; i<ak.size(); i++)
    {
        fprintf(fp,"%.1f-%.1f\t", ma[i], ma[i]+0.1);
        for(int j=0; j<ak[0].size(); j++)
            fprintf(fp, "%.4f\t", ak[i][j]);
        fprintf(fp, "\n");
    }
    return 0;
}

int PARAMETER_IDENTIFICATION::parameter_estimation(void)  //参数估计
{
    vector<double> Y_exp(y_num, 0);  //飞行器状态向量
    vector<vector<double>> exp_batch;  //状态向量缓冲区
    double head_flag = 1, head_Ma = 0;  //状态向量对应的马赫数
    int ck = -1;
    for(int i = 0;i < x.size()-1; i++)  //遍历样本, 按马赫数分段更新待辨识参数
    {
        {
            double dy = (y[i+1] - y[i]) / h;
            double dz = (z[i+1] - z[i]) / h;

            Y_exp[0] = v[i];  //v
            Y_exp[1] = asin(dy / v[i]);  //theta
            Y_exp[2] = asin(dz / v[i] / cos(Y_exp[1]));  //psi
            Y_exp[3] = x[i];  //x
            Y_exp[4] = y[i];  //y
            Y_exp[5] = z[i];  //z
            
            exp_batch.push_back(Y_exp);
        }
        double Ma = v[i] / atmosStd(y[i]).soundsp;  //计算样本的马赫数
        if(head_flag == 1)  //起始样本
        {
            head_flag = 0;
            auto it = find_if(ma.begin(), ma.end(), std::bind2nd(std::greater<double>(),Ma));
            head_Ma = *(it-1);  //样本区间:[it-1,it]
            int ck_1 = ck;
            ck = std::distance(ma.begin(), it)-1;
            if(ck_1 >= 0)
                ak[ck]=ak[ck_1];
        }
        else if(Ma < head_Ma || Ma > head_Ma + ma_num)  //当前样本马赫数超出了样本区间
        {
            head_flag = 1;
            i--;
            exp_batch.pop_back();
            double epsilon = 1, epsilon_1 = 1;
            while(epsilon > 0.0001)
            {
                if(method == CK)
                {
                    vector<double> B = chapman_kirk(exp_batch, ck);
                    epsilon = mul(B,B);
                    if(B.size() <= 0)
                    return -1;
                }
                else if(method == MLE)
                {
                    double J = max_likelihood_estimate(exp_batch, ck);
                    epsilon = fabs(1 - J / epsilon_1);
                    epsilon_1 = J;
                    if(epsilon_1 == 0)
                        epsilon_1 = 1;
                }
                printf("Ma = %.1f, epsilon = %lf, cd = %lf\n",ma[ck], epsilon, ak[ck][0]);
            }
            exp_batch.clear();  //清空上一批样本
        }
    }
    return 0;
}

vector<double> PARAMETER_IDENTIFICATION::chapman_kirk(const vector<vector<double>> exp_batch, const int ck)  //C-K法
{
    vector<double> B;
    int exp_l = exp_batch.size();
    {
        vector<vector<double>> exp = exp_batch, cal(exp_l, vector<double>(y_num, 0));
        vector<double> Y_missile = exp_batch[0];
        for(int j=0; j<exp_l; j++)
        {
            cal[j] = Y_missile;
            vector<double> dy = rk4(Y_missile);
            Y_missile = add(Y_missile, dy);
        }

        vector<vector<vector<double>>> P(exp_l, vector<vector<double>>(y_num, vector<double>(k_num, 0)));
        vector<vector<double>> dP(y_num, vector<double>(k_num, 0));  //维度1：观测量; 维度2：参数量
        for(int i=0; i<exp_l; i++)
        {
            vector<double> a = ak[ck];
            double v = cal[i][0], theta = cal[i][1], psi = cal[i][2], y = cal[i][4];
            const ATMOSPHERE atmosphere = atmosStd(y);
            double Ma = v / atmosphere.soundsp;  //计算样本的马赫数    0~5:v,theta,psi,x,y,z
            double pvs = atmosphere.density * Area * v;  //动压一阶导数因子
            double rho_dot = 0.5 * atmosphere.density_dot * Area * v * v;

            //求(P_jk)_i
            P[i] = dP;

            //维度1：观测量; 维度2：参数量
            /**********dv/da**********/
            dP[0][0] -= ((pvs*dP[0][0]*a[0] + rho_dot*dP[4][0]*a[0] + m0*g*dP[1][0]*cos(theta) + 0.5*pvs*v)/m0) * h;
            dP[0][1] -= ((pvs*dP[0][1]*a[0] + rho_dot*dP[4][1]*a[0] + m0*g*dP[1][1]*cos(theta))/m0) * h;
            
            /**********dtheta/da**********/
            dP[1][0] += (0.5*pvs/v*dP[0][0]*a[1]/m0 + rho_dot/v*dP[4][0]*a[1]/m0 + g*sin(theta)/v*dP[1][0] + g*cos(theta)*dP[0][0]/v/v) * h;
            dP[1][1] += (0.5*pvs/v*dP[0][1]*a[1]/m0 + rho_dot/v*dP[4][1]*a[1]/m0 + g*sin(theta)/v*dP[1][1] + g*cos(theta)*dP[0][1]/v/v + 0.5*pvs/m0) * h;
            
            /**********dpsi/da**********/
            dP[2][0] -= (0.5*pvs/v*dP[0][0]*a[1]/m0/cos(theta) + rho_dot/v*dP[4][0]*a[1]/m0/cos(theta) + 0.5*pvs*sin(theta)*dP[1][0]*a[1]/m0/cos(theta)/cos(theta)) * h;
            dP[2][1] -= (0.5*pvs/v*dP[0][1]*a[1]/m0/cos(theta) + rho_dot/v*dP[4][0]*a[1]/m0/cos(theta) + 0.5*pvs*sin(theta)*dP[1][1]*a[1]/m0/cos(theta)/cos(theta) + 0.5*pvs/m0/cos(theta)) * h;
            
            /**********dx/da**********/
            dP[3][0] += (dP[0][0]*cos(theta)*cos(psi) - dP[1][0]*sin(theta)*cos(psi)*v - dP[2][0]*cos(theta)*sin(psi)*v) * h;
            dP[3][1] += (dP[0][1]*cos(theta)*cos(psi) - dP[1][1]*sin(theta)*cos(psi)*v - dP[2][1]*cos(theta)*sin(psi)*v) * h;

            /**********dy/da**********/
            dP[4][0] += (dP[0][0]*sin(theta) + dP[1][0]*v*cos(theta)) * h;
            dP[4][1] += (dP[0][1]*sin(theta) + dP[1][1]*v*cos(theta)) * h;

            /**********dz/da**********/
            dP[5][0] -=(dP[0][0]*cos(theta)*sin(psi) - dP[1][0]*sin(theta)*sin(psi)*v + dP[2][0]*cos(theta)*cos(psi)*v) * h;
            dP[5][1] -=(dP[0][1]*cos(theta)*sin(psi) - dP[1][1]*sin(theta)*sin(psi)*v + dP[2][1]*cos(theta)*cos(psi)*v) * h;
        }


        vector<vector<double>> A(k_num, vector<double>(k_num, 0));  //敏感系数矩阵
        vector<double> C(k_num, 0);  //误差向量
        for(int r=0; r<k_num; r++)  //遍历行row
        {
            for(int c=0; c<k_num; c++)  //遍历列column, A为对角矩阵
            {
                for(int j=0; j<y_num; j++)  //遍历j, 外层累加
                    for(int i=0; i<P.size(); i++)  //遍历i, 内层累加
                        A[r][c] += P[i][j][c] * P[i][j][r];
            }
            for(int j=0; j<y_num; j++)  //遍历j, 外层累加
                for(int i=0; i<P.size(); i++)  //遍历i, 内层累加
                    C[r] += (exp[i][j] - cal[i][j]) * P[i][j][r]; 
        }

        vector<vector<double>> A_inv = LUP_solve_inverse(A);
        if(A_inv.size()<=0)
            return vector<double>(0);

        B = mul(A_inv, C);
        for(int i=0; i<k_num ;i++)
            ak[ck][i] += B[i];
    }
    return B;
}

double PARAMETER_IDENTIFICATION::max_likelihood_estimate(const vector<vector<double>> exp_batch, const int ck)  //最大似然估计法
{
    double J = 0;  //目标函数
    int exp_l = exp_batch.size();
    {
        vector<vector<double>> exp = exp_batch, cal(exp_l, vector<double>(y_num, 0));
        vector<double> Y_missile = exp_batch[0];
        for(int j=0; j<exp_l; j++)
        {
            cal[j] = Y_missile;
            vector<double> dy = rk4(Y_missile);
            Y_missile = add(Y_missile, dy);
        }

        vector<vector<vector<double>>> P(exp_l, vector<vector<double>>(y_num, vector<double>(k_num, 0)));
        vector<vector<double>> dP(y_num, vector<double>(k_num, 0));  //维度1：观测量; 维度2：参数量
        for(int i=0; i<exp_l; i++)
        {
            vector<double> a = ak[ck];
            double v = cal[i][0], theta = cal[i][1], psi = cal[i][2], y = cal[i][4];
            const ATMOSPHERE atmosphere = atmosStd(y);
            double Ma = v / atmosphere.soundsp;  //计算样本的马赫数    0~5:v,theta,psi,x,y,z
            double pvs = atmosphere.density * Area * v;  //动压一阶导数因子
            double rho_dot = 0.5 * atmosphere.density_dot * Area * v * v;

            //求(P_jk)_i
            P[i] = dP;

            //维度1：观测量; 维度2：参数量
            /**********dv/da**********/
            dP[0][0] -= ((pvs*dP[0][0]*a[0] + rho_dot*dP[4][0]*a[0] + m0*g*dP[1][0]*cos(theta) + 0.5*pvs*v)/m0) * h;
            dP[0][1] -= ((pvs*dP[0][1]*a[0] + rho_dot*dP[4][1]*a[0] + m0*g*dP[1][1]*cos(theta))/m0) * h;
            
            /**********dtheta/da**********/
            dP[1][0] += (0.5*pvs/v*dP[0][0]*a[1]/m0 + rho_dot/v*dP[4][0]*a[1]/m0 + g*sin(theta)/v*dP[1][0] + g*cos(theta)*dP[0][0]/v/v) * h;
            dP[1][1] += (0.5*pvs/v*dP[0][1]*a[1]/m0 + rho_dot/v*dP[4][1]*a[1]/m0 + g*sin(theta)/v*dP[1][1] + g*cos(theta)*dP[0][1]/v/v + 0.5*pvs/m0) * h;
            
            /**********dpsi/da**********/
            dP[2][0] -= (0.5*pvs/v*dP[0][0]*a[1]/m0/cos(theta) + rho_dot/v*dP[4][0]*a[1]/m0/cos(theta) + 0.5*pvs*sin(theta)*dP[1][0]*a[1]/m0/cos(theta)/cos(theta)) * h;
            dP[2][1] -= (0.5*pvs/v*dP[0][1]*a[1]/m0/cos(theta) + rho_dot/v*dP[4][0]*a[1]/m0/cos(theta) + 0.5*pvs*sin(theta)*dP[1][1]*a[1]/m0/cos(theta)/cos(theta) + 0.5*pvs/m0/cos(theta)) * h;
            
            /**********dx/da**********/
            dP[3][0] += (dP[0][0]*cos(theta)*cos(psi) - dP[1][0]*sin(theta)*cos(psi)*v - dP[2][0]*cos(theta)*sin(psi)*v) * h;
            dP[3][1] += (dP[0][1]*cos(theta)*cos(psi) - dP[1][1]*sin(theta)*cos(psi)*v - dP[2][1]*cos(theta)*sin(psi)*v) * h;

            /**********dy/da**********/
            dP[4][0] += (dP[0][0]*sin(theta) + dP[1][0]*v*cos(theta)) * h;
            dP[4][1] += (dP[0][1]*sin(theta) + dP[1][1]*v*cos(theta)) * h;

            /**********dz/da**********/
            dP[5][0] -=(dP[0][0]*cos(theta)*sin(psi) - dP[1][0]*sin(theta)*sin(psi)*v + dP[2][0]*cos(theta)*cos(psi)*v) * h;
            dP[5][1] -=(dP[0][1]*cos(theta)*sin(psi) - dP[1][1]*sin(theta)*sin(psi)*v + dP[2][1]*cos(theta)*cos(psi)*v) * h;
        }
    

        vector<vector<vector<double>>> V(exp_l, vector<vector<double>>(y_num, vector<double>(1, 0)));  //最大似然估计-新息, y_num*1列向量
        vector<vector<double>> R(y_num, vector<double>(y_num, 0));  //最大似然估计-协方差矩阵
        for(int i=0; i<exp_l; i++)
        {
            for(int j=0; j<y_num; j++)
                V[i][j][0] = exp[i][j] - cal[i][j];  //计算新息
        }

        R = eye(y_num, 1.);  //使用单位矩阵作为测量误差的协方差矩阵, 各参数均方差为1, 且相互独立(协方差=0)
        auto R_inv = LUP_solve_inverse(R);
        if(R_inv.size()<=0)
            return -1;

        vector<vector<double>> nable_J(1, vector<double>(k_num, 0));  //目标函数的一阶导数
        vector<vector<double>> nable2_J(k_num, vector<double>(k_num, 0));  //目标函数的二阶导数

        for(int i=0; i<exp_l; i++)
        {
            J += mul(mul(transpose(V[i]), R_inv), V[i])[0][0] + log(absdet(R));
            nable_J = add(nable_J, mul(mul(transpose(V[i]), R_inv), P[i]));
            nable2_J = add(nable2_J, mul(mul(transpose(P[i]), R_inv), P[i]));
        }
        nable_J = mul(-2./exp_l, nable_J);
        nable2_J = mul(2./exp_l, nable2_J);

        auto nable2_J_inv = LUP_solve_inverse(nable2_J);
        if(nable2_J_inv.size()<=0)
            return -1;

        vector<double> B = mul(nable2_J_inv, nable_J[0]);
        for(int i=0; i<k_num ;i++)
            ak[ck][i] -= B[i];
    }
    J /= exp_l;
    return J;
}

static ATMOSPHERE atmosStd(double height)  //获取气象参数
{
    double TH, PH;
    double T_dot;
    double height_V = height+1400;
    ATMOSPHERE atmosphere;

    if (height_V <= 11000.0)
    {
        TH = 288.15 - 0.0065*height_V;
        PH = 101325 * pow((TH / 288.15), 5.256);
        T_dot = (101325 * 5.256 * pow((TH / 288.15), 4.256) / 288.15 * -0.0065) / (TH * 287.053) + \
                0.0065 * 287.053 * PH / (TH*287.053) / (TH*287.053);
    }
    else if (height_V > 11000.0 && height_V <= 20000.0)
    {
        TH = 216.65;
        PH = 22632 * exp((11000 - height_V) / 6341.61557);
        T_dot = 22632 / -6341.61557 * exp((11000 - height_V) / 6341.61557) / (TH*287.053);
    }
    else if (height_V > 20000.0 && height_V <= 32000.0)
    {
        TH = 216.65 + 0.001*(height_V - 20000);
        PH = 5474.8*pow((TH / 216.65), -34.1632);
        T_dot = (5474.8 * -34.1632 * pow((TH / 216.65), -35.1632) / 216.65 * 0.001) / (TH * 287.053) - \
                0.001 * 287.053 * PH / (TH*287.053) / (TH*287.053);
    }
    else if (height_V > 32000.0 && height_V <= 47000.0)
    {
        TH = 228.65 + 0.0028*(height_V - 32000);
        PH = 868.01*pow((TH / 228.65), -12.20115);
        T_dot = (868.01 * -12.20115 * pow((TH / 228.65), -13.20115) / 228.65 * 0.0028) / (TH * 287.053) - \
                0.0028 * 287.053 * PH / (TH*287.053) / (TH*287.053);
    }
    else
    {
        TH = 270.65;
        PH = 110.9*exp((47000 - height_V) / 7922.2626);
        T_dot = 110.9 / -7922.2626 * exp((47000 - height_V) / 7922.2626) / (TH*287.053);
    }
    atmosphere.density = PH / (TH*287.053);
    atmosphere.soundsp = sqrt(1.4*287.053*TH);
    atmosphere.density_dot = T_dot;

    return atmosphere;
}

static double kalman_filter(const double z_measure, double &z_last, double &pz_last)  //卡尔曼滤波
{
    const double KALMAN_Q = 0.1;  //过程噪声
    const double KALMAN_R = 1.0;  //测量噪声

    double z_mid = z_last;  //z_last=x(k-1|k-1),z_mid=x(k|k-1)
    double z_now;

    double pz_mid = pz_last + KALMAN_Q;  //pz_mid=p(k|k-1),pz_last=p(k-1|k-1),Q=噪声
    double pz_now;

    double kg;

    /*卡尔曼滤波的五个公式*/
    kg=pz_mid/(pz_mid+KALMAN_R);  //kg为kalman filter，R 为噪声
    z_now=z_mid+kg*(z_measure-z_mid);  //估计出的最优值
    pz_now=(1-kg)*pz_mid;  //最优值对应的covariance
    pz_last = pz_now;  //更新covariance 值
    z_last = z_now;  //更新系统状态值

    // return z_now;  //启用滤波：辨识结果偏差很大
    return z_measure;  //禁用滤波
}
