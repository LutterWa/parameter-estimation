#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "algorithm.h"

int main()
{
    PARAMETER_IDENTIFICATION para(CK);
    if(para.load_data()>=0)//加载测试数据
        if(para.parameter_estimation()>=0)  //C-K法辨识模型参数
            if(para.save_data()>=0)
                return 0;
    system("pause");
    return -1;
}