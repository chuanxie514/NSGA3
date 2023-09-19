#ifndef GLOBAL_H_  // 全局头文件防止重复包含
#define GLOBAL_H_

#include <iostream>   // 输入输出流
#include <fstream>    // 文件输入输出流
#include <string>     // 字符串处理
#include <stdlib.h>   // 标准库函数，如内存分配和类型转换
#include <stdio.h>    // 标准输入输出
#include <math.h>     // 数学函数
#include <time.h>     // 时间函数
#include <memory.h>   // 内存处理
#include <vector>     // 向量容器
#include <limits>     // 数值极限
#include <algorithm>  // 算法库

#define MAX_INT 	numeric_limits<size_t>::max()  // 最大整数值
#define MAX_DOUBLE 	numeric_limits<double>::max()  // 最大双精度浮点数值
#define MIN_DOUBLE 	numeric_limits<double>::min()  // 最小双精度浮点数值

using namespace std;  // 使用标准命名空间

char strTestInstance[256];  // 测试实例字符串

// 变量的范围界限
double lowBound = 0, uppBound = 1;

// 变量和目标的维度
int numVariables;  // 变量维度
int numObjectives;  // 目标维度

// 最大迭代次数
int max_gen;

// 分布索引用于SBX交叉和多项式变异
int id_cx = 30;  // 交叉
int id_mu = 20;  // 变异

// 用于标准化的理想点
double* idealpoint;

// 用于标准化的最差点
double* nadirpoint;

// 参考点数量
size_t nrefpoints_2D = 99;

// 边界和内部层的划分点数
size_t p_boundary = 0, p_inside = 0;

// 随机数生成参数
int seed = 237;
long rnd_uni_init;

#include "random.h"  // 随机数生成函数

#endif /* GLOBAL_H_ */  