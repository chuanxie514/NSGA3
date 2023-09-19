#ifndef GLOBAL_H_  // ȫ��ͷ�ļ���ֹ�ظ�����
#define GLOBAL_H_

#include <iostream>   // ���������
#include <fstream>    // �ļ����������
#include <string>     // �ַ�������
#include <stdlib.h>   // ��׼�⺯�������ڴ���������ת��
#include <stdio.h>    // ��׼�������
#include <math.h>     // ��ѧ����
#include <time.h>     // ʱ�亯��
#include <memory.h>   // �ڴ洦��
#include <vector>     // ��������
#include <limits>     // ��ֵ����
#include <algorithm>  // �㷨��

#define MAX_INT 	numeric_limits<size_t>::max()  // �������ֵ
#define MAX_DOUBLE 	numeric_limits<double>::max()  // ���˫���ȸ�����ֵ
#define MIN_DOUBLE 	numeric_limits<double>::min()  // ��С˫���ȸ�����ֵ

using namespace std;  // ʹ�ñ�׼�����ռ�

char strTestInstance[256];  // ����ʵ���ַ���

// �����ķ�Χ����
double lowBound = 0, uppBound = 1;

// ������Ŀ���ά��
int numVariables;  // ����ά��
int numObjectives;  // Ŀ��ά��

// ����������
int max_gen;

// �ֲ���������SBX����Ͷ���ʽ����
int id_cx = 30;  // ����
int id_mu = 20;  // ����

// ���ڱ�׼���������
double* idealpoint;

// ���ڱ�׼��������
double* nadirpoint;

// �ο�������
size_t nrefpoints_2D = 99;

// �߽���ڲ���Ļ��ֵ���
size_t p_boundary = 0, p_inside = 0;

// ��������ɲ���
int seed = 237;
long rnd_uni_init;

#include "random.h"  // ��������ɺ���

#endif /* GLOBAL_H_ */  