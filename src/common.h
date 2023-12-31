
#ifndef COMMON_H_
#define COMMON_H_

#include "global.h"


inline double square(double n) {
	return n * n;
}

// 计算缩放函数 （ASF）
double ASF(const vector<double> y_obj, vector<double> w)
{
	double max_ratio = -MAX_DOUBLE;
	for (size_t i = 0; i < numObjectives; ++i)
		max_ratio = max(max_ratio, y_obj[i] / w[i]);


	return max_ratio;
}

double perpendicular_distance(vector<double> rpoint, vector<double> solution)
{
	double numerator = 0, denominator = 0;
	for (int i = 0; i < rpoint.size(); i++) {
		numerator += rpoint[i] * solution[i];
		denominator += square(rpoint[i]);
	}
	double k = numerator / denominator;

	double dist = 0;
	for (int i = 0; i < rpoint.size(); i++) {
		dist += square(k * rpoint[i] - solution[i]);
	}
	return sqrt(dist);
}
//计算垂直距离的函数，它接受一个参考点向量rpoint和一个解向量solution，并返回它们之间的欧氏距离
size_t find_refpoint(vector<TRefpoint> rpoints, vector<bool> rp_isable)
{
	size_t min_size = MAX_INT;
	for (size_t i = 0; i < rpoints.size(); ++i) {
		if (rp_isable[i]) // 如果参考点在前面有成员
			min_size = min(min_size, rpoints[i].member_size());
	}

	vector<size_t> min_rpoints;
	for (size_t i = 0; i < rpoints.size(); ++i)
		if (rpoints[i].member_size() == min_size)
			min_rpoints.push_back(i);

	// 返回一个随机参考点
	return min_rpoints[rand() % min_rpoints.size()];
}

int select_member(TRefpoint rpoint)
{
	int idx = -1;
	if (rpoint.has_candidate_members()) {
		if (rpoint.member_size() == 0)
			idx = rpoint.get_nearest_member();
		else
			idx = rpoint.get_random_member();
	}
	return idx;
}

#endif /* COMMON_H_ */
