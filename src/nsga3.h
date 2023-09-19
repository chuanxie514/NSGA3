#ifndef NSGA3_H_
#define NSGA3_H_

#include "global.h"
#include "individual.h"
#include "refpoint.h"
#include "recombination.h"
#include "common.h"
#include <iomanip> 

class TNSGA3 {

public:
	TNSGA3();
	virtual ~TNSGA3();

	int popsize;					// 种群大小
	vector<TIndividual> parent_pop;	// 父代种群
	vector<TIndividual> child_pop;	// 子代种群
	vector<TIndividual> mixed_pop;	// 组合种群
	vector<TRefpoint> refepoints;	// 参考点集合

	void init_population();			 	// 初始化种群
	void init_refpoints2D(); 			// 初始化二维参考点集合
	void init_refpointsND();			// 生成多维参考点集合
	void gen_refpointsND(vector<TRefpoint>* rps, TRefpoint* rp, size_t left, size_t count);
	void merge_populations();			// 合并父代和子代种群

	void update_idealpoint(TIndividual& indv);
	void run(size_t ngen, size_t run);
	void evolution();
	TParetofront::TPfronts nondominatedsort();

	void normalize(TParetofront::TPfronts Fronts);
	void associate(TParetofront::TPfronts Fronts);
	void niching();

	void save_front(char file_name[512]);
};

TNSGA3::TNSGA3() {
	popsize = 0;
	idealpoint = new double[numObjectives];
	nadirpoint = new double[numObjectives];
	for (int i = 0; i < numObjectives; i++) {
		idealpoint[i] = MAX_DOUBLE;
		nadirpoint[i] = MIN_DOUBLE;
	}
}

TNSGA3::~TNSGA3() {
	delete[] idealpoint;
}

// 为 2 个物镜初始化一组均匀分布的权重向量
void TNSGA3::init_refpoints2D()
{
	for (size_t i = 0; i <= nrefpoints_2D; i++) {
		TRefpoint rpoint;
		vector<size_t> array;
		array.push_back(i);
		array.push_back(nrefpoints_2D - i);
		for (size_t j = 0; j < array.size(); j++)
			rpoint.refpoint.push_back(1.0 * array[j] / nrefpoints_2D);
		refepoints.push_back(rpoint);
	}
}
// 为 大于2 个目标初始化一组均匀分布的权重向量
void TNSGA3::init_refpointsND()
{
	TRefpoint rpoint(numObjectives);
	gen_refpointsND(&refepoints, &rpoint, p_boundary, 0);

	if (p_inside > 0) {
		vector<TRefpoint> inside_rps;	// inside reference points
		gen_refpointsND(&refepoints, &rpoint, p_inside, 0);

		double center = 1.0 / numObjectives;
		for (size_t i = 0; i < inside_rps.size(); ++i) {
			for (size_t j = 0; j < numObjectives; ++j)
				inside_rps[i].refpoint[j] = (center + inside_rps[i].refpoint[j]) / 2;
			refepoints.push_back(inside_rps[i]);
		}
	}
}

// 为 <2 个生成一组均匀分布的权重向量
void TNSGA3::gen_refpointsND(vector<TRefpoint>* rps, TRefpoint* rp, size_t left, size_t count)
{
	if (count == numObjectives - 1) {
		rp->refpoint[count] = 1.0 * left / p_boundary;
		rps->push_back(*rp);
	}
	else {
		for (size_t i = 0; i <= left; ++i) {
			rp->refpoint[count] = 1.0 * i / p_boundary;
			gen_refpointsND(rps, rp, left - i, count + 1);
		}
	}
}

void TNSGA3::init_population()
{
	for (size_t i = 0; i < popsize; ++i) {
		TIndividual indv;
		indv.init_individual();
		indv.eval_objective();
		parent_pop.push_back(indv);
	}
}

void TNSGA3::merge_populations()
{
	for (size_t i = 0; i < popsize; ++i)
		mixed_pop.push_back(parent_pop[i]);

	for (size_t i = 0, j = popsize; i < popsize; ++i, j++)
		mixed_pop.push_back(child_pop[j]);

}

TParetofront::TPfronts TNSGA3::nondominatedsort()
{
	int rank = 1, count = 0;
	vector<int> ranks(mixed_pop.size(), 0);
	TParetofront::TPfronts fronts;

	while (count < mixed_pop.size()) {
		TParetofront::TPFront cfront;
		for (size_t i = 0; i < mixed_pop.size(); ++i) {
			if (ranks[i] > 0)
				continue;
			bool isdominated = false;
			for (size_t j = 0; j < cfront.size(); ++j) {
				if (mixed_pop[cfront[j]] < mixed_pop[i]) {
					isdominated = true;
					break;
				}
				else if (mixed_pop[i] < mixed_pop[cfront[j]]) {
					cfront.erase(cfront.begin() + j);
					j = j - 1;
				}
			}
			if (!isdominated) cfront.push_back(i);
		}
		for (size_t i = 0; i < cfront.size(); ++i) ranks[cfront[i]] = rank;
		fronts.push_back(cfront);
		count += cfront.size();
		rank++;
	}

	// 识别最后一个前端索引 （FIDX）
	int fidx = 0, block = 0;
	while (block < popsize) block += fronts[fidx++].size();

	// 删除无用的fronts
	fronts.erase(fronts.begin() + fidx, fronts.end());

	return fronts;
}

void TNSGA3::normalize(TParetofront::TPfronts Fronts) {
	// 计算理想点
	for (size_t i = 0; i < Fronts[0].size(); ++i) {
		size_t idx = Fronts[0][i];
		for (size_t j = 0; j < numObjectives; ++j) {
			if (mixed_pop[idx].y_obj[j] < idealpoint[j])
				idealpoint[j] = mixed_pop[idx].y_obj[j];
		}
	}
	// 计算最低点
	for (size_t i = 0; i < Fronts.size(); ++i) {
		for (size_t j = 0; j < Fronts[i].size(); ++j) {
			size_t idx = Fronts[i][j];
			for (size_t k = 0; k < numObjectives; ++k) {
				if (mixed_pop[idx].y_obj[k] > nadirpoint[k])
					nadirpoint[k] = mixed_pop[idx].y_obj[k];
			}
		}
	}
	// 归一化化目标
	for (size_t i = 0; i < Fronts.size(); ++i)
		for (size_t j = 0; j < Fronts[i].size(); ++j)
			mixed_pop[Fronts[i][j]].norm_objective();
}

void TNSGA3::associate(TParetofront::TPfronts Fronts)
{
	for (size_t i = 0; i < refepoints.size(); ++i)
		refepoints[i].clear_refpoint();

	for (size_t i = 0; i < Fronts.size(); ++i) {
		for (size_t j = 0; j < Fronts[i].size(); ++j) {

			double mindist = MAX_DOUBLE;
			size_t idx = refepoints.size();
			for (size_t r = 0; r < refepoints.size(); ++r) {
				double pdist = perpendicular_distance(refepoints[r].refpoint,
					mixed_pop[Fronts[i][j]].n_obj);
				if (pdist < mindist) {
					mindist = pdist;
					idx = r;
				}
			}
			if (i + 1 == Fronts.size())
				refepoints[idx].add_member(Fronts[i][j], mindist);
			else
				refepoints[idx].sum_member();
		}
	}
}

void TNSGA3::niching()
{
	vector<bool> rp_isable(refepoints.size(), true);
	while (parent_pop.size() < popsize) {

		size_t rp_idx = find_refpoint(refepoints, rp_isable);
		int ind_idx = select_member(refepoints[rp_idx]);

		if (ind_idx < 0) rp_isable[rp_idx] = false;
		else {
			refepoints[rp_idx].sum_member();
			refepoints[rp_idx].remove_member(ind_idx);
			parent_pop.push_back(mixed_pop[ind_idx]);
		}
	}
}

void TNSGA3::evolution()
{
	child_pop.clear(); // 清理子代种群
	for (size_t i = 0; i < popsize; ++i) {
		int p1 = int(popsize * rnd_uni(&rnd_uni_init));
		int p2 = int(popsize * rnd_uni(&rnd_uni_init));

		TIndividual child1, child2;
		realbinarycrossover(parent_pop[p1], parent_pop[p2], child1, child2);
		realmutation(child1, 1.0 / numVariables);
		realmutation(child2, 1.0 / numVariables);

		child1.eval_objective();
		child2.eval_objective();

		child_pop.push_back(child1);
		child_pop.push_back(child2);
	}

	// mixed_pop: 父子代种群组合
	mixed_pop.clear();
	merge_populations();

	// 非支配排序算法
	TParetofront::TPfronts Fronts = nondominatedsort();

	// 将前 F-1 个前沿的个体复制到下一代种群中
	parent_pop.clear();
	for (int i = 0; i < Fronts.size() - 1; ++i)
		for (int j = 0; j < Fronts[i].size(); ++j)
			parent_pop.push_back(mixed_pop[Fronts[i][j]]);

	// 如果种群已经完整，则结束
	if (parent_pop.size() == popsize) return;

	// 对目标进行归一化（算法2：normalize 函数）
	normalize(Fronts);

	// 成员关联（算法3：associate 函数）
	associate(Fronts);

	// 选择 k 个个体（算法4：niching 函数）
	niching();
}

void TNSGA3::run(size_t ngen, size_t run) {

	if (numObjectives == 2) init_refpoints2D();
	else init_refpointsND();

	popsize = refepoints.size();//将种群规模设置为与参考点数量相等
	while (popsize % 4) popsize++;	// 人口规模必须是 4 的倍数

	init_population();
	for (size_t g = 0; g < ngen; ++g) evolution();

	char file_name[512];
	sprintf_s(file_name, "ParetoFront/NSGA3_%s_%dD_R%lu.data", strTestInstance, numObjectives, run);
	save_front(file_name);

	parent_pop.clear();
}

void TNSGA3::save_front(char file_name[512])
{
	std::fstream fout;
	fout.open(file_name, std::ios::out);

	// 设置输出精度为小数点后四位
	fout << std::fixed << std::setprecision(4);

	for (size_t i = 0; i < parent_pop.size(); i++) {
		for (int j = 0; j < numObjectives; j++) {
			fout << parent_pop[i].y_obj[j] << "\t";
			if (j == numObjectives - 1) {
				fout << "\n";
			}
		}
	}

	fout.close();
}



#endif /* NSGA3_H_ */
