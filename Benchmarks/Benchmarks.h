#ifndef BENCHMARKS_H
#define BENCHMARKS_H

#include <string>
#include <vector>
#include "../util/json.hpp"
using json = nlohmann::json;
using std::string;
using std::vector;

class Benchmarks {
private:
	int group_num;
    json func_config;
	string funcID;
	template<typename T>T** read_data(string);
	
	// 第一类函数的变量
	int **group;
	double ***R, **xopt;
	double local_eva_for_global_solution(double*, int);
	bool overlap_grouping;
	double local_eva_type1(double* x, int groupIndex);
	// 第二类函数的变量
	double **R_global, **A, **W=nullptr;
	bool if_rotate;
	bool if_shift;
	double local_eva_type2(double* x, int groupIndex);
	double global_eva_type2(double* x);
	// 第三类函数的变量
	double local_eva_type3(double* x,int groupIndex);
	int funcClass;
	// string funcType;
	double weight;
	bool if_heterogeneous;
	// 第四类函数的变量
	double local_eva_type4(double* x,int groupIndex);
	int target_num;
	int source_num;
	int coordinate_dim;
	double **target, **source, **noisy_dis = nullptr;
	// 第五类函数的变量	
	double local_eva_type5(double* x,int groupIndex);
	bool dfs(int x);
	double KM();
	int **match_truth;
	// const int N=205;
	// double w[N][N];
	// double la[N],lb[N];
	// bool va[N],vb[N];
	// int match[N];
	// int n;
	double **w;
	double *Lx,*Ly;
	bool *VisX,*VisY;
	int *MatchY;
	double* Slack;

public:
	int max_eva_times;
	int eva_count;
	bool reach_max_eva_times;
	Benchmarks(string ID,int max_eva_times = 3000000, bool = true);
	~Benchmarks();
	double global_eva(double* x);
	double local_eva(double* x, int groupIndex);
	double getMinX();
	double getMaxX();
	int getGroupNum();
	int getDimension();
    vector<int> getOverlapDim(int g1,int g2);
    vector<int> getOverlapDimIndex(int g1,int g2);
    vector<int> getGroupDim(int g);
	vector<int> getGroupExcluDim(int groupIndex);
	vector<int> getOverlapGroup(int g);
	double getBestFitness();
	bool reachMaxEva();
	double getLocalOpt(int groupIndex);
	double** getNetworkGraph();
	double getMatchRes(double* x);
};

#endif