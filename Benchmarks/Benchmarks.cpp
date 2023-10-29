#include <sstream>
#include <fstream>
#include <iostream>
#include<cstring>
#include<cstdio>
#include<cmath>
#include<algorithm>
// #include<regex>
#include "../Benchmarks/Benchmarks.h"
#include "../Benchmarks/BaseFunction.h"
using namespace std;
using json = nlohmann::json;

Benchmarks::Benchmarks(string ID,int max_eva_times,bool overlap){
    
    string config_path = "../Benchmarks/default_config.json";
    string data_path = "../Benchmarks/data/";
    ifstream file_config(config_path);
    json config;
    file_config >> config;
    file_config.close();
    this->func_config = config["benchmarks"][ID];
    this->funcID = ID;
    this->group_num = func_config["group_num"];

    this->max_eva_times = max_eva_times;
    this->eva_count = 0;
    this->reach_max_eva_times = false;
    this->overlap_grouping = overlap;

    this->funcClass = 1;
    this->if_rotate = true;
    this->if_shift = false;
    this->if_heterogeneous = false;

    try{
        this->funcClass = func_config["funcClass"];
    }catch(...){
    }

    if(funcClass == 1){

        // 读取维度数据,并将取值范围由1-1000变为0-999
        group = read_data<int>(data_path+"group_"+ID+".txt");
        for (int i = 0; i < group_num; i++) {
            int size =func_config["subproblems"][i]["dimension"];
            for (int j = 0; j < size; j++) {
                group[i][j] -= 1;
            }
        }
        
        // 读取最优值
        xopt = read_data<double>(data_path+"xopt_"+ID+".txt");

        try
        {
            if_rotate = func_config["if_rotate"];
            if(if_rotate==false){
                return;
            }
        }
        catch(...)
        {
        }
        
        // 读取旋转矩阵
        R = new double**[group_num];
        double** testPtr = read_data<double>(data_path+"R"+to_string(1)+"_"+ID+".txt");
        if(testPtr == nullptr){
            // double** R100 = read_data<double>(data_path+ "R100_"+ID+".txt");
            // double** R250 = read_data<double>(data_path+ "R250_"+ID+".txt");
            // double** R500 = read_data<double>(data_path+ "R500_"+ID+".txt");
            
            // for(int i=0;i<group_num;i++){
            //     if(func_config["subproblems"][i]["dimension"] == 100)
            //         R[i] = R100;
            //     else if(func_config["subproblems"][i]["dimension"] == 250)
            //         R[i] = R250;
            //     else if(func_config["subproblems"][i]["dimension"] == 500)
            //         R[i] = R500;
            // }
            for(int i=0;i<group_num;i++){
                string d = to_string(func_config["subproblems"][i]["dimension"]);
                R[i] = read_data<double>(data_path+ "R"+d+"_"+ID+".txt");
            }
        }else{
            for(int i=0;i<group_num;i++){
                R[i] = read_data<double>(data_path+"R"+to_string(i+1)+"_"+ID+".txt");
            }
        }
    }else if (funcClass == 2){  
        // funcType = func_config["base_function"];
        A = read_data<double>(data_path+"A_"+ID);
        // W = read_data<double>(data_path+"W_"+ID);

        if_rotate = func_config["if_rotate"];
        if_heterogeneous = func_config["heterogeneous"];
        weight = func_config["weight"];
        if(if_rotate == true)
            R_global = read_data<double>(data_path+"R_"+ID);
        
        try
        {
            if_shift = func_config["if_shift"];
        }
        catch(...)
        {
        }

        if(if_shift == true)
            xopt = read_data<double>(data_path+"xopt_"+ID);
        else{
            int len = getDimension();
            xopt = new double*[1];
            xopt[0] = new double[len]{0};
        }

        getNetworkGraph();
        
    }else if (funcClass == 3){
        try
        {
            if_shift = func_config["if_shift"];
        }
        catch(...)
        {
        }
        if(if_shift == true)
            xopt = read_data<double>(data_path+"xopt_"+ID);
    }else if (funcClass == 4 || funcClass == 5){
        target_num = func_config["target_num"];
        source_num = func_config["group_num"];
        // int dimension = func_config["dimension"];
        coordinate_dim = func_config["coordinate"];
        int dimension = target_num * coordinate_dim;
        func_config["dimension"] = dimension;
        // target = read_data<double>(data_path+"target_WSN_location_1");
        // source = read_data<double>(data_path+"source_WSN_location_1");
        // noisy_dis = read_data<double>(data_path+"noisydis_WSN_location_1");
        // W=read_data<double>(data_path+"W_WSN_location_1");
        // noisy_dis = read_data<double>(data_path+"seq_dis_s"+to_string(source_num)+"_WSN_location");
        target = read_data<double>(data_path+"target_WSN_location");
        source = read_data<double>(data_path+"source_"+to_string(source_num)+"_WSN_location");
        W=read_data<double>(data_path+"W_"+to_string(source_num)+"_WSN_location");
        bool if_noisy = func_config["noisy"];
        double degree;
        if(if_noisy == true){
            degree = func_config["noisy_degree"];
        }
        noisy_dis = new double*[source_num];
        for(int i=0;i<source_num;i++){
            noisy_dis[i] = new double[target_num];
            for(int j=0;j<target_num;j++){
                if(coordinate_dim == 3)
                    noisy_dis[i][j] = sqrt(pow(target[j][0] - source[i][0],2) + pow(target[j][1]  - source[i][1],2) + pow(target[j][2] - source[i][2],2));
                else if(coordinate_dim == 2)
                    noisy_dis[i][j] = sqrt(pow(target[j][0] - source[i][0],2) + pow(target[j][1]  - source[i][1],2) );
                if(if_noisy == true)
                    noisy_dis[i][j] += (rand()*1.0/RAND_MAX-0.5) * degree;
            }
        }
        if (funcClass == 5){
            match_truth = new int*[source_num];      
            for(int i=0;i<source_num;i++){
                double* origin_dis = new double[target_num];
                for(int j=0;j<target_num;j++){
                    origin_dis[j] = noisy_dis[i][j];
                }
                sort(noisy_dis[i],noisy_dis[i]+target_num);
                match_truth[i] = new int[target_num];
                for(int j=0;j<target_num;j++){
                    for(int k=0;k<target_num;k++){
                        if(noisy_dis[i][j] == origin_dis[k]){
                            match_truth[i][j] = k;
                        }
                    }
                }
            }
            
            Lx = new double[target_num];
            Ly = new double[target_num];
            VisX = new bool[target_num];
            VisY = new bool[target_num];
            MatchY = new int[target_num];
            Slack = new double[target_num];
            w = new double*[target_num];
            for (size_t t = 0; t < target_num; t++){
                w[t] = new double[target_num];
            }
        }
    }
    // else if (funcClass == 5){
    //     // 读取目标的位置，读取传感器的位置
    //     target = read_data<double>(data_path+"target_"+ID);
    //     source = read_data<double>(data_path+"source_"+ID);
    //     noisy_dis = read_data<double>(data_path+"noisydis_"+ID);
    //     target_num = func_config["target_num"];
    //     source_num = func_config["group_num"];
        
    //     Lx = new double[target_num];
    //     Ly = new double[target_num];
    //     VisX = new bool[target_num];
    //     VisY = new bool[target_num];
    //     MatchY = new int[target_num];
    //     Slack = new double[target_num];
    //     w = new double*[target_num];
    //     for (size_t t = 0; t < target_num; t++){
    //         w[t] = new double[target_num];
    //     }
    // }
}

Benchmarks::~Benchmarks() {
    if(funcClass == 1){
        vector<double**> trash;
        for (int i = 0; i < group_num; i++) {
            delete[] group[i];
            delete[] xopt[i];
            int size = func_config["subproblems"][i]["dimension"];
            if(if_rotate == true){
                if(find(trash.begin(),trash.end(),R[i]) == trash.end()){
                    for(int j=0;j<size;j++){
                        delete[] R[i][j];
                    }
                    delete[] R[i];
                    trash.push_back(R[i]);
                }
            }
        }
        delete[] group;
        delete[] xopt;
        
        if(if_rotate == true)
            delete[] R;
    }else if (funcClass == 2){
        for (int i = 0; i < group_num; i++) {
            delete[] A[i];
            if(if_rotate == true)
                delete[] R_global[i];
        }
        delete[] A;
        if(if_rotate == true)
            delete[] R_global;
        if(if_shift == true){
            delete[] xopt[0];
            delete[] xopt;
        }
    }else if (funcClass == 3){
        if(if_shift == true){
            for (int i = 0; i < group_num; i++) {
                delete[] xopt[i];
            }
            delete[] xopt;
        }
    }else if (funcClass == 4 || funcClass == 5){        
        for (int i = 0; i < source_num; i++) {
            delete[] source[i];
            delete[] noisy_dis[i];
        }
        delete[] source;
        delete[] noisy_dis;
        
        for (int i = 0; i < target_num; i++) {
            delete[] target[i];
        }
        delete[] target;
    }
}

bool Benchmarks::reachMaxEva(){
    if(eva_count>=max_eva_times){
        if(!reach_max_eva_times){
            // cout<<"The time of evaluation has reached the maximum bound. Later evaluation results would not be recorded.\n"; 
            reach_max_eva_times = true;
        }
        return true;
    }
    return false;
}
double Benchmarks::local_eva(double* x, int groupIndex) {
    if(funcClass == 1){
        return local_eva_type1(x,groupIndex);
    }else if(funcClass == 2){
        return local_eva_type2(x,groupIndex);
    }
    else if(funcClass == 3){
        return local_eva_type3(x,groupIndex);
    }
    else if(funcClass == 4){
        return local_eva_type4(x,groupIndex);
    }
    else if(funcClass == 5){
        return local_eva_type5(x,groupIndex);
    }
    return 0;
}

double Benchmarks::local_eva_type4(double* x,int groupIndex){
    double res = 0;
    for (size_t t = 0; t < target_num; t++)
    {
        double est_dis;
        if(coordinate_dim == 3)
            est_dis = sqrt(pow(x[t*3] - source[groupIndex][0],2) + pow(x[t*3+1] - source[groupIndex][1],2) + pow(x[t*3+2] - source[groupIndex][2],2));
        else if(coordinate_dim == 2)
            est_dis = sqrt(pow(x[t*2] - source[groupIndex][0],2) + pow(x[t*2+1] - source[groupIndex][1],2) );
        res += pow(noisy_dis[groupIndex][t] - est_dis,2);
    }
    // res = res/target_num;
    eva_count += 1;
    return res;
}

double Benchmarks::local_eva_type5(double* x,int groupIndex){
    double res = 0;
    double tar_dis[target_num];
    // cout<<"local_eva"<<endl;
    for (size_t t = 0; t < target_num; t++)
    {
        if (coordinate_dim == 3)
            tar_dis[t] = sqrt(pow(x[t*3] - source[groupIndex][0],2) + pow(x[t*3+1] - source[groupIndex][1],2) + pow(x[t*3+2] - source[groupIndex][2],2));
        else if (coordinate_dim == 2)
            tar_dis[t] = sqrt(pow(x[t*2] - source[groupIndex][0],2) + pow(x[t*2+1] - source[groupIndex][1],2));
    }
    // sort(tar_dis,tar_dis+target_num);
    // for (size_t t = 0; t < target_num; t++){
    //     cout<<tar_dis[t]<<" "<<noisy_dis[groupIndex][t]<<endl;
    // }
    // double mtx[target_num][target_num] = {0};
    for (size_t t = 0; t < target_num; t++){
        for (size_t d = 0;  d< target_num; d++){
            w[t][d] = -((tar_dis[t] - noisy_dis[groupIndex][d])*(tar_dis[t] - noisy_dis[groupIndex][d]));
            // cout<<w[t][d]<<" ";
        }
        // cout<<endl;
    }
    res = -1*KM();
    // res = res/target_num;
    // 很迷，好像是加上这句之后就会影响算法效果。。。不太懂原因，所以我在draw-table里面做了这个除法
    eva_count += 1;
    return res;
}

// bool Benchmarks::dfs(int x) {
//     va[x]=1;
//     for(int y=0;y<target_num;y++)
//         if(!vb[y]){
//             // cout<<la[x]+lb[y]-w[x][y]<<endl;
//             if(fabs(la[x]+lb[y]-w[x][y])<0.001) {
//             // if(la[x]+lb[y]-w[x][y] == 0) {    
//                 vb[y]=1;
//                 if(match[y]<0||dfs(match[y])) {
//                     match[y]=x;
//                     return true;
//                 }
//             }
//             else{
//                 delta=min(delta,la[x]+lb[y]-w[x][y]);
//             }
//         }
//     return false;
// }

// double Benchmarks::KM() {
//     cout<<"KM"<<endl;
//     // memset(match,0,sizeof(match));
//     for(int i=0;i<target_num;i++) {
//         match[i] = -1;
//         la[i]=-(1<<30);
//         lb[i]=0;
//         for(int j=0;j<target_num;j++)
//             la[i]=max(la[i],w[i][j]);
//     }
//     for(int i=0;i<target_num;i++){
//         cout<<"KM"<<i<<endl;
//         while(true) {
//             memset(va,0,sizeof(va));
//             memset(vb,0,sizeof(vb));
//             delta=1<<30;
//             if(dfs(i))
//                 break;
//             for(int j=0;j<target_num;j++) {
//                 if(va[j])
//                     la[j]-=delta;
//                 if(vb[j])
//                     lb[j]+=delta;
//             }
//         }
//     }
//     cout<<"KM3"<<endl;
//     double ans=0;
//     for(int i=0;i<target_num;i++){
//         ans+=w[match[i]][i];
//         // cout<< i <<" "<<match[i]<<endl;
//     }
//     return ans;
// }

bool Benchmarks::dfs(int x)
{
    VisX[x] = true;
    for (int y = 0; y < target_num; y++)
    {
        if (VisY[y])
            continue;
        double t = Lx[x] + Ly[y] - w[x][y];
        // if (t == 0)
        if (fabs(t)<0.001)
        {
            VisY[y] = true;
            if (MatchY[y] == -1 || dfs(MatchY[y]))
            {
                MatchY[y] = x;
                return true;        //找到增广轨
            }
        }
        else if (Slack[y] > t)
            Slack[y] = t;
    }
    return false;                   //没有找到增广轨（说明顶点x没有对应的匹配，与完备匹配(相等子图的完备匹配)不符）
}

double Benchmarks::KM()                //返回最优匹配的值
{    
    // int INF = (1 << 31) - 1;
    double INF;
    int i, j;
    for(int i=0;i<target_num;i++){
        MatchY[i] = -1;
        Ly[i] = 0;
    }
    for (i = 0; i < target_num; i++)
        for (j = 0, Lx[i] = -INF; j < target_num; j++)
            if (w[i][j] > Lx[i])
                Lx[i] = w[i][j];
    for (int x = 0; x < target_num; x++)
    {
        // cout<<"km "<<x<<endl;
        for (i = 0; i < target_num; i++)
            Slack[i] = INF;
        while (true)
        {
            // memset(VisX, 0, sizeof(VisX));
            // memset(VisY, 0, sizeof(VisY));
            for(i=0;i<target_num;i++){
                VisX[i] = 0;
                VisY[i] = 0;
            }
            if (dfs(x))                     //找到增广轨，退出
                break;
            double d = INF;
            for (i = 0; i < target_num; i++)          //没找到，对l做调整(这会增加相等子图的边)，重新找
            {
                if (!VisY[i] && d > Slack[i])
                    d = Slack[i];
            }
            for (i = 0; i < target_num; i++)
            {
                if (VisX[i])
                    Lx[i] -= d;
            }
            for (i = 0; i < target_num; i++)
            {
                if (VisY[i])
                    Ly[i] += d;
                else
                    Slack[i] -= d;
            }
        }
    }
    double result = 0;
    for (i = 0; i < target_num; i++){
        if (MatchY[i] > -1)
            result += w[MatchY[i]][i];
    }
    return result;
}



double Benchmarks::getMatchRes(double* x)                //返回最优匹配的值
{    
    // for(int i=0;i<group_num;i++){
    //     for(int j=0;j<target_num;j++){
    //         cout<<match_truth[i][j]<<" ";
    //     }
    //     cout<<"match truth"<<endl;
    // }
    int ref_idx = rand()%source_num;
    int* match_ref = new int[target_num];
    local_eva(x,ref_idx);
    for(int j=0;j<target_num;j++){
        match_ref[MatchY[j]] = j;
    }
    int score = 0;
    for(int i=0;i<group_num;i++){
        double r = local_eva(x,i);
        // for(int j=0;j<target_num;j++){
        //     cout<<MatchY[j]<<" ";
        // }
        for(int j=0;j<target_num;j++){
            int measure_idx = 0;
            for(;measure_idx<target_num;measure_idx++){
                if(MatchY[measure_idx] == j)
                    break;
            }
            int origin_idx = match_truth[i][measure_idx];
            if(origin_idx == match_truth[ref_idx][match_ref[j]])
                score ++;
        }
        // cout<<score<<endl;
    }
    return score*1.0/group_num/target_num;
}

double Benchmarks::local_eva_type3(double* x,int groupIndex){
    
    if (groupIndex < 0 || groupIndex >= group_num) {
		cout << "groupIndex error\n";
		return 0;
	}
    int len = getDimension();
    
    double* shftx = new double[len];
    for (int j = 0; j < len; j++) {
        shftx[j] = x[j] - xopt[groupIndex][j];
    }
    string funcType=func_config["base_function"];

    double res = ellipsoid(shftx, len);
       
    delete[] shftx;
    eva_count += 1;
        
    return res;
}

double Benchmarks::local_eva_type2(double* x, int groupIndex) {
    
    if (groupIndex < 0 || groupIndex >= group_num) {
		cout << "groupIndex error\n";
		return 0;
	}
    // int len = func_config["subproblems"][groupIndex]["dimension"];
    int len = getDimension();
    // double ub = func_config["upper_bound"];
    // double lb = func_config["lower_bound"];
   
    
    double* shftx = new double[len];
    for (int j = 0; j < len; j++) {
        shftx[j] = x[j] - xopt[0][j];
    }

    double *input_x ;    
    if(if_rotate == true){
        input_x = multiply(shftx, R_global, len);
    }else{
        input_x = new double[len];
        memcpy(input_x,shftx,len*sizeof(double));
    }        

    string funcType;
    if(if_heterogeneous){
        funcType = func_config["base_function"][groupIndex%2];
    }else{
        funcType=func_config["base_function"];
    }

    double res = 0;
    if(funcType=="elliptic"){
        res = elliptic(input_x, len);
    }
    else if(funcType=="rastrigin"){
        res = rastrigin(input_x, len);
    }
    else if(funcType=="schwefel"){
        res = schwefel(input_x, len);
    }
    else if(funcType=="ackley"){
        res = ackley(input_x, len);
    }
    else if(funcType=="rosenbrock"){
        res = rosenbrock(input_x, len);
    }
    else if(funcType=="griewank"){
        res = griewank(input_x, len);
    }
    else if(funcType=="ellipsoid"){
        res = ellipsoid(input_x, len);
    }
    
    for(int j=0;j<len;j++){
        res += A[groupIndex][j] * input_x[j] * weight;
    }    
	delete[] input_x;
    delete[] shftx;
    
    eva_count += 1;
        
    return res;
}

double Benchmarks::global_eva_type2(double* x) {
    
    // int len = func_config["subproblems"][groupIndex]["dimension"];
    int len = getDimension();
    // double ub = func_config["upper_bound"];
    // double lb = func_config["lower_bound"];

    double* shftx = new double[len];
    for (int j = 0; j < len; j++) {
        shftx[j] = x[j] - xopt[0][j];
    }

    double *input_x;       

    if(if_rotate == true){
        input_x = multiply(shftx, R_global, len);
    }else{
        input_x  = new double[len];
        memcpy(input_x,shftx,len*sizeof(double));
    }        

    string funcType = func_config["base_function"];
    double res = 0;
    if(funcType=="elliptic"){
        res = elliptic(input_x, len);
    }
    else if(funcType=="rastrigin"){
        res = rastrigin(input_x, len);
    }
    else if(funcType=="schwefel"){
        res = schwefel(input_x, len);
    }
    else if(funcType=="ackley"){
        res = ackley(input_x, len);
    }
    else if(funcType=="rosenbrock"){
        res = rosenbrock(input_x, len);
    }
    else if(funcType=="griewank"){
        res = griewank(input_x, len);
    }
    
	delete[] input_x;
    delete[] shftx;
    
    eva_count += group_num;
        
    return res;
}

double Benchmarks::local_eva_type1(double* x, int groupIndex) {
    if (groupIndex < 0 || groupIndex >= group_num) {
		cout << "groupIndex error\n";
		return 0;
	}
    int len = func_config["subproblems"][groupIndex]["dimension"];
    double ub = func_config["upper_bound"];
    double lb = func_config["lower_bound"];

    double res = 0;
    string funcType = func_config["subproblems"][groupIndex]["base_function"];
    double *shftx;
    double *rotate_x;

    if(funcType=="quadratic"){
        shftx = new double[len];
        for (int j = 0; j < len; j++) {
            int index = group[groupIndex][j];
            // if (x[index] > ub || x[index] < lb) {
            //     cout << "solution out of range;" << endl;
            //     return 0;
            // }
            shftx[j] = x[index];
        }
        if(if_rotate==true)
            rotate_x = multiply(shftx, R[groupIndex], len);
        else
            rotate_x = shftx;
        
        for(int j=0;j<len;j++){
            res += rotate_x[j]*shftx[j] + xopt[groupIndex][j]*shftx[j];
        }
    }
    else{
        shftx = new double[len];
        for (int j = 0; j < len; j++) {
            int index = group[groupIndex][j];
            if (x[index] > ub || x[index] < lb) {
                cout << "solution out of range;" << endl;
                return 0;
            }
            shftx[j] = x[index] - xopt[groupIndex][j];
        }

        if(if_rotate==true)
            rotate_x = multiply(shftx, R[groupIndex], len);
        else
            rotate_x = shftx;

        if(funcType=="elliptic"){
            res = elliptic(rotate_x, len);
        }
        else if(funcType=="rastrigin"){
            res = rastrigin(rotate_x, len);
        }
        else if(funcType=="schwefel"){
            res = schwefel(rotate_x, len);
        }
        else if(funcType=="ackley"){
            res = ackley(rotate_x, len);
        }
        else if(funcType=="rosenbrock"){
            res = rosenbrock(rotate_x, len);
        }
        else if(funcType=="griewank"){
            res = griewank(rotate_x, len);
        }
    }

	delete[] shftx;
    if(if_rotate)
    	delete[] rotate_x;

    eva_count += 1;

    return res;

}

double Benchmarks::global_eva(double* x) {
	double res = 0;
    if(funcClass == 1 || funcClass == 3){
        for(int i=0;i<group_num;i++){
            double r = local_eva(x,i);
            if(r==0){
                res = 0;
                // cout<<"evaluate error!\n";
                break;
            }
            res += r;
        }
    }else if(funcClass == 2){
        if(if_heterogeneous){
            for(int i=0;i<group_num;i++){
                double r = local_eva(x,i);
                res += r;
            }
            res/=group_num;
        }else{
            res += global_eva_type2(x);
        }
    }else if(funcClass == 4 || funcClass == 5){
        for(int i=0;i<group_num;i++){
            double r = local_eva(x,i);
            res += r;
        }
        res/=group_num;
    } 
    try{
        double opt = func_config["opt"];
        res -= opt;
    }catch(...){}

	return res;
}

double Benchmarks::getLocalOpt(int groupIndex){
    double fopt = 0;
    try{
        fopt = func_config["subproblems"][groupIndex]["fopt"];
    }catch(...){
        if(funcClass == 2){
            int dim = getDimension();
            for(int i=0;i<dim;i++){
                fopt -= weight*fabs(A[groupIndex][i])*getMaxX();
            }
        }
    }
    
    return fopt;

}

double Benchmarks::getMinX() {
	return func_config["lower_bound"];
}
double Benchmarks::getMaxX() {
	return func_config["upper_bound"];
}
int Benchmarks::getGroupNum() {
	return func_config["group_num"];
}
int Benchmarks::getDimension(){
    return func_config["dimension"];
}

// template<typename T>
// T** Benchmarks::read_data(string fileName) {
//     int row,col;
// 	T **data = nullptr;
//     ifstream file(fileName);
// 	cout << fileName;
// 	if (file.is_open()) {
// 		cout << " is opened;\n";
// 		string a;
// 		data = new T*[row];
// 		for (int i = 0; i < row; i++) {
// 			data[i] = new T[col];
// 			for (int j = 0; j < col; j++) {
// 				file >> a;
// 				data[i][j] = stod(a);
// 			}
// 		}
// 		file.close();
// 	}
// 	else {
// 		cout << " can not be opened;\n";
// 	}
// 	return data;
// }

template<typename T>
T** Benchmarks::read_data(string fileName) {
	// cout << fileName<<endl;
    T** res = nullptr;
    ifstream file(fileName);
	if (file.is_open()) {
		// cout << " is opened;\n";
        vector<vector<T>> data;
		string line;
        while(getline(file,line)){
    // cout<<"1"<<endl;
            stringstream ss(line);
            vector<T> rowData;
            T x;
            while(ss>>x){
                rowData.push_back(x);
            }
            data.push_back(rowData);
        }
		file.close();
    // cout<<"1"<<endl;

        res = new T*[data.size()];
        int count = 0;
        for(auto rowData:data){
    // cout<<"1"<<endl;
            // cout<<rowData.size()<<endl;
            T* row = new T[rowData.size()];
            memcpy(row,&rowData[0],rowData.size()*sizeof(T));
            res[count]=row;
            count++;
        }
	}
	else {
		// cout << " can not be opened;\n";
	}
	return res;
}

vector<int> Benchmarks::getOverlapDim(int g1,int g2){
    vector<int> dim1 = getGroupDim(g1);
    vector<int> dim2 = getGroupDim(g2);
    vector<int> overlap;
    sort(dim1.begin(),dim1.end());
    sort(dim2.begin(),dim2.end());
    set_intersection(dim1.begin(),dim1.end(),dim2.begin(),dim2.end(),back_inserter(overlap));
    return overlap;
}


vector<int> Benchmarks::getOverlapDimIndex(int g1,int g2){
    vector<int> overlap = getOverlapDim(g1,g2);
    vector<int> groupDim = getGroupDim(g1);
    vector<int> dimIndex;
    dimIndex.resize(overlap.size());
    for(unsigned int i=0;i<overlap.size();i++){
        for(unsigned int j=0;j<groupDim.size();j++){
            if(overlap[i] == groupDim[j]){
                dimIndex[i] = j;
                break;
            }
        }
    }
    return dimIndex;
}

vector<int> Benchmarks::getGroupDim(int groupIndex){
    if(this->overlap_grouping){
        int size = func_config["subproblems"][groupIndex]["dimension"];
        vector<int> v(group[groupIndex],group[groupIndex]+size);
        return v;
    }else{
        int size = func_config["subproblems"][groupIndex]["dimension"];
        vector<int> v(group[groupIndex],group[groupIndex]+size);
        vector<int> overlap = getOverlapGroup(groupIndex);
        for(int g:overlap){
            if(g<groupIndex){
                int s = func_config["subproblems"][g]["dimension"];
                vector<int> neighbor(group[g],group[g]+s);
                for(vector<int>::iterator it = v.begin();it!=v.end();){
                    if(find(neighbor.begin(),neighbor.end(),*it)!=neighbor.end()){
                        it = v.erase(it);
                    }else{
                        it++;
                    }
                }
            }
        }
        return v;
    }
} 

vector<int> Benchmarks::getGroupExcluDim(int groupIndex){
    int size = func_config["subproblems"][groupIndex]["dimension"];
    vector<int> v(group[groupIndex],group[groupIndex]+size);
    vector<int> overlap = getOverlapGroup(groupIndex);
    for(int g:overlap){
        int s = func_config["subproblems"][g]["dimension"];
        vector<int> neighbor(group[g],group[g]+s);
        for(vector<int>::iterator it = v.begin();it!=v.end();){
            if(find(neighbor.begin(),neighbor.end(),*it)!=neighbor.end()){
                it = v.erase(it);
            }else{
                it++;
            }
        }
        
    }
    return v;
    
} 

vector<int> Benchmarks::getOverlapGroup(int groupIndex){
    vector<int> groups;
    if(W){
        for(int i=0;i<group_num;i++){
            if(i == groupIndex)
                continue;
            if(W[groupIndex][i]>0)
                groups.push_back(i);
        }
    }else{
        vector<int> overlap = func_config["subproblems"][groupIndex]["overlap"];
        for(int i=0;i<group_num;i++){
            if(i == groupIndex)
                continue;
            if(overlap[i]>0)
                groups.push_back(i);
        }
    }
    return groups;
}

double** Benchmarks::getNetworkGraph(){
    if(W==nullptr){
        string data_path = "../Benchmarks/data/";
        W=read_data<double>(data_path+"W_"+funcID);
    }
    return W;
}