#include "../framework/framework.h"
#include <iostream>
#include <mpi.h>
using namespace std;

void agent_function(Framework* handler){
    // Participants can only modify the code here

    cout<<"hello world "<<handler->get_self_id()<<endl;

    cout<<"weights ";
    vector<double> weight = handler->get_adjacent_weights();
    for(double w:weight){
        cout<<w<<" ";
    }
    cout<<endl;
    
    int myid = handler->get_self_id();
    
    vector<int> neighbor_list = handler->get_neighbor_id();
    MPI_Request *req = new MPI_Request[neighbor_list.size()];
    for(int j=0;j<neighbor_list.size();j++){
        int data = myid*10;
        handler->Message_Isend(&data,1,MPI_INT,neighbor_list[j],0,&req[j]);
    }

    for(int j=0;j<neighbor_list.size();j++){
        MPI_Status stat;
        int recv_data;
        handler->Message_Recv(&recv_data,1,MPI_INT,neighbor_list[j],0,&stat);
        printf("agent %d receive %d from agent %d\n",myid,recv_data,neighbor_list[j]);
    }

    int dim = handler->get_problem_dim();
    double* x = new double[dim]{0};
    x[0] = myid;
    cout<<myid<<" "<<handler->local_evaluation(x)<<endl;

    handler->submit_final_solution(x);
    return;
}