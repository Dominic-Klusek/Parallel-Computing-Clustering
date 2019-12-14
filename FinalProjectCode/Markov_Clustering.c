#include<stdio.h>
#include<mpi.h>
#include<math.h>

int main(int argc, char ** argv){
    // initialize mpi environment
    MPI_Init(&argc, &argv);
    
    // variable to hold number of processors
    int world_size;
    
    //variable to hold rank of processor
    int world_rank;
    
    
    // data array to hold intiial markov matrix (only processor 0 uses this)
    float *data;
    
    // variable to hold difference between expanded and inflated matrix values
    float difference = 1;
    
    // matrix values that are held by processor
    float matrix_val;
    float matrix_val_2;
    
    // variables to hold the calculated expanded and inflated values
    float expanded_val;
    float inflated_val;
    
    // variable to hold very temporary values
    float temp;
    
    // variable to hold the sum of the column
    float column_sum;
 
    // inflation factor, which can be changed to any float value > 1
    float inflation_factor = 2.0;
    
    // indices of the column group
    int column_indices[sqrt_p];

    
    MPI_Comm_size(MPI_COMM_WORLD, &world_size); //get total number of processors
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank); //get rank of processor
    
 
    // processor 0 must get data and add self loops, needs to be 1-D for easy data transfer to
    // other processors
    if(rank == 0){
        // processor 0 gets data from source
        data = new float[num_data];
    
        // add self loops
        for(i = 0; i < sqrt_p; i++){
            data[sqrt_p*i] = 1;
        }
    }   
 
 
 
    // create subsets of the communication world for easier operations
    MPI_Comm row_comm;
    MPI_Comm_Split(MPI_COMM_WORLD, rank / sqrt_p, rank, &row_comm);
 
 
 
    // construct column comminication
    MPI_Comm column_comm;
 
    // get first processor in column
    int start = rank - rank%sqrt_p;
    int end = start + sqrt_p;
 
    // create column indices for group
    for(i = start; i < end; i++){
        column_indices[i % sqrt_p] = i;
    }
 
    // Get the group of processes in MPI_COMM_WORLD
    MPI_Group world_group;
    MPI_Comm_group(MPI_COMM_WORLD, &world_group);
 
    // Construct a group containing all of the prime ranks in world_group
    MPI_Group prime_group;
    MPI_Group_incl(world_group, sqrt_p, column_indices, &prime_group);
 
    // Create a new communicator based on the group
    MPI_Comm_create_group(MPI_COMM_WORLD, prime_group, 0, &column_comm);
 
 
 
    // Scatter data to each processor
    MPI_Scatter(data, 1, MPI_FLOAT, matrix_val, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
 
 
 
    // create Markov Matrix
    // get sum of column, and store result in column_sum of each processor in column
    MPI_Allreduce(matrix_val, column_sum, 1, MPI_FLOAT, MPI_SUM, column_comm);
 
 
 
    // get final markov associative value
    matrix_val /= column_sum;
 
    // matrix_val_2 is equal to the same value as matrix_val
    matrix_val_2 = matrix_val;
 
 
 
    // every processor performs these actions
    while(difference != 0){
        // reset column sum
        column_sum = 0;
    
        // expansion
        for(i = 0; i < sqrt_p; i++){
            // non blocking send to the next row (cicular array)
            iSend(matrix_val, rank+1 % p, row_comm);
            recv(matrix_val);
 
            //non blocking send to next column (cicular array)
            iSend(matrix_val, rank+1 % p, column_comm_comm);
            recv(matrix_val_2);
        
            expanded_val += (matrix_val * matrix_val_2);
        }
    
        // inflation
        inflated_val = matrix_val ** inflation_factor;
 
        // get sum of column, and store result in column_sum of each processor in column
        MPI_Allreduce(inflated_val, column_sum, 1, MPI_FLOAT, MPI_SUM, column_comm);
    
        // normalize inflated matrix values
        inflated_val / column_sum;
    
        // calculate difference between expanded matrix, and inflated matrix
        difference = expanded_val - inflated_val;
    
        // get sum the total difference between expanded and inflated values
        MPI_Allreduce(difference, difference, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    }
 
 
 
    // finalize MPI session
    MPI_Finalize();
}
