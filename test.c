#include<stdio.h>
#include <stdlib.h> 
#include<mpi.h>
#include<math.h>
//mpirun -np 8 --use-hwthread-cpus ./parrallelSum

int main(int argc, char *argv[]){
    MPI_Init(&argc, &argv); //initialize MPI environment
        
    // loop counter variables
    int i, j;
    
    // variable to hold number of processors
    int world_size;
    
    //variable to hold rank of processor
    int world_rank;
    
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
    
    
    MPI_Comm_size(MPI_COMM_WORLD, &world_size); //get total number of processors
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank); //get rank of processor
    
    //hold the value for sqrt(p) since it is used so often
    int sqrt_p = sqrt(world_size);
    
    // indices of the column group
    int column_indices[sqrt_p];
    
    // data array to hold intiial markov matrix (only processor 0 uses this)
    // float data[sqrt_p]; final array code, but not being used for debugging
    float *data;

    
    // processor 0 must get data and add self loops, needs to be 1-D for easy data transfer to
    // other processors
    if(world_rank == 0){
        float temp;
        
        data = (float*)malloc(sizeof(float) * world_size);
        
        // file handler
        FILE *infile;
        
        // open file containing graph data
        infile = fopen("graph.txt","r");
        
        i = 0;
        while(fscanf(infile,"%f",&temp) != 0){
            //printf("Temp val: %f\n", temp);
            data[i] = temp;
            //printf("Data[%i] value: %f\n", i, &data);
            i++;
        }
        
        // close file
        fclose(infile);
        
    }
    
    // Scatter data to other processors, and copy data into second value
    MPI_Scatter(data, 1, MPI_FLOAT, &matrix_val, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
    matrix_val_2 = matrix_val;
    
    //printf("Processor %i: %f\n", world_rank, matrix_val);
    
    //free up dynamic memory
    if(world_rank == 0){
        free(data);
    }
    
    
    // create subsets of the communication world for easier operations
    MPI_Comm row_comm;
    MPI_Comm_split(MPI_COMM_WORLD, world_rank / sqrt_p, world_rank, &row_comm);
    
    int size;
    MPI_Comm_size(row_comm, &size);
    //printf("Size of row: %i\n", size);
    
    
    
    
    // construct column comminication
    MPI_Comm column_comm;
    
    
    //align rows and columns
    int row_number = floor(world_rank / sqrt_p);
    int column_number = (world_rank % sqrt_p);
    
    
    // get first processor in column
    int start = world_rank - row_number*sqrt_p;
    //printf("Start: %i\n", start);
    int end = start + sqrt_p;
    //printf("End: %i\n", end);
 
    // create column indices for group
    for(i = 0; i < sqrt_p; i++){
        //printf("%i: %i\n", i, start + i * sqrt_p);
        column_indices[i] = start + i * sqrt_p;
    }
 
    // Get the group of processes in MPI_COMM_WORLD
    MPI_Group world_group;
    MPI_Comm_group(MPI_COMM_WORLD, &world_group);
 
    // Construct a group containing all of the prime ranks in world_group
    MPI_Group prime_group;
    MPI_Group_incl(world_group, sqrt_p, column_indices, &prime_group);
 
    // Create a new communicator based on the group
    MPI_Comm_create_group(MPI_COMM_WORLD, prime_group, 0, &column_comm);
   
    
    MPI_Comm_size(column_comm, &size);
    //printf("Size of column: %i\n", size);
    
    MPI_Comm_rank(row_comm, &size);
    //printf("World rank: %i\tRow Rank: %i\n", world_rank, size);
    
    MPI_Comm_rank(column_comm, &size);
    //printf("World rank: %i\tColumn Rank: %i\n", world_rank, size);
    
    //printf("World rank: %i\tRow Number: %i\n", world_rank, row_number);
    //printf("World rank: %i\tColumn Number: %i\n", world_rank, column_number);
    
    // normalize matrix values before aligning
    MPI_Allreduce(&matrix_val, &matrix_val, 1, MPI_FLOAT, MPI_SUM, column_comm);
    MPI_Allreduce(&matrix_val_2, &matrix_val_2, 1, MPI_FLOAT, MPI_SUM, column_comm);
    
    // align rows and columns
    //MPI_ANY_TAG
    // get row and column ranks
    int row_rank, column_rank;
    MPI_Comm_rank(row_comm, &row_rank);
    MPI_Comm_rank(column_comm, &column_rank);
    
    MPI_Request req;
    printf("Processor %i sends data to left processor %i\n", row_rank, ((row_rank - row_number) + sqrt_p) % sqrt_p);
    MPI_Isend(&matrix_val, 1, MPI_FLOAT, ((row_rank - row_number) + sqrt_p) % sqrt_p, 0, row_comm, &req);
    MPI_Recv(&matrix_val, 1, MPI_FLOAT, MPI_ANY_SOURCE, 0, row_comm, MPI_STATUS_IGNORE);
    printf("Processor %i receives data from right processor\n", world_rank);
    MPI_Wait(&req, MPI_STATUS_IGNORE);
    
    
    printf("Processor %i sends data to up processor %i\n", row_rank, ((column_rank - column_number) + sqrt_p) % sqrt_p);
    MPI_Isend(&matrix_val, 1, MPI_FLOAT, ((column_rank - column_number) + sqrt_p) % sqrt_p, 1, column_comm, &req);
    MPI_Recv(&matrix_val_2, 1, MPI_FLOAT, MPI_ANY_SOURCE, 1, column_comm, MPI_STATUS_IGNORE);
    printf("Processor %i receives data from lower processor\n", world_rank);
    MPI_Wait(&req, MPI_STATUS_IGNORE);
    
    
    MPI_Finalize(); // close MPI environment
}
