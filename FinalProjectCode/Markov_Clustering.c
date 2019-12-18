#include <stdio.h>
#include <stdlib.h> 
#include <mpi.h>
#include <math.h>
// to compile 
// mpicc -o Markov Markov_Clustering.c -lm

int main(int argc, char *argv[]){
    MPI_Init(&argc, &argv); //initialize MPI environment
    
    // variable to hold number of processors
    int world_size;
    
    //variable to hold rank of processor
    int world_rank;
    
    // get number of processors and rank of processor
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    
    // get program start time
    double start_time = MPI_Wtime();
        
    // loop counter variables
    int i, j;
    
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
    
    //hold the value for sqrt(p) since it is used so often
    int sqrt_p = sqrt(world_size);
    
    // indices of the column group
    int column_indices[sqrt_p];
    
    // data array to hold intiial markov matrix (only processor 0 uses this), must be one dimensional for MPI_Scatter function
    float *data;
    
    //request variables to hold non-blocking send request data
    MPI_Request req1, req2;
    
    
    
    ////////////////////////////    Get Data    ////////////////////////////////////////
    // processor 0 must get data and add self loops, needs to be 1-D for easy data transfer to
    // other processors
    if(world_rank == 0){
        int node_1, node_2;
        float weight;
        
        data = (float*)malloc(sizeof(float) * world_size);
        for(i = 0; i < world_size; i++){
            data[i] = 0.0;
        }
        
        // file handler
        FILE *infile;
        
        // open file containing graph data
        infile = fopen("graph2.txt","r");
        
        while(fscanf(infile,"%d %d %8f",&node_1, &node_2, &weight) != 0){
            data[node_1 * sqrt_p + node_2] = weight;
            data[node_2 * sqrt_p + node_1] = weight;
        }
        
        // add self loops
        for(i = 0; i < world_size; i+=(sqrt_p + 1)){
            data[i] = 1.0;
        }
        
        // close file
        fclose(infile);
        
    }
    
    if(world_rank == 0){
        printf("Original Matrix\n");
            for(i = 0; i < sqrt_p; i++){
                for(j = 0; j < sqrt_p; j++){
                    printf("%8f ", data[sqrt_p * i + j]);
                }
            printf("\n");
            }
        printf("\n");
    }
    
    // Scatter data to other processors, and copy data into second value
    MPI_Scatter(data, 1, MPI_FLOAT, &matrix_val, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
    
    ////////////////////////////    Split processors   ////////////////////////////////////////
    // create subsets of the communication world for easier operations
    MPI_Comm row_comm;
    MPI_Comm_split(MPI_COMM_WORLD, world_rank / sqrt_p, world_rank, &row_comm);
    
    int size;
    MPI_Comm_size(row_comm, &size);
    
    
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
    
    // normalize matrix values before aligning
    MPI_Allreduce(&matrix_val, &column_sum, 1, MPI_FLOAT, MPI_SUM, column_comm);
    matrix_val = matrix_val / column_sum;
    
    // copy matrix val to second matrix variable
    matrix_val_2 = matrix_val;
    
    MPI_Gather(&matrix_val, 1, MPI_FLOAT, data, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
    if(world_rank == 0){
        printf("Markov Matrix\n");
        for(i = 0; i < sqrt_p; i++){
            for(j = 0; j < sqrt_p; j++){
                //printf("Final Data[%i][%i] value: %8f\n", i, j, data[i][j]);
                printf("%8f ", data[sqrt_p * i + j]);
            }
            printf("\n");
        }
        printf("\n");
    }

    ////////////////////////////    Perform Markov Clustering    ////////////////////////////////////////
    // get row and column ranks
    int row_rank, column_rank;
    MPI_Comm_rank(row_comm, &row_rank);
    MPI_Comm_rank(column_comm, &column_rank);
    
    while(difference > 0.000001){
        // reset column sum, and expanded val
        column_sum = 0;
        expanded_val = 0;
        
        // align matrices, this must be done after each value change
        MPI_Isend(&matrix_val, 1, MPI_FLOAT, ((row_rank - row_number) + sqrt_p) % sqrt_p, 0, row_comm, &req1);
        MPI_Recv(&matrix_val, 1, MPI_FLOAT, (row_rank + row_number) % sqrt_p, 0, row_comm, MPI_STATUS_IGNORE);
        MPI_Wait(&req1, MPI_STATUS_IGNORE);
    
        MPI_Isend(&matrix_val_2, 1, MPI_FLOAT, ((column_rank - column_number) + sqrt_p) % sqrt_p, 1, column_comm, &req2);
        MPI_Recv(&matrix_val_2, 1, MPI_FLOAT, (column_rank + column_number) % sqrt_p, 1, column_comm, MPI_STATUS_IGNORE);
        MPI_Wait(&req2, MPI_STATUS_IGNORE);
    
        /*
        // print out aligned matrices DEBUGGING
        MPI_Gather(&matrix_val, 1, MPI_FLOAT, data, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
        if(world_rank == 0){
            printf("Aligned Row Matrix\n");
            for(i = 0; i < sqrt_p; i++){
                for(j = 0; j < sqrt_p; j++){
                    //printf("Final Data[%i][%i] value: %8f\n", i, j, data[i][j]);
                    printf("%8f ", data[sqrt_p * i + j]);
                }
                printf("\n");
            }
        }
    
        MPI_Gather(&matrix_val_2, 1, MPI_FLOAT, data, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
        if(world_rank == 0){
            printf("Aligned Column Matrix\n");
            for(i = 0; i < sqrt_p; i++){
                for(j = 0; j < sqrt_p; j++){
                    //printf("Final Data[%i][%i] value: %8f\n", i, j, data[i][j]);
                    printf("%8f ", data[sqrt_p * i + j]);
                }
                printf("\n");
            }
        }
        */
        
        // perform expansion operations (Cannon's Algorithm)
        for(i = 0; i < sqrt_p; i++){
            // sum up product of matrix values
            expanded_val += (matrix_val * matrix_val_2);
            
            //non blocking send to next column (cicular array)
            MPI_Isend(&matrix_val, 1, MPI_FLOAT, (row_rank - 1 + sqrt_p) % sqrt_p, 0, row_comm, &req1);
            MPI_Recv(&matrix_val, 1, MPI_FLOAT, (row_rank + 1) % sqrt_p, 0, row_comm, MPI_STATUS_IGNORE);
            
            
            // non blocking send to the next row (cicular array)
            MPI_Isend(&matrix_val_2, 1, MPI_FLOAT, (column_rank - 1 + sqrt_p) % sqrt_p, 1, column_comm, &req2);
            MPI_Recv(&matrix_val_2, 1, MPI_FLOAT, (column_rank + 1) % sqrt_p, 1, column_comm, MPI_STATUS_IGNORE);
            MPI_Wait(&req1, MPI_STATUS_IGNORE);
            MPI_Wait(&req2, MPI_STATUS_IGNORE);
        }
        /*
        // output expanded matrix
        MPI_Gather(&expanded_val, 1, MPI_FLOAT, data, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
        if(world_rank == 0){
            printf("Expanded Matrix\n");
            for(i = 0; i < sqrt_p; i++){
                for(j = 0; j < sqrt_p; j++){
                    printf("%8f ", data[sqrt_p * i + j]);
                }
                printf("\n");
            }
            printf("\n");
        }*/
        
        // perform first inflation step
        inflated_val = pow(expanded_val, inflation_factor);
        
        // get sum of column, and store result in column_sum of each processor in column
        MPI_Allreduce(&inflated_val, &column_sum, 1, MPI_FLOAT, MPI_SUM, column_comm);
        
        // normalize inflated matrix values
        inflated_val /= column_sum;
    
        // calculate difference between expanded matrix, and inflated matrix
        difference = fabs(expanded_val - inflated_val);
        
        //store new values 
        matrix_val = inflated_val;
        matrix_val_2 = matrix_val;
    
        // get sum the total difference between expanded and inflated values
        MPI_Allreduce(&difference, &difference, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
        /*
        // output next iteration matrix
        MPI_Gather(&matrix_val, 1, MPI_FLOAT, data, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
        if(world_rank == 0){
            printf("Next iteration\n");
            for(i = 0; i < sqrt_p; i++){
                for(j = 0; j < sqrt_p; j++){
                    //printf("Final Data[%i][%i] value: %8f\n", i, j, data[i][j]);
                    printf("%8f ", data[sqrt_p * i + j]);
                }
            printf("\n");
            }
            printf("\n");
        }*/

    }
    
    
    ////////////////////////////    Finalize Program    ////////////////////////////////////////
    // output final matrix
    MPI_Gather(&matrix_val, 1, MPI_FLOAT, data, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
    if(world_rank == 0){
        printf("Final Matrix\n");
        for(i = 0; i < sqrt_p; i++){
            for(j = 0; j < sqrt_p; j++){
                //printf("Final Data[%i][%i] value: %8f\n", i, j, data[i][j]);
                printf("%8f ", data[sqrt_p * i + j]);
            }
        printf("\n");
        }
    }

    //free up dynamic memory
    if(world_rank == 0){
        free(data);
    }
    
    double end_time = MPI_Wtime();
    if(world_rank == 0){
        printf("Time difference: %8f\n", end_time - start_time);
    }
    MPI_Finalize(); // close MPI environment
    return 0;
}
