#include <stdio.h>
#include <stdlib.h> 
#include <mpi.h>
#include <math.h>
#include <time.h>
//mpirun -np 8 --use-hwthread-cpus ./parrallelSum

int main(int argc, char *argv[]){
    MPI_Init(&argc, &argv); //initialize MPI environment
    
    double start_time = MPI_Wtime();
    // loop counter variables
    int i, j, z, number_of_nodes;
    
    // variable to hold difference between expanded and inflated matrix values
    float difference = 1;
    
    // matrix values that are held by processor
    float matrix_val;
    float matrix_val_2;
    
    // variables to hold the calculated expanded and inflated values
    float **expanded_val;
    float **inflated_val;
    
    // variable to hold very temporary values
    float temp;
    
    // variable to hold the sum of the column
    float column_sum;
 
    // inflation factor, which can be changed to any float value > 1
    float inflation_factor = 2.0;
    
    // data array to hold intiial markov matrix (only processor 0 uses this)
    float **data;

    
    // processor 0 must get data and add self loops, needs to be 1-D for easy data transfer to
    // other processors
    int node_1, node_2;
    float weight;
        
    // file handler
    FILE *infile;
        
    // open file containing graph data
    infile = fopen("graph_sequentail.txt","r");
    
    fscanf(infile, "%d", &number_of_nodes);
    
    // allocate data array to hold node connection data
    data = (float**)malloc(number_of_nodes * sizeof(float *));
    
    for(i = 0; i < number_of_nodes; i++){
        data[i] = (float *)malloc(number_of_nodes * sizeof(float));
    }
    
    // get node weights
    while(fscanf(infile,"%d %d %f",&node_1, &node_2, &weight) != 0){
        data[node_1][node_2] = weight;
        data[node_2][node_1] = weight;
    }
        
    // add self loops
    for(i = 0; i < number_of_nodes; i++){
        data[i][i] = 1.0;
    }
        
    // close file
    fclose(infile);
    printf("Closed the file\n");
        
    for(i = 0; i < number_of_nodes; i++){
        for(j = 0; j < number_of_nodes; j++){
            //printf("Final Data[%i][%i] value: %f\n", i, j, data[i][j]);
            printf("%f ", data[i][j]);
        }
        printf("\n");
    }

    // create Markov Matrix
    for(i = 0; i < number_of_nodes; i++){
        column_sum = 0;
        for(j = 0; j < number_of_nodes; j++){
            column_sum += data[j][i];
        }
        for(j = 0; j < number_of_nodes; j++){
            data[j][i] /= column_sum;
        }
    }
    
    printf("Original Markov Matrix\n");
    for(i = 0; i < number_of_nodes; i++){
        for(j = 0; j < number_of_nodes; j++){
            //printf("Final Data[%i][%i] value: %f\n", i, j, data[i][j]);
            printf("%f ", data[i][j]);
        }
        printf("\n");
    }
    
    // allocate data array to hold node connection data
    expanded_val = (float**)malloc(number_of_nodes * sizeof(float *));
    
    for(i = 0; i < number_of_nodes; i++){
        expanded_val[i] = (float *)malloc(number_of_nodes * sizeof(float));
    }
    
    // allocate data array to hold node connection data
    inflated_val = (float**)malloc(number_of_nodes * sizeof(float *));
    
    for(i = 0; i < number_of_nodes; i++){
        inflated_val[i] = (float *)malloc(number_of_nodes * sizeof(float));
    }
    printf("Dynamically alocatted all memory\n");
    
    // every processor performs these actions
    while(difference >= 0.0000001){
        // perform expansion and first step of inflation operations in same loop
        for(i = 0; i < number_of_nodes; i++){
            for(j = 0; j < number_of_nodes; j++){
                expanded_val[i][j] = 0.0;
                for(z = 0; z < number_of_nodes; z++){
                    expanded_val[i][j] += (data[i][z] * data[z][j]);
                }
                inflated_val[i][j] = pow(expanded_val[i][j], inflation_factor);
            }
        }
        
        printf("Expanded iteration\n");
        for(i = 0; i < number_of_nodes; i++){
            for(j = 0; j < number_of_nodes; j++){
                //printf("Final Data[%i][%i] value: %f\n", i, j, data[i][j]);
                printf("%f ", expanded_val[i][j]);
            }
        printf("\n");
        }
        
        
        //perform second step of inflation operation
        for(i = 0; i < number_of_nodes; i++){
            column_sum = 0;
            for(j = 0; j < number_of_nodes; j++){
                column_sum += inflated_val[j][i];
            }
            for(j = 0; j < number_of_nodes; j++){
                inflated_val[j][i] /= column_sum;
            }
        }
        
        // reset difference variable
        difference = 0;
        
        // find difference and replace data values with inflated values
        for(i = 0; i < number_of_nodes; i++){
            for(j = 0; j < number_of_nodes; j++){
                difference += fabs((expanded_val[i][j] - inflated_val[i][j]));
                data[i][j] = inflated_val[i][j];
            }
        }
        
        printf("Next iteration\n");
        for(i = 0; i < number_of_nodes; i++){
            for(j = 0; j < number_of_nodes; j++){
                //printf("Final Data[%i][%i] value: %f\n", i, j, data[i][j]);
                printf("%f ", data[i][j]);
            }
        printf("\n");
        }
    }
    
    printf("Final Array\n");
    for(i = 0; i < number_of_nodes; i++){
        for(j = 0; j < number_of_nodes; j++){
            //printf("Final Data[%i][%i] value: %f\n", i, j, data[i][j]);
            printf("%f ", data[i][j]);
        }
        printf("\n");
    }
    

    for(i = 0; i < number_of_nodes; i++){
        free(inflated_val[i]);
        free(expanded_val[i]);
        free(data[i]);
    }
    
    free(inflated_val);
    free(expanded_val);
    free(data);
    
    double end_time = MPI_Wtime();
    printf("Time difference: %f\n", end_time - start_time);
    MPI_Finalize(); // close MPI environment
}
