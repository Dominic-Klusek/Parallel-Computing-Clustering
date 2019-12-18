#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>
#include<math.h>

// global variables
int num_features;
int world_size = 5;
int num_data;
float **data;
 
// function prototypes
float calc_euclidean_distance(float *coord1, float *coord2);
int check_for_duplicates(float **found_coords, float coord[]);
void find_min_feature_values(int *min_features);
void find_max_feature_values(int *max_features);
void mean_shift_clustering(float *new_coord, float window_size);


int main(int argc, char **argv){
    // seed random number generator
    srand(time(0));
    
    // intialize MPI
    MPI_Init(&argc, &argv);
    
    //get start time of program
    double start_time = MPI_Wtime();
    
    int i = 0, j = 0, mod = 0;
    float temp;
    float **generated_coordinate;
    
    
    // file handler
    FILE *infile;

    // open file containing graph data
    infile = fopen("iris.data","r");
    
    // get number of data points
    fscanf(infile,"%d", &num_data);
    
    // get number of data features
    fscanf(infile,"%d", &num_features);

    // allocate  space for rows of data points
    data = (float **)malloc(num_data * sizeof(float*));
    
    // load data set
    while(fscanf(infile,"%f", &temp) != 0){
        // allocate space for row of data
        if(j == 0){
            data[i] = (float *)malloc(num_features * sizeof(float));
        }
        
        // store retrieved data points
        data[i][j] = temp;
        j++;
        
        // when all data features are loaded for a single data point, reset variables
        if(j == num_features){
            j = 0;
            i++;
        }
    }
    
    // close file
    fclose(infile);
    
    
    // create array to hold generated coordinates
    generated_coordinate = (float **)malloc(world_size * num_features * sizeof(float)); //make array to hold coordinates in row major fashion
        
    // get minimum of data points
    int *min_features = (int *)malloc(num_features * sizeof(int));
    for(i = 0; i < num_features; i++){
        min_features[i] = 0;
    }
    find_min_feature_values(min_features);

    // get maximum of data points
    int *max_features = (int *)malloc(num_features * sizeof(int));
    for(i = 0; i < num_features; i++){
        max_features[i] = 0;
    }
    find_max_feature_values(max_features);

        
    // j will be used to keep track of the feature that we are working on
    j = 0;
    // generate random coordinates, and then scatter coordinates to processors
    for(i = 0; i <= world_size; i++){
        generated_coordinate[i] = (float *)malloc(num_features * sizeof(float));
        
        for(j = 0; j < num_features; j++){
            generated_coordinate[i][j] = rand() % max_features[j] + min_features[j];
        }
        
        if(j == num_features){
            j = 0;
        }
    }
    
    // free memory of min and max features
    free(min_features);
    free(max_features);
    printf("Created Coordinates\n");
    
    for(i = 0; i < world_size; i++){
        //perform Mean Shift CLustering
        mean_shift_clustering(generated_coordinate[i], 10.0);
    }
    printf("Completed Clustering\n");
    
    
    // some variables to hold data before processing
    int z, valid_coords = 0;
    float **unique_centroids = (float **)malloc(world_size * sizeof(float *));
        
    // allocate memory to contain a max of world_size uniqye centroids
    for(i = 0; i < world_size; i++){
            unique_centroids[i] = (float *)malloc(num_features * sizeof(float));
    }
        
    // place the centroid found by processor 0
    for(i = 0; i < num_features; i++){
        unique_centroids[0][i] = generated_coordinate[0][i];
    }
    valid_coords++;
        
    printf("\nFound Clusters\n");
    // retrieve coordinates from other processors
    for(i = 1; i < world_size; i++){
        //printf("Processor 0 received coordinates from Processor %i\n", i);
        for(j = 0; j < num_features; j++){
                printf("%f ", generated_coordinate[i][j]);
        }
        printf("\n");
        
        // check to see it this is a unique coordinate, store if unique
        if(check_for_duplicates(unique_centroids, generated_coordinate[i]) == 0){
            for(j = 0; j < num_features; j++){
                unique_centroids[valid_coords][j] = generated_coordinate[i][j];
            }
            valid_coords++;
        }
    }
        
    // print the total number of unique centroids
    printf("\n\nNumber of unique centroids: %i\n", valid_coords);
        
    // print the unique centroids and free up dynamic memory
    for(i = 0; i < valid_coords; i++){
        printf("UNique Centroids %i: ", i);
        for(j = 0; j < num_features; j++){
            printf("%f ", unique_centroids[i][j]);
        }
        printf("\n");
        free(unique_centroids[i]);
    }
    free(unique_centroids);
    
    
    // free dynamically allocatted memory
    for(i = 0; i < num_data; i++){
        free(data[i]);
    }
    
    for(i = 0; i < world_size; i++){
        free(generated_coordinate[i]);
    }
    
    double end_time = MPI_Wtime();
    printf("\nTime difference: %f\n", end_time - start_time);
    MPI_Finalize(); // close MPI environment
    
    return 0;
}

// function to calculate the euclidean distance between two points
float calc_euclidean_distance(float *coord1, float *coord2){
    // function to calculate the euclidean distance
    // distance variable 
    float distance = 0.0;
    
    // for each data feature calculate the 
    for(int i = 0; i < num_features; i++){
        distance += pow(coord1[i] - coord2[i], 2);
    }
    
    // return the total distance
    return distance;
}

// function to check if coord is already in found_coords 
int check_for_duplicates(float **found_coords, float coord[]){
    int i, j, similar;
    
    // assume that the coordinates are the same, then look for differences
    for(i = 0; i < world_size; i++){
        similar = 1;
        for(j = 0; j < num_features; j++){
            similar = similar & (found_coords[i][j] == coord[j]);
            //printf("%f\n", found_coords[i][j]);
        }
        if(similar == 1){
            return 1;
        }
    }
    
    
    // if no match is found, return 0 to add the element to the list
    return 0;
}


// find minimum features for data
void find_min_feature_values(int *min_features){
    // array to store minimum values for each data feature
    int i, j;
    
    // go through each point and data feature
    for(i = 0; i < num_data; i++){
        for(j = 0; j < num_features; j++){
            // if a smaller value is found for the data feauture, store in min array
            if(data[i][j] < min_features[j]){
                min_features[j] = floor(data[i][j]);
            }
        }
    }
}

// find maximum features for data
void find_max_feature_values(int *max_features){
    int i, j;
    
    // go through each point and data feature
    for(i = 0; i < num_data; i++){
        for(j = 0; j < num_features; j++){
            // if a smaller value is found for the data feauture, store in min array
            if(data[i][j] > max_features[j]){
                max_features[j] = ceil(data[i][j]);
            }
        }
    }
}

// function to perform mean shift clustering
void mean_shift_clustering(float *new_coord, float window_size){
    int i, j, valid_points;
    float *mean_coord = (float *)malloc(num_features * sizeof(float));
    float *previous_coord = (float *)malloc(num_features * sizeof(float));
    
    do{
        valid_points = 0;
        // set value of mean coord to 0
        for(i = 0; i < num_features; i++){
            mean_coord[i] = 0;
        }
        
        // find valid points
        for(i = 0; i < num_data; i++){
            //printf("# Valid Points: %d\n", valid_points);
            if(calc_euclidean_distance(new_coord, data[i]) <= window_size){
                valid_points++;
                for(j = 0; j < num_features; j++){
                    mean_coord[j] += data[i][j];
                }
            }
        }
        
        //printf("# Valid Points: %d\n", valid_points);
        if(valid_points == 0){
            window_size *= 2;
        }
        else{
            // calculate mean coordinate, and store in new_coord
            for(i = 0; i < num_features; i++){
                previous_coord[i] = new_coord[i];
                new_coord[i] = mean_coord[i] / valid_points;
            }
        }
        
    } while(calc_euclidean_distance(previous_coord, new_coord) > 0.01);
    //printf("Finished Algorithm\n");
    
    
    free(mean_coord);
    free(previous_coord);
}
