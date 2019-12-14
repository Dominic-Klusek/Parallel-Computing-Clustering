#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>
#include<math.h>

// global variables
int num_features;
int world_size;
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
    
    int world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size); //get total number of processors
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank); //get rank of processor
    
    int i = 0, j = 0, mod = 0;
    float temp;
    float *generated_coordinate;
    
    
    // file handler
    FILE *infile;

    // open file containing graph data
    infile = fopen("iris.data","r");
    
    // get number of data points
    fscanf(infile,"%d", &num_data);
    
    // get number of data features
    fscanf(infile,"%d", &num_features);
    
    float *processor_coordinates = (float *)malloc(num_features * sizeof(float)); // coordinates that the processor will operate on
    
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
    
    if(world_rank == 0){
        // create array to hold generated coordinates
        generated_coordinate = (float *)malloc(world_size * num_features * sizeof(float)); //make array to hold coordinates in row major fashion
        
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
        for(i = 0; i <= world_size * num_features; i++){
            generated_coordinate[i] = rand() % max_features[j] + min_features[j];
            if(j == num_features){
                j = 0;
            }
        }
        
        free(min_features);
        free(max_features);
    }
    
    MPI_Scatter(generated_coordinate, num_features, MPI_FLOAT, processor_coordinates, num_features, MPI_FLOAT, 0, MPI_COMM_WORLD);
    
    printf("Gotten Coordinate: %f, %f, %f, %f\n", processor_coordinates[0], processor_coordinates[1], processor_coordinates[2], processor_coordinates[3]);
    
    //perform Mean Shift CLustering
    mean_shift_clustering(processor_coordinates, 100.0);
    
    printf("New Coordinate: %f, %f, %f, %f\n", processor_coordinates[0], processor_coordinates[1], processor_coordinates[2], processor_coordinates[3]);
    
    // send final coordinates to Processor 0
    if(world_rank == 0){
        float *temp_coord = (float *)malloc(num_features * sizeof(float));
        int z, valid_coords = 0;
        float **unique_centroids = (float **)malloc(world_size * sizeof(float *));
        for(i = 0; i < world_size; i++){
                unique_centroids[i] = (float *)malloc(num_features * sizeof(float));
        }
        // retrieve coordinates from other processors
        for(i = 1; i < world_size; i++){
            // recieve coordinate
            MPI_Recv(temp_coord, num_features, MPI_FLOAT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            printf("Processor 0 received coordinates from Processor %i\n", i);
            for(j = 0; j < num_features; j++){
                    printf("%d ", temp_coord[j]);
            }
            printf("\n");
            // check to see it this is a unique coordinate, store if unique
            if(check_for_duplicates(unique_centroids, temp_coord) == 1){
                for(j = 0; j < num_features; j++){
                    unique_centroids[valid_coords][j] = temp_coord[j];
                }
                valid_coords++;
            }
        }
        
        printf("Number of uniques: %i\n", valid_coords);
        
        for(i = 0; i < valid_coords; i++){
            printf("UNique Centroids %i: ", i);
            for(j = 0; j < num_features; j++){
                printf("%d ", unique_centroids[i][j]);
            }
            printf("\n");
            free(unique_centroids[i]);
        }
        free(temp_coord);
        //free(unique_centroids);
    } else {
        MPI_Send(processor_coordinates, num_features, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
        printf("Processor %i sent coordinates to Processor 0\n", world_rank);
    }
    
    // free dynamically allocatted memory
    for(i = 0; i < num_data; i++){
        free(data[i]);
    }
    
    printf("Got here with error\n");
    
    free(processor_coordinates);
    printf("Processor is fucked\n");
    
    if(world_rank == 0){
        free(generated_coordinate);
    }
    // finalize mpi
    MPI_Finalize();
    printf("Got here\n");
    return 0;
}

// good 
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

// good 
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


// good
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

// good
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

// good
void mean_shift_clustering(float *new_coord, float window_size){
    int i, j, valid_points= 100;
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
        
        printf("# Valid Points: %d\n", valid_points);
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
    printf("Finished Algorithm\n");
    
    
    free(mean_coord);
    free(previous_coord);
}
