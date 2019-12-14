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
float calc_euclidean_distance(float coord1[], float coord2[]);
int check_for_duplicates(float **found_coords, float coord[]);
void find_min_feature_values(int *min_features);
void find_max_feature_values(int *max_features);


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
        printf("World_size * NUmber of features = %i\n", world_size * num_features);
        
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
        printf("Size of Generated_Coords: %i\n", sizeof(generated_coordinate)/sizeof(float));
        for(i = 0; i < world_size * num_features; i++){
            printf("%f\n", generated_coordinate[i]);
        }
        
        free(min_features);
        free(max_features);
    }
    
    MPI_Scatter(generated_coordinate, num_features, MPI_FLOAT, processor_coordinates, num_features, MPI_FLOAT, 0, MPI_COMM_WORLD);
    
    printf("Gotten Coordinate: %f, %f, %f, %f\n", processor_coordinates[0], processor_coordinates[1], processor_coordinates[2], processor_coordinates[3]);
    
    //perform Mean Shift CLustering
    
    // send final coordinates to Processor 0
    
    // if(rank == 0), check if similar
    
    // free dynamically allocatted memory
    for(i = 0; i < num_data; i++){
        free(data[i]);
    }
    free(processor_coordinates);
    
    if(world_rank == 0){
        free(generated_coordinate);
    }
    
    // finalize mpi
    MPI_Finalize();
    printf("Got here\n");
    return 0;
}

// good 
float calc_euclidean_distance(float coord1[], float coord2[]){
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
            printf("%f\n", found_coords[i][j]);
        }
        if(similar == 1){
            return 1;
        }
        printf("Look at me go\n");
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
