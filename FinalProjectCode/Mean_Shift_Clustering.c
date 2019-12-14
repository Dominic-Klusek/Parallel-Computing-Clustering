/ global variables
int num_features;
 
// function prototypes
float calc_euclidean_distance(float[] coord1, float[] coord2);
bool check_for_duplicates(List<int[]> found_coords, float[] coord);
int[] find_min_feature_values(int[][] data);
int[] find_max_feature_values(int[][] data);
 
int main(int argc, char **argv){
    // variable to hold the data points
    float data_points[num_data][num_features];
    
    // variable to hold the number of centers to create
    int num_centers = 10;
    
    // variable to hold the radius of the window
    float radius = 10;
    
    // variable to hold the difference of the previous and current centroid coordinates
    float difference = 1;
    
    // variable to hold the number of valid points in window
    float valid_points;
    
    // float scaling factor
    float scaling_factor = 1.5;
    
    // array to hold the sum of all varible points
    float point_sums[num_features];
    
    // array to hold the centroid coordinates
    float prev_coordinates[num_features];
    float coordinates[num_features];
    
    // arrays to hold data about the data features
    float min_array[num_features];
    float max_array[num_features];
    
    // list to hold the final coordinates
    List<int[]> final_coordinationes;
    
    
    // initialize the MPI environment
    MPI_Init(&argc, &argv);
    
    
    
    // each processor reads data from a filelength
    data_points = load_data("data_file_name");
    
    
    
    // only processor 0 calculates the centers
    if(rank == 0){
        // find minimum and max values for each data feature
        min_array = find_min_feature_values(data_points);
        max_array = find_max_feature_values(data_points);
        
        // create centroid coordinates, and send them to each processor
        for(i=0; i < num_centers; i++){
            for(j=0; j < num_features; j++){
                coordinates[j] = rand() % max_array[j] + min_array[j];
            }
            send(coordinates, num_features, MPI_FLOAT, i, MPI_COMM_WORLD);
        }
    } else{
        // receive initial center coordinates from processor 0
        recv(coordinates, num_features, MPI_FLOAT, 0, MPI_COMM_WORLD);
    }
    
    
    
    // while point difference change is larger than alpha
    while(difference > alpha){
        valid_points = 0;
    
        // go through all the points and sum up points withing the radius
        for(i = 0; i < num_data; i++){
            difference = calc_euclidean_distance(data_points[i], coordinates);
        
            if(difference <= radius){
                valid_points += 1;
                point_sums += data_points[i];
            }
        }
        
        // if there were no valid points found, increase radius, and try again
        // else change store the current coordinates, recalculate the coordinates using the sum of the points, 
            // and then calculate difference to see if convergence is achieved
        if(valid_points == 0){
            radius *= scaling_factor;
            difference = 1;
        } else{
            prev_coordinates = coordinates;
            coordinates = point_sums/valid_points;
            difference = coordinates - prev_coordinates;
        }
    
    }
    
    
    
    // rank 0 collects final centers 
    if(rank == 0){
        // append the coordinates at processor 0 so that the list is not empty
        final_coordinationes.append(coordinates);
        
        // get center points, and store in a linked list
        for(i = 1; i < p; i++){
            recv(coordinates, num_features, MPI_FLOAT, i, MPI_COMM_WORLD);
        
            // check if coordinates are already in the list
            if(check_for_duplicates(final_coordinationes, coordinates)){
                final_coordinationes.append(coordinates);
            }
        }
    
        // output final cetners
    
    }
    
    
    
    //finalize the MPI session
    MPI_Finalize();
}
 
float calc_euclidean_distance(float[] coord1, float[] coord2){
    // function to calculate the euclidean distance
    // distance variable 
    float distance = 0.0;
    
    // for each data feature calculate the 
    for(int i = 0; i < (sizeof(coord1)/sizeof(*coord1); i++){
        distance += (coord1[i] - coord2[i]) ** 2;
    }
    
    // return the total distance
    return distance;
}
 
bool check_for_duplicates(List<int[]> found_coords, float[] coord){
    // function uses an iterator to go through entire list, and checks if the list contains the coordinate in itself
    for (std::list<int>::iterator it=found_coords.begin(); it != found_coords.end(); ++it){
        // if a match is found, return false
        if(*it == coord){
            return false;
        }
    }
    
    // if no match is found, return true to add the element to the list
    return true;
}
 
int[] find_min_feature_values(int[][] data){
    // array to store minimum values for each data feature
    int min[(sizeof(data[0])/sizeof(*data[0])] = 0;
    
    // go through each point and data feature
    for(i = 0; i < (sizeof(data)/sizeof(*data); i++){
        for(j = 0; j < num_features; j++){
            // if a smaller value is found for the data feauture, store in min array
            if(data[i][j] < min[j]){
                min[j] = data[i][j];
            }
        }
    }
}
 
int[] find_max_feature_values(int[][] data){
    // array to store maximum values for each data feature
    int max[(sizeof(data[0])/sizeof(*data[0])] = 0;
    
    // go through each point and data feature
    for(i = 0; i < (sizeof(data)/sizeof(*data); i++){
        for(j = 0; j < num_features; j++){
            // if a smaller value is found for the data feauture, store in min array
            if(data[i][j] > max[j]){
                max[j] = data[i][j];
            }
        }
    }
