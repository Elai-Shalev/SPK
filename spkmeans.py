import sys
import spkmeansmodule as spkm
import numpy as np


# Initializes K centroids of input matrix. 
# Returns initial centroids and their indicies 
def initialize_centroids(vectors, K):
    # Initialize centroid and centroid indices lists
    centroids = [0]*K
    centroids_indices = [0]*K

    # Randomly select first initial centroid
    centroids_indices[0] = np.random.choice(len(vectors))
    centroids[0] = vectors[centroids_indices[0]]

    # Initialize vector distances array
    vector_distances = np.full(len(vectors), sys.float_info.max)

    # Main Loop for selecting all other initial centroids
    i = 0
    probability_distribution = [0]*len(vectors)
    while (i<K-1):
        # For all 0<=l<=N-1  (N=number_of_vectors): 
        # set D[l] as vector with minimum distance
        for l in range(len(vectors)):
            for j in range(i+1):
                vec = (vectors[l] - centroids[j]) ** 2
                vector_distances[l] = min(vector_distances[l], np.sum(vec))

        # Set probability distribution as speficied in requirements
        sum_distances = float(np.sum(vector_distances))
        
        for l in range(len(vectors)):
            probability_distribution[l] = vector_distances[l] / sum_distances

        i = i+1
        
        # Set current initial centroid
        centroids_indices[i] = np.random.choice(len(vectors), 
                                                p=probability_distribution)
        centroids[i] = vectors[centroids_indices[i]]
    
    return centroids_indices, centroids
    

# Implementation of K-Means++ Algorithm
def run_kmeans(K, EPSILON, max_iter, vector_data):
    # Calculate initial centroids
    centroids_indices, centroids = initialize_centroids(vector_data, K)
    centroids = list(np.array(centroids).flatten())

    number_of_vectors =  len(vector_data)
    dimension = len(vector_data[0])

    vector_data = list(vector_data.flatten())
    # Run C library's kms.fit
    new_centroids = spkm.fit(vector_data, centroids, number_of_vectors, 
                             dimension, K, max_iter, float(EPSILON))
    return centroids_indices, new_centroids


def dimension_reduction(file_name, K):
    vectors_list = spkm.reduce(file_name, K)
    return vectors_list


if __name__ == '__main__':
    try:
        # Set random seed as 0 as specified in requirements
        np.random.seed(0)

        max_iter = 300
        EPSILON = 0
        if len(sys.argv) != 4:
            print('Invalid Input!')
            quit()

        try:
            K = int(sys.argv[1])
            if K < 0 or '.' in sys.argv[1]:
                print('Invalid Input!')
                quit()
        except ValueError:
            print('Invalid Input!')
            quit()

        operation = sys.argv[2]
        file_name = sys.argv[3]

        if operation == "spk":
            reduced_data = dimension_reduction(file_name, K)
            K = spkm.get_K()
            reduced_data = np.array(reduced_data)
            initial_centroid_indices, centroids = run_kmeans(K, 
                                                             EPSILON, 
                                                             max_iter, 
                                                             reduced_data)

            if initial_centroid_indices is not None:
                print(','.join([str(i) for i in initial_centroid_indices]))
                for vec in centroids:
                    print(','.join(["%.4f" % i for i in vec]))

        elif operation == "wam":
            spkm.wam(file_name, K)
        elif operation == "ddg":
            spkm.ddg(file_name, K)
        elif operation == "lnorm":
            spkm.lnorm(file_name, K)
        elif operation == "jacobi":
            spkm.jacobi(file_name, K)
        else:
            print('Invalid Input!')
            quit()

    except Exception:
        print('An Error Has Occurred')
        quit()
