import sys
import spkmeansmodule as spkm
import pandas as pd
import numpy as np


def initialize_centroids(vector_dataframe, K):
    # Set random seed as 0 as specified in requirements
    np.random.seed(0)
    
    # Initialize centroid and centroid indices lists
    centroids = [0]*K
    centroids_indices = [0]*K
    vectors = vector_dataframe.to_numpy()

    # Randomly select first initial centroid
    centroids_indices[0] = np.random.choice(len(vectors))
    centroids[0] = vectors[centroids_indices[0]]

    # Initialize vector distances array
    vector_distances = np.full(len(vectors), sys.float_info.max)

    # Main Loop for selecting all other initial centroids
    i = 0
    probability_distribution = [0]*len(vectors)
    while (i<K-1):
        # For all 0<=l<=N-1  (N=number_of_vectors), set D[l] as vector with minimum distance
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
        centroids_indices[i] = np.random.choice(len(vectors), p=probability_distribution)
        centroids[i] = vectors[centroids_indices[i]]
    
    return centroids_indices, centroids
    

def run_kmeans(K, EPSILON, max_iter, vector_dataframe):
    # Sort dataframe and calculate initial centroids
    vector_dataframe.sort_index(inplace=True)
    centroids_indices, centroids = initialize_centroids(vector_dataframe, K)
    centroids = list(np.array(centroids).flatten())

    vector_dataframe = vector_dataframe.to_numpy()
    number_of_vectors =  len(vector_dataframe)
    dimension = len(vector_dataframe[0])

    vector_dataframe = list(vector_dataframe.flatten())
    # Run C library's kms.fit
    new_centroids = spkm.fit(vector_dataframe, centroids, number_of_vectors, dimension, K, max_iter, float(EPSILON))
    return centroids_indices, new_centroids


def dimension_reduction(vector_dataframe, K):
    pass


def read_data(file_name):
    # Read data from input files
    vector_dataframe = pd.read_csv(file_name, header=None).set_index(0)
    return vector_dataframe


if __name__ == '__main__':    
    try:
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

        vector_dataframe = read_data(file_name)
        if operation == "spk":
            reduced_data, K = dimension_reduction(vector_dataframe, K)
            initial_centroid_indices, centroids = run_kmeans(K, EPSILON, max_iter, reduced_data)

            if initial_centroid_indices is not None:
                print(','.join([str(i) for i in initial_centroid_indices]))
                for vec in centroids:
                    print(','.join(["%.4f" % i for i in vec]))

        elif operation == "wam":
            pass
        elif operation == "ddg":
            pass
        elif operation == "lnorm":
            pass
        elif operation == "jacobi":
            pass 

    except Exception:
        print('An Error Has Occurred')
