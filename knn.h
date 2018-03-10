#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>
#include "common.h" 
#include "distance.h"
#include "sort.h"
#include <string.h>


void printArray( double ** array,int dimension);
void printArray3d( double *** array, int numElements,int dimension);
void writeToFile();

double avgRunTimeEuclideanDistance;
double avgRunTimeManhattanDistance;
double avgRunTimeQuickSort; 
double avgRunTimeInsertionSort; 
double avgRunTimeSelectionSort; 


void KNN(record_t** qiArray, record_t ** pjArray,int numberOfElementsInQuerySet, int numberOfElementsInReferenceSet, int k,FILE *fileToWriteTo,int numMaxElements,int dimension, int numberOfTimesRepeated){
    double startTimeEuclideanDistance, runTimeEuclideanDistance=0;
    double startTimeManhattanDistance, runTimeManhattanDistance=0;
    double startTimeQuickSort, runTimeQuickSort=0;
    double startTimeInsertionSort, runTimeInsertionSort=0;
    double startTimeSelectionSort, runTimeSelectionSort=0;

    if(fileToWriteTo){
        fprintf(fileToWriteTo,"\n %s \n","***********************************************************************************");
        fprintf(fileToWriteTo,"NumElements: %d NumDimensions: %d \n",numMaxElements, dimension);
        //fprintf(fileToWriteTo,"Run number: %d \n", i );

        fprintf(fileToWriteTo,"%s %d \n","Reference Set Size:",numberOfElementsInReferenceSet );
        fprintf(fileToWriteTo,"%s %d \n \n","Query Set Size:",numberOfElementsInQuerySet);
    }

    
    //printStruct(qiArray,numberOfElementsInQuerySet,numberOfElementsInReferenceSet);
    //Step 1 compute all the distances between qi and pj , 0 ≤ j ≤ m,0 ≤ i < n − 1,
    //Euclidean Distance
    distanceFormulae = Euclidean;
    startTimeEuclideanDistance = omp_get_wtime();
	double*** distanceArrayEuclidean = getDistancesSerial(qiArray,pjArray,numberOfElementsInQuerySet,numberOfElementsInReferenceSet, dimension);
	runTimeEuclideanDistance += omp_get_wtime() - startTimeEuclideanDistance;
    avgRunTimeEuclideanDistance += runTimeEuclideanDistance;
    if(fileToWriteTo != NULL){
        fprintf(fileToWriteTo,"Average run time( %d iterations) for distance step (Euclidean Distance): %f seconds\n\n",numberOfTimesRepeated,avgRunTimeEuclideanDistance/numberOfTimesRepeated);
    }

    //Manhattan Distance
    distanceFormulae = Manhattan;
    startTimeManhattanDistance = omp_get_wtime();
	double*** distanceArrayManhattan = getDistancesSerial(qiArray,pjArray,numberOfElementsInQuerySet,numberOfElementsInReferenceSet, dimension);
	runTimeManhattanDistance += omp_get_wtime() - startTimeManhattanDistance;
    avgRunTimeManhattanDistance += runTimeEuclideanDistance;
    if(fileToWriteTo != NULL){
        fprintf(fileToWriteTo,"Average run time( %d iterations) for distance step (Manhattan Distance): %f seconds\n\n",numberOfTimesRepeated,avgRunTimeManhattanDistance/numberOfTimesRepeated);
    }

    // //Step 2 
    // //Quick Sort
    // if(sortFormulae == QuickSort){
    //     startTimeQuickSort = omp_get_wtime();
	//     sortElements(distanceArrayEuclidean,numberOfElementsInQuerySet,numberOfElementsInReferenceSet,QuickSort);
	//     runTimeQuickSort += omp_get_wtime() - startTimeQuickSort;
    //     if(fileToWriteTo != NULL){
    //         fprintf(fileToWriteTo,"Run time for sort step (Quick Sort): %f seconds\n\n",runTimeQuickSort);
    //     }
    // }
    // //Insertion Sort
    // if(sortFormulae == InsertionSort){
    //     startTimeInsertionSort = omp_get_wtime();
	//     sortElements(distanceArrayManhattan,numberOfElementsInQuerySet,numberOfElementsInReferenceSet,InsertionSort);
	//     runTimeInsertionSort += omp_get_wtime() - startTimeInsertionSort;
    //     if(fileToWriteTo != NULL){
    //         fprintf(fileToWriteTo,"Run time for sort step (Insertion Sort): %f seconds\n\n",runTimeInsertionSort);
    //     }
    // }
    // //Selection Sort
    // if(sortFormulae == SelectionSort){
    //     startTimeSelectionSort = omp_get_wtime();
	//     sortElements(distanceArrayManhattan,numberOfElementsInQuerySet,numberOfElementsInReferenceSet,SelectionSort);
	//     runTimeSelectionSort += omp_get_wtime() - startTimeSelectionSort;
    //     if(fileToWriteTo != NULL){
    //         fprintf(fileToWriteTo,"Run time for sort step (Selection Sort): %f seconds\n\n",runTimeSelectionSort);
    //     }
    // }
    // //Step 3 select the k reference points corresponding to the k smallest distances;
    // findNeighboursAllData(distanceArrayEuclidean,qiArray,k,numberOfElementsInQuerySet);

}



//https://stackoverflow.com/questions/15094834/check-if-a-value-exist-in-a-array
int indexOfElement(int val, double **arr, int size){
    int i;
    for (i=0; i < size; i++) {
        if (arr[i][0] == val)
            return i;
    }
    return -1;
}

int indexOfLargestElement( double ** arr, int size){
    double largest = arr[0][1];
    int index =0;
    for (int i = 1; i < size; i++) 
    {
	    if (largest< arr[i][1]){
	        largest = arr[i][1];
            index = i;
        }
	}
    return index;
}

void findNeighboursAllData(double ***array, record_t** qiArray,int k,int size){
    for(int i= 0; i < size;i++){
        findNeighbours(array[i],k,qiArray[i]->output);
    }
}

void findNeighbours(double **array, int k, double actualOutput){
    //printf("%s %d \n","Finding neighbours with k: ",k);
    double predictedOutputs[k];
    double **countOccurences = (double **)malloc(k * sizeof(double *));
    for (int i=0; i<k; i++)
         countOccurences[i] = (double *)malloc(2 * sizeof(int));

    int indexInOccurenceArray;
    for(int i=0; i< k;i++){
        predictedOutputs[i] = array[i][1];
        //printf("%s %d %s %lf","Predicted output, neighbour ",i,":" ,predictedOutputs[i]);
        //printf(" \n");
        
        indexInOccurenceArray =  indexOfElement(predictedOutputs[i],countOccurences,k);
        if(indexInOccurenceArray == -1){
            countOccurences[i][0] = predictedOutputs[i];
            countOccurences[i][1] = 1;
        } else {
            countOccurences[indexInOccurenceArray][0] = predictedOutputs[i];
            countOccurences[indexInOccurenceArray][1] =  countOccurences[indexInOccurenceArray][1] +1;
        }
        
    }
    int indexLargest = indexOfLargestElement(countOccurences,k);
    //printf("\nAssigned output:  %lf", countOccurences[indexLargest][0]);
    //printf("\nActual output:  %lf", actualOutput);
    //printf("%s","\n");

    //free(countOccurences);
}