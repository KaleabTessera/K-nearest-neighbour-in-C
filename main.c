#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <omp.h>

int dimension = 0;
int numColumns = 0;
int numberOfElementsInQuerySet = 0;
int numberOfElementsInReferenceSet = 0;
int k;

typedef struct Element {
    double output;
    double inputs[];
} record_t;
 

//void printArray( double ** pj, int numElements,int dimension);
void printArray( double ** array,int dimension);
double*** getDistances(record_t** qiArray, record_t ** pjArray);
double** getDistancePerQueryPoint(record_t * qi, record_t ** pjArray, int pjArrayLength);
double euclideanDistance(double * qi, double * pj);
void quickSort(double *** array, int arraySize);
void KNN(record_t** qiArray, record_t ** pjArray);
//int indexOfElement(int val, int *arr, int size);


//https://www.tutorialspoint.com/c_standard_library/c_function_qsort.htm
int cmpfunc (const void ** a, const void ** b) {
   return ( *(double*)a[0] - *(double*)b[0] );
}

int main() {

    FILE* stream = fopen("mnist_train.csv", "r");
    char line[2048];
    
    int lineNumber = 0;
    int maxNumberOfElements = 5000;
    int indexOfOutput = 0;
    k = 5;

    //+1 to cater for first line with feature headings
    int maxNumLinesRead = maxNumberOfElements + 1;
    //Percentage of data placed in reference set
    double percentageReferenceSet = 0.8;

    numColumns = 10000;
    //dimension excludes output column
    dimension = numColumns-1;
    numberOfElementsInReferenceSet = (int)(percentageReferenceSet * maxNumberOfElements);
    numberOfElementsInQuerySet = maxNumberOfElements - numberOfElementsInReferenceSet;

    record_t** qiArray= (record_t**)malloc(sizeof(record_t*)*numberOfElementsInQuerySet);
    record_t** pjArray=(record_t**)malloc(sizeof(record_t*)*numberOfElementsInReferenceSet);
    
    size_t count = 0;

    const char delimeter[2] = ",";
    
   
    int indexReferenceArray = 0;
    int indexQueryArray = 0;

    //https://stackoverflow.com/questions/12911299/read-csv-file-in-c
    while (fgets(line, 2048, stream) && (lineNumber < maxNumLinesRead))
    {
        if(lineNumber > 0){
            char* tmp = strdup(line);
            int column = 0;
            record_t* record =(record_t*)malloc(sizeof(record_t)*dimension);
            double* element = (double*)malloc(sizeof(double)*dimension);
            char* token;
            token = strtok(tmp, delimeter);
            // https://www.tutorialspoint.com/c_standard_library/c_function_strtok.htm
            while( token != NULL && column < numColumns) {
                if(lineNumber > 0){
                    if(column == indexOfOutput){
                        record->output = strtod (token, NULL);    
                    }
                    else record->inputs[column] = strtod (token, NULL);
                    //element[column] = strtod (token, NULL); 
                }
                column=column+1;
                token = strtok(NULL,delimeter);
            }
            
            if(lineNumber <= numberOfElementsInReferenceSet){
                pjArray[indexReferenceArray] = record;
                indexReferenceArray = indexReferenceArray+1;
            }
            else {
                qiArray[indexQueryArray] = record;
                indexQueryArray = indexQueryArray+1;
            } 
            free(tmp);
            //free(token);
           // free(record);
        }
        lineNumber = lineNumber+1;
    }
    printf("%s %d \n","Reference Set Size:",numberOfElementsInReferenceSet );
    //printStruct(pjArray,numberOfElementsInReferenceSet,dimension);
    printf("%s %d \n \n","Query Set Size:",numberOfElementsInQuerySet);
    //printStruct(qiArray,numberOfElementsInQuerySet,dimension);

    KNN(qiArray,pjArray);

    return 0;
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
}

void printArray( double ** array,int dimension){
        printf( "[" );
        for(int x = 0; x < dimension; x++){
            printf( "%lf", array[x][0] );
            if(x != dimension -1){
                printf(",");
            }
        }
        printf( "] \n" );
}

void printStruct( record_t ** array, int numElements,int dimension){
    printf( "[");
    for(int x = 0; x < numElements; x++){
        
        printf( "Inputs: [" );
        for(int y = 0; y < dimension; y++){
            printf( "%lf",array[x]->inputs[y]);
            if(y != dimension -1){
                printf(",");
            }
        }
        printf( "]," );
        printf("%s %f","Output: [",array[x]->output);
        printf("] \n");
        if(x != numElements -1){
            printf(",");
        }
    }
    printf("] \n");
}
void quickSort(double *** array, int arraySize){
    for(int index = 0 ; index < arraySize; index++){
        qsort(array[index], numberOfElementsInReferenceSet, sizeof(double*), cmpfunc);
    }
}

void KNN(record_t** qiArray, record_t ** pjArray){
    double startTimeDistance, runTimeDistance=0;
    double startTimeSort, runTimeSort=0;
    
    //Step 1 compute all the distances between qi and pj , 0 ≤ j ≤ m,0 ≤ i < n − 1,
    startTimeDistance = omp_get_wtime();
	double*** distanceArray = getDistances(qiArray,pjArray);
	runTimeDistance += omp_get_wtime() - startTimeDistance;
    printf("Run time for distance step: %f seconds\n\n",runTimeDistance);

    //Step 2 
    startTimeSort = omp_get_wtime();
	quickSort(distanceArray,numberOfElementsInQuerySet);
	runTimeSort += omp_get_wtime() - startTimeSort;
    printf("Run time for sort step: %f seconds\n\n",runTimeSort);

    //Step 3 select the k reference points corresponding to the k smallest distances;
    findNeighboursAllData(distanceArray,qiArray,k,numberOfElementsInQuerySet);


    // double*** euclideanDistanceArray= (double***)malloc(sizeof(double**)*numberOfElementsInReferenceSet);
    // //Step 0 & 4. repeat steps 1 to 3 for q i+1 .
    // for(int index = 0 ; index < numberOfElementsInQuerySet; index++){
    //     printf("%s","\n");
    //     //Step 1 compute all the distances between q i and p j , 0 ≤ j ≤ m;
    //     //printf("%s \n","Step 1 compute all the distances between q i and p j , 0 ≤ j ≤ m;");
    //     euclideanDistanceArray[index] = getDistancePerQueryPoint(qiArray[index],pjArray,numberOfElementsInReferenceSet);
    //     //printf("%s","Distance array: ");
    //     //printArray(euclideanDistanceArray[index],numberOfElementsInReferenceSet);
        
    //     //Step 2 sort the computed distances;
    //     //printf("\n%s","Step 2 sort the computed distances; \n");
    //     qsort(euclideanDistanceArray[index], numberOfElementsInReferenceSet, sizeof(double*), cmpfunc);
    //     //printf("%s","Sorted distance array: ");
    //     //printArray(euclideanDistanceArray[index],numberOfElementsInReferenceSet);
        
    //     //printf("\n%s","Step 3 select the k reference points corresponding to the k smallest distances \n");
    //     //Step 3 select the k reference points corresponding to the k smallest distances;
    //     findNeighbours(euclideanDistanceArray[index],k,qiArray[index]->output);
    // }
    // /return void;
}

double*** getDistances(record_t** qiArray, record_t ** pjArray){
    double*** euclideanDistanceArray= (double***)malloc(sizeof(double**)*numberOfElementsInReferenceSet);
    for(int index = 0 ; index < numberOfElementsInQuerySet; index++){
        euclideanDistanceArray[index] = getDistancePerQueryPoint(qiArray[index],pjArray,numberOfElementsInReferenceSet);
    }
    return euclideanDistanceArray;
}


double** getDistancePerQueryPoint(record_t * qi, record_t ** pjArray, int pjArrayLength){
    double **distanceArray = (double **)malloc(pjArrayLength * sizeof(double *));
    for (int i=0; i<pjArrayLength; i++)
         distanceArray[i] = (double *)malloc(2 * sizeof(double));

    for(int index = 0 ; index < pjArrayLength; index++){
        distanceArray[index][0] = euclideanDistance(qi->inputs,pjArray[index]->inputs);
        distanceArray[index][1] = pjArray[index]->output; 
    }
    return distanceArray;
}


double euclideanDistance(double * qi, double * pj){
    double distance = 0;
    for(int index = 0; index < dimension; index++){
            distance += pow((qi[index] - pj[index]),2);
    }
    return sqrt(distance);
}




