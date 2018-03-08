#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>

typedef enum {Euclidean, Manhattan} DistanceFormulae;
typedef enum {QuickSort, InsertionSort, SelectionSort } SortFormulae;

int dimension = 0;
int numColumns = 0;
int numberOfElementsInQuerySet = 0;
int numberOfElementsInReferenceSet = 0;
int k;
DistanceFormulae distanceFormulae;
SortFormulae sortFormulae;
FILE *fileToWrite;

typedef struct Element {
    double output;
    double inputs[];
} record_t;

void KNN(record_t** qiArray, record_t ** pjArray);
double*** getDistances(record_t** qiArray, record_t ** pjArray);
double** getDistancePerQueryPoint(record_t * qi, record_t ** pjArray, int pjArrayLength);
double euclideanDistance(double * qi, double * pj);
double manhattanDistance(double * qi, double * pj);

void sortElements(double *** array, int arraySize);

void printArray( double ** array,int dimension);
void printArray3d( double *** array, int numElements,int dimension);

int main() {

    int numMaxElementsConfig = 3;
    int numDimensionConfig = 6;
    int numMaxElements[4] = {100,250,1000,600};//,2000,4000,8000,16000};
    int dimensions[6] = {32,40,128,256,512,600};
    int numberOfTimesRepeated = 8;

    //https://en.wikipedia.org/wiki/C_date_and_time_functions
     /* Obtain current time. */
    time_t current_time;
    char* c_time_string;
    current_time = time(NULL);

    if (current_time == ((time_t)-1))
    {
        (void) fprintf(stderr, "Failure to obtain the current time.\n");
        exit(EXIT_FAILURE);
    }

    /* Convert to local time format. */
    c_time_string = ctime(&current_time);

    if (c_time_string == NULL)
    {
        (void) fprintf(stderr, "Failure to convert the current time.\n");
        exit(EXIT_FAILURE);
    }


    //strcat(c_time_string,".txt");

    char str[2000] ;
    strcat(str, "./results/ ");
    strcat(str, c_time_string);
    strcat(str, ".txt");
    fileToWrite = fopen(str, "w");
    if (fileToWrite == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }
  
    FILE* stream = fopen("mnist_train.csv", "r");
    for(int total = 0; total < numMaxElementsConfig; total++){
        for(int d =0; d< numDimensionConfig; d++){
            for(int i=0; i<numberOfTimesRepeated; i++ ){
                fprintf(fileToWrite,"\n %s \n","***********************************************************************************");
                fprintf(fileToWrite,"NumElmennts: %d NumDimensions: %d \n",numMaxElements[total], dimensions[d] );
                fprintf(fileToWrite,"Run number: %d \n", i );
                char line[2048];

                int lineNumber = 0;
                int maxNumberOfElements = numMaxElements[total];
                int indexOfOutput = 0;
                k = 5;

                //+1 to cater for first line with feature headings
                int maxNumLinesRead = maxNumberOfElements + 1;
                //Percentage of data placed in reference set
                double percentageReferenceSet = 0.8;

                //dimension excludes output column
                dimension = dimensions[d];
                numColumns = dimension +1;

                numberOfElementsInReferenceSet = (int)(percentageReferenceSet * maxNumberOfElements);
                numberOfElementsInQuerySet = maxNumberOfElements - numberOfElementsInReferenceSet;

                record_t** qiArray= (record_t**)malloc(sizeof(record_t*)*numberOfElementsInQuerySet);
                record_t** pjArray=(record_t**)malloc(sizeof(record_t*)*numberOfElementsInReferenceSet);

                size_t count = 0;

                const char delimeter[2] = ",";


                int indexReferenceArray = 0;
                int indexQueryArray = 0;

                //Reading in CSV
                //https://stackoverflow.com/questions/12911299/read-csv-file-in-c
                while (fgets(line, 2048, stream) && (lineNumber < maxNumLinesRead))
                {
                    //if(lineNumber > 0){
                        char* tmp = strdup(line);
                        int column = 0;
                        record_t* record =(record_t*)malloc(sizeof(record_t)*numColumns+1);
                        char* token;
                        token = strtok(tmp, delimeter);
                        // https://www.tutorialspoint.com/c_standard_library/c_function_strtok.htm
                        while( token != NULL && column < numColumns) {
                            if(lineNumber > 0){
                                if(column == indexOfOutput){
                                    record->output = strtod (token, NULL);    
                                }
                                else record->inputs[column] = strtod (token, NULL);
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
                        //free(record);
                    //}
                    lineNumber = lineNumber+1;

                }
                fprintf(fileToWrite,"%s %d \n","Reference Set Size:",numberOfElementsInReferenceSet );
                //printStruct(pjArray,numberOfElementsInReferenceSet,dimension);
                fprintf(fileToWrite,"%s %d \n \n","Query Set Size:",numberOfElementsInQuerySet);
                //printStruct(qiArray,numberOfElementsInQuerySet,dimension);


                KNN(qiArray,pjArray);
                free(qiArray);
                free(pjArray);
                }
        }
    }

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

void printArray3d( double *** array, int numElements,int dimension){
    printf( "[");
    for(int x = 0; x < numElements; x++){
        printf( "[" );
        for(int y = 0; y < dimension; y++){
            printf( "%lf", array[x][y][0] );
            if(y != dimension -1){
                printf(",");
            }
        }
        printf( "] \n" );
        if(x != numElements -1){
            printf(",");
        }
    }
    printf("] \n");
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


void KNN(record_t** qiArray, record_t ** pjArray){
    double startTimeEuclideanDistance, runTimeEuclideanDistance=0;
    double startTimeManhattanDistance, runTimeManhattanDistance=0;
    double startTimeQuickSort, runTimeQuickSort=0;
    double startTimeInsertionSort, runTimeInsertionSort=0;
    double startTimeSelectionSort, runTimeSelectionSort=0;
    
    
    //Step 1 compute all the distances between qi and pj , 0 ≤ j ≤ m,0 ≤ i < n − 1,
    //Euclidean Distance
    distanceFormulae = Euclidean;
    startTimeEuclideanDistance = omp_get_wtime();
	double*** distanceArrayEuclidean = getDistances(qiArray,pjArray);
	runTimeEuclideanDistance += omp_get_wtime() - startTimeEuclideanDistance;
    fprintf(fileToWrite,"Run time for distance step (Euclidean Distance): %f seconds\n\n",runTimeEuclideanDistance);

    //Manhattan Distance
    distanceFormulae = Manhattan;
    startTimeManhattanDistance = omp_get_wtime();
	double*** distanceArrayManhattan = getDistances(qiArray,pjArray);
	runTimeManhattanDistance += omp_get_wtime() - startTimeManhattanDistance;
    fprintf(fileToWrite,"Run time for distance step (Manhattan Distance): %f seconds\n\n",runTimeManhattanDistance);

    //Step 2 
    //Quick Sort
    if(sortFormulae == QuickSort){
        startTimeQuickSort = omp_get_wtime();
	    sortElements(distanceArrayEuclidean,numberOfElementsInQuerySet);
	    runTimeQuickSort += omp_get_wtime() - startTimeQuickSort;
        fprintf(fileToWrite,"Run time for sort step (Quick Sort): %f seconds\n\n",runTimeQuickSort);
    }
    //Insertion Sort
    if(sortFormulae == InsertionSort){
        startTimeInsertionSort = omp_get_wtime();
	    sortElements(distanceArrayManhattan,numberOfElementsInQuerySet);
	    runTimeInsertionSort += omp_get_wtime() - startTimeInsertionSort;
        fprintf(fileToWrite,"Run time for sort step (Insertion Sort): %f seconds\n\n",runTimeInsertionSort);
    }
    //Selection Sort
    if(sortFormulae == SelectionSort){
        startTimeSelectionSort = omp_get_wtime();
	    sortElements(distanceArrayManhattan,numberOfElementsInQuerySet);
	    runTimeSelectionSort += omp_get_wtime() - startTimeSelectionSort;
        fprintf(fileToWrite,"Run time for sort step (Selection Sort): %f seconds\n\n",runTimeSelectionSort);
    }
    //Step 3 select the k reference points corresponding to the k smallest distances;
    findNeighboursAllData(distanceArrayEuclidean,qiArray,k,numberOfElementsInQuerySet);

    // free(distanceArrayEuclidean);
    // free(distanceArrayManhattan);
    // free(qiArray);
    // free(pjArray);
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

    if(distanceFormulae == Euclidean){
        for(int index = 0 ; index < pjArrayLength; index++){
            distanceArray[index][0] = euclideanDistance(qi->inputs,pjArray[index]->inputs);
            distanceArray[index][1] = pjArray[index]->output; 
        }
    }
    else if(distanceFormulae == Manhattan){
        for(int index = 0 ; index < pjArrayLength; index++){
            distanceArray[index][0] = manhattanDistance(qi->inputs,pjArray[index]->inputs);
            distanceArray[index][1] = pjArray[index]->output; 
        }
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

double manhattanDistance(double * qi, double * pj){
    double distance = 0;
    for(int index = 0; index < dimension; index++){
            distance += abs(qi[index] - pj[index]);
    }
    return distance;
}

//https://www.tutorialspoint.com/c_standard_library/c_function_qsort.htm
//https://stackoverflow.com/questions/6103636/c-qsort-not-working-correctly
int cmpfunc (const void ** a, const void ** b) {
    double i = *(double*)a[0];
    double r = *(double*)b[0];
    return (i > r) - (i < r);
}

void sortElements(double *** array, int arraySize){
    if(sortFormulae == QuickSort){
        for(int index = 0 ; index < arraySize; index++){
            qsort(array[index], numberOfElementsInReferenceSet, sizeof(double**), cmpfunc);
        }
    }
    else if( sortFormulae == InsertionSort){
        for(int index = 0 ; index < arraySize; index++){
            insertionSort(array[index], numberOfElementsInReferenceSet);
        }
    }
    else if(sortFormulae == SelectionSort){
        for(int index = 0 ; index < arraySize; index++){
            selectionSort(array[index], numberOfElementsInReferenceSet);
        }
    }
}

//https://www.geeksforgeeks.org/insertion-sort/
void insertionSort(double** arr, int n)
{
   double  key;
   int i, j;
   for (i = 1; i < n; i++)
   {
       key = arr[i][0];
       j = i-1;
 
       /* Move elements of arr[0..i-1], that are
          greater than key, to one position ahead
          of their current position */
       while (j >= 0 && arr[j][0] > key)
       {
           arr[j+1][0] = arr[j][0];
           j = j-1;
       }
       arr[j+1][0] = key;
   }
}

 
//https://www.geeksforgeeks.org/selection-sort/
void swap(double *xp, double *yp)
{
    double temp = *xp;
    *xp = *yp;
    *yp = temp;
}
 
void selectionSort(double** arr, int n)
{
    int i, j, min_idx;
 
    // One by one move boundary of unsorted subarray
    for (i = 0; i < n-1; i++)
    {
        // Find the minimum element in unsorted array
        min_idx = i;
        for (j = i+1; j < n; j++)
          if (arr[j][0] < arr[min_idx][0])
            min_idx = j;
 
        // Swap the found minimum element with the first element
        swap(&arr[min_idx][0], &arr[i][0]);
    }
}