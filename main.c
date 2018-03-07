#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int dimension = 0;
int numColumns = 0;
int numberOfElementsInQuerySet = 0;
int numberOfElementsInReferenceSet = 0;

enum ElementType { et_str, et_int, et_dbl };
typedef struct Element {
    double output;
    enum ElementType type;
    double inputs[];
} record_t;
    

void printArray( double * pj[], int numElements,int dimension);
double** getDistances(record_t** qiArray, record_t ** pjArray);
double* getDistancePerQueryPoint(record_t * qi, record_t ** pjArray, int pjArrayLength);
double euclideanDistance(double * qi, double * pj);
void quickSort(double ** array, int arraySize);


const char* getfield(char* line, int num)
{
    const char* tok;
    for (tok = strtok(line, ";");
            tok && *tok;
            tok = strtok(NULL, ";\n"))
    {
        if (!--num)
            return tok;
    }
    return NULL;
}

//https://www.tutorialspoint.com/c_standard_library/c_function_qsort.htm
int cmpfunc (const void * a, const void * b) {
   return ( *(double*)a - *(double*)b );
}

int main() {

    FILE* stream = fopen("winequality-red-large.csv", "r");
    char line[1024];
    
    int lineNumber = 0;
    int maxNumberOfElements = 10;

    //+1 to cater for first line with feature headings
    int maxNumLinesRead = maxNumberOfElements + 1;
    //Percentage of data placed in reference set
    double percentageReferenceSet = 0.8;

    numColumns = 12;
    //dimension excludes output column
    dimension = numColumns-1;
    numberOfElementsInReferenceSet = (int)(percentageReferenceSet * maxNumberOfElements);
    numberOfElementsInQuerySet = maxNumberOfElements - numberOfElementsInReferenceSet;

    record_t** qiArray= (record_t**)malloc(sizeof(record_t*)*numberOfElementsInQuerySet);
    record_t** pjArray=(record_t**)malloc(sizeof(record_t*)*numberOfElementsInReferenceSet);
    
    size_t count = 0;

    const char delimeter[2] = ";";
    
   
    int indexReferenceArray = 0;
    int indexQueryArray = 0;

    //https://stackoverflow.com/questions/12911299/read-csv-file-in-c
    while (fgets(line, 1024, stream) && (lineNumber < maxNumLinesRead))
    {
        if(lineNumber > 0){
            char* tmp = strdup(line);
            int column = 0;
            record_t* record =(record_t*)malloc(sizeof(record_t)*dimension);
            double* element = (double*)malloc(sizeof(double)*dimension);
            char* token;
            token = strtok(tmp, delimeter);
            // https://www.tutorialspoint.com/c_standard_library/c_function_strtok.htm
            while( token != NULL ) {
                if(lineNumber > 0){
                    if(column == numColumns-1){
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
            free(token);
           // free(record);
        }
        lineNumber = lineNumber+1;
    }
    printf("%s %d \n","Reference Set Size:",numberOfElementsInReferenceSet );
    printStruct(pjArray,numberOfElementsInReferenceSet,dimension);
    printf("%s %d \n","Query Set Size:",numberOfElementsInQuerySet);
    printStruct(qiArray,numberOfElementsInQuerySet,dimension);

     double** array = getDistances(qiArray,pjArray);
     printf("%s","Distance array: ");
     printArray(array,numberOfElementsInQuerySet,numberOfElementsInReferenceSet);
     quickSort(array,numberOfElementsInQuerySet);
     printf("%s","Sorted distance array: ");
     printArray(array,numberOfElementsInQuerySet,numberOfElementsInReferenceSet);
    
     //findNeighbours(array,3);
     return 0;
}

void printArray( double ** array, int numElements,int dimension){
    printf( "[");
    for(int x = 0; x < numElements; x++){
        printf( "[" );
        for(int y = 0; y < dimension; y++){
            printf( "%lf", array[x][y] );
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
void quickSort(double ** array, int arraySize){
    for(int index = 0 ; index < arraySize; index++){
        qsort(array[index], numberOfElementsInReferenceSet, sizeof(double*), cmpfunc);
    }
}

double** getDistances(record_t** qiArray, record_t ** pjArray){
    double** euclideanDistanceArray= (double**)malloc(sizeof(double*)*numberOfElementsInReferenceSet);
    for(int index = 0 ; index < numberOfElementsInQuerySet; index++){
        euclideanDistanceArray[index] = getDistancePerQueryPoint(qiArray[index],pjArray,numberOfElementsInReferenceSet);
    }
    return euclideanDistanceArray;
}

double* getDistancePerQueryPoint(record_t * qi, record_t ** pjArray, int pjArrayLength){

    double *distanceArray = (double*)malloc(sizeof(double)*pjArrayLength);
    for(int index = 0 ; index < pjArrayLength; index++){
        distanceArray[index] = euclideanDistance(qi->inputs,pjArray[index]->inputs);
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




