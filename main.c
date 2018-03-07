#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int dimension = 0;
int numElementsQueryArray = 0;
int numElementsReferenceArray = 0;

void printArray( double * pj[], int numElements,int dimension);
double** getDistances(double* qiArray[], double * pjArray[]);
double* getDistancePerQueryPoint(double * qi, double * pjArray[], int pjArrayLength);
double euclideanDistance(double * qi, double * pj);

enum ElementType { et_str, et_int, et_dbl };
typedef struct Element {
    double output;
    enum ElementType type;
    float inputs[10];
} record_t;
    
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

int main() {

   
    double** qiArray= (double**)malloc(sizeof(double*)*1024);
    double ** pjArray=(double**)malloc(sizeof(double*)*1024);
    record_t records[100];
    size_t count = 0;

    FILE* stream = fopen("winequality-red-large.csv", "r");
    char line[1024];
    
    int lineNumber = 0;
    int maxNumberOfElements = 1000;

    //+1 to cater for first line with feature headings
    int maxNumLinesRead = maxNumberOfElements + 1;

    const char delimeter[2] = ";";
    //Percentage of data placed in reference set
    double percentageReferenceSet = 0.8;
    int numberOfElementsInReferenceSet = (int)(percentageReferenceSet * maxNumberOfElements);
    int numberOfElementsInQuerySet = maxNumberOfElements - numberOfElementsInReferenceSet;
    int indexReferenceArray = 0;
    int indexQueryArray = 0;

    //https://stackoverflow.com/questions/12911299/read-csv-file-in-c
    while (fgets(line, 1024, stream) && (lineNumber < maxNumLinesRead))
    {
        if(lineNumber > 0){
        char* tmp = strdup(line);
        int column = 0;
        double* element = (double*)malloc(sizeof(double)*12);
        char* token;
        token = strtok(tmp, delimeter);
        // https://www.tutorialspoint.com/c_standard_library/c_function_strtok.htm
        while( token != NULL ) {
            if(lineNumber > 0){
                element[column] = strtod (token, NULL); 
            }
            column=column+1;
            token = strtok(NULL,delimeter);
        }
        
        if(lineNumber <= numberOfElementsInReferenceSet){
            pjArray[indexReferenceArray] = element;
            indexReferenceArray = indexReferenceArray+1;
        }
        else {
            qiArray[indexQueryArray] = element;
            indexQueryArray = indexQueryArray+1;
        } 
        free(tmp);
        free(token);
        }
        lineNumber = lineNumber+1;
    }
    printf("%s %d \n","Reference Set Size:",numberOfElementsInReferenceSet );
    printArray(pjArray,numberOfElementsInReferenceSet,12);
    printf("%s %d \n","Query Set",numberOfElementsInQuerySet);
    printArray(qiArray,numberOfElementsInQuerySet,12);

     //double** array = getDistances(qiArray,pjArray);
     //printArray(array,2,2);
    
     return 0;
}

void printArray( double * array[], int numElements,int dimension){
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

double** getDistances(double* qiArray[], double * pjArray[]){
    double** euclideanDistanceArray= (double**)malloc(sizeof(double*)*numElementsReferenceArray);
    for(int index = 0 ; index < numElementsQueryArray; index++){
        euclideanDistanceArray[index] = getDistancePerQueryPoint(qiArray[index],pjArray,numElementsReferenceArray);
    }
    return euclideanDistanceArray;
}

double* getDistancePerQueryPoint(double * qi, double * pjArray[], int pjArrayLength){

    double *distanceArray = (double*)malloc(sizeof(double)*pjArrayLength);
    for(int index = 0 ; index < pjArrayLength; index++){
        distanceArray[index] = euclideanDistance(qi,pjArray[index]);
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




