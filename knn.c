#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>
#include "common.h" 
#include <stdbool.h>
#include "knn.h"

int dimension = 0;
int numColumns = 0;
int numberOfElementsInQuerySet = 0;
int numberOfElementsInReferenceSet = 0;
int k;
DistanceFormulae distanceFormulae;
SortFormulae sortFormulae;
FILE *fileToWriteSort;
FILE *fileToWriteDistance;

void KNN(record_t** qiArray, record_t ** pjArray,int numberOfElementsInQuerySet, int numberOfElementsInReferenceSet, int k,FILE *fileToWriteTo,int numMaxElements,int dimension, int numberOfTimesRepeated);
double*** getDistances(record_t** qiArray, record_t ** pjArray,int x, int y, int dimension);
void sortElements(double *** array, int arraySize, int numberOfElementsInReferenceSet,SortFormulae sortFormulae );

void printArray( double ** array,int dimension);
void printArray3d( double *** array, int numElements,int dimension);

int main() {

    int numMaxElementsConfig = 3;
    int numDimensionConfig = 3;
    int numMaxElements[] = {100,250,500};
    int dimensions[] = {32,40,128};

    int numberOfTimesRepeated = 20;
    
    
  
    FILE* stream = fopen("mnist_train.csv", "r");
    writeToFile();
    for(int total = 0; total < numMaxElementsConfig; total++){
        for(int d =0; d< numDimensionConfig; d++){
            for(int i=0; i<numberOfTimesRepeated; i++ ){
                //fprintf(fileToWrite,"\n %s \n","***********************************************************************************");
                //fprintf(fileToWrite,"NumElements: %d NumDimensions: %d \n",numMaxElements[total], dimensions[d] );
                //fprintf(fileToWrite,"Run number: %d \n", i );
                char line[4096];

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
                while (fgets(line, 4096, stream) && (lineNumber < maxNumLinesRead))
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
                                double numToAdd = strtod (token, NULL);
                                if(isnan(numToAdd)){
                                    numToAdd=0;
                                }
                                if(column == indexOfOutput){
                                    record->output = numToAdd;    
                                }
                                else record->inputs[column] =numToAdd;
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
                    lineNumber = lineNumber+1;

                }
                //fprintf(fileToWrite,"%s %d \n","Reference Set Size:",numberOfElementsInReferenceSet );
                //printStruct(pjArray,numberOfElementsInReferenceSet,dimension);
                //fprintf(fileToWrite,"%s %d \n \n","Query Set Size:",numberOfElementsInQuerySet);
                //printStruct(qiArray,numberOfElementsInQuerySet,dimension);

                //Only write to file on last run
                if(i == numberOfTimesRepeated -1){
                    KNN(qiArray,pjArray,numberOfElementsInQuerySet,numberOfElementsInReferenceSet, k,fileToWriteDistance,numMaxElements[total], dimensions[d],numberOfTimesRepeated);
                }
                else  KNN(qiArray,pjArray,numberOfElementsInQuerySet,numberOfElementsInReferenceSet, k,NULL,numMaxElements[total], dimensions[d],numberOfTimesRepeated);
                }
        }
    }

    return 0;
}

void writeToFile(){
    //Writing to file
    time_t rawtime;
    struct tm * timeinfo;

    time ( &rawtime );
    timeinfo = localtime ( &rawtime );

    //distance results
    char fileNameDistance[100] ;
    strcat(fileNameDistance, "./results/distance_");
    strcat(fileNameDistance, asctime (timeinfo));
    strcat(fileNameDistance, ".txt");
    fileToWriteDistance = fopen(fileNameDistance, "w");
    if (fileToWriteDistance == NULL)
    {
        fileToWriteDistance = fopen("/results/distance_test.txt", "w");
    }

    //sort results
    char fileNameSort[2000] ;
    strcat(fileNameSort, "./results/sort_");
    strcat(fileNameSort, asctime (timeinfo));
    strcat(fileNameSort, ".txt");
    fileToWriteSort = fopen(fileNameSort, "w");
    if (fileToWriteDistance == NULL)
    {
        fileToWriteDistance = fopen("/results/sort_test.txt", "w");
    }
}

