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
int main() {
    
    numElementsQueryArray = 2;
    numElementsReferenceArray = 2;
    dimension = 2; 

    double * qiArray[numElementsQueryArray];
    double * pjArray[numElementsReferenceArray];

    double qi[2] = {5,6};
    double qi2[2] = {7,8};

    double pj[2] = {10,1.5};
    double pj2[2] = {11,2.5};

    qiArray[0] = qi;
    qiArray[1] = qi2;
    pjArray[0] = pj;
    pjArray[1] = pj2;

    double** array = getDistances(qiArray,pjArray);
     printArray(qiArray,2,2);
     printArray(pjArray,2,2);
     printArray(array,2,2);
    
    return 0;
}

void printArray( double * array[], int numElements,int dimension){
    printf( "[");
    for(int x = 0; x < numElements; x++){
        printf( "[" );
        for(int y = 0; y < dimension; y++){
            printf( "%lf", array[x][y] );
            if(x != numElements -1){
                printf(",");
            }
        }
        printf( "]" );
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




