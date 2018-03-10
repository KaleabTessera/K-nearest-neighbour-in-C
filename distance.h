#ifndef DISTANCE_H_  
#define DISTANCE_H_
#include "common.h" 
#include <stdlib.h>
#include <math.h>

DistanceFormulae distanceFormulae;
SortFormulae sortFormulae;
int dimension;

double manhattanDistance(double * qi, double * pj);
double euclideanDistance(double * qi, double * pj);
double** getDistancePerQueryPoint(record_t * qi, record_t ** pjArray, int pjArrayLength);


double*** getDistancesParallel(record_t** qiArray, record_t ** pjArray,int x, int y, int dimension){
    dimension = dimension;
    double*** euclideanDistanceArray= (double***)malloc(sizeof(double**)*y);
    #pragma omp parallel
    #pragma omp for
    for(int index = 0 ; index < x; index++){
        euclideanDistanceArray[index] = getDistancePerQueryPoint(qiArray[index],pjArray,y);
    }
    return euclideanDistanceArray;
}

double*** getDistancesSerial(record_t** qiArray, record_t ** pjArray,int x, int y, int dimension){
    dimension = dimension;
    double*** euclideanDistanceArray= (double***)malloc(sizeof(double**)*y);
    for(int index = 0 ; index < x; index++){
        euclideanDistanceArray[index] = getDistancePerQueryPoint(qiArray[index],pjArray,y);
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
#endif 