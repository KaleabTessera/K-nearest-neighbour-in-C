#include <stdio.h>
#include <math.h>

double * getDistances(int qi, double * pj){
    return pj;
}
    
void printArray( double * pj, int numElements){
    for(int index = 0; index < numElements; index++){
        printf( "%lf", pj[index] );
        printf( "\n" );
    }
}

double euclideanDistance(int dimention, double * qi, double * pj){
    int distance = 0;
    for(int index = 0; index < dimention; index++){
            distance += powf((qi[index] - pj[index]),2);
    }
    return sqrt(distance);
}

int main() {
    int numElements;
    numElements = 100;

    double pj[numElements];
    double qi2[2] = {5,7};
    double pj2[2] = {7,10};
    for(int i = 0; i < numElements; i++){
        pj[i] = 2;
    }
    double distance = euclideanDistance(2,qi2,pj2);
    printf("%lf", distance );  
    //printArray(pj,numElements) ;
    return 0;
}

