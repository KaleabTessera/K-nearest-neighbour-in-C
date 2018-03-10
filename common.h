#ifndef _common_h   
#define _common_h
typedef enum {Euclidean, Manhattan} DistanceFormulae;
typedef enum {InsertionSort,QuickSort , SelectionSort } SortFormulae;
typedef enum {Serial, Parallel } Computation;

typedef struct Element {
    double output;
    double inputs[];
} record_t;

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
        
        printf("Inputs: [" );
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
#endif 