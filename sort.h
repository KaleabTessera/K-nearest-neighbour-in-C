#ifndef _sort_h_  
#define _sort_h_
#include "common.h" 

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>
#include "common.h" 

//https://www.tutorialspoint.com/c_standard_library/c_function_qsort.htm
//https://stackoverflow.com/questions/6103636/c-qsort-not-working-correctly

int cmpfunc (const void ** a, const void ** b) ;

int cmpfunc (const void ** a, const void ** b) {
    double i = *(double*)a[0];
    double r = *(double*)b[0];
    return (i > r) - (i < r);
}

void sortElements(double *** array, int arraySize, int numberOfElementsInReferenceSet,SortFormulae sortFormulae ){
    if(sortFormulae == QuickSort){
        for(int index = 0 ; index < arraySize; index++){
            qsort(array[index], 10, sizeof(double**), cmpfunc);
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
#endif 