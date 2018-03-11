#ifndef _sort_h_  
#define _sort_h_
#include "common.h" 

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>
#include "common.h" 

#include <string.h>

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
    else if( sortFormulae == RadixSort){
        for(int index = 0 ; index < arraySize; index++){
            radixsort(array[index], numberOfElementsInReferenceSet);
        }
    }
    else if(sortFormulae == SelectionSort){
        for(int index = 0 ; index < arraySize; index++){
            selectionSort(array[index], numberOfElementsInReferenceSet);
        }
    }
}

void sortElementsParallel(double *** array, int arraySize, int numberOfElementsInReferenceSet,SortFormulae sortFormulae ){
    if(sortFormulae == QuickSort){
        for(int index = 0 ; index < arraySize; index++){
            qsort(array[index], 10, sizeof(double**), cmpfunc);
        }
    }
    else if( sortFormulae == RadixSort){
        for(int index = 0 ; index < arraySize; index++){
            radixsortParallel(array[index], numberOfElementsInReferenceSet);
        }
    }
    else if(sortFormulae == SelectionSort){
        for(int index = 0 ; index < arraySize; index++){
            selectionSort(array[index], numberOfElementsInReferenceSet);
        }
    }
}

//https://www.geeksforgeeks.org/radix-sort/
// A utility function to get maximum value in arr[]
double getMax(double** arr, int n)
{
    double mx = arr[0][0];
    for (int i = 1; i < n; i++)
        if (arr[i][0] > mx)
            mx = arr[i][0];
    return mx;
}
 
void countSortParallel(double** arr, int n, int exp)
{
    double output[n]; // output array
    int i, count[10] = {0},count1[10] = {0},count2[10] = {0},count3[10] = {0},count4[10] = {0};
    int numSplits = 4;
    double **split1 = (double **)malloc(n/numSplits * sizeof(double *));
    double **split2 = (double **)malloc(n/numSplits * sizeof(double *));
    double **split3 = (double **)malloc(n/numSplits * sizeof(double *));
    double **split4 = (double **)malloc(n/numSplits * sizeof(double *));
    memcpy(split1, arr                          , n/numSplits * sizeof(double *));
    memcpy(split2, arr+(n/numSplits)            ,n/numSplits* sizeof(double *));
    memcpy(split3, arr +(2 * n/numSplits)       , n/numSplits * sizeof(double *));
    memcpy(split4, arr+ (3 * n/numSplits)       ,n/numSplits * sizeof(double *));
    // Store count of occurrences in count[]
    
    #pragma omp parallel sections num_threads(4)
    {
        #pragma omp section
    for (int i = 0; i < n/numSplits; i++)
        count1[((int)split1[i][0]/exp)%10 ]++;
        #pragma omp section
     for (int i = 0; i < n/numSplits; i++)
        count2[((int)split2[i][0]/exp)%10 ]++;
        #pragma omp section
    for (int i = 0; i < n/numSplits; i++)
        count3[((int)split3[i][0]/exp)%10 ]++;
        #pragma omp section
     for (int i = 0; i < n/numSplits; i++)
        count4[((int)split4[i][0]/exp)%10 ]++;
    }

    for (i = 0; i < 10; i++)
        count[i] =  count1[i] + count2[i] + count3[i]+ count4[i];


    for (i = 1; i < 10; i++)
        count[i] += count[i - 1];
 
    for (int i = n - 1; i >= 0; i--)
    {
        output[count[ ((int)arr[i][0]/exp)%10 ] - 1] = arr[i][0];
        count[ ((int)arr[i][0]/exp)%10 ]--;
    }

    for (i = 0; i < n; i++)
        arr[i][0] = (int)output[i];
}

void countSort(double** arr, int n, int exp)
{
    double output[n]; // output array
    int i, count[10] = {0};
 
    // Store count of occurrences in count[]
    for (i = 0; i < n; i++)
        count[((int)arr[i][0]/exp)%10 ]++;
 
    // Change count[i] so that count[i] now contains actual
    //  position of this digit in output[]
    for (i = 1; i < 10; i++)
        count[i] += count[i - 1];
 
    // Build the output array
    for (i = n - 1; i >= 0; i--)
    {
        output[count[ ((int)arr[i][0]/exp)%10 ] - 1] = arr[i][0];
        count[ ((int)arr[i][0]/exp)%10 ]--;
    }
 
    // Copy the output array to arr[], so that arr[] now
    // contains sorted numbers according to current digit
    for (i = 0; i < n; i++)
        arr[i][0] = (int)output[i];
}

//https://www.geeksforgeeks.org/radix-sort/

void radixsort(double** arr, int n)
{
    // Find the maximum number to know number of digits
    double m = getMax(arr, n);
 
    // Do counting sort for every digit. Note that instead
    // of passing digit number, exp is passed. exp is 10^i
    // where i is current digit number
    
    for (int exp = 1; m/exp > 0; exp *= 10)
        countSort(arr, n, exp);
}

void radixsortParallel(double** arr, int n)
{
    // Find the maximum number to know number of digits
    double m = getMax(arr, n);
 
    // Do counting sort for every digit. Note that instead
    // of passing digit number, exp is passed. exp is 10^i
    // where i is current digit number
    
    for (int exp = 1; m/exp > 0; exp *= 10)
        countSortParallel(arr, n, exp);
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