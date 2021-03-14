#include <stdio.h>
#include <time.h>
#include <windows.h>
#include <winbase.h>
#include <psapi.h>
#include <Windows.h>
#include <exception>
#include <string.h>
#include <iostream>
#include <fstream>
using namespace std;

HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
double get_memory_used_MB()
{
    PROCESS_MEMORY_COUNTERS info;
    GetProcessMemoryInfo(GetCurrentProcess(), &info, sizeof(info));
    return (double)info.PeakWorkingSetSize / (1024 * 1024);
}
typedef struct array_list
{
    double* data;
    int length;
    int block;
} array_list;
#define block_size 100
array_list create_array_list()
{
    array_list ret = { (double*)malloc(block_size * sizeof(double)), 0, block_size, };
    return ret;
}

void increment(array_list* a)
{
    (*a).length++;
    if ((*a).length > (*a).block)
    {
        (*a).block += block_size;
        (*a).data = (double*)realloc((*a).data, (*a).block * sizeof(double));
    }
}
void insert(array_list* a, int p, double value)
{
    increment(a);
    for (int i = (*a).length;i >= p;i--)
        (*a).data[i + 1] = (*a).data[i];
    (*a).data[p] = value;

}
void add(array_list* a, double value)
{
    increment(a);
    (*a).data[(*a).length - 1] = value;
}
int append(array_list* list, double value)
{
    int l = 0;
    int r = (*list).length - 1;

    int c = 0;

    double com = 0;
    while (l <= r)
    {
        c = (l + r) / 2;
        com = (*list).data[c] - value;
        if (com > 0) r = c - 1;
        else if (com < 0) l = c + 1;
        else
        {
            insert(list, c, value);
            return c;
        }

    }
    if (c >= (*list).length)
        add(list, value);

    else
    {
        if (com < 0) c++;
        insert(list, c, value);
    }

    return c;
}
int find(array_list* list, double value)
{
    int l = 0;
    int r = (*list).length - 1;

    int c = 0;

    double com = 0;
    while (l <= r)
    {
        c = (l + r) / 2;
        com = (*list).data[c] - value;
        if (com > 0) r = c - 1;
        else if (com < 0) l = c + 1;
        else
        {
            return c;
        }

    }
    return -1;
}

//int trocas = 0;
//int comparacoes = 0;
// A utility function to swap two elements
void swap(int* a, int* b)
{
    int t = *a;
    *a = *b;
    *b = t;
}

/* This function takes last element as pivot, places
   the pivot element at its correct position in sorted
    array, and places all smaller (smaller than pivot)
   to left of pivot and all greater elements to right
   of pivot */
int partition(int arr[], int low, int high)
{
    int pivot = arr[high];    // pivot
    int i = (low - 1);  // Index of smaller element

    for (int j = low; j <= high - 1; j++)
    {
        // If current element is smaller than or
        // equal to pivot
        //comparacoes++;
        if (arr[j] <= pivot)
        {
            i++;    // increment index of smaller element
            swap(&arr[i], &arr[j]);
        }
    }
    swap(&arr[i + 1], &arr[high]);
    return (i + 1);
}
/* Function to print an array */
void printArray(int arr[], int size)
{
    int i;
    for (i = 0; i < size; i++)
        printf("%d ", arr[i]);
    printf("\n");
}

void quickSortPart(int vetor[], int inicio, int fim) {

    int pivo, aux, i, j, meio;

    i = inicio;
    j = fim;

    meio = (int)((i + j) / 2);
    pivo = vetor[meio];

    do {
        while (vetor[i] < pivo) i = i + 1;
        while (vetor[j] > pivo) j = j - 1;

        if (i <= j) {
            aux = vetor[i];
            vetor[i] = vetor[j];
            vetor[j] = aux;
            i = i + 1;
            j = j - 1;
        }
    } while (j > i);

    if (inicio < j) quickSortPart(vetor, inicio, j);
    if (i < fim) quickSortPart(vetor, i, fim);

}
void quickSort(int arr[], int n)
{
    quickSortPart(arr, 0, n - 1);
}
void selectionSort(int arr[], int n)
{
    int i, j, min_idx;

    // One by one move boundary of unsorted subarray  
    for (i = 0; i < n - 1; i++)
    {
        // Find the minimum element in unsorted array  
        min_idx = i;
        for (j = i + 1; j < n; j++)
        {
            if (arr[j] < arr[min_idx])
            {
                min_idx = j;
            }
        }
            

        // Swap the found minimum element with the first element  
        swap(&arr[min_idx], &arr[i]);
    }
}
void insertionSort(int arr[], int n)
{
    int i, key, j;
    for (i = 1; i < n; i++)
    {
        key = arr[i];
        j = i - 1;

        /* Move elements of arr[0..i-1], that are
        greater than key, to one position ahead
        of their current position */
        while (j >= 0 && arr[j] > key)
        {
            arr[j + 1] = arr[j];
            j = j - 1;
        }
        arr[j + 1] = key;
    }
}
void bubbleSort(int arr[], int n)
{
    int i, j;
    for (i = 0; i < n - 1; i++)

        // Last i elements are already in place  
        for (j = 0; j < n - i - 1; j++)
        {
            if (arr[j] > arr[j + 1])
            {
                swap(&arr[j], &arr[j + 1]);
            }
        }
    
}
void heapify(int arr[], int n, int i)
{
    int largest = i; // Initialize largest as root
    int l = 2 * i + 1; // left = 2*i + 1
    int r = 2 * i + 2; // right = 2*i + 2

    // If left child is larger than root
    //comparacoes += 2;
    if (l < n && arr[l] > arr[largest])
    {
        largest = l;
        
    }

    // If right child is larger than largest so far

    if (r < n && arr[r] > arr[largest])
    {
        largest = r;
    }

    // If largest is not root
   // comparacoes++;
    if (largest != i) {
        swap(&arr[i], &arr[largest]);

        // Recursively heapify the affected sub-tree
        heapify(arr, n, largest);
    }
}

// main function to do heap sort
void heapSort(int arr[], int n)
{
    // Build heap (rearrange array)
    for (int i = n / 2 - 1; i >= 0; i--)
        heapify(arr, n, i);

    // One by one extract an element from heap
    for (int i = n - 1; i > 0; i--) {
        // Move current root to end
        swap(&arr[0], &arr[i]);

        // call max heapify on the reduced heap
        heapify(arr, i, 0);
    }
}
void merge(int arr[], int l, int m, int r)
{
    int n1 = m - l + 1;
    int n2 = r - m;

    // Create temp arrays
    int* L = (int*)malloc(n1 * sizeof(int));
    int* R= (int*)malloc(n2 * sizeof(int));

    // Copy data to temp arrays L[] and R[]
    for (int i = 0; i < n1; i++)
        L[i] = arr[l + i];
    for (int j = 0; j < n2; j++)
        R[j] = arr[m + 1 + j];

    // Merge the temp arrays back into arr[l..r]

    // Initial index of first subarray
    int i = 0;

    // Initial index of second subarray
    int j = 0;

    // Initial index of merged subarray
    int k = l;

    while (i < n1 && j < n2) {
        //comparacoes++;
        if (L[i] <= R[j]) {
            arr[k] = L[i];
            i++;
        }
        else {
            arr[k] = R[j];
            j++;
        }
        k++;
    }

    // Copy the remaining elements of
    // L[], if there are any
    while (i < n1) {
        arr[k] = L[i];
        i++;
        k++;
    }

    // Copy the remaining elements of
    // R[], if there are any
    while (j < n2) {
        arr[k] = R[j];
        j++;
        k++;
    }
}

// l is for left index and r is
// right index of the sub-array
// of arr to be sorted */
void mergeSortPart(int arr[], int l, int r) {
    //comparacoes++;
    if (l >= r) {
        return;//returns recursively
    }
    int m = l + (r - l) / 2;
    mergeSortPart(arr, l, m);
    mergeSortPart(arr, m + 1, r);
    merge(arr, l, m, r);
}
void mergeSort(int arr[], int l)
{
    mergeSortPart(arr, 0, l - 1);
}
// A utility function to get maximum value in arr[] 
int getMax(int arr[], int n)
{
    int mx = arr[0];
    for (int i = 1; i < n; i++)
    {
        if (arr[i] > mx)
            mx = arr[i];
    }
        
    return mx;
}

// A function to do counting sort of arr[] according to 
// the digit represented by exp. 
void countSort(int arr[], int n, int exp)
{
    int * output=(int*)malloc(n*sizeof(int)); // output array 
    int i, count[10] = { 0 };

    // Store count of occurrences in count[] 
    for (i = 0; i < n; i++)
        count[(arr[i] / exp) % 10]++;

    // Change count[i] so that count[i] now contains actual 
    //  position of this digit in output[] 
    for (i = 1; i < 10; i++)
        count[i] += count[i - 1];

    // Build the output array 
    for (i = n - 1; i >= 0; i--) {
        output[count[(arr[i] / exp) % 10] - 1] = arr[i];
        count[(arr[i] / exp) % 10]--;
    }

    // Copy the output array to arr[], so that arr[] now 
    // contains sorted numbers according to current digit 
    for (i = 0; i < n; i++)
        arr[i] = output[i];
}

// The main function to that sorts arr[] of size n using 
// Radix Sort 
void radixSort(int arr[], int n)
{
    // Find the maximum number to know number of digits 
    int m = getMax(arr, n);

    // Do counting sort for every digit. Note that instead 
    // of passing digit number, exp is passed. exp is 10^i 
    // where i is current digit number 
    for (int exp = 1; m / exp > 0; exp *= 10)
        countSort(arr, n, exp);
}
//metodo criado combinando dois algoritimos
void insert(int arr[], int n,int p, int v)
{
    for (int i = n - 1; i >= p;i--)
    {
        arr[i+1] = arr[i];
    }
    arr[p] = v;
}
int append_vet(int arr[],int n, int value)
{
    int l = 0;
    int r = n;

    int c = 0;

    int com = 0;
    while (l <= r)
    {
        c = (l + r) / 2;
        com = arr[c] - value;
        if (com > 0) r = c - 1;
        else if (com < 0) l = c + 1;
        else
        {
            insert(arr, n,c, value);
            return c;
        }

    }
        if (com < 0) c++;
        insert(arr,n,c, value);

    return c;
}

void betterInsert(int arr[],int n)
{
    for (size_t i = 0; i < n; i++)
    {
        append_vet(arr,i,arr[i]);
    }
}

void fill(int* v, int l)
{
    for (size_t i = 0; i < l;i++)
    {
        v[i] = rand() ;
    }
}
//shell sort
void shellSort(int arr[], int n)
{
    // Start with a big gap, then reduce the gap 
    for (int gap = n / 2; gap > 0; gap /= 2)
    {
        // Do a gapped insertion sort for this gap size. 
        // The first gap elements a[0..gap-1] are already in gapped order 
        // keep adding one more element until the entire array is 
        // gap sorted  
        for (int i = gap; i < n; i += 1)
        {
            // add a[i] to the elements that have been gap sorted 
            // save a[i] in temp and make a hole at position i 
            int temp = arr[i];

            // shift earlier gap-sorted elements up until the correct  
            // location for a[i] is found 
            int j;
            for (j = i; j >= gap && arr[j - gap] > temp; j -= gap)
                arr[j] = arr[j - gap];

            //  put temp (the original a[i]) in its correct location 
            arr[j] = temp;
           // trocas++;
            
        }
    }
}


// l is for left index and r is
// right index of the sub-array
// of arr to be sorted */
void betterInsertPart(int arr[], int l,int r)
{
    for (size_t i = l; i < r; i++)
    {
        append_vet(arr, i, arr[i]);
        //trocas++;
    }
}
void mergeInsertPart(int arr[], int l, int r) {
    //comparacoes++;
    if (l >= r) {
        return;//returns recursively
    }
    int m = l + (r - l) / 2;
    if (m-l<=2||r-m<=2)
    {
        betterInsertPart(arr, l,r);
    }
    else
    {
        mergeInsertPart(arr, l, m);
        mergeInsertPart(arr, m + 1, r);
        merge(arr, l, m, r);
    }
    
}
void mergeInsert(int arr[], int l)
{
    mergeSortPart(arr, 0, l - 1);
}
typedef void (*SortAlgorithim) (int arr[], int n);
SortAlgorithim function[] = { &quickSort, &selectionSort, &insertionSort, &bubbleSort, &heapSort, &mergeSort, &shellSort, &radixSort,&betterInsert ,&mergeInsert };
char name[][100] = { "Quick","Selection","Insertion","Bubble","Heap","Merge","ShellSort","Radix","BetterInsert","MergeInsert" };

void testAlgorithim(SortAlgorithim s, char * name, int arr[], int n,int d,ofstream * toutput, ofstream* moutput)
{
    
    double time = 0;
    double memory = 0;
    array_list times = create_array_list();
    for (size_t j = 0; j < 5; j++)
    {
        int* c = (int*)malloc(n * sizeof(int));
        for (int k = 0; k < n; k++)
            c[k] = arr[k];
        clock_t start_time, end_time;
        size_t peakSize;
        start_time = clock();

        s(c, n);

        end_time = clock();
        memory = get_memory_used_MB();
        append(&times, (end_time - (double)start_time) / CLOCKS_PER_SEC);
        free(c);

    }
    for (size_t j = 1; j < 4; j++)
    {
        time += times.data[j];
    }
    time /= 3;
    free(times.data);

    


    printf("Algoritimo: %s, Tempo: %f segundos, Memoria: %f MB\n", name, time, memory);
    *toutput << "\"" << time << "\", ";
    *moutput << "\"" << memory << "\", ";
}
ofstream trandom;
ofstream tasc;
ofstream tdesc;
ofstream mrandom;
ofstream masc;
ofstream mdesc;

void handle() {};
void testAlgorithims(int arr[],int n,int d, ofstream * toutput,ofstream* moutput)
{
    *toutput << n <<", ";
    *moutput << n <<", ";
    for (size_t j = 0; j < sizeof(name)/sizeof(name[0]); j++)
    {
        testAlgorithim(function[j], name[j], arr, n,d, toutput,moutput);
    }
    *toutput << "\n";
    *moutput << "\n";

}
// Driver program to test above functions
int main()
{
    srand(time(NULL));
    trandom.open("time/random.csv");
    tasc.open("time/asc.csv");
    tdesc.open("time/desc.csv");
    mrandom.open("memory/random.csv");
    masc.open("memory/asc.csv");
    mdesc.open("memory/desc.csv");
    trandom << "N, ";
    tasc << "N, ";
    tdesc << "N, ";
    mrandom << "N, ";
    masc << "N, ";
    mdesc << "N, ";
    for (size_t i = 0; i < sizeof(name) / sizeof(name[0]); i++)
    {
        trandom << name[i] << ", ";
        tasc << name[i] << ", ";
        tdesc << name[i] << ", ";
        mrandom << name[i] << ", ";
        masc << name[i] << ", ";
        mdesc << name[i] << ", ";
    }
    trandom << "\n";
    tasc << "\n";
    tdesc << "\n";
    mrandom << "\n";
    masc << "\n";
    mdesc << "\n";

    for (size_t i = 2000; i <= 128000; i*=2)
    {
        int *a =(int*)malloc(i*sizeof(int));
        fill(a, i);


        SetConsoleTextAttribute(hConsole, 11);
        printf("N = %d\n",i);
        printf("Random\n");
        SetConsoleTextAttribute(hConsole, 7);
        testAlgorithims(a,i,0, &trandom,&mrandom);
        SetConsoleTextAttribute(hConsole, 11);

        radixSort(a, i);
        printf("ASC\n");

        SetConsoleTextAttribute(hConsole, 7);
        testAlgorithims(a, i,1,&tasc, &masc);
        

        
        SetConsoleTextAttribute(hConsole, 11);
        int* b = (int*)malloc(i * sizeof(int));
        for (int j = 0; j < i; j++)
        {
            b[j] = a[i - j - 1];
        }
        printf("DESC\n");
        SetConsoleTextAttribute(hConsole, 7);
        testAlgorithims(b, i,2, &tdesc, &mdesc);

        free(a);
    }
    trandom.close();
    tasc.close();
    tdesc.close();
    mrandom.close();
    masc.close();
    mdesc.close();

    return 0;
}