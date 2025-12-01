#include <iostream>
#include <fstream>
#include <atomic>
#include <vector>
#include <algorithm>
#include <functional> 
#include <chrono>     
#include <string>     
#include "pcm/src/cpucounters.h"
#include <bits/stdc++.h>
#include <omp.h>     

using namespace pcm;

std::vector<CoreCounterState> cstates1, cstates2;
std::vector<SocketCounterState> sktstate1, sktstate2;
SystemCounterState sstate1, sstate2;
int EXECS_PER_ALG = 10;
const int RUN = 32;
std::string ARRAY_FILE = "medio.bin";
std::string CSV_FILE = "./resultados.csv";

/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////


void bubbleSortMain(std::vector<int>& arr, int n) {
    bool swapped = false;
    for (int i = 0; i < n - 1; ++i){
        swapped = false;
        for (int j = 0; j < n - 1 - i; ++j){ //o último elemento vai estar certo, não precisa chegar até ele
            if(arr[j] > arr[j+1]){
                std::swap(arr[j], arr[j+1]);
                swapped = true;
            }
        }
        if(!swapped)
            break;
    }
}

/////////////////////////////////

static void heapify(std::vector<int>& arr, int n, int i){
    while (true) {
        int largest = i;

        int l = 2 * i + 1;
        int r = 2 * i + 2;
        if (l < n && arr[l] > arr[largest])
            largest = l;

        if (r < n && arr[r] > arr[largest])
            largest = r;

        if (largest == i)
            break; 

        std::swap(arr[i], arr[largest]);
        i = largest; 
    }
}

void heapSortMain(std::vector<int>& arr, int n){

    for (int i = n / 2 - 1; i >= 0; i--)
        heapify(arr, n, i);

    for (int i = n - 1; i > 0; i--) {

        std::swap(arr[0], arr[i]);

        heapify(arr, i, 0);
    }
}

/////////////////////////////////////////////////////////////////

void insertionSortMain(std::vector<int>& arr, int n)
{
    for (int i = 1; i < n; ++i) {
        int key = arr[i];
        int j = i - 1;

        while (j >= 0 && arr[j] > key) {
            arr[j + 1] = arr[j];
            j--;
        }
        arr[j + 1] = key;
    }
}

///////////////////////////////////////////////////////////////////

static void iterativeMerge(std::vector<int>& arr, std::vector<int>& temp, int left, int mid, int right) {
    int i = left;    
    int j = mid + 1; 
    int k = left;    

    for (int idx = left; idx <= right; ++idx) {
        temp[idx] = arr[idx];
    }

    while (i <= mid && j <= right) {
        if (temp[i] <= temp[j]) {
            arr[k++] = temp[i++];
        } else {
            arr[k++] = temp[j++];
        }
    }
    while (i <= mid) {
        arr[k++] = temp[i++];
    }
}

void iterativeMergeSortMain(std::vector<int>& arr, int n) {
    std::vector<int> temp(n);

    for (int currSize = 1; currSize < n; currSize = 2*currSize) {
        for (int leftStart = 0; leftStart < n; leftStart += 2*currSize) {

            int mid = std::min(leftStart + currSize - 1, n-1);
            int rightEnd = std::min(leftStart + 2*currSize - 1, n-1);

            if (mid < rightEnd) {
                iterativeMerge(arr, temp, leftStart, mid, rightEnd);
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////

//Hoare's Partition Algorithm
static int partition(std::vector<int> &arr, int low, int high) {
    int pivot = arr[low];

    int i = low -1, j = high + 1;
    while (true) {

        do {
            i++;
        } while (i <= high && arr[i] < pivot);

        do {
            j--;
        } while (j >= low && arr[j] > pivot);

        if (i >= j) return j;

        std::swap(arr[i], arr[j]);
    }
}

static void midOfThree(std::vector<int>& array, int low, int high){
    int mid = (high + low) / 2;

    if (array[mid] < array[low]) {
        std::swap(array[mid], array[low]);
    }
    if (array[high] < array[low]) {
        std::swap(array[high], array[low]);
    }
    if (array[high] < array[mid]) {
        std::swap(array[high], array[mid]);
    }

    std::swap(array[mid], array[low]);
}

static void quickSortIterative(std::vector<int>& arr, int l, int h)
{
    std::stack<int> s;

    s.push(l);
    s.push(h);

	while (!s.empty()) {
        h = s.top();
        s.pop();
        l = s.top();
        s.pop();


        midOfThree(arr, l, h);

		int p = partition(arr, l, h);
		if (p > l) {
            s.push(l);
            s.push(p);
		}

		if (p + 1 < h) {
            s.push(p + 1);
            s.push(h);
		}
	}
}

void iterativeQuickSortMain(std::vector<int>& arr, int n)
{
	quickSortIterative(arr, 0, n - 1);
}

//////////////////////////////////////////////////////////////////////

static void mergeSortMerge(std::vector<int>& arr, int left, int mid, int right)
{
    int n1 = mid - left + 1;
    int n2 = right - mid;

    std::vector<int> L(n1), R(n2);

    for (int i = 0; i < n1; i++)
        L[i] = arr[left + i];
    for (int j = 0; j < n2; j++)
        R[j] = arr[mid + 1 + j];

    int i = 0, j = 0;
    int k = left;

    while (i < n1 && j < n2) {
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

    while (i < n1) {
        arr[k] = L[i];
        i++;
        k++;
    }

    while (j < n2) {
        arr[k] = R[j];
        j++;
        k++;
    }
}

static void mergeSort(std::vector<int>& arr, int left, int right)
{
    if (left >= right)
        return;

    int mid = left + (right - left) / 2;
    mergeSort(arr, left, mid);
    mergeSort(arr, mid + 1, right);
    mergeSortMerge(arr, left, mid, right);
}

void mergeSortMain(std::vector<int>& arr, int n) {
    mergeSort(arr, 0, n - 1);
}

/////////////////////////////////////////////////////////////////////////

void oddEvenSortMain(std::vector<int>& arr, int n) {
    bool isSorted = false;

    while (!isSorted)
    {
        isSorted = true;

        // Bubble sort nos elementos de index ímpar
        for (int i=1; i<=n-2; i=i+2)
        {
            if (arr[i] > arr[i+1])
            {
                std::swap(arr[i], arr[i+1]);
                isSorted = false;
            }
        }

        // Bubble sort nos elementos de index par
        for (int i=0; i<=n-2; i=i+2)
        {
            if (arr[i] > arr[i+1])
            {
                std::swap(arr[i], arr[i+1]);
                isSorted = false;
            }
        }
    }

    return;
}

////////////////////////////////////////////////////////////////////////


static void quickSort(std::vector<int>& array, int low, int high) {
    if (low >= 0 && high>= 0 && low < high) {

        midOfThree(array, low, high);

        int partition_index = partition(array, low, high);

        quickSort(array, low, partition_index);

        quickSort(array, partition_index + 1, high);
    }
}

void quickSortMain(std::vector<int>& arr, int n) {

    quickSort(arr, 0, n - 1);

}

////////////////////////////////////////////////////////////////////////

static int radixGetMax(std::vector<int>& arr, int n)
{
    int max = arr[0];
    for (int i = 1; i < n; i++)
        if (arr[i] > max)
            max = arr[i];
    return max;
}

static void radixCountSort(std::vector<int>& arr, int n, int exp)
{
    std::vector<int> output(n);
    int i, count[10] = { 0 };

    for (i = 0; i < n; i++)
        count[(arr[i] / exp) % 10]++;

    for (i = 1; i < 10; i++)
        count[i] += count[i - 1];

    for (i = n - 1; i >= 0; i--) {
        output[count[(arr[i] / exp) % 10] - 1] = arr[i];
        count[(arr[i] / exp) % 10]--;
    }

    for (i = 0; i < n; i++)
        arr[i] = output[i];
}

void radixSortMain(std::vector<int>& arr, int n)
{
    int max = radixGetMax(arr, n);

    for (int exp = 1; max / exp > 0; exp *= 10)
        radixCountSort(arr, n, exp);
}

/////////////////////////////////////////////////////////////////////////////////

void selectionSortMain(std::vector<int>& arr, int n)
{
    for (int i = 0; i < n-1; ++i) {
        int min_value_index = i;

        for (int j = i+1; j < n; j++) {
            if (arr[j] < arr[min_value_index]) {
                min_value_index = j;
            }
        }

        std::swap(arr[i], arr[min_value_index]);
    }
}

////////////////////////////////////////////////////////////////////////////

void shellSortMain(std::vector<int>& arr, int n) {
    int gap = 1;

    while (gap * 3 + 1 < n)
        gap = gap * 3 + 1;

    while (gap >= 1)
    {
        for (int i = 0; i < gap; i++)
        {
            for (int j = i; j < n - gap; j += gap)
            {
                int k = j;
                while (k >= i && arr[k] > arr[k + gap])
                {
                    std::swap(arr[k], arr[k + gap]);
                    k -= gap;
                }
            }
        }

        gap = gap / 3;
    }
}

///////////////////////////////////////////////////////////////////////////

static void timSortInsertionSort(std::vector<int>& arr, int left, int right)
{
    for (int i = left + 1; i <= right; i++) {
        int temp = arr[i];
        int j = i - 1;
        while (j >= left && arr[j] > temp) {
            arr[j + 1] = arr[j];
            j--;
        }
        arr[j + 1] = temp;
    }
}

static void timSortMerge(std::vector<int>& arr, int l, int m, int r)
{
    int len1 = m - l + 1, len2 = r - m;

    std::vector<int> left(len1);
    std::vector<int> right(len2);

    for (int i = 0; i < len1; i++)
        left[i] = arr[l + i];
    for (int i = 0; i < len2; i++)
        right[i] = arr[m + 1 + i];

    int i = 0;
    int j = 0;
    int k = l;

    while (i < len1 && j < len2) {
        if (left[i] <= right[j]) {
            arr[k] = left[i];
            i++;
        }
        else {
            arr[k] = right[j];
            j++;
        }
        k++;
    }

    while (i < len1) {
        arr[k] = left[i];
        k++;
        i++;
    }

    while (j < len2) {
        arr[k] = right[j];
        k++;
        j++;
    }
}

void timSortMain(std::vector<int>& arr, int n)
{
    for (int i = 0; i < n; i += RUN)
        timSortInsertionSort(arr, i, std::min((i + RUN - 1), (n - 1)));

    for (int size = RUN; size < n; size = 2 * size) {
        for (int left = 0; left < n; left += 2 * size) {

            int mid = left + size - 1;
            int right = std::min((left + 2 * size - 1), (n - 1));

            if (mid < right)
                timSortMerge(arr, left, mid, right);
        }
    }
}

void shakerSortMain(std::vector<int> &arr, int n)
{
    bool swapped = true;
    int start = 0;
    int end = n - 1;

    while (swapped)
    {
        swapped = false;

        for (int i = start; i < end; ++i)
        {
            if (arr[i] > arr[i + 1])
            {
                std::swap(arr[i], arr[i + 1]);
                swapped = true;
            }
        }

        if (!swapped)
            break;

        swapped = false;

        --end;

        for (int i = end - 1; i >= start; --i)
        {
            if (arr[i] > arr[i + 1])
            {
                std::swap(arr[i], arr[i + 1]);
                swapped = true;
            }
        }
        ++start;
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////

void parallelMergeSort(std::vector<int> &arr, int left, int right)
{
    if (left >= right)
        return;

    const int PARALLEL_CUTOFF = 2000;
    int size = right - left + 1;

    if (size <= PARALLEL_CUTOFF)
    {
        mergeSort(arr, left, right);
    }
    else
    {
        int mid = left + (right - left) / 2;
        
        #pragma omp task shared(arr)
        {
            parallelMergeSort(arr, left, mid);
        }

        #pragma omp task shared(arr)
        {
            parallelMergeSort(arr, mid + 1, right);
        }
.
        #pragma omp taskwait

        mergeSortMerge(arr, left, mid, right);
    }
}

void parallelMergeSortMain(std::vector<int> &arr, int n, int max_threads)
{
#pragma omp parallel num_threads(max_threads)
    {
#pragma omp single 
        parallelMergeSort(arr, 0, n - 1);
    }
}

void ParallelOddEvenMain(std::vector<int> &arr, int n, int max_threads)
{
    int phase;
    int tmp, i;

#pragma omp parallel num_threads(max_threads) default(none) shared(arr, n) private(i, tmp, phase)
    for (phase = 0; phase < n; phase++)
    {
        if (phase % 2 == 0)
        {
#pragma omp for // omp for tem uma barreira implícita no fim do loop
            for (i = 1; i < n; i += 2)
            {
                if (arr[i - 1] > arr[i])
                {
                    std::swap(arr[i - 1], arr[i]);
                }
            }
        }
        else
        {
#pragma omp for
            for (i = 1; i < n - 1; i += 2)
            {
                if (arr[i] > arr[i + 1])
                {
                    std::swap(arr[i], arr[i + 1]);
                }
            }
        }
    }
}


static void quickSortSequential(std::vector<int> &array, int low, int high)
{
    if (low < high)
    {
        midOfThree(array, low, high);
        int part = partition(array, low, high);

        quickSortSequential(array, low, part);

        quickSortSequential(array, part + 1, high);
    }
}

void quickSortOmpTask(std::vector<int> &arr, int low, int high, int deep)
{
    if (low < high)
    {
        if (deep >= 0)
        {            
            midOfThree(arr, low, high);
            int part = partition(arr, low, high);

#pragma omp task shared(arr)
            quickSortOmpTask(arr, low, part, deep - 1);

#pragma omp task shared(arr)
            quickSortOmpTask(arr, part + 1, high, deep - 1);

        }
        else
        {
            midOfThree(arr, low, high);
            int part = partition(arr, low, high);
            quickSortSequential(arr, low, part);
            quickSortSequential(arr, part + 1, high);
        }
    }
}

void parallelQuickMain(std::vector<int> &arr, int n, int max_threads)
{
#pragma omp parallel num_threads(max_threads) shared(arr)
#pragma omp single nowait
    quickSortOmpTask(arr, 0, n - 1, 10);
}


void ParallelShellMain(std::vector<int> &arr, int n, int max_threads)
{
    int gap = 1;

    while (gap * 3 + 1 < n)
        gap = gap * 3 + 1;

    while (gap >= 1)
    {
#pragma omp parallel for num_threads(max_threads)
        for (int i = 0; i < gap; i++)
        {
            for (int j = i; j < n - gap; j += gap)
            {
                int k = j;
                while (k >= i && arr[k] > arr[k + gap])
                {
                    std::swap(arr[k], arr[k + gap]);
                    k -= gap;
                }
            }
        }

        gap = gap / 3;
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////


void atualiza_csv(std::vector<std::string> linhas){
    std::cout <<"Abrindo arquivo \n";
    std::ofstream ofs;
    ofs.open (CSV_FILE, std::ofstream::out | std::ofstream::app);

    for(int i = 0; i < linhas.size(); ++i){
        ofs << linhas[i] << std::endl;
    }

    ofs << std::endl;

    ofs.close();
    std::cout <<"Arquivo fechado \n";
}

void read_array_from_binary_file(const std::string& filename, std::vector<int>& target_array) {
    std::ifstream inputFile(filename, std::ios::binary | std::ios::in);

    if (!inputFile.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for reading in binary mode." << std::endl;
        return;
    }

    inputFile.seekg(0, std::ios::end); // Move to the end of the file
    long long file_size = inputFile.tellg(); // Get current position (file size in bytes)
    inputFile.seekg(0, std::ios::beg); // Move back to the beginning of the file

    if (file_size % sizeof(int) != 0) {
        std::cerr << "Error: Binary file size is not a multiple of sizeof(int). File might be corrupt." << std::endl;
        inputFile.close();
        return;
    }

    const size_t num_elements = file_size / sizeof(int);

    target_array.resize(num_elements);

    std::cout << "Reading " << num_elements << " integers from binary file " << filename << " into vector..." << std::endl;

    inputFile.read(reinterpret_cast<char*>(target_array.data()), num_elements * sizeof(int));

    inputFile.close();

    std::cout << "Successfully read " << target_array.size() << " integers from " << filename << std::endl;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////

void benchmark_algorithm(PCM * m,
    const std::string& algorithm_name,
    std::function<void(std::vector<int>&, int)> sort_function,
    const std::vector<int>& original_array
){

    bool is_pre_ordenado = std::is_sorted(original_array.begin(), original_array.end());
    bool is_pre_ordenado_dec = std::is_sorted(original_array.begin(), original_array.end(), std::greater<int>());
    std::vector<std::string> resultados(EXECS_PER_ALG);

    std::cout << "Executando " << algorithm_name << "\n";

    for (int j = 0; j < EXECS_PER_ALG; ++j) {
        std::vector<int> copied_array = original_array;

        m->getAllCounterStates(sstate1, sktstate1, cstates1);
        auto start = std::chrono::high_resolution_clock::now();

        sort_function(copied_array, copied_array.size());

        auto stop = std::chrono::high_resolution_clock::now();
        m->getAllCounterStates(sstate2, sktstate2, cstates2);

        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
        std::string duration_str = std::to_string(duration.count());
        std::string consumed_joules = std::to_string(getConsumedJoules(sstate1, sstate2));
        resultados[j] = algorithm_name + "," + duration_str + "," + consumed_joules + "," + std::to_string(copied_array.size());
        resultados[j] += "," + std::to_string(is_pre_ordenado) + "," + std::to_string(is_pre_ordenado_dec) + ",";

        /*for (uint32 i = 0; i < m->getNumCores(); ++i)
        {
            for (int s = 0; s <= PCM::MAX_C_STATE; ++s)
            {

                if (m->isCoreCStateResidencySupported(s)){
                    int coreCStateResidency = getCoreCStateResidency(s, cstates1[i], cstates2[i]) * 100;
                    resultados[j] += "," + std::to_string(coreCStateResidency);
                }
            }
        }*/

        if(j == 0){
            if(std::is_sorted(copied_array.begin(), copied_array.end()))
                std::cout << "Está ordenado! \n";
            else
                std::cout << "Não está ordenado \n";
        }
    }

    atualiza_csv(resultados);
}


void benchmark_parallel_algorithm(PCM * m,
    const std::string& algorithm_name,
    std::function<void(std::vector<int>&, int, int)> sort_function,
    const std::vector<int>& original_array,
    int max_threads
){

    bool is_pre_ordenado = std::is_sorted(original_array.begin(), original_array.end());
    bool is_pre_ordenado_dec = std::is_sorted(original_array.begin(), original_array.end(), std::greater<int>());
    std::vector<std::string> resultados(EXECS_PER_ALG);

    std::cout << "Executando " << algorithm_name << "\n";

    for (int j = 0; j < EXECS_PER_ALG; ++j) {
        std::vector<int> copied_array = original_array;

        m->getAllCounterStates(sstate1, sktstate1, cstates1);
        auto start = std::chrono::high_resolution_clock::now();

        sort_function(copied_array, copied_array.size(), max_threads);

        auto stop = std::chrono::high_resolution_clock::now();
        m->getAllCounterStates(sstate2, sktstate2, cstates2);

        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
        std::string duration_str = std::to_string(duration.count());
        std::string consumed_joules = std::to_string(getConsumedJoules(sstate1, sstate2));
        resultados[j] = algorithm_name + "," + duration_str + "," + consumed_joules + "," + std::to_string(copied_array.size());
        resultados[j] += "," + std::to_string(is_pre_ordenado) + "," + std::to_string(is_pre_ordenado_dec);
        resultados[j] += "," + std::to_string(max_threads);

        if(j == 0){
            if(std::is_sorted(copied_array.begin(), copied_array.end()))
                std::cout << "Está ordenado! \n";
            else
                std::cout << "Não está ordenado \n";
        }
    }

    atualiza_csv(resultados);
}

void call_algorithms(PCM * m, const std::vector<int>& original_array, int max_threads){
    benchmark_algorithm(m, "Bubble Sort", bubbleSortMain, original_array);
    benchmark_algorithm(m, "Heap Sort", heapSortMain, original_array);
    benchmark_algorithm(m, "Insertion Sort", insertionSortMain, original_array);
    benchmark_algorithm(m, "Merge Sort", mergeSortMain, original_array);
    benchmark_algorithm(m, "Odd-Even Sort", oddEvenSortMain, original_array);
    benchmark_algorithm(m, "Radix Sort", radixSortMain, original_array);
    benchmark_algorithm(m, "Selection Sort", selectionSortMain, original_array);
    benchmark_algorithm(m, "Tim Sort", timSortMain, original_array);
    benchmark_algorithm(m, "Iterative Quick Sort", iterativeQuickSortMain, original_array);
    benchmark_algorithm(m, "Iterative Merge Sort", iterativeMergeSortMain, original_array);
    benchmark_algorithm(m, "Shell Sort", shellSortMain, original_array);
    benchmark_algorithm(m, "Quick Sort", quickSortMain, original_array);
    benchmark_algorithm(m, "Shaker Sort", shakerSortMain, original_array);
    

    benchmark_parallel_algorithm(m, "Parallel Odd-Even Sort", ParallelOddEvenMain, original_array, max_threads);
    benchmark_parallel_algorithm(m, "Parallel Merge Sort", parallelMergeSortMain, original_array, max_threads);    
    benchmark_parallel_algorithm(m, "Parallel Shell Sort", ParallelShellMain, original_array, max_threads);
    benchmark_parallel_algorithm(m, "Parallel Quick Sort", parallelQuickMain, original_array, max_threads);
}

int main() {
    PCM *m; // = PCM::getInstance();
    m = PCM::getInstance();
    PCM::ErrorCode returnResult = m->program();
    if (returnResult != PCM::Success) {
        std::cerr << "PCM couldn't start" << std::endl;
        std::cerr << "Error code: " << returnResult << std::endl;
        exit(1);
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////

    std::vector<int> original_array;
    read_array_from_binary_file(ARRAY_FILE, original_array);

    std::cout << "Array gerado \n";

    call_algorithms(m, original_array, 2);
    call_algorithms(m, original_array, 4);
    call_algorithms(m, original_array, 8);
    call_algorithms(m, original_array, 12);
    call_algorithms(m, original_array, 16);
    call_algorithms(m, original_array, 20);

    std::sort(original_array.begin(), original_array.end());
    std::cout << "Array ordenado pelo std \n";

    call_algorithms(m, original_array, 2);
    call_algorithms(m, original_array, 4);
    call_algorithms(m, original_array, 8);
    call_algorithms(m, original_array, 12);
    call_algorithms(m, original_array, 16);
    call_algorithms(m, original_array, 20);

    std::sort(original_array.begin(), original_array.end(), std::greater<int>());
    std::cout << "Array ordenado em ordem decrescente pelo std \n";

    call_algorithms(m, original_array, 2);
    call_algorithms(m, original_array, 4);
    call_algorithms(m, original_array, 8);
    call_algorithms(m, original_array, 12);
    call_algorithms(m, original_array, 16);
    call_algorithms(m, original_array, 20);

    /////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////

    m->cleanup();

    return 0;
}
