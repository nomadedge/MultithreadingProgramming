#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <omp.h>
#include <algorithm>
#include <stdlib.h>
#include <chrono>
#include <ctime>
#include <cmath>

using namespace std;

class Timer
{
public:
    void start()
    {
        m_StartTime = std::chrono::system_clock::now();
        m_bRunning = true;
    }

    void stop()
    {
        m_EndTime = std::chrono::system_clock::now();
        m_bRunning = false;
    }

    double elapsedMilliseconds()
    {
        std::chrono::time_point<std::chrono::system_clock> endTime;

        if (m_bRunning)
        {
            endTime = std::chrono::system_clock::now();
        }
        else
        {
            endTime = m_EndTime;
        }

        return std::chrono::duration_cast<std::chrono::milliseconds>(endTime - m_StartTime).count();
    }

    double elapsedSeconds()
    {
        return elapsedMilliseconds() / 1000.0;
    }

private:
    std::chrono::time_point<std::chrono::system_clock> m_StartTime;
    std::chrono::time_point<std::chrono::system_clock> m_EndTime;
    bool                                               m_bRunning = false;
};

bool canParseInt(string line)
{
    char* p;
    strtol(line.c_str(), &p, 10);
    return *p == 0;
}

bool canParseDouble(string line)
{
    char* p;
    strtod(line.c_str(), &p);
    return *p == 0;
}

int factorial(int number)
{
    int res = 1;
    for (int i = 2; i <= number; i++)
    {
        res *= i;
    }
    return res;
}

int merge(int* arr, int* temp, int left, int mid, int right)
{
    int i, j, k;
    int inv_count = 0;

    i = left;
    j = mid;
    k = left;
    while ((i <= mid - 1) && (j <= right)) {
        if (arr[i] <= arr[j]) {
            temp[k++] = arr[i++];
        }
        else {
            temp[k++] = arr[j++];
            inv_count = inv_count + (mid - i);
        }
    }

    while (i <= mid - 1)
        temp[k++] = arr[i++];

    while (j <= right)
        temp[k++] = arr[j++];

    for (i = left; i <= right; i++)
        arr[i] = temp[i];

    return inv_count;
}

int mergeSort(int* arr, int* temp, int left, int right)
{
    int mid, inv_count = 0;
    if (right > left) {
        mid = (right + left) / 2;

        inv_count += mergeSort(arr, temp, left, mid);
        inv_count += mergeSort(arr, temp, mid + 1, right);

        inv_count += merge(arr, temp, left, mid + 1, right);
    }
    return inv_count;
}

int countInversions(int* arr, int array_size)
{
    int* buf = new int[array_size];
    for (int i = 0; i < array_size; i++)
    {
        buf[i] = arr[i];
    }
    int* temp = new int[array_size];
    int result = mergeSort(buf, temp, 0, array_size - 1);
    delete[] temp;
    delete[] buf;
    return result;
}

int main(int argc, char** argv)
{
    if (argc != 3)
    {
        cout << "Pass the file path and threads count to the program.";
        return 1;
    }
    int threadsCount = atoi(argv[2]);
    if (argv[2] == "" || !canParseInt(argv[2]) || threadsCount < -1)
    {
        cout << "Threads count must be equal to or more than -1.";
        return 1;
    }
    ifstream matrixFile;
    matrixFile.open(argv[1]);
    if (!matrixFile.is_open())
    {
        cout << "File either does not exist or you do not have enough permissions.";
        return 1;
    }
    vector<string> lines = {};
    string s;
    while (getline(matrixFile, s))
        lines.push_back(s);
    matrixFile.close();

    int matrixSize;
    double** matrix;

    try
    {
        matrixSize = atoi(lines[0].c_str());
        matrix = new double* [matrixSize];
        for (int i = 0; i < matrixSize; i++)
        {
            matrix[i] = new double[matrixSize];
        }

        for (int i = 0; i < matrixSize; i++)
        {
            istringstream iss(lines[i + 1]);
            for (int j = 0; j < matrixSize; j++)
            {
                string subs;
                iss >> subs;
                if (subs == "" || !canParseDouble(subs))
                {
                    throw exception();
                }
                matrix[i][j] = atof(subs.c_str());
            }
        }
    }
    catch (const exception&)
    {
        cout << "File content does not correspond required format.";
        return 1;
    }

    int factor = factorial(matrixSize);
    int** permutatuions = new int*[factor];
    for (int i = 0; i < factor; i++)
    {
        permutatuions[i] = new int[matrixSize];
    }
    int* original = new int[matrixSize];
    for (int i = 0; i < matrixSize; i++)
    {
        original[i] = i;
    }
    int k = 0;
    do {
        for (int i = 0; i < matrixSize; i++)
        {
            permutatuions[k][i] = original[i];
        }
        k++;
    } while (next_permutation(original, original + matrixSize));

    /*string scheduleTypes[3] = {"static", "dynamic", "guided"};
    int chunksCount = matrixSize / threadsCount;
    int* chunkSizes = new int[chunksCount];
    for (int i = 0; i < chunksCount; i++)
    {
        chunkSizes[i] = i + 1;
    }
    int paramsCount = chunksCount * 3;
    string* scheduleParams = new string[paramsCount];
    int index = 0;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < chunksCount; j++)
        {
            scheduleParams[index] = scheduleTypes[i] + "," + to_string(chunkSizes[j]);
            index++;
        }
    }*/
    double sum;
    omp_set_dynamic(0);
    double averageTime = 0;
    for (int i = 0; i < 5; i++)
    {
        //_putenv_s("OMP_SCHEDULE", scheduleParams[i].c_str());
        Timer timer;
        timer.start();
        sum = 0;
        if (threadsCount == -1)
        {
            for (int i = 0; i < factor; i++)
            {
                double number = pow(-1, countInversions(permutatuions[i], matrixSize));
                for (int j = 0; j < matrixSize; j++)
                {
                    number *= matrix[j][permutatuions[i][j]];
                }
                sum += number;
            }
        }
        else if (threadsCount == 0)
        {
            #pragma omp parallel
            {
                double localSum = 0;
                #pragma omp for
                for (int i = 0; i < factor; i++)
                {
                    double number = pow(-1, countInversions(permutatuions[i], matrixSize));
                    for (int j = 0; j < matrixSize; j++)
                    {
                        number *= matrix[j][permutatuions[i][j]];
                    }
                    localSum += number;
                }
                #pragma omp atomic
                sum += localSum;
            }
        }
        else
        {
            #pragma omp parallel num_threads(threadsCount)
            {
                double localSum = 0;
                #pragma omp for
                for (int i = 0; i < factor; i++)
                {
                    double number = pow(-1, countInversions(permutatuions[i], matrixSize));
                    for (int j = 0; j < matrixSize; j++)
                    {
                        number *= matrix[j][permutatuions[i][j]];
                    }
                    localSum += number;
                }
                #pragma omp atomic
                sum += localSum;
            }
        }
        
        timer.stop();
        averageTime += timer.elapsedMilliseconds();
    }
    averageTime /= 5;
    printf("\nTime (%i thread(s)): %f ms\n", threadsCount, averageTime);

    for (int i = 0; i < matrixSize; i++)
    {
        delete[] matrix[i];
    }
    delete[] matrix;
    for (int i = 0; i < factor; i++)
    {
        delete[] permutatuions[i];
    }
    delete[] permutatuions;

    printf("Determinant: %f\n", sum);

    return 0;
}