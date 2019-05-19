#include <vector>
#include <stdexcept>
#include <string>
#include <omp.h>
#include <cstdlib>
#include <math.h>
#include <iostream>
#include <ctime>
#include <functional>

using namespace std;

const int block_size = 200;
template <class T>
class Row
{
public:
    Row(T *data, size_t size) : data(data), size(size) {};
    T &operator[](size_t index)
    {
        if (index >= size)
        {
            throw out_of_range("Mat does not have column with index: " + to_string(index));
        }
        return data[index];
    }
    const T &operator[](size_t index) const
    {
        if (index >= size)
        {
            throw out_of_range("Mat does not have column with index: " + to_string(index));
        }
        return data[index];
    }
    T *data;
    size_t size;
};

template <class T>
class Matrix
{
public:
    Matrix(T *data, size_t rows, size_t columns)
    {
        for (size_t i = 0; i < rows; i++)
        {
            Row<T> tmp(data + i * columns, columns);
            mat.push_back(tmp);
        }
        this->rows = rows;
        this->columns = columns;
        this->data = data;
    }
    Row<T>& operator[](size_t index)
    {
        if (index >= rows)
        {
            throw out_of_range("Mat does not have row with index: " + to_string(index));
        }
        return mat[index];
    }
    const Row<T>& operator[](size_t index) const
    {
        if (index >= rows)
        {
            throw out_of_range("Mat does not have row with index: " + to_string(index));
        }
        return mat[index];
    }

    void copy_from(const Matrix<T>& m, std::function<size_t(size_t)> r, std::function<size_t(size_t)>c) {
#pragma omp parallel for
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < columns; j++)
                this->operator[](i)[j] = m[r(i)][c(j)];
    }

    void release()
    {
        mat.clear();
        delete[] data;
        data = nullptr;
        rows = 0;
    }

    T *data;
    vector<Row<T> > mat;
    size_t rows;
    size_t columns;
};


//template<class T>
//class LU_decompositor{
//    static void decompose() {
//
//    }
//};

template <class T>
void LU_Simple_Decomposition(Matrix<T>& A, Matrix<T>& L, Matrix<T>& U)
{
    size_t N = L.rows;
#pragma omp parallel for 
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            U[i][j] = A[i][j];

#pragma omp parallel
    for (int i = 0; i < N; i++)
    {
        L[i][i] = 1;
#pragma omp parallel for 
        for (int k = i + 1; k < N; k++)
        {
            double mu = U[k][i] / U[i][i];
            for (int j = i; j < N; j++)
            {
                U[k][j] -= mu * U[i][j];
                L[k][i] = mu;
                L[i][k] = 0;
            }
        }
    }
#pragma omp parallel for 
    for (int i = 1; i < N; i++)
        for (int j = 0; j < i; j++)
            U[i][j] = 0;
}
template <class T>
void Right_Up_System(Matrix<T>& L, Matrix<T>& X, Matrix<T>& B, int N, int K)
{
#pragma omp parallel for
    for (int k = 0; k < K; k++)
        for (int i = 0; i < N; i++)
        {
            X[i][k] = B[i][k];
            for (int j = 0; j < i; j++)
                X[i][k] -= L[i][j] * X[j][k];
            X[i][k] /= L[i][i];
        }
}
template <class T>
void Left_Down_System(Matrix<T>& U, Matrix<T>& X, Matrix<T>& B, int N, int K)
{
#pragma omp parallel for
    for (int k = 0; k < K; k++)
        for (int i = 0; i < N; i++)
        {
            X[k][i] = B[k][i];
            for (int j = 0; j < i; j++)
                X[k][i] -= U[j][i] * X[k][j];
            X[k][i] /= U[i][i];
        }
}




template <class T>
void LU_Decomposition_(Matrix<T>& A, Matrix<T>& L, Matrix<T>& U, int N)
{
    if (N <= block_size)
    {
        LU_Simple_Decomposition(A, L, U);
        return;
    }

    Matrix<double> A11(new double[block_size * block_size], block_size, block_size);
    Matrix<double> A12(new double[block_size * (N - block_size)], block_size, N - block_size);
    Matrix<double> A21(new double[(N - block_size) * block_size], N - block_size, block_size);
    Matrix<double> A22(new double[(N - block_size) * (N - block_size)], N - block_size, N - block_size);

    Matrix<double> L11(new double[block_size * block_size], block_size, block_size);

    Matrix<double> U11(new double[block_size * block_size], block_size, block_size);

    A11.copy_from(A, [](size_t i) { return i; }, [](size_t j) { return j; });
    A12.copy_from(A, [](size_t i) { return i; }, [](size_t j) { return j + block_size; });
    A21.copy_from(A, [](size_t i) { return i + block_size; }, [](size_t j) { return j; });
    A22.copy_from(A, [](size_t i) { return i + block_size; }, [](size_t j) { return j + block_size; });

/*
#pragma omp parallel for
    for (int i = 0; i < block_size; i++)
        for (int j = 0; j < block_size; j++)
            A11[i][j] = A[i][j];

#pragma omp parallel for
    for (int i = 0; i < block_size; i++)
        for (int j = 0; j < N - block_size; j++)
            A12[i][j] = A[i][j + block_size];

#pragma omp parallel for
    for (int i = 0; i < N - block_size; i++)
        for (int j = 0; j < block_size; j++)
            A21[i][j] = A[i + block_size][j];

#pragma omp parallel for
    for (int i = 0; i < N - block_size; i++)
        for (int j = 0; j < N - block_size; j++)
            A22[i][j] = A[i + block_size][j + block_size];
*/
    LU_Simple_Decomposition(A11, L11, U11);
    A11.release();
    Matrix<double> U12(new double[block_size * (N - block_size)], block_size, N - block_size);
    Right_Up_System(L11, U12, A12, block_size, N - block_size);
    A12.release();

    Matrix<double> L21(new double[(N - block_size) * block_size], N - block_size, block_size);
    Left_Down_System(U11, L21, A21, block_size, N - block_size);
    A21.release();
    Matrix<double> L22(new double[(N - block_size) * (N - block_size)], N - block_size, N - block_size);

    Matrix<double> U22(new double[(N - block_size) * (N - block_size)], N - block_size, N - block_size);
#pragma omp parallel for
    for (int i = 0; i < N - block_size; i++)
        for (int j = 0; j < N - block_size; j++)
            for (int k = 0; k < block_size; k++)
                A22[i][j] -= L21[i][k] * U12[k][j];
    LU_Decomposition_(A22, L22, U22, N - block_size);
    A22.release();
#pragma omp parallel for 
    for (int i = 0; i < block_size; i++)
        for (int j = 0; j < block_size; j++)
        {
            L[i][j] = L11[i][j];
            U[i][j] = U11[i][j];
        }
    L11.release();
    U11.release();
#pragma omp parallel for 
    for (int i = 0; i < block_size; i++)
        for (int j = 0; j < N - block_size; j++)
        {
            L[i][j + block_size] = 0;
            U[i][j + block_size] = U12[i][j];
        }
    U12.release();
#pragma omp parallel for 
    for (int i = 0; i < N - block_size; i++)
        for (int j = 0; j < block_size; j++)
        {
            L[(i + block_size)][j] = L21[i][j];
            U[(i + block_size)][j] = 0;
        }
    L21.release();
#pragma omp parallel for 
    for (int i = 0; i < N - block_size; i++)
        for (int j = 0; j < N - block_size; j++)
        {
            L[(i + block_size)][j + block_size] = L22[i][j];
            U[(i + block_size)][j + block_size] = U22[i][j];
        }
    L22.release();
    U22.release();
}

void LU_Decomposition(double *A_, double *L_, double *U_, int N)
{
    Matrix<double> A(A_, N, N);
    Matrix<double> L(L_, N, N);
    Matrix<double> U(U_, N, N);
    LU_Decomposition_(A, L, U, N);
}

int main()
{
    

    return 0;
}