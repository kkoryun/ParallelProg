#include <omp.h>
#include <iostream>
#include <math.h>
#include <cstdlib>
#include <ctime>

#include <vector>
#include <string>
#include <stdexcept>

#define r 100

using namespace std;

template <class T>
class Row
{
public:
	Row(T *data, size_t size) : data(data), size(size){};
	T &operator[](size_t index)
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
	}
	Row<T> operator[](size_t index)
	{
		if (index >= rows)
		{
			throw out_of_range("Mat does not have row with index: " + to_string(index));
		}
		return mat[index];
	}
    
    //void copy() {

    //}
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
void LU_Simple_Decomposition(Matrix<T> A, Matrix<T> L, Matrix<T> U)
{
	size_t N = L.rows;
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			U[i][j] = A[i][j];

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
	for (int i = 1; i < N; i++)
		for (int j = 0; j < i; j++)
			U[i][j] = 0;
}
template <class T>
void Left_Triangle_Many_Systems_Solve(Matrix<T> L, Matrix<T> X, Matrix<T> B, int N, int K)
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
void Right_Triangle_Many_Systems_Solve(Matrix<T> U, Matrix<T> X, Matrix<T> B, int N, int K)
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
void LU_Decomposition_(Matrix<T> A, Matrix<T> L, Matrix<T> U, int N)
{
	if (N <= r)
	{
		LU_Simple_Decomposition(A, L, U);
		return;
	}

	Matrix<double> A11(new double[r * r], r, r);
	Matrix<double> A12(new double[r * (N - r)], r, N - r);
	Matrix<double> A21(new double[(N - r) * r], N - r, r);
	Matrix<double> A22(new double[(N - r) * (N - r)], N - r, N - r);

	Matrix<double> L11(new double[r * r], r, r);
	Matrix<double> L21(new double[(N - r) * r], N - r, r);
	Matrix<double> L22(new double[(N - r) * (N - r)], N - r, N - r);

	Matrix<double> U11(new double[r * r], r, r);
	Matrix<double> U12(new double[r * (N - r)], r, N - r);
	Matrix<double> U22(new double[(N - r) * (N - r)], N - r, N - r);

#pragma omp parallel for
	for (int i = 0; i < r; i++)
		for (int j = 0; j < r; j++)
			A11[i][j] = A[i][j];

#pragma omp parallel for
	for (int i = 0; i < r; i++)
		for (int j = 0; j < N - r; j++)
			A12[i][j] = A[i][j + r];

#pragma omp parallel for
	for (int i = 0; i < N - r; i++)
		for (int j = 0; j < r; j++)
			A21[i][j] = A[(i + r)][j];

#pragma omp parallel for
	for (int i = 0; i < N - r; i++)
		for (int j = 0; j < N - r; j++)
			A22[i][j] = A[(i + r)][j + r];

	LU_Simple_Decomposition(A11, L11, U11);
	Left_Triangle_Many_Systems_Solve(L11, U12, A12, r, N - r);
	Right_Triangle_Many_Systems_Solve(U11, L21, A21, r, N - r);

#pragma omp parallel for
	for (int i = 0; i < N - r; i++)
		for (int j = 0; j < N - r; j++)
			for (int k = 0; k < r; k++)
				A22[i][j] -= L21[i][k] * U12[k][j];
	LU_Decomposition_(A22, L22, U22, N - r);
#pragma omp parallel for
	for (int i = 0; i < r; i++)
		for (int j = 0; j < r; j++)
		{
			L[i][j] = L11[i][j];
			U[i][j] = U11[i][j];
		}
#pragma omp parallel for
	for (int i = 0; i < r; i++)
		for (int j = 0; j < N - r; j++)
		{
			L[i][j + r] = 0;
			U[i][j + r] = U12[i][j];
		}
#pragma omp parallel for
	for (int i = 0; i < N - r; i++)
		for (int j = 0; j < r; j++)
		{
			L[(i + r)][j] = L21[i][j];
			U[(i + r)][j] = 0;
		}
#pragma omp parallel for
	for (int i = 0; i < N - r; i++)
		for (int j = 0; j < N - r; j++)
		{
			L[(i + r)][j + r] = L22[i][j];
			U[(i + r)][j + r] = U22[i][j];
		}
    A11.release();
    A12.release();
    A21.release();
	A22.release();
	L21.release();
	L11.release();
	L22.release();
	U12.release();
	U11.release();
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
    int N = 10;
    double* A = new double[N*N];
    srand(clock());
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            //A[i*N + j] = 1 + (1.1 + i * j / (i + j + 1) - i + 2 * j) / (1.1 + N*N / (2*N + 1) + N);
            A[i*N + j] = (double)rand() / RAND_MAX;
        }
    }
    double* L = new double[N * N];
    double* U = new double[N * N];
    omp_set_num_threads(1);
    int t1 = clock();
    LU_Decomposition(A, L, U, N);
    int t2 = clock();
    cout << t2 - t1 << " ms \n";
    int diff1 = abs(t2 - t1);
    omp_set_num_threads(4);
    t1 = clock();
    //LU_Simple_Decomposition(A, L, U, N);
    t2 = clock();
    cout << t2 - t1 << " ms \n";
    int diff2 = abs(t2 - t1);
 
    double max_err = 0;
    double max_A = 0;
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            double lu = 0;
            for (int k = 0; k < N; k++)
            {
                lu += L[N * i + k] * U[N * k + j];
            }
            double err = fabs(lu - A[N * i + j]);
            if (err > max_err) max_err = err;
            if (A[N * i + j] > max_A) max_A = A[N * i + j];
        }
    }
    cout << max_err / max_A;
    //cin >> N;
    return 0;
}