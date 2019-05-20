#include <omp.h>
#include <iostream>
#include <math.h>
#include <cstdlib>
#include <ctime>
#include <vector>

using namespace std;

// struct CRSMatrix
// {
//     int n;
//     int m;
//     int nz;
//     vector<double> val;
//     vector<int> colIndex;
//     vector<int> rowPtr;
// };

using vec = vector<double>;
vec operator*(double a, const vec& v) {
    vec res(v.size());
    #pragma omp parallel for
    for (int i = 0; i < v.size(); i++)
    {
        res[i] = a * v[i];
    }
    return res;
}
vec operator*(const vec& v, double a) {
    return operator*(a, v);
}
double operator*(const vec& a, const vec& b) {
    int size = a.size();
    double sum = 0;
    #pragma omp parallel for shared(sum, a) reduction(+: sum)
    for (int i = 0; i < size; i++)
    {
        sum += a[i] * b[i];
    }
    return sum;
}
vec operator-(const vec& a, const vec& b) {
    int size = a.size();
    vec c(size);
    #pragma omp parallel for
    for (int i = 0; i < size; i++)
    {
        c[i] = a[i] - b[i];
    }
    return c;
}
vec operator+(const vec& a, const vec& b) {
    int size = a.size();
    vec c(size);
    #pragma omp parallel for
    for (int i = 0; i < size; i++)
    {
        c[i] = a[i] + b[i];
    }
    return c;
}

vec operator*(CRSMatrix &A, vec &x)
{
    int size = A.n;
    vec res(A.n);

    #pragma omp parallel for
    for (int i = 0; i < size; i++)
    {
        res[i] = 0;
        if (i < A.n - 1)
        {
            for (int j = A.rowPtr[i]; (j < A.nz) && (j < A.rowPtr[i + 1]); j++)
                res[i] += A.val[j] * x[A.colIndex[j]];
        }
        else
        {
            for (int j = A.rowPtr[i]; j < A.nz; j++)
                res[i] += A.val[j] * x[A.colIndex[j]];
        }
    }

    #pragma omp parallel for
    for (int i = 0; i < size; i++)
    {
        if (i < A.n - 1)
        {
            for (int j = A.rowPtr[i]; (j < A.nz) && (j < A.rowPtr[i + 1]); j++)
                if (A.colIndex[j] != i)
                    res[A.colIndex[j]] += A.val[j] * x[i];
        }
        else
        {
            for (int j = A.rowPtr[i]; j < A.nz; j++)
                if (A.colIndex[j] != i)
                    res[A.colIndex[j]] += A.val[j] * x[i];
        }
    }
    return res;
}

void SLE_Solver_CRS(CRSMatrix &A, vec &b, double eps, int max_iter, vec &x, int &count)
{
    int size = A.n;
    x = b;

    vec r(A.n);
    vec p(A.n);
    vec Ax = A * x;


    r = b - Ax;
    p = r;

    double b_norm = sqrt(b * b);
    
    for (count = 1; count <= max_iter; count++)
    {
        vec Ap = A * p;

        double App = Ap * p;
        if (App == 0)
            break;
        double alpha = r * r / App;
        //  double x_norm = 0;
        // TODO

         vec tmp = alpha * p;
        x = x +  tmp;
        double x_norm = tmp * tmp;

        x_norm = sqrt(x_norm);
        if (x_norm / b_norm < eps)
        {
            return;
        }
        vec new_r(size);
        #pragma omp parallel for
        for (int i = 0; i < A.n; i++)
        {
            new_r[i] = r[i] - alpha * Ap[i];
        }
        double rr = r * r;
        if (rr == 0)
            break;
        double beta = (new_r * new_r) / rr;

        r = new_r;
        p = r + beta * p;
     
    }
}

void SLE_Solver_CRS(CRSMatrix &A, double *b, double eps, int max_iter, double *x, int &count)
{
    int size = A.n;
    vec b_ = vec(b, b + size);
    vec x_= vec(x, x + size);
    SLE_Solver_CRS(A, b_, eps, max_iter, x_ , count);
    for (int j = 0; j < size; j++)
    x[j] = x_[j];
}

int main(){
    return 0;
}
