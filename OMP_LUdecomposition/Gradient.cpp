#include <omp.h>
#include <iostream>
#include <math.h>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <algorithm>

using namespace std;

using vec = vector<double>;

struct CRSMatrix
{
    int n;
    int m;
    int nz;
    vec val;
    vector<int> colIndex;
    vector<int> rowPtr;
    vec operator*(const vec& v) const {
        vec res;
        //omp
        for (size_t i = 0; i < n; i++)
        {
            int first = rowPtr[i];
            int last;
            if(i == n-1)
                last = n;
            else
                last = rowPtr[i + 1];
            double sum = 0;
            for (int j = first;j < last;j++)
            {
                sum += val[j] * v[colIndex[j]];
            }
            for (int j = 0; j < first; j++)
            {
                if (colIndex[j] == colIndex[first])
                {
                    int pos = binary_search(rowPtr.begin(), rowPtr.end(), j);
                    sum += val[j] * v[pos];
                }
            }
            res.push_back(sum);
        }
        return res;
    }
};

struct Decision {
    vec x;
    vec p;
    vec r;
};

vec operator*(double a, vec& v);
vec operator*(vec& v, double a);
vec operator-(vec& a, vec& b);
vec operator+(vec& a, vec& b);
double operator*(vec& a, vec& b);

void fill_first_decision(const CRSMatrix & A, const double* const b, Decision& d);

void solve(CRSMatrix & A, Decision& last_d, Decision& curr_d) {
    double r_dot = last_d.r * last_d.r;
    double p_dot = (A*last_d.p) * last_d.p;
    double alpha = r_dot / p_dot;

    curr_d.x = last_d.x + alpha * last_d.p;
    curr_d.r = last_d.r - alpha * (A*last_d.p);
    curr_d.p = last_d.r + ((curr_d.r * curr_d.r) / (last_d.r * last_d.r))*last_d.p;
}

void SLE_Solver_CRS(CRSMatrix & A, double * b, double eps, int max_iter, double * x, int & count) {
    count = 0;
    double err;
    Decision last_d;
    Decision curr_d;
    fill_first_decision(A, b, last_d);
    while (err > eps && count < max_iter)
    {
        count++;
        solve(A, last_d, curr_d);
        auto tmp = (last_d.x - curr_d.x);
        err = sqrt(tmp * tmp);
    }
    //copy to x
}

int main()
{
    int n = 5;
    CRSMatrix A;
    A.m = n;
    A.n = n;
    A.nz = 8;
    double val[] = { 1, 1, 1, 1, 1, 1, 1, 1 };
    vec val1(val, val + 8);
    int col[] = { 0, 2, 4, 2, 4, 3, 3, 4 };
    vector<int> col1(col, col + 8);
    int row[] = { 0, 3, 5, 6, 7 };
    vector<int> row1(row, row + 5);
    A.val = val1;
    A.colIndex = col1;
    A.rowPtr = row1;
    vec x{ 1, 1, 1, 1, 1 };
    vec y = A * x;
    for (int i = 0; i < n; i++)
        cout << y[i] << endl;
    cout << endl;
    double *X = new double[n];
    int count;
    SLE_Solver_CRS(A, y.data(), 0.0001, 10, X, count);
    for (int i = 0; i < n; i++)
        cout << X[i] << endl;
    cout << "Count: " << count << endl;
    cin >> n;
    return 0;
}


vec operator*(double a, vec& v) {
    vec res(v.size());
    //#omp 
    for (size_t i = 0; i < v.size(); i++)
    {
        res[i] = a * v[i];
    }
    return res;
}
vec operator*(vec& v, double a) {
    return operator*(a, v);
}
double operator*(vec& a, vec& b) {
    int size = a.size();
    double sum = 0;
    //#omp share
    for (size_t i = 0; i < size; i++)
    {
        sum += a[i] * b[i];
    }
    return sum;
}
vec operator-(vec& a, vec& b) {
    int size = a.size();
    vec c(size);
    //#omp
    for (size_t i = 0; i < size; i++)
    {
        c[i] = a[i] - b[i];
    }
    return c;
}
vec operator+(vec& a, vec& b) {
    int size = a.size();
    vec c(size);
    //#omp
    for (size_t i = 0; i < size; i++)
    {
        c[i] = a[i] + b[i];
    }
    return c;
}

void fill_first_decision(const CRSMatrix & A, const double* const b, Decision& d) {
    d.x.resize(A.n);
    for (size_t i = 0; i < A.n; i++)
        d.x[i] = b[i];
    d.r = vec(b, b + A.n) - A * d.x;
    d.p = d.r;
}