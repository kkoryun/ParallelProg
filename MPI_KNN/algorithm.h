#pragma once
#include"struct.h"
#include <math.h>
template<class T>void Swap(T * left, T * right) {
	T tmp = *left;
	*left = *right;
	*right = tmp;
}
//METRIC
template<class T> T metric(point<T> & x, point<T> & y, int size) {
	T res = 0;
	for (size_t i = 0; i < size; i++)
		res+=(x.data[i] - y.data[i])*(x.data[i] - y.data[i]);

	res = sqrt(res);
	return res;
}
//findKMin
template<class T>  void findKMin(T * vec, int size, int K,int * minIndex) {


	T min;
	T * tmp = new T[size];
	int * ind = new int[size];
	for (size_t i = 0; i < size; i++)
		ind[i] = i;
	memcpy(tmp, vec, size * sizeof(T));
	

	int minInd = -1;
	for (size_t j = 0; j < K; j++)
	{
		min = tmp[j];
		minInd = j;
		for (size_t i = j; i < size; i++) 
		{
			if (tmp[i] < min) 
			{
				min = tmp[i];
				minInd = i;
			}
		}
		Swap<int>(&ind[j], &ind[minInd]);
		Swap<T>(&tmp[j], &tmp[minInd]);
	}
	memcpy(minIndex, ind, K*sizeof(int));
		delete[] tmp;
	delete[] ind;
}
//caclWeigth
template<class T> void Weight(T * metrics, int size, T * weight ){

}


