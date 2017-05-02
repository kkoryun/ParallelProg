#pragma once
#include<vector>
#include<omp.h>
#include<math.h>

template<class T> class BatcherSort
{
public:

	BatcherSort() {
		data = 0;
		size = 0;
	}

	BatcherSort(T * _data, int _size) {
		Set(_data, _size);
	}

	BatcherSort(std::vector<T> & _data) {
		Set(_data);
	}

	void Sort() {
		if (log2(size) > (int)log2(size))
			t = (int)log2(size) + 1;
		else
			t = (int)log2(size);
		h = 1 << (t - 1);
		while (h > 0)
		{
			int tmp;
			r = 0;
			z = h;
			q = 1 << (t - 1);
			while (q >= h)
			{

				for (int i = 0; (i < size - z); i++)
				{
					comparers(i, i + z);
				}
				z = q - h;
				q = q >> 1;
				r = h;

			}
			h = h >> 1;
		}
	}

	void Sort(int thread_count) {

		if (log2(size) > (int)log2(size))
			t = (int)log2(size) + 1;
		else
			t = (int)log2(size);
		h = 1 << (t - 1);

		while (h > 0)
		{
			int tmp;
			r = 0;
			z = h;
			q = 1 << (t - 1);
			while (q >= h)
			{

omp_set_num_threads(thread_count);
#pragma omp parallel for
				for (int i = 0; (i < size - z); i++)
				{
					comparers(i, i + z);
				}
				z = q - h;
				q = q >> 1;
				r = h;

			}
			h = h >> 1;
		}
	}


	void Set(T * _data, int _size) {
		data = = _data;
		size = _size;
	}

	void Set(std::vector<T> & _data) {
		data = _data.data();
		size = _data.size();
	}

private:

	T * data;
	int size;
	int h;
	int q;
	int r;
	int z;
	int t;
	
	inline void comparers(int index1, int index2) 
	{
		if ((index1 & h) == r)
		{
			if (data[index1] > data[index2])
			{
				swap(index1,index2);
			}
		}
	}

	inline void swap(int index1, int index2) 
	{
		int tmp = data[index1];
		data[index1] = data[index2];
		data[index2] = tmp;
	}

};


