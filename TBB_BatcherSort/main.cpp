#include<vector>
#include<time.h>
#include<iostream>
#include<algorithm>
#include<math.h>

#include "tbb/task_scheduler_init.h"
#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"
using namespace std;

void fill_random(std::vector<int>& data) {
	for (size_t i = 0; i < data.size(); i++)
	{
		data[i] = rand() % 100;
	}
}

template<class T> class BatcherSorttbb
{

public:
	BatcherSorttbb(T* data, int h, int r, int z) :
		data(data), h(h), r(r), z(z) {

	}

	void operator()(const tbb::blocked_range<int> & range) const {
		for (int i = range.begin(); i != range.end(); i++)
		{
			if ((i & h) == r)
			{
				if (data[i] > data[i + z])
				{
					int tmp = data[i];
					data[i] = data[i + z];
					data[i + z] = tmp;
				}
			}
		}
	}
private:
	T * data;
	int h;
	int r;
	int z;


};

int main() {
	double t1, t2;
	size_t data_size = 1000000;
	std::vector<int> vec(data_size), vec_copy(data_size);

	//random fill time
	t1 = clock();
	fill_random(vec);
	vec_copy = vec;
	if (vec_copy != vec) return -1;
	t2 = clock();
	std::cout << "fill and copy " << t2 - t1 << std::endl;


	//stl qsort timetest
	t1 = clock();
	std::sort(vec_copy.begin(), vec_copy.end());
	t2 = clock();
	std::cout << "slt sort " << t2 - t1 << std::endl;

	tbb::task_scheduler_init init;

	int t;
	if (log2(data_size) > (int)log2(data_size))
		t = (int)log2(data_size) + 1;
	else
		t = (int)log2(data_size);

	int h = 1 << (t - 1);
	int q;
	int r;
	int dist;

	int count = 0;
	int * data = vec.data();

	t1 = clock();

	while (h > 0)
	{
		int tmp;
		r = 0;
		dist = h;
		q = 1 << (t - 1);
		while (q >= h)
		{
			tbb::parallel_for(tbb::blocked_range<int>(0, data_size - dist), BatcherSorttbb<int>(data, h, r, dist));
			dist = q - h;
			q = q >> 1;
			r = h;
			count++;
		}
		h = h >> 1;
	}

	t2 = clock();
	std::cout << "betcher sort " << t2 - t1 << " iter count " << count << std::endl;

	if (vec_copy != vec) {
		std::cout << "no sorted" << std::endl;
		return -1;
	}
	char c;
	cin >> c;
	return 0;
}
