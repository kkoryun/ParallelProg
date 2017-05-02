#include<math.h>
#include<iostream>
//#include<iostream>
#include<random>
#include<vector>
#include<algorithm>
#include<time.h>
#include<omp.h>

#include"BatcherSort.h"

void swap(int * data, int index1, int index2) {
	int tmp = data[index1];
	data[index1] = data[index2];
	data[index2] = tmp;
}

void fill_random(std::vector<int>& data) {
	//srand(1);
	for (size_t i = 0; i < data.size(); i++)
	{
		data[i] = rand() % 100;
	}
}

int main() {
	double t1, t2, t3, t4;
	size_t data_size = 8876665;
	size_t time = 0;
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

			t3 = clock();

//omp_set_num_threads(4);
#pragma omp parallel for private(tmp) 
				for (int i = 0; (i < data_size - dist); i++)
				{
					if ((i & h) == r) {
						if (data[i] > data[i + dist]) {
							 tmp = data[i];
							data[i] = data[i + dist];
							data[i + dist] = tmp;
						}
					}
				}
			

			t4 = clock();

			time += t4 - t3;
			//std::cout << "sort step " << t4 - t3 << std::endl;
			dist = q - h;
			q = q >> 1;
			r = h;
			count++;
		} 
		h = h >> 1;
		
	}

	t2 = clock();
	std::cout << "betcher sort "<<t2-t1<< std::endl;

	
	//BatcherSort<int> bsort(vec);

	/*t1 = clock();
	bsort.Sort();
	t2 = clock();
	std::cout << "batcher sort " << t2 - t1 << std::endl;*/
	
	/*t1 = clock();
	bsort.Sort(4);
	t2 = clock();
	std::cout << "batcher sort " << t2 - t1 << std::endl;*/
	
	
	char c;
	if (vec_copy != vec) {
		std::cout << "no sorted" << std::endl;
		std::cin >> c;
		return -1;
	}

	std::cin >> c;
	return 0;
}