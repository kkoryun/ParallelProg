#pragma once
template<typename TYPE>
class mat
{
public:
	TYPE * data;
		   size_t size, remaind ;
    TYPE ** datapnt;

	mat(size_t w = 1, size_t h = 1) {
		size = w*h;
		data = new TYPE [size];

		for (size_t i = 0; i < size; i++)
		{
			srand(5);
			data[i] =rand()%10;
		}
		datapnt = new TYPE * [h];
		for (size_t i = 0; i < h; i++)
			datapnt[i] = data + w*i;

		remaind = 0;
	};
	void resize(size_t pn){
		if (data != 0 && datapnt!=0) {
			delete[] datapnt;
			pnt = new TYPE*[pn];
			size_t len = size / pn;
			remaind = size % pn;

			for (size_t i = 0; i < pn; i++)
				pnt[i] = m + len*i;
		}
	}
	size_t getRemaund() {
		return remaind;
	}
	void showmatr(int** m, size_t w, size_t h)
	{
		for (int i = 0; i < h; i++)
		{
			for (int j = 0; j < w; j++)
				cout << m[i][j] << " ";
			cout << endl;
		}
	}
	~mat() {
		delete[] data;
		delete[] pnt;
	}
};

