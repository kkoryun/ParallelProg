#pragma once
#include<string>
#include<vector>
#include<assert.h>

using namespace std;
template<typename T>struct TPoint {
	T * data;
    int size;
	int Class;
};
template<class T> class point
{
public:
	T * data;
	int size;
	int Class;
	point(int _size = 0) {
		if (_size != 0) {
			size = _size;
			data = new T[size];
			Class = -1;
		}
		else {
			data = 0;
			Class = -1;
		}
	}
	point(T * datapnt, int pointsize,int _Class) {
		size = pointsize;
		if (size != 0) {
			data = new T[size];
			for (size_t i = 0; i < size; i++)
				data[i] = datapnt[i];
			Class = _Class;
		}
	}
	~point() {
		if (data != 0)
			delete[] data;
	}
	point<T> & operator=(point<T> & p) {
		Class = p.Class;
		size = p.size;
		if (data != 0)
			delete[] data;
		data = new T[size];
		memcpy(data, p.data, size * sizeof(T));
			
		return *this;
	}
	friend ostream & operator<<(ostream & os, point<T> & p) {
		for (size_t i = 0; i < p.size; i++)
			os << p.data[i] << " " ;
		os << "class: " << p.Class << endl;
		return os;
	}
	
};
/*template<class T> class set   {
public:
	set(int _size = 0, int _pointsize=0) {
		size = _size;
		pointsize = _pointsize;
		if (size != 0) {
			data = new point<T>[size];
			for (size_t i = 0; i < size; i++)
				data[i] = point<T>(pointsize);
		}
		else
		{
			data = 0;
		}
	}
	set(int _size, int _pointsize, T pointdata) {
		size = _size;
		pointsize = _pointsize;
		data = new point<T>[size];
		
		for (size_t i = 0; i < size; i++) 
			data[i] = point<T>(pointdata + (pointSize + 1)*i, pointSize);
		
	}
	~set() {
		if (data != 0)
			delete[] data;
	}
	

	point<T> * data;
	int size, pointsize;
};*/
typedef point<double> pointd;

