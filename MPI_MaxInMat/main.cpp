#include<mpi.h>
#include<iostream>
#include<conio.h>
#include"mat.h"
using namespace std;
void showmatr(int** m, size_t w, size_t h)
{
	for (int i = 0; i < h; i++)
	{
		for (int j = 0; j < w; j++)
			cout << m[i][j] << " ";
		cout << endl;
	}
}

int FindMax(int * v, size_t count) {
	int max = v[0];
	for (size_t i = 0; i < count; i++)
		if (v[i] >= max) max = v[i];

	return max;
}

int main(int argc, char **argv) {
	int rank, size;
	size_t h = 5000, w = 5000;
	double time1, time2;
	int MAX;
	size_t remaind;
	int * sendbuf;
	int * sendcount;
	int * displs;
	int * recvbuf;
	mat<int> *m;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	time1 = MPI_Wtime();

	if (rank == 0) {
		m = new mat<int>(w, h);
		sendbuf = m->data;
	}
		sendcount = new int[size];
		displs = new int[size];
		int siz = w*h / size, rem = w*h %size;

		for (int i = 0; i < size - 1; i++)
		{
			sendcount[size - i - 1] = siz;
			displs[size - i - 1] = i * siz;
		}
		sendcount[0] = siz + rem;
		displs[0] = (size - 1) * siz;


	recvbuf = new int[sendcount[rank]];
	time2 = MPI_Wtime();
	MPI_Barrier(MPI_COMM_WORLD);
	if (rank == 0) {
		
		cout << time2 - time1 << endl;
	}
	
	
	MPI_Scatterv(sendbuf, sendcount, displs, MPI_INT, recvbuf, sendcount[rank], MPI_INT, 0, MPI_COMM_WORLD);
	int max = FindMax(recvbuf, sendcount[rank]);
	MPI_Reduce(&max, &MAX, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

	time2 = MPI_Wtime();

	if (rank == 0) {
		cout << MAX << endl;
		cout << time2 - time1 << endl;
	}


	MPI_Finalize();


	return 0;
}