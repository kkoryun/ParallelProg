#define _CRT_SECURE_NO_WARNINGS
#define MSMPI_NO_DEPRECATE_20
#include<stdlib.h>
#include<stdio.h>
#include<mpi.h>
#include"struct.h"
#include<iostream>
#include<fstream>
#include<string>
#include"algorithm.h"

void fileReader(string filename, double * buf, int size) {
	FILE * file = fopen(filename.c_str(), "r");
	string filebuf = "";
	char c = 0;
	int i = 0;

	while (/*!feof(file) && */i < size)
	{
		c = 0;
		while (c != ',' /*&& !feof(file)*/)
		{
			fread(&c, sizeof(char), 1, file);
			//		switch (c)
			//		{
			//		case '\0':break;
			//		case '\n':break;
			//		case ',': break;
			//		case ' ': break;
			//	
			//		default: {filebuf += c;
			//			break;}
			//		}
		}
		//	if(filebuf[0]=='"'){
		//		string s = filebuf;
		//		filebuf = "";
		//		int i = 1;
		//		while (s[i]!='"')
		//		{
		//			filebuf += s[i];
		//			i++;
		//		}
		//	}

		//	buf[i] = stod(filebuf);
		i++;
		//	filebuf = "";
	}
	//fclose(file);
}
void partition(int size, int part, int blocksize, int * counts, int * displs) {
	int blockcount = size / blocksize;
	int partblockcount = blockcount / part;
	//cout << "partblockcount: " << partblockcount << endl;
	int remaind = blockcount % part;
	//cout << "remaind: " << remaind << endl;
	for (size_t i = 0; i < part - 1; i++)
	{
		counts[i] = partblockcount * blocksize;
		displs[i] = i * partblockcount * blocksize;
	}
	counts[part - 1] = (partblockcount + remaind)*blocksize;
	displs[part - 1] = (part - 1) * partblockcount * blocksize;

}


const int pointSize = 50, setSize = 100000, K = 300, POINTSIZE = pointSize + 1;
int main(int argc, char ** argv) {
#pragma region DATA
	double * sendbuf = 0;
	double * recvbuf = 0;
	double * newpoint = 0;
	double * metrics = 0;
	double * weight = 0;
	double   time1 = 0;
	double   time2 = 0;
	int    * sendcount = 0;
	int    * displs = 0;
	int    * minIndex = 0;
	pointd * trainingSet = 0;
	pointd   newpnt;

#pragma endregion

	MPI_Init(&argc, &argv);


	int size, rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	newpoint = new double[POINTSIZE];
#pragma region RANK0

	if (rank == 0)
	{

		sendbuf = new double[setSize * POINTSIZE];
		for (size_t i = 0; i < setSize * POINTSIZE; i++)
		{
			sendbuf[i] = ((double)rand()) / RAND_MAX;
		}
		for (size_t i = 0; i < pointSize; i++)
		{
			newpoint[i] = ((double)rand()) / RAND_MAX;
		}
		time1 = MPI_Wtime();
	}
#pragma endregion


#pragma region ALLOPERATIONS

	MPI_Bcast(newpoint, POINTSIZE, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	newpnt = pointd(newpoint, pointSize, -1);

	//get training point buf for all proc
	sendcount = new int[size];
	displs = new int[size];
	partition(setSize * POINTSIZE, size, POINTSIZE, sendcount, displs);
	recvbuf = new double[sendcount[rank]];
	MPI_Scatterv(sendbuf, sendcount, displs, MPI_DOUBLE, recvbuf, sendcount[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
	//
	//create training point on proc
	int blockcount = sendcount[rank] / POINTSIZE;
	trainingSet = new pointd[blockcount];
	for (size_t i = 0; i < blockcount; i++) {
		trainingSet[i] = pointd(recvbuf + POINTSIZE * i, pointSize, *(recvbuf + pointSize) + 0.5);
	}
	//

	//calc metrics for all point on proc 
	metrics = new double[blockcount];

	for (size_t i = 0; i < blockcount; i++) {
		metrics[i] = metric<double>(trainingSet[i], newpnt, pointSize);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	//find K point's index with minimum metric
	minIndex = new int[K];
	findKMin(metrics, blockcount, K, minIndex);


	double * tmprecvbuf = 0;
	if (rank == 0)
	{
		tmprecvbuf = new double[POINTSIZE*K*size];
	}

	double * tmpsendbuf = new double[POINTSIZE*K];
	for (size_t i = 0; i < K; i++)
	{
		memcpy(tmpsendbuf + POINTSIZE*i, trainingSet[minIndex[i]].data, pointSize * sizeof(double));
		tmpsendbuf[pointSize + POINTSIZE*i] = trainingSet[minIndex[i]].Class;
	}
	MPI_Gather(tmpsendbuf, POINTSIZE*K, MPI_DOUBLE, tmprecvbuf, POINTSIZE*K, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	delete[] tmpsendbuf;
#pragma endregion

#pragma region RANK0
	if (rank == 0)
	{
		pointd * applicants = new pointd[K*size];
		for (size_t i = 0; i < K*size; i++) {
			applicants[i] = pointd(tmprecvbuf + POINTSIZE * i, pointSize, *(tmprecvbuf + pointSize + POINTSIZE * i));
		}
		////
		////
		cout << endl;
		double * applicantsmetric = new double[K*size];
		for (size_t i = 0; i < K*size; i++) {
			applicantsmetric[i] = metric<double>(applicants[i], newpnt, pointSize);
		}
		////
		int * applicantminindex = new int[K];
		findKMin(applicantsmetric, K*size, K, applicantminindex);
		////
		weight = new double[K];
		for (size_t i = 0; i < K; i++) {
			weight[i] = 1.0 / (applicantsmetric[applicantminindex[i]] * applicantsmetric[applicantminindex[i]]);
		}
		time2 = MPI_Wtime();
		
		cout<< "TIME"<< time2 - time1 << endl;
		delete[] tmprecvbuf;
		delete[] applicants;
		delete[] applicantsmetric;
		delete[] applicantminindex;
		delete[] weight;
	}
#pragma endregion



	MPI_Finalize();

	return 0;
}
