#include "mpi.h"
#include<math.h>
#include <string>
#include<iostream>
#include<stdio.h>
using namespace std;

int main_BCast(int argc, char ** argv) {
	void * sendbuf = new int[160], *recvbuf = new int[10];

	int sendcount = 10, recvcount = 10, root = 0;
	int procNum, procRank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &procNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
	if (procRank == 0) {
		for (size_t i = 0; i < 160; i++)
		{
			((int*)sendbuf)[i] = i + 1;
			std::cout << ((int*)sendbuf)[i] << " ";
		}
		std::cout << std::endl;
	}


	int dataBlockCount = procNum;


	float l = log2(procNum);
	for (size_t i = 0; i < l; i++)
	{
		MPI_Barrier(MPI_COMM_WORLD);
		for (size_t j = 0; j < (1 << (i + 1)); j++)
		{
			if (procRank == j)
			{
				if (procRank < (1 << i) && procRank + (1 << i) < procNum)
				{
					std::cout << " send" << " ";
					std::cout << "i= " << i << " ";
					std::cout << "j= " << j << " ";
					std::cout << "procRank " << procRank << " ";
					std::cout << "procRank+(1<<i) " << procRank + (1 << i) << " ";
					std::cout << "dataBlockCount " << dataBlockCount << " ";

					int dataBlockCountSend = dataBlockCount / 2 + dataBlockCount % 2;
					MPI_Send(&dataBlockCountSend, 1, MPI_INT, procRank + (1 << i), 1, MPI_COMM_WORLD);
					std::cout << "dataBlockCountSend " << dataBlockCountSend << " ";
					MPI_Send((int*)sendbuf + (dataBlockCount / 2 + dataBlockCount % 2) * sendcount,
						(dataBlockCount / 2) * sendcount,
						MPI_INT, procRank + (1 << i), 0, MPI_COMM_WORLD);
					dataBlockCount = dataBlockCount / 2 + dataBlockCount % 2;

					std::cout << "dataBlockCount " << dataBlockCount << " ";
					std::cout << std::endl;
				}
				if (procRank >= (1 << i))
				{
					std::cout << " recv" << " ";
					std::cout << "i= " << i << " ";
					std::cout << "j= " << j << " ";
					std::cout << "procRank " << procRank << " ";
					std::cout << "procRank-(1<<i) " << procRank - (1 << i) << " ";
					std::cout << "dataBlockCount " << dataBlockCount << " ";
					MPI_Recv(&dataBlockCount, 1, MPI_INT, procRank - (1 << i), 1, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
					std::cout << "dataBlockCount " << dataBlockCount << " ";
					MPI_Recv(sendbuf, dataBlockCount  * recvcount, MPI_INT, procRank - (1 << i), 0,
						MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					//dataBlockCount = dataBlockCount / 2;
					//std::cout << "dataBlockCount " << dataBlockCount << " ";
					std::cout << std::endl;
				}
			}
		}


		//std::cout << "barrier ";
		MPI_Barrier(MPI_COMM_WORLD);
		std::cout << "circle finished " << procRank;
		std::cout << std::endl;
		std::cout << std::endl;
	}
	std::cout << std::endl;

	memcpy(recvbuf, sendbuf, recvcount);

	for (size_t j = 0; j < procNum; j++)
	{
		if (procRank == j) {
			std::cout << " proc:" << procRank << " recvbuf= ";
			for (size_t i = 0; i < recvcount; i++)
				std::cout << ((int*)sendbuf)[i] << " ";
		}
	}

	if (procRank != 0)
		delete sendbuf;
	MPI_Finalize();
	//getchar();
	return 0;
}


void Scatter(const void * sendbuf, int sendcount, MPI_Datatype sendtype,
	void * recvbuf, int recvcount, MPI_Datatype recvtype,
	int root, MPI_Comm comm)
{
	int procRank, procNum;
	MPI_Comm_size(comm, &procNum);
	MPI_Comm_rank(comm, &procRank);

	int sendtypesize, recvtypesize;
	MPI_Type_size(sendtype, &sendtypesize);
	MPI_Type_size(recvtype, &recvtypesize);


	if (procRank % 2) 
	{
		MPI_Recv(recvbuf, recvcount, recvtype, procRank - 1, procRank, comm, MPI_STATUSES_IGNORE);
	}
	else 
	{
		int i = 2;
		while (i < procNum) i = i << 1;
			
		int tmp = procRank;
		while (tmp % i != 0) { tmp = tmp % i; i = i >> 1; }
		//cout << "rank: " << procRank << " i: " << i << endl;
		
		char *buf;
		if (procRank != root) {
			buf = new char[i*recvtypesize*recvcount];
			//cout << "rank: " << rank << " i*recvtypesize*recvcount: " << i*recvtypesize*recvcount << endl;
			MPI_Recv(buf, i*recvcount, recvtype, procRank - i, procRank, comm, MPI_STATUSES_IGNORE);				
		}
		else 
		{
			buf = (char*)sendbuf;
		}
		
		
		memcpy(recvbuf, buf, sendcount*sendtypesize);

		//cout << "rank: " << rank << " sendcount * sendtypesize: " << sendcount * sendtypesize << endl;
		
		i = i >> 1;
		//cout << "rank: " << procRank << " i= " << i << endl;
		while (i >= 1)
		{
			cout << "rank: " << procRank << " i= " << i << endl;
			if (procRank + i < procNum)
				MPI_Send(buf + i * sendcount * sendtypesize, i*sendcount, sendtype, procRank + i, procRank + i, comm);
			//cout << "rank: " << rank << "i * sendcount * sendtypesize= " << i * sendcount * sendtypesize <<" i * sendcount ="<< i * sendcount<<" rank+i ="<<rank+1<< endl;
			i = i >> 1;
		}

		if (procRank != root)
			delete[] buf;
	}
}

void linScatter(const void * sendbuf, int sendcount, MPI_Datatype sendtype,
	void * recvbuf, int recvcount, MPI_Datatype recvtype,
	int root, MPI_Comm comm) 
{
	int procRank, procNum;
	MPI_Comm_size(comm, &procNum);
	MPI_Comm_rank(comm, &procRank);
	if (procRank == 0) {
		for (size_t i = 1; i < procNum; i++)
			MPI_Send(sendbuf, sendcount, MPI_BYTE, i, 0, comm);
	}
	else {
		MPI_Recv(recvbuf, recvcount, MPI_BYTE, 0, 0, comm, MPI_STATUSES_IGNORE);
	}
	
}

int main(int argc, char ** argv) {
	void * sendbuf = new char[40], *recvbuf = new char[10];
	double time1, time2;
	int sendcount = 10, recvcount = 10, root = 0;
	int procNum, procRank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &procNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
	if (procRank == 0) {
		for (size_t i = 0; i < 40; i++)
		{
			((char*)sendbuf)[i] = i + 48;
			std::cout << ((char*)sendbuf)[i] << " ";
		}
		std::cout << std::endl;
	}
	time1 = MPI_Wtime();
	Scatter(sendbuf, sendcount, MPI_BYTE, recvbuf, recvcount, MPI_BYTE, 0, MPI_COMM_WORLD);
	time2 = MPI_Wtime();
	for (size_t j = 0; j < procNum; j++)
	{
		if (procRank == j) {
			std::cout << " proc:" << procRank << " recvbuf= ";
			for (size_t i = 0; i < recvcount; i++)
				std::cout << ((char*)recvbuf)[i] << " ";
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	std::cout << std::endl;
	if (procRank == 0) std::cout << time2 - time1 << std::endl;
	MPI_Finalize();
	return 0;
}