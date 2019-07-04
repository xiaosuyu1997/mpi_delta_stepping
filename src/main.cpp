/*
	Course project for parallel and distributed system in Peking University
	A MPI implementation of Delta stepping algorithm
    Copyright (C) <year>  <name of author>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/


#include <stdlib.h>       // for atoi
#include <iostream>
#include <stdio.h>        // for printf
#include <string.h>
#include"nodearc.h"
#include<math.h>
#include<set>
#include<queue>
#include<vector>
#include<time.h>
#include<mpi.h>
#include <cassert>
#include <list>

using namespace std;

//#define TIME

#define SYN_FREC 1000
#define MAX_BUCKET_CNT 999999
#define MODUL ((long long) 1 << 62)
#define VERY_FAR            9223372036854775807LL // LLONG_MAX

extern double timer();            // in timer.cc: tells time use
extern int parse_gr(int* n_ad, int* m_ad, Node** nodes_ad, Arc** arcs_ad,
	int* node_min_ad, char* problem_name);

extern int parse_ss(int* sN_ad, int** source_array, char* aName);

extern void calculate_bucket(int bucket[], Node* nodes, int n, int Delta);

// Given function
void ArcLen(long cNodes, Node* nodes,
	long long* pMin /* = NULL */, long long* pMax /* = NULL */)
	// finds the max arc length and min arc length of a graph.
	// Used to init for buckets.  pMax or pMin can be NULL.  For fun,
	// and utility, returns MaxArcLen.
{
	Arc* lastArc, * arc;
	long long maxLen = 0, minLen = VERY_FAR;

	// arcs are stored sequentially.  The last arc overall
	// is one before the first arc of the sentinel node.
	lastArc = (nodes + cNodes)->first - 1;
	for (arc = nodes->first; arc <= lastArc; arc++)
	{
		if (arc->len > maxLen)
			maxLen = arc->len;
		if (arc->len < minLen)
			minLen = arc->len;
	}
	if (pMin)* pMin = minLen;
	if (pMax)* pMax = maxLen;
}

// Request datatype
typedef pair<int, long long> req_type;





int main(int argc, char** argv)
{
	int comm_size = 0, comm_id = 0;
	double tm = 0;                               // time recorder

	Arc* arcs = NULL;
	Node* nodes = NULL;
	int source = 0;

	int n = 0, m = 0, nmin = 0, nQ = 0;          // edge num, arc num,  
	int* source_array = NULL;
	long long dist = 0;

	char gName[100], aName[100], oName[100];     // graph file name, auxiliary file name, output file name
	FILE* oFile = NULL;
	long long maxArcLen = 0, minArcLen = 0;
	long long dDist;                             // max tentative path length


	long long Delta = 1024;                      // most important hyperparameter, delta



	// start mpi part
	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &comm_id);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_size);


	if (comm_id == 0)
	{
		if ((argc < 4) || (argc > 5)) {
			fprintf(stderr,
				"Usage: \"%s <graph file> <aux file> <out file>\"\n    or \"%s <graph file> <aux file> <out file> [<Delta>] \"\n", argv[0], argv[0]);
			exit(0);
		}
	}


	strcpy(gName, argv[1]);
	strcpy(aName, argv[2]);
	strcpy(oName, argv[3]);
	oFile = fopen(oName, "a");

	if (argc == 5) {
		Delta = atoi(argv[4]);
	}

	if(comm_id == 0)
	{
	fprintf(stderr, "c ---------------------------------------------------\n");
	fprintf(stderr, "c Chen Ziyu implemented version \n");
	fprintf(stderr, "c ---------------------------------------------------\n");
	}

	parse_gr(&n, &m, &nodes, &arcs, &nmin, gName);

#ifdef DEBUG
	printf("number of nodes: %lld, nmin: %lld\n", n, nmin);
#endif

	parse_ss(&nQ, &source_array, aName);

	ArcLen(n, nodes, &minArcLen, &maxArcLen);      // other useful stats

	// 申请空间，并且initialize其他节点(此处实现的为通过主进程读入图信息，并且传递给其他进程，但是由于
	// 发现存在进程基地址移动问题，node和arc之间的偶联会存在问题，所以弃用，而采用上述各个进程均读入的方式)
	/*
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&minArcLen, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&maxArcLen, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&nQ, 1, MPI_INT, 0, MPI_COMM_WORLD);

	if (comm_id != 0)
	{
		nodes = new Node[n + 1];
		arcs = new Arc[m + 2];
		source_array = new int[nQ];
	}

	MPI_Bcast(nodes, sizeof(Node) * n, MPI_BYTE, 0, MPI_COMM_WORLD);
	MPI_Bcast(arcs, sizeof(Arc) * m, MPI_BYTE, 0, MPI_COMM_WORLD);

	MPI_Bcast(source_array, nQ, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
	*/
	dDist = maxArcLen * (double)(n - 1);

	// sanity check
	if (comm_id == 0)
	{
		if (dDist > VERY_FAR) {
			fprintf(stderr, "Warning: distances may overflow\n");
			fprintf(stderr, "         proceed at your own risk!\n");
		}

		fprintf(stderr, "c Nodes: %24ld       Arcs: %22ld\n", n, m);
		fprintf(stderr, "c MinArcLen: %20lld       MaxArcLen: %17lld\n",
			minArcLen, maxArcLen);
		fprintf(stderr, "c Trials: %23ld\n", nQ);
	}




	//start delta stepping
	int lb = comm_id * n / comm_size;
	int ub = (comm_id + 1) * n / comm_size;

	// 数组记录每个process下的大小，还有displace，用于后续同步D数组
	int* recv_cnt = new int[comm_size];
	int* recv_disp = new int[comm_size + 1];
	recv_disp[0] = 0;
	for (int i = 0; i < comm_size; ++i)
	{
		recv_cnt[i] = (i + 1) * n / comm_size - i * n / comm_size;
		recv_disp[i + 1] = recv_disp[i] + recv_cnt[i];
	}

#ifdef DEBUG
	for (int i = 0; i < comm_size; ++i)
	{
		printf("recv_cnt: %d, recv_disp: %d\n", recv_cnt[i], recv_disp[i]);
	}
#endif

	// 对于每个节点当前的tentative distance通过一个大小为n（节点数）的数组D存储，并且对于每个进程，只保证$i \in [lb, ub)$
	// 的D[i]是在任何更新完成后都正确的。后续会有同步。
	long long* D = new long long[n + 2];
	long long* D_recv = new long long[n + 2];

	Request R(n, comm_size, comm_id);

	list<int>* B = new list<int>[MAX_BUCKET_CNT];     //  一共 dDist/Delta 个bucket

	//set<int> S;     // left out to record nodes with heavy arcs(however, in this implementation we do not distinguish heavy/light nodes)
	int v = 0;

#ifdef TIME
	double find_bucket_time = 0;
	double cal_time = 0;
	double cal_R_time = 0;
	double release_R_time = 0;
	double syn_R_time = 0;
	double syn_D_time = 0;
	double tmp1 = 0;
	double tmp2 = 0;
#endif

	int cur_bucket = 0;
	int max_process_bucket = 0;
	int dest = 0;

	//cout << comm_id << ", flag0" << endl;




	// start timing
	if (comm_id == 0)
		tm = timer();          


	for (int ques_id = 0; ques_id < nQ; ques_id++)
	{
		// init problem
		cur_bucket = 0;
		max_process_bucket = 0;
		dist = 0;
		source = source_array[ques_id] - 1;

		// 负责当前 source 的进程将其插入bucket[0]中
		if (comm_id == (comm_size * (source + 1) - 1) / n)
		{
			B[0].push_back(source);
		}

		// init tentative distance
		for (int i = 0; i < n; ++i)
			D[i] = VERY_FAR;
		D[source] = 0;

		//cout << comm_id << ", flag1" << endl;
#ifdef TIME
		find_bucket_time = 0;
		cal_time = 0;
		cal_R_time = 0;
		release_R_time = 0;
		syn_R_time = 0;
		syn_D_time = 0;
		tmp1 = 0;
		tmp2 = 0;
#endif

		// outer loop
		while (true)
		{

#ifdef TIME
			tmp1 = -clock();
#endif

			//cout << comm_id << ", flag2" << endl;

			// 此处为report中提到的通过数组存储bucket结构的实现，性能低
			/*
			long long tmp_min = VERY_FAR;
			for (int i = 0; i < n; ++i)
			{
				if (D[i] < tmp_min && D[i] >= cur_bucket * Delta)
				{
					tmp_min = D[i];
					if (tmp_min / Delta == cur_bucket) break;
				}

			}

			if (tmp_min == VERY_FAR) break;
			else cur_bucket = tmp_min / Delta;

			for (int i = lb; i < ub; ++i)
			{
				if (D[i] / Delta == cur_bucket)
				{
					B.push(i);
				}

			}
			*/

			// find next bucket, global information, need to synchronize
			bool flag = false;
			while (true)
			{
				bool all_empty = B[cur_bucket].empty();
				MPI_Allreduce(MPI_IN_PLACE, &all_empty, 1, MPI_C_BOOL, MPI_LAND, MPI_COMM_WORLD);
				if (all_empty)    // there is no process with nodes in current buckeet
				{
					cur_bucket += 1;
					// 如何确认终止（曾经记录的插入过的最大bucket < cur_bucket），每个进程分别记录
					bool finish = max_process_bucket < cur_bucket;
					MPI_Allreduce(MPI_IN_PLACE, &finish, 1, MPI_C_BOOL, MPI_LAND, MPI_COMM_WORLD);
					if (finish)
					{
						flag = true;
						break;
					}
				}
				else break;
			}

			if (flag) break;

#ifdef TIME
			find_bucket_time += tmp1 + clock();
			tmp1 = -clock();
#endif

			// 和上述类似的判断是否完成phase loop内部的bucket计算
			bool all_current_bucket_clear = true;
			bool process_bucket_clear = true;


			// inner loop(phase loop)
			while (true)
			{

#ifdef TIME
				tmp2 = -clock();
#endif

				process_bucket_clear = B[cur_bucket].empty();
				//cout << "local: " << comm_id << " " << process_bucket_clear << endl;

				MPI_Allreduce(&process_bucket_clear, &all_current_bucket_clear, 1, MPI_C_BOOL, MPI_LAND, MPI_COMM_WORLD);
				//cout << "all:" << comm_id << " "<< all_current_bucket_clear << endl;
				if (all_current_bucket_clear) break;

				R.clear();

				while (!B[cur_bucket].empty())
				{
					int pop_node = B[cur_bucket].front();
					//S.insert(pop_node);
					B[cur_bucket].pop_front();

					// generate requests
					for (Arc* arc = (nodes + pop_node)->first; arc < (nodes + pop_node + 1)->first; ++arc)
					{
						//cout << comm_id << ", tmpflag " << pop_node << " " << arc << " " << arcs << endl;
						dest = arc->head - nodes;
						//cout << comm_id << ", tmpflag " << dest << endl;
						if (D[pop_node] + arc->len < D[dest])
						{
							R.insert(dest, D[pop_node] + arc->len);
						}
					}
				}

#ifdef TIME
				cal_R_time += tmp2 + clock();
				tmp2 = -clock();
#endif
				//cout << comm_id << ", flag3" << endl;
				//synchronization
				req_type* recv_linear = NULL;
				int cnt;

				R.synchronize(&cnt, &recv_linear);

#ifdef TIME
				syn_R_time += tmp2 + clock();
				tmp2 = -clock();
#endif
				
				//cout << comm_id << ", flag4" << endl;
				// relax process
				for (int i = 0; i < cnt; ++i)
				{
					if (D[recv_linear[i].first] > recv_linear[i].second)
					{

						//删除节点
						v = recv_linear[i].first;
						if (D[v] != VERY_FAR)   // has been inserted into one bucket
						{
							B[D[v] / Delta].remove(v);     // erase if existed (maybe not, but only because it's stored in the current working bucket, or have been once inserted in this loop)
						}

						// 更新权重
						D[v] = recv_linear[i].second;

						// 插入节点
						B[D[v] / Delta].push_back(v);

						max_process_bucket = max((long long)max_process_bucket, D[v] / Delta);
					}

				}

#ifdef TIME
				release_R_time += tmp2 + clock();
#endif

			}


#ifdef TIME
			cal_time += tmp1 + clock();
			tmp1 = -clock();
#endif

			// 但是为了减少生成的request数量，我们在MAX_BUCKET_CNT这个超参控制的外层循环次数后，
			// 就对所有进程内的D数组进行同步（一次Allgatherv操作，每个数组提供[lb, ub)，放置到新的数组中）。
			if (cur_bucket % SYN_FREC == 0)
			{
				MPI_Allgatherv(D + lb, ub - lb, MPI_LONG_LONG, D_recv, recv_cnt, recv_disp, MPI_LONG_LONG, MPI_COMM_WORLD);

				swap(D, D_recv);
			}


#ifdef TIME
			syn_D_time += tmp1 + clock();
#endif

			cur_bucket++;
		}
		MPI_Allgatherv(D + lb, ub - lb, MPI_LONG_LONG, D_recv, recv_cnt, recv_disp, MPI_LONG_LONG, MPI_COMM_WORLD);

		swap(D, D_recv);
		//cout << comm_id << ", flag6" << endl;
#ifdef CHECKSUM
		if (comm_id == 0)
		{
			for (int i = 0; i < n; ++i)
				if (D[i] != VERY_FAR)
				{
					//printf("%d\n", i + 1);
					dist = (dist + (D[i] % MODUL)) % MODUL;
				}
			//cout << dist << endl;
			fprintf(oFile, "d %lld\n", dist);
			//cout << comm_id << ", flag7" << endl;
		}
#endif

	}

#ifndef CHECKSUM
	if (comm_id == 0)
	{
		fprintf(oFile, "f %s %s\n", gName, aName);
		tm = (timer() - tm);   // finish timing, ms
		/* *round* the time to the nearest .01 of ms */
		fprintf(stderr, "c Time (ave, ms): %18.2f\n",
			1000.0 * tm / (float)nQ);

		fprintf(oFile, "g %ld %ld %lld %lld\n",
			n, m, minArcLen, maxArcLen);
		fprintf(oFile, "t %f\n", 1000.0 * tm / (float)nQ);
	}
#endif




	fclose(oFile);

	delete[] recv_cnt;
	delete[] recv_disp;
	delete[] D;
	delete[] D_recv;
	delete[] B;

	MPI_Finalize();
	return 0;
}
