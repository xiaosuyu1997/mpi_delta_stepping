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

# include"nodearc.h"

using namespace std;



Request::Request(int _n, int _size, int _id) :
	n(_n), comm_size(_size), comm_id(_id)
{
	send_counter = new int[comm_size];
	recv_counter = new int[comm_size];
	req_vec = new vector<req_type>[comm_size];
	send_linear = new req_type[1000];    // 临时先用着，如果不够再重新申请
	send_linear_len = 1000;

	recv_linear = new req_type[1000];
	recv_linear_len = 1000;

	MPI_Type_contiguous(16, MPI_BYTE, &mytype);
	MPI_Type_commit(&mytype);


}

Request::~Request()
{
	delete[] send_counter;
	delete[] recv_counter;
	delete[] send_linear;
	MPI_Type_free(&mytype);
}

void Request::insert(int node_id, long long dist)
{
	// 记录存储此node的进程
	int process = (comm_size * (node_id + 1) - 1) / n;
	send_counter[process] += 1;
	req_vec[process].push_back(req_type(node_id, dist));
}


// clear is done every inner loop
void Request::clear()
{
	for (int i = 0; i < comm_size; ++i)
	{
		send_counter[i] = 0;
		req_vec[i].clear();
	}
}



void Request::synchronize(int* cnt, req_type** reqs)
{
	int cnt_recorded = 0;
	for (int i = 0; i < comm_size; ++i)
		cnt_recorded += send_counter[i];

	if (cnt_recorded > send_linear_len)
	{
#ifdef DEBUG
		printf("req length is not enough, from %d to %d\n", current_req_linear_len, cnt_recorded);
#endif
		delete[] send_linear;
		send_linear_len = cnt_recorded + 100;
		send_linear = new req_type[send_linear_len];
	}


	int cnt_in = 0;
	for (int i = 0; i < comm_size; ++i)
	{
		for (vector<req_type>::iterator it = req_vec[i].begin(); it < req_vec[i].end(); ++it)
			send_linear[cnt_in++] = *it;
	}

	MPI_Alltoall(send_counter, 1, MPI_INT, recv_counter, 1, MPI_INT, MPI_COMM_WORLD);

	// calculate displace
	int* senddis = new int[comm_size];
	int* recvdis = new int[comm_size];
	senddis[0] = 0;
	recvdis[0] = 0;
	for (int i = 1; i < comm_size; ++i)
	{
		senddis[i] = senddis[i - 1] + send_counter[i - 1];
		recvdis[i] = recvdis[i - 1] + recv_counter[i - 1];
	}

	// affirm recv array is enough
	if (recv_linear_len < recvdis[comm_size - 1] + recv_counter[comm_size - 1])
	{
		delete[] recv_linear;
		recv_linear_len = recvdis[comm_size - 1] + recv_counter[comm_size - 1] + 100;
		recv_linear = new req_type[recv_linear_len];
	}

	MPI_Alltoallv(send_linear, send_counter, senddis, mytype, recv_linear, recv_counter, recvdis, mytype, MPI_COMM_WORLD);


	*reqs = recv_linear;
	*cnt = recvdis[comm_size - 1] + recv_counter[comm_size - 1];
	//cout << "process" << comm_id << "cnts:" << *cnt << endl;

	delete[] senddis;
	delete[] recvdis;
}
