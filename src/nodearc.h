/* nodearc.h  by Andrew Goldberg, started 5/25/01
 *     Contains the definition of the node and arc class needed
 *     for all the code.
 */
# include<cassert>
# include<mpi.h>
# include<vector>


#ifndef NODEARC_H
#define NODEARC_H
using namespace std;

typedef pair<int, long long> req_type;
//#define NOT_SCANNED       0   // possible data structures the node can be in
//#define IN_HEAP       1
//#define IN_F          2   // in smartq F
//#define IN_BUCKETS    4
//#define IN_SCANNED    5

typedef struct Node;      // so the arc knows about it

typedef struct Arc {
  long long len;         // arc length; long long may be an overkill
  struct Node *head;     // where the arc ends up
} Arc;

typedef struct Node {
  long long dist;        // tentative shortest path length to some node
  Arc   *first;          // first outgoing arc
  
  Arc* firstlight;
  Arc* firstheavy;
  
  //struct Node *parent;   // parent on the heap
  //char state;   // what data structure we're in:  IN_* (above)
  
  //unsigned int tStamp;

//  struct {
//    struct Node *next;          // next in bucket
//    struct Node *prev;          // prev in bucket
//    void *bucket;               // you should cast this
//#ifndef MLB
//    long long caliber;          // minimum incoming arc length
//#endif
//  } sBckInfo;                   // FOR SMART BUCKETS
} Node;

class Request
{
public:
	void insert(int node_id, long long dist);
	void clear();
	void synchronize(int* cnt, req_type** reqs);

	int n;
	int comm_size, comm_id;

	int* send_counter;  // 记录需要分到每个comm_id 进程的request数目
	int* recv_counter;  // 记录每个进程从其他进程接受的request数目
	vector<req_type>* req_vec;  // 记录产生的所有request，等待线性化后交换
	req_type* send_linear;      // 从req_vec线性化后等待发送的buffer
	int send_linear_len;        // 总长度

	req_type* recv_linear;      // 接收的buffer
	int recv_linear_len;        // 总长度

	MPI_Datatype mytype;

	Request(int _n, int _size, int _id);
	~Request();
};



#endif
