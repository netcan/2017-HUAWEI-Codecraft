#ifndef _MCMF_SCALING_
#define _MCMF_SCALING_

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <iomanip>
#include <math.h>
#include <deque>
#include "mcmf.h"
#include "deploy.h"
#include <sys/time.h>

using namespace std;
class MCMF_SCALING {
	private:
		static const int MAXN = 200000 + 5;
		int N; // 总节点数
		int A; // 弧数
		int superSource;

		vector<int> startnode;
		vector<int> endnode;
		vector<int> d;
		vector<int> c;
		vector<int> b;
		vector<int> degree;
		int maxdeg;
		int maxdegnode;
		int networkNum, consumerNum, needFlow = 0;

		vector<vector<int>> arcout;
		vector<vector<int>> arcin;

		vector<int> numarcout;
		vector<int> numarcin;
		vector<double> u;
		vector<int> x;
		vector<double> r;
		vector<int> s;

	public:
		double minpot1( vector<int> &v, vector<double> &r ){
			double m;
			m = -r[v.back()];
			v.pop_back();

			while(v.empty() != true){
				m = min(m, -r[v.back()]);
				v.pop_back();
			}
			return m;
		}

		double minpot2( vector<int> &v, vector<double> &r ){
			double m;

			m = r[v.back()];
			v.pop_back();

			while(v.empty() != true){
				m = min(m, r[v.back()]);
				v.pop_back();
			}
			return m;
		}

		void initflow( vector<int> &x, vector<double> const &r, vector<int> const &c, int A){
			int i;
			for(i=0; i < A; i++){
				if( r[i] >= 0){
					x[i] = 0;
				}else{
					x[i] = c[i];
				}
			}
		}

		double cost( vector<int> const &x, vector<int> const &d, int A){
			int i;
			double totalcost = 0;
			for(i = 0; i < A; i++){
				totalcost += double(x[i])*double(d[i]);
			}
			return totalcost;
		}

		bool checkfeas( vector<int> const &x, vector<int> const &c, vector<int> const &s, int A, int N){
			int i;
			for(i=0; i < A; i++){
				if( (x[i] > c[i]) || (x[i]<0) ){ return false; }
			}
			for(i=0; i < N; i++){
				if(s[i] != 0){ return false; }
			}
			return true;
		}

		void printflow( vector<int> &x, int A){
			int i;
			for(i=0; i < A; i++){
				cout<<"x["<<i<<"] = "<<x[i]<<endl;
			}
		}

		double diff_sec(timeval t1, timeval t2){
			return (t2.tv_sec - t1.tv_sec +  (t2.tv_usec - t1.tv_usec)/1000000.0 );
		}

		void AddEdge(int from, int to, int cap, int cost) {
			startnode.push_back(from);
			endnode.push_back(to);
			degree[from] += 1;
			degree[to] += 1;
			c.push_back(cap);
			d.push_back(cost);
			arcin[to].push_back(endnode.size() - 1);
			arcout[from].push_back(startnode.size() - 1);
		}

		void loadGraph(char * topo[MAX_EDGE_NUM], int line_num) {
			int aa, bb, cc, dd;
			sscanf(topo[0], "%d%d%d", &networkNum, &bb, &consumerNum); // 网络节点数量 网络链路数量 消费节点数量
			superSource = networkNum + consumerNum;
			A = bb * 2 + consumerNum;
			N = superSource + 1; // 网络节点+消费节点+1

			b.resize(N, 0);
			degree.resize(N, 0);
			arcin.resize(N);
			arcout.resize(N);
			numarcin.resize(N, 0);
			numarcout.resize(N, 0);

			int i;
			for(i = 2; i < line_num && !isspace(topo[i][0]); ++i);
			for(++i; i < line_num && !isspace(topo[i][0]); ++i);

			for(++i; i < line_num && !isspace(topo[i][0]); ++i) {
				sscanf(topo[i], "%d%d%d%d", &aa, &bb, &cc, &dd); // 链路起始节点ID 链路终止节点ID 总带宽大小 单位网络租用费
				AddEdge(aa, bb, cc, dd);
				AddEdge(bb, aa, cc, dd);
				// printf("u: %d v: %d bandwidth: %d cost: %d\n", a, b, c, d);
			}

			for(++i; i < line_num; ++i) {
				sscanf(topo[i], "%d%d%d", &aa, &bb, &cc); // 消费节点ID 相连网络节点ID 视频带宽消耗需求
				AddEdge(bb, aa + networkNum, cc, 0); // 与网络节点相连
				b[aa + networkNum] = -cc;
				needFlow += cc;
			}
			b[superSource] = needFlow;

			for(int u = 0; u < N; ++u) {
				numarcin[u] = arcin[u].size();
				numarcout[u] = arcout[u].size();
			}
		}

		void showData() {
			// printf("N=%d A=%d\n", N, A);
			// puts("numarcin");
			// for(int u = 0; u < N; ++u)
				// printf("%d ", numarcin[u]);
			// puts("");
			// puts("numarcout");
			// for(int u = 0; u < N; ++u)
				// printf("%d ", numarcout[u]);
			// puts("");
			// puts("d");
			// for(int u = 0; u < N; ++u)
				// printf("%d ", d[u]);
			// puts("");
			// puts("c");
			// for(int u = 0; u < N; ++u)
				// printf("%d ", c[u]);
			// puts("");
			// puts("b");
			// for(int u = 0; u < N; ++u)
				// printf("%d ", b[u]);
			// puts("");
			// puts("degree");
			// for(int u = 0; u < N; ++u)
				// printf("%d ", degree[u]);
			// puts("");
			// puts("startnode");
			// for(int u = 0; u < A; ++u)
				// printf("%d ", startnode[u]);
			// puts("");
			// puts("endnode");
			// for(int u = 0; u < A; ++u)
				// printf("%d ", endnode[u]);
			// puts("");
			// puts("arcin");
			// for(int u = 0; u < N; ++u) {
				// printf("u: %d\n", u);
				// for(size_t i = 0; i < arcin[u].size(); ++i)
					// printf("%d ", arcin[u][i]);
				// puts("");
			// }
			// puts("");
			// puts("arcout");
			// for(int u = 0; u < N; ++u)
				// for(size_t i = 0; i < arcout[u].size(); ++i)
					// printf("%d ", arcout[u][i]);

			// puts("");
		}

		void auction() {

			int i; int j;
			int m; int n;
			double d_max;


			// readdata(filename, N, A, startnode, endnode, d, c, b, degree);

			maxdeg=0;
			for(j=0; j<N; j++){
				if(maxdeg < degree[j]){
					maxdegnode = j;
					maxdeg = degree[j];
				}
			}


			for(i=0; i < A; i++){
				m = startnode[i];
				n = endnode[i];
				arcin[n][numarcin[n]] = i;
				numarcin[n] += 1;
				arcout[m][numarcout[m]] = i;
				numarcout[m] += 1;
			}

			d_max = fabs(double(d[0]));
			for(i=0; i < A; i++){
				if(fabs(double(d[i]))> d_max){
					d_max = fabs(double(d[i]));
				}
			}

			double epsilon = d_max;
			deque<int> nodeq;
			int ibar;

			for(i=0; i < A; i++){
				r[i] = double(d[i]) + u[startnode[i]] - u[endnode[i]];
			}

			bool proceed;
			int k;
			double alpha;
			double alpha1;
			double alpha2;
			int beta;
			vector<int> temp;
			int arc;

			while(epsilon >= 1/( double(N) )){

				initflow(x,r,c,A);

				for(i=0; i < N; i++){
					s[i] = -b[i];
					for(j=0; j < numarcout[i]; j++){
						k = arcout[i][j];
						s[i] += x[k];
					}
					for(j=0; j < numarcin[i]; j++){
						k = arcin[i][j];
						s[i] -= x[k];
					}
					if(s[i]>0){
						nodeq.push_back(i);
					}
				}

				while(nodeq.empty() == false){
					ibar = nodeq.back();
					nodeq.pop_back();
					while(s[ibar] > 0){
						proceed = true;
						temp.clear();

						for(j=0; j<numarcout[ibar] && proceed == true; j++){
							k = arcout[ibar][j];
							if( (r[k] <= epsilon) && (r[k] >= (epsilon/2)) && (x[k] > 0) ){
								arc = k;
								proceed = false;
							}
						}

						if(proceed == false ){
							k = arc;
							beta = min( x[k], s[ibar] );
							x[k] -= beta;
							s[ibar] -= beta;
							s[endnode[k]] += beta;

							if( (s[endnode[k]] > 0) && (s[endnode[k]] - beta <= 0) ){
								nodeq.push_front(endnode[k]);
							}
						}else{
							for( j=0; j<numarcin[ibar] && proceed==true; j++ ){
								k = arcin[ibar][j];
								if( (r[k] >= (-epsilon)) && (r[k] <= (-epsilon/2)) && (x[k] < c[k]) ){
									arc = k;
									proceed = false;
								}
							}

							if(proceed == false ){
								k = arc;
								beta = min( c[k] - x[k], s[ibar]);
								x[k] += beta;
								s[ibar] -= beta;
								s[startnode[k]] += beta;

								if( (s[startnode[k]] > 0) && (s[startnode[k]] - beta) <= 0 ){
									nodeq.push_front(startnode[k]);
								}
							}
						}
						if( proceed==true ){
							for(j=0; j < numarcout[ibar]; j++){
								k = arcout[ibar][j];
								if( x[k] > 0 ){ temp.push_back(k); }
							}
							if(temp.empty() == false){
								alpha1 = minpot1( temp, r );
								temp.clear();

								for(j=0; j < numarcin[ibar]; j++){
									k = arcin[ibar][j];
									if( x[k] < c[k] ){ temp.push_back(k); }
								}
								if(temp.empty() == false){
									alpha2 = minpot2( temp, r );
									temp.clear();
									alpha = min(alpha1, alpha2) + epsilon;
								}else{
									alpha = alpha1 + epsilon;
								}
							}else{
								for(j=0; j < numarcin[ibar]; j++){
									k = arcin[ibar][j];
									if( x[k] < c[k] ){ temp.push_back(k); }
								}
								alpha = minpot2( temp, r) + epsilon;
								temp.clear();
							}

							u[ibar] = u[ibar] + alpha;

							for(j=0; j < numarcout[ibar]; j++){
								k = arcout[ibar][j];
								r[k] += alpha;
							}
							for(j=0; j < numarcin[ibar]; j++){
								k = arcin[ibar][j];
								r[k] -= alpha;
							}
						}
					}
				}
				epsilon = epsilon/2;
			}


			if( checkfeas(x,c,s,A,N) ){
				cout<<"Current flow is feasible"<<endl;
				cout<<"Minimum cost is: "<<cost(x,d,A)<<endl;
			}else{
				cout<<"Current flow not feasible"<<endl;
			}
		}
};


// int main( int argc, char *argv[] ){
// mcmf_scaling ac;
// ac.auction(argc, argv);
// return 0;
// }
//
#endif
