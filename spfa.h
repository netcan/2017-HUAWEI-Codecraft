/*************************************************************************
	> File Name: spfa.h
	  > Author: Netcan
	  > Blog: http://www.netcan666.com
	  > Mail: 1469709759@qq.com
	  > Created Time: 2017-03-22 Wed 12:33:02 CST
 ************************************************************************/

#include <cstdio>
#include <cstring>
#include <algorithm>
#include <cctype>
#include <queue>
#include <cstdlib>
#include <vector>
#include <unordered_set>
#include "deploy.h"
using namespace std;

#ifndef __MCMF__
#define __MCMF__

class MCMF{
	private:
		struct Edge{
			int from, to, cap, flow ,cost;
			Edge() {}
			Edge(int from, int to, int cap, int flow, int cost):from(from), to(to), cap(cap), flow(flow), cost(cost) {}
		};
		static const int N = 1500+5;
		static char topo[50000*1000*6];

		int superSource, superSink; // 总节点数，超级源点/汇点，需要的流量
		int d[N], f[N], p[N]; // 最小费用，当前流量，父节点（增广路径），中间变量
		bool vis[N]; // 标记指针
		pair<int, vector<vector<int>>> solutionPath; // 当前最优费用，路径

		bool BellmanFord(int s, int t, int &flow, int &cost);

		inline void reset() { // 还原初始状态，删除源点
			G[superSource].clear();
			edges = oldEdges;
		}

		int findPath(vector<int> & tmpPath, int u, int minFlow, int totalFlow);
		void getPath(int cost);
	public:
		vector<int> G[N]; // 图
		vector<Edge> edges, oldEdges; // 边集，边集备份
		int networkNum, edgeNum, consumerNum, costPerCDN, needFlow;
		static const int INF = 0x3f3f3f3f;

		MCMF() {
			needFlow = 0;
		};
		MCMF(int superSource, int superSink, int networkNum, int edgeNum, int consumerNum ,int costPerCDN, int needFlow = 0):
			superSource(superSource), superSink(superSink), networkNum(networkNum), edgeNum(edgeNum),
			consumerNum(consumerNum), costPerCDN(costPerCDN), needFlow(needFlow) {}
		void AddEdge(int from, int to, int cap, int cost);
		void showSolution() const;
		void loadGraph();
		void loadGraph(char * topo[MAX_EDGE_NUM], int line_num);
		const char* outputPath();

		inline int minCost() { // 调用setCDN后再调用minCost!! 注意不能连续调用多次minCost!!!
			int flow = 0, cost = 0;
			while (BellmanFord(superSource, superSink, flow, cost));
			cost += G[superSource].size() * costPerCDN;

			if(flow < needFlow) return -1;
			else {
				if(cost < solutionPath.first) getPath(cost); // 更新方案

				// for(int i=0; i < networkNum + consumerNum +2; ++i) {
					// for(size_t j = 0; j < G[i].size(); ++j)
						// printf("%d->%d flow: %d\n", edges[G[i][j]].from, edges[G[i][j]].to, edges[G[i][j]].flow);
				// }
				return cost;
			}
		}
		inline void setCdn(const unordered_set<int> & cdn) {
			reset();
			for(int x: cdn)
				AddEdge(superSource, x, MCMF::INF, 0);
		}

};
extern MCMF mcmf;

#endif
