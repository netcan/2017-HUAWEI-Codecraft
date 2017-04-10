/*************************************************************************
	> File Name: spfa.h
	  > Author: Netcan
	  > Blog: http://www.netcan666.com
	  > Mail: 1469709759@qq.com
	  > Created Time: 2017-03-22 Wed 12:33:02 CST
 ************************************************************************/

#include <cstdio>
#include <assert.h>
#include <cstring>
#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <deque>
#include <vector>
#include <unordered_set>
#include "deploy.h"
using namespace std;

#ifndef __MCMF__
#define __MCMF__

class MCMF{
	private:
		struct Edge{
			int from, to, cap, flow ,cost, oldCost;
			Edge() {}
			Edge(int from, int to, int cap, int flow, int cost):from(from), to(to), cap(cap), flow(flow), cost(cost), oldCost(cost) {}
		};
		struct Server {
			int level, outFlow, cost;
			Server(int level, int outFlow, int cost): level(level), outFlow(outFlow), cost(cost) {}
			Server() {
				level = outFlow = cost = 0;
			}
			bool operator<(const Server &b) const {
				if(this->outFlow != b.outFlow) return this->outFlow < b.outFlow;
				else return this->cost < b.cost;
			}
			bool operator<(int flow) const {
				return this->outFlow < flow;
			}
		};

		static const int N = 30000+5;
		static char topo[50000*1000*6]; // 网络路径数量不得超过300000条, 单条路径的节点数量不得超过10000个, 所有数值必须为大于等于0的整数，数值大小不得超过1000000。

		struct {
			int q[N];
			int tail, head;
			void reset() {
				tail = head = 0;
			}
			bool empty() {
				return tail == head;
			}
			void push(int x) {
				q[tail] = x;
				tail = (tail + 1) % N;
			}
			void pop() {
				if(empty()) return;
				head = (head + 1) % N;
			}
			int front() {
				return q[head];
			}

		} queue;

		int Vn, superSource, superSink; // 总节点数，超级源点/汇点，需要的流量
		int d[N], f[N], p[N]; // 最小费用，当前流量，父节点（增广路径），中间变量
		vector<Server> servers; // 服务器
		Server maxFlowServer;
		int deployCost[10000+5]; // 节点部署费用
		int flow[10000+5]; // 节点流出的流量

		bool vis[N]; // 标记指针
		bool mcmfMethod = 1; // 0为BellmanFord最小费用流，1为ZKW最小费用流算法
		pair<int, vector<vector<int>>> solutionPath; // 当前可行解，路径

		bool BellmanFord(int s, int t, int &flow, int &cost);
		// ZKW算法
		int aug(int u, int minFlow, int &tmpCost, int &cost);
		bool modLabel(int &tmpCost);

		inline void reset() { // 还原初始状态，删除源点
			for(int i = G[superSource].size() - 1; i >= 0; --i) {
				Edge &eS = edges[G[superSource][i]]; // 获取指向cdn的边
				Edge &eV = edges[G[eS.to].back()]; // 获取cdn指向虚拟节点的边
				G[eS.to].pop_back(); // 删除虚拟节点的指针
				for(int j = G[eV.to].size() - 1;  j >= 0; --j)
					G[edges[G[eV.to][j]].to].pop_back();
				G[eV.to].clear();
			}

			G[superSource].clear(); // 清空超源的指针
			for(int i=edges.size(); i > edgeNum; --i) edges.pop_back(); // 删除超源的边
			for(size_t i = 0; i < edges.size(); ++i) {
				if(mcmfMethod) edges[i].cost = edges[i].oldCost; // 恢复费用
				edges[i].flow = 0; // 重置流量
			}
		}

		int findPath(vector<int> & tmpPath, int u, int minFlow, int totalFlow);
		void getPath(int cost);

		inline int minCost(const unordered_set<int> &cdn) { // 调用setCDN后再调用minCost!! 注意不能连续调用多次minCost!!!
			int cost = 0, flow = 0;

			if(mcmfMethod == 0) { // BellmanFord算法
				while (BellmanFord(superSource, superSink, flow, cost));
			} else { // zkw算法
				int tmpCost = 0;

				do
					do
						bzero(vis, sizeof(vis));
					while(aug(superSource, INF, tmpCost, cost));
				while(modLabel(tmpCost));

				// SLF优化
				/*
				while(modLabel(tmpCost))
					do bzero(vis, sizeof(vis));
					while(aug(superSource, INF, tmpCost, cost));
				*/

				for (size_t i = 0; i < G[superSource].size(); i++)
					flow += edges[G[superSource][i]].flow;
			}

			if(flow < needFlow) return -1;
			for(auto c: cdn) cost += deployCost[c];
			cost += G[superSource].size() * costPerCDN;
			if(cost < solutionPath.first) getPath(cost); // 更新方案
			return cost;
		}


		inline void setCdn(const unordered_set<int> & cdn) {
			reset();

			for(int u: cdn) {
				AddEdge(superSource, u, MCMF::INF, 0);
				AddEdge(u, u+Vn, maxFlowServer.outFlow, 0); // 拆点
				for(size_t i = 0; i < G[u].size(); ++i) {
					const Edge &e = edges[G[u][i]];
					if(e.to == superSource || e.to == u+Vn || e.cost < 0) continue;
					AddEdge(u+Vn, e.to, e.cap, e.cost);
					edges[G[u][i]].cost = 1000000; // 调大些
					edges[G[u][i] ^ 1].cost = 1000000;
				}
			}
		}
	public:
		inline bool isConsumer(int u) {
			return u >= networkNum && u < superSource;
		}

		vector<int> G[N]; // 图
		vector<Edge> edges; // 边集
		int networkNum, edgeNum, consumerNum, needFlow, costPerCDN = 0;
		static const int INF = 0x3f3f3f3f;

		MCMF() {
			needFlow = 0;
		};
		void AddEdge(int from, int to, int cap, int cost);
		void showSolution() const;
		void loadGraph(char * topo[MAX_EDGE_NUM], int line_num);
		const char* outputPath();


		inline int minCost_Set(const unordered_set<int> &cdn) {
			setCdn(cdn);
			return minCost(cdn);
		}

};
extern MCMF mcmf;

#endif
