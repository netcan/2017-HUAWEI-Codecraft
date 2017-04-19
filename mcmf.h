/*************************************************************************
	> File Name: spfa.h
	  > Author: Netcan
	  > Blog: http://www.netcan666.com
	  > Mail: 1469709759@qq.com
	  > Created Time: 2017-03-22 Wed 12:33:02 CST
 ************************************************************************/
#ifndef __MCMF__
#define __MCMF__

#include <cstdio>
#include <assert.h>
#include <cstring>
#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <deque>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include "deploy.h"
using namespace std;


extern bool runing;

class MCMF{
	private:
		struct Edge{
			int to, cap, flow ,cost, oldCost;
			Edge() {}
			Edge(int to, int cap, int flow, int cost): to(to), cap(cap), flow(flow), cost(cost), oldCost(cost) {}
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

		static const int N = 20000+5;
		static char topo[50000*1000*6]; // 网络路径数量不得超过300000条, 单条路径的节点数量不得超过10000个, 所有数值必须为大于等于0的整数，数值大小不得超过1000000。

		int Vn, superSource, superSink; // 总节点数，超级源点/汇点，需要的流量
		size_t minCostCdnGap = 50; // 当cdn小于这个数的时候，进行贪心降档
		int d[N];
		bool vis[N]; // 标记数组
		vector<Server> servers; // 服务器
		Server maxFlowServer;
#ifdef _DEBUG
		int realMinCost = INF; // 保存真实的最小费用，最后打印，调试用
#endif

		pair<int, vector<vector<int>>> solutionPath; // 当前可行解，路径

		// ZKW算法
		int aug(int u, int minFlow, int &tmpCost, int &cost);
		bool modLabel(int &tmpCost);

		inline void reset() { // 还原初始状态，删除源点
			for(size_t i = 0; i < G[superSource].size(); ++i)
				G[edges[G[superSource][i]].to].pop_back(); // 删除链接超源的边
			for(int i=G[superSource].size() * 2; i > 0; --i) edges.pop_back(); // 删除超源的边
			for(size_t i = 0; i < edges.size(); ++i) {
				edges[i].cost = edges[i].oldCost;
				edges[i].flow = 0; // 重置流量
			}
			G[superSource].clear();
		}
		inline void calcEvaluation() { // 评估函数
			for(int u = 0; u < networkNum; ++u)
				nodes[u].evaluation =  nodes[u].nodeFlow * 1000.0 / nodes[u].deployCost;
		}

		int findPath(vector<int> & tmpPath, int u, int minFlow, int totalFlow);
		void getPath(int cost);
		inline int pathFlowCost() { // 路径流量费
			int cost = 0, flow = 0;
			int tmpCost = 0;
			do {
				int f;
				do {
					memset(vis, 0, sizeof(vis[0]) * Vn);
					f = aug(superSource, INF, tmpCost, cost);
					flow += f;
				}
				while(f);
			}
			while(modLabel(tmpCost));
			if(flow < needFlow) return -1;
			return cost;

			// SLF优化
			/*
			   while(modLabel(tmpCost))
			   do bzero(vis, sizeof(vis));
			   while(aug(superSource, INF, tmpCost, cost));
			*/
		}

		inline int minCost(const unordered_set<int> &cdn) { // 调用setCDN后再调用minCost!! 注意不能连续调用多次minCost!!!
			int cost = pathFlowCost(), cdnCost = 0, cdnFlow = 0, minCdnFlowCost = INF;
			if(cost == -1) return -1;

			vector<pair<int, int>> diff; // 存放差/节点名
			unordered_map<int, int> eId; // 存放超源到cdn的边的下标
			bool downShift = false;
			if(cdn.size() < minCostCdnGap) downShift = true;

			for (size_t i = 0; i < G[superSource].size(); i++) { // 降档
				Edge &e = edges[G[superSource][i]];

				vector<Server>::iterator it;
				if( (it = lower_bound(servers.begin(), servers.end(), e.flow))  != servers.end()) // >= 降档
					nodes[e.to].bestCdnId = it - servers.begin(); // 存放下标，nodes输出路径的时候用
				else nodes[e.to].bestCdnId = servers.size() - 1; // 最大的level

				if(downShift) {
					eId[e.to] = G[superSource][i];
					e.cap = servers[nodes[e.to].bestCdnId].outFlow;
					cdnFlow += servers[nodes[e.to].bestCdnId].outFlow; // 计算总费用
				}
				cdnCost += servers[nodes[e.to].bestCdnId].cost; // 计算总费用
				// printf("%d e.flow: %d/%d(%d)\n", e.to, e.flow, servers[nodes[e.to].bestCdnId].outFlow, servers[nodes[e.to].bestCdnId].level);
			}
			minCdnFlowCost = min(minCdnFlowCost, cost + cdnCost); // 更新总费用
			// printf("minCdnFlowCost %d\n", minCdnFlowCost);

			if(downShift) {
				bool exit = false;

				while(! exit) {
					diff.clear();
					for (size_t i = 0; i < G[superSource].size(); i++) {
						Edge &e = edges[G[superSource][i]];
						if(nodes[e.to].bestCdnId  == 0) diff.push_back(make_pair(0, e.to));
						else diff.push_back(make_pair(
									(servers[nodes[e.to].bestCdnId].cost - servers[nodes[e.to].bestCdnId - 1].cost) /
									(e.flow - servers[nodes[e.to].bestCdnId - 1].outFlow)
									, e.to));
					}
					sort(diff.begin(), diff.end(), greater<pair<int, int>>());

					for(size_t i = 0; i < diff.size(); ++i) {
						int u = diff[i].second;
						if(nodes[u].bestCdnId == 0) continue;
						for(size_t j = 0; j < edges.size(); ++j) {
							edges[j].cost = edges[j].oldCost;
							edges[j].flow = 0; // 重置流量
						}

						int dc = servers[nodes[u].bestCdnId].cost - servers[nodes[u].bestCdnId - 1].cost;
						int df = servers[nodes[u].bestCdnId].outFlow - servers[nodes[u].bestCdnId - 1].outFlow;
						if(cdnFlow - df >= needFlow) { // 预判
							cdnCost -= dc; // 更新Cdn费用
							cdnFlow -= df;
							edges[eId[u]].cap = servers[nodes[u].bestCdnId - 1].outFlow;
							--nodes[u].bestCdnId;
							cost = pathFlowCost();

							if(cost == -1 || cost + cdnCost > minCdnFlowCost) { // 降档失败
								cdnCost += dc;
								cdnFlow += df;
								edges[eId[u]].cap = servers[nodes[u].bestCdnId + 1].outFlow;
								++nodes[u].bestCdnId;
								if(cost == -1){
									exit = true;
									break;
								}
							} else {// 降档成功，有解
								minCdnFlowCost = cost + cdnCost;
								break;
								// printf("success: %d\n", u);
							}
						}  else {
							exit = true;
							break;
						}
					}

					// printf("minCdnFlowCost %d\n", minCdnFlowCost);
					// printf("%d diff: %d\n", u, diff[i].first);
					// printf("cdnFlow: %d cdnCost: %d needFlow: %d minCdnFlowCost: %d\n", cdnFlow, cdnCost, needFlow, minCdnFlowCost);
				}
			}


			cost = minCdnFlowCost;


			// 计算部署费用
#ifdef _DEBUG
			int realCost = cost;
#endif
			for(auto c: cdn) {
				cost += nodes[c].deployCost;
#ifdef _DEBUG
				realCost += nodes[c].deployCost;
#endif
			}

#ifdef _DEBUG
			realMinCost = min(realMinCost, realCost);
#endif

			if(cost < solutionPath.first) {
				if(downShift) {
					for(size_t j = 0; j < edges.size(); ++j) {
						edges[j].cost = edges[j].oldCost;
						edges[j].flow = 0; // 重置流量
					}
					pathFlowCost();
				}

				// 打印档次
				/*
				vector<pair<int,int>> v;
				for (size_t i = 0; i < G[superSource].size(); i++) { // 降档
					Edge &e = edges[G[superSource][i]];
					v.push_back(make_pair(e.to, G[superSource][i]));
				}
				sort(v.begin(), v.end());

				for(size_t i = 0; i < v.size(); ++i) {
					Edge &e =  edges[v[i].second];
					// printf("%d e.flow: %d/%d(%d)\n", e.to, e.flow, servers[nodes[e.to].bestCdnId].outFlow, servers[nodes[e.to].bestCdnId].level);
					printf("%d\t%d\n", e.to, servers[nodes[e.to].bestCdnId].level);
				}
				*/
				getPath(cost); // 更新方案
			}
			return cost;
		}


		inline void setCdn(const unordered_set<int> & cdn) {
			reset();
			for(int x: cdn)
				AddEdge(superSource, x, maxFlowServer.outFlow, 0);
		}
	public:
		struct Node{
			int deployCost; // 节点部署费用
			int nodeFlow; // 每个节点的流量
			int bestCdnId; // 存放每个节点最适合的服务器档次（下标）
			double evaluation; // 每个节点的评估值
			Node() {
				deployCost = bestCdnId = nodeFlow = evaluation = 0;
			}
		} nodes[10000 + 5];

		vector<int> G[N]; // 图
		vector<Edge> edges; // 边集
		int networkNum, edgeNum, consumerNum, needFlow, costPerCDN = 0;
		static const int INF = 0x3f3f3f3f;
		friend class MCMF_SCALING;

		void inline showRealMinCost() {
#ifdef _DEBUG
			printf("\x1B[31mReal minCost: %d/%d\x1B[0m\n", realMinCost, consumerNum * costPerCDN);
#endif
		}

		inline bool isConsumer(int u) {
			return u >= networkNum && u < superSource;
		}


		MCMF() {
			needFlow = 0;
		};
		inline void setCostCdnGap(int x) {
			minCostCdnGap = x;
		}
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
