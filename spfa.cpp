/*************************************************************************
	> File Name: spfa.cpp
	  > Author: Netcan
	  > Blog: http://www.netcan666.com
	  > Mail: 1469709759@qq.com
	  > Created Time: 2017-03-21 Tue 21:25:20 CST
 ************************************************************************/

#include <cstdio>
#include <cstring>
#include <algorithm>
#include <queue>
#include <vector>
using namespace std;


class MCMF{
	private:
		struct Edge{
			int from, to, cap, flow ,cost;
			Edge() {}
			Edge(int from, int to, int cap, int flow, int cost):from(from), to(to), cap(cap), flow(flow), cost(cost) {}
		};
		static const int N = 1500+5;
		int n, superSource, superSink; // 点数，超级源点/汇点，需要的流量
		int networkNode, edgeNum, consumerNode, costPerCDN, needFlow;
		vector<Edge> edges, oldEdges; // 边集，边集备份
		vector<int> G[N]; // 图
		int d[N], f[N], p[N]; // 最小费用，当前流量，父节点（增广路径），中间变量
		bool vis[N]; // 标记指针
		pair<pair<int, int>, vector<pair<vector<int>, int>>> path; // 费用/流量，路径/路径上的流量


		bool BellmanFord(int s, int t, int &flow, int &cost) {
			memset(d, 0x3f, sizeof(d));
			memset(vis, 0, sizeof(vis));
			vis[s] = 1; d[s] = 0; f[s] = MCMF::INF; p[s] = 0;

			queue<int> Q;
			Q.push(s);

			while (!Q.empty()) {
				int u = Q.front();
				Q.pop();
				vis[u] = 0;

				for (size_t i = 0; i < G[u].size(); i++) {
					const Edge &e = edges[G[u][i]];
					if (e.cap > e.flow && d[e.to] > d[u] + e.cost) {
						d[e.to] = d[u] + e.cost;
						p[e.to] = G[u][i];
						f[e.to] = min(f[u], e.cap - e.flow);
						if (!vis[e.to]) {
							vis[e.to] = true;
							Q.push(e.to);
						}
					}
				}
			}

			if (d[t] == INF)
				return false;

			flow += f[t];
			cost += d[t] * f[t];

			int u = t;
			vector<int> argumentPath;
			while (u != s) {
				edges[p[u]].flow += f[t];
				edges[p[u] ^ 1].flow -= f[t];
				u = edges[p[u]].from;
				if(u != s) argumentPath.push_back(u);
			}
			path.second.push_back(make_pair(move(argumentPath), f[t]));

			return true;
		}
		void reset() { // 还原初始状态，删除源点
			G[superSource].clear();
			path.second.clear();
			if(! oldEdges.empty())
				edges = move(oldEdges);
		}
	public:
		static const int INF = 0x3f3f3f3f;
		MCMF(int n, int superSource, int superSink, int networkNode, int edgeNum, int consumerNode ,int costPerCDN, int needFlow = 0):
			n(n), superSource(superSource), superSink(superSink), networkNode(networkNode), edgeNum(edgeNum),
			consumerNode(consumerNode), costPerCDN(costPerCDN), needFlow(needFlow) {}

		void AddEdge(int from, int to, int cap, int cost) {
			edges.push_back(Edge(from, to, cap, 0, cost));
			if(to != superSink && from != superSource)
				edges.push_back(Edge(to, from, 0, 0, -cost));

			int m = edges.size();
			if(from == superSource || to == superSink)
				G[from].push_back(m - 1);
			else {
				G[from].push_back(m - 2);
				G[to].push_back(m - 1);
			}

			if(from < networkNode && to >= networkNode) // 网络节点直连消费节点，计算需要的总共流量
				needFlow += cap;
		}

		void showPath() {
			int totalFlow = 0;
			for(auto &x : path.second) {
				for(vector<int>::const_reverse_iterator i = x.first.rbegin(); i != x.first.rend(); ++i) {
					if(i == x.first.rbegin()) printf("%d", *i);
					else printf("->%d", *i);
				}
				totalFlow += x.second;
				printf(" flow: %d\n", x.second);
			}
			printf("Flow :%d/%d Cost: %d\n", totalFlow, needFlow, path.first.first);
		}

		int minCost() {
			if(oldEdges.empty()) oldEdges = edges; // 备份边集
			else return path.first.first;

			int flow = 0, cost = 0;
			while (BellmanFord(superSource, superSink, flow, cost));
			path.first.first = cost + G[superSource].size() * costPerCDN; // 加上服务器费用
			path.first.second = flow;
			return flow < needFlow ?-1: cost;
		}
		void setCdn(const vector <int> & cdn) {
			reset();
			for(size_t i=0; i<cdn.size(); ++i)
				AddEdge(superSource, cdn[i], MCMF::INF, 0);
		}
};

int main() {
	freopen("case_example/case0.txt", "r", stdin);
	vector<int> cdn = {7, 13, 15, 22, 37, 38, 43};
	// load graph
	int networkNode, edgeNum, consumerNode, costPerCDN;
	scanf("%d%d%d", &networkNode, &edgeNum, &consumerNode);
	scanf("%d", &costPerCDN);

	printf("networkNode: %d edgeNum: %d cosumerNode: %d\n",
			networkNode, edgeNum, consumerNode);
	printf("costPerCDN: %d\n", costPerCDN);

	int superSource = consumerNode + networkNode; // 超级源点、汇点
	int superSink = consumerNode + networkNode + 1;

	MCMF mcmf(superSink + 1, superSource, superSink,
			networkNode, edgeNum, consumerNode, costPerCDN);

	int from, to, bandwidth, cpb;

	for(int i=0; i < edgeNum; ++i) {
		scanf("%d%d%d%d", &from, &to, &bandwidth, &cpb);
		mcmf.AddEdge(from, to, bandwidth, cpb);
		mcmf.AddEdge(to, from, bandwidth, cpb);
		// printf("from: %d to: %d bandwidth: %d cpb: %d\n",
				// from, to, bandwidth, cpb);
	}

	for(int i=0; i < consumerNode; ++i) {
		scanf("%d%d%d", &from, &to, &bandwidth);
		// mcmf.AddEdge(from + networkNode, to, bandwidth, 0);
		mcmf.AddEdge(to, from + networkNode, bandwidth, 0);
		mcmf.AddEdge(from + networkNode, superSink, bandwidth, 0);
		// printf("consumer: from: %d to: %d bandwidth: %d\n",
				// from, to, bandwidth);
	}

	// set CDN and build super Source
	puts("======Output=======");
	mcmf.setCdn(cdn);
	mcmf.minCost();
	mcmf.showPath();



    return 0;
}
