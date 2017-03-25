/*************************************************************************
  > File Name: spfa.cpp
  > Author: Netcan
  > Blog: http://www.netcan666.com
  > Mail: 1469709759@qq.com
  > Created Time: 2017-03-21 Tue 21:25:20 CST
 ************************************************************************/

#include "spfa.h"

char MCMF::topo[50000*1000*6];

void MCMF::getPath(int cost) {
	if(cost < solutionPath.first) {
		solutionPath.second.clear(); // 记得清理
		solutionPath.first = cost;
		vector<int> tmpPath;
		bzero(vis, sizeof(vis));
		solutionPath.second.clear(); // 清空可行解
		findPath(tmpPath, superSource, INF, INF);
	}
}

int MCMF::findPath(vector<int> & tmpPath, int u, int minFlow, int totalFlow) { // dfs，深搜路径，路径上的最小流量，总流量
	if(vis[u]) return 0;
	else if(u >= networkNum && u < superSource) { // 到达消费节点，找到一条路径
		solutionPath.second.push_back(tmpPath);
		solutionPath.second.back().push_back(u - networkNum); // 转换为消费节点的id
		solutionPath.second.back().push_back(minFlow);
		return minFlow;
	}

	vis[u] = true;
	if(u != superSource) tmpPath.push_back(u);

	int tf = totalFlow;
	for(size_t i = 0; i < G[u].size(); ++i) {
		Edge &e = edges[G[u][i]];
		// printf("%d->%d flow: %d\n", e.from, e.to, e.flow);
		if(e.flow > 0) { // 流过的流量>0
			int v = e.to;
			if(!vis [v]) {
				if(totalFlow > 0) {
					int t = findPath(tmpPath, v,
							min(minFlow, min(totalFlow, e.flow)),
							min(totalFlow, e.flow));
					e.flow -= t;
					totalFlow -= t;
				}
				else break;
			}
		}
	}

	vis[u] = false;
	tmpPath.pop_back();
	return tf;
}

bool MCMF::BellmanFord(int s, int t, int &flow, int &cost) {
	memset(d, 0x3f, sizeof(d));
	memset(vis, 0, sizeof(vis));
	vis[s] = 1; d[s] = 0; f[s] = MCMF::INF; p[s] = 0;

	queue<int> Q; Q.push(s);

	while (!Q.empty()) {
		int u = Q.front(); Q.pop();
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

	if (d[t] == INF) return false;

	flow += f[t];
	cost += d[t] * f[t];

	int u = t;
	while (u != s) {
		edges[p[u]].flow += f[t];
		edges[p[u] ^ 1].flow -= f[t];
		u = edges[p[u]].from;
	}
	return true;
}

void MCMF::AddEdge(int from, int to, int cap, int cost) {
	edges.push_back(Edge(from, to, cap, 0, cost));
	edges.push_back(Edge(to, from, 0, 0, -cost));

	int m = edges.size();
	G[from].push_back(m - 2);
	G[to].push_back(m - 1);

	if(from < networkNum && to >= networkNum)  // 网络节点直连消费节点，计算需要的总共流量
		needFlow += cap;
}


void MCMF::showSolution() const{
	int totalFlow = 0;
	for(const auto &x : solutionPath.second) {
		for(vector<int>::const_iterator i = x.begin(); i != x.end() - 1; ++i) {
			if(i == x.begin()) printf("%d", *i);
			else printf("->%d", *i);
		}
		totalFlow += x.back();
		printf(" flow: %d\n", x.back());
	}
	printf("Flow :%d/%d Cost: %d/%d\n", totalFlow, needFlow, solutionPath.first, costPerCDN * consumerNum);
}

void MCMF::loadGraph() {
	scanf("%d%d%d", &networkNum, &edgeNum, &consumerNum);
	scanf("%d", &costPerCDN);

	solutionPath.first = consumerNum * costPerCDN;

	printf("networkNode: %d edgeNum: %d cosumerNode: %d\n",
			networkNum, edgeNum, consumerNum);
	printf("costPerCDN: %d\n", costPerCDN);

	superSource = consumerNum + networkNum; // 超级源点、汇点
	superSink = consumerNum + networkNum + 1;

	int from, to, bandwidth, cpb;

	for(int i=0; i < edgeNum; ++i) {
		scanf("%d%d%d%d", &from, &to, &bandwidth, &cpb);
		AddEdge(from, to, bandwidth, cpb);
		AddEdge(to, from, bandwidth, cpb);
	}

	for(int i=0; i < consumerNum; ++i) {
		scanf("%d%d%d", &from, &to, &bandwidth);
		AddEdge(to, from + networkNum, bandwidth, 0);
		AddEdge(from + networkNum, superSink, bandwidth, 0);

		vector<int> path{to, from, bandwidth}; // 直连策略
		solutionPath.second.push_back(move(path));
	}

	oldEdges = edges;
}

void MCMF::loadGraph(char * topo[MAX_EDGE_NUM], int line_num) {
	sscanf(topo[0], "%d%d%d", &networkNum, &edgeNum, &consumerNum);
	sscanf(topo[2], "%d", &costPerCDN);

	solutionPath.first = consumerNum * costPerCDN; // asdf


	superSource = consumerNum + networkNum; // 超级源点、汇点
	superSink = consumerNum + networkNum + 1;

	int from, to, bandwidth, cpb;

	int i;
	for(i=0; i < edgeNum; ++i) {
		sscanf(topo[i+4], "%d%d%d%d", &from, &to, &bandwidth, &cpb);
		AddEdge(from, to, bandwidth, cpb);
		AddEdge(to, from, bandwidth, cpb);
	}

	i += 4;
	for(++i; i < line_num; ++i) {
		sscanf(topo[i], "%d%d%d", &from, &to, &bandwidth);
		AddEdge(to, from + networkNum, bandwidth, 0);
		AddEdge(from + networkNum, superSink, bandwidth, 0);

		vector<int> path{to, from, bandwidth}; // 直连策略
		solutionPath.second.push_back(move(path));
	}

	oldEdges = edges;
}

const char* MCMF::outputPath() {
	char buffer[10];
	char *pt = topo, *pb = buffer;
	snprintf(buffer, sizeof(buffer), "%ld\n\n", solutionPath.second.size());
	while(*pb && (*pt++ = *pb++));
	for(auto &x: solutionPath.second) {
		for(auto it = x.begin(); it != x.end(); ++it) {
			snprintf(buffer, sizeof(buffer), it == x.begin() ? "%d":" %d", *it);
			pb = buffer;
			while(*pb && (*pt++ = *pb++));
		}
		*pt++ = '\n';
	}
	*--pt = 0;
	return topo;
}

MCMF mcmf;

// int main() {
	// freopen("case_example/case4.txt", "r", stdin);
	// unordered_set<int> cdn = {
		// 12, 15, 20, 22, 26, 37, 48
	// };
	// // load graph
	// // set CDN and build super Source
	// mcmf.loadGraph();

	// puts("======Output=======");
	// mcmf.setCdn(cdn);
	// mcmf.minCost();
	// mcmf.showPath();

	// return 0;
// }
