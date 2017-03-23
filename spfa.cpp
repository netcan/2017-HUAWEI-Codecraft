/*************************************************************************
  > File Name: spfa.cpp
  > Author: Netcan
  > Blog: http://www.netcan666.com
  > Mail: 1469709759@qq.com
  > Created Time: 2017-03-21 Tue 21:25:20 CST
 ************************************************************************/

#include "spfa.h"

char MCMF::topo[50000*1000*6];

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
	vector<int> argumentPath;
	argumentPath.push_back(f[t]);
	while (u != s) {
		edges[p[u]].flow += f[t];
		edges[p[u] ^ 1].flow -= f[t];
		u = edges[p[u]].from;
		if(u != s) argumentPath.push_back(u < networkNum? u: u - networkNum);
	}
	path.second.push_back(move(argumentPath));

	return true;
}

void MCMF::AddEdge(int from, int to, int cap, int cost) {
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

	if(from < networkNum && to >= networkNum)  // 网络节点直连消费节点，计算需要的总共流量
		needFlow += cap;
}

void MCMF::showPath() const {
	if(path.second.empty()) return;
	int totalFlow = 0;
	for(auto &x : path.second) {
		for(vector<int>::const_reverse_iterator i = x.rbegin(); i != x.rend() - 1; ++i) {
			if(i == x.rbegin()) printf("%d", *i);
			else printf("->%d", *i);
		}
		totalFlow += x[0];
		printf(" flow: %d\n", x[0]);
	}
	printf("Flow :%d/%d Cost: %d\n", totalFlow, needFlow, path.first.first);
}

void MCMF::showSolution() const{
	int totalFlow = 0;
	for(const auto &x : solutionPath.second) {
		for(vector<int>::const_reverse_iterator i = x.rbegin(); i != x.rend() - 1; ++i) {
			if(i == x.rbegin()) printf("%d", *i);
			else printf("->%d", *i);
		}
		totalFlow += x[0];
		printf(" flow: %d\n", x[0]);
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

		vector<int> path{bandwidth, from, to};
		solutionPath.second.push_back(move(path));
	}

}

void MCMF::loadGraph(char * topo[MAX_EDGE_NUM], int line_num) {
	sscanf(topo[0], "%d%d%d", &networkNum, &edgeNum, &consumerNum);
	sscanf(topo[2], "%d", &costPerCDN);

	solutionPath.first = consumerNum * costPerCDN;


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

		vector<int> path{bandwidth, from, to};
		solutionPath.second.push_back(move(path));
	}

}

const char* MCMF::outputPath() {
	char buffer[10];
	char *pt = topo, *pb = buffer;
	snprintf(buffer, sizeof(buffer), "%ld\n\n", solutionPath.second.size());
	while(*pb && (*pt++ = *pb++));
	for(auto &x: solutionPath.second) {
		for(auto it = x.rbegin(); it != x.rend(); ++it) {
			snprintf(buffer, sizeof(buffer), it == x.rbegin() ? "%d":" %d", *it);
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
