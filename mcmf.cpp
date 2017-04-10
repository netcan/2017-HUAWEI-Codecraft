/*************************************************************************
  > File Name: spfa.cpp
  > Author: Netcan
  > Blog: http://www.netcan666.com
  > Mail: 1469709759@qq.com
  > Created Time: 2017-03-21 Tue 21:25:20 CST
 ************************************************************************/

#include "mcmf.h"

char MCMF::topo[50000*1000*6];

void MCMF::getPath(int cost) {
	if(cost != -1 && cost < solutionPath.first) {
		solutionPath.first = cost;
		solutionPath.second.clear(); // 记得清理
		vector<int> tmpPath;
		bzero(vis, sizeof(vis));
		findPath(tmpPath, superSource, INF, INF);
	}

}

int MCMF::findPath(vector<int> & tmpPath, int u, int minFlow, int totalFlow) { // dfs，深搜路径，路径上的最小流量，总流量
	if(vis[u]) return 0;
	else if(isConsumer(u)) { // 到达消费节点，找到一条路径
		solutionPath.second.push_back(tmpPath);
		solutionPath.second.back().push_back(u - networkNum); // 转换为消费节点的id
		solutionPath.second.back().push_back(minFlow);
		solutionPath.second.back().push_back(maxFlowServer.level); // 档次
		return minFlow;
	}

	vis[u] = true;
	if(u < superSource) tmpPath.push_back(u);

	int tf = totalFlow;
	for(size_t i = 0; i < G[u].size(); ++i) {
		Edge &e = edges[G[u][i]];
		if(e.flow > 0) { // 流过的流量>0
			// printf("%d->%d flow: %d\n", e.from, e.to, e.flow);
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
	if(u < superSource) tmpPath.pop_back();
	return tf;
}

bool MCMF::BellmanFord(int s, int t, int &flow, int &cost) {
	memset(d, 0x3f, sizeof(d));
	bzero(vis, sizeof(vis));
	vis[s] = 1; d[s] = 0; f[s] = MCMF::INF; p[s] = 0;

	queue.reset();
	queue.push(s);

	while (!queue.empty()) {
		int u = queue.front(); queue.pop();
		vis[u] = 0;

		for (size_t i = 0; i < G[u].size(); i++) {
			const Edge &e = edges[G[u][i]];
			if (e.cap > e.flow && d[e.to] > d[u] + e.cost) {
				d[e.to] = d[u] + e.cost;
				p[e.to] = G[u][i];
				f[e.to] = min(f[u], e.cap - e.flow);
				if (!vis[e.to]) {
					vis[e.to] = true;
					queue.push(e.to);
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

int MCMF::aug(int u, int minFlow, int &tmpCost, int &cost) {
	if(u == superSink) { // 到达终点
		cost += tmpCost * minFlow;
		return minFlow;
	}
	vis[u] = true;
	int tf = minFlow;
	for(size_t i = 0; i < G[u].size(); ++i) {
		Edge &e = edges[G[u][i]];
		if( (e.cap - e.flow) && !e.cost && !vis[e.to]) {
			int d = aug(e.to, min(tf, (e.cap - e.flow)), tmpCost, cost);
			e.flow += d;
			edges[G[u][i] ^ 1].flow -= d;
			tf -= d;
			if(! tf) return minFlow;
		}
	}
	return minFlow - tf;
}

bool MCMF::modLabel(int &tmpCost) {
	int d = INF;
	for(int u=0; u< Vn + networkNum; ++u) // 遍历完全部节点
		if(vis[u]) {
			for(size_t i = 0; i < G[u].size(); ++i) {
				Edge &e = edges[G[u][i]];
				if( (e.cap - e.flow) && !vis[e.to] && e.cost < d)
					d = e.cost;
			}
		}
	if(d == INF) return false;

	for(int u=0; u<Vn + networkNum; ++u)
		if(vis[u]) {
			for(size_t i = 0; i < G[u].size(); ++i) {
				edges[G[u][i]].cost -= d;
				edges[G[u][i] ^ 1].cost += d;
			}
		}
	tmpCost += d;
	return true;

	// SLF优化
	/*
	memset(d, 0x3f, sizeof(d));
	d[superSink] = 0;
	static deque<int> que; que.push_back(superSink);
	while(que.size())
	{
		int dt, u = que.front(); que.pop_front();
		for(size_t i = 0; i < G[u].size(); ++i) {
			Edge &e = edges[G[u][i]], &re = edges[G[u][i] ^ 1];
			if( (re.cap - re.flow) && (dt = d[u] - e.cost) < d[e.to] )
				(d[e.to] = dt) <= d[que.size() ? que.front() : 0]
					? que.push_front(e.to) : que.push_back(e.to);
		}
	}
	for(int u=0; u<=superSink; ++u)
		for(size_t i = 0; i < G[u].size(); ++i) {
			Edge &e = edges[G[u][i]];
			e.cost += d[e.to] - d[u];
		}

	tmpCost += d[superSource];
	return d[superSource] < INF;
	*/
}



void MCMF::AddEdge(int from, int to, int cap, int cost) {
	edges.push_back(Edge(from, to, cap, 0, cost));
	edges.push_back(Edge(to, from, 0, 0, -cost));

	int m = edges.size();
	G[from].push_back(m - 2);
	G[to].push_back(m - 1);

	if(from < networkNum && isConsumer(to))  // 网络节点直连消费节点，计算需要的总共流量
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

void MCMF::loadGraph(char * topo[MAX_EDGE_NUM], int line_num) {
	sscanf(topo[0], "%d%d%d", &networkNum, &edgeNum, &consumerNum); // 网络节点数量 网络链路数量 消费节点数量

	// solutionPath.first = consumerNum * costPerCDN; // 待求

	superSource = consumerNum + networkNum; // 超级源点、汇点
	superSink = consumerNum + networkNum + 1;
	Vn = superSink + 1;

	int a, b, c, d;
	int i;
	for(i = 2; i < line_num && !isspace(topo[i][0]); ++i) {
		sscanf(topo[i], "%d%d%d", &a, &b, &c); // 服务器硬件档次ID 输出能力 硬件成本
		servers.push_back(Server(a, b, c));
		if(b > maxFlowServer.outFlow) maxFlowServer = servers.back();
		// printf("level: %d outFlow: %d cost: %d\n", a, b, c);
	}
	// printf("maxFlowServer level: %d outFlow: %d cost: %d\n", maxFlowServer.level, maxFlowServer.outFlow, maxFlowServer.cost);

	for(++i; i < line_num && !isspace(topo[i][0]); ++i) {
		sscanf(topo[i], "%d%d", &a, &b); // 网络节点ID 部署成本
		deployCost[a] = b;
		// printf("node: %d cost: %d\n", a, b);
	}

	for(++i; i < line_num && !isspace(topo[i][0]); ++i) {
		sscanf(topo[i], "%d%d%d%d", &a, &b, &c, &d); // 链路起始节点ID 链路终止节点ID 总带宽大小 单位网络租用费
		AddEdge(a, b, c, d);
		AddEdge(b, a, c, d);
		// printf("u: %d v: %d bandwidth: %d cost: %d\n", a, b, c, d);
	}

	for(++i; i < line_num; ++i) {
		sscanf(topo[i], "%d%d%d", &a, &b, &c); // 消费节点ID 相连网络节点ID 视频带宽消耗需求
		AddEdge(b, a + networkNum, c, 0); // 与网络节点相连
		AddEdge(a + networkNum, superSink, c, 0); // 与汇点相连
		// printf("consumer: %d connect: %d need: %d\n", a, b, c);

		// vector<int> path{to, from, bandwidth}; // 直连策略
		// solutionPath.second.push_back(move(path));
	}
	edgeNum = edges.size(); // 边数
	costPerCDN = maxFlowServer.cost; // 以最大档次的费用为准
	solutionPath.first = INF;
}

const char* MCMF::outputPath() {
	// getPath(solutionPath.first, true); // 放到最后才遍历路径，提高性能
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
