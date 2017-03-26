#include "deploy.h"
#include <stdio.h>

typedef void (sigFunc)(int);
bool runing = true;

sigFunc *
Signal(int signo, sigFunc *func) {
	struct sigaction	act, oact;
	act.sa_handler = func;
	sigemptyset(&act.sa_mask);
	act.sa_flags = 0;
	if (sigaction(signo, &act, &oact) < 0)
		return(SIG_ERR);
	return(oact.sa_handler);
}
/* end signal */

void timeOutHandler(int signo) {
	runing = false;
	return;
}

void SA(unordered_set<int>init = {}, double T = 20.0, double delta = 0.99999) { // 模拟退火，初始温度，迭代系数
	// double T = 20.0, delta = 0.99999; // 初始温度20, 0.999-0.999999

	unordered_set<int> backup, cur;

	if(init.empty()) {
		for(int u=0; u < mcmf.consumerNum; ++u)  // 初始位置
			backup.insert(mcmf.edges[mcmf.G[u + mcmf.networkNum][0]].to);
	} else backup = move(init);

	int minCost = MCMF::INF, backCost = MCMF::INF, curCost = MCMF::INF;
	backCost = mcmf.minCost_Set(backup);
	minCost = min(minCost, backCost);

	int iterationCnt = 0;
	while(T > 0.1 && runing) {
		int u = -1;
		do {
			for(auto x: backup) {
				if(rand() < RAND_MAX * 1.0 / mcmf.networkNum) {
					u = x;
					break;
				}
			}
		} while(u == -1);

		int selectEdge = 0, v; // (u, v)随机选点

		do {
			for(selectEdge = 0; (selectEdge < (int)mcmf.G[u].size() - 1) &&
					rand() > RAND_MAX * 1.0 / mcmf.G[u].size(); ++selectEdge);
		}
		while( (v = mcmf.edges[mcmf.G[u][selectEdge]].to) >= mcmf.networkNum);

		for(int x: backup) {
			if(x == u) cur.insert(v);
			else cur.insert(x);
		}

		curCost = mcmf.minCost_Set(cur);
		++iterationCnt;

		if(curCost == -1)  {// 无解
			cur.clear();
		}
		else {
			int dC = curCost - backCost;
			// printf("dC: %d\n", dC);
			if(dC < 0 || exp(-dC / T) * RAND_MAX > rand())  {// 接受
				backup = move(cur);
				backCost = curCost;
			} else {
				cur.clear();
			}

			minCost = min(minCost, backCost);
		}
		T *= delta;
	}

	printf("T=%lf iterationCnt=%d\n", T, iterationCnt);
	// printf("Deploy CDN(%ld):\n", backup.size());
	// for(int x: backup)
		// printf("%d ", x);
	// puts("\n=====Solution======");
	// mcmf.showSolution();
	printf("minCost: %d/%d cdnNum: %ld\n\n", minCost, mcmf.consumerNum * mcmf.costPerCDN, backup.size());
}


void Tabu(unordered_set<int>init = {}) { // 禁忌搜索
	typedef unordered_set<int> X;
	list<int> H; // 禁忌表，队列

	X x_best;
	int minCost;
	X x_now;
	if(init.empty()) {
		for(int u=0; u < mcmf.consumerNum; ++u)  // 初始位置
			x_now.insert(mcmf.edges[mcmf.G[u + mcmf.networkNum][0]].to);
	} else x_now = move(init);

	pair<int, X> x_next{MCMF::INF, {}}; // 转移
	H.push_back(minCost = mcmf.minCost_Set(x_now));

	// for(int x: x_now)
		// printf("%d ", x);
	// puts("");

	int iterationCnt = 0;
	while(runing) {
		int Len = 0;
		for(int u: x_now) {
			for(size_t i = 0; i < mcmf.G[u].size(); i+=2) {
				++Len;
				int v = mcmf.edges[mcmf.G[u][i]].to; // u->v
				if(v < mcmf.networkNum) {
					X tmp{}; // 邻居
					for(int uu: x_now) {
						if(uu != u) tmp.insert(uu);
						else tmp.insert(v);
					}

					// puts("==============");
					// for(int x: tmp)
						// printf("%d ", x);
					// puts("\n==============");

					int cost = mcmf.minCost_Set(tmp);
					if(find(H.begin(), H.end(), cost) == H.end() && cost < x_next.first) {
						x_next.first = cost;
						x_next.second = move(tmp);
						if(minCost > cost) {
							minCost = cost;
							x_best = x_next.second;
						}
					}
				}
			}
		}
		H.push_back(x_next.first); // 入队
		x_next.first = MCMF::INF;
		x_now = move(x_next.second);
		++iterationCnt;
		while(H.size() > sqrt(Len)) H.pop_front();
	}

	mcmf.showSolution();
	printf("iterationCnt = %d\n", iterationCnt);
}


void deploy_server(char * topo[MAX_EDGE_NUM], int line_num,char * filename)
{
	srand(time(0));
	Signal(SIGALRM, timeOutHandler);
	// 启动计时器
	alarm(88);
	mcmf.loadGraph(topo, line_num);
	// Tabu();
	SA();

	//- test
	/*
	double T = 1.0, delta = 0.99999;
	double bestT = T, bestDelta = delta;
	int minCost = MCMF::INF;
	int cost = 0;
	for(; T <= 100.0; T+=1) {
		alarm(88);
		if( (cost = SA(T,delta)) < minCost) {
			minCost = cost;
			bestT = T;
			bestDelta = delta;
		}
		printf("bestT = %lf bestDelta = %lf minCost = %d\n", bestT, bestDelta, minCost);
		runing = true;
	}
	*/
	//- test End

	// 开始计算
	write_result(mcmf.outputPath(), filename);

}
