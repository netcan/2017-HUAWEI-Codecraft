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

void SA() {
	unordered_set<int> backup, cur;

	for(int u=0; u < mcmf.consumerNum; ++u)
		backup.insert(mcmf.edges[mcmf.G[u + mcmf.networkNum][0]].to);

	int minCost = MCMF::INF, backCost = MCMF::INF, curCost = MCMF::INF;
	mcmf.setCdn(backup);
	backCost = mcmf.minCost();
	minCost = min(minCost, backCost);

	double T = 20.0, delta = 0.9999; // 初始温度

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

		mcmf.setCdn(cur);
		curCost = mcmf.minCost();

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

	printf("Deploy CDN:\n");
	for(int x: backup)
		printf("%d ", x);
	puts("");
	printf("minCost: %d/%d cdnNum: %ld\n", minCost, mcmf.consumerNum * mcmf.costPerCDN, backup.size());
}

void deploy_server(char * topo[MAX_EDGE_NUM], int line_num,char * filename)
{
	srand(time(0));
	Signal(SIGALRM, timeOutHandler);
	alarm(85);
	// 启动计时器
	mcmf.loadGraph(topo, line_num);
	SA();

	// 开始计算
	write_result(mcmf.outputPath(), filename);

}
