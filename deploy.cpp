#include "deploy.h"
#include <stdio.h>

typedef void Sigfunc(int);
char *fileName;
bool runing = true;

Sigfunc *
Signal(int signo, Sigfunc *func) {
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
	write_result(mcmf.outputPath(), fileName);
	runing = false;
	return;
}

void deploy_server(char * topo[MAX_EDGE_NUM], int line_num,char * filename)
{
	fileName = filename;
	srand(time(0));
	Signal(SIGALRM, timeOutHandler);
	alarm(85);
	// 启动计时器

	// 开始计算
	mcmf.loadGraph(topo, line_num);
	unordered_set<int> backup, cur;
	// load graph
	// mcmf.setCdn(cdn);
	// mcmf.minCost();
	// mcmf.showPath();

	for(int u=0; u < mcmf.consumerNum; ++u)
		backup.insert(mcmf.edges[mcmf.G[u + mcmf.networkNum][0]].to);

	int minCost = MCMF::INF, backCost = MCMF::INF, curCost = MCMF::INF;
	mcmf.setCdn(backup);
	backCost = mcmf.minCost();
	minCost = min(minCost, backCost);

	double T = 20.0;
	// for(auto x: backup)
		// printf("%d\n", x);

	while(T > 10e-5 && runing) {

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
		// puts("=====BACKUP=====");
		// for(auto x: backup)
			// printf("%d ", x);
		// printf("\nu->v: %d->%d\n", u, v);
		// puts("=====CURRENT====");
		// for(auto x: cur)
			// printf("%d ", x);
		// puts("");

		mcmf.setCdn(cur);
		curCost = mcmf.minCost();

		if(curCost == -1)  {// 无解
			// puts("ERROR\n");
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

		T *= 0.9999;
		// printf("T=%lf\n", T);
		// puts("");
	}

	printf("minCost: %d/%d cdnNum: %ld\n", minCost, mcmf.consumerNum * mcmf.costPerCDN, backup.size());
}
