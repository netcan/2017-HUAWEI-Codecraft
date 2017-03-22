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
	pid_t pid;
	int stat;
	pid = wait(&stat);
	write_result(mcmf.outputPath(), fileName);
	runing = false;
	return;
}

void deploy_server(char * topo[MAX_EDGE_NUM], int line_num,char * filename)
{
	fileName = filename;
	srand(time(0));
	Signal(SIGALRM, timeOutHandler);
	alarm(45);
	// 启动计时器

	// 开始计算
	mcmf.loadGraph(topo, line_num);
	unordered_set<int> cdn;
	// load graph
	// mcmf.setCdn(cdn);
	// mcmf.minCost();
	// mcmf.showPath();

	size_t cnt = 1;
	int minCost = MCMF::INF, cost = MCMF::INF;
	bool find = false;

	while(runing) {
		while(cdn.size() < cnt)
			cdn.insert(rand() % mcmf.networkNum);
		mcmf.setCdn(cdn);
		cost = mcmf.minCost();
		if(cost != -1) {
			minCost = min(minCost, cost);
			find = true;
		}
		else if(!find) cnt = min<int>(cnt+1, mcmf.consumerNum);

		cdn.clear();
	}
	printf("minCost: %d\n", minCost);

}
