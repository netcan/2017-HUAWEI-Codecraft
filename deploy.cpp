#include "deploy.h"
#include <stdio.h>
#include "random.h"
#include "gene.h"

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

//- GA begin
int fitness(const Gene &p) { // 适应性
	int cost = mcmf.minCost_Set(p.to_Set());
	// printf("cost = %d\n", cost);
	int Total = mcmf.networkNum * mcmf.costPerCDN;
	if(cost == -1) return 1;
	else return max(1, Total - cost);
}

// 返回一个选中基因的下标
int select(const vector<Gene> & genes) {
	double R = Rand.Random_Real(0, 1);
	double s = 0.0;
	size_t genes_lens = genes.size();
	for(size_t i = 0; i < genes_lens; ++i) {
		s += genes[i].P;
		// printf("%f/%f\n", s, R);
		if(s >= R) {
			// printf("select %d\n", i);
			return i;
		}
	}
	return 0;
}

void GA(int geneCnt = 50, double crossP = 0.95, double mutationP = 0.15) { // 遗传算法
	int iterationCnt = 0;
	int minCost = MCMF::INF;

	vector<Gene> genes(geneCnt);
	vector<Gene> next_genes(geneCnt);
	priority_queue<Gene> que; // 最大堆选出最强的那20条染色体
	unordered_set<int> inital;
	// 初始化基因
	for(int u=0; u < mcmf.consumerNum; ++u)  // 初始位置
		inital.insert(mcmf.edges[mcmf.G[u + mcmf.networkNum][0]].to);
	genes[0].set(inital, mcmf.networkNum);

	for(int i = 1; i < geneCnt; ++i)
		genes[i].reset(mcmf.networkNum);


	while(runing && iterationCnt < 300) {

		// for(int i = 0; i < geneCnt; ++i) {
			// printf("基因型%d: ", i);
			// genes[i].show();
		// }

		// 适应度计算
		int sum = 0;
		for(int i = 0; i < geneCnt; ++i) {
			genes[i].fitness = fitness(genes[i]);
			sum += genes[i].fitness;
			minCost = min(minCost, mcmf.networkNum * mcmf.costPerCDN - genes[i].fitness);
			que.push(genes[i]); // 最大堆
		}

		for(int i = 0; i < geneCnt; ++i)
			genes[i].P = genes[i].fitness*1.0 / sum;

		next_genes.clear();

		// 选择
		for(int i = 0; i < geneCnt; ++i) {
			if(que.size() > geneCnt * 0.6) next_genes[i] = que.top();
			else next_genes[i] = genes[select(genes)];
			que.pop();
		}

		for(int i = 0; i < geneCnt; ++i) // 复制
			genes[i] = next_genes[i];

		// XXOO
		for(int i = 0; i < geneCnt; i+=2)
			if(Rand.Random_Real(0, 1) < crossP)
				genes[i] * genes[i+1];

		// 突变
		for(int i = 0; i < geneCnt; ++i)
			if(Rand.Random_Real(0, 1) < mutationP)
				genes[i].mutation();

		++iterationCnt;
		// printf("iterationCnt: %d minCost = %d\n", iterationCnt, minCost);
		// break;
	}

	// mcmf.showSolution();
	printf("iterationCnt=%d\n", iterationCnt);
	printf("minCost: %d/%d\n\n", minCost, mcmf.consumerNum * mcmf.costPerCDN);
}


//- GA end

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
	while(runing && T > 0.1) {
		int u = -1;
		do {
			for(auto x: backup) {
				if(Rand.Random_Real(0, 1) <  1.0 / mcmf.networkNum) {
					u = x;
					break;
				}
			}
		} while(u == -1);

		int selectEdge = 0, v; // (u, v)随机选点

		do {
			for(selectEdge = 0; (selectEdge < (int)mcmf.G[u].size() - 1) &&
					Rand.Random_Real(0, 1) > 1.0 / mcmf.G[u].size(); ++selectEdge);
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
			if(dC < 0 || exp(-dC / T) > Rand.Random_Real(0, 1))  {// 接受
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


unordered_set<int> Tabu(unordered_set<int>init = {}, int times = MCMF::INF) { // 禁忌搜索
	typedef unordered_set<int> X;
	list<int> H; // 禁忌表，队列

	pair<int, X> x_best;
	X x_now;
	if(init.empty()) {
		for(int u=0; u < mcmf.consumerNum; ++u)  // 初始位置
			x_now.insert(mcmf.edges[mcmf.G[u + mcmf.networkNum][0]].to);
	} else x_now = move(init);

	pair<int, X> x_next{MCMF::INF, {}}; // 转移
	H.push_back(x_best.first = mcmf.minCost_Set(x_now));

	// for(int x: x_now)
		// printf("%d ", x);
	// puts("");

	int iterationCnt = 0;
	while(runing && iterationCnt < times) {
		int Len = 0;
		for(int u: x_now) {
			size_t mcmf_Gu_size = mcmf.G[u].size();
			for(size_t i = 0; i < mcmf_Gu_size && runing; i+=2) {
				++Len;
				int v = mcmf.edges[mcmf.G[u][i]].to; // u->v
				if(v < mcmf.networkNum) {
					X tmp{}; // 邻居
					for(int uu: x_now) {
						if(uu != u) tmp.insert(uu);
						else tmp.insert(v);
					}

					int cost = mcmf.minCost_Set(tmp);
					if(cost == -1) continue;

					if(find(H.begin(), H.end(), cost) == H.end() && cost < x_next.first) {
						x_next.first = cost;
						x_next.second = move(tmp);
						if(x_best.first > cost) {
							x_best.first = cost;
							x_best.second = x_next.second;
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

	printf("iterationCnt = %d\n", iterationCnt);
	printf("minCost: %d/%d cdnNum: %ld\n\n", x_best.first, mcmf.consumerNum * mcmf.costPerCDN, x_best.second.size());
	return x_best.second;
}


void deploy_server(char * topo[MAX_EDGE_NUM], int line_num,char * filename)
{
	Signal(SIGALRM, timeOutHandler);
	// 启动计时器
	alarm(88);
	mcmf.loadGraph(topo, line_num);
	// SA(Tabu({}, 20));
	// SA();
	// GA();
	if(mcmf.networkNum < 200) GA();
	else SA();

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
