#include "deploy.h"
#include <stdio.h>
#include "random.h"
#include "mcmf.h"
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
	for(size_t i = 0; i < genes.size(); ++i) {
		s += genes[i].P;
		// printf("%f/%f\n", s, R);
		if(s >= R) {
			// printf("select %d\n", i);
			return i;
		}
	}
	return 0;
}

void GA(int geneCnt = 20, double retain = 12, double crossP = 0.95, double mutationP = 0.15) { // 遗传算法
	// 初始基因数，精英保留(geneCnt-retain)，交叉率，变异率
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


	while(runing && iterationCnt < 800) {

		// for(int i = 0; i < geneCnt; ++i) {
			// printf("基因型%d: ", i);
			// genes[i].show();
		// }

		// 适应度计算
		int sum = 0;
		for(int i = 0; i < geneCnt; ++i) {
			genes[i].fitness = fitness(genes[i]);
			sum += genes[i].fitness;
			minCost = min<double>(minCost, mcmf.networkNum * mcmf.costPerCDN - genes[i].fitness);
			que.push(genes[i]); // 最大堆
		}

		for(int i = 0; i < geneCnt; ++i)
			genes[i].P = genes[i].fitness*1.0 / sum;

		next_genes.clear();

		// 选择
		for(int i = 0; i < geneCnt; ++i) {
			if(que.size() > retain) next_genes[i] = que.top();
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

int SA(unordered_set<int>init = {}, double T = 20.0, double delta = 0.99999, double poi = 0.02) { // 模拟退火，初始温度，迭代系数，0.15的增点概率
	// double T = 20.0, delta = 0.99999; // 初始温度20, 0.999-0.999999

	unordered_set<int> backup, cur;

	if(init.empty()) {
		for(int u=0; u < mcmf.consumerNum; ++u)  // 初始位置
			backup.insert(mcmf.edges[mcmf.G[u + mcmf.networkNum][0]].to);
	} else backup = move(init);

	int minCost = MCMF::INF, backCost = MCMF::INF, curCost = MCMF::INF;
	backCost = mcmf.minCost_Set(backup);
	minCost = min(minCost, backCost);


	int iterationCnt = 0, poiCnt = 0;
	while(runing && T > 0.1) {
		//- 随机选点u
		int u = -1;
		int i = Rand.Random_Int(0, backup.size() - 1);
		auto it = backup.begin();
		for(; it != backup.end() && i; ++it, --i);
		u = *it;
		// - 选完了

		// 随机选u->v
		int v = -1;
		do {
			v = mcmf.edges[mcmf.G[u][Rand.Random_Int(0, mcmf.G[u].size() - 1)]].to; // (u, v)随机选点
		} while(v >= mcmf.networkNum); // 防止移动到消费节点
		// - 选完v了

		for(int x: backup) {
			if(x == u) cur.insert(v);
			else cur.insert(x);
		}

		if(Rand.Random_Real(0, 1) < poi) {
			++poiCnt;
			cur.insert(Rand.Random_Int(0, mcmf.networkNum - 1)); // 增加一个点
		}

		curCost = mcmf.minCost_Set(cur);
		++iterationCnt;

		if(curCost == -1)  {// 无解
			cur.clear();
		}
		else {
			int dC = curCost - backCost;
			// printf("dC: %d\n", dC);
			if(min(1.0, exp(-dC / T)) > Rand.Random_Real(0, 1))  {// 接受
				backup = move(cur);
				backCost = curCost;
			} else {
				cur.clear();
			}
			minCost = min(minCost, backCost);
	 	}
		T *= delta;

		// printf("T=%lf iterationCnt=%d minCost = %d\n", T, iterationCnt, minCost);
	}

	printf("T=%lf iterationCnt=%d poiCnt=%d\n", T, iterationCnt, poiCnt);
	// printf("Deploy CDN(%ld):\n", backup.size());
	// for(int x: backup)
		// printf("%d ", x);
	// puts("\n=====Solution======");
	// mcmf.showSolution();
	printf("minCost: %d/%d cdnNum: %ld\n\n", minCost, mcmf.consumerNum * mcmf.costPerCDN, backup.size());
	return minCost;
}

//- SAGA
void SAGA(unordered_set<int>init = {}, double T = 20.0, double poi = 0.05, double delta = 0.999, int geneCnt = 26, double crossP = 0.95, double mutationP = 0.15) { // 模拟退火，初始温度，迭代系数
	// double T = 20.0, delta = 0.99999; // 初始温度20, 0.999-0.999999

	unordered_set<int> initial;
	vector<Gene> genes(geneCnt);
	vector<Gene> next_genes(geneCnt);

	if(init.empty()) {
		for(int u=0; u < mcmf.consumerNum; ++u)  // 初始位置，直连
			initial.insert(mcmf.edges[mcmf.G[u + mcmf.networkNum][0]].to);
	} else initial = move(init);

	int minCost = MCMF::INF;

	// for(int i = 0; i < geneCnt; ++i)
		// genes[i].set(initial, mcmf.networkNum);
	genes[0].set(initial, mcmf.networkNum);
	for(int i = 1; i < geneCnt; ++i)
		genes[i].reset(mcmf.networkNum);


	int iterationCnt = 0;
	Gene elite; // 精英基因
	while(T > 0.1) {
		next_genes.clear();

		int fmin = MCMF::INF;
		for(int idx = 0; runing && idx < geneCnt; ++idx) {
			unordered_set<int> s = genes[idx].to_Set(); // 每条染色体
			int fi = mcmf.minCost_Set(s), fj;
			unordered_set<int> cur; // 邻域
			// 计算领域

			//- 随机选点u
			int u = -1;
			int i = Rand.Random_Int(0, s.size() - 1);
			auto it = s.begin();
			for(; it != s.end() && i; ++it, --i);
			u = *it;
			// - 选完了

			// 随机选u->v
			int v = -1;
			do {
				v = mcmf.edges[mcmf.G[u][Rand.Random_Int(0, mcmf.G[u].size() - 1)]].to; // (u, v)随机选点
			} while(v >= mcmf.networkNum); // 防止移动到消费节点
			// - 选完v了

			for(int x: s) {
				if(x == u) cur.insert(v);
				else cur.insert(x);
			}

			if(Rand.Random_Real(0, 1) < poi)
				cur.insert(Rand.Random_Int(0, mcmf.networkNum - 1)); // 增加一个点
			// 邻域计算完毕

			fj = mcmf.minCost_Set(cur);

			if(fj != -1)  {// 有解
				int dC = fj - fi;
				// printf("dC: %d\n", dC);
				if(fi == -1 || min<double>(1, exp(-dC / T)) > Rand.Random_Real(0, 1)) {// 接受
					genes[idx].set(cur, mcmf.networkNum);
					genes[idx].fitness = fj;
				} else
					genes[idx].fitness = fi;

				if(fmin > fj) {
					fmin = fj;
					elite = genes[idx];
				}
			} else { // 无解，不接受
				genes[idx].fitness = (fi == -1?mcmf.networkNum * mcmf.costPerCDN:fi);
			}
		}

		// 计算适应度
		double sum = 0.0;
		for(int idx = 0; idx < geneCnt; ++idx) {
			int dC = genes[idx].fitness - fmin;
			genes[idx].fitness = exp(-dC / T);
			sum += genes[idx].fitness;
		}

		for(int idx = 0; idx < geneCnt; ++idx)
			genes[idx].P = genes[idx].fitness / sum;

		// 轮盘赌选择
		next_genes[0] = elite; // 精英
		for(int idx = 1; idx < geneCnt; ++idx)
			next_genes[idx] = genes[select(genes)];


		for(int idx = 0; idx < geneCnt; ++idx)
			genes[idx] = next_genes[idx];

		// XXOO
		for(int i = 0; i < geneCnt; i+=2)
			if(Rand.Random_Real(0, 1) < crossP)
				genes[i] * genes[i+1];

		// 突变
		for(int i = 0; i < geneCnt; ++i)
			if(Rand.Random_Real(0, 1) < mutationP)
				genes[i].mutation();

		minCost = min(minCost, fmin);
		T *= delta;

		++iterationCnt;
		// printf("T=%lf iterationCnt=%d minCost = %d\n", T, iterationCnt, minCost);
	}

	printf("T=%lf iterationCnt=%d\n", T, iterationCnt);
	// mcmf.showSolution();
	printf("minCost: %d/%d\n\n", minCost, mcmf.consumerNum * mcmf.costPerCDN);
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
			for(size_t i = 0; i < mcmf.G[u].size() && runing; i+=2) {
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
	// SA({}, 20, 0.99999, 0.02);
	// GA();
	// SAGA();

	// 初始解{}，初始温度，增点概率，迭代系数，基因数，交叉率，变异率
	if(mcmf.networkNum < 200)
		SAGA({}, 20, 0.01, 0.99, 30, 0.95, 0.15);
	else if(mcmf.networkNum < 500)
		SAGA({}, 20, 0.01, 0.999, 26, 0.95, 0.15);
	else
		SAGA({}, 20, 0.01, 0.999, 6, 0.95, 0.15);

	// unordered_set<int> cdn{0, 3, 22};
	// printf("cost = %d\n", mcmf.minCost_Set(cdn));


	//- test
	/*
	double T = 20.0, delta = 0.99999, poi = 0.02;
	double bestT = T, bestDelta = delta, bestPoi = poi;
	int minCost = MCMF::INF;
	int cost = 0;
	// for(; T <= 100.0; T+=1) {
	for(poi = 0.01; poi <= 1; poi += 0.01) {
		alarm(88);
		if( (cost = SA({}, T,delta, poi)) < minCost && cost != -1) {
			minCost = cost;
			bestT = T;
			bestDelta = delta;
			bestPoi = poi;
		}
		puts("--------------------");
		printf("bestT = %lf/%lf bestDelta = %lf/%lf bestPoi = %lf/%lf minCost = %d\n", bestT, T, bestDelta, delta, bestPoi, poi, minCost);
		puts("--------------------");
		runing = true;
	}
	*/
	//- test End

	// 开始计算
	write_result(mcmf.outputPath(), filename);

}
