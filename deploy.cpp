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

// 直连状态
unordered_set<int> directConn() {
	static unordered_set<int> direct;
	if(direct.empty()) {
		for(int u=0; u < mcmf.consumerNum; ++u)  // 初始位置，直连
			direct.insert(mcmf.edges[mcmf.G[u + mcmf.networkNum][0]].to);
	}
	return direct;
}
// XJBS
bool cmp(int u1, int u2) { // 比较函数，消费降低需要的流量越低，排在越前
	int u1cap = 0, u2cap = 0;
	for(size_t i = 0; i < mcmf.G[u1].size(); ++i)
		if(mcmf.isConsumer(mcmf.edges[mcmf.G[u1][i]].to)) {
			u1cap = mcmf.edges[mcmf.G[u1][i]].cap;
			break;
		}

	for(size_t i = 0; i < mcmf.G[u2].size(); ++i)
		if(mcmf.isConsumer(mcmf.edges[mcmf.G[u2][i]].to)) {
			u2cap = mcmf.edges[mcmf.G[u2][i]].cap;
			break;
		}

	return u1cap < u2cap;
}

unordered_set<int> XJBS(bool sorted = false) {
	unordered_set<int> init = directConn();

	vector<int> tmp(init.begin(), init.end());
	if(sorted) sort(tmp.begin(), tmp.end(), cmp);
	list<int> cdn(tmp.begin(), tmp.end());

	int minCost = mcmf.minCost_Set(unordered_set<int>(cdn.begin(), cdn.end()));

	// 删点
	int iterationCnt = 0;
	for(auto itr = cdn.begin(); itr != cdn.end(); ) {
		int node = *itr;
		int cost = -1;
		itr = cdn.erase(itr);
		++iterationCnt;
		if( (cost = mcmf.minCost_Set(unordered_set<int>(cdn.begin(), cdn.end()))) < minCost && cost != -1) {
			minCost = cost;
			// printf("deleted: %d\n", node);
		}
		else {
			itr = cdn.insert(itr, node); // 恢复
			++itr;
		}
		// printf("cost: %d\n", minCost);
	}

	// 替换
	/*
	for(auto itr = cdn.begin(); itr != cdn.end(); ++itr) {
		int u = *itr;
		for(size_t i = 0; i < mcmf.G[u].size(); ++i) {
			int cost = -1;
			int v = mcmf.edges[mcmf.G[u][i]].to;
			if(v < mcmf.networkNum) {
				*itr = v; // 替换
				if( (cost = mcmf.minCost_Set(unordered_set<int>(cdn.begin(), cdn.end()))) < minCost && cost != -1) {
					minCost = cost;
					printf("replace %d with %d\n", u, v);
				}
				else {
					*itr = u;
				}
			}

			size_t next = mcmf.G[u][i] + 1;
			if( next < mcmf.G[u].size() && mcmf.edges[next].to == v) ++i;
		}
		// printf("minCost: %d/%d\n", minCost, mcmf.consumerNum * mcmf.costPerCDN);
	}
	*/

	// printf("minCost: %d/%d iterationCnt: %d\n", minCost, mcmf.consumerNum * mcmf.costPerCDN, iterationCnt);
	// printf("cdn sz: %ld\n", cdn.size());
	return unordered_set<int>(cdn.begin(), cdn.end());
}

// 按评估值从高到低选址
unordered_set<int> evaluationSelect() {
	vector<pair<double, int>> evaluation;
	unordered_set<int> cdn{};
	for(int u = 0; u < mcmf.networkNum; ++u)
		evaluation.push_back(make_pair(mcmf.nodes[u].evaluation, u));
	sort(evaluation.begin(), evaluation.end(), greater_equal<pair<double, int>>());

	for(int i = 0; i < mcmf.networkNum; ++i) {
		cdn.insert(evaluation[i].second);
		// printf("%d: %lf\n", evaluation[i].second, evaluation[i].first);
		if(mcmf.minCost_Set(cdn) != -1) break;
	}
	mcmf.showRealMinCost();
	return cdn;
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

void GA(unordered_set<int> init = {}, int geneCnt = 20, double retain = 12, double crossP = 0.95, double mutationP = 0.25) { // 遗传算法
	// 初始基因数，精英保留(geneCnt-retain)，交叉率，变异率
	int iterationCnt = 0;
	int minCost = MCMF::INF;

	vector<Gene> genes(geneCnt);
	vector<Gene> next_genes(geneCnt);
	priority_queue<Gene> que; // 最大堆选出最强的那20条染色体
	unordered_set<int> initial;
	if(init.empty()) initial = directConn();
	else initial = move(init);

	// 初始化基因
	genes[0].set(initial, mcmf.networkNum);

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
	// printf("minCost: %d/%d\n\n", minCost, mcmf.consumerNum * mcmf.costPerCDN);
	mcmf.showRealMinCost();
}


//- GA end

//- 模拟退火 begin
int select(const vector<pair<int, double>> & cdn) {
	double R = Rand.Random_Real(0, 1);
	double s = 0.0;
	for(size_t i = 0; i < cdn.size(); ++i) {
		s += cdn[i].second;
		// printf("%f/%f\n", s, R);
		if(s >= R) {
			// printf("select %d\n", i);
			return i;
		}
	}
	return 0;
}
unordered_set<int> SA(unordered_set<int>init = {}, int innerLoop = 10, double T = 20.0, double delta = 0.99999, double poi = 0.02) { // 模拟退火，初始温度，迭代系数，0.15的增点概率
	// double T = 20.0, delta = 0.99999; // 初始温度20, 0.999-0.999999

	unordered_set<int> backup, cur, best;

	if(init.empty()) backup = directConn();
	else backup = move(init);

	int minCost = MCMF::INF, backCost = MCMF::INF, curCost = MCMF::INF;
	backCost = mcmf.minCost_Set(backup);
	minCost = min(minCost, backCost);

	int iterationCnt = 0;
	while(runing && T > 0.1) {

		for(int loop = 0; loop < innerLoop && runing; ++loop) {
			vector<pair<int, double>> cdn; // cdn选中的概率，概率越大，越容易被选中
			double sum = 0.0;
			int u = -1, v = -1;
			// 随机选点u->v
			for(auto x: backup)
				sum += mcmf.nodes[x].evaluation;
			for(auto x: backup)
				cdn.push_back(make_pair(x, mcmf.nodes[x].evaluation / sum));
			u = cdn[select(cdn)].first;

			cdn.clear();
			sum = 0.0;
			for(size_t i = 0; i < mcmf.G[u].size(); ++i) {
				int t = mcmf.edges[mcmf.G[u][i]].to;
				if(t < mcmf.networkNum)
					sum += 100000.0 / mcmf.nodes[t].evaluation;
				// printf("%d %lf\n", t, mcmf.nodes[t].evaluation);
			}
			for(size_t i = 0; i < mcmf.G[u].size(); ++i) {
				int t = mcmf.edges[mcmf.G[u][i]].to;
				if(t < mcmf.networkNum)
					cdn.push_back(make_pair(t, (100000.0 / mcmf.nodes[t].evaluation) / sum));
			}
			v = cdn[select(cdn)].first;


			for(int x: backup) {
				if(x == u) { // u->v，新点v
					// printf("%d->%d\n", u, v);
					if(cur.count(v)) { // 合并uv
						int flow = mcmf.servers[mcmf.nodes[v].bestCdnId].outFlow + mcmf.servers[mcmf.nodes[u].bestCdnId].outFlow;
						vector<MCMF::Server>::iterator it;
						if( (it = lower_bound(mcmf.servers.begin(), mcmf.servers.end(), flow))  != mcmf.servers.end()) // >= 升档
							mcmf.nodes[v].bestCdnId = it - mcmf.servers.begin(); // 存放下标，nodes输出路径的时候用
						else mcmf.nodes[v].bestCdnId = mcmf.servers.size() - 1;

					} else cur.insert(v); // 新点
				}
				else cur.insert(x);
			}

			if(Rand.Random_Real(0, 1) < poi)
				cur.insert(Rand.Random_Int(0, mcmf.networkNum - 1)); // 增加一个点

			curCost = mcmf.minCost_Set(cur);

			if(curCost == -1) { // 无解继续
				cur.clear();
				continue;
			}

			// 对新点进行降档处理
			for(mcmf.nodes[v].bestCdnId -= 1; mcmf.nodes[v].bestCdnId >= 0; --mcmf.nodes[v].bestCdnId) {
				int tmp = mcmf.minCost_Set(cur);
				if(tmp == -1 || tmp >= curCost) break;
				curCost = tmp;
				// printf("%d\n", curCost);
			}
			++mcmf.nodes[v].bestCdnId;

			++iterationCnt;

			int dC = curCost - backCost;
			// printf("dC: %d ratio: %lf probability: %lf\n", dC, curCost * 1.0 / minCost, exp(-dC / T));
			if(min(1.0, exp(-dC / T)) > Rand.Random_Real(0, 1))  {// 接受
				// printf("T: %lf dC: %d ratio: %lf\n", T, dC, curCost * 1.0 / minCost);
				backup = move(cur);
				backCost = curCost;
			} else {
				cur.clear();
			}

			if(minCost > backCost) {
				minCost = backCost;
				best = backup;
				// mcmf.showRealMinCost();
			}
		}
		T *= delta;

		// printf("T=%lf iterationCnt=%d minCost = %d\n", T, iterationCnt, minCost);
	}

#ifdef _DEBUG
	printf("T=%lf iterationCnt=%d\n", T, iterationCnt);
	mcmf.showRealMinCost();
#endif
	// printf("Deploy CDN(%ld):\n", backup.size());
	// for(int x: backup)
		// printf("%d ", x);
	// puts("\n=====Solution======");
	// mcmf.showSolution();
	return best;
}
//- 模拟退火 end

//- SAGA begin
void SAGA(unordered_set<int>init = {}, double T = 20.0, double poi = 0.05, double delta = 0.999, int geneCnt = 26, double crossP = 0.95, double mutationP = 0.15) { // 模拟退火，初始温度，迭代系数
	// double T = 20.0, delta = 0.99999; // 初始温度20, 0.999-0.999999

	unordered_set<int> initial;
	vector<Gene> genes(geneCnt);
	vector<Gene> next_genes(geneCnt);

	if(init.empty()) initial = directConn();
	else initial = move(init);

	int minCost = MCMF::INF;

	// for(int i = 0; i < geneCnt; ++i)
		// genes[i].set(initial, mcmf.networkNum);
	genes[0].set(initial, mcmf.networkNum);
	unordered_set<int> direct = directConn();
	for(int i = 1; i < geneCnt; ++i)
		genes[i].reset(mcmf.networkNum);


	int iterationCnt = 0;
	// 忘记初始化了！导致段错误！！
	Gene elite{mcmf.networkNum}; // 精英基因
	while(runing && T > 0.1) {
		next_genes.clear();
		int fmin = MCMF::INF;

		for(int idx = 0; runing && idx < geneCnt; ++idx) {
			unordered_set<int> s = genes[idx].to_Set(); // 每条染色体
			if(s.empty()) continue; // 空集的时候需要跳过

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
					if(fmin > fj) {
						fmin = fj;
						elite = genes[idx];
					}
				} else {
					genes[idx].fitness = fi; // 不接收
					if(fmin > fi) {
						fmin = fi;
						elite = genes[idx];
					}
				}
			} else { // 无解，不接受
				if(fmin > fi && fi != -1) {
					fmin = fi;
					elite = genes[idx];
				}
				genes[idx].fitness = (fi == -1?mcmf.networkNum * mcmf.costPerCDN:fi);
			}
		}

		// 计算适应度
		double sum = 0.0;
		for(int idx = 0; runing && idx < geneCnt; ++idx) {
			if(fmin == MCMF::INF) fmin = 0;
			int dC = genes[idx].fitness - fmin;
			genes[idx].fitness = exp(-dC / T);
			sum += genes[idx].fitness;
		}

		for(int idx = 0; runing && idx < geneCnt; ++idx)
			genes[idx].P = genes[idx].fitness / sum;

		// 轮盘赌选择
		next_genes[0] = elite; // 精英

		for(int idx = 1; runing && idx < geneCnt; ++idx)
			next_genes[idx] = genes[select(genes)];

		for(int idx = 0; runing && idx < geneCnt; ++idx)
			genes[idx] = next_genes[idx];

		// 洗牌，打乱顺序，考虑是否必要
		// random_shuffle(genes.begin(), genes.end());
		// XXOO
		for(int i = 0; runing && i < geneCnt; i+=2)
			if(Rand.Random_Real(0, 1) < crossP)
				genes[i] * genes[i+1];

		// 突变
		for(int i = 0; runing && i < geneCnt; ++i)
			if(Rand.Random_Real(0, 1) < mutationP)
				genes[i].mutation();

		if(fmin != 0) minCost = min(minCost, fmin);
		T *= delta;

		++iterationCnt;
		// printf("minCost: %d/%d\n\n", minCost, mcmf.consumerNum * mcmf.costPerCDN);
	}

	printf("T=%lf iterationCnt=%d\n", T, iterationCnt);
	// mcmf.showSolution();
	// printf("minCost: %d/%d\n\n", minCost, mcmf.consumerNum * mcmf.costPerCDN);
	mcmf.showRealMinCost();
}
//- SAGA end

//- 禁忌搜索 begin
// 这块没写好，效果太差
unordered_set<int> Tabu(unordered_set<int>init = {}, int times = MCMF::INF) { // 禁忌搜索
	typedef unordered_set<int> X;
	list<int> H; // 禁忌表，队列

	pair<int, X> x_best;
	X x_now;
	if(init.empty()) x_now = directConn();
	else x_now = move(init);

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
//- 禁忌搜索 end

//- BPSO begin
// 这个BPSO有坑，效果没想象中的好
double sig(double v, double Vmax, double Vmin) { // v->[0, 1]
	return 1/(1+ pow((Vmax - v)/(v- Vmin), 2));
}

void BPSO(unordered_set<int> init = {}, int particleCnt = 10, double Vmin = 0.0, double Vmax = 10.0, double c1 = 1.0, double c2 = 1.0) {
	vector<Particle> particles(particleCnt);
	Particle pBest, gBest;
	int fpBest = -1, fgBest = -1, fCur; // 当代最小费用，全局最小费用
	vector<double> v[particleCnt]; // Vij
	unordered_set<int> initial;

	if(init.empty()) initial = directConn(); // 初始状态
	else initial = move(init);

	for(int i = 0; i < particleCnt; ++i) {
		if(i) particles[i].reset(mcmf.networkNum);
		else particles[i].set(initial, mcmf.networkNum); // 直连状态

		for(int j = 0; j < mcmf.networkNum; ++j) // 初始化速度
			v[i].push_back(Vmin + (Vmax - Vmin) * Rand.Random_Real(0, 1));
	}

	// for(int i = 0; i < particleCnt; ++i)
		// for(int j = 0; j < mcmf.networkNum; ++j)
			// printf("%lf\n", v[i][j]);


	int iterationCnt = 0;
	while(runing) {
		for(int i = 0; i < particleCnt; ++i) {
			if( (fCur = mcmf.minCost_Set(particles[i].to_Set())) != -1) {
				if(fpBest == -1 || fpBest > fCur) {
					fpBest = fCur;
					pBest = particles[i];
				}
				if(fgBest == -1 || fgBest > fpBest) {
					fgBest = fpBest;
					gBest = pBest;
				}
			}
		}

		for(int i = 0; i < particleCnt; ++i) {
			// puts("------------------------");
			// printf("p[%d]: \n", i);
			// particles[i].show();
			for(int j = 0; j < mcmf.networkNum; ++j) {
				v[i][j] = v[i][j] +
					c1 * Rand.Random_Real(0, 1) * (pBest.getBit(j) - particles[i].getBit(j)) +
					c2 * Rand.Random_Real(0, 1) * (gBest.getBit(j) - particles[i].getBit(j));
				if(v[i][j] > Vmax) v[i][j] = Vmax;
				else if(v[i][j] < Vmin) v[i][j] = Vmin;
				// printf("sig(%lf) = %lf\n", v[i][j], sig(v[i][j], Vmax, Vmin));

				if(Rand.Random_Real(0, 1) < sig(v[i][j], Vmax, Vmin)) particles[i].setBit(j, 1);
				else particles[i].setBit(j, 0);
			}
			// particles[i].show();
		}

		fpBest = -1;
		++iterationCnt;
		// printf("iterationCnt = %d\n", iterationCnt);
		// printf("minCost: %d/%d\n", fgBest, mcmf.consumerNum * mcmf.costPerCDN);
		// break;
	}

	printf("iterationCnt = %d\n", iterationCnt);
	printf("minCost: %d/%d\n", fgBest, mcmf.consumerNum * mcmf.costPerCDN);
}
//- BPSO end


void deploy_server(char * topo[MAX_EDGE_NUM], int line_num,char * filename)
{
	Signal(SIGALRM, timeOutHandler);
	// 启动计时器
	alarm(88);
	mcmf.loadGraph(topo, line_num);

	if(mcmf.networkNum < 800){
		unordered_set<int> s = SA({}, 1, 200, 0.99999, 0.00);
		mcmf.showRealMinCost();
		// GA(XJBS(true));
		// SAGA(XJBS(true), 200, 0.00, 0.99, 20, 0.95, 0.05);
	} else {
		unordered_set<int> s = SA({}, 1, 500, 0.9999, 0.00);
		mcmf.showRealMinCost();
	}


	// SA(Tabu({}, 20));
	// GA(XJBS(true));
	// SAGA();
	// BPSO(XJBS(true));
	// XJBS();

	// 初始解{}，初始温度，增点概率，迭代系数，基因数，交叉率，变异率
	// if(mcmf.networkNum < 200) {
		// mcmf.setCostPerCdnMethod(false); // 动态变动
		// SAGA(XJBS(), 2000, 0.00, 0.99, 30, 0.8, 0.05);
	// }
	// else if(mcmf.networkNum < 500) {
		// mcmf.setCostPerCdnMethod(false); // 服务器费用固定
		// SAGA(XJBS(false), 2000, 0.00, 0.99, 50, 0.8, 0.05);
	// }
	// else {
		// mcmf.setCostPerCdnMethod(false); // 服务器费用固定
		// SAGA(XJBS(true), 20, 0.00, 0.999, 6, 0.8, 0.05);
	// }

	// unordered_set<int> cdn{
		// 0, 45, 55, 56, 60, 78, 105, 107, 133, 134, 142, 152, 161, 177, 236, 242, 245, 274, 278, 290, 291, 296, 314, 333, 343, 359, 373, 389, 390, 394, 409, 411, 416, 445, 458, 460, 467, 470, 495, 497, 515, 518, 526, 527, 538, 556, 557, 570, 577, 582, 586, 597, 615, 617, 625, 640, 641, 650, 656, 657, 666, 669, 683, 688, 697, 700, 714, 724, 751, 767, 804, 835, 847, 872, 883, 894, 920, 934, 940, 952, 970, 984, 991, 993, 1002, 1017, 1019, 1029, 1031, 1032, 1034, 1053, 1056, 1070, 1076, 1080, 1090, 1103, 1109, 1110, 1118, 1119, 1184, 1187
	// };
	// mcmf.setCostCdnGap(1000);
	// mcmf.minCost_Set(cdn);
	// mcmf.showRealMinCost();


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
