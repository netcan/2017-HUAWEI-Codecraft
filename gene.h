/*************************************************************************
	> File Name: gene.cpp
	  > Author: Netcan
	  > Blog: http://www.netcan666.com
	  > Mail: 1469709759@qq.com
	  > Created Time: 2017-03-26 Sun 20:31:48 CST
 ************************************************************************/

#ifndef __GENE__
#define __GENE__
#include <cstdio>
#include <vector>
#include <algorithm>
#include <ctime>
#include <bitset>
#include <unordered_set>
#include "random.h"
using namespace std;

class Gene {
	private:
		int len; // 长度，0-1200
		bitset<10000+5> code;
	public:
		double fitness; // 适应度/费用
		double P; // 选中概率
		Gene(): len(0), code(), fitness(0), P(0) {}
		inline void reset(int len) { // 重置，亦即随机
			this->len = len;
			this->P = this->fitness = 0;
			for(int i = 0; i < len; ++i)
				code[i] = Rand.Random_Int(0, 1);
		}

		Gene(int len): len(len), fitness(0), P(0) {}

		inline void operator*(Gene &b) { // 交叉，同时改变2条染色体
			int end = Rand.Random_Int(1, len); // 交换的位置，交换一边就行了，因为另一边不动，这里交换两边
			int begin = Rand.Random_Int(0, end - 1);
			// printf("begin = %d end = %d len = %d\n", begin, end, len);
			for(int i = begin; i < end; ++i)
				if(code[i] != b.code[i]) {
					code[i] = !code[i];
					b.code[i] = !b.code[i];
				}
		}

		inline void set(unordered_set<int> &s, int len) {
			this->len = len;
			std::vector<int> ss(s.begin(), s.end());
			sort(ss.begin(), ss.end()); // 排序
			size_t j = 0;
			size_t ss_size = ss.size();
			for(int i = 0; i < len && j < ss_size; ++i) {
				if(ss[j] == i) {
					code[i] = 1;
					++j;
				} else // 忘记置0了
					code[i] = 0;
			}
		}

		inline bool operator<(const Gene &b) const { // 最小堆用
			return this->fitness < b.fitness;
		}

		inline void operator=(const Gene &b) { // 赋值
			this->len = b.len;
			this->fitness = b.fitness;
			this->P = b.P;
			code = b.code;
		}
		inline bool operator==(const Gene &b)const { // 判断序列是否相等
			return code == b.code;
		}

		inline void mutation(int loc = -1) { // 突变，[0, len)
			if(loc == -1) loc = Rand.Random_Int(0, len - 1);
			else if(loc >= len) return;
			// printf("loc = %d\n", loc);
			code[loc] = !code[loc];
		}

		inline void show() const {
			for(int i = 0; i < len; ++i) {
				if(i != 0 && i % 8 == 0) printf(",");
				printf(code[i]?"1":"0");
			}
			puts("");
		}

		inline unordered_set<int> to_Set() const {
			unordered_set<int> S;
			for(int i = 0; i < len; ++i)
				if(code[i]) S.insert(i);
			return S;
		}
		inline bool getBit(int loc) {
			return code[loc];
		}
		inline void setBit(int loc, bool x) {
			code[loc] = x;
		}
};

typedef Gene Particle;


#endif


// int main(void) {
	// Gene ga(10);
	// Gene gb(10);

	// ga.show();
	// gb.show();
	// ga * gb;
	// ga.show();
	// gb.show();

	// ga.show();
	// ga.mutation();
	// ga.show();

	// for(auto x: ga.to_Set()) {
		// printf("%d\n", x);
	// }

	// return 0;
// }
