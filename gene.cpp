/*************************************************************************
	> File Name: gene.cpp
	  > Author: Netcan
	  > Blog: http://www.netcan666.com
	  > Mail: 1469709759@qq.com
	  > Created Time: 2017-03-26 Sun 20:31:48 CST
 ************************************************************************/

#include <cstdio>
#include <ctime>
#include <cstdlib>
#include <algorithm>
using namespace std;

template<class T>
class Gene {
	private:
		const int len; // 长度，0-1200
		const int codeSize = sizeof(T) * 8; // T的位数
		T code[1200 / sizeof(T) + 1];
	public:
		Gene(int len) {
			for(int i = 0; i*codeSize < len; ++i)
				code[i] = rand() % ((1 << codeSize) - 1);
		}
		Gene(T *c, int len): len(len) {
			for(int i = 0; i*codeSize < len; ++i)
				code[i] = c[i];
		}
		void operator*(Gene &b) {
			int loc = rand() % len + 1; // 交换的位置，交换一边就行了，因为另一边不动，这里交换左边
			// int loc = 12;
			// printf("loc = %d\n", loc);
			int i;
			for(i = 1; i*codeSize < loc; ++i)
				swap<T>(code[i-1], b.code[i-1]);
			--i;
			for(int j = 0; j < codeSize && i * codeSize + j < loc; ++j) {
				if( (code[i] & (1 << (codeSize - j - 1))) == (b.code[i] & (1 << (codeSize - j - 1))) )
					continue;
				else {
					code[i] ^= 1 << (codeSize - j - 1);
					b.code[i] ^= 1 << (codeSize - j - 1);
				}
			}

		}

		void mutation(int loc = -1) { // 突变，[0, len)
			if(loc == -1) loc = rand() % len;
			else if(loc >= len) return;
			printf("loc = %d\n", loc);

			code[loc / codeSize] ^= 1 << (codeSize - loc % codeSize - 1);
		}

		void show() {
			for(int i = 0; i*codeSize < len; ++i) {
				if(i != 0) printf(",");
				for(int j = 0; j < codeSize && i * codeSize + j < len; ++j)
					printf( (code[i] & (1 << (codeSize - j - 1)) )?"1":"0");
			}
			puts("");
		}

};


int main(void) {
	typedef unsigned char T;
	srand(time(0));
	T a[] = {197, 51};
	T b[] = {255, 160};

	Gene<T> ga(a, 11);
	Gene<T> gb(b, 11);
	// ga.show();
	// gb.show();
	// ga * gb;
	// ga.show();
	// gb.show();

	ga.show();
	ga.mutation(10);
	ga.show();

	return 0;
}
