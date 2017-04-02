/*************************************************************************
	> File Name: random.cpp
	  > Author: Netcan
	  > Blog: http://www.netcan666.com
	  > Mail: 1469709759@qq.com
	  > Created Time: 2017-03-27 Mon 07:40:54 CST
 ************************************************************************/

#ifndef __RANDOM__
#define __RANDOM__
#include <cstdio>
#include <random>

class Random {
	private:
		std::default_random_engine generator;

	public:
		Random(time_t t = time(NULL)): generator(t) {}
		// Random(time_t t = 1491132293): generator(t) {}
		inline uint32_t Random_Int(uint32_t min, uint32_t max) { // [min, max]
			return std::uniform_int_distribution<uint32_t>{min, max}(generator);
		}
		inline double Random_Real(double min, double max) { // [min, max)
			return std::uniform_real_distribution<double>{min, max}(generator);
		}
};
Random Rand;

#endif
