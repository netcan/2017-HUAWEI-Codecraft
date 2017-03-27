#include<cstring>
#include<cstdio>
#include<cmath>
#define rand_01 ((float)rand() / (float)RAND_MAX)

const int numofdims = 30;
const int numofparticles = 50;

void fitnessfunc(float X[numofparticles][numofdims], float fitnesses[numofparticles])
{
	memset(fitnesses, 0, sizeof (float) * numofparticles);
	for(int i = 0; i < numofparticles; i++)
	{
		for(int j = 0; j < numofdims; j++)
		{
			fitnesses[i] += X[i][j] * X[i][j]; 
		}
	}
}
