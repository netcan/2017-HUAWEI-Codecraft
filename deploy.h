#ifndef __ROUTE_H__
#define __ROUTE_H__

#include "lib_io.h"
#include "spfa.h"
#include <signal.h>
#include <unistd.h>
#include <cmath>
#include <sys/types.h>
#include <sys/wait.h>

void deploy_server(char * graph[MAX_EDGE_NUM], int edge_num, char * filename);



#endif
