#ifndef __GRAPH_TOOLS_H__
#define __GRAPH_TOOLS_H__

#include "graphGenetic.h"

int** createAdjMatrix(int nSommets);
void showMatrix(graph_genetic_t* graphe);
void freeMatrix(graph_genetic_t* g);


#endif
