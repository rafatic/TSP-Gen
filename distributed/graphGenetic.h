#ifndef __GRAPHE_GENETIC_H__
#define __GRAPHE_GENETIC_H__

#include <stdio.h>

typedef struct graph{
	int nSommets;
	int nArretes;
	int oriente;
	int value; // valué
	int complet;
	int** matriceAdj;

}graph_genetic_t;

graph_genetic_t* createGraph(char* nomFichier);
void fillMatrix(FILE* fichier, char* line, graph_genetic_t* graphe);
graph_genetic_t* createGraphFromTSP(char* nomFichier);
void fillFromEuclideanDistances(graph_genetic_t* graphe, FILE* file, char* line);

#endif
