#include "graphTools.h"
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <limits.h>


int** createAdjMatrix(int nSommets)
{
	int** matrice = (int**)calloc(nSommets, sizeof(int**));
	int i;
	for(i=0; i < nSommets; i++)
	{
		matrice[i] = (int*)calloc(nSommets, sizeof(int*));
	}

	return matrice;
}


void showMatrix(graph_genetic_t* graphe)
{
	printf("Matrice d'adjacence du graphe :\n");

	printf("    ");

	int i,j,k;
	for(i = 0; i< graphe->nSommets; i++)
	{
		printf("%d  ", i);
	}
	printf("\n");
	for(i = 0; i< graphe->nSommets + 2 ; i++)
	{
		printf("---");
	}
	printf("\n");
	for(j = 0; j< graphe->nSommets; j++)
	{
		if(j < 10)
		 	printf("%d   |", j);
		else
			printf("%d  |", j);
		for(k = 0; k < graphe->nSommets; k++)
		{
			printf("%d  ", graphe->matriceAdj[j][k]);

		}
		printf("\n");
	}
}

void freeMatrix(graph_genetic_t* g)
{
	int i;

	for (i=0; i < g->nSommets; i++)
		free(g->matriceAdj[i]);

	free(g->matriceAdj);
}
