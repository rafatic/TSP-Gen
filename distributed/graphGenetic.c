#define _GNU_SOURCE
#include "graphGenetic.h"
#include "graphTools.h"
#include "debug.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

graph_genetic_t* createGraph(char* nomFichier)
{
	FILE*fichier;
	fichier = fopen(nomFichier, "r");
	char* line = NULL;
    size_t len = 0;
    int stop = 0;
    char* token;



    graph_genetic_t* graphe = (graph_genetic_t*)malloc(sizeof(graph_genetic_t));

	while (stop == 0 && getline(&line, &len, fichier) != -1 )
	{
		token = strtok(line, " ");


		if(strcmp(token, "nSommets") == 0)
		{
			token = strtok(NULL, " ");
			graphe->nSommets = atoi(token);
		}
		else if(strcmp(token, "oriente") == 0)
		{
			token = strtok(NULL, " ");
			graphe->oriente = atoi(token);
		}
		else if(strcmp(token, "value") == 0)
		{
			token = strtok(NULL, " ");
			graphe->value = atoi(token);
		}
		else if(strcmp(token, "complet") == 0)
		{
			token = strtok(NULL, " ");
			graphe->complet = atoi(token);
		}
		else if(strcmp(token, "debutDefAretes\n") == 0)
		{



			fillMatrix(fichier, line, graphe);

			fclose(fichier);
			

			
			stop = 1;
		}

    }
    return graphe;
}


void fillMatrix(FILE* fichier, char* line, graph_genetic_t* graphe)
{
	int indexI, indexJ;
	char* token;
	size_t len = 0;
	int** matrice = createAdjMatrix(graphe->nSommets);



	while ((getline(&line, &len, fichier)) != -1)
	{
		
		token = strtok(line, " ");

		indexI = atoi(token);
		while(token != NULL)
		{
			token = strtok(NULL, " ");

			if(token != NULL)
			{
				indexJ = atoi(token);
				if(graphe->value == 1)
					matrice[indexI][indexJ] = atoi(strtok(NULL, " "));
				else
					matrice[indexI][indexJ] = 1;

				if(graphe->oriente == 0)
					matrice[indexJ][indexI] = matrice[indexI][indexJ];
			}
		}
	}
	graphe->matriceAdj = matrice;
}

graph_genetic_t* createGraphFromTSP(char* nomFichier)
{
	FILE* fichier;
	fichier = fopen(nomFichier, "r");
	char* line = NULL;
    size_t len = 0;
    int stop = 0;
    char* token;
	graph_genetic_t* graphe;

    graphe = (graph_genetic_t*) malloc(sizeof(graph_genetic_t));

	graphe->oriente = 0;
	graphe->value = 1;
	graphe->complet = 1;

	LOG("Loading problem...");

	while ((getline(&line, &len, fichier)) != -1 && stop == 0)
	{
		token = strtok(line, " :");

		if(strcmp(token, "NAME") == 0)
		{
			token = strtok(NULL, ": \n");
			LOG("Name : %s", token);
		}
		else if(strcmp(token, "COMMENT") == 0)
		{
			token = strtok(NULL, ":\n");
			LOG("Description :%s", token);
		}
		else if(strcmp(token, "TYPE") == 0)
		{
			token = strtok(NULL, ": \n");
			LOG("Type : %s", token);
		}
		else if(strcmp(token, "DIMENSION") == 0)
		{
			token = strtok(NULL, ": \n");
			graphe->nSommets = atoi(token);
			LOG("Dimension : %d", atoi(token));
		}
		else if(strcmp(token, "EDGE_WEIGHT_TYPE") == 0)
		{
			token = strtok(NULL, ": \n");
			if(strcmp(token, "EUC_2D") == 0)
			{
				LOG("Edge weight type : %s", token);
				fillFromEuclideanDistances(graphe, fichier, line);
			}
			else
			{
				ERROR("Unknown edge weight type ! Stopping...");
				fclose(fichier);
				exit(2);
			}
		}
    }
    fclose(fichier);
	free(line);
    return graphe;
}

void fillFromEuclideanDistances(graph_genetic_t* graphe, FILE* file, char* line)
{
	int i, j;
	double *x, *y, xd, yd;
	char* token;
	size_t len = 0;
	int** matrice = createAdjMatrix(graphe->nSommets);



	// tabArretes = (arrete_t**)malloc((graphe->nSommets - 1) * (graphe->nSommets - 1) * sizeof(arrete_t));
	x = (double*) malloc(graphe->nSommets * sizeof(double));
	y = (double*) malloc(graphe->nSommets * sizeof(double));

	/* Skipping "NODE_COORD_SECTION" */
	if(getline(&line, &len, file) == -1);

	while (getline(&line, &len, file) != -1 && strcmp(line, "EOF\n") != 0)
	{
		/* List */
		token = strtok(line, " ");
		i = atoi(token) - 1;

		token = strtok(NULL, " ");
		x[i] = atof(token);

		token = strtok(NULL, " ");
		y[i] = atof(token);
	}

	// #pragma omp parallel for collapse(2)
	for (i = 0; i < graphe->nSommets; i++)
	{
		for (j = 0; j < graphe->nSommets; j++)
		{
			xd = x[i] - x[j];
			yd = y[i] - y[j];
			matrice[i][j] = (int)(sqrt(xd*xd + yd*yd)+0.5f);
		}
	}

	graphe->matriceAdj = matrice;

	free(x);
	free(y);
	free(line);
}
