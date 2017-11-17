#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <limits.h>
#include <string.h>
#include "debug.h"
#include "genetique.h"
#include "graphGenetic.h"
#include "graphTools.h"

#define NGENERATIONS 5
#define NPERSON	10
#define NPARENTS 4
#define MUTATION_RATE 40

int main(int argc, char *argv[])
{
	int *selectedParents = NULL, i = 0, j, bestFitnessValue = INT_MAX;
	Population* population;
	graph_genetic_t* graph = NULL;
	Person* bestPerson;
	char* ext;
	Genetic* genetic;
	char fileName[80];

	srand(time(NULL));
	

    if (argc < 2)
	{
		//ERROR("Invalid arguments. The command must be under the following form : commandName <program> <graph_file_path>.txt|.tsp <population_size> <generation_count> <parents_count> <mutation_rate>");
		ERROR("Invalid arguments. The command must be under the following form : commandName <program> <graph_file_path> (*.txt|*.tsp");
		exit(1);
	}
	ext = strrchr(argv[1], '.');
	LOG("%s file detected.", ext);

	strcpy(fileName, "param.cfg");

	genetic = configureAlgorithm(fileName);
	printf("nGenerations : %d\n", genetic->nGenerations);
	printf("nPersons : %d\n", genetic->nPersons);
	printf("nParents : %d\n", genetic->nParents);
	printf("mutationRate : %d\n", genetic->mutationRate);

	

	if (strcmp(ext, ".txt") == 0)
	{
    	graph = createGraph(argv[1]);
	}	
	else
	{
		graph = createGraphFromTSP(argv[1]);
	}

	if (DEBUG)
		showMatrix(graph);

	
	selectedParents = (int*) malloc(genetic->nParents * sizeof(int));
	population = populate(genetic->nPersons, graph->nSommets, graph->matriceAdj);

	
	for(i=0; i < genetic->nGenerations; i ++)
	{
		selection(population, genetic->nParents, selectedParents);

		createNewPopulation(population, selectedParents, genetic->nParents, graph->matriceAdj, genetic->mutationRate);

		for (j = 0; j < population->size; j++)
		{
			if (population->persons[j].fitnessValue < bestFitnessValue)
			{
				bestPerson = &population->persons[j];
				bestFitnessValue = bestPerson->fitnessValue;
			}
		}
	}
	
	

	printf("Best result :\n");
	showPerson(*bestPerson);

	//freePopulation(population);
	//free(selectedParents);

	exit(0);
}
