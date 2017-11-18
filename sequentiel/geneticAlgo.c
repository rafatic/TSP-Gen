#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <limits.h>
#include <errno.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>
#include <time.h>

#include "debug.h"
#include "genetique.h"
#include "graphGenetic.h"
#include "graphTools.h"





int main(int argc, char *argv[])
{
	int *selectedParents = NULL, i = 0, j, bestFitnessValue = INT_MAX;
	Population* population;
	graph_genetic_t* graph = NULL;
	Person* bestPerson;
	char* ext;
	Genetic* genetic;
	char fileName[255];
	clock_t startTime, endTime;

	srand(time(NULL) + getpid());
	

    if (argc < 3)
	{
		//ERROR("Invalid arguments. The command must be under the following form : commandName <program> <graph_file_path>.txt|.tsp <population_size> <generation_count> <parents_count> <mutation_rate>");
		ERROR("Invalid arguments. The command must be under the following form : commandName <program> <graph_file_path> (*.txt|*.tsp <configuration file>");
		exit(1);
	}
	ext = strrchr(argv[1], '.');
	LOG("%s file detected.", ext);

	LOG("%s configuration detected.", argv[2]);
	strcpy(fileName, argv[2]);

	genetic = configureAlgorithm(fileName);
	/*printf("nGenerations : %d\n", genetic->nGenerations);
	printf("nPersons : %d\n", genetic->nPersons);
	printf("nParents : %d\n", genetic->nParents);
	printf("mutationRate : %d\n", genetic->mutationRate);*/

	

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

	startTime = clock();

	
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
	
	endTime = clock();

	//printf("Best result :\n");
	//showPerson(*bestPerson);

	//printf("Exec time : %ldms\n", endTime - startTime);

	//csv format
	printf("%d, %ld\n", bestPerson->fitnessValue, endTime - startTime);

	freePopulation(population);
	free(selectedParents);
	

	exit(0);
}
