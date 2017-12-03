#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <limits.h>
#include <errno.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>
#include <time.h>
#include <mpi.h>
#include "debug.h"
#include "genetique.h"
#include "graphGenetic.h"
#include "graphTools.h"

#define ROOT 0
#define MSG_BESTPERSON 1
#define MSG_BESTPERSON_WAY 2


int main(int argc, char *argv[])
{
	int *selectedParents = NULL, i = 0, j, bestFitnessValue = INT_MAX, rank, p, indexBestPerson;
	Population* population;
	graph_genetic_t* graph = NULL;
	Person *bestPerson, *reducedPersons;
	char* ext;
	Genetic* genetic;
	char fileName[255];
	clock_t startTime, endTime;
	MPI_Status status;
	

	

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
	if(genetic == NULL)
	{
		exit(1);
	}
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

	startTime = clock();

	MPI_Init(NULL, NULL);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	defineMPI_TypePerson(&MPI_PersonType);
	if(rank == ROOT)
	{
		reducedPersons = (Person*)malloc(p * sizeof(Person));
		for(i = 0; i < p; i++)
		{
			reducedPersons[i].hamiltonianWay = (int*)malloc(graph->nSommets * sizeof(int));
		}
		

	}

	srand(time(NULL) + getpid() + rank);
	population = populate(genetic->nPersons, graph->nSommets, graph->matriceAdj);


	
	for(i=0; i < genetic->nGenerations; i ++)
	{
		selection(population, genetic->nParents, selectedParents);

		createNewPopulation(population, selectedParents, graph->matriceAdj, genetic, i);

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
	//printf("%d, %d, %ld\n", rank, bestPerson->fitnessValue, endTime - startTime);

	if(rank != ROOT)
	{
		MPI_Send(bestPerson, 1, MPI_PersonType, ROOT, MSG_BESTPERSON, MPI_COMM_WORLD);
		MPI_Send(bestPerson->hamiltonianWay, bestPerson->townCount, MPI_INT, 0, MSG_BESTPERSON_WAY, MPI_COMM_WORLD);
	}
	else
	{
		reducedPersons[0] = *bestPerson;
		if(p > 1)
		{
			for(i = 1; i < p; i++)
			{
				MPI_Recv(&reducedPersons[i], 1, MPI_PersonType, i, MSG_BESTPERSON, MPI_COMM_WORLD, &status);
				MPI_Recv(reducedPersons[i].hamiltonianWay, reducedPersons[i].townCount, MPI_INT, i, MSG_BESTPERSON_WAY, MPI_COMM_WORLD, &status);
			}	
		}

		indexBestPerson = 0;

		for(i = 1; i < p; i++)
		{
			if(reducedPersons[i].fitnessValue < reducedPersons[indexBestPerson].fitnessValue)
			{
				indexBestPerson = i;
			}
		}

		printf("%d, %d\n", indexBestPerson, reducedPersons[indexBestPerson].fitnessValue);

		/*for(i = 0; i < p; i++)
		{
			printf("%d, %d, ", i, reducedPersons[i].fitnessValue);

			for (j = 0; j < bestPerson->townCount; j++)
			{
				printf("%d ", reducedPersons[i].hamiltonianWay[j]);
			}
			printf("\n");
		}*/
		
	}

	

	freePopulation(population);
	free(selectedParents);

	MPI_Type_free(&MPI_PersonType);
	MPI_Finalize();


	exit(0);
}
