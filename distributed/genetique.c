#define _GNU_SOURCE
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <string.h>
#include <mpi.h>
#include "debug.h"
#include "genetique.h"

#define MSG_MIGRANT 3
#define MSG_MIGRANT_WAY 4

Genetic* configureAlgorithm(char* fileName)
{
	FILE* configFile;
	configFile = fopen(fileName, "r");
	char* line = NULL;
	size_t len = 0;
	int stop = 0;
	char* token;
	Genetic* configuration;

	configuration = (Genetic*)malloc(sizeof(Genetic));

	while((getline(&line, &len, configFile)) != - 1 && stop == 0)
	{
		token = strtok(line, " ");
		if(strcmp(token, "nGenerations") == 0)
		{
			token = strtok(NULL, " ");
			configuration->nGenerations = atoi(token);
		}
		else if(strcmp(token, "nPersons") == 0)
		{
			token = strtok(NULL, " ");
			configuration->nPersons = atoi(token);	
		}
		else if(strcmp(token, "nParents") == 0)
		{
			token = strtok(NULL, " ");
			configuration->nParents = atoi(token);	
			
		}
		else if(strcmp(token, "mutationRate") == 0)
		{
			token = strtok(NULL, " ");
			configuration->mutationRate = atoi(token);
			if(configuration->mutationRate > 100)
			{
				configuration->mutationRate = 100;
			}
		}
		else if(strcmp(token, "migrationPeriodicity") == 0)
		{
			token = strtok(NULL, " ");
			configuration->migrationPeriodicity = atoi(token);
			if(configuration->migrationPeriodicity < 0)
			{
				configuration->migrationPeriodicity = 0;
			}
		}
		else if(strcmp(token, "nMigrants") == 0)
		{
			token = strtok(NULL, " ");
			configuration->nMigrants = atoi(token);
			if(configuration->nMigrants < 0)
			{
				configuration->nMigrants = 0;
			}
			fclose(configFile);
			stop = 1;
		}
		else
		{
			ERROR("Unexpected token : '%s'\n", token);
			ERROR("Please check your configuration file. Aborting execution\n");
			fclose(configFile);
			stop = 1;
			return NULL;
		}
	}

	if(configuration->nParents * 2 > configuration->nPersons)
	{
		WARNING("Warning : Invalid parameter : nParents = %d", configuration->nParents);
		WARNING("Parents should be less than the half of the entire population");
		WARNING("Setting nParents to nPersons / 2");
		configuration->nParents = configuration->nPersons / 2;
	}

	if(configuration->migrationPeriodicity > configuration->nGenerations)
	{
		WARNING("WARNING : Number of generation must be greater than the migration periodicity");
		WARNING("Please check your configuration file :");
		WARNING("Migration parameters will be ignored");
		configuration->migrationPeriodicity = 0;
		configuration->nMigrants = 0;
	}
	if(configuration->nMigrants > configuration->nPersons)
	{
		WARNING("WARNING : Population size must be greater than the number of migrants");
		WARNING("Please check your configuration file !");
		WARNING("Migration parameters will be ignored");
		configuration->migrationPeriodicity = 0;
		configuration->nMigrants = 0;
	}
	return configuration;
}

void createPerson(Person* buffer, int townCount)
{
	buffer->townCount = townCount;
	buffer->hamiltonianWay = (int*) calloc(townCount, sizeof(int));
	buffer->fitnessValue = 0;
}

void copyPerson(Person* p, Person* buffer, int cutPoint)
{
	int i;

	buffer->townCount = p->townCount;
	for (i = 0; i < p->townCount; i++)
		buffer->hamiltonianWay[i] = p->hamiltonianWay[i];
	buffer->fitnessValue = p->fitnessValue;
}

void setRandomHamiltonianWay(Person* p)
{
	int i, k, tmp;


	/* Fill the way array with 0..n-1 */
	
    for (i = 0; i < p->townCount; i++)
	{
		p->hamiltonianWay[i] = i;
	}

	/* Then shuffle the array (Fisher-Yates method) */
	for (i = p->townCount - 1; i > 1; i--)
	{
		k = rand() % i;

		tmp = p->hamiltonianWay[i];
		p->hamiltonianWay[i] = p->hamiltonianWay[k];
		p->hamiltonianWay[k] = tmp;
	}

}
/**
 * Function populate
 * Creates a population, giving their persons a random way
 * @param int	popSize	The number of persons in the population
 * @param int	nEdges	The number of edges in the graph
 * @return	The pointer to the created population
 */
Population* populate(int popSize, int nEdges, int** matriceAdj)
{
	int i;
	Population* population;

	/* Allocating the population structure */
	population = (Population*) malloc(sizeof(Population));

	population->size = popSize;
	population->persons = (Person*) calloc(popSize, sizeof(Person));

	/* For each people, we're putting a random hamiltonian way */
	for (i = 0; i < popSize; i++)
	{
		// population->persons[i].townCount = nEdges;
		// population->persons[i].hamiltonianWay = (int*) malloc(nEdges * sizeof(int));

		createPerson(&population->persons[i], nEdges);

		setRandomHamiltonianWay(&population->persons[i]);

		population->persons[i].fitnessValue = evaluate(&population->persons[i], matriceAdj);
	}

	SUCCESS("Succesfully generated a population with %d persons.", population->size);

	if (DEBUG == 1)
		showPopulation(*population);

	return population;
}

void showPerson(Person person)
{
	int i;

	printf("\tWay : ");
	for (i = 0; i < person.townCount; i++)
	{
		printf("%d ", person.hamiltonianWay[i]);
	}
	printf("\tScore : %d\n", person.fitnessValue);
}

/**
 * Function showPopulation
 * Display the whole population
 * @param	Population	population	A copy of the population structure
 * @return	NULL
 */
void showPopulation(Population population)
{
	int i;

	PRINT("Displaying the %d-persons population :", population.size);

	for (i = 0; i < population.size; i++)
	{
		PRINT("Person %d:", i);
		showPerson(population.persons[i]);
	}
}

/**
 * Function evaluate
 *
 * @param	Person*		person
 * @param	int**		matriceAdj
 * @return	int			The score
 */
int evaluate(Person* person, int** matriceAdj)
{
	int i, score = 0;

	for (i = 0; i < person->townCount - 1; i++)
	{
		score += matriceAdj[person->hamiltonianWay[i]][person->hamiltonianWay[i+1]];
	}
	score += matriceAdj[person->hamiltonianWay[person->townCount - 1]][person->hamiltonianWay[0]];

	return score;
}

void selection(Population* population, int nParents, int* buffer)
{
	int i, k;
	double *probability, reverseScoreTotal = 0.0, random;




	/* Here, we're working with indexes, each selected person will get it index put in the array */
	probability = (double*) malloc(population->size * sizeof(double));


	/* Computing 1/score */
	for (i = 0; i < population->size; i++)
	{
		probability[i] = 1.0 / population->persons[i].fitnessValue;
		reverseScoreTotal += probability[i];
	}


	/* Dividing each one by s (reverseScoreTotal) */
	for (i = 0; i < population->size; i++)
	{
		probability[i] /= reverseScoreTotal;
	}

	/* Random number between 0 and 1 */
	random = (double)rand() / (double)RAND_MAX;
	if(DEBUG == 1)
	{
		LOG("Randomly-selected number : %lf", random);
	}
	

	/* Doing partial sum */
	for (i = 1; i < population->size; i++)
	{
		probability[i] += probability[i-1];
		if(DEBUG == 1)
			LOG("Partial sum for person %d : %lf", i, probability[i]);
	}

	/* Deciding */
	for (i = 0; i < nParents; i++)
	{
		k = 0;

		while(k < population->size - 1 && probability[k] < random)
		{
			k++;
		}
		/* If i is already selected, we have to roll back until another person */
		while(probability[k] == -1.0 && i > 0)
		{
			k--;
		}
		buffer[i] = k;

		probability[k] = -1.0;
	}


}

void copyCutPerson(Person* p, Person* buffer, int cutPoint)
{
	int i;

	buffer->townCount = p->townCount;
	for (i = 0; i < cutPoint; i++)
		buffer->hamiltonianWay[i] = p->hamiltonianWay[i];
	buffer->fitnessValue = p->fitnessValue;
}

void reproduce(Person* p1, Person* p2, Person* firstChild, Person* secondChild)
{
	int i, j, currentIndex, cutPoint;


	cutPoint = p1->townCount / 2;

	copyCutPerson(p1, firstChild, cutPoint);

	i = cutPoint;

	for (j = 0; j < firstChild->townCount; j++)
	{

		currentIndex = p2->hamiltonianWay[(j + cutPoint) % p2->townCount];

		if(isUnique(firstChild->hamiltonianWay, cutPoint, currentIndex))
		{
			firstChild->hamiltonianWay[i] = currentIndex;
			i++;
		}
	}
	i = cutPoint;
	copyCutPerson(p2, secondChild, cutPoint);

	for (j = 0; j < secondChild->townCount; j++)
	{
		currentIndex = p1->hamiltonianWay[(j + cutPoint) % p1->townCount];

		if(isUnique(secondChild->hamiltonianWay, i, currentIndex))
		{
			secondChild->hamiltonianWay[i] = currentIndex;
			i++;
		}
	}




}

void freePerson(Person* p)
{
	free(p->hamiltonianWay);
	free(p);
}

void freePopulation(Population* p)
{
	int i;
	for (i = 0; i < p->size; i++)
		free(p->persons[i].hamiltonianWay);
	free(p->persons);
	free(p);
}

int isUnique(int* t, int size, int index)
{
	int i = 0, unique = 1;

	while (i < size && unique == 1)
	{
		if (t[i] == index)
			unique = 0;
		else
			i++;
	}

	return unique;
}

void createNewPopulation(Population* population, int* selectedParents, int** matriceAdj, Genetic* configuration, int generationNumber)
{
	
	int i, nChildrenToAdd = configuration->nParents, mutationRoll, numberOfMutations, rank, commSize;
	int *worstPersons, *bestPersons;
	Person* children = NULL;
	Person* migrants;
	MPI_Status status;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &commSize);

	// TEMPORARY, NEED TO BE CHECKED DURING INITIALIZATION
	if(configuration->nParents % 2 != 0)
	{
		configuration->nParents--;
	}

	children = (Person*) malloc(nChildrenToAdd * sizeof(Person));
	for (i = 0; i < nChildrenToAdd; i++)
	{
		createPerson(&children[i], population->persons[0].townCount);
	}

	worstPersons = getWorstPersons(population, nChildrenToAdd);

	for(i = 0; i < configuration->nParents; i = i + 2)
	{
		reproduce(&population->persons[selectedParents[i]], &population->persons[selectedParents[i+1]], &children[i], &children[i+1]);

		population->persons[worstPersons[i]] = children[i];
		population->persons[worstPersons[i+1]] = children[i+1];
		
	}

	


	mutationRoll = rand() % 100;
	if(mutationRoll < configuration->mutationRate)
	{
		numberOfMutations = rand() % population->size;
		
		for(i=0; i < numberOfMutations; i++)
		{
			mutate(&population->persons[rand()% (population->size - 1)], matriceAdj);
		}
	}

	if(configuration->migrationPeriodicity > 0)
	{
		if(generationNumber % configuration->migrationPeriodicity == 0)
		{
			bestPersons = getBestPersons(population, configuration->nMigrants);
			
			/*printf("%d : SENDING PERSONS\n", rank);
			for(i = 0; i < configuration->nMigrants; i++)
			{
				showPerson(population->persons[bestPersons[i]]);
			}*/

			
			
			
			//printf("rank : %d\n", rank);

			for(i = 0; i < configuration->nMigrants; i++)
			{
				if(rank == commSize - 1)
				{
					//printf("%d sending to : %d\n", rank, 0);
					MPI_Send(&population->persons[bestPersons[i]], 1, MPI_PersonType, 0, MSG_MIGRANT, MPI_COMM_WORLD);
					MPI_Send(population->persons[bestPersons[i]].hamiltonianWay, population->persons[bestPersons[i]].townCount, MPI_INT, 0, MSG_MIGRANT_WAY, MPI_COMM_WORLD);
				}
				else
				{
					//printf("%d sending to : %d\n", rank, rank + 1);
					MPI_Send(&population->persons[bestPersons[i]], 1, MPI_PersonType, rank + 1, MSG_MIGRANT, MPI_COMM_WORLD);
					MPI_Send(population->persons[bestPersons[i]].hamiltonianWay, population->persons[bestPersons[i]].townCount, MPI_INT, rank + 1, MSG_MIGRANT_WAY, MPI_COMM_WORLD);
				}

				
			}

			worstPersons = getWorstPersons(population, configuration->nMigrants);
			migrants = (Person*)malloc(configuration->nMigrants * sizeof(Person));
			for(i = 0; i < configuration->nMigrants; i++)
			{
				migrants[i].hamiltonianWay = (int*)malloc(population->persons[0].townCount * sizeof(int));
			}
			/*printf("%d : REPLACING PERSONS\n", rank);
			for(i = 0; i < configuration->nMigrants; i++)
			{
				showPerson(population->persons[worstPersons[i]]);
			}*/
			for(i = 0; i < configuration->nMigrants; i++)
			{
				// if(rank == 0)
				// {

				// 	printf("%d receiving from %d\n", rank, commSize - 1);
				// 	MPI_Recv(&population->persons[worstPersons[i]], 1, MPI_PersonType, commSize - 1, MSG_MIGRANT, MPI_COMM_WORLD, &status);

				// 	//printf("size : %ld\n", sizeof(population->persons[worstPersons[i]].hamiltonianWay));
				// 	MPI_Recv(population->persons[worstPersons[i]].hamiltonianWay, population->persons[worstPersons[i]].townCount, MPI_INT, commSize - 1, MSG_MIGRANT_WAY, MPI_COMM_WORLD, &status);
				// }
				// else
				// {
				// 	printf("%d receiving from %d\n", rank, rank - 1 );
				// 	MPI_Recv(&population->persons[worstPersons[i]], 1, MPI_PersonType, rank - 1, MSG_MIGRANT, MPI_COMM_WORLD, &status);	

				// 	//printf("size : %ld\n", sizeof(population->persons[worstPersons[i]].hamiltonianWay));
				// 	MPI_Recv(population->persons[worstPersons[i]].hamiltonianWay, population->persons[worstPersons[i]].townCount, MPI_INT, rank - 1, MSG_MIGRANT_WAY, MPI_COMM_WORLD, &status);
				// }
				if(rank == 0)
				{
					MPI_Recv(&migrants[i], 1, MPI_PersonType, commSize - 1, MSG_MIGRANT, MPI_COMM_WORLD, &status);

					//printf("size : %ld\n", sizeof(population->persons[worstPersons[i]].hamiltonianWay));
					MPI_Recv(migrants[i].hamiltonianWay, migrants[i].townCount, MPI_INT, commSize - 1, MSG_MIGRANT_WAY, MPI_COMM_WORLD, &status);
				}
				else
				{
					MPI_Recv(&migrants[i], 1, MPI_PersonType, rank - 1, MSG_MIGRANT, MPI_COMM_WORLD, &status);	

					//printf("size : %ld\n", sizeof(population->persons[worstPersons[i]].hamiltonianWay));
					MPI_Recv(migrants[i].hamiltonianWay, migrants[i].townCount, MPI_INT, rank - 1, MSG_MIGRANT_WAY, MPI_COMM_WORLD, &status);
				}
			}

			for(i = 0; i < configuration->nMigrants; i++)
			{
				memcpy(&population->persons[worstPersons[i]], &migrants[i], sizeof(Person));
			}
		}
	}





}


int* getWorstPersons(Population* population, int nPersons)
{
	int* worstPersons = (int*)malloc(sizeof(int) * nPersons);
	int i, j, min, indexMin= 0;

	min = INT_MAX;
	for (i = 0; i < nPersons; i++)
	{
		worstPersons[i] = i;
		if(population->persons[worstPersons[i]].fitnessValue < min)
		{
			min = population->persons[worstPersons[i]].fitnessValue;
			indexMin = i;
		}		
	}

	for (i = nPersons; i < population->size; i++)
	{
		if(population->persons[i].fitnessValue > min)
		{
			worstPersons[indexMin] = i;
			min = INT_MAX;
			for(j = 0; j < nPersons; j++)
			{
				if(population->persons[worstPersons[j]].fitnessValue < min)
				{
					min = population->persons[worstPersons[j]].fitnessValue;
					indexMin = j;
				}
			}	
		}
	}

	return worstPersons;
}

int* getBestPersons(Population* population, int nPersons)
{
	int* bestPersons = (int*) malloc(sizeof(int)* nPersons);
	int i, j, max = 0, indexMax = 0;

	for (i = 0; i < population->size; i++)
	{
		bestPersons[i] = i;
		if(population->persons[bestPersons[i]].fitnessValue > max)
		{
			max = population->persons[bestPersons[i]].fitnessValue;
			indexMax = i;
		}
	}

	for (i = population->size; i < nPersons; i++)
	{
		if (population->persons[i].fitnessValue < max)
		{
			bestPersons[indexMax] = i;
			max = 0;
			for(j = 0; j < population->size; j++)
			{
				if (population->persons[bestPersons[j]].fitnessValue > max)
				{
					max = population->persons[bestPersons[j]].fitnessValue;
					indexMax = j;
				}
			}
		}
	}
	return bestPersons;

}

/*int* getBestPersons(Person* p, int size, int nPersons)
{
	int* bestPersons = (int*) malloc(sizeof(int)* nPersons);
	int i, j, max, indexMax = 0;

	max = 0;

	for (i = 0; i < nPersons; i++)
	{
		bestPersons[i] = i;
		if(p[bestPersons[i]].fitnessValue > max)
		{
			max = p[bestPersons[i]].fitnessValue;
			indexMax = i;
		}
	}

	for (i = nPersons; i < size; i++)
	{
		if (p[i].fitnessValue < max)
		{
			bestPersons[indexMax] = i;
			max = 0;
			for(j = 0; j < nPersons; j++)
			{
				if (p[bestPersons[j]].fitnessValue > max)
				{
					max = p[bestPersons[j]].fitnessValue;
					indexMax = j;
				}
			}
		}
	}

	return bestPersons;
}*/

void mutate(Person* person, int** matriceAdj)
{
	int indexA, indexB, tmp;

	indexA = rand() % (person->townCount - 1);
	indexB = rand() % (person->townCount - 1);


	tmp = person->hamiltonianWay[indexA];
	person->hamiltonianWay[indexA] = person->hamiltonianWay[indexB];
	person->hamiltonianWay[indexB] = tmp;

	person->fitnessValue = evaluate(person, matriceAdj);
}

void defineMPI_TypePerson(MPI_Datatype* newType)
{
	Person tmp[2];
	MPI_Aint extent = &tmp[1] - &tmp[0];

	MPI_Type_create_resized(MPI_2INT, 0, extent, newType);
	MPI_Type_commit(newType);
}
