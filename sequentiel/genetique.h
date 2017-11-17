#ifndef GENETIQUE__H
#define GENETIQUE__H

#include "graphGenetic.h"

typedef struct __geneticAlgorithm {
	int nPersons;
	int nGenerations;
	int nParents;
	int mutationRate;
} Genetic;

typedef struct __person {
	int townCount;
	int fitnessValue;
	int* hamiltonianWay;
} Person;

typedef struct __population {
	int size;
	Person* persons;
} Population;

Genetic* configureAlgorithm(char* fileName);
Population* populate(int popSize, int nEdges, int** matriceAdj);
void createPerson(Person* buffer, int townCount);
void copyPerson(Person* p, Person* buffer, int cutPoint);
void copyCutPerson(Person* p, Person* buffer, int cutPoint);
void setRandomHamiltonianWay(Person* p);
void showPerson(Person person);
void showPopulation(Population population);
int evaluate(Person* person, int** matriceAdj);
void selection(Population* population, int nParents, int* buffer);
void reproduce(Person* p1, Person* p2, Person* firstChild, Person* secondChild);

void freePerson(Person* p);
void freePopulation(Population* p);

int isUnique(int* t, int size, int index);

void createNewPopulation(Population* population, int* selectedParents, int nParents, int** matriceAdj, int mutationRate);
int* getWorstPersons(Population* population, int nPersons);
int* getBestPersons(Person* p, int size, int nPersons);
int* sortPersonsByWorst(Population* population);

void mutate(Person* person, int** matriceAdj);

#endif
