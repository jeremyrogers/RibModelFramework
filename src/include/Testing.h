#ifndef Testing_H
#define Testing_H


#include "SequenceSummary.h"
#include "Gene.h"
#include "Genome.h"


void testSequenceSummary();
void testGene();
void testGenome(std::string testFileDir);
void testCodonTable();
void initMutation(std::vector<double> mutationValues, unsigned mixtureElement, std::string aa);



//Blank header
#endif // Testing_H