#ifndef GENOME_H
#define GENOME_H

#include <vector>
#include <string>
#include <map>
#include <fstream>
#include <sstream>
#include <cmath>

#include "Gene.h"
#include "CodonTable.h"

#ifndef STANDALONE
#include <Rcpp.h>
#endif

#include "Gene.h"

class Model;
class Genome
{
	private:

		std::vector<Gene> genes;
		std::vector<Gene> simulatedGenes;
		std::vector <unsigned> numGenesWithPhi;

		static std::vector <std::string> defaultVector;

	public:

		//Constructors & destructors:
		explicit Genome();
		explicit Genome(unsigned codonTableId, std::string model, bool splitAA);
		explicit Genome(unsigned codonTableId, std::string model, bool splitAA, std::vector <std::string> groupList);
		Genome& operator=(const Genome& other);
		bool operator==(const Genome& other) const;
		virtual ~Genome();


		//File I/O Functions:
		void readFasta(std::string filename, bool Append = false);
		void writeFasta(std::string filename, bool simulated = false);
		void readRFPFile(std::string filename);
		void writeRFPFile(std::string filename, bool simulated = false);
		void readObservedPhiValues(std::string filename, bool byId = true);


		//Gene Functions:
		void addGene(const Gene& gene, bool simulated = false);
		std::vector <Gene> getGenes(bool simulated = false);
		unsigned getNumGenesWithPhi(unsigned index);
		Gene& getGene(unsigned index, bool simulated = false);
		Gene& getGene(std::string id, bool simulated = false);


		//Other Functions:
		unsigned getGenomeSize();
		void clear();
		Genome getGenomeForGeneIndicies(std::vector <unsigned> indicies, bool simulated = false); //NOTE: If simulated is true, it will return a genome with the simulated genes, but the returned genome's genes vector will contain the simulated genes.
		std::vector<unsigned> getCodonCountsPerGene(std::string codon);



		//R Section:

#ifndef STANDALONE


		bool checkIndex(unsigned index, unsigned lowerbound, unsigned upperbound);
		Gene& getGeneByIndex(unsigned index, bool simulated = false);
		Gene& getGeneById(std::string ID, bool simulated = false);
		Genome getGenomeForGeneIndiciesR(std::vector <unsigned> indicies, bool simulated = false);
        std::vector <std::string> getGroupListFromGenomeR();
        std::vector <std::string> AAToCodonFromGenomeR(std::string, bool withoutReference);

#endif //STANDALONE

	protected:
};

#endif // GENOME_H
