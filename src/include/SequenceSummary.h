#ifndef SequenceSummary_H
#define SequenceSummary_H


#include <string>
#include <map>
#include <algorithm>
#include <cctype>
#include <vector>
#include <array>
#include <iostream>
#include "CodonTable.h"

#ifndef STANDALONE
#include <Rcpp.h>
#endif

class SequenceSummary
{
	private:

		unsigned ncodons[64];
		unsigned RFPObserved[64];
		unsigned naa[26];
		std::vector <std::vector <unsigned>> codonPositions;

	public:

		//Constructors & destructors:
		explicit SequenceSummary();
		SequenceSummary(const std::string& sequence);
		SequenceSummary(const SequenceSummary& other);
		SequenceSummary& operator=(const SequenceSummary& other);
		virtual ~SequenceSummary(); //TODO:Why is this virtual????



		//Data Manipulation Functions:
		unsigned getAACountForAA(std::string aa);
		unsigned getAACountForAA(unsigned aaIndex);
		unsigned getCodonCountForCodon(std::string& codon);
		unsigned getCodonCountForCodon(unsigned codonIndex);
		unsigned getRFPObserved(std::string codon);
		unsigned getRFPObserved(unsigned codonIndex);
		void setRFPObserved(unsigned codonIndex, unsigned value);
		std::vector <unsigned> *getCodonPositions(std::string codon);
		std::vector <unsigned> *getCodonPositions(unsigned index);



		//Other Functions:
		void clear();
		bool processSequence(const std::string& sequence);  //TODO: WHY return a bool


		//Static functions:
		static char complimentNucleotide(char ch);



		//R Section:

#ifndef STANDALONE

		//Data Manipulation Functions:
		unsigned getAACountForAAR(std::string aa);
		unsigned getAACountForAAIndexR(unsigned aaIndex); //TEST THAT ONLY!
		unsigned getCodonCountForCodonR(std::string& codon);
		unsigned getCodonCountForCodonIndexR(unsigned codonIndex); //TEST THAT ONLY!
		unsigned getRFPObservedForCodonR(std::string codon);
		unsigned getRFPObservedForCodonIndexR(unsigned codonIndex); //TEST THAT ONLY!
		std::vector <unsigned> getCodonPositionsForCodonR(std::string codon);
		std::vector <unsigned> getCodonPositionsForCodonIndexR(unsigned codonIndex); //TEST THAT ONLY!

#endif //STANDALONE


	protected:
};

#endif // SequenceSummary_H

/*--------------------------------------------------------------------------------------------------
 *                                   !!!RCPP NOTE!!!
 * All functions are exposed to R. However, some functions are only exposed for the purpose of
<<<<<<< HEAD:ribModel/src/include/SequenceSummary.h
 * unit testing (Test That). These functions can be accessed in R, but do not check indices or
 * potential problems with user input (such as codon strings being lower case). The two functions
 * that have R wrappers that only call the C++ function had to be implemented because RCPP cannot
 * deal with overloaded functions.
 -------------------------------------------------------------------------------------------------*/
