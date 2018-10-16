#ifndef BRUTE_FORCE
#define BRUTE_FORCE

#include<stdio.h>
#include<iostream>
#include <map>
#include <vector>

typedef std::map<std::string, int> annotationCounter;

struct ThreePrimeModification {
 
  int max_nts = 5;
  int max_nts_losses = 5;
  
};

struct FivePrimeModification {
 
  int max_nts = 3;
  int max_nts_losses = 3;

};

struct SNPmodification {
 
  int max_nts = 1;

};

struct InDelModification {
 
  int max_nts = 1;

};

void bruteForce(std::string input, std::ofstream& output, std::string file_type, int option, int times);

std::vector<std::pair<std::string, std::string>> sequenceModification(int option, std::string sequence, std::string name, std::string mimat);

void fivePrimeAddition(std::string sequence, std::string alphabet);

void threePrimeAddition(std::string sequence, std::string alphabet);

void SNPs(std::string sequence, int max_nts, int position);
 
void indel(std::string sequence, int max_nts, int position);
 
void threePrimeLoss(std::string input_sequence, std::string sequence, int max_nts);
 
void fivePrimeLoss(std::string input_sequence, std::string sequence, int max_nts);

void storeSequences(std::vector<std::pair<std::string, std::string>> modified_sequences, 
		    std::vector<std::pair<std::string, std::string>> & sequences);
  
void writeSequences(std::vector<std::pair<std::string, std::string>> sequences, std::ofstream& output_file, 
		    std::string read_id_token, int times, std::map<std::string, int> &annotation_counter);
 
#endif
