#include<stdio.h>
#include<iostream>
#include<fstream>
#include <string>
#include <regex>
#include <vector> 
#include <map>
#include "project/brute_force.h"
#include "project/simulator.h"
#include "project/utils.h"

struct ThreePrimeModification three_prime_modification;

struct FivePrimeModification five_prime_modification;

struct SNPmodification snp_modification;

struct InDelModification indel_modification;

std::vector<std::pair<std::string, std::string>> modified_sequences;

std::string nucleotides [4] = {"A", "T", "G", "C"};

typedef std::map<std::string, int> annotationCounter;

int read_id_counter;

annotationCounter annotation_counter;

std::vector<std::pair<std::string, std::string>> sequenceModification(int option, std::string sequence, std::string name, std::string mimat){
    
  std::vector<std::pair<std::string, std::string>> sequences;
  
  std::ifstream input;
  
  int flag = 0;
  int start = 1;
  int counter = 1;
  std::string alphabet = "";
  
  switch(option){
    
    case 0: //Canonical miRNA sequence.
      
      modified_sequences.emplace_back(sequence,"nm"); //Store plain sequence without modification.
      storeSequences(modified_sequences, sequences);
      modified_sequences.clear();
      
    break;
    
    case 1: //5' templated addition.
      
        input.open(fixed_hairpin_species_spec.c_str());
 	
 	for( std::string line; getline( input, line ); )
           {
	     std::string mi;
	     for(auto const& item:mi_mimat){
	       if(item.second.compare(mimat) == 0){
		 mi = item.first;
	       }
	     }
	  if (line.find(mi) != std::string::npos) flag = 1;
	    if(counter % 2 == 0 && flag == 1){ //Sequence line
	      std::size_t found = line.find(sequence);
	      if(found!=std::string::npos){
		if(found >= five_prime_modification.max_nts)
		  alphabet = line.substr(found-five_prime_modification.max_nts, five_prime_modification.max_nts);
		else if(found != 0)
		  alphabet = line.substr(0, found);
		reverseStr(alphabet);
		fivePrimeAddition(sequence, alphabet);
		storeSequences(modified_sequences, sequences);
		modified_sequences.clear();
	      } 
	      flag = 0;
	    }
	  counter++;
	   }
    break;
    
    case 2: //3' templated addition.
          input.open(fixed_hairpin_species_spec.c_str());
 	
 	for( std::string line; getline( input, line ); )
           {
	     std::string mi;
	     for(auto const& item:mi_mimat){
	       if(item.second.compare(mimat) == 0){
		 mi = item.first;
	       }
	     }
	  if (line.find(mi) != std::string::npos) flag = 1;
	  
	  if(counter % 2 == 0 && flag == 1){ //Sequence line
	    
	    std::size_t found = line.find(sequence);
	    if(found!=std::string::npos){
	    
	    if(line.size() - found >= three_prime_modification.max_nts)
	      alphabet = line.substr(found+sequence.size(), three_prime_modification.max_nts);
	    else if(found != line.size())
	      alphabet = line.substr(found,line.size() - found);
	      threePrimeAddition(sequence, alphabet);
	      storeSequences(modified_sequences, sequences);
	      modified_sequences.clear();
	    } 
	    flag = 0;
	  }
	  
	  counter++;
	   }
    break;
    
    case 3: //A and T 5' addition.
      alphabet = "";
      for (int i = 0; i < five_prime_modification.max_nts; ++i){
	alphabet = alphabet + "A";
	
      }
      fivePrimeAddition(sequence, alphabet);
      storeSequences(modified_sequences, sequences);
      modified_sequences.clear();
      alphabet = "";
      for (int i = 0; i < five_prime_modification.max_nts; ++i){
	alphabet = alphabet + "T";
      }
      fivePrimeAddition(sequence, alphabet);
      storeSequences(modified_sequences, sequences);
      modified_sequences.clear();
    break;
    
    case 4: //A and T 3' addition.
      alphabet = "";
      for (int i = 0; i < three_prime_modification.max_nts; ++i){
	alphabet = alphabet + "A";
	
      }
      threePrimeAddition(sequence, alphabet);
      storeSequences(modified_sequences, sequences);
      modified_sequences.clear();
      alphabet = "";
      for (int i = 0; i < three_prime_modification.max_nts; ++i){
	alphabet = alphabet + "T";
      }
      threePrimeAddition(sequence, alphabet);
      storeSequences(modified_sequences, sequences);
      modified_sequences.clear();
    break;
    
    case 5: //SNP modification.
      SNPs(sequence, snp_modification.max_nts, 0); //Start at position 0 (1) of the sequence.
      storeSequences(modified_sequences, sequences);
      modified_sequences.clear();
    break;
    
    case 6: //Indel modification.
      indel(sequence, indel_modification.max_nts, 1); //Start at position 1 (2) of the sequence.
      storeSequences(modified_sequences, sequences);
      modified_sequences.clear();
    break;
     
    case 7: //5' loss/
      fivePrimeLoss(sequence, sequence, five_prime_modification.max_nts_losses); //We start with position 1 always.
      storeSequences(modified_sequences, sequences);
      modified_sequences.clear();
    break;
    
    case 8: //3' loss.
      threePrimeLoss(sequence, sequence, three_prime_modification.max_nts_losses); //We start with position 1 always.
      storeSequences(modified_sequences, sequences);
      modified_sequences.clear();
    break;
    
  };

  return sequences;
  
}

void fivePrimeAddition(std::string sequence, std::string alphabet){
  
  std::string modified_sequence = "";
  
  for (int i = 0; i < alphabet.size(); ++i){
    
    if(modified_sequence == "")
      modified_sequence = alphabet[i] + sequence; //Add nt to 5'
    else
      modified_sequence = alphabet[i] + modified_sequence; //Add nt to 5'

    std::string tag = alphabet.substr(0,i+1);
    std::replace( tag.begin(), tag.end(), 'U', 'T' ); //Replace all Us with Ts.
    reverseStr(tag); //Reverse tag.
    modified_sequences.emplace_back(modified_sequence, tag);
    
  }
  
}

void threePrimeAddition(std::string sequence, std::string alphabet){ 
  
  std::string modified_sequence = "";
  
  for (int i = 0; i < alphabet.size(); ++i){
    
    if(modified_sequence == "")
      modified_sequence = sequence + alphabet[i]; //Add nt to 5'
    else
      modified_sequence = modified_sequence + alphabet[i]; //Add nt to 5'
      
    std::string tag = alphabet.substr(0,i+1);
    std::replace( tag.begin(), tag.end(), 'U', 'T' ); //Replace all Us with Ts.
    modified_sequences.emplace_back(modified_sequence, tag);
  }
  
}

void SNPs(std::string sequence, int max_nts, int position){
  
  std::string modified_sequence;
  
      if(position < sequence.size()){
	
	for (int i = 0; i < sizeof(nucleotides)/sizeof(nucleotides[0]); ++i){
	  
	  modified_sequence = sequence;
	  std::string ms = modified_sequence.substr(position,1);
	 if( ms != nucleotides[i]){
	    modified_sequence.replace(position,1,nucleotides[i]);
	    std::string tag = std::to_string(position+1) + nucleotides[i];
	    std::replace( tag.begin(), tag.end(), 'U', 'T' ); //Replace all Us with Ts.
            modified_sequences.emplace_back(modified_sequence, tag);
	    //modified_sequences.push_back(modified_sequence); //Store sequence to the global sequence vector
// 	    if(max_nts > 1){
// 	      SNPs(modified_sequence, max_nts-1, position + 1);
// 	    }
	 }
	}
	SNPs(sequence, max_nts, position + 1);
    }
}


void indel(std::string sequence, int max_nts, int position){
  
  std::string modified_sequence;
  
      if(position < sequence.size()){
	
	for (int i = 0; i < sizeof(nucleotides)/sizeof(nucleotides[0]); ++i){
	  modified_sequence = sequence;
	  modified_sequence.insert(position, nucleotides[i]);
	  
	  //modified_sequences.push_back(modified_sequence); //Store sequence to the global sequence vector
	  
	  std::string tag = std::to_string(position+1) + nucleotides[i];
	  std::replace( tag.begin(), tag.end(), 'U', 'T' ); //Replace all Us with Ts.
          modified_sequences.emplace_back(modified_sequence, tag);
	  //if(i == 0)//We don't want any repeptitive sequences!
	  //indel(sequence, max_nts, position + 1);
	}
	indel(sequence, max_nts, position + 1);
    }
}


void threePrimeLoss(std::string input_sequence, std::string sequence, int max_nts){ 
  
    if(max_nts > 0){
       
       sequence.replace(sequence.size()-1,1,""); //Remove nt from 3'

       std::string tag = input_sequence.substr(input_sequence.size()-
       (abs(max_nts-three_prime_modification.max_nts_losses)+1),
					       (abs(max_nts-three_prime_modification.max_nts_losses)+1));
       std::replace( tag.begin(), tag.end(), 'U', 'T' ); //Replace all Us with Ts.
       modified_sequences.emplace_back(sequence, tag);
       threePrimeLoss(input_sequence, sequence, max_nts - 1); //Recursive call
    }
  
}

void fivePrimeLoss(std::string input_sequence, std::string sequence, int max_nts){
  
    if(max_nts > 0){
      
       sequence.replace(0,1,""); //Remove nt from 5'
       
       std::string tag = input_sequence.substr(0,(abs(max_nts-five_prime_modification.max_nts_losses)+1));
       std::replace( tag.begin(), tag.end(), 'U', 'T' ); //Replace all Us with Ts.
       modified_sequences.emplace_back(sequence, tag); //Store sequence to the global sequence vector
       fivePrimeLoss(input_sequence, sequence, max_nts - 1); //Recursive call
    }
}

//void storeSequences(std::vector<std::string> modified_sequences, std::vector<std::string> & sequences){
void storeSequences(std::vector<std::pair<std::string, std::string>> modified_sequences, 
		    std::vector<std::pair<std::string, std::string>> & sequences){
//   for (auto i = modified_sequences.begin(); i != modified_sequences.end(); ++i) //Read all the generated sequences obtained from equenceModification function.
// 	    {
// 	      sequences.push_back(*i);
// 	    }
    for (auto const& item:modified_sequences) //Read all the generated sequences obtained from equenceModification function.
	    {
	      sequences.emplace_back(item.first, item.second);
	    }
  
}

void writeSequences(std::vector<std::pair<std::string, std::string>> sequences, std::ofstream& output_file, 
		    std::string read_id_token, int times, std::map<std::string, int> &annotation_counter)
{
  int counter = 1;
  std::string sequence; 
  
  //for (auto i = sequences.begin(); i != sequences.end(); ++i) //Read all the generated sequences obtained from equenceModification function.
  for(auto const& item:sequences)
  {
      
      sequence = item.first;
      std::replace( sequence.begin(), sequence.end(), 'U', 'T' ); //Replace all Us with Ts.
      for (int j = 0; j < times; ++j)
      {
	output_file << read_id_token << "_" << modification_type <<  "_" << item.second << "_" << counter << "\n";
	output_file << sequence << "\n";
	++counter;
	annotation_counter[read_id_token]++;
      }
    }
}

void bruteForce(std::string input, std::ofstream& output, std::string file_type, int option, int times){
 
  //Reading files
  read_id_counter = 1;

  int counter = 1;
  
  int flag = 0;
  
  std::string read_id;
  //Stream class to read from input file.
  std::ifstream input_file;
  input_file.open(input.c_str());
  //Stream class to write on output file.
  
//   std::ofstream output_file;
//   output_file.open(output.c_str());
  
  std::string transcript_counts_file = output_path + "simulated_miRNA_BF_counts.tsv";
  std::ofstream transcript_file;
  transcript_file.open(transcript_counts_file.c_str());
  
  for( std::string line; getline( input_file, line ); )
  {
      //if (std::regex_match (line, std::regex("(>hsa)(.*)"))) { //FIXME regex should be replaced with >
	std::size_t found = line.find(">"+species);
	if(found!=std::string::npos){
      //if (std::regex_match (line, std::regex("(>species)(.*)"))) { //FIXME regex should be replaced with >
	
	flag = 1;
	read_id = line;
      }
      
      if(counter % 2 == 0 && flag == 1){ //Sequence line

	std::string read_id_token = read_id.substr(0, read_id.find(" ")); //Obtain the read id from the input file.
	
	std::vector<std::string> read_id_splt = split(read_id, ' ');
	
	//std::vector<std::string> sequences;
	std::vector<std::pair<std::string, std::string>> sequences;
	
        std::string name = read_id_splt[0];
 	std::string mimat = read_id_splt[1];
	//sequences = sequenceModification(option, line); //Call sequenceModification function.
	sequences = sequenceModification(option, line, name, mimat); //Call sequenceModification function.
	
	writeSequences(sequences, output, read_id_token, times, annotation_counter); //Write sequences to fasta file
	//and store transcript couts.
	
	flag = 0;
      }
      ++counter;
  }
  
  //Write transcript counts to output file.
  transcript_file << "Transcript" << "\t" << "Count" << "\n"; //Write header.
  for (annotationCounter::iterator p = annotation_counter.begin(); p != annotation_counter.end(); ++p) {
          std::string id = p->first;
          id.replace(0,1,"");
          transcript_file << id << "\t" << p->second << "\n";
   }
   
   //output_file.close();
}

