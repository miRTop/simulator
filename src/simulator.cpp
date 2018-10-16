#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<string>
#include<fstream>
#include <regex>
#include <map>
#include <getopt.h>
#include "project/simulator.h"
#include "project/brute_force.h"
#include "project/utils.h"

std::string modification_type;
 
std::string mode;

std::string output_path;

std::string fixed_hairpin = "hairpin.fixed.fa";

std::string fixed_hairpin_species_spec = "hairpin.fixed.species.spec.fa";

int times; /* Number of generated sequences of the same type */

std::string species = "hsa"; /* Default species for which the simulation will be performed */

std::vector <std::pair<string, string>> mi_mimat;

void helpMessage (std::string program_name) { //FIXME
         printf("\nUsage: %s [OPTIONS]\n", program_name.c_str());
	 printf("  -a <file>	: Path to mature.fa file.                   \n");
	 printf("  -b <file>	: Path to hairpin.fa file.                   \n");
	 printf("  -g <file>	: Path to miRBase gtf file.                   \n");
	 printf("  -o <path>	: Output path.                   \n");
	 printf("  -s <string>	: Species for the generated sequences (e.g. hsa).                   \n");
         printf("  -m <int>	: Mode for the  brute force execution.                \n");
	 printf("  -t <int>	: Number of identical modified sequences sequences. \n");
	 printf("  Available modes:\n");
	 printf("   	0: Canonical miRNAs \n");   
	 printf("	1: 5' templated addition \n");
         printf("	2: 3' templated addition \n");
         printf("	3: 5' non-templated addition (A and Ts) \n");
	 printf("	4: 3' non-templated addition (A and Ts)\n");
         printf("	5: SNPs \n");
	 printf("	6: Indels \n");
	 printf("	7: 3' templated losses \n");
         printf("	8: 5' templated losses \n");
         printf("	9: Options 0-8 \n");
         printf("\n");
         exit(1);
	 
}

/* Concatenate sequences lines in hairpin file. */
void fixHairpinFile (std::string input_file, std::string output_file) {
  
  std::ifstream input;
  std::ofstream output;
  
  input.open(input_file.c_str());
  output.open(output_file.c_str());
  int start = 0;
  std::string hairpin;

	for( std::string line; getline( input, line ); )
          {
	    if( std::regex_match (line, std::regex("(>)(.*)")) ){
	      
	      if(start) output << hairpin << "\n";
	      
	      output << line << "\n";

	      hairpin = "";
	      
	    }else{
	      
	      hairpin = hairpin + line;
	      start = 1;
	      
	    }
	    
	  }
	  
	  output << hairpin << "\n"; //Write the last line
	  output.close();
}

//Create mapping of MI and MIMAT IDs form GTF file.
void mapMItoMIMAT(std::string input, std::vector <std::pair<std::string, std::string>> &mi_mimat){
  
  std::ifstream fp;
  fp.open(input.c_str());
  std::smatch match;
   	for( std::string line; getline( fp, line ); )
           {
	      if (std::regex_match (line, std::regex("(.*)(Derives_from)(.*)"))) 
	      {
		std::vector<std::string> tab = split(line, '\t');
		std::vector<std::string> semicolon = split(tab[8], ';');
		std::string mimat = semicolon[0].substr(3,semicolon[0].size());
		std::string mi = semicolon[3].substr(13,semicolon[3].size());
		mi_mimat.emplace_back(mi, mimat);
	      }
	  
	   }
  
  
}

void line_by_species(std::string input_file, std::string output_file){
  
	int flag = 0;
	int counter = 1;
        std::ifstream input;
	std::ofstream output;
        input.open(input_file.c_str());
 	output.open(output_file.c_str());
	
 	for( std::string line; getline( input, line ); )
           {
	  if (line.find('>'+species) != std::string::npos) {
	    output << line << "\n";
	    flag = 1;
	  }
	    if(counter % 2 == 0 && flag == 1){ //Sequence line
	      output << line << "\n";
	      flag = 0;
	    }
	  counter++;
	   }
        output.close();
  
}

int main(int argc, char **argv)
{
  
  int aflag = 0;
  int bflag = 0;
  char *cvalue = NULL;
  std::string fasta;
  int index;

   //int             opt;
   const char    * short_opt = "hm:t:a:b:g:o:s:";
   struct option   long_opt[] =
   {
      {"help",          no_argument,       nullptr, 'h'},
      {"mature",         required_argument, nullptr, 'a'},
      {"hairpin",         required_argument, nullptr, 'b'},
      {"gtf",         required_argument, nullptr, 'g'},
      {"output",         required_argument, nullptr, 'o'},
      {"mode",        	required_argument, nullptr, 'm'},
      {"times",        	required_argument, nullptr, 't'},
      {nullptr,            no_argument, nullptr, 0  }
   };

   if(argc == 1) helpMessage(argv[0]); //If no arguments have been passed 
   
   std::string input;
   std::string mature;
   std::string hairpin;
   std::string gtf;

   int mode;   
   
   
//    while((opt = getopt_long(argc, argv, short_opt, long_opt, NULL)) != -1)
//    {
   int flag = 0;
   while(true)
   {
      const auto opt = getopt_long(argc, argv, short_opt, long_opt, NULL);
      
      if (-1 == opt)
      break;
      
      switch(opt) //FIXME Add modes for BF, MonteCarlo +++
      {
//          case -1:       /* no more arguments */
//          case 0:        /* long options toggles */
//          break;

         case 'a':
	 mature = optarg;
	 flag++;
         break;
	 
	 case 'b':
	 hairpin = optarg;
	 flag++;
         break;
	 
	 case 'g':
	 gtf = optarg;
	 flag++;
         break;
	 
	 case 'o':
	 output_path = optarg;
	 flag++;
         break;
	 
	 case 'm':
	 mode = atoi(optarg);
	 flag++;
         break;
	 
	 case 't':
	 times = atoi(optarg);
	 flag++;
         break;
	 
         case 'h':
         helpMessage(argv[0]);
	 break;

         //case ':':
         case '?':
	 //  if(optopt
         fprintf(stderr, "Try `%s --help' for more information.\n", argv[0]);
         return(-2);

         //default:
//          fprintf(stderr, "%s: invalid option -- %c\n", argv[0], opt);
//          fprintf(stderr, "Try `%s --help' for more information.\n", argv[0]);
//          return(-2);
       }
     // };
   //};
   }

  //FIXME write better way to communicate with user.
  if(flag != 6) helpMessage(argv[0]); //If not enough arguments have been passed. 
  
  path_check(output_path);
  
  fixHairpinFile(hairpin, fixed_hairpin);
  
  mapMItoMIMAT(gtf, mi_mimat);
  
  std::string mature_species_spec = "mature.species.spec.fa";
  line_by_species(mature, mature_species_spec); //Select lines with specific species.
  
  line_by_species(hairpin, fixed_hairpin_species_spec); //Select lines with specific species.
  
  std::string output = output_path + "simulated_miRNA_BF.fa";
  
  std::ofstream output_fp; //Output FASTA file.
  output_fp.open(output.c_str());
  
  input = mature_species_spec;
  
  if(mode == 1){//FIXME Better calls?
    modification_type = "5prime_templated_addition";
  }else if (mode == 2){
    modification_type = "3prime_templated_addition";
  }else if (mode == 3){
    modification_type = "5prime_addition";
  }else if (mode == 4){
    modification_type = "3prime_addition";
  }else if (mode == 5){
    modification_type = "SNP";
  }else if (mode == 6){
    modification_type = "indel";
  }else if (mode == 7){
    modification_type = "3prime_loss";
  }else if (mode == 8){
    modification_type = "5prime_loss";
  }
  
  if(mode == 9){
    modification_type = "no_modification";
    bruteForce(input, output_fp, "fasta", 0, times);
    modification_type = "5prime_templated_addition";
    bruteForce(input, output_fp, "fasta", 1, times);
    modification_type = "3prime_templated_addition";
    bruteForce(input, output_fp, "fasta", 2, times);
    modification_type = "5prime_addition";
    bruteForce(input, output_fp, "fasta", 3, times);
    modification_type = "3prime_addition";
    bruteForce(input, output_fp, "fasta", 4, times);
    modification_type = "SNP";
    bruteForce(input, output_fp, "fasta", 5, times);
    modification_type = "indel";
    bruteForce(input, output_fp, "fasta", 6, times);
    modification_type = "5prime_loss";
    bruteForce(input, output_fp, "fasta", 7, times);
    modification_type = "3prime_loss";
    bruteForce(input, output_fp, "fasta", 8, times);
  }else{
    bruteForce(input, output_fp, "fasta", mode, times);
  }
 output_fp.close();

 delete_file(mature_species_spec);
 delete_file(fixed_hairpin);
 delete_file(fixed_hairpin_species_spec);
 
}
