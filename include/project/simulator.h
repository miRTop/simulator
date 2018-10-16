#ifndef SIMULATOR
#define SIMULATOR

#include<stdio.h>
#include<iostream>
#include <vector>
#include <map>

extern std::vector <std::pair<std::string, std::string>> mi_mimat;

extern std::string modification_type;

extern std::string output_path;

extern std::string fixed_hairpin_species_spec;

extern std::string species;

void helpMessage(std::string program_name);

void fixHairpinFile (std::string input_file, std::string output_file);

void mapMItoMIMAT(std::string input, std::vector <std::pair<std::string, std::string>> &mi_mimat);

void line_by_species(std::string input_file, std::string output_file);

#endif
