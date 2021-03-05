#ifndef JSON_MANAGER_H
#define JSON_MANAGER_H

#include <stdio.h>
#include <stdlib.h>
#include <cstdlib> //as stdlib.h
#include <cstdio>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>      // std::istringstream ; to read array of numbers from a line in a file
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <utility>      // std::pair
#include <iomanip> //for input/output manipulators
#include <boost/functional/hash.hpp>

#include "TROOT.h"

using namespace std;

//std::unordered_map< UInt_t, std::vector< std::pair<UInt_t,UInt_t> > > theJsonMap;
// std::unordered_map< UInt_t, std::vector< std::pair<UInt_t,UInt_t> > > jsonMap_preVFP;
// std::unordered_map< UInt_t, std::vector< std::pair<UInt_t,UInt_t> > > jsonMap_postVFP;
// std::unordered_map< UInt_t, std::vector< std::pair<UInt_t,UInt_t> > > jsonMap_all;
std::unordered_map< UInt_t, std::vector< std::pair<UInt_t,UInt_t> > > theJsonMap;

//==========================================================

std::unordered_map< UInt_t, std::vector< std::pair<UInt_t,UInt_t> > > makeMapFromJson(const std::string& myJsonFile = "") {

  // format is 
  // run: [ls1,ls2] [ls3,ls4] [...]  note that spaces are important
  // given a json, this format can be obtained with python/plotter/myFormatJson.py
  std::string run;
  std::string lumiBlocks;
  std::unordered_map< UInt_t, std::vector< std::pair<UInt_t,UInt_t> > > runsAndLumiBlocks;

  ifstream inputFile(myJsonFile.c_str());

  //cout << "Printing content of " << myJsonFile << endl;

  if (inputFile.is_open()) {

    while (inputFile >> run) {

      // read line without first object that was put in run (the space separates objects in the line)
      getline(inputFile, lumiBlocks);  
      run.assign(run,0,6); // run has 6 digits
      //cout << run << " --> ";

      std::vector< pair<UInt_t,UInt_t> > blocks;
      stringstream ss(lumiBlocks);     
      std::string block;

      while (ss >> block) {

    	//cout << block << " ";
    	size_t pos = block.find(",");
    	UInt_t LSin, LSfin;
    	string num1, num2;
    	// we have block = "[a,b]" where a and b are integers. We want to get a and b
    	num1.assign(block,1,pos);
    	num2.assign(block,pos+1,block.size()-1);
    	LSin  = static_cast<UInt_t>(std::stoi(num1));
    	LSfin = static_cast<UInt_t>(std::stoi(num2));
    	// cout << "LSin,LSfin = " << LSin << "," << LSfin << endl;
    	blocks.push_back(std::make_pair(LSin,LSfin));

      }

      runsAndLumiBlocks.insert ( std::pair< UInt_t, std::vector< std::pair<UInt_t,UInt_t > > >(static_cast<UInt_t>(stoi(run)), blocks) );
      
    }

    /////////////////////////
    // check that it works
    // cout << "printing map for " << myJsonFile << " ..." << endl;
    // for (std::unordered_map<UInt_t, vector< pair<UInt_t,UInt_t> > >::iterator it = runsAndLumiBlocks.begin(); it != runsAndLumiBlocks.end(); ++it) {
    //   cout << it->first << " --> "; 
    //   for (UInt_t i = 0; i < it->second.size(); i++) {
    // 	cout << "[" << it->second.at(i).first << "," << it->second.at(i).second << "]  ";
    //   } 
    //   cout << endl;
    // }
    // cout << endl;
    ///////////////////////////

  } else {
    
    std::cout << "Error in makeMapFromJson(): could not open file " << myJsonFile << std::endl;
    exit(EXIT_FAILURE);

  }

  return runsAndLumiBlocks;

}

//==========================================================

void initializeJson() {

  // format is 
  // run: [ls1,ls2] [ls3,ls4] [...]  note that spaces are important
  // given a json, this format can be obtained with python/plotter/myFormatJson.py
  // jsonMap_preVFP  = makeMapFromJson("./pileupStuff/json_preVFP_myFormatJson.txt");
  // jsonMap_postVFP = makeMapFromJson("./pileupStuff/json_postVFP_myFormatJson.txt");
  // checked that there are no duplicate keys, i.e. no common runs between the two eras
  // jsonMap_all = jsonMap_preVFP;
  // jsonMap_all.insert(jsonMap_postVFP.begin(), jsonMap_postVFP.end());
  // jsonMap_all = makeMapFromJson("./pileupStuff/json_all_myFormatJson.txt");

  theJsonMap = makeMapFromJson("./pileupStuff/json_all_myFormatJson.txt");
  
}

//==========================================================

Bool_t isGoodRunLS(const Bool_t& isData, const UInt_t& run, const UInt_t& lumis) {
  
  // for MC this function always returns true
  if (not isData) return true;

  // std::unordered_map< UInt_t, std::vector< std::pair<UInt_t,UInt_t> > > theJsonMap; // would it be better as a pointer?
  // if       (era == 1) theJsonMap = jsonMap_preVFP;
  // else if  (era == 2) theJsonMap = jsonMap_postVFP;
  // else                theJsonMap = jsonMap_all;

  //std::cout << "run,lumi = " << run << "," << lumis << " --> ";
  if ( theJsonMap.find(run) == theJsonMap.end() ) {
    //std::cout << "false" << std::endl;
    return false; // run not found
  }
  
  Bool_t LSfound = false;

  for (UInt_t i = 0; i < theJsonMap.at(run).size() && !LSfound; ++i) {
    
    // evaluate second value, skip if lumis is bigger (block does not contain it)
    if (lumis > theJsonMap.at(run).at(i).second) continue;  
    // if arrive here, check lower boundary
    if (lumis >= theJsonMap.at(run).at(i).first ) LSfound = true;

  }
  //std::cout << ((LSfound) ? "true" : "false")  << std::endl;
  
  // if (not LSfound) {
  //   std::cout << "NOT IN JSON: run,lumi = " << run << "," << lumis << std::endl;
  // } else {
  //   std::cout << "RUN/LUMI FOUND IN JSON: run,lumi = " << run << "," << lumis << std::endl;
  // }
  return LSfound;
    
}

//==========================================================

#endif
