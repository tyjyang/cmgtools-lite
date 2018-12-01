#include "TSystem.h"
#include "TROOT.h"

#include <iostream>
#include <string>
#include <cstdlib> //as stdlib.h                    
#include <cstdio>

#include "../interface/utility.h"

using namespace std;

// root -l -b -q makeLoopTree.C++

void makeLoopTree(const bool isMuon = false, 
		  const string& outfileName = "wmass_varhists.root") 
{

  // load source codes with ++, so that they are always compiled (you never know ...)

  string cmssw_base = getEnvVariable("CMSSW_BASE");
  cout << "CMSSW_BASE = " << cmssw_base << endl;

  string host_name = getEnvVariable("HOSTNAME");
  cout << "HOSTNAME = " << host_name << endl;

 
  cout << "Loading functions.cc" << endl;
  gROOT->ProcessLine(Form(".L %s/src/CMGTools/WMass/python/plotter/functions.cc+",cmssw_base.c_str())); 
  cout << "Loading functionsWMass.cc" << endl;
  gROOT->ProcessLine(Form(".L %s/src/CMGTools/WMass/python/plotter/w-helicity-13TeV/functionsWMass.cc+",cmssw_base.c_str()));
  cout << "Loading loopNtuplesSkeleton.cc" << endl;
  gROOT->ProcessLine(".L loopNtuplesSkeleton.C++");

  //string command = "loopNtuplesSkeleton(\"/eos/cms/store/cmst3/group/wmass/mciprian/TREES_1LEP_80X_V3_WSKIM_NEW/\",\"./\",\"" + outfileName + "\")";
  string command = "";
  if (isMuon) {

    if (host_name.find("lxplus") != string::npos) 
      command = "loopNtuplesSkeleton(\"/eos/cms/store/cmst3/group/wmass/mciprian/TREES_1LEP_80X_V3_SIGSKIM_WMUNU_FULLSEL_NOMT_V2/\",\"./\",\"" + outfileName + "\",true)";
    else if (host_name.find("pccmsrm") != string::npos) 
      command = "loopNtuplesSkeleton(\"/u2/mciprian/TREES_13TeV/muon/signalSkim/\",\"./\",\"" + outfileName + "\",true)";

  } else {

    if (host_name.find("lxplus") != string::npos) 
      command = "loopNtuplesSkeleton(\"/eos/cms/store/cmst3/group/wmass/mciprian/TREES_1LEP_80X_V3_SIGSKIM_WENU_FULLSEL_NOMT/\",\"./\",\"" + outfileName + "\",false)";
    else if (host_name.find("pccmsrm") != string::npos) 
      command = "loopNtuplesSkeleton(\"/u2/mciprian/TREES_13TeV/electron/signalSkim/\",\"./\",\"" + outfileName + "\",false)";
  }

  //string command = "loopNtuplesSkeleton(\"/eos/cms/store/cmst3/group/wmass/w-helicity-13TeV/trees/TREES_electrons_1l_V6_TINY/\",\"./\",\"" + outfileName + "\")";
  cout << "Executing " << command << endl;
  gROOT->ProcessLine(command.c_str());
  cout << endl;                                          
  cout << "===========================" << endl;                    
  cout << " THE END!" << endl;                          
  cout << "===========================" << endl;         



}
