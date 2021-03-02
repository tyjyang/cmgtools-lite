#ifndef TnPNtuplesBase_cxx
#define TnPNtuplesBase_cxx

#include "TnPNtuplesBase.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TLorentzVector.h>
#include <TH1F.h>

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
#include <utility>      // std::pair
#include <iomanip> //for input/output manipulators


// #include <TROOT.h>
// #include <TChain.h>
// #include <TFile.h>
// #include <TString.h>
// #include <TH1F.h>

using namespace std;

void TnPNtuplesBase::setFlavor(int flavor){
  fFlavor = flavor;
}

void TnPNtuplesBase::setOutfile(TString outfilepath){
  fOutfile = outfilepath;
}

bool TnPNtuplesBase::isTagLepton(int jj){
  return true;
}

bool TnPNtuplesBase::getStatusFlag(Int_t flags, int index) { 
  return ((flags >> index) & 1);
}

//std::map< UInt_t, std::vector< std::pair<UInt_t,UInt_t> > > makeMapFromJson(const string myJsonFile = "") {
void TnPNtuplesBase::makeMapFromJson(const string myJsonFile = "") {

  // format is 
  // run: [ls1,ls2] [ls3,ls4] [...]  note that spaces are important
  // given a json, this format can be obtained with myFormatJson.py
  string run;
  string lumiBlocks;
  std::map< UInt_t, std::vector< std::pair<UInt_t,UInt_t> > > runsAndLumiBlocks;

  string jsonFile = Form(myJsonFile.c_str());

  ifstream inputFile(jsonFile.c_str());

  //cout << "Printing content of " << myJsonFile << endl;

  if (inputFile.is_open()) {

    while (inputFile >> run) {

      // read line without first object that was put in run (the space separates objects in the line)
      getline(inputFile, lumiBlocks);  

      run.assign(run,0,6); // run has 6 digits
      //cout << run << " --> ";

      vector< pair<UInt_t,UInt_t> > blocks;
      stringstream ss(lumiBlocks);     
      string block;

      while (ss >> block) {

    	//cout << block << " ";

    	size_t pos = block.find(",");
    	UInt_t LSin, LSfin;
    	string num1, num2;
    	// we have block = "[a,b]" where a and b are integers. We want to get a and b
    	num1.assign(block,1,pos);
    	num2.assign(block,pos+1,block.size()-1);
    	LSin  = (UInt_t) std::stoi(num1);
    	LSfin = (UInt_t) std::stoi(num2);
    	// cout << "LSin,LSfin = " << LSin << "," << LSfin << endl;
    	blocks.push_back(std::make_pair(LSin,LSfin));

      }

      //runsAndLumiBlocks[run] = blocks;
      //cout << endl;
      runsAndLumiBlocks.insert ( std::pair< UInt_t, std::vector< std::pair<UInt_t,UInt_t > > >((UInt_t) stoi(run), blocks) );
      //      runsAndLumiBlocks.at(stoi(run)) = blocks;

    }

    /////////////////////////
    // check that it works
    cout << "printing map ..." << endl;
    for (std::map<UInt_t, vector< pair<UInt_t,UInt_t> > >::iterator it = runsAndLumiBlocks.begin(); it != runsAndLumiBlocks.end(); ++it) {
      cout << it->first << " --> "; 
      for (UInt_t i = 0; i < it->second.size(); i++) {
    	cout << "[" << it->second.at(i).first << "," << it->second.at(i).second << "]  ";
      } 
      cout << endl;
    }
    cout << endl;
    ///////////////////////////

  } else {
    
    cout << "Error in makeMapFromJson(): could not open file " << jsonFile << endl;
    exit(EXIT_FAILURE);

  }


  fJsonMap = runsAndLumiBlocks;
  //return runsAndLumiBlocks;

}

Bool_t TnPNtuplesBase::isGoodRunLS(Bool_t isData, UInt_t run, UInt_t lumis) {
  
  // for MC thsi fucntion always return true
  if (not isData) return true;

  // I should make it load the json somewhere else, something like loading the FR file
  // if (theJsonMap.empty()) theJsonMap = makeMapFromJson(formattedJson);

  // if (theJsonMap.empty()) {
  //   cout << "Warning in isGoodRunLS(): mymap is empty. Returning false, but please check what's happening!" << endl;
  //   return false;
  // }

  if ( fJsonMap.find(run) == fJsonMap.end() ) return false; // run not found

  Bool_t LSfound = false;

  for (UInt_t i = 0; i < fJsonMap.at(run).size() && !LSfound; ++i) {
    
    // evaluate second value, skip if lumis is bigger (block does not contain it)
    if (lumis >  fJsonMap.at(run).at(i).second) continue;  
    // if arrive here, check lower boundary
    if (lumis >= fJsonMap.at(run).at(i).first ) LSfound = true;

  }

  return LSfound;
    
}

void TnPNtuplesBase::Loop(int maxentries) {} // Loop method

// Histos booking
void TnPNtuplesBase::bookOutputTree() 
{
  outFile_ = new TFile(fOutfile, "RECREATE");    
  outFile_->cd();

  //cddir = outFile_->mkdir("IDIsoToHLT");
  //cddir->cd();    

  std::cout << "Booking output tree" << endl;
  outTree_ = new TTree("fitter_tree", "fitter_tree");

  outTree_->Branch("tag_pt"             , &tag_pt             , "tag_pt/F");
  outTree_->Branch("tag_eta"            , &tag_eta            , "tag_eta/F");
  outTree_->Branch("tag_matchMC"        , &tag_matchMC        , "tag_matchMC/I");
  outTree_->Branch("probe_truept"       , &probe_truept       , "probe_truept/F");
  outTree_->Branch("probe_trueeta"      , &probe_trueeta      , "probe_trueeta/F");
  outTree_->Branch("probe_pt"           , &probe_pt           , "probe_pt/F");
  outTree_->Branch("probe_eta"          , &probe_eta          , "probe_eta/F");
  outTree_->Branch("probe_dxy"          , &probe_dxy          , "probe_dxy/F");
  outTree_->Branch("probe_dz"           , &probe_dz           , "probe_dz/F");
  outTree_->Branch("probe_iso"          , &probe_iso          , "probe_iso/F");
  outTree_->Branch("probe_iso03"        , &probe_iso03        , "probe_iso03/F");
  outTree_->Branch("probe_chgiso03"     , &probe_chgiso03     , "probe_chgiso03/F");
  outTree_->Branch("probe_mediumId"     , &probe_mediumId     , "probe_mediumId/I");
  outTree_->Branch("probe_sc_eta"       , &probe_sc_eta       , "probe_sc_eta/F");
  outTree_->Branch("probe_phi"          , &probe_phi          , "probe_phi/F");
  outTree_->Branch("probe_charge"       , &probe_charge       , "probe_charge/F");
  outTree_->Branch("probe_pdgId"        , &probe_pdgId        , "probe_pdgId/I");

  outTree_->Branch("probe_triggerMatch" , &probe_triggerMatch , "probe_triggerMatch/I");
  outTree_->Branch("probe_isGlobal"     , &probe_isGlobal     , "probe_isGlobal/I");
  outTree_->Branch("probe_isTracker"    , &probe_isTracker    , "probe_isTracker/I");
  outTree_->Branch("probe_isMuon"       , &probe_isMuon       , "probe_isMuon/I");

  outTree_->Branch("probe_matchMC"      , &probe_matchMC      , "probe_matchMC/I");
  outTree_->Branch("probe_tightCharge"  , &probe_tightCharge  , "probe_tightCharge/I");
  outTree_->Branch("probe_fullLepId"    , &probe_fullLepId    , "probe_fullLepId/I");
  outTree_->Branch("probe_alsoTag"      , &probe_alsoTag      , "probe_alsoTag/I");
  outTree_->Branch("pair_mass"          , &pair_mass          , "pair_mass/F");

  // some event variables
  outTree_->Branch("nvtx"               , &nvtx               , "nvtx/I");
  outTree_->Branch("run"                , &_run               , "run/I");
  outTree_->Branch("thisEntry"          , &thisEntry          , "thisEntry/I");
  outTree_->Branch("puw_inc"            , &puw_inc            , "puw_inc/F");
  outTree_->Branch("puw_bf"             , &puw_bf             , "puw_bf/F");
  outTree_->Branch("puw_gh"             , &puw_gh             , "puw_gh/F");
  outTree_->Branch("totWeight"          , &totWeight          , "totWeight/F");
  outTree_->Branch("mcTrue"             , &mcTrue             , "mcTrue/I");

  std::cout << "Booking output histos for event breakdown" << endl;
  h_entries   = new TH1F("h_entries"  , "h_entries"  , 10,   0., 10.);
  h_selection = new TH1F("h_selection", "h_selection",  6, -0.5, 5.5);
  h_entries->Sumw2();
  h_selection->Sumw2();
}

// To compute the pileup weight
float TnPNtuplesBase::puw_2016(int nTrueInt, int period=0)
{ 
// inclusive
  float _puw2016_nTrueInt_36fb[100] = {0.6580788810596797, 0.4048055020122533, 0.8615727417125882, 0.78011444438815, 0.7646792294649241, 0.4428435234418468, 0.21362966755761287, 0.1860690038399781, 0.26922404664693694, 0.3391088547725875, 0.4659907891982563, 0.6239096711946419, 0.7379729344299205, 0.8023612182327676, 0.8463429262728119, 0.8977203788568622, 0.9426445582538644, 0.9730606480643851, 0.9876071912898365, 0.9967809848428801, 1.0101217492530405, 1.0297182331139065, 1.0495800245321394, 1.062993596017779, 1.0719371123276982, 1.0807330103942534, 1.0874070624201613, 1.092162599537252, 1.1002282592549013, 1.1087757206073987, 1.1145692898936854, 1.1177137330250873, 1.1227674788957978, 1.127920444541154, 1.1304392102417709, 1.1363304087643262, 1.1484550576573724, 1.1638015070216376, 1.1833470331674814, 1.2079950937281854, 1.2274834937176624, 1.2444596857655403, 1.269207386640866, 1.2896345429784277, 1.3441092120892817, 1.318565807546743, 1.3245909554149022, 1.2438485820227307, 1.1693858006040387, 0.9959473986900877, 0.8219068296672665, 0.6541529992036196, 0.4693498259046716, 0.32997121176496014, 0.31310467423384075, 0.23393084684715088, 0.2693460933531134, 0.2514663792161372, 0.3548113113715822, 0.5892264892080126, 0.6554126787320194, 0.7565044847439312, 1.0033865725393603, 1.4135667773217673, 0.8064239843371759, 1.096753894690768, 1.0, 1.0, 1.0, 1.0, 1.0, 1.409161148159025, 2.436803424843253, 1.0, 1.0, 1.5123015283558174, 1.0, 1.0, 1.0, 2.184264934077207, 0.44722379540985835, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

// preVFP
  float _puw2016_nTrueInt_BtoF[100] = {0.017867211276858426, 0.25019181678340197, 0.7305862075344697, 0.8203306355218927, 0.8568507658433443, 0.5021252568573517, 0.26598094658249777, 0.2525594126097894, 0.43091072564779115, 0.5816151668719733, 0.7400719694832786, 0.8635841784600747, 0.9396947027327512, 0.9867166485834531, 1.020291825978338, 1.058677896425578, 1.0951610451340548, 1.1270599475402712, 1.149199627966607, 1.161062844833281, 1.1570323009561523, 1.1408277403854514, 1.1237795117739837, 1.1080566938230716, 1.0897761375745778, 1.064304242937247, 1.0256371301024838, 0.974592047960344, 0.919208105883519, 0.8614029922190752, 0.80250010071568, 0.7447070188841107, 0.691304388080654, 0.6402835441510876, 0.5891723690900388, 0.5400972141379156, 0.49313582436893444, 0.4463956167604181, 0.400690015469447, 0.3570993263445747, 0.313769141057574, 0.272973076482914, 0.23756151770204517, 0.2052480930611912, 0.18165295449088267, 0.15154171165962477, 0.13020442654298184, 0.1060542785336088, 0.08931962151840153, 0.07320668426537356, 0.06699868236915886, 0.07342602448799117, 0.09022522247731071, 0.1221948811200131, 0.21751127388127575, 0.26185148964193156, 0.3993515054053398, 0.4259481886684673, 0.6348999885572649, 1.0768669803515931, 1.2073680188570837, 1.3977538479690192, 1.8560151439282793, 2.6159138791249634, 1.492621348422977, 2.030155370346692, 1.0, 1.0, 1.0, 1.0, 1.0, 2.608617469665849, 4.5109759375776735, 1.0, 1.0, 2.7995528359963338, 1.0, 1.0, 1.0, 4.043483083032464, 0.8278949263309181, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

// postVFP
  float _puw2016_nTrueInt_GtoH[100] = {1.5094769667755739, 0.5678507573322055, 1.0105884486266141, 0.7259960323564019, 0.6619453439299897, 0.37363396015230216, 0.15239333007202369, 0.10791755923476429, 0.07918472234473897, 0.05419619726622047, 0.14361053792634917, 0.34283985083406093, 0.5010734484221103, 0.5859542596067007, 0.6417013874657507, 0.7081679228576429, 0.7630928885638684, 0.7927011423655624, 0.7977355318365416, 0.804189334635704, 0.8379574350304874, 0.8991008934966008, 0.9625166892724227, 1.0100064501247878, 1.0501599728932576, 1.0991686204042044, 1.159807646806213, 1.230197853416741, 1.3121394462817406, 1.3989162630341463, 1.4799356011552773, 1.5557799985845535, 1.6310266488456922, 1.7003552120673373, 1.7676611780885065, 1.8405616412980672, 1.9250530644650858, 2.010419966610567, 2.10422235818776, 2.2195187847594804, 2.2997062399954857, 2.395010762430828, 2.4803953120930684, 2.5918893616537275, 2.7047108284644694, 2.71073345533101, 2.6519181014415363, 2.6277548910588657, 2.3664596415368297, 2.1655717980383105, 1.8274047657200394, 1.3313115741327561, 0.9113129353116672, 0.5869561379336159, 0.4284308024180917, 0.20595003106437168, 0.10999572232736675, 0.049536718778780756, 0.035006854789963834, 0.012922619333388894, 0.009214760376410974, 0.009219346474718641, 0.0015371306555079958, 0.001800361515424057, 0.0003153024223498953, 0.00010484990280832637, 1.0, 9.561632218479775e-05, 4.4263424947054785e-05, 1.0, 1.0, 8.750970959266416e-06, 1.9850705244600614e-06, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0235991128869821e-08, 1.3636783092893852e-09, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};


  if (nTrueInt<100) {
    if      (period==0) return _puw2016_nTrueInt_36fb[nTrueInt]; 
    else if (period==1) return _puw2016_nTrueInt_BtoF[nTrueInt]; 
    else if (period==2) return _puw2016_nTrueInt_GtoH[nTrueInt]; 
    else { std::cout<< " you are giving a wrong period to the pileup weight function, returning 0 for the events!!!!" << std::endl; return 0.;}
  }
  else  return 1.;
}
#endif
