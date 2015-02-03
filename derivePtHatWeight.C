#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TMath.h"

#include "TH1F.h"
#include "TCanvas.h"


Int_t ptHatCuts_PYTH[7] = {15, 30, 50, 80, 120, 170, 1000000};
Float_t crossSections_PYTH[7] = {.2034, .01075, .001025, .00009865, .00001129, .000001465, 0.000000000};

Int_t ptHatCuts_PYTHPawan[7] = {15, 30, 50, 80, 120, 170, 10000000};
Float_t crossSections_PYTHPawan[7] = {.20340, .01075, .001025, .00009865, .00001129, .000001465, 0.000000000};

Int_t ptHatCuts_PYTHHYD[11] = {15, 30, 50, 80, 100, 120, 170, 220, 280, 370, 1000000};
Float_t crossSections_PYTHHYD[11] = {.20340, .01075, .001025, .00009865, .00003069, .00001129, .000001465, .0000002837, .00000005323, .000000005934, .0000000000};

Int_t ptHatCuts_PYTH_HITrk[11] = {15, 30, 50, 80, 120, 170, 220, 280, 370, 10000000};
Float_t crossSections_PYTH_HITrk[10] = {.20340, .01075, .001025, .00009865, .00001129, .000001465, .0000002837, .00000005323, .000000005934};

void derivePtHatWeights(const Int_t numCut, Int_t ptHatCuts[], Float_t crossSect[], std::string fList = "")
{
  std::string buffer;
  std::vector<std::string> listOfFiles;
  int nLines = 0;
  ifstream inFile(fList.data());

  std::cout << fList << std::endl;
  std::cout << inFile.is_open() << std::endl;

  if(!inFile.is_open()){
    std::cout << "Error opening file. Exiting." <<std::endl;
    return;
  }
  else{
    while(true){
      inFile >> buffer;
      if(inFile.eof()) break;
      listOfFiles.push_back(buffer);
      nLines++;
    }
  }

  std::cout << "FileList Loaded" << std::endl;

  TChain* ptHatChain_p = new TChain("akVs3CaloJetAnalyzer/t");

  for(Int_t iter = 0; iter < (Int_t)(listOfFiles.size()); iter++){
    ptHatChain_p->Add(listOfFiles[iter].c_str());

    std::cout << listOfFiles[iter] << std::endl;
  }

  Float_t ptHat_ = 0;

  ptHatChain_p->SetBranchStatus("*", 0);
  ptHatChain_p->SetBranchStatus("pthat", 1);
  ptHatChain_p->SetBranchAddress("pthat", &ptHat_);

  Int_t nEntries = ptHatChain_p->GetEntries();
  std::cout << nEntries << std::endl;

  Int_t hatEntries[numCut];
  Float_t hatWeight[numCut];

  for(Int_t iter = 0; iter < numCut; iter++){
    hatEntries[iter] = 0;
    hatWeight[iter] = 0;
  }

  for(Int_t evtIter = 0; evtIter < nEntries; evtIter++){
    ptHatChain_p->GetEntry(evtIter);

    for(Int_t hatIter = 0; hatIter < numCut; hatIter++){
      if(ptHat_ > ptHatCuts[hatIter] && ptHat_ < ptHatCuts[hatIter+1]){
	hatEntries[hatIter]++;
	break;
      }
    }
  }
  
  for(Int_t hatIter = 0; hatIter < numCut; hatIter++){
    std::cout << hatIter << std::endl;
    std::cout << "  hatEntries: " << hatEntries[hatIter] << std::endl;
    hatWeight[hatIter] = (crossSect[hatIter] - crossSect[hatIter+1])/hatEntries[hatIter];
    std::cout << "  hatWeight: " << hatWeight[hatIter] << std::endl;    
  }

  TH1F* testHist_NoWeight_p = new TH1F("testHist_NoWeight_h", "testHist_NoWeight_h", 100, 0, 500);
  TH1F* testHist_Weight_p = new TH1F("testHist_Weight_h", "testHist_Weight_h", 100, 0, 500);

  for(Int_t evtIter = 0; evtIter < nEntries; evtIter++){
    ptHatChain_p->GetEntry(evtIter);

    testHist_NoWeight_p->Fill(ptHat_);

    for(Int_t hatIter = 0; hatIter < numCut; hatIter++){
      if(ptHat_ > ptHatCuts[hatIter] && ptHat_ < ptHatCuts[hatIter+1]){
	testHist_Weight_p->Fill(ptHat_, hatWeight[hatIter]);
	break;
      }
    }
  }

  TCanvas* quickCanv_p = new TCanvas("quickCanv", "quickCanv", 300*2, 350*1);
  quickCanv_p->Divide(2, 1, 0.0, 0.0);

  quickCanv_p->cd(1);
  gPad->SetLogy();
  testHist_NoWeight_p->DrawCopy("HIST");

  quickCanv_p->cd(2);
  gPad->SetLogy();
  testHist_Weight_p->DrawCopy("HIST");

  return;
}
