#include "TChain.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

TChain *akChain_p;
TChain *hiChain_p;
TChain *skChain_p;


void getLogBins(const Float_t lower, const Float_t higher, const Int_t nBins, Float_t bins[])
{
  Float_t logBins[nBins+1];
  bins[0] = lower;
  bins[nBins] = higher;

  logBins[0] = TMath::Log10(lower);
  logBins[nBins] = TMath::Log10(higher);

  Float_t interval = (logBins[nBins] - logBins[0])/nBins;

  for(Int_t iter = 1; iter < nBins; iter++){
    logBins[iter] = logBins[0] + iter*interval;
    bins[iter] = TMath::Power(10, logBins[iter]);
  }

  return;
}



void quickRun(const std::string fList, Int_t nBins = 20)
{
  std::string buffer;
  std::vector<std::string> listOfFiles;
  int nLines = 0;
  ifstream inFile(fList.data());

  std::cout << fList << std::endl;
  std::cout << inFile.is_open() << std::endl;

  if(!inFile.is_open()){
    std::cout << "Error opening file. Exiting." << std::endl;
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

  akChain_p = new TChain("akVs5CaloJetAnalyzer/t");
  hiChain_p = new TChain("hiEvtAnalyzer/HiTree");
  skChain_p = new TChain("skimanalysis/HltTree");


  for(Int_t iter = 0; iter < (Int_t)listOfFiles.size(); iter++){
    std::cout << listOfFiles.at(iter) << "." << std::endl;

    akChain_p->Add(listOfFiles[iter].c_str());
    hiChain_p->Add(listOfFiles[iter].c_str());
    skChain_p->Add(listOfFiles[iter].c_str());
  }

  akChain_p->AddFriend(hiChain_p);
  akChain_p->AddFriend(skChain_p);

  const Float_t ptLower = 20.;
  const Float_t ptHigher = 700.;
  Float_t ptBins[nBins+1];

  getLogBins(ptLower, ptHigher, nBins, ptBins);

  for(Int_t iter = 0; iter < nBins; iter++){
    std::cout << "iter, cuts: " << iter << ", " << ptBins[iter] << ", " << ptBins[iter+1] << std::endl;
  }

  return;
}

