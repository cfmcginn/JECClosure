#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include "TChain.h"
#include "TCut.h"

#include "TCanvas.h"
#include "TLine.h"

#include "commonSetup.h"

#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>

#include "commonUtility.h"
#include "TMath.h"
#include "TLatex.h"

#include <vector>
#include "cfmVectFunc.h"

TChain* getChain_p[3] = {0, 0, 0};

Int_t ptHatCuts_PYTHHI[10] = {15, 30, 50, 80, 120, 170, 220, 280, 370, 1000000};
Int_t ptHatCuts_PYTHPP[7] = {15, 30, 50, 80, 120, 170, 1000000};
Int_t ptHatCuts_PYTHHYD[9] = {15, 30, 50, 80, 120, 220, 280, 370, 1000000};

Float_t ptHatWeights_PYTHHYD[8] = {.611066, .0399951, .00243874, .000241009, .0000273228, .00000147976, .000000618337, .000000250369};

Int_t centCutArray[9] = {0, 10, 20, 40, 60, 80, 100, 140, 200};
Float_t centCutArray_F[9] = {0, 10, 20, 40, 60, 80, 100, 140, 200};

Int_t pcollisionEventSelection_;
Int_t pPAcollisionEventSelectionPA_;
Float_t vz_;
Float_t pthat_;
Int_t hiBin_;

Int_t nref_;
Float_t jtpt_[maxJets];
Float_t jteta_[maxJets];
Float_t rawpt_[maxJets];
Float_t refpt_[maxJets];
Float_t refeta_[maxJets];
Float_t refphi_[maxJets];
Int_t subid_[maxJets];

Int_t ngen_;
Float_t genpt_[maxJets];
Float_t geneta_[maxJets];
Float_t genphi_[maxJets];
Int_t genmatchindex_[maxJets];
Int_t gensubid_[maxJets];

const Float_t jtPtCut = 20.0;

const std::string PtEtaString[2] = {"Pt", "Eta"};
const std::string MeanResString[2] = {"Mean", "Res"};
const std::string PFCaloString[2] = {"PF", "Calo"};
const std::string PuVsString[3] = {"", "Pu", "Vs"};

void BookChain(Bool_t isPbPb)
{
  getChain_p[0]->SetBranchStatus("*", 0);

  if(isPbPb) getChain_p[0]->SetBranchStatus("pcollisionEventSelection", 1);
  if(!isPbPb) getChain_p[0]->SetBranchStatus("pPAcollisionEventSelectionPA", 1);

  getChain_p[0]->SetBranchStatus("vz", 1);
  getChain_p[0]->SetBranchStatus("pthat", 1);
  getChain_p[0]->SetBranchStatus("hiBin", 1);

  getChain_p[0]->SetBranchStatus("nref", 1);
  getChain_p[0]->SetBranchStatus("jtpt", 1);
  getChain_p[0]->SetBranchStatus("jteta", 1);
  getChain_p[0]->SetBranchStatus("rawpt", 1);
  getChain_p[0]->SetBranchStatus("refpt", 1);
  getChain_p[0]->SetBranchStatus("refeta", 1);
  getChain_p[0]->SetBranchStatus("refphi", 1);
  getChain_p[0]->SetBranchStatus("subid", 1);

  getChain_p[0]->SetBranchStatus("ngen", 1);
  getChain_p[0]->SetBranchStatus("genpt", 1);
  getChain_p[0]->SetBranchStatus("geneta", 1);
  getChain_p[0]->SetBranchStatus("genphi", 1);
  getChain_p[0]->SetBranchStatus("genmatchindex", 1);
  getChain_p[0]->SetBranchStatus("gensubid", 1);

  if(isPbPb) getChain_p[0]->SetBranchAddress("pcollisionEventSelection", &pcollisionEventSelection_);
  if(!isPbPb) getChain_p[0]->SetBranchAddress("pPAcollisionEventSelectionPA", &pPAcollisionEventSelectionPA_);

  getChain_p[0]->SetBranchAddress("vz", &vz_);
  getChain_p[0]->SetBranchAddress("pthat", &pthat_);
  getChain_p[0]->SetBranchAddress("hiBin", &hiBin_);

  getChain_p[0]->SetBranchAddress("nref", &nref_);
  getChain_p[0]->SetBranchAddress("jtpt", jtpt_);
  getChain_p[0]->SetBranchAddress("jteta", jteta_);
  getChain_p[0]->SetBranchAddress("rawpt", rawpt_);
  getChain_p[0]->SetBranchAddress("refpt", refpt_);
  getChain_p[0]->SetBranchAddress("refeta", refeta_);
  getChain_p[0]->SetBranchAddress("refphi", refphi_);
  getChain_p[0]->SetBranchAddress("subid", subid_);

  getChain_p[0]->SetBranchAddress("ngen", &ngen_);
  getChain_p[0]->SetBranchAddress("genpt", genpt_);
  getChain_p[0]->SetBranchAddress("geneta", geneta_);
  getChain_p[0]->SetBranchAddress("genphi", genphi_);
  getChain_p[0]->SetBranchAddress("genmatchindex", genmatchindex_);
  getChain_p[0]->SetBranchAddress("gensubid", gensubid_);

  return;
}

void GetChain(std::vector<std::string> inList, const std::string alg, Bool_t isPbPb)
{
  getChain_p[0] = new TChain(Form("%sJetAnalyzer/t", alg.c_str()));
  getChain_p[1] = new TChain("skimanalysis/HltTree");
  getChain_p[2] = new TChain("hiEvtAnalyzer/HiTree");

  for(Int_t iter = 0; iter < (Int_t)(inList.size()); iter++){
    getChain_p[0]->Add(inList[iter].c_str());
    getChain_p[1]->Add(inList[iter].c_str());
    getChain_p[2]->Add(inList[iter].c_str());
  }
  getChain_p[0]->AddFriend(getChain_p[1]);
  getChain_p[0]->AddFriend(getChain_p[2]);

  BookChain(isPbPb);

  return;
}


void CleanChain()
{
  if(getChain_p[0] != 0){
    delete getChain_p[0];
    getChain_p[0] = 0;
  }

  if(getChain_p[1] != 0){
    delete getChain_p[1];
    getChain_p[1] = 0;
  }

  if(getChain_p[2] != 0){
    delete getChain_p[2];
    getChain_p[2] = 0;
  }
  return;
}


Float_t getHatWeight(Float_t inHat, Bool_t isPbPb, Bool_t isHITrk)
{
  if(isPbPb){
    for(Int_t iter = 0; iter < 8; iter++){
      if(inHat > ptHatCuts_PYTHHYD[iter] && inHat < ptHatCuts_PYTHHYD[iter+1])
	return ptHatWeights_PYTHHYD[iter];
    }
  }

  std::cout << inHat << std::endl;
  std::cout << "No weight assigned; check for error." << std::endl;
  return 0;
}


Bool_t jtPasses(Float_t inPt, Float_t inHat, Bool_t isPbPb, Bool_t isHITrk)
{
  if(isPbPb){
    for(Int_t iter = 0; iter < 8; iter++){
      if(inHat >= ptHatCuts_PYTHHYD[iter] && inHat < ptHatCuts_PYTHHYD[iter+1]){
        if(inPt >= ptHatCuts_PYTHHYD[iter] && inPt < ptHatCuts_PYTHHYD[iter+1])
	  return true;
	else
	  return false;
      }
    }
  }
  else if(isHITrk){
    for(Int_t iter = 0; iter < 9; iter++){
      if(inHat >= ptHatCuts_PYTHHI[iter] && inHat < ptHatCuts_PYTHHI[iter+1]){
        if(inPt >= ptHatCuts_PYTHHI[iter] && inPt < ptHatCuts_PYTHHI[iter+1])
          return true;
        else
          return false;
      }
    }
  }
  else{
    for(Int_t iter = 0; iter < 6; iter++){
      if(inHat >= ptHatCuts_PYTHPP[iter] && inHat < ptHatCuts_PYTHPP[iter+1]){
        if(inPt >= ptHatCuts_PYTHPP[iter] && inPt < ptHatCuts_PYTHPP[iter+1])
          return true;
        else
          return false;
      }
    }
  }

  std::cout << inHat << std::endl;
  std::cout << "No Bool assigned; check for error." << std::endl;
  return false;
}


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


void getEtaBins(const Float_t lower, const Float_t higher, const Int_t nBins, Float_t bins[])
{
  Float_t interval = (higher - lower)/nBins;

  for(Int_t iter = 0; iter < nBins+1; iter++){
    bins[iter] = lower + interval*iter;
  }
  bins[0] = -1.99999999999;
  bins[nBins] = 1.99999999999;

  return;
}

void getPhiBins(const Float_t lower, const Float_t higher, const Int_t nBins, Float_t bins[])
{
  Float_t interval = (higher - lower)/nBins;

  for(Int_t iter = 0; iter < nBins+1; iter++){
    bins[iter] = lower + interval*iter;
  }
  bins[0] = -TMath::Pi();
  bins[nBins] = TMath::Pi();

  return;
}


Int_t posSearch(Float_t val, Int_t arrSize, Float_t arr[])
{
  Int_t pos = arrSize/2;
  Int_t low = 0;
  Int_t high = arrSize-1;

  Int_t iter = 0;

  if(val > arr[high]) return high - 1;

  while(val < arr[pos] || val >= arr[pos+1]){
    iter++;

    if(val < arr[pos]) high = pos;
    else low = pos;
    pos = (high+low)/2;
  }

  return pos;
}


Int_t posSearch(Float_t val, Int_t arrSize, Int_t arr[])
{
  Int_t pos = arrSize/2;
  Int_t low = 0;
  Int_t high = arrSize-1;

  Int_t iter = 0;

  if(val > arr[high]) return high - 1;

  while(val < arr[pos] || val >= arr[pos+1]){
    iter++;

    if(val < arr[pos]) high = pos;
    else low = pos;
    pos = (high+low)/2;
  }

  return pos;
}


std::string getCentString(Bool_t isPbPb, Int_t centLow, Int_t centHi)
{
  if(isPbPb) return Form("%d%d", (Int_t)(centLow*.5), (Int_t)((centHi)*.5));
  else return "PP";
}


void InitHist(TH1F* inHist_p, const std::string PtEta, const std::string MeanRes)
{
  handsomeTH1(inHist_p);
  inHist_p->GetXaxis()->SetTitleOffset(.75);
  if(!strcmp(PtEta.c_str(), "Pt")) inHist_p->SetXTitle("p_{T}^{gen}");
  else if(!strcmp(PtEta.c_str(), "Eta")) inHist_p->SetXTitle("#eta_{gen}");
  else inHist_p->SetXTitle("#phi_{gen}");

  if(!strcmp(MeanRes.c_str(), "Mean")) inHist_p->SetYTitle("<p_{T}^{jet}/p_{T}^{gen}>");
  else if(!strcmp(MeanRes.c_str(), "Res")) inHist_p->SetYTitle("#sigma_{reco/gen}");
  else inHist_p->SetYTitle("Efficiency");

  return;
}


void CleanTH1(TH1F* cleanHist_p)
{
  delete cleanHist_p;
  cleanHist_p = 0;
}


void FitGauss(TH1F* inHist_p, Float_t &mean, Float_t &meanError, Float_t &res, Float_t &resError)
{
  inHist_p->Fit("gaus", "Q L", "");

  mean = inHist_p->GetFunction("gaus")->GetParameter(1);
  meanError = inHist_p->GetFunction("gaus")->GetParError(1);
  res = inHist_p->GetFunction("gaus")->GetParameter(2);
  resError = inHist_p->GetFunction("gaus")->GetParError(2);

  Float_t prob = inHist_p->GetFunction("gaus")->GetProb();

  //  if(TMath::Abs(1.00 - mean) < .01) return;

  Int_t meanBin = inHist_p->FindBin(mean);
  Float_t meanRMS = -1;  
  Float_t total = inHist_p->Integral();

  for(Int_t iter = 0; iter < inHist_p->GetNbinsX(); iter++){
    Int_t lowBound = 0;
    if(meanBin - iter > 0) lowBound = meanBin - iter; 

    if(inHist_p->Integral(lowBound, meanBin + iter)/total > .95 || lowBound == 0 || inHist_p->GetBinContent(lowBound) < .01){
      meanRMS = inHist_p->GetBinCenter(meanBin + iter) - inHist_p->GetBinCenter(meanBin);
      break;
    }
  }
  
  double minPt = inHist_p->GetBinCenter(meanBin) - meanRMS;
  double maxPt = inHist_p->GetBinCenter(meanBin) + meanRMS;

  minPt = std::max(std::max(minPt, 0.0), inHist_p->GetXaxis()->GetXmin());
  maxPt = std::min(maxPt, inHist_p->GetXaxis()->GetXmax());

  inHist_p->Fit("gaus", "Q L", "", minPt, maxPt);
  
  if(TMath::Abs(1.00 - mean) < TMath::Abs(1.00 - inHist_p->GetFunction("gaus")->GetParameter(1)) && prob > 0.0001){
    inHist_p->Fit("gaus", "Q L", "");
    return;
  }

  mean = inHist_p->GetFunction("gaus")->GetParameter(1);
  meanError = inHist_p->GetFunction("gaus")->GetParError(1);
  res = inHist_p->GetFunction("gaus")->GetParameter(2);
  resError = inHist_p->GetFunction("gaus")->GetParError(2);

  return;
}


void makeDiJetIniHist(std::vector<std::string> inList, const std::string outName, const std::string alg = "akVs3Calo", Bool_t isPbPb = false, Bool_t isHITrk = false)
{
  TH1::SetDefaultSumw2();

  std::cout << "Init: " << alg << std::endl;

  GetChain(inList, alg, isPbPb);

  const Int_t nPtBins = 20;
  const Float_t ptLower = 20.;
  const Float_t ptHigher = 700.;

  const Int_t nEtaBins = 20;
  const Float_t etaLower = -2.0;
  const Float_t etaHigher = 2.0;

  const Int_t nPhiBins = 20;
  const Float_t phiLower = -TMath::Pi();
  const Float_t phiHigher = TMath::Pi();

  Int_t nCentBins;

  if(isPbPb) nCentBins = 8;
  else nCentBins = 1;

  std::string centString[nCentBins];

  if(isPbPb){
    for(Int_t centIter = 0; centIter < nCentBins; centIter++){
      centString[centIter] = getCentString(isPbPb, centCutArray[centIter], centCutArray[centIter+1]);
    }
  }
  else centString[0] = getCentString(isPbPb, 0, 0);

  //Edit Here

  Float_t ptBins[nPtBins+1];
  Float_t etaBins[nEtaBins+1];
  Float_t phiBins[nPhiBins+1];

  getLogBins(ptLower, ptHigher, nPtBins, ptBins);
  getEtaBins(etaLower, etaHigher, nEtaBins, etaBins);
  getPhiBins(phiLower, phiHigher, nPhiBins, phiBins);

  Float_t effBinning[3] = {-.5, .5, 1.5};

  TH1F* pthatNoWeight_p = new TH1F("pthatNoWeight_h", "pthatNoWeight_h", 50, 0, 500);
  TH1F* pthatWeight_p = new TH1F("pthatWeight_h", "pthatWeight_h", 50, 0, 500);
  TH1F* temp_jetOverGen_Pt_p[nCentBins][nPtBins];
  TH1F* temp_jetOverGen_Eta_p[nCentBins][nEtaBins];
  TH1F* temp_jetOverGen_Phi_p[nCentBins][nPhiBins];
  TH1F* temp_genEff_Pt_p[nCentBins][nPtBins];
  TH1F* temp_genEff_Eta_p[nCentBins][nEtaBins];
  TH1F* temp_genEff_Phi_p[nCentBins][nPhiBins];

  for(Int_t centIter = 0; centIter < nCentBins; centIter++){
    for(Int_t iter = 0; iter < nPtBins; iter++){
      temp_jetOverGen_Pt_p[centIter][iter] = new TH1F(Form("%s_tempJetHist_Pt_%s_%d", alg.c_str(), centString[centIter].c_str(), iter), Form("%s_tempJetHist_Pt_%s_%d", alg.c_str(), centString[centIter].c_str(), iter), 150, -0.01, 2.99);
      temp_genEff_Pt_p[centIter][iter] = new TH1F(Form("%s_tempGenEffHist_Pt_%s_%d", alg.c_str(), centString[centIter].c_str(), iter), Form("%s_tempGenEffHist_Pt_%s_%d", alg.c_str(), centString[centIter].c_str(), iter), 2, effBinning);      
    }

    for(Int_t iter = 0; iter < nEtaBins; iter++){
      temp_jetOverGen_Eta_p[centIter][iter] = new TH1F(Form("%s_tempJetHist_Eta_%s_%d", alg.c_str(), centString[centIter].c_str(), iter), Form("%s_tempJetHist_Eta_%s_%d", alg.c_str(), centString[centIter].c_str(), iter), 150, -0.01, 2.99);
      temp_genEff_Eta_p[centIter][iter] = new TH1F(Form("%s_tempGenEffHist_Eta_%s_%d", alg.c_str(), centString[centIter].c_str(), iter), Form("%s_tempGenEffHist_Eta_%s_%d", alg.c_str(), centString[centIter].c_str(), iter), 2, effBinning);
    }

    for(Int_t iter = 0; iter < nPhiBins; iter++){
      temp_jetOverGen_Phi_p[centIter][iter] = new TH1F(Form("%s_tempJetHist_Phi_%s_%d", alg.c_str(), centString[centIter].c_str(), iter), Form("%s_tempJetHist_Phi_%s_%d", alg.c_str(), centString[centIter].c_str(), iter), 150, -0.01, 2.99);
      temp_genEff_Phi_p[centIter][iter] = new TH1F(Form("%s_tempGenEffHist_Phi_%s_%d", alg.c_str(), centString[centIter].c_str(), iter), Form("%s_tempGenEffHist_Phi_%s_%d", alg.c_str(), centString[centIter].c_str(), iter), 2, effBinning);
    }
  }

  TH1F* jetOverGenMean_Pt_p[nCentBins];
  TH1F* jetOverGenMean_Eta_p[nCentBins];
  TH1F* jetOverGenMean_Phi_p[nCentBins];

  TH1F* jetOverGenRes_Pt_p[nCentBins];
  TH1F* jetOverGenRes_Eta_p[nCentBins];
  TH1F* jetOverGenRes_Phi_p[nCentBins];

  TH1F* genEff_Pt_p[nCentBins];
  TH1F* genEff_Eta_p[nCentBins];
  TH1F* genEff_Phi_p[nCentBins];

  for(Int_t centIter = 0; centIter < nCentBins; centIter++){
    jetOverGenMean_Pt_p[centIter] = new TH1F(Form("%s_Mean_Pt_%s_h", alg.c_str(), centString[centIter].c_str()), Form("%s_Mean_Pt_%s_h", alg.c_str(), centString[centIter].c_str()), nPtBins, ptBins);
    jetOverGenMean_Eta_p[centIter] = new TH1F(Form("%s_Mean_Eta_%s_h", alg.c_str(), centString[centIter].c_str()), Form("%s_Mean_Eta_%s_h", alg.c_str(), centString[centIter].c_str()), nEtaBins, etaBins);
    jetOverGenMean_Phi_p[centIter] = new TH1F(Form("%s_Mean_Phi_%s_h", alg.c_str(), centString[centIter].c_str()), Form("%s_Mean_Phi_%s_h", alg.c_str(), centString[centIter].c_str()), nPhiBins, phiBins);

    jetOverGenRes_Pt_p[centIter] = new TH1F(Form("%s_Res_Pt_%s_h", alg.c_str(), centString[centIter].c_str()), Form("%s_Res_Pt_%s_h", alg.c_str(), centString[centIter].c_str()), nPtBins, ptBins);
    jetOverGenRes_Eta_p[centIter] = new TH1F(Form("%s_Res_Eta_%s_h", alg.c_str(), centString[centIter].c_str()), Form("%s_Res_Eta_%s_h", alg.c_str(), centString[centIter].c_str()), nEtaBins, etaBins);
    jetOverGenRes_Phi_p[centIter] = new TH1F(Form("%s_Res_Phi_%s_h", alg.c_str(), centString[centIter].c_str()), Form("%s_Res_Phi_%s_h", alg.c_str(), centString[centIter].c_str()), nPhiBins, phiBins);

    genEff_Pt_p[centIter] = new TH1F(Form("%s_Eff_Pt_%s_h", alg.c_str(), centString[centIter].c_str()), Form("%s_Eff_Pt_%s_h", alg.c_str(), centString[centIter].c_str()), nPtBins, ptBins);
    genEff_Eta_p[centIter] = new TH1F(Form("%s_Eff_Eta_%s_h", alg.c_str(), centString[centIter].c_str()), Form("%s_Eff_Eta_%s_h", alg.c_str(), centString[centIter].c_str()), nEtaBins, etaBins);
    genEff_Phi_p[centIter] = new TH1F(Form("%s_Eff_Phi_%s_h", alg.c_str(), centString[centIter].c_str()), Form("%s_Eff_Phi_%s_h", alg.c_str(), centString[centIter].c_str()), nPhiBins, phiBins);

    InitHist(jetOverGenMean_Pt_p[centIter], "Pt", "Mean");
    InitHist(jetOverGenMean_Eta_p[centIter], "Eta", "Mean");
    InitHist(jetOverGenMean_Phi_p[centIter], "Phi", "Mean");

    InitHist(jetOverGenRes_Pt_p[centIter], "Pt", "Res");
    InitHist(jetOverGenRes_Eta_p[centIter], "Eta", "Res");
    InitHist(jetOverGenRes_Phi_p[centIter], "Phi", "Res");

    InitHist(genEff_Pt_p[centIter], "Pt", "Eff");
    InitHist(genEff_Eta_p[centIter], "Eta", "Eff");
    InitHist(genEff_Phi_p[centIter], "Phi", "Eff");
  }

  TH1F* genRefRawJtPt_p[4];
  const std::string genRefRawJt[4] = {"Gen", "Ref", "Raw", "Jt"};

  for(Int_t iter = 0; iter < 4; iter++){
    genRefRawJtPt_p[iter] = new TH1F(Form("%s_%sPt_h", alg.c_str(), genRefRawJt[iter].c_str()), Form("%s_%sPt_h", alg.c_str(), genRefRawJt[iter].c_str()), 10, 0, 30);
  }

  Int_t nEntries = getChain_p[0]->GetEntries();
  std::cout << "# Entries: " << nEntries << std::endl;

  for(Int_t jEntry = 0; jEntry < nEntries; jEntry++){
    getChain_p[0]->GetEntry(jEntry);

    if(jEntry%10000 == 0) std::cout << jEntry << std::endl;

    if(!pcollisionEventSelection_ && isPbPb) continue;

    if(!pPAcollisionEventSelectionPA_ && !isPbPb) continue;

    if(TMath::Abs(vz_) > 15) continue;

    Int_t centPos = 0;
    if(isPbPb) centPos = posSearch(hiBin_, 9, centCutArray);

    Float_t weight = getHatWeight(pthat_, isPbPb, isHITrk);

    pthatNoWeight_p->Fill(pthat_);
    pthatWeight_p->Fill(pthat_, weight);

    for(Int_t jtIter = 0; jtIter < nref_; jtIter++){
      if(!jtPasses(refpt_[jtIter], pthat_, isPbPb, isHITrk)) continue;

      if(TMath::Abs(refeta_[jtIter]) >= 2.0) continue;
      if(subid_[jtIter] != 0) continue;

      if(refpt_[jtIter] > 0 && refpt_[jtIter] < 30) genRefRawJtPt_p[1]->Fill(refpt_[jtIter]);
      if(rawpt_[jtIter] > 0 && rawpt_[jtIter] < 30) genRefRawJtPt_p[2]->Fill(rawpt_[jtIter]);
      if(jtpt_[jtIter] > 0 && jtpt_[jtIter] < 30) genRefRawJtPt_p[3]->Fill(jtpt_[jtIter]);

      if(refpt_[jtIter] < jtPtCut) continue;

      Int_t ptPos = posSearch(refpt_[jtIter], nPtBins+1, ptBins);
      temp_jetOverGen_Pt_p[centPos][ptPos]->Fill(jtpt_[jtIter]/refpt_[jtIter]);

      if(refpt_[jtIter] > 50){
	Int_t etaPos = posSearch(refeta_[jtIter], nEtaBins+1, etaBins);
	temp_jetOverGen_Eta_p[centPos][etaPos]->Fill(jtpt_[jtIter]/refpt_[jtIter]);

	Int_t phiPos = posSearch(refphi_[jtIter], nPhiBins+1, phiBins);
	temp_jetOverGen_Phi_p[centPos][phiPos]->Fill(jtpt_[jtIter]/refpt_[jtIter]);
      }
    }

    for(Int_t genIter = 0; genIter < ngen_; genIter++){
      if(!jtPasses(genpt_[genIter], pthat_, isPbPb, isHITrk)) continue;

      if(TMath::Abs(geneta_[genIter]) >= 2.0) continue;
      if(gensubid_[genIter] != 0) continue;

      if(genpt_[genIter] > 0 && genpt_[genIter] < 30) genRefRawJtPt_p[0]->Fill(genpt_[genIter]);

      if(genpt_[genIter] < jtPtCut) continue;

      Int_t ptPos = posSearch(genpt_[genIter], nPtBins+1, ptBins);
      if(genmatchindex_[genIter] >= 0) temp_genEff_Pt_p[centPos][ptPos]->Fill(1.0);
      else temp_genEff_Pt_p[centPos][ptPos]->Fill(0.0);
      
      if(genpt_[genIter] > 50){
	Int_t etaPos = posSearch(geneta_[genIter], nEtaBins+1, etaBins);
	if(genmatchindex_[genIter] >= 0) temp_genEff_Eta_p[centPos][etaPos]->Fill(1.0);
	else temp_genEff_Eta_p[centPos][etaPos]->Fill(0.0);

	Int_t phiPos = posSearch(genphi_[genIter], nPhiBins+1, phiBins);
	if(genmatchindex_[genIter] >= 0) temp_genEff_Phi_p[centPos][phiPos]->Fill(1.0);
	else temp_genEff_Phi_p[centPos][phiPos]->Fill(0.0);
      }
    }
  }

  for(Int_t centIter = 0; centIter < nCentBins; centIter++){
    for(Int_t rawIter = 0; rawIter < nPtBins; rawIter++){

      Float_t mean = 0;
      Float_t meanError = 0;

      Float_t res = 0;
      Float_t resError = 0;
      
      FitGauss(temp_jetOverGen_Pt_p[centIter][rawIter], mean, meanError, res, resError);
      
      jetOverGenMean_Pt_p[centIter]->SetBinContent(rawIter+1, mean);
      jetOverGenRes_Pt_p[centIter]->SetBinContent(rawIter+1, res);
      jetOverGenMean_Pt_p[centIter]->SetBinError(rawIter+1, meanError);
      jetOverGenRes_Pt_p[centIter]->SetBinError(rawIter+1, resError);	       


      genEff_Pt_p[centIter]->SetBinContent(rawIter+1, temp_genEff_Pt_p[centIter][rawIter]->GetMean());
      genEff_Pt_p[centIter]->SetBinError(rawIter+1, temp_genEff_Pt_p[centIter][rawIter]->GetMeanError());
    }

    for(Int_t rawIter = 0; rawIter < nEtaBins; rawIter++){

      Float_t mean = 0;
      Float_t meanError = 0;

      Float_t res = 0;
      Float_t resError = 0;

      FitGauss(temp_jetOverGen_Eta_p[centIter][rawIter], mean, meanError, res, resError);      

      jetOverGenMean_Eta_p[centIter]->SetBinContent(rawIter+1, mean);
      jetOverGenRes_Eta_p[centIter]->SetBinContent(rawIter+1, res);	
      jetOverGenMean_Eta_p[centIter]->SetBinError(rawIter+1, meanError);
      jetOverGenRes_Eta_p[centIter]->SetBinError(rawIter+1, resError);
      
      genEff_Eta_p[centIter]->SetBinContent(rawIter+1, temp_genEff_Eta_p[centIter][rawIter]->GetMean());
      genEff_Eta_p[centIter]->SetBinError(rawIter+1, temp_genEff_Eta_p[centIter][rawIter]->GetMeanError());
    }

    for(Int_t rawIter = 0; rawIter < nPhiBins; rawIter++){

      Float_t mean = 0;
      Float_t meanError = 0;

      Float_t res = 0;
      Float_t resError = 0;

      FitGauss(temp_jetOverGen_Phi_p[centIter][rawIter], mean, meanError, res, resError);      

      jetOverGenMean_Phi_p[centIter]->SetBinContent(rawIter+1, mean);
      jetOverGenRes_Phi_p[centIter]->SetBinContent(rawIter+1, res);	
      jetOverGenMean_Phi_p[centIter]->SetBinError(rawIter+1, meanError);
      jetOverGenRes_Phi_p[centIter]->SetBinError(rawIter+1, resError);
      
      genEff_Phi_p[centIter]->SetBinContent(rawIter+1, temp_genEff_Phi_p[centIter][rawIter]->GetMean());
      genEff_Phi_p[centIter]->SetBinError(rawIter+1, temp_genEff_Phi_p[centIter][rawIter]->GetMeanError());
    }
  }

  TFile* out = new TFile(outName.c_str(), "UPDATE");
  std::cout << outName << std::endl;
  if(!strcmp(alg.c_str(), "akVs3PF")){
    pthatWeight_p->Scale(1./pthatWeight_p->Integral());
    handsomeTH1(pthatWeight_p);
    pthatWeight_p->GetXaxis()->SetTitleOffset(1.15);
    pthatWeight_p->SetXTitle("p_{T}^{#hat}");
    pthatWeight_p->GetXaxis()->SetTitleOffset(1.15);
    pthatWeight_p->SetYTitle("Event Fraction");
    pthatWeight_p->Write("", TObject::kOverwrite);

    pthatNoWeight_p->Scale(1./pthatNoWeight_p->Integral());
    handsomeTH1(pthatNoWeight_p);
    pthatNoWeight_p->GetXaxis()->SetTitleOffset(1.15);
    pthatNoWeight_p->SetXTitle("p_{T}^{#hat}");
    pthatNoWeight_p->GetXaxis()->SetTitleOffset(1.15);
    pthatNoWeight_p->SetYTitle("Event Fraction");
    pthatNoWeight_p->Write("", TObject::kOverwrite);
  }
  for(Int_t centIter = 0; centIter < nCentBins; centIter++){
    jetOverGenMean_Pt_p[centIter]->Write("", TObject::kOverwrite);
    jetOverGenMean_Eta_p[centIter]->Write("", TObject::kOverwrite);
    jetOverGenMean_Phi_p[centIter]->Write("", TObject::kOverwrite);

    jetOverGenRes_Pt_p[centIter]->Write("", TObject::kOverwrite);
    jetOverGenRes_Eta_p[centIter]->Write("", TObject::kOverwrite);
    jetOverGenRes_Phi_p[centIter]->Write("", TObject::kOverwrite);

    genEff_Pt_p[centIter]->Write("", TObject::kOverwrite);
    genEff_Eta_p[centIter]->Write("", TObject::kOverwrite);
    genEff_Phi_p[centIter]->Write("", TObject::kOverwrite);

    for(Int_t iter = 0; iter < nPtBins; iter++){
      temp_jetOverGen_Pt_p[centIter][iter]->Write("", TObject::kOverwrite);
      //      temp_jetOverGen_Eta_p[centIter][iter]->Write("", TObject::kOverwrite);
      //      temp_genEff_Pt_p[centIter][iter]->Write("", TObject::kOverwrite);
      //      temp_genEff_Eta_p[centIter][iter]->Write("", TObject::kOverwrite);
    }
  }
  for(Int_t iter = 0; iter < 4; iter++){
    genRefRawJtPt_p[iter]->Write("", TObject::kOverwrite);
  }
  out->Close();
  delete out;

  for(Int_t centIter = 0; centIter < nCentBins; centIter++){
    CleanTH1(jetOverGenMean_Pt_p[centIter]);
    CleanTH1(jetOverGenMean_Eta_p[centIter]);
    CleanTH1(jetOverGenMean_Phi_p[centIter]);

    CleanTH1(jetOverGenRes_Pt_p[centIter]);
    CleanTH1(jetOverGenRes_Eta_p[centIter]);
    CleanTH1(jetOverGenRes_Phi_p[centIter]);

    CleanTH1(genEff_Pt_p[centIter]);
    CleanTH1(genEff_Eta_p[centIter]);
    CleanTH1(genEff_Phi_p[centIter]);

    for(Int_t rawIter = 0; rawIter < nPtBins; rawIter++){
      CleanTH1(temp_jetOverGen_Pt_p[centIter][rawIter]);
      CleanTH1(temp_genEff_Pt_p[centIter][rawIter]);
    }

    for(Int_t rawIter = 0; rawIter < nEtaBins; rawIter++){
      CleanTH1(temp_jetOverGen_Eta_p[centIter][rawIter]);
      CleanTH1(temp_genEff_Eta_p[centIter][rawIter]);
    }

    for(Int_t rawIter = 0; rawIter < nPhiBins; rawIter++){
      CleanTH1(temp_jetOverGen_Phi_p[centIter][rawIter]);
      CleanTH1(temp_genEff_Phi_p[centIter][rawIter]);
    }
  }

  for(Int_t iter = 0; iter < 4; iter++){
    CleanTH1(genRefRawJtPt_p[iter]);
  }

  CleanTH1(pthatWeight_p);
  CleanTH1(pthatNoWeight_p);

  CleanChain();

  return;
}

int runMakeDiJetIniHist(std::string fList = "", const char* outFileName = "raw_rawOverGen", Bool_t isPbPb = false, Bool_t isHITrk = false)
{
  TH1::SetDefaultSumw2();

  std::string buffer;
  std::vector<std::string> listOfFiles;
  int nLines = 0;
  ifstream inFile(fList.data());

  std::cout << fList << std::endl;
  std::cout << inFile.is_open() << std::endl;

  if(!inFile.is_open()){
    std::cout << "Error opening file. Exiting." << std::endl;
    return 1;
  }
  else{
    while(true){
      inFile >> buffer;
      if(inFile.eof()) break;
      listOfFiles.push_back(buffer);
      nLines++;
    }
  }
  
  for(Int_t iter = 0; iter < (Int_t)listOfFiles.size(); iter++){
    std::cout << listOfFiles.at(iter) << "." << std::endl;
  }
  if(isPbPb){
    for(Int_t iter = 0; iter < 9; iter++){
      std::cout << iter << ": " << ptHatCuts_PYTHHYD[iter] << std::endl;
    }
  }
  else if(isHITrk){
    for(Int_t iter = 0; iter < 10; iter++){
      std::cout << iter << ": " << ptHatCuts_PYTHHI[iter] << std::endl;
    }
  }
  else{
    for(Int_t iter = 0; iter < 7; iter++){
      std::cout << iter << ": " << ptHatCuts_PYTHPP[iter] << std::endl;
    }
  }

  Int_t pfCaloStart = 0;
  Int_t numStart = 1;
  Int_t numEnd = 5;

  if(isHITrk){
    pfCaloStart = 1;
    numStart = 2;
    numEnd = 5;
  }

  for(Int_t numIter = numStart; numIter < numEnd; numIter++){
    for(Int_t pfCaloIter = pfCaloStart; pfCaloIter < 2; pfCaloIter++){

      if(!isPbPb) makeDiJetIniHist(listOfFiles, Form("%s.root", outFileName), Form("ak%d%s", numIter+1, PFCaloString[pfCaloIter].c_str()), isPbPb, isHITrk);

      if(isPbPb || isHITrk || numIter != 5)
	makeDiJetIniHist(listOfFiles, Form("%s.root", outFileName), Form("akVs%d%s", numIter+1, PFCaloString[pfCaloIter].c_str()), isPbPb, isHITrk);
      
      //      makeDiJetIniHist(listOfFiles, Form("%s.root", outFileName), Form("akPu%d%s", numIter+1, PFCaloString[pfCaloIter].c_str()), isPbPb, isHITrk);
      
    }    
  }

  return(0);
}


int main(int argc, char* argv[])
{
  if(argc != 5){
    std::cout << "Usage: runMakeDiJetIniHist <inputList> <outFileName> <isPbPb> <isHiTrk>" << std::endl;
    return 1;
  }

  int rStatus = -1;

  rStatus = runMakeDiJetIniHist(argv[1], argv[2], (Bool_t)(atoi(argv[3])), (Bool_t)(atoi(argv[4])));

  return rStatus;
}
