#include "TFile.h"
#include "TH1F.h"

#include "TCanvas.h"
#include "TLine.h"

#include <iostream>
#include <fstream>
#include <string>

#include "commonUtility.h"
#include "TMath.h"
#include "TLatex.h"

Int_t centCutArray[9] = {0, 10, 20, 40, 60, 80, 100, 140, 200};
Float_t centCutArray_F[9] = {0, 10, 20, 40, 60, 80, 100, 140, 200};

const std::string PtEtaString[3] = {"Pt", "Eta", "Phi"};
const std::string MeanResString[4] = {"Mean", "Res", "Eff", "Fake"};
const std::string PFCaloString[2] = {"PF", "Calo"};
const std::string PuVsString[3] = {"", "Pu", "Vs"};
const std::string QGString[3] = {"", "_Q", "_G"};

std::string getCentString(Bool_t isPbPb, Int_t centLow, Int_t centHi)
{
  if(isPbPb) return Form("%d%d", (Int_t)(centLow*.5), (Int_t)((centHi)*.5));
  else return "PP";
}


void CleanTH1(TH1F* cleanHist_p)
{
  delete cleanHist_p;
  cleanHist_p = 0;
}


void drawLine(const std::string MeanRes, const std::string PtEta)
{
  if(!strcmp(MeanRes.c_str(), "Mean") || !strcmp(MeanRes.c_str(), "Eff")){
    TLine* oneLine_p;
    if(!strcmp(PtEta.c_str(), "Pt")) oneLine_p = new TLine(20.0, 1.0, 700.0, 1.0);
    else if(!strcmp(PtEta.c_str(), "Eta")) oneLine_p = new TLine(-2.0, 1.0, 2.0, 1.0);
    else oneLine_p = new TLine(-TMath::Pi(), 1.0, TMath::Pi(), 1.0);

    oneLine_p->SetLineColor(1);
    oneLine_p->SetLineStyle(2);
    oneLine_p->Draw("SAME");

    if(!strcmp(MeanRes.c_str(), "Mean")){
      TLine* oneUpLine_p;
      if(!strcmp(PtEta.c_str(), "Pt")) oneUpLine_p = new TLine(20.0, 1.01, 700.0, 1.01);
      else if(!strcmp(PtEta.c_str(), "Eta")) oneUpLine_p = new TLine(-2.0, 1.01, 2.0, 1.01);
      else oneUpLine_p = new TLine(-TMath::Pi(), 1.01, TMath::Pi(), 1.01);

      oneUpLine_p->SetLineColor(1);
      oneUpLine_p->SetLineStyle(2);
      oneUpLine_p->Draw("SAME");
      
      TLine* oneDownLine_p;
      if(!strcmp(PtEta.c_str(), "Pt")) oneDownLine_p = new TLine(20.0, 0.99, 700.0, 0.99);
      else if(!strcmp(PtEta.c_str(), "Eta")) oneDownLine_p = new TLine(-2.0, 0.99, 2.0, 0.99);
      else oneDownLine_p = new TLine(-TMath::Pi(), 0.99, TMath::Pi(), 0.99);


      oneDownLine_p->SetLineColor(1);
      oneDownLine_p->SetLineStyle(2);
      oneDownLine_p->Draw("SAME");
    }
  }    

  if(!strcmp(PtEta.c_str(), "Pt")){
    TLine* fiftyLine_p;
    if(!strcmp(MeanRes.c_str(), "Mean")) fiftyLine_p = new TLine(50.0, 0.90, 50.0, 1.10);
    else if(!strcmp(MeanRes.c_str(), "Eff")) fiftyLine_p = new TLine(50.0, 0.80, 50.0, 1.10);
    else fiftyLine_p = new TLine(50.0, 0.00, 50.0, 0.50);
    
    fiftyLine_p->SetLineColor(1);
    fiftyLine_p->SetLineStyle(2);
    fiftyLine_p->Draw("SAME");
  }

  return;
}


void plotInitHist(TH1F* inHist_p)
{
  inHist_p->GetXaxis()->SetTitleFont(43);
  inHist_p->GetXaxis()->SetLabelFont(43);

  inHist_p->GetYaxis()->SetTitleFont(43);
  inHist_p->GetYaxis()->SetLabelFont(43);

  inHist_p->GetXaxis()->SetTitleSize(24);
  inHist_p->GetXaxis()->SetTitleOffset(1.8);
  inHist_p->GetXaxis()->SetLabelSize(16);
  inHist_p->GetYaxis()->SetTitleSize(24);
  inHist_p->GetYaxis()->SetTitleOffset(2.6);
  inHist_p->GetYaxis()->SetLabelSize(16);
}


Float_t getHistMax(const std::string MeanRes, const std::string PtEta)
{
  if(!strcmp(MeanRes.c_str(), "Mean")) return 1.09999;
  else if(!strcmp(MeanRes.c_str(), "Eff") && !strcmp(PtEta.c_str(), "Pt")) return 1.09999;
  else if(!strcmp(MeanRes.c_str(), "Eff") && !strcmp(PtEta.c_str(), "Eta")) return 1.01999;
  else return 0.34999;
}


Float_t getHistMin(const std::string MeanRes, const std::string PtEta)
{
  if(!strcmp(MeanRes.c_str(), "Mean")) return 0.90001;
  else if(!strcmp(MeanRes.c_str(), "Eff") && !strcmp(PtEta.c_str(), "Pt")) return .80001;
  else if(!strcmp(MeanRes.c_str(), "Eff") && !strcmp(PtEta.c_str(), "Eta")) return .98001;
  else return 0.00001;
}


void drawGenHist(TH1F* inHist_p, const std::string MeanRes, const std::string PtEta)
{
  if(!strcmp(PtEta.c_str(), "Pt")) gPad->SetLogx();
  else if(!strcmp(PtEta.c_str(), "Phi")) inHist_p->SetXTitle("#phi_{gen}");
  inHist_p->DrawCopy("P");
  drawLine(MeanRes, PtEta);
  return;
}


void plotDiJetIniHist_PYTH(const std::string histFileName, const std::string VsPu, const std::string PFCalo, const std::string MeanRes, const std::string qg, const std::string PtEta, Bool_t isFrag = false, Bool_t isRes = false)
{
  TH1::SetDefaultSumw2();
  TFile *f = new TFile(histFileName.c_str(), "UPDATE");

  Float_t max = getHistMax(MeanRes, PtEta);
  Float_t min = getHistMin(MeanRes, PtEta);

  TH1F* jecHist_p[6];

  std::string corrTitle = "";
  if(isRes) corrTitle = "RES";
  else if(isFrag) corrTitle = "FRAG";

  for(Int_t iter = 0; iter < 6; iter++){
    std::cout << Form("ak%s%d%s%s_%s%s_%s_PP_h", VsPu.c_str(), iter+1, PFCalo.c_str(), corrTitle.c_str(), MeanRes.c_str(), qg.c_str(), PtEta.c_str()) << std::endl;

    jecHist_p[iter] = (TH1F*)f->Get(Form("ak%s%d%s%s_%s%s_%s_PP_h", VsPu.c_str(), iter+1, PFCalo.c_str(), corrTitle.c_str(), MeanRes.c_str(), qg.c_str(), PtEta.c_str()));
    jecHist_p[iter]->SetMaximum(max);
    jecHist_p[iter]->SetMinimum(min);
    plotInitHist(jecHist_p[iter]);
  }

  TCanvas* plotCanv_p = new TCanvas(Form("ak%s%s%s_%s%s_%s_PP_c", VsPu.c_str(), PFCalo.c_str(), corrTitle.c_str(), MeanRes.c_str(), qg.c_str(), PtEta.c_str()), Form("ak%s%s%s_%s%s_%s_PP_c", VsPu.c_str(), PFCalo.c_str(), corrTitle.c_str(), MeanRes.c_str(), qg.c_str(), PtEta.c_str()), 3*300, 2*350);
  plotCanv_p->Divide(3, 2, 0.0, 0.0);

  TLatex* label_p = new TLatex();
  label_p->SetNDC();
  label_p->SetTextFont(43);
  label_p->SetTextSizePixels(23);

  for(Int_t iter = 0; iter < 6; iter++){
        plotCanv_p->cd(iter+1);

    if(!strcmp(VsPu.c_str(), "Vs") && iter == 5) continue;

    drawGenHist(jecHist_p[iter], MeanRes, PtEta);
    label_p->DrawLatex(.5, .9, Form("ak%s%d%s", VsPu.c_str(), iter+1, PFCalo.c_str()));

    if(iter == 0){
      if(!strcmp(PtEta.c_str(), "Pt")) label_p->DrawLatex(.5, .80, Form("p_{T}^{gen} > 20 GeV/c"));
      else label_p->DrawLatex(.5, .80, Form("p_{T}^{gen} > 50 GeV/c"));
    }

    if(iter == 1) label_p->DrawLatex(.5, .80, Form("|#eta| < 2.0"));
    if(iter == 2) label_p->DrawLatex(.5, .80, "PYTHIA");
  }

  plotCanv_p->Write("", TObject::kOverwrite);
  claverCanvasSaving(plotCanv_p, Form("pdfDir/ak%s%s%s_%s%s_%s_PP", VsPu.c_str(), PFCalo.c_str(), corrTitle.c_str(), MeanRes.c_str(), qg.c_str(), PtEta.c_str()), "pdf");

  delete label_p;
  delete plotCanv_p;
  f->Close();
  delete f;
}

void plotDiJetIniHist_PYTHHYD(const std::string histFileName, const std::string alg, const std::string qg, const std::string MeanRes, const std::string PtEta, Bool_t isFrag = false, Bool_t isRes = false)
{
  TH1::SetDefaultSumw2();
  TFile *f = new TFile(histFileName.c_str(), "UPDATE");

  Float_t max = getHistMax(MeanRes, PtEta);
  Float_t min = getHistMin(MeanRes, PtEta);

  TH1F* jecHist_p[8];
  TH1F* jecPointHist_p[8][20];

  std::cout << "A" << std::endl;

  std::string resFragStr = "";
  if(isRes) resFragStr = "RES";
  else if(isFrag) resFragStr = "FRAG";

  for(Int_t iter = 0; iter < 8; iter++){
    std::cout << Form("%s%s_%s%s_%s_%d%d_h", alg.c_str(), resFragStr.c_str(), MeanRes.c_str(), qg.c_str(), PtEta.c_str(), centCutArray[iter]/2, centCutArray[iter+1]/2) << std::endl;

    jecHist_p[iter] = (TH1F*)f->Get(Form("%s%s_%s%s_%s_%d%d_h", alg.c_str(), resFragStr.c_str(), MeanRes.c_str(), qg.c_str(), PtEta.c_str(), centCutArray[iter]/2, centCutArray[iter+1]/2));
    jecHist_p[iter]->SetMaximum(max);
    jecHist_p[iter]->SetMinimum(min);

    if(!strcmp(PtEta.c_str(), "Pt")){
      for(Int_t subIter = 0; subIter < 20; subIter++){
	jecPointHist_p[iter][subIter] = (TH1F*)f->Get(Form("%s%s_tempJetHist%s_%s_%d%d_%d", alg.c_str(), resFragStr.c_str(), qg.c_str(), PtEta.c_str(), centCutArray[iter]/2, centCutArray[iter+1]/2, subIter));
      }					      
    }
  }

  std::cout << "B" << std::endl;

  TCanvas* plotCanv_p = new TCanvas(Form("%s%s_%s%s_%s_PbPb_c", alg.c_str(), resFragStr.c_str(), MeanRes.c_str(), qg.c_str(), PtEta.c_str()), Form("%s%s_%s_%s_PbPb_c", alg.c_str(), resFragStr.c_str(), MeanRes.c_str(), PtEta.c_str()), 4*300, 2*350);
  plotCanv_p->Divide(4, 2, 0.0, 0.0);

  TCanvas* plotPointCanv_p[8];
  if(!strcmp(PtEta.c_str(), "Pt")){
    for(Int_t iter = 0; iter < 8; iter++){
      plotPointCanv_p[iter] = new TCanvas(Form("%s%s_tempJetHist%s_%s_PbPb_%d%d", alg.c_str(), resFragStr.c_str(), qg.c_str(), PtEta.c_str(), centCutArray[iter]/2, centCutArray[iter+1]/2), Form("%s%s_tempJetHist%s_%s_PbPb_%d%d", alg.c_str(), resFragStr.c_str(), qg.c_str(), PtEta.c_str(), centCutArray[iter]/2, centCutArray[iter+1]/2), 5*300, 4*350);
      plotPointCanv_p[iter]->Divide(5, 4); 
    }
  }

  std::cout << "C" << std::endl;

  TLatex* label_p = new TLatex();
  label_p->SetNDC();
  label_p->SetTextFont(43);
  label_p->SetTextSizePixels(23);

  for(Int_t iter = 0; iter < 8; iter++){
    plotInitHist(jecHist_p[7-iter]);
    plotCanv_p->cd(iter+1);
    drawGenHist(jecHist_p[7-iter], MeanRes, PtEta);
    label_p->DrawLatex(.5, .90, Form("%s%s", alg.c_str(), resFragStr.c_str()));
    label_p->DrawLatex(.5, .80, Form("%d-%d%%", centCutArray[7-iter]/2, centCutArray[8-iter]/2));

    if(iter == 0){
      if(!strcmp(PtEta.c_str(), "Pt")) label_p->DrawLatex(.5, .70, Form("p_{T}^{gen} > 20 GeV/c"));
      else label_p->DrawLatex(.5, .70, Form("p_{T}^{gen} > 50 GeV/c"));
    }

    if(iter == 1) label_p->DrawLatex(.5, .70, Form("|#eta| < 2.0"));
    if(iter == 2) label_p->DrawLatex(.5, .70, "PYT+HYD");

    if(!strcmp(PtEta.c_str(), "Pt")){
      for(Int_t subIter = 0; subIter < 20; subIter++){
	plotPointCanv_p[iter]->cd(subIter+1);
	jecPointHist_p[iter][subIter]->DrawCopy();
      }
    }
  }

  std::cout << "D" << std::endl;

  plotCanv_p->Write("", TObject::kOverwrite);
  claverCanvasSaving(plotCanv_p, Form("pdfDir/%s%s_%s%s_%s_PbPb", alg.c_str(), resFragStr.c_str(), MeanRes.c_str(), qg.c_str(), PtEta.c_str()), "pdf");

  if(!strcmp(PtEta.c_str(), "Pt")){
    for(Int_t iter = 0; iter < 8; iter++){
      plotPointCanv_p[iter]->Write("", TObject::kOverwrite);
      claverCanvasSaving(plotPointCanv_p[iter], Form("pdfDir/%s%s_tempJetHist%s_%s_PbPb_%d%d", alg.c_str(), resFragStr.c_str(), qg.c_str(), PtEta.c_str(), centCutArray[iter]/2, centCutArray[iter+1]/2), "pdf");
    }
  }

  std::cout << "E" << std::endl;

  delete label_p;
  delete plotCanv_p;

  if(!strcmp(PtEta.c_str(), "Pt")){
    for(Int_t iter = 0; iter < 8; iter++){
      delete plotPointCanv_p[iter];
    }
  }

  f->Close();
  delete f;

}


void runPlotDiJetIniHist_PYTH(const std::string histFileName)
{
  for(Int_t pfCaloIter = 1; pfCaloIter < 2; pfCaloIter++){
    for(Int_t ptEtaIter = 0; ptEtaIter < 2; ptEtaIter++){
      for(Int_t meanResIter = 0; meanResIter < 3; meanResIter++){
	for(Int_t puVsIter = 0; puVsIter < 1; puVsIter++){

	  for(Int_t qgIter = 0; qgIter < 3; qgIter++){
	    if(qgIter > 0 && meanResIter == 2) continue;

	    //	    plotDiJetIniHist_PYTH(histFileName, PuVsString[puVsIter], PFCaloString[pfCaloIter], MeanResString[meanResIter], QGString[qgIter], PtEtaString[ptEtaIter], false, false);
	    //plotDiJetIniHist_PYTH(histFileName, PuVsString[puVsIter], PFCaloString[pfCaloIter], MeanResString[meanResIter], QGString[qgIter], PtEtaString[ptEtaIter], true, false);
	    plotDiJetIniHist_PYTH(histFileName, PuVsString[puVsIter], PFCaloString[pfCaloIter], MeanResString[meanResIter], QGString[qgIter], PtEtaString[ptEtaIter], false, false);
	  }

	}
      }
    }
  }
  return;
}


void runPlotDiJetIniHist_PYTHHYD(const std::string histFileName)
{
  for(Int_t ptEtaIter = 0; ptEtaIter < 3; ptEtaIter++){
    for(Int_t meanResIter = 0; meanResIter < 4; meanResIter++){
      for(Int_t puVsIter = 1; puVsIter < 2; puVsIter++){
	for(Int_t pfCaloIter = 1; pfCaloIter < 2; pfCaloIter++){
	  for(Int_t numIter = 1; numIter < 5; numIter++){

	    for(Int_t qgIter = 0; qgIter < 3; qgIter++){
	      if(meanResIter > 1 && qgIter > 0) continue;

	      std::cout << qgIter << ", " << numIter << ", " << pfCaloIter << ", " << puVsIter << ", " << meanResIter << ", " << ptEtaIter << std::endl;	      
	      plotDiJetIniHist_PYTHHYD(histFileName, Form("ak%s%d%s", PuVsString[puVsIter+1].c_str(), numIter+1, PFCaloString[pfCaloIter].c_str()), QGString[qgIter], MeanResString[meanResIter], PtEtaString[ptEtaIter]);

	      //	      plotDiJetIniHist_PYTHHYD(histFileName, Form("ak%s%d%s", PuVsString[puVsIter+1].c_str(), numIter+1, PFCaloString[pfCaloIter].c_str()), QGString[qgIter], MeanResString[meanResIter], PtEtaString[ptEtaIter], true, false);
	      
	      //	      plotDiJetIniHist_PYTHHYD(histFileName, Form("ak%s%d%s", PuVsString[puVsIter+1].c_str(), numIter+1, PFCaloString[pfCaloIter].c_str()), QGString[qgIter], MeanResString[meanResIter], PtEtaString[ptEtaIter], true, true);
	    }

	  }
	}
      }    
    }
  }
  return;
}

