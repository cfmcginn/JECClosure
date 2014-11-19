#include <vector>
#include <iostream>
#include "TMath.h"

Float_t getMean(std::vector<Float_t>* inVect_p)
{
  if(inVect_p->size() == 0){
    std::cout << "Passed empty vector. Return 0." << std::endl;
    return 0;
  }

  Float_t sum = 0;
  for(Int_t jtIter = 0; jtIter < (Int_t)inVect_p->size(); jtIter++){
    sum += inVect_p->at(jtIter);
  }

  return sum/((Float_t)inVect_p->size());
}


Float_t getError(std::vector<Float_t>* inVect_p, Float_t mean)
{
  if(inVect_p->size() < 2){
    std::cout << "Passed vector of size 1 or less. Return 0." << std::endl;
    return 0;
  }

  Float_t error = 0;
  for(Int_t jtIter = 0; jtIter < (Int_t)inVect_p->size(); jtIter++){
    error += TMath::Power((inVect_p->at(jtIter) - mean), 2);
  }

  error = TMath::Sqrt(error/((Float_t)inVect_p->size() - 1.0));
  return error/TMath::Sqrt((Float_t)inVect_p->size());
}



Float_t getMeanWeighted(std::vector<Float_t>* inVect_p, std::vector<Float_t>* weightVect_p)
{
  if(inVect_p->size() == 0){
    std::cout << "Passed empty vector. Return 0." << std::endl;
    return 0;
  }

  Float_t sum = 0;
  Float_t denom = 0;
  for(Int_t jtIter = 0; jtIter < (Int_t)inVect_p->size(); jtIter++){
    sum += inVect_p->at(jtIter)*weightVect_p->at(jtIter);
    denom += weightVect_p->at(jtIter);
  }

  return sum/denom;
}


Float_t getErrorWeighted(std::vector<Float_t>* inVect_p, Float_t mean, std::vector<Float_t>* weightVect_p)
{
  if(inVect_p->size() < 2){
    std::cout << "Passed vector of size 1 or less. Return 0." << std::endl;
    return 0;
  }

  Float_t error = 0;
  Float_t denom = 0;
  for(Int_t jtIter = 0; jtIter < (Int_t)inVect_p->size(); jtIter++){
    error += TMath::Power((inVect_p->at(jtIter) - mean), 2)*weightVect_p->at(jtIter);
    denom += weightVect_p->at(jtIter);
  }

  error = TMath::Sqrt(error/(denom - 1.0));
  return error/TMath::Sqrt(denom);
}

