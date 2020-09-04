#ifndef ROODOUBLECB
#define ROODOUBLECB

#include <TROOT.h>

#include <RooStats/RooStatsUtils.h>
#include <RooStats/HLFactory.h>
#include <RooRealVar.h>
#include <RooArgSet.h>
#include <RooDataHist.h>
#include <RooFitResult.h>
#include <RooWorkspace.h>
#include <RooPlot.h>
#include <RooHistPdf.h>
#include <RooAddPdf.h>
#include <RooAbsPdf.h>
#include <RooRandom.h>
#include <RooBinning.h>
#include <RooMinimizer.h>

#include <RooRealProxy.h>
#include <RooAbsReal.h>


//#include "RooAbsPdf.h"
//#include "RooRealProxy.h"
//#include "RooAbsReal.h"

class RooDoubleCBFast : public RooAbsPdf {
public:
  RooDoubleCBFast();
  RooDoubleCBFast(const char *name, const char *title,
              RooAbsReal& _x,
              RooAbsReal& _mean,
              RooAbsReal& _width,
              RooAbsReal& _alpha1,
              RooAbsReal& _n1,
              RooAbsReal& _alpha2,
              RooAbsReal& _n2
           );
  RooDoubleCBFast(const RooDoubleCBFast& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooDoubleCBFast(*this,newname); }
  inline virtual ~RooDoubleCBFast() { }
  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;

protected:

  RooRealProxy x ;
  RooRealProxy mean;
  RooRealProxy width;
  RooRealProxy alpha1;
  RooRealProxy n1;
  RooRealProxy alpha2;
  RooRealProxy n2;
  
  Double_t evaluate() const ;

private:

  ClassDef(RooDoubleCBFast,1)
};
#endif
