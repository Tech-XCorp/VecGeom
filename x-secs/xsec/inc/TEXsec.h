// @(#)root/base:$Id: $
// Author: Federico Carminati   27/05/13

/*************************************************************************
 * Copyright (C) 1995-2000, fca                                          *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TEXsec
#define ROOT_TEXsec


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TEXSec                                                               //
//                                                                      //
// X-section for G5 per material                                        //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <RTypes.h>
#include <TNamed.h>
class TPartIndex;
class TPXsec;
class TGraph;

class TEXsec : public TNamed {

public:
   TEXsec();
   TEXsec(Int_t z, Int_t a, Float_t emin, Float_t emax, Int_t nen, Int_t np);
   ~TEXsec();
   Bool_t AddPart(Int_t kpart, Int_t pdg, Int_t nen, Int_t nxsec, Float_t emin, Float_t emax);
   Bool_t AddPartXS(Int_t kpart, const Float_t xsec[], const Short_t dict[]);
   Bool_t AddPartIon(Int_t kpart, const Float_t dedx[]);
   Bool_t AddPartMS(Int_t kpart, const Float_t angle[], const Float_t ansig[],
		    const Float_t length[], const Float_t lensig[]);
   
   Int_t Mat() const {return fMat;}
   Float_t Emin() const {return fEmin;}
   Float_t Emax() const {return fEmax;}
   Int_t NEbins() const {return fNEbins;}
   Double_t ElDelta() const {return fElDelta;}
   Float_t XSPDG(Int_t pdg, Short_t rcode, Float_t en) const;
   Float_t XS(Int_t pindex, Short_t rindex, Float_t en) const;
   Float_t DEdxPDG(Int_t pdg, Float_t en) const;
   Float_t DEdx(Int_t pindex, Float_t en) const;
   Bool_t MSPDG(Int_t pdg, Float_t en, Float_t &ang, Float_t &asig, 
	     Float_t &len, Float_t &lsig) const;
   Bool_t MS(Int_t index, Float_t en, Float_t &ang, Float_t &asig, 
		  Float_t &len, Float_t &lsig) const;
   TGraph *XSGraph(const char* part, const char *reac, 
		   Float_t emin, Float_t emax, Int_t nbin) const;
   TGraph *DEdxGraph(const char* part, 
		   Float_t emin, Float_t emax, Int_t nbin) const;
   TGraph *MSGraph(const char* part, const char *what,
		   Float_t emin, Float_t emax, Int_t nbin) const;

   Float_t LambdaPDG(Int_t pdg, Double_t en) const;
   Float_t Lambda(Int_t pindex, Double_t en) const;
   Int_t SampleReacPDG(Int_t pdg, Double_t en) const;
   Int_t SampleReac(Int_t pindex, Double_t en) const;

   const Double_t *Cuts() const {return fCuts;}

   void DumpPointers() const;
   void Draw(Option_t *option);

private:
   Int_t          fMat;     // Material code Z*10000+A*10+metastable level
   Double_t       fAtcm3;   // Atoms per cubic cm unit density
   Float_t        fEmin;    // Min en in GeV
   Float_t        fEmax;    // Max en in Gev
   Short_t        fNEbins;  // Number of log steps in energy
   Double_t       fElDelta; // Log energy step
   Int_t          fNpart;   // Number of particles
   TPXsec        *fPXsec;   // [fNpart] Cross section table per particle
   Double_t      *fCuts;    // [fNpart] Just a placeholder for the moment

   ClassDef(TEXsec,1)  // Element X-secs

};


#endif
