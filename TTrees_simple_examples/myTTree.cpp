#include "TSystem.h"
#include "TH1F.h"
#include "TClonesArray.h"
#include "TPythia8.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TClonesArray.h"
#include "MyPart.h"
#include "TFile.h"

void myTree(Int_t nev = 10000, Int_t ndeb = 1, Int_t energia = 14000)
{
  gSystem->Load("libEg");
  gSystem->Load("libEGPythia8");

  TFile *outfile = TFile::Open("output_TTree_test.root", "RECREATE");

  TClonesArray *fparticles =  new TClonesArray("MyPart");
  TTree *tree = new TTree("partTree", "Final stage particles TTree");
  tree->Branch("particles", &fparticles, 0, 0);

  TClonesArray *particles = new TClonesArray("TParticle", 1000);
  TPythia8 *pythia8 = new TPythia8();
  pythia8->ReadString("HardQCD:all = on");
  pythia8->ReadString("Random:setSeed = on");
  pythia8->ReadString("Random:seed = 42");

  pythia8->Initialize(2212, 2212, energia /* TeV */);
  Float_t eta, phi = 0, px, py, pz, pt = 0;
  Float_t mass, energy = 0;
  Int_t pdg = 0;

  for ( Int_t iev = 0; iev < nev; iev++)
  {
    pythia8->GenerateEvent();
    if (iev < ndeb) pythia8->EventListing();
    pythia8->ImportParticles(particles, "All");
    Int_t np = particles->GetEntriesFast();
    Int_t nfp = 0;

    for (Int_t ip = 0; ip < np; ip++)
    {
      TParticle *part = (TParticle*) particles->At(ip);
      Int_t ist = part->GetStatusCode();
      Int_t pdg =part->GetPdgCode();

      Float_t eta = part->Eta();
      Float_t pt  = part->Pt();
      Int_t iMom = part->GetFirstMother();

      if (ist > 0)
      {
        MyPart *mp = static_cast<MyPart*>(fparticles->New(nfp++));
        mp->fPt   = pt;
   	    mp->fEta  = eta;
   	    mp->fPhi  = part->Phi();
   	    mp->fMass = part->GetMass();
   	    mp->fPdg  = pdg;
   	    mp->fPx   = part->Px();
   	    mp->fPy   = part->Py();
   	    mp->fPz   = part->Pz();
   	    mp->fE    = part->Energy();
        mp->fFirstMother = iMom;
      }
    }
    tree->Fill();
    fparticles->Clear();
  }
    pythia8->PrintStatistics();
    tree->Write();
    outfile->Close();
}