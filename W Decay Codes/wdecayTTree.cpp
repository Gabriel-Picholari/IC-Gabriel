//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
/* 

  Macro with the main objective of creating a TTree containing specific data from a pp collision that necessarily results in the formation of a W boson. Note that the collected data is filtered to 
include only the final particles (necessary for jet reconstruction) and the last quark (refer to the function "findLastQuark" for further details) of the showering process resulting from the 
specified decay of the W boson. The specific decay mode is either s-cbar or c-sbar.                                            

*/                                                                                                  
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#include "TH1F.h"
#include "MyJet.h"
#include "TMath.h"
#include "TFile.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "MyQuark.h"
#include "TPythia8.h"
#include "TParticle.h"
#include "TClonesArray.h"
#include "TDatabasePDG.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

Int_t findLastQuark(TClonesArray* particles, Int_t index = -1);

void wdecayTTree(Int_t nev = 1000, Int_t ndeb = 1 /* Listagem */ )
{

  gSystem->Load("libEG");
  gSystem->Load("libEGPythia8");

  //---------------------------------------------------------------------------------------------------------
  // Inicializacoes e configuracoes da TTrees
  //---------------------------------------------------------------------------------------------------------

  TClonesArray *jets_array =  new TClonesArray("MyJet");
  TClonesArray *quarks = new TClonesArray("MyQuark");

  TFile *outfile = new TFile("wdecay1K_QCD_OFF.root", "RECREATE");
  TTree *ttree = new TTree("W decay TTree", "Fast_Jet TTree");

  ttree->Branch("jets_array", &jets_array);
  ttree->Branch("quarks", &quarks);

  //---------------------------------------------------------------------------------------------------------
  // Inicializacoes e configuracoes do Pythia:
  //---------------------------------------------------------------------------------------------------------

  TPythia8 *pythia8 = new TPythia8();
  pythia8->ReadString("HardQCD:all = off");
  pythia8->ReadString("Random:setSeed = on");
  pythia8->ReadString("Random:seed = 1");

  pythia8->ReadString("WeakSingleBoson:ffbar2W = on");
  
  pythia8->ReadString("24:onMode = off");
  pythia8->ReadString("24:onIfMatch = 3 -4");
  
  pythia8->ReadString("-24:onMode = off");
  pythia8->ReadString("-24:onIfMatch = -3 4");

  pythia8->Initialize(2212 /* Proton */, 2212 /* Proton */, 14000 /* TeV */); /* 14000 TeV = 14000000 GeV */

  TClonesArray *particles = new TClonesArray("TParticle", 1000);

  //Loop de eventos
  for ( Int_t iev = 0; iev < nev; iev++)
  {
    pythia8->GenerateEvent();
    if (iev < ndeb) pythia8->EventListing();
    pythia8->ImportParticles(particles, "All");
    Int_t np = particles->GetEntriesFast();
    Int_t nfp = 0;
    Int_t nfp2 = 0;


    //Loop de particulas
    for (Int_t ip = 0; ip < np; ip++)
    {
      TParticle *part = (TParticle*) particles->At(ip);
      Int_t ist = part->GetStatusCode();
      Int_t pdg = part->GetPdgCode();

      if (TMath::Abs(pdg) == 24)
      {

        Int_t fd = part->GetFirstDaughter();
        Int_t ld = part->GetLastDaughter();

        for ( Int_t id = fd; id <= ld; id++ )
        {
          Int_t iLastQuark = findLastQuark(particles, id);
          TParticle *idPart = (TParticle*) particles->At(iLastQuark);
          Int_t idPdg = idPart->GetPdgCode();
          
          MyQuark *mq = static_cast<MyQuark *>(quarks->New(nfp++));
          
          mq->qPdg  = idPdg;
          mq->qpT = idPart->Pt();
          mq->qEta = idPart->Eta();
          mq->qPhi = idPart->Phi();
        }
      }

      if (ist == 1)
      {
        TParticle *fpart = (TParticle*) particles->At(ip);
        MyJet *fp = static_cast<MyJet*>(jets_array->New(nfp2++));

        fp->fPt   = fpart->Pt();
        fp->fEta  = fpart->Eta();
        fp->fPhi  = fpart->Phi();
        fp->fMass = fpart->GetMass();
        fp->fPx   = fpart->Px();
        fp->fPy   = fpart->Py();
        fp->fPz   = fpart->Pz();
        fp->fE    = fpart->Energy();
      }
    }
    
    //cout << count_c << " e " << count_sbar << endl;
    //Com os jatos criados, podemos colher suas informacoes

    ttree->Fill();

    particles->Clear();
    jets_array->Clear();
    quarks->Clear();
  }

  pythia8->PrintStatistics();
  ttree->Write();
  outfile->Close();
}

Int_t findLastQuark(TClonesArray* particles, Int_t index=-1)
{
  if(!particles)
    return -1;
    //cout << "index=" << index << endl;

  TParticle *part = (TParticle*)particles->At(index);

  if(!part)
    return 0;

  Int_t pdgi = part->GetPdgCode();
  Int_t pdgf = pdgi;
  Int_t imom = index;
  TString PdgList;
  Int_t iLastQuark = index;
  //cout << "pdgi= " << pdgi << endl;
  bool samepdg = false;

  while(pdgi == pdgf)
  {
    part = (TParticle*)particles->At(iLastQuark);
    Int_t fd = part->GetFirstDaughter();
    Int_t ld = part->GetLastDaughter();
    //cout << "fd,ld =" << fd << "," << ld << endl;

    for(int i = fd; i <= ld; i++)
    {
      TParticle *partd = (TParticle*)particles->At(i);
      Int_t thisPdg = partd->GetPdgCode();
      //PdgList.Append(Form(";%d",thisPdg));
      //cout << "thispdg = " << thisPdg << endl;

      if(thisPdg == pdgi)
      {
	       iLastQuark = i;
	       samepdg= true;
	    }
    }
    if(!samepdg)
      pdgf = 0;
      //cout << pdgf << "," << iLastQuark << ";";
    samepdg = false;
  }
  //cout << endl;
  return iLastQuark;
}