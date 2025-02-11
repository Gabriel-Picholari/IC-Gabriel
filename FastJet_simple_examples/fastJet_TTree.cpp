#include "TSystem.h"
#include "TH1F.h"
#include "TClonesArray.h"
#include "TPythia8.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TClonesArray.h"
#include "MyJet.h"
#include "MyQuark.h"
#include "TFile.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

Int_t findLastQuark(TClonesArray* particles, Int_t index = -1);

void template_FastJet_TTree(Int_t nev = 1000000, Int_t ndeb = 1 /* Listagem */ )
{

  gSystem->Load("libEg");
  gSystem->Load("libEGPythia8");

  //---------------------------------------------------------------------------------------------------------
  // Inicializacoes e configuracoes da TTree:s
  //---------------------------------------------------------------------------------------------------------

  TFile *outfile = TFile::Open("testeTrack.root", "RECREATE");

  TClonesArray *jets_array =  new TClonesArray("MyJet");
  TClonesArray *quarks = new TClonesArray("MyQuark");

  TTree *ttree = new TTree("TTree", "Fast_Jet TTree");

  ttree->Branch("jets_array", jets_array, 0, 0);
  ttree->Branch("quarks", quarks, 0, 0);

  //---------------------------------------------------------------------------------------------------------
  // Inicializacoes e configuracoes do Pythia:
  //---------------------------------------------------------------------------------------------------------

  TPythia8 *pythia8 = new TPythia8();
  pythia8->ReadString("HardQCD:all = on");
  pythia8->ReadString("Random:setSeed = on");
  pythia8->ReadString("Random:seed = 42");


  pythia8->Initialize(2212 /* Proton */, 2212 /* Proton */, 14000 /* TeV */); /* 14000 TeV = 14000000 GeV */

  TClonesArray *particles = new TClonesArray("TParticle", 1000);

  //---------------------------------------------------------------------------------------------------------
  // Inicializacoes e configuracoes do FastJet:
  //---------------------------------------------------------------------------------------------------------

  Double_t R = 0.4;

  fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R);

  std::vector<fastjet::PseudoJet> particles_fastjet;
  std::vector<fastjet::PseudoJet> jets;

  //---------------------------------------------------------------------------------------------------------

  //Loop de eventos
  for ( Int_t iev = 0; iev < nev; iev++)
  {
    pythia8->GenerateEvent();
    if (iev < ndeb) pythia8->EventListing();
    pythia8->ImportParticles(particles, "All");
    Int_t np = particles->GetEntriesFast();
    Int_t nfp = 0;
    Int_t nfp2 = 0;

    /* particles_fastjet.clear(); ----- wrong location? */

    //Loop de particulas
    for (Int_t ip = 0; ip < np; ip++)
    {
      TParticle *part = (TParticle*) particles->At(ip);
      Int_t ist = part->GetStatusCode();
      Int_t pdg =part->GetPdgCode();

      /* const TParticlePDG* pdgInfo = TDatabasePDG::Instance()->GetParticle(pdg); ----- perhaps unnecessary */

      if (ip == 4 || ip == 5)
      {
        Int_t iLastQuark = findLastQuark(particles, ip);

        TParticle *lqpart = (TParticle*)particles->At(iLastQuark);

        MyQuark *mq = static_cast<MyQuark *>(quarks->New(nfp++));

        mq->lqpT = lqpart->Pt();
        mq->lqEta  = lqpart->Eta();
        mq->lqPhi  = lqpart->Phi();
      }

      if (ist > 0)
      {
        Double_t px = part->Px();
        Double_t py = part->Py();
        Double_t pz = part->Pz();
        Double_t e = part->Energy();

        fastjet::PseudoJet particle(px, py, pz, e);
        particles_fastjet.push_back(particle);
      }
    }

    fastjet::ClusterSequence clusterSeq(particles_fastjet, jet_def);
    jets = clusterSeq.inclusive_jets();

    for (const fastjet::PseudoJet& jet : jets)
    {
      MyJet *mj = static_cast<MyJet*>(jets_array->New(nfp2++));
      mj->jetPt   = jet.pt();
      mj->jetEta  = jet.eta();
      mj->jetPhi  = jet.phi_std();
      mj->jetMass = jet.m();
      mj->jetPx   = jet.px();
      mj->jetPy   = jet.py();
      mj->jetPz   = jet.pz();
      mj->jetE    = jet.E();
      mj->nConstituent = jet.constituents().size();
      mj->pTLeadConstituent = (mj->nConstituent > 0) ? jet.constituents()[0].pt() : 0.0;
    }

    //Com os jatos criados, podemos colher suas informacoes

    ttree->Fill();

    particles_fastjet.clear();
    jets.clear();
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