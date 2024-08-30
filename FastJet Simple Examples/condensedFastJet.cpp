#include <string>
#include <iostream>
#include <TString.h>

#include "TH1F.h"
#include "TFile.h"
#include "TMath.h"
#include "MyJet.h"
#include "MyQuark.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TPythia8.h"
#include "TParticle.h"
#include "TClonesArray.h"
#include "TDatabasePDG.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

Int_t findLastQuark(TClonesArray* particles, Int_t index = -1);

void condensado(Int_t nev = 500000, Int_t ndeb = 1)
{
  gSystem->Load("libEG");
  gSystem->Load("libEGPythia8");

  /*
    Inicializacoes
  */

  Double_t jetEta, jetPhi, jetPt = 0;
  Double_t lastQuark_pT, lastQuark_Phi, lastQuark_Eta = 0;

  TPythia8 *pythia8 = new TPythia8();
  pythia8->ReadString("HardQCD:all = on");
  pythia8->ReadString("Random:setSeed = on");
  pythia8->ReadString("Random:seed = 1");


  pythia8->Initialize(2212 /* Proton */, 2212 /* Proton */, 14000 /* TeV */); /* 14000 TeV = 14000000 GeV */

  TClonesArray *particles = new TClonesArray("TParticle", 1000);

  /*
    Histograma(s)
  */

  TH2F *quarks_and_Jets_pTs = new TH2F ("h1", "Graph of quarks by jets with equal pTs", 100, 0, 0, 100, 0, 0);

  /*
    Configuracoes do jato
  */

  Double_t R = 0.4;

  fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R);

  std::vector<fastjet::PseudoJet> particles_fastjet;
  std::vector<fastjet::PseudoJet> jets;

  /*
    Loop de eventos
  */
  for (Int_t iev = 0; iev < nev; iev++)
  {
    pythia8->GenerateEvent();
    if (iev < ndeb) pythia8->EventListing();
    pythia8->ImportParticles(particles, "All");
    Int_t np = particles->GetEntriesFast();

    for (Int_t ip = 0; ip < np; ip++)
    {
      TParticle *part = (TParticle*) particles->At(ip);
      Int_t ist = part->GetStatusCode();

      if (ip == 4 || ip == 5)
      {
        Int_t ilq = findLastQuark(particles, ip);

        TParticle *lqpart = (TParticle*) particles->At(ilq);

        lastQuark_pT = lqpart->Pt();
        lastQuark_Eta = lqpart->Eta();
        lastQuark_Phi = lqpart->Phi();
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
    } // Particle loop

    fastjet::ClusterSequence clusterSeq(particles_fastjet, jet_def);
    jets = clusterSeq.inclusive_jets();

    for (const fastjet::PseudoJet& jet : jets)
    {
      jetEta = jet.eta();
      jetPhi = jet.phi_std();
      jetPt = jet.pt();

      Double_t distancia = TMath::Sqrt(TMath::Power(jetPhi - lastQuark_Phi, 2) + TMath::Power(jetEta - lastQuark_Eta, 2));

      if (distancia < 0.1)
      {
        // cout << distancia << lastQuark_pT << jetPt << endl;
        quarks_and_Jets_pTs->Fill(lastQuark_pT, jetPt);
      }
    }

    particles_fastjet.clear();
    jets.clear();
    particles->Clear();

  } // Event loop

  TCanvas *c1 = new TCanvas("c1", "Histograms and distributions", 2500, 2500);
  c1->Divide(1, 1);

  c1->cd(1);
  quarks_and_Jets_pTs->SetTitle("Quarks pTs by jets pTs");
  quarks_and_Jets_pTs->GetXaxis()->SetTitle("LastQuark pT");
  quarks_and_Jets_pTs->GetYaxis()->SetTitle("Jet pT");
  quarks_and_Jets_pTs->Draw();


} // Main

Int_t findLastQuark(TClonesArray* particles, Int_t index=-1)
{

  if(!particles)
    return -1;

  TParticle *part = (TParticle*)particles->At(index);

  if(!part)
    return 0;

  Int_t pdgi = part->GetPdgCode();
  Int_t pdgf = pdgi;
  Int_t imom = index;
  TString pdgList;
  Int_t iLastQuark = index;

  while(pdgi == pdgf)
  {
    part = (TParticle*)particles->At(iLastQuark);
    Int_t fd = part->GetFirstDaughter();
    Int_t ld = part->GetLastDaughter();
    for(int i=fd;i<=ld;i++)
    {
      TParticle *partd = (TParticle*)particles->At(i);
      Int_t thisPdg = partd->GetPdgCode();
      //pdgList.Append(Form(";%d",thisPdg));
      if(thisPdg == pdgi) iLastQuark = i;
    }
    //if(pdgList.Contains(Form("%d",pdgi)))
      pdgf = 0;
  }
  return iLastQuark;
}