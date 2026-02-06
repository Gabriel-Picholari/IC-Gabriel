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
#include "TClonesArray.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

void wdecayTTree2(Int_t nev = 30000, Int_t ndeb = 1 /* Listing */ )
{
  Long_t count = 0;
  gSystem->Load("libEG");
  gSystem->Load("libEGPythia8");

  //---------------------------------------------------------------------------------------------------------
  // TTree initializations and configurations
  //---------------------------------------------------------------------------------------------------------

  TClonesArray *jets_array =  new TClonesArray("MyJet");
  TClonesArray *quarks = new TClonesArray("MyQuark");

  TFile *outfile = new TFile("wdecay2_seed_1_30K_hardQCD_all_off.root", "RECREATE"); 

  TTree *ttree = new TTree("W decay TTree 2", "Fast_Jet TTree");

  ttree->Branch("jets_array", &jets_array);
  ttree->Branch("quarks", &quarks);

  //---------------------------------------------------------------------------------------------------------
  // Initialization of histograms
  //---------------------------------------------------------------------------------------------------------

  TH1F *distanciaAngular = new TH1F("h1", "Angular separation between quarks cbar(c) and s(sbar)", 100, 0, 10);
  TH1F *bosonW_rapidity_distribution = new TH1F("bosonW_rapidity", "Initial W^{+-} boson rapidity distribution", 100, 0, 100);
  TH1F *bosonW_pT_distribution = new TH1F("bosonW_pT", "Initial W^{+-} boson transverse momentum distribution", 100, 0, 100);


  //---------------------------------------------------------------------------------------------------------
  // Pythia initializations and configurations:
  //---------------------------------------------------------------------------------------------------------

  TPythia8 pythia8 = new TPythia8();
  pythia8.ReadString("HardQCD:all = off");
  pythia8.ReadString("Random:setSeed = on");
  pythia8.ReadString("Random:seed = 1");

  pythia8.ReadString("WeakSingleBoson:ffbar2W = on");
  
  pythia8.ReadString("24:onMode = off");
  pythia8.ReadString("24:onIfMatch = 3 -4");
  
  pythia8.ReadString("-24:onMode = off");
  pythia8.ReadString("-24:onIfMatch = -3 4");
  

  pythia8.Initialize(2212 /* Proton */, 2212 /* Proton */, 14000 /* TeV */); /* 14000 TeV = 14000000 GeV */

  TClonesArray *particles = new TClonesArray("TParticle", 1000);

  Int_t wPtFlag = 0;

  // Event loop
  for ( Int_t iev = 0; iev < nev; iev++)
  {

    wPtFlag = 0;    // By default, we set the flag to 0 at the beginning of each event, that is, we assume the W boson pT is <= 10 GeV/c

    pythia8.GenerateEvent();
    if (iev == 0) pythia8.EventListing();
    pythia8.ImportParticles(particles, "All");


    Int_t np = particles->GetEntriesFast();
    Int_t nfp = 0;
    Int_t nfp2 = 0;
    
    bool foundW = false;

    // Particle loop
    for (Int_t ip = 0; ip < np; ip++)
    {
      TParticle *part = (TParticle*) particles->At(ip);
      Int_t ist = part->GetStatusCode();
      Int_t partPdg = part->GetPdgCode();

      if (!foundW && abs(partPdg) == 24)
      {
        TLorentzVector vec;
        vec.SetPxPyPzE(part->Px(), part->Py(), part->Pz(), part->Energy());
        Float_t bosonW_pT = vec.Pt();

        if (bosonW_pT > 0.0)
        {

        foundW = true;

          Float_t bosonW_rapidity = vec.Rapidity();
          bosonW_rapidity_distribution->Fill(bosonW_rapidity);

          bosonW_pT_distribution->Fill(bosonW_pT);

          if (bosonW_pT > 10)
          { 
            wPtFlag = 1;
          }
        }
      }

      if (ist > 0)
      {
        
        MyJet *fp = static_cast<MyJet*>(jets_array->New(nfp++));

        fp->fPt   = part->Pt();
        fp->fEta  = part->Eta();
        fp->fPhi  = part->Phi();
        fp->fMass = part->GetMass();
        fp->fPx   = part->Px();
        fp->fPy   = part->Py();
        fp->fPz   = part->Pz();
        fp->fE    = part->Energy();
        fp->fVx   = part->Vx();
        fp->fVy   = part->Vy();
        fp->fVz   = part->Vz();
        fp->wPtFlag = wPtFlag; // Even thoug it is rewritten for each particle, it is the same for all particles within the event

        fp->finalParticlePdg = partPdg; 
        Int_t motherIndex = part->GetFirstMother();
        TParticle *motherPart = (TParticle*)particles->At(motherIndex);
        Int_t partMotherPdg = motherPart->GetPdgCode();
        fp->finalParticleMotherPdg = partMotherPdg; // Final particle mother PDG

        Int_t secondMotherIndex = motherPart->GetFirstMother();
        TParticle *secondMotherPart = (TParticle*)particles->At(secondMotherIndex);
        Int_t secondMotherPdg = secondMotherPart->GetPdgCode();
        fp->finalParticleSecondMotherPdg = secondMotherPdg; // Final particle second mother PDG
        
        Int_t thirdMotherPdg = 0;
        Int_t thirdMotherIndex = secondMotherPart->GetFirstMother();
        if (thirdMotherIndex >= 0 && thirdMotherIndex < np) 
        { 
          TParticle *thirdMotherPart = (TParticle*)particles->At(thirdMotherIndex);
          thirdMotherPdg = thirdMotherPart->GetPdgCode();
        }
        
        Int_t index = ip;

        while (index > 1)                                                      // The loop continues until the particle reaches the proton
        {
          TParticle *ipPart = (TParticle*)particles->At(index);
          Int_t ipPdg = ipPart->GetPdgCode();                                  // Takes the pdg of the current particle 
          Int_t abs_ipPdg = abs(ipPdg);

          Int_t motherIdx1st = ipPart->GetFirstMother();                       // Takes the index of the particle's mother
          TParticle *motherPart = (TParticle*)particles->At(motherIdx1st);     // Creates a pointer to the mother
          Int_t motherPdg= motherPart->GetPdgCode();                           // Takes the pdg of the mother               

          Double_t deltaR_s, deltaR_c = 0;
          Double_t daughterEta_c, daughterPhi_c = 0;
          Double_t daughterEta_s, daughterPhi_s = 0;

          if(abs(motherPdg) == 24)                                              // Checks if we have reached a W boson
          {
            Int_t fd = motherPart->GetFirstDaughter();
            Int_t ld = motherPart->GetLastDaughter();

            for (Int_t id = fd; id <= ld; id++)
            {
              TParticle *daughterPart = (TParticle*)particles->At(id);
              Int_t daughterPdg = daughterPart->GetPdgCode();
              
              if (abs(daughterPdg) == 4 || abs(daughterPdg) == 3)               // Obtaining quark raw data
              {
                MyQuark *mq = static_cast<MyQuark *>(quarks->New(nfp2++));
          
                mq->qPdg  = daughterPdg;
                mq->qpT = daughterPart->Pt();
                mq->qEta = daughterPart->Eta();
                mq->qPhi = daughterPart->Phi();
              }

              if (abs(daughterPdg) == 4)          // Geometric match for c
              {
                daughterEta_c = daughterPart->Eta();
                daughterPhi_c = daughterPart->Phi();

                deltaR_c = TMath::Sqrt( TMath::Power(part->Eta() - daughterEta_c, 2) + TMath::Power(part->Phi() - daughterPhi_c, 2) );
              }
              else if (abs(daughterPdg) == 3)         // Geometric match for s
              {
                daughterEta_s = daughterPart->Eta();
                daughterPhi_s = daughterPart->Phi();

                deltaR_s = TMath::Sqrt( TMath::Power(part->Eta() - daughterEta_s, 2) + TMath::Power(part->Phi() - daughterPhi_s, 2) );
              }
              else
              {
                index = motherIdx1st;
              }
            }

            Float_t R_quarks = TMath::Sqrt( TMath::Power(daughterEta_c - daughterEta_s, 2) + TMath::Power(daughterPhi_c - daughterPhi_s, 2) );
            distanciaAngular->Fill(R_quarks);
            ///std::cout << deltaR_c << std::endl;
            //std::cout << deltaR_s << std::endl;

            if ( deltaR_c < deltaR_s)
            {
              fp->signalType = "charm";
              //std::cout << "Saved a charm" << std::endl;
            }
            else
            {
              fp->signalType = "strange";
              //std::cout << "Saved a strange" << std::endl;
            }
            break;

          }
          else
          {
            index = motherIdx1st;
          }
        }
      }
    }

    ttree->Fill();

    particles->Clear();
    jets_array->Clear();
    quarks->Clear();

  }


  pythia8.PrintStatistics();
  ttree->Write();

  
  TCanvas *c1 = new TCanvas("c1", "Angular distance distribution", 2500, 2500);
  c1->Divide(1, 1);

  c1->cd(1);
  distanciaAngular->SetTitle("Quarks angular separation distribution");
  distanciaAngular->GetXaxis()->SetTitle("Angular distance");
  distanciaAngular->GetYaxis()->SetTitle("Frequency");
  distanciaAngular->DrawCopy();
  

  TCanvas *c2 = new TCanvas("c2", "W^{+-} rapidity distribution", 2500, 2500);
  c2->Divide(1, 1);

  c2->cd(1);
  bosonW_rapidity_distribution->SetTitle("W^{+-} boson rapidity distribution");
  bosonW_rapidity_distribution->GetXaxis()->SetTitle("Rapidity");
  bosonW_rapidity_distribution->GetYaxis()->SetTitle("Frequency");
  bosonW_rapidity_distribution->DrawCopy();

  TCanvas *c3 = new TCanvas("c3", "W^{+-} p_{T} distribution", 2500, 2500);
  c3->Divide(1, 1);

  c3->cd(1);
  bosonW_pT_distribution->SetTitle("W^{+-} boson p_{T} distribution");
  bosonW_pT_distribution->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  bosonW_pT_distribution->GetYaxis()->SetTitle("Frequency");
  bosonW_pT_distribution->DrawCopy();

  outfile->Close();
}
