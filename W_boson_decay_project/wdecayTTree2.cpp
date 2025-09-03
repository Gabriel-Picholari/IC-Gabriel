//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

/*

The primary difference between this macro and its predecessor is that, in this version, data from all quarks—c, c̅, s, and s̅—are being saved, whereas previously we were only interested in the 
quarks from the appropriate W boson decay channel.

 */

//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

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

void wdecayTTree2(Int_t nev = 10000, Int_t ndeb = 1 /* Listing */ )
{
  Long_t count = 0;
  gSystem->Load("libEG");
  gSystem->Load("libEGPythia8");

  //---------------------------------------------------------------------------------------------------------
  // TTree initializations and configurations
  //---------------------------------------------------------------------------------------------------------

  TClonesArray *jets_array =  new TClonesArray("MyJet");
  TClonesArray *quarks = new TClonesArray("MyQuark");

  TFile *outfile = new TFile("wdecay2_seed_1_10K_hardQCD_all_off.root", "RECREATE");
  TTree *ttree = new TTree("W decay TTree 2", "Fast_Jet TTree");

  ttree->Branch("jets_array", &jets_array);
  ttree->Branch("quarks", &quarks);

  //---------------------------------------------------------------------------------------------------------
  // Initialization of histograms
  //---------------------------------------------------------------------------------------------------------

  TH1F *distanciaAngular = new TH1F("h1", "Distância angular entre os quarks cbar(c) e s(sbar)", 100, 0, 10);
  TH1F *bosonW_rapidity_distribution = new TH1F("bosonW_rapidity", "Boson W^{+-} rapidity distribution", 100, 0, 100);

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

  // Event loop
  for ( Int_t iev = 0; iev < nev; iev++)
  {
    pythia8.GenerateEvent();
    if (iev == 0) pythia8.EventListing();
    pythia8.ImportParticles(particles, "All");


    Int_t np = particles->GetEntriesFast();
    Int_t nfp = 0;
    Int_t nfp2 = 0;
    
    //std::cout << std::endl;
    //std::cout << "NEW EVENT ITERATION" << std::endl;
    //std::cout << std::endl;


    // Particle loop
    for (Int_t ip = 0; ip < np; ip++)
    {
      TParticle *part = (TParticle*) particles->At(ip);
      Int_t ist = part->GetStatusCode();
      Int_t partPdg = part->GetPdgCode();

      /*
      std::cout << "ip = " << ip << "\t PDG = " << partPdg << "\t firstDaughter = " << part->GetFirstDaughter() << std::endl;
      Int_t fd = 0, ld = 0;
      if (abs(partPdg) == 24) 
      {
        fd = part->GetFirstDaughter();
        ld = part->GetLastDaughter();
        TParticle *partFD = (TParticle*) particles->At(fd);
        TParticle *partLD = (TParticle*) particles->At(ld);
        std::cout << "Stack number: " << ip << " ; First daughter: " << fd << " ; Last daughter: " << ld << std::endl;
        std::cout << "PDG first daughter: " << partFD->GetPdgCode() << "; PDG last daughter: " << partLD->GetPdgCode() << std::endl;
        std::cout << "------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
      }
      */

      if (abs(partPdg) == 24)
      {
        TLorentzVector vec;
        vec.SetPxPyPzE(part->Px(), part->Py(), part->Pz(), part->Energy());

        Float_t bosonW_rapidity = vec.Rapidity();
        bosonW_rapidity_distribution->Fill(bosonW_rapidity);

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
          
          //if (thirdMotherPart) 
          //{
            thirdMotherPdg = thirdMotherPart->GetPdgCode();
          //}
        }

        const std::unordered_set<int> charmPdgSet = {411, 421, 413, 423, 415, 425, 431, 433, 435};
        const std::unordered_set<int> strangePdgSet = {130, 310, 311, 321, 313, 323, 315, 325, 317, 327, 319, 329};
        
        Int_t index = ip;

        Bool_t hasCharmedHadron, hasStrangeHadron = false;

        while (index > 1)                                                      // The loop continues until the particle reaches the proton
        {
          TParticle *ipPart = (TParticle*)particles->At(index);
          Int_t ipPdg = ipPart->GetPdgCode();                                  // Takes the pdg of the current particle 
          Int_t abs_ipPdg = abs(ipPdg);
          
          if (charmPdgSet.count(abs_ipPdg)) hasCharmedHadron = true;           // If true, then the decay chain passed through a hadron of interest (check sets above)
          if (strangePdgSet.count(abs_ipPdg)) hasStrangeHadron = true;

          Int_t motherIdx1st = ipPart->GetFirstMother();                       // Takes the index of the particle's mother
          TParticle *motherPart = (TParticle*)particles->At(motherIdx1st);     // Creates a pointer to the mother
          Int_t motherPdg= motherPart->GetPdgCode();                           // Takes the pdg of the mother

          //std::cout << "--------------------------------------------------------------------------------------------------------------------------" << std::endl;
          //std::cout << "Particle Index:" << index << ";" << "Particle PDG: " << ipPdg << ";" << "Mother Index: " <<  motherIdx1st << "; " << "Mother PDG: " << motherPdg << std::endl;                   

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
  distanciaAngular->SetTitle("Quarks angular distance distribution");
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

  outfile->Close();
}

/*

while (index > 1)                                                    // O loop continua até que a partícula 
{
  TParticle *ipPart = (TParticle*)particles->At(index);
  Int_t ipPdg = ipPart->GetPdgCode();                                // Toma o pdg da particula atual 
  Int_t motherIdx1st = ipPart->GetFirstMother();                     // Toma o índice da mãe da particula
  TParticle *motherPart = (TParticle*)particles->At(motherIdx1st);   // Cria um ponteiro na mãe
  Int_t motherPdg= motherPart->GetPdgCode();                         // Toma o pdg da mãe

  std::cout << "--------------------------------------------------------------------------------------------------------------------------" << std::endl;
  std::cout << "Particle Index:" << index << ";" << "Particle PDG: " << ipPdg << ";" << "Mother Index: " <<  motherIdx1st << "; " << "Mother PDG: " << motherPdg << std::endl;

  if (abs(ipPdg) != 4 && abs(ipPdg) != 3)                          // Verifica se a particula atual é tipo sinal
  {
    index = motherIdx1st;
    continue;
  }                         

  if(abs(motherPdg) == 24)                                         // Verifica se a particula é de fato sinal
  {
    MyQuark *mq = static_cast<MyQuark *>(quarks->New(nfp2++));

    mq->qPdg  = ipPdg;
    mq->qpT = ipPart->Pt();
    mq->qEta = ipPart->Eta();
    mq->qPhi = ipPart->Phi();

    count++;
  
    if (abs(ipPdg) == 3) fp->signalType = "strange";
    if (abs(ipPdg) == 4) fp->signalType = "charm";
  }
  else
  {
    index = motherIdx1st;
  }
}

*/

/*
while (index > 1) 
{

  TParticle *ipPart = (TParticle*)particles->At(index);
  Int_t ipPdg = ipPart->GetPdgCode();

  std::vector<int> motherList = pythia8.event[index]->motherList();      // Vetor das partículas mães
  Int_t bestMotherIdx = -1;                                             // Inicialização da mãe sinal            
  Double_t minDeltaR = std::numeric_limits<double>::max();              // Inicialização do Delta

  for (Int_t motherIdx : motherList) 
  {
    TParticle *motherPart = (TParticle*)particles->At(motherIdx);
    Int_t motherPdg = motherPart->GetPdgCode();

    Int_t grandmotherIdx = motherPart->GetFirstMother();                        //Verifica se a mãe do quark de interesse é um bóson W (o que configura o quark como sinal)
    TParticle *grandmotherPart = (TParticle*)particles->At(grandmotherIdx);
    Int_t grandmotherPdg = grandmotherPart->GetPdgCode();

    if ( (abs(motherPdg) == 4 || abs(motherPdg) == 3) && abs(grandmotherPdg) == 24 )
    {
      Double_t deltaR = sqrt( pow(ipPart->Eta() - motherPart->Eta(), 2) + pow(ipPart->Phi() - motherPart->Phi(), 2) );        // Matching geométrico

      if (deltaR < minDeltaR) 
      {
        minDeltaR = deltaR;
        bestMotherIdx = motherIdx;
      }
    }
  }

  if (bestMotherIdx != -1)          // Se a partícula passou pela verificação e é sinal, então fazemos o tagging
  {
    TParticle *bestMother = (TParticle*)particles->At(bestMotherIdx);
    Int_t bestMotherPdg = bestMother->GetPdgCode();

    MyQuark *mq = static_cast<MyQuark *>(quarks->New(nfp2++));
    mq->qPdg  = bestMotherPdg;
    mq->qpT = bestMother->Pt();
    mq->qEta = bestMother->Eta();
    mq->qPhi = bestMother->Phi();

    if (abs(bestMotherPdg) == 3) fp->signalType = "strange";
    if (abs(bestMotherPdg) == 4) fp->signalType = "charm";
    break;

  } 
  else          // Se o índice da mãe não foi modificado (permanece -1, como foi inicializado) então a partícula não é sinal e devemos subir
  {
    Int_t motherIdx1st = ipPart->GetFirstMother();
    index = motherIdx1st;
  }
}
*/