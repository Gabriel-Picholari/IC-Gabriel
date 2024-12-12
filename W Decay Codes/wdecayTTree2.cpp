//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
/*

The primary difference between this macro and its predecessor is that, in this version, data from all quarks—c, c̅, s, and s̅—are being saved, whereas previously we were only interested in the 
quarks from the appropriate W boson decay channel.

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
#include "TClonesArray.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

void wdecayTTree2(Int_t nev = 1, Int_t ndeb = 1 /* Listagem */ )
{
  Long_t count = 0;
  gSystem->Load("libEG");
  gSystem->Load("libEGPythia8");

  //---------------------------------------------------------------------------------------------------------
  // Inicializacoes e configuracoes da TTrees
  //---------------------------------------------------------------------------------------------------------

  TClonesArray *jets_array =  new TClonesArray("MyJet");
  TClonesArray *quarks = new TClonesArray("MyQuark");

  TFile *outfile = new TFile("wdecay2.root", "RECREATE");
  TTree *ttree = new TTree("W decay TTree 2", "Fast_Jet TTree");

  ttree->Branch("jets_array", &jets_array);
  ttree->Branch("quarks", &quarks);

  //---------------------------------------------------------------------------------------------------------
  // Inicializacoes e configuracoes do Pythia:
  //---------------------------------------------------------------------------------------------------------

  TPythia8 pythia8 = new TPythia8();
  pythia8.ReadString("HardQCD:all = off");
  pythia8.ReadString("Random:setSeed = on");
  pythia8.ReadString("Random:seed = 3");

  pythia8.ReadString("WeakSingleBoson:ffbar2W = on");

  pythia8.ReadString("24:onMode = off");
  pythia8.ReadString("24:onIfMatch = 3 -4");
  
  pythia8.ReadString("-24:onMode = off");
  pythia8.ReadString("-24:onIfMatch = -3 4");

  pythia8.Initialize(2212 /* Proton */, 2212 /* Proton */, 14000 /* TeV */); /* 14000 TeV = 14000000 GeV */

  TClonesArray *particles = new TClonesArray("TParticle", 1000);

  //Loop de eventos
  for ( Int_t iev = 0; iev < nev; iev++)
  {
    pythia8.GenerateEvent();
    if (iev == 0) pythia8.EventListing();
    pythia8.ImportParticles(particles, "All");


    Int_t np = particles->GetEntriesFast();
    Int_t nfp = 0;
    Int_t nfp2 = 0;
    
    std::cout << std::endl;
    std::cout << "NEW EVENT ITERATION" << std::endl;
    std::cout << std::endl;


    //Loop de particulas
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

      Float_t eta = part->Eta();
      //if (eta < -2 || eta > 2) continue;

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

        Int_t index = ip;

        while (index > 1)                                                      // O loop continua até que a partícula chegue ao próton
        {
          TParticle *ipPart = (TParticle*)particles->At(index);
          Int_t ipPdg = ipPart->GetPdgCode();                                  // Toma o pdg da particula atual 
          Int_t motherIdx1st = ipPart->GetFirstMother();                       // Toma o índice da mãe da particula
          TParticle *motherPart = (TParticle*)particles->At(motherIdx1st);     // Cria um ponteiro na mãe
          Int_t motherPdg= motherPart->GetPdgCode();                           // Toma o pdg da mãe

          //std::cout << "--------------------------------------------------------------------------------------------------------------------------" << std::endl;
          //std::cout << "Particle Index:" << index << ";" << "Particle PDG: " << ipPdg << ";" << "Mother Index: " <<  motherIdx1st << "; " << "Mother PDG: " << motherPdg << std::endl;                   

          Double_t deltaR_s, deltaR_c = 0;
          Double_t daughterEta, daughterPhi = 0;

          if(abs(motherPdg) == 24)          // Verifica se chegamos a um bóson W
          {
            Int_t fd = motherPart->GetFirstDaughter();
            Int_t ld = motherPart->GetLastDaughter();

            for (Int_t id = fd; id <= ld; id++)
            {
              TParticle *daughterPart = (TParticle*)particles->At(id);
              Int_t daughterPdg = daughterPart->GetPdgCode();

              if (abs(daughterPdg) == 4 || abs(daughterPdg) == 3)         // Apenas para teste
              {
                MyQuark *mq = static_cast<MyQuark *>(quarks->New(nfp2++));

                mq->qPdg  = ipPdg;
                mq->qpT = ipPart->Pt();
                mq->qEta = ipPart->Eta();
                mq->qPhi = ipPart->Phi();
              }
              if (abs(daughterPdg) == 4)          // Match geométrico para o c
              {
                daughterEta = daughterPart->Eta();
                daughterPhi = daughterPart->Phi();

                deltaR_c = TMath::Sqrt( TMath::Power(part->Eta() - daughterEta, 2) + TMath::Power(part->Phi() - daughterPhi, 2) );
              }
              else if (abs(daughterPdg) == 3)         // Match geométrico para o s
              {
                daughterEta = daughterPart->Eta();
                daughterPhi = daughterPart->Phi();

                deltaR_s = TMath::Sqrt( TMath::Power(part->Eta() - daughterEta, 2) + TMath::Power(part->Phi() - daughterPhi, 2) );
              }
              else
              {
                index = motherIdx1st;
              }
            }
            std::cout << deltaR_c << std::endl;
            std::cout << deltaR_s << std::endl;

            if ( deltaR_c > deltaR_s)
            {
              fp->signalType = "charm";
              std::cout << "Salvamos um charm" << std::endl;
            }
            else
            {
              fp->signalType = "strange";
              std::cout << "Salvamos um strange" << std::endl;
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