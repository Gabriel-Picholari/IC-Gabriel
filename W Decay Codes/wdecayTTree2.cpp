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


Int_t findLastQuark(TClonesArray* particles, Int_t index = -1);

void wdecayTTree2(Int_t nev = 1000, Int_t ndeb = 1 /* Listagem */ )
{

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
      Float_t eta = part->Eta();

      if (eta < -2 || eta > 2) continue;

      if (ist == 1)
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

        Int_t idx = ip;
        TParticle *ipPart = part;

        while (idx >= 0) 
        {

          Int_t ipPdg = ipPart->GetPdgCode();
          Int_t motherIdx = ipPart->GetMother(idx);
          

          if (motherIdx >= 0) 
          {
            
            TParticle *motherPart = (TParticle*) particles->At(motherIdx);
            motherPart->GetPdgCode();

            Int_t motherPdg = motherPart->GetPdgCode();

            if (abs(ipPdg) == 4 || abs(ipPdg) == 3) 
            {
              if (abs(motherPdg) == 24) 
              {
                MyQuark *mq = static_cast<MyQuark *>(quarks->New(nfp2++)); // Partículas de sinal

                mq->qPdg  = ipPdg;
                mq->qpT = ipPart->Pt();
                mq->qEta = ipPart->Eta();
                mq->qPhi = ipPart->Phi();

                if (abs(ipPdg) == 3) fp->signalType = "strange";
                if (abs(ipPdg) == 4) fp->signalType = "charm";

                break;
              }
            }

            idx = motherIdx;
            ipPart = motherPart;  // Atualiza o ponteiro da partícula para a mãe, repetimos o loop com a primeira mãe
          }
          else
          {
            fp->signalType= "background";
            break;  // Não há mais mãe, saímos do loop, i.e., o índice da mãe é menor que 0
          }

        }
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
