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

void wdecayTTree2(Int_t nev = 9, Int_t ndeb = 1 /* Listagem */ )
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

  TPythia8 *pythia8 = new TPythia8();
  pythia8->ReadString("HardQCD:all = off");
  pythia8->ReadString("Random:setSeed = on");
  pythia8->ReadString("Random:seed = 3");

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
    if (iev == 8) pythia8->EventListing();
    pythia8->ImportParticles(particles, "All");
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

      //std::cout << "ip = " << ip << "\t PDG = " << partPdg << "\t firstDaughter = " << part->GetFirstDaughter() << std::endl;
      
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

        while (index > 1)                                                    // O loop continua até que a partícula 
        {
          
          TParticle *ipPart = (TParticle*)particles->At(index);
          Int_t ipPdg = ipPart->GetPdgCode();                                // Toma o pdg da particula atual 
          Int_t motherIdx1st = ipPart->GetFirstMother();                        // Toma o índice da mãe da particula
          Int_t motherIdx2nd = ipPart->GetSecondMother();                        // Toma o índice da segunda mãe da particula
          TParticle *motherPart = (TParticle*)particles->At(motherIdx1st);      // Cria um ponteiro na mãe
          TParticle *mother2ndPart = 0;   // Cria um ponteiro na segunda mãe
          if(motherIdx2nd>0)
            mother2ndPart = (TParticle*)particles->At(motherIdx2nd);
          Int_t motherPdg= motherPart->GetPdgCode();                          // Toma o pdg da mãe
          Int_t mother2Pdg=-99;
          if(mother2ndPart)                      
            mother2Pdg= mother2ndPart->GetPdgCode();                         // Toma o pdg da segunda mãe
          // if(iev==9){
          //   cout<<"index="<<index<<endl;
          //   std::cout << "Particle Index:" << index << ";" << "Particle PDG: " << ipPdg << ";" << "Mother Index: " <<  motherIdx1st << "; " << "Mother PDG: " << motherPdg << std::endl;
          //   std::cout << "Contador: " << count << " ; Quark encontrado: " << ipPdg << " ; Índice do quark: " << index << " ; Partícula final associada: " << partPdg << " ; Índice: " << ip << std::endl;
          //   std::cout << "==========================================================================================================================" << std::endl;
          // }
          if (abs(ipPdg) == 4 || abs(ipPdg) == 3)                          // Verifica se a particula atual é tipo sinal
          {
            if(abs(motherPdg) == 24)                                       // Verifica se a particula é de fato sinal
            {
              MyQuark *mq = static_cast<MyQuark *>(quarks->New(nfp2++));

              mq->qPdg  = ipPdg;
              mq->qpT = ipPart->Pt();
              mq->qEta = ipPart->Eta();
              mq->qPhi = ipPart->Phi();

              count++;
              if(iev==8){
              std::cout << "Particle Index:" << index << ";" << "Particle PDG: " << ipPdg << ";" << "1st Mother Index: " <<  motherIdx1st << "; " << "Mother PDG: " << motherPdg << std::endl;
              std::cout << "Particle Index:" << index << ";" << "Particle PDG: " << ipPdg << ";" << "2nd Mother Index: " <<  motherIdx2nd << "; " << "Mother PDG: " << mother2Pdg << std::endl;
              std::cout << "Contador: " << count << " ; Quark encontrado: " << ipPdg << " ; Índice do quark: " << index << " ; Partícula final associada: " << partPdg << " ; Índice: " << ip << std::endl;
              std::cout << "==========================================================================================================================" << std::endl;
              }
            
              if (abs(ipPdg) == 3) fp->signalType = "strange";
              if (abs(ipPdg) == 4) fp->signalType = "charm";
              break;
            }
            else
            {
              index = motherIdx1st;
            }
          }
          else
          {
            index = motherIdx1st;
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

/*

while (index > 1)                                                       // O loop continua até que a partícula 
{
  TParticle *ipPart = (TParticle*)particles->At(index);
  Int_t ipPdg = ipPart->GetPdgCode();                                // Toma o pdg da particula atual 
  Int_t motherIdx1st = ipPart->GetFirstMother();                        // Toma o índice da mãe da particula
  TParticle *motherPart = (TParticle*)particles->At(motherIdx1st);      // Cria um ponteiro na mãe
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
