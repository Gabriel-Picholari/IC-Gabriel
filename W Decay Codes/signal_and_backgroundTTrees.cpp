//---------------------------------------------------------------------------------------------------------//
//                                                                                                         //
// Macro with the objective of reading the TTree created in the "wdecayTTree.cpp" macro. It performs       //
// jet reconstruction and then matches it with the information of the last quarks to correlate the jet to  //
// the quark that hadronized into the jet. Data corresponding to the matched quark is then organized into  //
// a signal tree. Non-matching quarks have their data organized in background trees. Note that these trees //
// do not distinguish between testing and training data; this responsibility remains under the control of  //
// another macro.                                                                                          //
//                                                                                                         //
//---------------------------------------------------------------------------------------------------------//

#include "TFile.h"
#include "MyJet.h"
#include "TMath.h"
#include "MyQuark.h"
#include "TSystem.h"
#include "TPythia8.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TClonesArray.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

void signal_and_backgroundTTrees(const char *fileName)
{

  gSystem->Load("libEG");
  gSystem->Load("libEGPythia8");

  //---------------------------------------------------------------------------------------------------------
  // Inicializacao das variaveis
  //---------------------------------------------------------------------------------------------------------

  Float_t fpPt, fpEta, fpPhi, fpE, fpPx, fpPy, fpPz, fpMass = 0;
  Float_t jetPt, jetEta, jetPhi, jetE, jetPx, jetPy, jetPz, jetMass, jetNConst, pT_LeadConst = 0;
  Float_t c_pT, c_Eta, c_Phi = 0;
  Float_t sbar_pT, sbar_Eta, sbar_Phi = 0;
  Float_t cbar_pT, cbar_Eta, cbar_Phi = 0;
  Float_t match_R = 0.1;

  //---------------------------------------------------------------------------------------------------------
  // Inicializacoes e configuracoes do FastJet:
  //---------------------------------------------------------------------------------------------------------

  Float_t jetR = 0.5;

  fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, jetR);

  std::vector<fastjet::PseudoJet> particles_fastjet;
  std::vector<fastjet::PseudoJet> jets;

  //---------------------------------------------------------------------------------------------------------
  // Inicializacao do arquivo.root e das TTrees
  //---------------------------------------------------------------------------------------------------------

  TFile *file = new TFile(fileName, "READ");
  TTree *ttree = dynamic_cast<TTree *>(file->Get("W decay TTree"));

  TClonesArray *jets_array = new TClonesArray("MyJet");
  TClonesArray *quarks = new TClonesArray("MyQuark");

  ttree->SetBranchAddress("jets_array", &jets_array);
  ttree->SetBranchAddress("quarks", &quarks);

  Long64_t ne = ttree->GetEntries();

  //---------------------------------------------------------------------------------------------------------
  // Inicializacao das TTrees TMVA
  //---------------------------------------------------------------------------------------------------------

  Float_t pT_c, pT_lConst_c, label_c, nConst_c, averageAng_c, sigmaKT_c = 0;

  TFile *filteredDataFile = new TFile("filteredOutput_3var.root", "RECREATE");

  TTree *signalTree_c = new TTree("SignalTree_c", "Tree with signal data from c quark");
  
  signalTree_c->Branch("pT_c", &pT_c);
  signalTree_c->Branch("label_c", &label_c);
  signalTree_c->Branch("nConst_c", &nConst_c);
  signalTree_c->Branch("sigmaKT_c", &sigmaKT_c);
  signalTree_c->Branch("pT_lConst_c", &pT_lConst_c);
  signalTree_c->Branch("averageAng_c", &averageAng_c);


  TTree *backgroundTree_c = new TTree("BackgroundTree_c", "Tree with background data from c quark");
  backgroundTree_c->Branch("pT_c", &pT_c);
  backgroundTree_c->Branch("label_c", &label_c);
  backgroundTree_c->Branch("nConst_c", &nConst_c);
  backgroundTree_c->Branch("sigmaKT_c", &sigmaKT_c);
  backgroundTree_c->Branch("pT_lConst_c", &pT_lConst_c);
  backgroundTree_c->Branch("averageAng_c", &averageAng_c);


  //---------------------------------------------------------------------------------------------------------
  // Loop equivalente ao loop de eventos
  //---------------------------------------------------------------------------------------------------------

  for (Long64_t ni = 0; ni < ne; ni++)
  {
    ttree->GetEntry(ni);

    //---------------------------------------------------------------------------------------------------------
    // Loop equivalente ao loop de particulas
    //---------------------------------------------------------------------------------------------------------

    for (Int_t nj = 0; nj < jets_array->GetEntries(); nj++)
    {
      MyJet *fp = static_cast<MyJet *>(jets_array->At(nj));

      fpPx = fp->fPx;
      fpPy = fp->fPy;
      fpPz = fp->fPz;
      fpE = fp->fE;

      fastjet::PseudoJet particle(fpPx, fpPy, fpPz, fpE);
      particles_fastjet.push_back(particle);
    }

    fastjet::ClusterSequence clusterSeq(particles_fastjet, jet_def);
    jets = clusterSeq.inclusive_jets();

    for (const fastjet::PseudoJet &jet : jets)
    {
      jetPt = jet.pt();
      jetEta = jet.eta();
      jetPhi = jet.phi_std();
      jetMass = jet.m(); // Invariant mass
      jetPx = jet.px();
      jetPy = jet.py();
      jetPz = jet.pz();
      jetE = jet.E();
      jetNConst = jet.constituents().size();
      pT_LeadConst = 0.0;
      for (const fastjet::PseudoJet &constituent : jet.constituents())
      {
        if (constituent.pt() > pT_LeadConst)
        {
          pT_LeadConst = constituent.pt();
        }
      }

      Float_t averAng, sigmaKT = 0;

      for (Int_t i = 0; i < jetNConst; ++i) 
      {
        Float_t pt_constituentes = jet.constituents()[i].pt();
        Float_t eta_constituentes = jet.constituents()[i].eta();
        Float_t phi_constituentes = jet.constituents()[i].phi();

        Float_t partAng = TMath::Sqrt(TMath::Power(TMath::Abs(phi_constituentes) - TMath::Abs(jetPhi), 2) + TMath::Power(eta_constituentes - jetEta, 2));

        averAng = averAng + partAng;
        sigmaKT = sigmaKT + (TMath::Power(pt_constituentes - (jetPt / jetNConst), 2));
      }

      averAng = averAng / jetNConst;
      sigmaKT = sigmaKT / jetNConst;

      //---------------------------------------------------------------------------------------------------------
      // Comparacao dos jatos com os quarks -> Match entre ambos
      //---------------------------------------------------------------------------------------------------------

      for (Int_t nj2 = 0; nj2 < quarks->GetEntries(); nj2++)
      {

        MyQuark *quarkPdg = static_cast<MyQuark *>(quarks->At(nj2));
        Int_t nj2Pdg = quarkPdg->qPdg;

        if (nj2Pdg == 4)
        {
          MyQuark *cq = static_cast<MyQuark *>(quarks->At(nj2));

          c_pT = cq->qpT;
          c_Eta = cq->qEta;
          c_Phi = cq->qPhi;

          Float_t distancia_c = TMath::Sqrt(TMath::Power(jetPhi - c_Phi, 2) + TMath::Power(jetEta - c_Eta, 2));

          //if (jetPt < 10) continue; SHALL BE DONE USING TCUT THROUGH TMVA
          
          if (distancia_c <= match_R) // Matched ---> SIGNAL
          {
            label_c = 1;
            pT_c = jetPt;
            sigmaKT_c = sigmaKT;
            nConst_c = jetNConst;
            averageAng_c = averAng;
            pT_lConst_c = pT_LeadConst;

            signalTree_c->Fill();
          }

          else// Not matched ---> BACKGROUND
          {
            label_c = 0;
            pT_c = jetPt;
            sigmaKT_c = sigmaKT;
            nConst_c = jetNConst;
            averageAng_c = averAng;
            pT_lConst_c = pT_LeadConst;
            
            backgroundTree_c->Fill();
          }
        }

      } // End of particles and jets particular matches

    } // End of individual jet creation

    particles_fastjet.clear();
    jets.clear();
    jets_array->Clear();
    quarks->Clear();

  } // End of event loop equivalent

  signalTree_c->Write();
  backgroundTree_c->Write();

  file->Close();
  filteredDataFile->Close();

  delete jets_array;
  delete quarks;

}
