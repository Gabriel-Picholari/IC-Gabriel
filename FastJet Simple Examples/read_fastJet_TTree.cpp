#include <iostream>
#include <string>
#include "TSystem.h"
#include "TH1F.h"
#include "TPythia8.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TClonesArray.h"
#include "MyJet.h"
#include "MyQuark.h"
#include "TFile.h"
#include <TString.h>
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

void template_read_FastJet_TTree(const char* fileName)
{

  gSystem->Load("libEG");
  gSystem->Load("libEGPythia8");

  //---------------------------------------------------------------------------------------------------------
  //Inicializacao das variaveis
  //---------------------------------------------------------------------------------------------------------

  Double_t jetMass, jetEta, jetPhi = 0;
  Double_t jetPx, jetPy, jetPz, jetPt, quarkPt, jetE = 0;
  Int_t nJet, nConstituent = 0;
  Double_t pTLeadConstituent = 0;
  Double_t lastQuark_pT, lastQuark_Phi, lastQuark_Eta = 0;

  //---------------------------------------------------------------------------------------------------------
  // Inicializacao dos histogramas
  //---------------------------------------------------------------------------------------------------------

  TH1F *nJets_event = new TH1F("h1", "Number of jets per event", 100, 0, 0);
  TH1F *pT_Jet = new TH1F("h2", "Jet's transverse momentum", 100, 0, 0);
  TH1F *rapidity_Jet = new TH1F("h3", "Jet's rapidity", 100, 0, 0);
  TH1F *phi_Jet = new TH1F("h4", "Jet's phi", 100, 0, 0);
  TH1F *mass_Jet = new TH1F("h5", "Jet's mass", 100, 0, 0);
  TH1F *nConstituent_Jet = new TH1F("h6", "Constituents number per Jet", 100, 0, 0);
  TH2F *pTLeadConstituent_pTJet = new TH2F ("h7", "Graph of lead constituent's pT by jet's pT", 100, 0, 0, 100, 0, 0);
  TH2F *nConstituent_Mass = new TH2F ("h8", "Graph of lead constituent's pT by jet's Mass", 20, 0, 20, 100, 0, 0);
  TH2F *quarks_and_Jets_pTs = new TH2F ("h9", "Graph of quarks by jets with equal pTs", 100, 0, 20, 100, 0, 20);

  //---------------------------------------------------------------------------------------------------------
  // Inicializacao do arquivo.root e das TTrees com prevencao de erros
  //---------------------------------------------------------------------------------------------------------

  TFile *file = TFile::Open(fileName, "READ");

  if (!file || !file->IsOpen())
  {
      std::cerr << "Error: Could not open the file." << std::endl;
      return;
  }

  TTree *ttree = dynamic_cast<TTree *>(file->Get("TTree"));


  if (!ttree)
  {
      std::cerr << "Error: Could not find the TTree in the file." << std::endl;
      file->Close();
      return;
  }

  TClonesArray *jets_array = new TClonesArray("MyJet");
  TClonesArray *quarks = new TClonesArray("MyQuark");

  ttree->SetBranchAddress("jets_array", &jets_array);
  ttree->SetBranchAddress("quarks", &quarks);

  Long64_t ne = ttree->GetEntries();

  // Loop equivalente ao loop de eventos
  for ( Long64_t ni = 0; ni < ne; ni++)
  {
    ttree->GetEntry(ni);

    for (Int_t nj = 0; nj < jets_array->GetEntries(); nj++)
    {

      //---------------------------------------------------------------------------------------------------------
      // Processamento das informacoes referentes aos jatos usando MyJet
      //---------------------------------------------------------------------------------------------------------

      MyJet *jet = static_cast<MyJet *> (jets_array->At(nj));

      jetPt = jet->jetPt;
      jetEta = jet->jetEta;
      jetPhi = jet->jetPhi;
      jetMass = jet->jetMass;
      jetPx = jet->jetPx;
      jetPy = jet->jetPy;
      jetPz = jet->jetPz;
      jetE = jet->jetE;
      nConstituent = jet->nConstituent;
      pTLeadConstituent = jet->pTLeadConstituent;

      // Falta o numero de jatos por evento

      pT_Jet->Fill(jetPt);
      rapidity_Jet->Fill(jetEta);
      phi_Jet->Fill(jetPhi);
      mass_Jet->Fill(jetMass);
      nConstituent_Jet->Fill(nConstituent);
      pTLeadConstituent_pTJet->Fill(jetPt, pTLeadConstituent);
      nConstituent_Mass->Fill(nConstituent, jetMass);

      //---------------------------------------------------------------------------------------------------------
      // Comparacao dos jatos com os quarks
      //---------------------------------------------------------------------------------------------------------

      for (Int_t nj2 = 0; nj2 < quarks->GetEntries(); nj2++)
      {
        MyQuark *lq = static_cast<MyQuark *>(quarks->At(nj2));

        lastQuark_pT = lq->lqpT;
        lastQuark_Eta = lq->lqEta;
        lastQuark_Phi = lq->lqPhi;

        Double_t distancia = TMath::Sqrt(TMath::Power(jetPhi - lastQuark_Phi, 2) + TMath::Power(jetEta - lastQuark_Eta, 2));

        if (distancia <= 0.1 && jetPt > 1)
        {
          // cout << distancia << lastQuark_pT << jetPt << endl;
          quarks_and_Jets_pTs->Fill(lastQuark_pT, jetPt);
        }

      }
    }
    jets_array->Clear();
    quarks->Clear();
  }

  TCanvas *c1 = new TCanvas("c1", "Histograms and distributions", 2500, 2500);
  c1->Divide(3, 2);
  /*
  c1->cd(7);
  nJets_event ->SetTitle("Number of Jets per Event");
  nJets_event ->GetXaxis()->SetTitle("Number of Jets");
  nJets_event ->GetYaxis()->SetTitle("Frequency");
  nJets_event ->Draw();
  */
  c1->cd(1);
  pT_Jet->SetTitle("Jet's Transverse Momentum");
  pT_Jet->GetXaxis()->SetTitle("Transverse Momentum [GeV]");
  pT_Jet->GetYaxis()->SetTitle("Frequency");
  pT_Jet->Draw();

  c1->cd(2);
  rapidity_Jet->SetTitle("Jet's Rapidity");
  rapidity_Jet->GetXaxis()->SetTitle("Rapidity");
  rapidity_Jet->GetYaxis()->SetTitle("Frequency");
  rapidity_Jet->Draw();

  c1->cd(3);
  phi_Jet->SetTitle("Jet's Phi");
  phi_Jet->GetXaxis()->SetTitle("Phi");
  phi_Jet->GetYaxis()->SetTitle("Frequency");
  phi_Jet->Draw();

  c1->cd(4);
  mass_Jet->SetTitle("Jet's Mass");
  mass_Jet->GetXaxis()->SetTitle("Mass");
  mass_Jet->GetYaxis()->SetTitle("Frequency");
  mass_Jet->Draw();

  c1->cd(5);
  nConstituent_Jet->SetTitle("Constituents' Number per Jet");
  nConstituent_Jet->GetXaxis()->SetTitle("Number of Constituents");
  nConstituent_Jet->GetYaxis()->SetTitle("Frequency");
  nConstituent_Jet->Draw();

  c1->cd(6);
  pTLeadConstituent_pTJet->SetTitle("Graph of lead constituent's pT by jet's pT");
  pTLeadConstituent_pTJet->GetXaxis()->SetTitle("Jet's pT [GeV]");
  pTLeadConstituent_pTJet->GetYaxis()->SetTitle("Lead Constituent pT [GeV]");
  pTLeadConstituent_pTJet->Draw();

  TCanvas *c2 = new TCanvas("c2", "Histograms and distributions", 2500, 2500);
  c2->Divide(1, 2);

  c2->cd(1);
  nConstituent_Mass->SetTitle("Constituents' Number by Mass");
  nConstituent_Mass->GetXaxis()->SetTitle("Number of Constituents");
  nConstituent_Mass->GetYaxis()->SetTitle("Mass");
  nConstituent_Mass->Draw();

  c2->cd(2);
  quarks_and_Jets_pTs->SetTitle("Quarks pTs by jets pTs");
  quarks_and_Jets_pTs->GetXaxis()->SetTitle("Quarks pTs");
  quarks_and_Jets_pTs->GetYaxis()->SetTitle("Jets pTs");
  quarks_and_Jets_pTs->Draw();

  file->Close();
}

//-------------------------------------------------------------------------------------------------------------------------
// Funcao que determina o ultimo quark da mesma cor do quark que prescindiu os decaimentos sucessivos antes da hadronizacao
//-------------------------------------------------------------------------------------------------------------------------

Int_t findLastQuark(TClonesArray* particles, Int_t index = -1)
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
  Int_t iLastQuark=index;

  while(pdgi == pdgf)
  {
    part = (TParticle*)particles->At(iLastQuark);

    Int_t fd = part->GetFirstDaughter();
    Int_t ld = part->GetLastDaughter();

    for(int i=fd; i<=ld; i++)
    {
      TParticle *partd = (TParticle*)particles->At(i);

      if (!partd)
        continue;

      Int_t thisPdg = partd->GetPdgCode();
      //pdgList.Append(Form(";%d",thisPdg));
      if(thisPdg == pdgi)
        iLastQuark = i;
    }
        //if(pdgList.Contains(Form(";%d",pdgi)))
          pdgf = 0;
  }
  return iLastQuark;
}