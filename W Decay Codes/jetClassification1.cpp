//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
/*

  This code determines the cuts for multiple potential discriminating variables through cumulative histograms. The differences from the previous version include the addition of a new discriminating 
variable and the creation of cumulative histograms for each of them. This macro is intended to be used in conjunction with wdecayTTree2.cpp.

*/
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#include <cmath>
#include <string>
#include "TH1.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TMath.h"
#include "MyJet.h"
#include "TFile.h"
#include <iostream>
#include "TSystem.h"
#include "MyQuark.h"
#include "TCanvas.h"
#include "TPythia8.h"
#include "TParticle.h"
#include "TClonesArray.h"
#include <TLorentzVector.h>
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

class JetInfo : public fastjet::PseudoJet::UserInfoBase
{
    public:
        JetInfo(const TString& type = "") : signalType(type) {}

        void setSignalType(const TString& type) { signalType = type; }
        TString getSignalType() const { return signalType; }

    private:
    TString signalType;
};

void jetClassification1(const char* fileName)
{

    gSystem->Load("libEG");
    gSystem->Load("libEGPythia8");

    //---------------------------------------------------------------------------------------------------------
    //Inicializacao das variaveis
    //---------------------------------------------------------------------------------------------------------

    Float_t fpPt, fpEta, fpPhi, fpE, fpPx, fpPy, fpPz, fpMass = 0;
    Float_t jetPt, jetEta, jetPhi, jetE, jetPx, jetPy, jetPz, jetMass, jetNConst, pT_LeadConst = 0;
    Float_t maxRho, nVert;
    TString signalType = "";

    // Criacao dos TLorentzVector

    TLorentzVector vec_s(0,0,0,0);
    TLorentzVector vec_c(0,0,0,0);
    //TLorentzVector vec_sbar(0,0,0,0); MAY BE USED LATER, FOR NOW, WE DO NOT NEED TO DISTINGUISH BETWEEN C AND CBAR or S AND SBAR
    //TLorentzVector vec_cbar(0,0,0,0);

    //---------------------------------------------------------------------------------------------------------
    // Inicializacao dos histogramas
    //---------------------------------------------------------------------------------------------------------

    TH1F *invariantMass = new TH1F("h1", "W^{+-} invariant mass spectrum [GeV/c^{2}]", 100, 0, 100);

    //---------------------------------------------------------------------------------------------------------
    // Inicializacoes e configuracoes do FastJet:
    //---------------------------------------------------------------------------------------------------------

    Float_t jetR = 0.4;

    fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, jetR);

    std::vector<fastjet::PseudoJet> particles_fastjet;
    std::vector<fastjet::PseudoJet> jets;

    //---------------------------------------------------------------------------------------------------------
    // Inicializacao do arquivo.root e das TTrees
    //---------------------------------------------------------------------------------------------------------

    TFile *file = TFile::Open(fileName, "READ");
    TTree *ttree = dynamic_cast<TTree *>(file->Get("W decay TTree 2"));

    TClonesArray *jets_array = new TClonesArray("MyJet");
    TClonesArray *quarks = new TClonesArray("MyQuark");

    ttree->SetBranchAddress("jets_array", &jets_array);
    ttree->SetBranchAddress("quarks", &quarks);

    Long64_t ne = ttree->GetEntries();

    //---------------------------------------------------------------------------------------------------------
    // Loop equivalente ao loop de eventos
    //---------------------------------------------------------------------------------------------------------

    for ( Long64_t ni = 0; ni < ne; ni++)
    {
        ttree->GetEntry(ni);

        vec_s.SetPtEtaPhiM(10,0,0,(3.141592));
        vec_c.SetPtEtaPhiM(10,0,0,(3.141592));
        //vec_sbar.SetPtEtaPhiM(10,0,0,(3.141592)); REASON FOR COMMENT ABOVE
        //vec_cbar.SetPtEtaPhiM(10,0,0,(3.141592));


        //---------------------------------------------------------------------------------------------------------
        // Loop equivalente ao loop de eventos
        //---------------------------------------------------------------------------------------------------------

        for (Int_t nj = 0; nj < jets_array->GetEntries(); nj++)
        {
            MyJet *fp = static_cast<MyJet *>(jets_array->At(nj));

            fpPx = fp->fPx;
            fpPy = fp->fPy;
            fpPz = fp->fPz;
            fpE = fp->fE;
            signalType = fp->signalType; // Retrieving the string label

            fastjet::PseudoJet particle(fpPx, fpPy, fpPz, fpE);
            
            particle.set_user_info(new JetInfo(signalType));

            particles_fastjet.push_back(particle);
        }

        fastjet::ClusterSequence clusterSeq(particles_fastjet, jet_def);
        jets = clusterSeq.inclusive_jets();

        for (const fastjet::PseudoJet& jet : jets)
        {
            jetPt = jet.pt();
            jetEta = jet.eta();
            jetPhi = jet.phi();
            jetMass = jet.m();  // Jet individual invariant mass
            jetPx = jet.px();
            jetPy = jet.py();
            jetPz = jet.pz();
            jetE = jet.E();
            jetNConst = jet.constituents().size();
            
            pT_LeadConst = 0;
            for (const fastjet::PseudoJet &constituent : jet.constituents())
            {
                if (constituent.pt() > pT_LeadConst)
                {
                    pT_LeadConst = constituent.pt();
                }
            }           

            Float_t angAve, sigmaKT = 0;

            for (Int_t i = 0; i < jetNConst; ++i) 
            {
                Double_t pt_constituentes = jet.constituents()[i].pt();
                Double_t eta_constituentes = jet.constituents()[i].eta();
                Double_t phi_constituentes = jet.constituents()[i].phi();

                Double_t angPart = TMath::Sqrt(TMath::Power(TMath::Abs(phi_constituentes) - TMath::Abs(jetPhi), 2) + TMath::Power(eta_constituentes - jetEta, 2));

                angAve = angAve + angPart;
                sigmaKT = sigmaKT + (TMath::Power(pt_constituentes - (jetPt / jetNConst), 2));
            }

            angAve = angAve / jetNConst;
            sigmaKT = sigmaKT / jetNConst;               

            //---------------------------------------------------------------------------------------------------------
            // Jet classification block (based on constituents info)
            //---------------------------------------------------------------------------------------------------------

            for (const fastjet::PseudoJet &constituent : jet.constituents())
            {
                TString signalType_jet = constituent.user_info<JetInfo>().getSignalType();

                if (signalType_jet == "charm") // Then one of the final particles that is part of this jet has a distant mother in a c quark
                {
                    vec_c = TLorentzVector(jetPx, jetPy, jetPz, jetE);
                }
                else if (signalType_jet == "strange") // Then it's a strange jet
                {
                    vec_s = TLorentzVector(jetPx, jetPy, jetPz, jetE);
                }
            }
        
            //---------------------------------------------------------------------------------------------------------

        } // End of individual jet creation

        const Float_t tolerance = 1e-5;
        
        if ( fabs(vec_c.M() - 3.141592) > tolerance && fabs(vec_s.M() - 3.141592) > tolerance ) 
        {
            TLorentzVector vec_W = vec_c + vec_s;
            invariantMass->Fill( vec_W.M() ); 
        }

        particles_fastjet.clear();
        jets.clear();
        jets_array->Clear();
        quarks->Clear();

    
    } // End of event loop equivalent

    //---------------------------------------------------------------------------------------------------------
    // Creation of the cumulative histograms regarding the classical analysis target object
    //---------------------------------------------------------------------------------------------------------

    // MISSING
    
    //---------------------------------------------------------------------------------------------------------
    // Plotting histograms
    //---------------------------------------------------------------------------------------------------------

    TCanvas *c1 = new TCanvas("c1", "Invariant mass distribution", 2500, 2500);
    c1->Divide(1, 2);

    c1->cd(1);
    invariantMass->SetTitle("Jet's invariant mass spectrum");
    invariantMass->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
    invariantMass->GetYaxis()->SetTitle("Frequency");
    invariantMass->Draw();



    file->Close();
}