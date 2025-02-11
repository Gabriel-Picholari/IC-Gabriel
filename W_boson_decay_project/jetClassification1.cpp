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
        JetInfo(const TString& type = "", const Int_t& pdg = 0, const Int_t& motherPdg = 0) : signalType(type), finalParticlePdg(pdg), finalParticleMotherPdg(motherPdg) {}

        void setSignalType(const TString& type) { signalType = type; }
        TString getSignalType() const { return signalType; }

        void setFinalParticlePdg(Int_t pdg) { finalParticlePdg = pdg; }
        Int_t getFinalParticlePdg() const { return finalParticlePdg; }

        void setFinalParticleMotherPdg(Int_t motherPdg) { finalParticleMotherPdg = motherPdg; }
        Int_t getFinalParticleMotherPdg() const { return finalParticleMotherPdg; }

    private:
        TString signalType;
        Int_t finalParticlePdg;
        Int_t finalParticleMotherPdg;
};

void jetClassification1(const char* fileName)
{

    gSystem->Load("libEG");
    gSystem->Load("libEGPythia8");

    //---------------------------------------------------------------------------------------------------------
    // Initialization of variables
    //---------------------------------------------------------------------------------------------------------

    Float_t fpPt, fpEta, fpPhi, fpE, fpPx, fpPy, fpPz, fpMass = 0;
    Float_t jetPt, jetEta, jetPhi, jetE, jetPx, jetPy, jetPz, jetMass, jetNConst, pT_LeadConst = 0;
    Float_t maxRho, nVert;
    TString signalType = "";
    Int_t finalParticlePdg = 0;
    Int_t finalParticleMotherPdg = 0;
    
    TLorentzVector vec_s(0,0,0,0);
    TLorentzVector vec_c(0,0,0,0);

    //---------------------------------------------------------------------------------------------------------
    // Initialization of histograms
    //---------------------------------------------------------------------------------------------------------

    TH1F *invariantMass = new TH1F("h1", "W^{+-} invariant mass spectrum [GeV/c^{2}]", 100, 0, 100);
    TH1F *primary_CharmRatioHist = new TH1F("primaryCharmRatio", "Charm: Ratio between tagged jets (with correct hadron) and quark p_{T}", 50, 0, 2);
    TH1F *secondary_CharmRatioHist = new TH1F("secondaryCharmRatio", "Charm: Ratio between tagged jets (whitout correct hadron) and quark p_{T}", 50, 0, 2);

    TH1F *primary_StrangeRatioHist = new TH1F("primaryStrangeRatio", "Strange: Ratio between tagged jets (with correct hadron) and quark p_{T}", 50, 0, 2);
    TH1F *secondary_StrangeRatioHist = new TH1F("secondaryStrangeRatio", "Strange: Ratio between tagged jets (whitout correct hadron) and quark p_{T}", 50, 0, 2);

    //---------------------------------------------------------------------------------------------------------
    // Initializations and FastJet configurations:
    //---------------------------------------------------------------------------------------------------------

    Float_t jetR = 0.4;

    fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, jetR);

    std::vector<fastjet::PseudoJet> particles_fastjet;
    std::vector<fastjet::PseudoJet> jets;

    //---------------------------------------------------------------------------------------------------------
    // Initialization of the .root file and TTrees
    //---------------------------------------------------------------------------------------------------------

    TFile *file = TFile::Open(fileName, "READ");
    TTree *ttree = dynamic_cast<TTree *>(file->Get("W decay TTree 2"));

    TClonesArray *jets_array = new TClonesArray("MyJet");
    TClonesArray *quarks = new TClonesArray("MyQuark");

    ttree->SetBranchAddress("jets_array", &jets_array);
    ttree->SetBranchAddress("quarks", &quarks);

    Long64_t ne = ttree->GetEntries();

    //---------------------------------------------------------------------------------------------------------
    // Event loop equivalent
    //---------------------------------------------------------------------------------------------------------

    for ( Long64_t ni = 0; ni < ne; ni++)
    {

        //std::cout << std::endl;
        //std::cout << "NEW EVENT ITERATION" << std::endl;
        //std::cout << std::endl;

        ttree->GetEntry(ni);

        vec_s.SetPtEtaPhiM(10,0,0,(3.141592));
        vec_c.SetPtEtaPhiM(10,0,0,(3.141592));

        std::vector<fastjet::PseudoJet> tagged_c_jets;
        std::vector<fastjet::PseudoJet>tagged_s_jets;

        Int_t count = 0;

        //---------------------------------------------------------------------------------------------------------
        // Quark information ( it's unique per event - a single pair of cbar(c) - s(sbar) quarks )
        //---------------------------------------------------------------------------------------------------------

        Float_t charmPt = 0;
        Float_t strangePt = 0;

        //std::cout << quarks->GetEntries() << std::endl;
        for (Int_t nk = 0; nk < quarks->GetEntries(); nk++)
        {
            MyQuark *mq = static_cast<MyQuark *>(quarks->At(nk));
            Int_t quarkPdg = mq->qPdg;
            
            if ( abs(quarkPdg) == 4) charmPt = mq->qpT;
            //std::cout << charmPt << std::endl;
            if ( abs(quarkPdg) == 3) strangePt = mq->qpT;
            //std::cout << strangePt << std::endl;                  // Check over that later (not a serious problem I suppose)

        }

        //---------------------------------------------------------------------------------------------------------
        // Particle loop equivalent
        //---------------------------------------------------------------------------------------------------------

        for (Int_t nj = 0; nj < jets_array->GetEntries(); nj++)
        {
            MyJet *fp = static_cast<MyJet *>(jets_array->At(nj));

            fpPx = fp->fPx;
            fpPy = fp->fPy;
            fpPz = fp->fPz;
            fpE = fp->fE;
            signalType = fp->signalType;
            finalParticlePdg = fp->finalParticlePdg;
            finalParticleMotherPdg = fp->finalParticleMotherPdg;

            fastjet::PseudoJet particle(fpPx, fpPy, fpPz, fpE);
            
            JetInfo* jetInfo = new JetInfo(signalType, finalParticlePdg, finalParticleMotherPdg);
            particle.set_user_info(jetInfo);

            particles_fastjet.push_back(particle);
        }

        fastjet::ClusterSequence clusterSeq(particles_fastjet, jet_def);
        jets = clusterSeq.inclusive_jets();

        for (const fastjet::PseudoJet& jet : jets)
        {
            //---------------------------------------------------------------------------------------------------------
            // Collecting jet's general information
            //---------------------------------------------------------------------------------------------------------

            //std::cout << "Jet number: " << count << std::endl;
            count++;

            jetPt = jet.pt();
            jetEta = jet.eta();
            jetPhi = jet.phi();
            jetMass = jet.m();
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

            TLorentzVector currentJet(jetPx, jetPy, jetPz, jetE);
            if (currentJet.M() < 0) continue;
            
            Bool_t isCharmTagged = false;
            Bool_t isStrangeTagged = false;

            for (const fastjet::PseudoJet &constituent : jet.constituents())
            {
                TString signalType_jet = constituent.user_info<JetInfo>().getSignalType();

                if (signalType_jet == "charm" && !isCharmTagged) // Then one of the final particles that is part of this jet has a distant mother in a c quark
                {   
                    tagged_c_jets.push_back(jet);
                    vec_c = TLorentzVector(jetPx, jetPy, jetPz, jetE); // Just to check consistency in the number of entries (will be reconsidered)
                    //std::cout << "Added charm jet with mass: " << currentJet.M() << std::endl;
                    isCharmTagged = true;
                }
                else if (signalType_jet == "strange" && !isStrangeTagged) // Then it's a strange jet
                {
                    tagged_s_jets.push_back(jet);
                    vec_s = TLorentzVector(jetPx, jetPy, jetPz, jetE); // Just to check consistency in the number of entries (will be reconsidered)
                    //std::cout << "Added strange jet with mass: " << currentJet.M() << std::endl;
                    isStrangeTagged = true;
                }
            }

        } // End of individual jet creation

        std::cout << "\n--- c-tagged Jets ---" << std::endl;
        for (const fastjet::PseudoJet  &jet : tagged_c_jets)
        {
            TLorentzVector cJet(jet.px(), jet.py(), jet.pz(), jet.E());

            std::cout << "New Jet" << std::endl;
            std::cout << "Jet Mass: " << cJet.M() << std::endl;
            std::cout << "Jet pT: " << cJet.Pt() << std::endl;
            Float_t charmRatio = jet.pt() / charmPt;
            std::cout << "Charm ratio: " << charmRatio << std::endl;

            for (const fastjet::PseudoJet &constituent : jet.constituents())
            {
                Int_t constituentPdg = constituent.user_info<JetInfo>().getFinalParticlePdg();
                Int_t abs_constituentPdg = abs(constituentPdg);

                if (abs_constituentPdg == 130 || abs_constituentPdg == 310 || abs_constituentPdg == 311 || abs_constituentPdg == 321 || abs_constituentPdg == 313 || abs_constituentPdg ==323 || abs_constituentPdg == 315 || abs_constituentPdg == 325 || abs_constituentPdg == 317 || abs_constituentPdg == 327 || abs_constituentPdg == 319 || abs_constituentPdg == 329 )
                {
                    primary_CharmRatioHist->Fill(charmRatio);
                }
                else
                {
                    secondary_CharmRatioHist->Fill(charmRatio);
                }

                std::cout << constituentPdg << std::endl;
            }
        }

        std::cout << "\n--- s-tagged Jets ---" << std::endl;
        for (const fastjet::PseudoJet &jet : tagged_s_jets)
        {
            TLorentzVector sJet(jet.px(), jet.py(), jet.pz(), jet.E());

            std::cout << "----- New Jet -----" << std::endl;
            std::cout << "Jet Mass: " << sJet.M() << std::endl;
            std::cout << "Jet pT: " << sJet.Pt() << std::endl;

            Float_t strangeRatio = jet.pt() / strangePt;
            std::cout << "Stange ratio: " << strangeRatio << std::endl;

            for (const fastjet::PseudoJet &constituent : jet.constituents())
            {
                Int_t constituentPdg = constituent.user_info<JetInfo>().getFinalParticlePdg();
                Int_t abs_constituentPdg = abs(constituentPdg);

                if (abs_constituentPdg == 130 || abs_constituentPdg == 310 || abs_constituentPdg == 311 || abs_constituentPdg == 321 || abs_constituentPdg == 313 || abs_constituentPdg == 323 || abs_constituentPdg == 315 || abs_constituentPdg == 325 || abs_constituentPdg == 317 || abs_constituentPdg == 327 || abs_constituentPdg == 319 || abs_constituentPdg == 329 )
                {
                    primary_StrangeRatioHist->Fill(strangeRatio);
                }
                else
                {
                    secondary_StrangeRatioHist->Fill(strangeRatio);
                }
                
                std::cout << constituentPdg << std::endl;
            }
        }

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
    // Plotting histograms
    //---------------------------------------------------------------------------------------------------------

    TCanvas *c1 = new TCanvas("c1", "Invariant mass distribution", 2500, 2500);
    c1->Divide(1, 1);

    c1->cd(1);
    invariantMass->SetTitle("Jet's invariant mass spectrum");
    invariantMass->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
    invariantMass->GetYaxis()->SetTitle("Frequency");
    invariantMass->Draw();

    TCanvas *c2 = new TCanvas("c2", "Charm", 2500, 2500);
    c2->Divide(1, 1);

    c2->cd(1);
    primary_CharmRatioHist->SetTitle("Charm jet and quark p_{T} ratio");
    primary_CharmRatioHist->GetXaxis()->SetTitle("Ratio");
    primary_CharmRatioHist->GetYaxis()->SetTitle("Frequency");
    primary_CharmRatioHist->SetLineColor(kRed);
    primary_CharmRatioHist->Draw();

    secondary_CharmRatioHist->SetLineColor(kBlue);
    secondary_CharmRatioHist->Draw("same");

    TCanvas *c3 = new TCanvas("c3", "Strange", 2500, 2500);
    c3->Divide(1, 1);

    c3->cd(1);
    primary_StrangeRatioHist->SetTitle("Strange jet and quark p_{T} ratio");
    primary_StrangeRatioHist->GetXaxis()->SetTitle("Ratio");
    primary_StrangeRatioHist->GetYaxis()->SetTitle("Frequency");
    primary_StrangeRatioHist->SetLineColor(kRed);
    primary_StrangeRatioHist->Draw();
    
    secondary_StrangeRatioHist->SetLineColor(kBlue);
    secondary_StrangeRatioHist->Draw("same");

    file->Close();
}