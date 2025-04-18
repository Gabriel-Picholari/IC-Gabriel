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
#include <unordered_set>
#include "TClonesArray.h"
#include <TLorentzVector.h>
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

class JetInfo : public fastjet::PseudoJet::UserInfoBase
{
    public:
        JetInfo(const TString& type = "", const Int_t& pdg = 0, const Int_t& motherPdg = 0, const Int_t& secondMotherPdg = 0, const Int_t& thirdMotherPdg = 0) : signalType(type), finalParticlePdg(pdg), finalParticleMotherPdg(motherPdg), finalParticleSecondMotherPdg(secondMotherPdg), finalParticleThirdMotherPdg(thirdMotherPdg){}

        void setSignalType(const TString& type) { signalType = type; }
        TString getSignalType() const { return signalType; }

        void setFinalParticlePdg(Int_t pdg) { finalParticlePdg = pdg; }
        Int_t getFinalParticlePdg() const { return finalParticlePdg; }

        void setFinalParticleMotherPdg(Int_t motherPdg) { finalParticleMotherPdg = motherPdg; }
        Int_t getFinalParticleMotherPdg() const { return finalParticleMotherPdg; }

        void setFinalParticleSecondMotherPdg(Int_t secondMotherPdg) { finalParticleSecondMotherPdg = secondMotherPdg; }
        Int_t getFinalParticleSecondMotherPdg() const { return finalParticleSecondMotherPdg; }

        void setFinalParticleThirdMotherPdg(Int_t thirdMotherPdg) { finalParticleThirdMotherPdg = thirdMotherPdg; }
        Int_t getFinalParticleThirdMotherPdg() const { return finalParticleThirdMotherPdg; }

    private:
        TString signalType;
        Int_t finalParticlePdg;
        Int_t finalParticleMotherPdg;
        Int_t finalParticleSecondMotherPdg;
        Int_t finalParticleThirdMotherPdg;
};

void jetClassification2_strange(const char* fileName)
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
    Int_t finalParticleSecondMotherPdg = 0;
    Int_t finalParticleThirdMotherPdg = 0;
    
    TLorentzVector vec_s(0,0,0,0);
    TLorentzVector vec_c(0,0,0,0);

    //---------------------------------------------------------------------------------------------------------
    // Initializations and FastJet configurations:
    //---------------------------------------------------------------------------------------------------------

    Float_t jetR = 0.7;

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
    // Inicializacao das TTrees TMVA
    //---------------------------------------------------------------------------------------------------------

    Float_t eventID_s, pT_s, label_s, nConst_s, eta_s, phi_s, mass_s = 0;

    TFile *filteredDataFile = new TFile("filteredOutput_2var_modelApplication_strange.root", "RECREATE");

    TTree *signalTree_s = new TTree("SignalTree_s", "Tree with signal data from s quark");
    signalTree_s->Branch("pT_s", &pT_s);
    signalTree_s->Branch("eta_s", &eta_s);
    signalTree_s->Branch("phi_s", &phi_s);
    signalTree_s->Branch("mass_s", &mass_s);
    signalTree_s->Branch("label_s", &label_s);
    signalTree_s->Branch("nConst_s", &nConst_s);
    signalTree_s->Branch("eventID_s", &eventID_s);

    TTree *backgroundTree_s = new TTree("BackgroundTree_s", "Tree with background data from s quark");
    backgroundTree_s->Branch("pT_s", &pT_s);
    backgroundTree_s->Branch("eta_s", &eta_s);
    backgroundTree_s->Branch("phi_s", &phi_s);
    backgroundTree_s->Branch("mass_s", &mass_s);
    backgroundTree_s->Branch("label_s", &label_s);
    backgroundTree_s->Branch("nConst_s", &nConst_s);
    backgroundTree_s->Branch("eventID_s", &eventID_s);
    

    //---------------------------------------------------------------------------------------------------------
    // Event loop equivalent
    //---------------------------------------------------------------------------------------------------------

    for ( Long64_t ni = 0; ni < ne; ni++)
    {
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
            finalParticleSecondMotherPdg = fp->finalParticleSecondMotherPdg;
            finalParticleThirdMotherPdg = fp->finalParticleThirdMotherPdg;

            fastjet::PseudoJet particle(fpPx, fpPy, fpPz, fpE);
            
            JetInfo* jetInfo = new JetInfo(signalType, finalParticlePdg, finalParticleMotherPdg, finalParticleSecondMotherPdg, finalParticleThirdMotherPdg);
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
                    isCharmTagged = true;
                }
                else if (signalType_jet == "strange" && !isStrangeTagged) // Then it's a strange jet
                {
                    tagged_s_jets.push_back(jet);
                    vec_s = TLorentzVector(jetPx, jetPy, jetPz, jetE); // Just to check consistency in the number of entries (will be reconsidered)
                    isStrangeTagged = true;
                }
            }
            
            if (!isStrangeTagged) // Then it' a background jet, that is, not a strange jet (can be charmed or any other jet)
            {
                label_s = 0;
                eventID_s = ni;
                pT_s = jetPt;
                eta_s = jetEta;
                phi_s = jetPhi;
                mass_s = jetMass;
                nConst_s = jetNConst;
                backgroundTree_s->Fill();
            }

        } // End of individual jet creation

        const std::unordered_set<int> strangePdgSet = {130, 310, 311, 321, 313, 323, 315, 325, 317, 327, 319, 329, 9000311, 9000321, 10311, 10321, 100311, 100321, 9010311, 9010321, 9020311, 9020321, 313, 323, 10313, 10323, 20313, 20323, 100313, 100323, 9000313, 9000323, 30313, 30323, 315, 325, 9000315, 9000325, 10315, 10325, 20315, 20325, 9010315, 9010325, 9020315, 9020325, 317, 327, 9010317, 9010327, 319, 329, 3122, 3222, 3212, 3112, 3224, 3214, 3114, 3322, 3312, 3324, 3314, 3334};

        //---------------------------------------------------------------------------------------------------------
        // Strange jets tagging
        //---------------------------------------------------------------------------------------------------------
        
        for (const fastjet::PseudoJet &jet : tagged_s_jets)
        {
            TLorentzVector sJet(jet.px(), jet.py(), jet.pz(), jet.E());
            Float_t strangeRatio = jet.pt() / strangePt;

           Bool_t hasStrangeConstituent = false;
           Int_t missingStrangePdg = 0;

            for (const fastjet::PseudoJet &constituent : jet.constituents())
            {
                Int_t constituentPdg = constituent.user_info<JetInfo>().getFinalParticlePdg();
                Int_t constituentMotherPdg = constituent.user_info<JetInfo>().getFinalParticleMotherPdg();
                Int_t constituentSecondMotherPdg = constituent.user_info<JetInfo>().getFinalParticleSecondMotherPdg();
                
                Int_t abs_constituentPdg = abs(constituentPdg);
                Int_t abs_constituentSecondMotherPdg = abs(constituentSecondMotherPdg);


                if (strangePdgSet.count(abs_constituentPdg) || strangePdgSet.count(abs_constituentSecondMotherPdg))
                {
                    hasStrangeConstituent = true;
                    break;
                }
            }
            if (hasStrangeConstituent)
            {
                label_s = 1;
                eventID_s = ni;
                pT_s = jetPt;
                eta_s = jetEta;
                phi_s = jetPhi;
                mass_s = jetMass;
                nConst_s = jetNConst;
                signalTree_s->Fill();
            }
        }
        

        particles_fastjet.clear();
        jets.clear();
        jets_array->Clear();
        quarks->Clear();

    
    } // End of event loop equivalent

    signalTree_s->Write();
    backgroundTree_s->Write();

    file->Close();
    filteredDataFile->Close();
}

