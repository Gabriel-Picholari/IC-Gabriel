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
        JetInfo(const Float_t& vx = 0, const Float_t& vy = 0, const Float_t& vz = 0, const TString& type = "", const Int_t& pdg = 0, const Int_t& motherPdg = 0, const Int_t& secondMotherPdg = 0, const Int_t& thirdMotherPdg = 0) : fpVx(vx), fpVy(vy), fpVz(vz), signalType(type), finalParticlePdg(pdg), finalParticleMotherPdg(motherPdg), finalParticleSecondMotherPdg(secondMotherPdg), finalParticleThirdMotherPdg(thirdMotherPdg){}

        void setVx(Float_t vx) { fpVx = vx; }
        Float_t getVx() const { return fpVx; }

        void setVy(Float_t vy) { fpVy = vy; }
        Float_t getVy() const { return fpVy; }

        void setVz(Float_t vz) { fpVz = vz; }
        Float_t getVz() const { return fpVz; }
        
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
        Float_t fpVx;
        Float_t fpVy;
        Float_t fpVz;
        TString signalType;
        Int_t finalParticlePdg;
        Int_t finalParticleMotherPdg;
        Int_t finalParticleSecondMotherPdg;
        Int_t finalParticleThirdMotherPdg;
};

void jetClassification2_charm(const char* fileName)
{

    gSystem->Load("libEG");
    gSystem->Load("libEGPythia8");

    //---------------------------------------------------------------------------------------------------------
    // Initialization of variables
    //---------------------------------------------------------------------------------------------------------

    Float_t fpPt, fpEta, fpPhi, fpE, fpPx, fpPy, fpPz, fpMass, fpVx, fpVy, fpVz = 0;
    Float_t jetPt, jetEta, jetPhi, jetE, jetPx, jetPy, jetPz, jetMass, jetNConst, pT_LeadConst = 0;
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
    // Initialization of TMVA TTrees
    //---------------------------------------------------------------------------------------------------------

    Float_t eventID_c, pT_c, label_c, nConst_c, eta_c, phi_c, mass_c = 0;

    TFile *filteredDataFile = new TFile("filteredOutput_2var_modelTraining_charm.root", "RECREATE");

    TTree *signalTree_c = new TTree("SignalTree_c", "Tree with signal data from c quark");
    signalTree_c->Branch("pT_c", &pT_c);
    signalTree_c->Branch("eta_c", &eta_c);
    signalTree_c->Branch("phi_c", &phi_c);
    signalTree_c->Branch("mass_c", &mass_c);
    signalTree_c->Branch("label_c", &label_c);
    signalTree_c->Branch("nConst_c", &nConst_c);
    signalTree_c->Branch("eventID_c", &eventID_c);

    TTree *backgroundTree_c = new TTree("BackgroundTree_c", "Tree with background data from c quark");
    backgroundTree_c->Branch("pT_c", &pT_c);
    backgroundTree_c->Branch("eta_c", &eta_c);
    backgroundTree_c->Branch("phi_c", &phi_c);
    backgroundTree_c->Branch("mass_c", &mass_c);
    backgroundTree_c->Branch("label_c", &label_c);
    backgroundTree_c->Branch("nConst_c", &nConst_c);
    backgroundTree_c->Branch("eventID_c", &eventID_c);

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

            fpVx = fp->fVx;
            fpVy = fp->fVy;
            fpVz = fp->fVz; 
            fpPx = fp->fPx;
            fpPy = fp->fPy;
            fpPz = fp->fPz;
            fpE = fp->fE;
            
            //---------------------------------------------------------------------------------------------------------
            // Adding new variables to each jet
            //---------------------------------------------------------------------------------------------------------

            signalType = fp->signalType;
            finalParticlePdg = fp->finalParticlePdg;
            finalParticleMotherPdg = fp->finalParticleMotherPdg;
            finalParticleSecondMotherPdg = fp->finalParticleSecondMotherPdg;
            finalParticleThirdMotherPdg = fp->finalParticleThirdMotherPdg;

            fastjet::PseudoJet particle(fpPx, fpPy, fpPz, fpE);
            
            JetInfo* jetInfo = new JetInfo(fpVx, fpVy, fpVz, signalType, finalParticlePdg, finalParticleMotherPdg, finalParticleSecondMotherPdg, finalParticleThirdMotherPdg);
            particle.set_user_info(jetInfo);

            //---------------------------------------------------------------------------------------------------------

            particles_fastjet.push_back(particle);
        }

        fastjet::ClusterSequence clusterSeq(particles_fastjet, jet_def);
        jets = clusterSeq.inclusive_jets();

        for (const fastjet::PseudoJet& jet : jets)
        {
            //---------------------------------------------------------------------------------------------------------
            // Collecting jet's general information (individually and as soon as they are created)
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
            // Jet classification block (based on the extra constituents info we previously added)
            //---------------------------------------------------------------------------------------------------------

            TLorentzVector currentJet(jetPx, jetPy, jetPz, jetE);
            if (currentJet.M() < 0) continue; 
            // I observed some abberrations while looking at jets: some originated by FastJet had negative mass and they were taken as physically irrelevant
            
            Bool_t isCharmTagged = false;   // One constituent being tagged as either charmed or strange is sufficient for tagging the jet
            Bool_t isStrangeTagged = false; // Notice that since we do it per jet, the bools do not get redundant and are automatically reseted

            for (const fastjet::PseudoJet &constituent : jet.constituents())
            {
                TString signalType_jet = constituent.user_info<JetInfo>().getSignalType();

                if (signalType_jet == "charm" && !isCharmTagged) // Then one of the final particles that is part of this jet has a distant mother in a c quark. We tag it!
                {   
                    tagged_c_jets.push_back(jet);
                    vec_c = TLorentzVector(jetPx, jetPy, jetPz, jetE);
                    isCharmTagged = true; // Tagged - we've got a label
                }
                else if (signalType_jet == "strange" && !isStrangeTagged) // Then it's a strange jet
                {
                    tagged_s_jets.push_back(jet);
                    vec_s = TLorentzVector(jetPx, jetPy, jetPz, jetE);
                    isStrangeTagged = true; // Analogous
                }
            }

            if (!isCharmTagged) // Then it' a background jet, that is, not a charmed jet (can be strange or any other jet)
            { 
                label_c = 0;
                eventID_c = ni;
                pT_c = jetPt;
                eta_c = jetEta;
                phi_c = jetPhi;
                mass_c = jetMass;
                nConst_c = jetNConst;
                backgroundTree_c->Fill();

                // Since we're dealing in this specific macro only with charmed jets, for a model dedicated to predicting charmed jets, anything else, that is, stranged jets
                // or any other kind of jet is to be considered background.
            }

        } // End of individual jet creation

        const std::unordered_set<int> charmPdgSet = {411, 421, 413, 423, 415, 425, 431, 433, 435, 10411, 10421, 413, 423, 10413, 10423, 20413, 20423, 415, 425, 431, 10431, 433, 10433, 20433, 435, 4122, 4222, 4212, 4112, 4224, 4214, 4114, 4232, 4132, 4322, 4312, 4324, 4314, 4332, 4334, 4412, 4422, 4414, 4424, 4432, 4434, 4444};

        //---------------------------------------------------------------------------------------------------------
        // Charmed jets tagging
        //---------------------------------------------------------------------------------------------------------

        for (const fastjet::PseudoJet  &jet : tagged_c_jets) // Opening the jet vector: the analysis object is a jet with a constituent found to have ancestry in a charm
        {
            TLorentzVector cJet(jet.px(), jet.py(), jet.pz(), jet.E());
            Float_t charmRatio = jet.pt() / charmPt;
            
            Bool_t hasCharmConstituent = false;

            for (const fastjet::PseudoJet &constituent : jet.constituents()) // Opening the jet itself: the analysis object is a constituent of the jet
            {
                Int_t constituentPdg                = constituent.user_info<JetInfo>().getFinalParticlePdg();
                Int_t constituentMotherPdg          = constituent.user_info<JetInfo>().getFinalParticleMotherPdg();
                Int_t constituentSecondMotherPdg    = constituent.user_info<JetInfo>().getFinalParticleSecondMotherPdg();
                Int_t constituentThirdMotherPdg     = constituent.user_info<JetInfo>().getFinalParticleThirdMotherPdg();

                Int_t abs_constituentPdg                = abs(constituentPdg);
                Int_t abs_constituentMotherPdg          = abs(constituentMotherPdg);
                Int_t abs_constituentSecondMotherPdg    = abs(constituentSecondMotherPdg);
                Int_t abs_constituentThirdMotherPdg     = abs(constituentThirdMotherPdg);

                // We found it to be necessary to go back to the third consecutive mother of the final particle because some intermediate decays made the presence of
                // hadrons on the list to be "hard to see". We are here trying to make sure that the final particle we said to have ancestry on a charm taht comes from 
                // the decay of W bósons is indeed a particle that comes from the decay, in first ( or till third ) instance of a confirmed charmed hadron.

                if (charmPdgSet.count(abs_constituentPdg) || charmPdgSet.count(abs_constituentMotherPdg) || charmPdgSet.count(abs_constituentSecondMotherPdg) || charmPdgSet.count(abs_constituentThirdMotherPdg))
                {
                    hasCharmConstituent = true; // If thats the case, then we welcome this jets to the list of charmed jets that will label the trainment
                    break;
                }
            }
            if (hasCharmConstituent) // Signal data
            {
                label_c = 1;
                eventID_c = ni;

                //Kinematics directly from PseudoJet
                pT_c = jet.pt();
                eta_c = jet.eta();
                phi_c = jet.phi();
                mass_c = jet.m();

                nConst_c = (Int_t)jet.constituents().size();

                /*  Just uncomment if leadingPt is required
                if (constituent.pt() > pT_LeadConst)
                {
                    pT_LeadConst = constituent.pt();
                }
                */

                signalTree_c->Fill();
            }
            // If our second confirmation fails, then I suppose we can throw this jet into background. I'll get confirmation over that supposition and then execute this idea
        }

        particles_fastjet.clear();
        jets.clear();
        jets_array->Clear();
        quarks->Clear();

    
    } // End of event loop equivalent
    signalTree_c->Write();
    backgroundTree_c->Write();

    //signalTree_s->Write();
    //backgroundTree_s->Write();

    file->Close();
    filteredDataFile->Close();
    //filteredDataFile_s->Close();
}

