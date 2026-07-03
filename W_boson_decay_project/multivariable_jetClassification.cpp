#include <map>
#include <cmath>
#include <string>
#include "TH1.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TMath.h"
#include "MyJet.h"
#include "TFile.h"
#include <iostream>
#include <algorithm>
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
        JetInfo(const Float_t& vx, const Float_t& vy, const Float_t& vz, const TString& type = "", const Int_t& pdg = 0, const Int_t& motherPdg = 0, const Int_t& secondMotherPdg = 0, const Int_t& thirdMotherPdg = 0) : fpVx(vx), fpVy(vy), fpVz(vz), signalType(type), finalParticlePdg(pdg), finalParticleMotherPdg(motherPdg), finalParticleSecondMotherPdg(secondMotherPdg), finalParticleThirdMotherPdg(thirdMotherPdg){}

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

struct JetTagResult
{
    Bool_t isCharmTagged = false;
    Bool_t isStrangeTagged = false;
};

JetTagResult classifyJet(
    const fastjet::PseudoJet& jet, 
    std::vector<fastjet::PseudoJet>& tagged_c_jets, 
    std::vector<fastjet::PseudoJet>& tagged_s_jets, 
    TLorentzVector& vec_c, 
    TLorentzVector& vec_s
)
{
    JetTagResult result;

    TLorentzVector currentJet(jet.px(), jet.py(), jet.pz(), jet.E());

    if (currentJet.M() < 0) return result;

    for (const fastjet::PseudoJet &constituent : jet.constituents())
    {
        TString signalType = constituent.user_info<JetInfo>().getSignalType();

        if (signalType == "charm" && !result.isCharmTagged)
        {
            tagged_c_jets.push_back(jet);
            vec_c = currentJet;
            result.isCharmTagged = true;
        }

        else if (signalType == "strange" && !result.isStrangeTagged)
        {
            tagged_s_jets.push_back(jet);
            vec_s = currentJet;
            result.isStrangeTagged = true;
        }

        if (result.isCharmTagged && result.isStrangeTagged)
        {
            break;
        }
    }

    return result;
}

struct JetVariables
{
    Float_t first_nRho = 0;
    Float_t second_nRho = 0;
    Float_t third_nRho = 0;

    Float_t leadingPt = 0;

    std::vector<Float_t> vertices_masses;
};

JetVariables extractJetVariables(
    const fastjet::PseudoJet& jet,
    Float_t rhoLowerBound,
    Float_t firstRhoUpperBound,
    Float_t secondRhoUpperBound,
    Float_t thirdRhoUpperBound
)
{
    JetVariables vars;

    std::map<Float_t, std::vector<fastjet::PseudoJet>> grupos_por_rho;

    for (const fastjet::PseudoJet &constituent : jet.constituents())
    {
        Float_t vx = constituent.user_info<JetInfo>().getVx();
        Float_t vy = constituent.user_info<JetInfo>().getVy();

        Float_t rho = TMath::Sqrt(vx*vx + vy*vy);
        grupos_por_rho[rho].push_back(constituent);

        if (rho >= rhoLowerBound && rho < firstRhoUpperBound) vars.first_nRho++;

        if (rho >= rhoLowerBound && rho < secondRhoUpperBound) vars.second_nRho++;

        if (rho >= rhoLowerBound && rho < thirdRhoUpperBound)
        {
            vars.third_nRho++;
            //grupos_por_rho[rho].push_back(constituent);  // See May 7th log for more information about this comment
        }

        if (constituent.pt() > vars.leadingPt) vars.leadingPt = constituent.pt();
    }

    for (auto const& [rho, constituintes] : grupos_por_rho) 
    {
        TLorentzVector vertexMomentum(0, 0, 0, 0);
        for (const auto& c : constituintes) 
        {
            vertexMomentum += TLorentzVector(c.px(), c.py(), c.pz(), c.E());
        }
        vars.vertices_masses.push_back(vertexMomentum.M());
    }

    return vars;
}

/*
    Parameters:
    - switch_string: either 'strange' or 'charm', designates the type of signal being analysed
    - mode_string: either 'training' or 'testing', designates whether the function is being used to filter data for model training or model testing
    - contaminatingGluonMode: either 'include' or 'exclude', designates whether the output data shouls include or exclude the contaminating gluon soft jets
*/

void multivariable_jetClassification(const char* fileName, std::string switch_string, std::string mode_string, std::string contaminatingGluonMode)
{

    gSystem->Load("libEG");
    gSystem->Load("libEGPythia8");

    //---------------------------------------------------------------------------------------------------------
    // Initialization of variables
    //---------------------------------------------------------------------------------------------------------

    Float_t fpPt, fpEta, fpPhi, fpE, fpPx, fpPy, fpPz, fpMass, fpVx, fpVy, fpVz;
    Float_t jetPt, jetEta, jetPhi, jetE, jetPx, jetPy, jetPz, jetMass, jetNConst, pT_LeadConst;
    Float_t rhoLowerBound = 0.01;
    Float_t firstRhoUpperBound = 1;
    Float_t secondRhoUpperBound = 1.5;
    Float_t thirdRhoUpperBound = 2;
    TString signalType = "";
    Int_t finalParticlePdg = 0;
    Int_t finalParticleMotherPdg = 0;
    Int_t finalParticleSecondMotherPdg = 0;
    Int_t finalParticleThirdMotherPdg = 0;
    
    TLorentzVector vec_s(0,0,0,0);
    TLorentzVector vec_c(0,0,0,0);

    //---------------------------------------------------------------------------------------------------------
    // Histogramas
    //---------------------------------------------------------------------------------------------------------

    TH1F* first_signal_nRho_Distribution = new TH1F("firstNRho_signal_hist", ("Signal jets:"+ switch_string +" N_{*rho} distribution (upper bound = 1)").c_str(), 20, 0, 20);
    TH1F* second_signal_nRho_Distribution = new TH1F("secondNRho_signal_hist", ("Signal jets:"+ switch_string +" N_{*rho} distribution (upper bound = 1.5)").c_str(), 20, 0, 20);
    TH1F* third_signal_nRho_Distribution = new TH1F("thirdNRho_signal_hist", ("Signal jets:"+ switch_string +" N_{*rho} distribution (upper bound = 2)").c_str(), 20, 0, 20);

    TH1F* first_background_nRho_Distribution = new TH1F("firstNRho_background_hist", ("Background jets:"+ switch_string +" N_{*rho} distribution (upper bound = 1)").c_str(), 20, 0, 20);
    TH1F* second_background_nRho_Distribution = new TH1F("secondNRho_background_hist", ("Background jets:"+ switch_string +" N_{*rho} distribution (upper bound = 1.5)").c_str(), 20, 0, 20);
    TH1F* third_background_nRho_Distribution = new TH1F("thirdNRho_background_hist", ("Background jets:"+ switch_string +" N_{*rho} distribution (upper bound = 2)").c_str(), 20, 0, 20);

    TH1F* signal_pTDistribution = new TH1F("signal_pTDistribution", ("Signal jets:"+ switch_string +" p_{T} distribution").c_str(), 100, 0, 100);
    TH1F* background_pTDistribution = new TH1F("background_pTDistribution", ("Background jets:"+ switch_string +" p_{T} distribution").c_str(), 100, 0, 100);

    TH1F* signal_nConstDistribution = new TH1F("signal_nConstDistribution", ("Signal jets:"+ switch_string +" N_{const}").c_str(), 100, 0, 100);
    TH1F* background_nConstDistribution = new TH1F("background_nConstDistribution", ("Background jets:"+ switch_string +" N_{const}").c_str(), 100, 0, 100);

    TH1F* signal_jetConstituentsVertex_invariantMassDistribution = new TH1F("signal_jetConstituentsVertex_invariantMassDistribution", ("Signal jets:"+ switch_string +" invariant mass distribution of the jet constituents' vertices").c_str(), 100, 0, 50);
    TH1F* background_jetConstituentsVertex_invariantMassDistribution = new TH1F("background_jetConstituentsVertex_invariantMassDistribution", ("Background jets:"+ switch_string +" invariant mass distribution of the jet constituents' vertices").c_str(), 100, 0, 50);

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

    Float_t eventID, pT, label, nConst, eta, phi, nRho, first_nRho, second_nRho, third_nRho;
    std::vector<Float_t> jetVerticesInvariantMasses;

    std::string sfx;
    std::string long_sfx;
    std::string short_sfx;
    std::string uppercase_switch_string;

    if (switch_string == "strange")
    {
        sfx = "_s";
        long_sfx = "_strange";
        short_sfx = "s";
        uppercase_switch_string = "Strange";
    } 
    else if (switch_string == "charm")
    {
        sfx = "_c";
        long_sfx = "_charm";
        short_sfx = "c";
        uppercase_switch_string = "Charm";
    }

    TFile *filteredDataFile = nullptr;

    if (switch_string == "strange" && mode_string == "training")
    {
        if (contaminatingGluonMode == "include") filteredDataFile = new TFile("multivariable_filteredOutput_modelTraining_gluonJetsIncluded_strange.root", "RECREATE");
        if (contaminatingGluonMode == "exclude") filteredDataFile = new TFile("multivariable_filteredOutput_modelTraining_gluonJetsExcluded_strange.root", "RECREATE");
    } 
    if (switch_string == "strange" && mode_string == "testing") filteredDataFile = new TFile("multivariable_filteredOutput_modelTesting_strange.root", "RECREATE");

    if (switch_string == "charm" && mode_string == "training")
    {
        if (contaminatingGluonMode == "include") filteredDataFile = new TFile("multivariable_filteredOutput_modelTraining_gluonJetsIncluded_charm.root", "RECREATE");
        if (contaminatingGluonMode == "exclude") filteredDataFile = new TFile("multivariable_filteredOutput_modelTraining_gluonJetsExcluded_charm.root", "RECREATE");
    }
    if (switch_string == "charm" && mode_string == "testing") filteredDataFile = new TFile("multivariable_filteredOutput_modelTesting_charm.root", "RECREATE");

    TTree *signalTree = new TTree(("SignalTree" + sfx).c_str(), ("TTree with signal data from " + short_sfx + " quark").c_str());
    signalTree->Branch(("pT" + sfx).c_str(), &pT);
    signalTree->Branch(("eta" + sfx).c_str(), &eta);
    signalTree->Branch(("phi" + sfx).c_str(), &phi);
    signalTree->Branch(("label" + sfx).c_str(), &label);
    signalTree->Branch(("nConst" + sfx).c_str(), &nConst);
    signalTree->Branch(("nRho" + sfx).c_str(), &nRho);
    signalTree->Branch(("eventID" + sfx).c_str(), &eventID);
    signalTree->Branch(("jetVerticesInvariantMasses" + sfx).c_str(), &jetVerticesInvariantMasses);

    TTree *backgroundTree = new TTree(("BackgroundTree" + sfx).c_str(), ("TTree with background data from " + short_sfx + " quark").c_str());
    backgroundTree->Branch(("pT" + sfx).c_str(), &pT);
    backgroundTree->Branch(("eta" + sfx).c_str(), &eta);
    backgroundTree->Branch(("phi" + sfx).c_str(), &phi);
    backgroundTree->Branch(("label" + sfx).c_str(), &label);
    backgroundTree->Branch(("nConst" + sfx).c_str(), &nConst);
    backgroundTree->Branch(("nRho" + sfx).c_str(), &nRho);
    backgroundTree->Branch(("eventID" + sfx).c_str(), &eventID);
    backgroundTree->Branch(("jetVerticesInvariantMasses" + sfx).c_str(), &jetVerticesInvariantMasses);

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
            //std::cout << strangePt << std::endl;

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
            signalType = fp->signalType;
            finalParticlePdg = fp->finalParticlePdg;
            finalParticleMotherPdg = fp->finalParticleMotherPdg;
            finalParticleSecondMotherPdg = fp->finalParticleSecondMotherPdg;
            finalParticleThirdMotherPdg = fp->finalParticleThirdMotherPdg;

            fastjet::PseudoJet particle(fpPx, fpPy, fpPz, fpE);
            
            JetInfo* jetInfo = new JetInfo(fpVx, fpVy, fpVz, signalType, finalParticlePdg, finalParticleMotherPdg, finalParticleSecondMotherPdg, finalParticleThirdMotherPdg);
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
            if (jetPt < 5) continue;

            jetEta = jet.eta();

            Float_t absEta = TMath::Abs(jetEta);
            if (absEta > 2) continue;

            jetPhi = jet.phi();
            jetMass = jet.m();
            jetPx = jet.px();
            jetPy = jet.py();
            jetPz = jet.pz();
            jetE = jet.E();
            jetNConst = jet.constituents().size();
            
            first_nRho = 0; // Counter that will serve as the new discriminatory variable: it measures the amount of constituents, per jet, that have a vertex in a given interval [0, rhoUpperBound]
            second_nRho = 0;
            third_nRho = 0;
            pT_LeadConst = 0;
            
            JetVariables vars = extractJetVariables(
                jet,
                rhoLowerBound,
                firstRhoUpperBound,
                secondRhoUpperBound,
                thirdRhoUpperBound
            );

            first_nRho = vars.first_nRho;
            second_nRho = vars.second_nRho;
            third_nRho = vars.third_nRho;
            pT_LeadConst = vars.leadingPt;
            jetVerticesInvariantMasses = vars.vertices_masses;

            // In general, the variables are initialized, the function is then executed and the updated struct is returned and the variables then assume their true values

            //---------------------------------------------------------------------------------------------------------
            // Jet classification block (based on constituents info)
            //---------------------------------------------------------------------------------------------------------

            JetTagResult tagResult = classifyJet(jet, tagged_c_jets, tagged_s_jets, vec_c, vec_s);

            Bool_t isSignalTagged;
            Bool_t isCharmTagged = tagResult.isCharmTagged;
            Bool_t isStrangeTagged = tagResult.isStrangeTagged;

            if (switch_string == "strange")
            {
                isSignalTagged = isStrangeTagged;
            }
            else if (switch_string == "charm")
            {
                isSignalTagged = isCharmTagged;
            }
            
            if (!isSignalTagged) // Then it' a background jet
            {
                label = 0;
                eventID = ni;
                pT = jetPt;
                eta = jetEta;
                phi = jetPhi;
                nConst = jetNConst;
                nRho = third_nRho;
                
                jetVerticesInvariantMasses = vars.vertices_masses;
                background_pTDistribution->Fill(pT);
                background_nConstDistribution->Fill(nConst);
                first_background_nRho_Distribution->Fill(first_nRho);
                second_background_nRho_Distribution->Fill(second_nRho);
                third_background_nRho_Distribution->Fill(third_nRho);

                backgroundTree->Fill();
            }

        } // End of individual jet creation

        std::unordered_set<int> pdgSet;
        std::vector<fastjet::PseudoJet> tagged_signal_jets;
        Float_t signal_quark_pT = 0;

        if (switch_string == "strange") 
        {
            pdgSet = {130, 310, 311, 321, 313, 323, 315, 325, 317, 327, 319, 329, 9000311, 9000321, 10311, 10321, 100311, 100321, 9010311, 9010321, 9020311, 9020321, 313, 323, 10313, 10323, 20313, 20323, 100313, 100323, 9000313, 9000323, 30313, 30323, 315, 325, 9000315, 9000325, 10315, 10325, 20315, 20325, 9010315, 9010325, 9020315, 9020325, 317, 327, 9010317, 9010327, 319, 329, 3122, 3222, 3212, 3112, 3224, 3214, 3114, 3322, 3312, 3324, 3314, 3334};
            tagged_signal_jets = tagged_s_jets;
            signal_quark_pT = strangePt;
        }
        if (switch_string == "charm") 
        {
            pdgSet = {411, 421, 413, 423, 415, 425, 431, 433, 435, 10411, 10421, 413, 423, 10413, 10423, 20413, 20423, 415, 425, 431, 10431, 433, 10433, 20433, 435, 4122, 4222, 4212, 4112, 4224, 4214, 4114, 4232, 4132, 4322, 4312, 4324, 4314, 4332, 4334, 4412, 4422, 4414, 4424, 4432, 4434, 4444};
            tagged_signal_jets = tagged_c_jets;
            signal_quark_pT = charmPt;
        }

        //---------------------------------------------------------------------------------------------------------
        // Signal jets tagging
        //---------------------------------------------------------------------------------------------------------
        

        for (const fastjet::PseudoJet &jet : tagged_signal_jets)
        {
            TLorentzVector sJet(jet.px(), jet.py(), jet.pz(), jet.E());
            Float_t signal_pTRatio = jet.pt() / signal_quark_pT;

           Bool_t hasSignalConstituent = false;
           Int_t missingSignalPdg = 0;

            for (const fastjet::PseudoJet &constituent : jet.constituents())
            {
                Int_t constituentPdg                = constituent.user_info<JetInfo>().getFinalParticlePdg();
                Int_t constituentMotherPdg          = constituent.user_info<JetInfo>().getFinalParticleMotherPdg();
                Int_t constituentSecondMotherPdg    = constituent.user_info<JetInfo>().getFinalParticleSecondMotherPdg();
                Int_t constituentThirdMotherPdg     = constituent.user_info<JetInfo>().getFinalParticleThirdMotherPdg();

                Int_t abs_constituentPdg                = abs(constituentPdg);
                Int_t abs_constituentMotherPdg          = abs(constituentMotherPdg);
                Int_t abs_constituentSecondMotherPdg    = abs(constituentSecondMotherPdg);
                Int_t abs_constituentThirdMotherPdg     = abs(constituentThirdMotherPdg);

                if (pdgSet.count(abs_constituentPdg) || pdgSet.count(abs_constituentMotherPdg) || pdgSet.count(abs_constituentSecondMotherPdg) || pdgSet.count(abs_constituentThirdMotherPdg))
                {
                    hasSignalConstituent = true;
                    break;
                }
            }
            if (hasSignalConstituent) // Data that got through the signal filters
            {
                eventID = ni;

                //Kinematics directly from PseudoJet
                pT = jet.pt();
                eta = jet.eta();
                phi = jet.phi();
                nConst = (Int_t)jet.constituents().size();

                //nRho from JetUserInfo data

                first_nRho = 0;
                second_nRho = 0;
                third_nRho = 0;
                pT_LeadConst = 0;

                JetVariables vars = extractJetVariables(
                    jet,
                    rhoLowerBound,
                    firstRhoUpperBound,
                    secondRhoUpperBound,
                    thirdRhoUpperBound
                );

                first_nRho = vars.first_nRho;
                second_nRho = vars.second_nRho;
                third_nRho = vars.third_nRho;
                nRho = third_nRho;

                pT_LeadConst = vars.leadingPt;
                jetVerticesInvariantMasses = vars.vertices_masses;

                if (signal_pTRatio > 0.7)   // Indeed a signal jet
                {
                    label = 1;
                    
                    signal_pTDistribution->Fill(pT);
                    signal_nConstDistribution->Fill(nConst);
                    first_signal_nRho_Distribution->Fill(first_nRho);
                    second_signal_nRho_Distribution->Fill(second_nRho);
                    third_signal_nRho_Distribution->Fill(third_nRho);

                    signalTree->Fill();
                }
                else if (contaminatingGluonMode == "include")
                {
                    label = 0;

                    background_pTDistribution->Fill(pT);
                    background_nConstDistribution->Fill(nConst);
                    first_background_nRho_Distribution->Fill(first_nRho);
                    second_background_nRho_Distribution->Fill(second_nRho);
                    third_background_nRho_Distribution->Fill(third_nRho);

                    backgroundTree->Fill();
                }
                // else if (contaminatingGluonMode == "exclude") do nothing, i.e., don't fill the background tree with this jet -> No coding needed (mode justified only for the file name)
            }
            
        }
        
        particles_fastjet.clear();
        jets.clear();
        jets_array->Clear();
        quarks->Clear();

    
    } // End of event loop equivalent

    signalTree->Write();
    backgroundTree->Write();

    TCanvas *c1 = new TCanvas("c1", "Strange-jet classifier: N_{#rho} distributions", 2500, 2500);
    c1->Divide(3, 2);

    // Signal (upper bound = 1, 1.5, 2)
    c1->cd(1);
    first_signal_nRho_Distribution->SetTitle((uppercase_switch_string+" signal jets: N_{#rho} distribution (upper bound = 1)").c_str());
    first_signal_nRho_Distribution->GetXaxis()->SetTitle("N_{#rho}");
    first_signal_nRho_Distribution->GetYaxis()->SetTitle("Entries");
    first_signal_nRho_Distribution->DrawCopy();

    c1->cd(2);
    second_signal_nRho_Distribution->SetTitle((uppercase_switch_string+" signal jets: N_{#rho} distribution (upper bound = 1.5)").c_str());
    second_signal_nRho_Distribution->GetXaxis()->SetTitle("N_{#rho}");
    second_signal_nRho_Distribution->GetYaxis()->SetTitle("Entries");
    second_signal_nRho_Distribution->DrawCopy();

    c1->cd(3);
    third_signal_nRho_Distribution->SetTitle((uppercase_switch_string+" signal jets: N_{#rho} distribution (upper bound = 2)").c_str());
    third_signal_nRho_Distribution->GetXaxis()->SetTitle("N_{#rho}");
    third_signal_nRho_Distribution->GetYaxis()->SetTitle("Entries");
    third_signal_nRho_Distribution->DrawCopy();

    // Background (upper bound = 1, 1.5, 2)
    c1->cd(4);
    first_background_nRho_Distribution->SetTitle((uppercase_switch_string+" background jets: N_{#rho} distribution (upper bound = 1)").c_str());
    first_background_nRho_Distribution->GetXaxis()->SetTitle("N_{#rho}");
    first_background_nRho_Distribution->GetYaxis()->SetTitle("Entries");
    first_background_nRho_Distribution->DrawCopy();

    c1->cd(5);
    second_background_nRho_Distribution->SetTitle((uppercase_switch_string+" background jets: N_{#rho} distribution (upper bound = 1.5)").c_str());
    second_background_nRho_Distribution->GetXaxis()->SetTitle("N_{#rho}");
    second_background_nRho_Distribution->GetYaxis()->SetTitle("Entries");
    second_background_nRho_Distribution->DrawCopy();

    c1->cd(6);
    third_background_nRho_Distribution->SetTitle((uppercase_switch_string+" background jets: N_{#rho} distribution (upper bound = 2)").c_str());
    third_background_nRho_Distribution->GetXaxis()->SetTitle("N_{*rho}");
    third_background_nRho_Distribution->GetYaxis()->SetTitle("Entries");
    third_background_nRho_Distribution->DrawCopy();


    TCanvas *c2 = new TCanvas("c2", "Transverse momentum (signal vs background)", 2500, 2500);
    c2->Divide(2, 1);

    c2->cd(1);
    signal_pTDistribution->SetTitle((uppercase_switch_string+" signal jets: p_{T} distribution").c_str());
    signal_pTDistribution->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    signal_pTDistribution->GetYaxis()->SetTitle("Entries");
    signal_pTDistribution->DrawCopy();

    c2->cd(2);
    background_pTDistribution->SetTitle((uppercase_switch_string+" background jets: p_{T} distribution").c_str());
    background_pTDistribution->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    background_pTDistribution->GetYaxis()->SetTitle("Entries");
    background_pTDistribution->DrawCopy();


    TCanvas *c3 = new TCanvas("c3", "Number of constituents (signal vs background)", 2500, 2500);
    c3->Divide(2, 1);

    c3->cd(1);
    signal_nConstDistribution->SetTitle((uppercase_switch_string+" signal jets: N_{const} distribution").c_str());
    signal_nConstDistribution->GetXaxis()->SetTitle("N_{const}");
    signal_nConstDistribution->GetYaxis()->SetTitle("Entries");
    signal_nConstDistribution->DrawCopy();

    c3->cd(2);
    background_nConstDistribution->SetTitle((uppercase_switch_string+" background jets: N_{const} distribution").c_str());
    background_nConstDistribution->GetXaxis()->SetTitle("N_{const}");
    background_nConstDistribution->GetYaxis()->SetTitle("Entries");
    background_nConstDistribution->DrawCopy();

    file->Close();
    filteredDataFile->Close();
}

