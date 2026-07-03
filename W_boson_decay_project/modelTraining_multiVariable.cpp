#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TMVA/Reader.h>
#include <TMVA/Factory.h>
#include <TMVA/DataLoader.h>

void modelTraining_multiVariable(const char *fileName, std::string switch_string, std::string contaminatingGluonMode) 
{

    std::string sfx;
    std::string datasetName;

    if (switch_string == "strange") sfx = "_s";
    else if (switch_string == "charm") sfx = "_c";

    TFile *outputFile = nullptr;
    if (switch_string == "charm")  
    {
        if (contaminatingGluonMode == "include")
        {
            outputFile = TFile::Open("TMVAOutput_multiVariable_gluonJetsIncluded_charm.root", "RECREATE");
            datasetName = "dataset_c_multiVariable_gluonJetsIncluded";
        }
        if (contaminatingGluonMode == "exclude")
        {
            outputFile = TFile::Open("TMVAOutput_multiVariable_gluonJetsExcluded_charm.root", "RECREATE");
            datasetName = "dataset_c_multiVariable_gluonJetsExcluded";
        }
    }
    if (switch_string == "strange")
    {
        if (contaminatingGluonMode == "include")
        {
            outputFile = TFile::Open("TMVAOutput_multiVariable_gluonJetsIncluded_strange.root", "RECREATE");
            datasetName = "dataset_s_multiVariable_gluonJetsIncluded";
        }
        if (contaminatingGluonMode == "exclude")
        {
            outputFile = TFile::Open("TMVAOutput_multiVariable_gluonJetsExcluded_strange.root", "RECREATE");
            datasetName = "dataset_s_multiVariable_gluonJetsExcluded";
        }
    }


    TFile* inputFile = TFile::Open(fileName, "READ");

    TTree* signalTree = dynamic_cast<TTree *>(inputFile->Get(("SignalTree" + sfx).c_str()));
    TTree* backgroundTree = dynamic_cast<TTree *>(inputFile->Get(("BackgroundTree" + sfx).c_str()));

    TMVA::Factory factory("TMVAClassification", outputFile, "AnalysisType=Classification");
    TMVA::DataLoader loader(datasetName.c_str());

    // Variáveis de entrada
    loader.AddVariable("pT" + sfx, 'F');
    loader.AddVariable("nRho" + sfx, 'I');
    loader.AddVariable("nConst" + sfx, 'F');

    // Variáveis espectadoras
    loader.AddSpectator("eta" + sfx, "F");
    loader.AddSpectator("phi" + sfx, "F");
    loader.AddSpectator("label" + sfx, "F");
    loader.AddSpectator("eventID" + sfx, "F");
    loader.AddSpectator("flavor" + sfx, "I");

    // Trees de sinal e fundo
    loader.AddSignalTree(signalTree, 1.0);
    loader.AddBackgroundTree(backgroundTree, 1.0);

    //---------------------------------------------------------------------------------------------------------
    // Divisão manual das Trees de sinal e fundo em treinamento e teste
    //---------------------------------------------------------------------------------------------------------

    TCut mycut = "";

    Int_t nTest_Signal = signalTree->GetEntries() * 0.2;
    Int_t nTrain_Signal = signalTree->GetEntries() * 0.8;

    Int_t nTest_Background = backgroundTree->GetEntries() * 0.2;
    Int_t nTrain_Background = backgroundTree->GetEntries() * 0.8;

    TString options = TString::Format("SplitMode=Random:SplitSeed=0:NormMode=None:nTrain_Signal=%d:nTrain_Background=%d:nTest_Signal=%d:nTest_Background=%d:", nTrain_Signal, nTrain_Background, nTest_Signal, nTest_Background);

    loader.PrepareTrainingAndTestTree(mycut, options.Data());

    //---------------------------------------------------------------------------------------------------------
    // Treinamento com GradBoost (mesma configuração do charm)
    //---------------------------------------------------------------------------------------------------------

    factory.BookMethod(&loader, TMVA::Types::kBDT, "GradBoost", "BoostType=Grad:Shrinkage=0.1:MaxDepth=3:nCuts=150");

    factory.TrainAllMethods();
    factory.TestAllMethods();
    factory.EvaluateAllMethods();

    inputFile->Close();
    outputFile->Close();    
}
