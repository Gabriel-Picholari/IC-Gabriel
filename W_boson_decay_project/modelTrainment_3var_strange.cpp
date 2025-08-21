#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TMVA/Reader.h>
#include <TMVA/Factory.h>
#include <TMVA/DataLoader.h>

void modelTrainment_3var_strange(const char *fileName) 
{
    
    TFile* outputFile = TFile::Open("TMVAOutput_3var_strange.root", "RECREATE");
    TFile* inputFile = TFile::Open(fileName, "READ");

    TTree* signalTree_s = dynamic_cast<TTree *>(inputFile->Get("SignalTree_s"));
    TTree* backgroundTree_s = dynamic_cast<TTree *>(inputFile->Get("BackgroundTree_s"));

    TMVA::Factory factory("TMVAClassification", outputFile, "AnalysisType=Classification");
    TMVA::DataLoader loader("dataset_s_3var");

    // Variáveis de entrada
    loader.AddVariable("pT_s", 'F');
    loader.AddVariable("nRho_s", 'I');
    loader.AddVariable("nConst_s", 'F');
    //loader.AddVariable("maxRho_s", 'F'); Discontinued

    // Variáveis espectadoras
    loader.AddSpectator("eta_s", "F");
    loader.AddSpectator("phi_s", "F");
    loader.AddSpectator("mass_s", "F");
    loader.AddSpectator("label_s", "F");
    loader.AddSpectator("eventID_s", "F");

    // Trees de sinal e fundo
    loader.AddSignalTree(signalTree_s, 1.0);
    loader.AddBackgroundTree(backgroundTree_s, 1.0);

    //---------------------------------------------------------------------------------------------------------
    // Divisão manual das Trees de sinal e fundo em treinamento e teste
    //---------------------------------------------------------------------------------------------------------

    TCut mycut = "";

    Int_t nTest_Signal = signalTree_s->GetEntries() * 0.2;
    Int_t nTrain_Signal = signalTree_s->GetEntries() * 0.8;

    Int_t nTest_Background = backgroundTree_s->GetEntries() * 0.2;
    Int_t nTrain_Background = backgroundTree_s->GetEntries() * 0.8;

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
