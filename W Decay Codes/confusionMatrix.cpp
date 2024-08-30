//---------------------------------------------------------------------------------------------------------//
//                                                                                                         //
// This code retrieves the trees of testing and verifies their compatibility with the training trees       //
// in order to determine true and false positives and negatives. Then, a confusion matrix is built.        //
// The code utilizes a loop counter strategy. Important to mention that the file read by this macro is     //
// the one that creates both the signal and bgr Trees and NOT the modelTraining.cpp output.                //
//                                                                                                         //
//---------------------------------------------------------------------------------------------------------//

#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <TROOT.h>
#include <TMVA/Reader.h>
#include <TMVA/Results.h>
#include <TMVA/DataLoader.h>

void confusionMatrix(const char *fileName)
{
    //---------------------------------------------------------------------------------------------------------
    // Recuperacao da Tree de teste
    //---------------------------------------------------------------------------------------------------------

    TFile* inputFile = TFile::Open(fileName, "READ");
    TTree* testTree = dynamic_cast<TTree*>(inputFile->Get("dataset/TestTree"));

    //---------------------------------------------------------------------------------------------------------
    // Inicializacao do objeto reader recuperando o arquivo xml gerado pela macro modelTrainment.cpp
    //---------------------------------------------------------------------------------------------------------

    TMVA::Reader reader;

    Float_t pT_c, nConst_c, pT_lConst_c, label_c;

    reader.AddVariable("pT_c", &pT_c);
    reader.AddVariable("nConst_c", &nConst_c);
    reader.AddVariable("pT_lConst_c", &pT_lConst_c);

    reader.AddSpectator("label_c", &label_c);

    reader.BookMVA("LogisticRegression", "dataset/weights/TMVARegression_LogisticRegression.weights.xml");

    //---------------------------------------------------------------------------------------------------------
    // Loop que contabiliza os verdadeiros positivos (TP), verdadeiros negativos (TN), falsos positivos (TP) e
    // falsos negativos (FN) usando o objeto reader
    //---------------------------------------------------------------------------------------------------------

    int TP = 0, TN = 0, FP = 0, FN = 0;

    testTree->SetBranchAddress("pT_c", &pT_c);
    testTree->SetBranchAddress("label_c", &label_c);
    testTree->SetBranchAddress("nConst_c", &nConst_c);
    testTree->SetBranchAddress("pT_lConst_c", &pT_lConst_c);

    Long64_t nEntries = testTree->GetEntries();

    for (Long64_t i = 0; i < nEntries; ++i) 
    {
        testTree->GetEntry(i);
        double response = reader.EvaluateMVA("LogisticRegression");

        int predictedLabel;
        if ( response >= 0.5 )
        {
            predictedLabel = 1;
        }
        else
        {
            predictedLabel = 0;
        }

        if (predictedLabel == 1 && label_c == 1) 
        {
            TP++;
        } 
        else if (predictedLabel == 0 && label_c == 0) 
        {
            TN++;
        } 
        else if (predictedLabel == 1 && label_c == 0) 
        {
            FP++;
        } 
        else if (predictedLabel == 0 && label_c == 1) 
        {
            FN++;
        }
        
    }

    //---------------------------------------------------------------------------------------------------------
    // Matriz de confusao
    //---------------------------------------------------------------------------------------------------------
    
    std::cout << "Confusion Matrix:" << std::endl;
    std::cout << "TP: " << TP << "  FP: " << FP << std::endl;
    std::cout << "FN: " << FN << "  TN: " << TN << std::endl;

    inputFile->Close();
}