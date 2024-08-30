#include "TSystem.h"
#include "TH1F.h"
#include "TClonesArray.h"
#include "TPythia8.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TClonesArray.h"
#include "MyPart.h"
#include "TFile.h"

void read_myTree(const char* fileName)
{

  Float_t mass, eta = 0, light_speed = 3e8;
  Float_t px, py, pz = 0;
  Float_t eX, eY, eZ, eTot, eT = 0;
  Float_t sum_eX, sum_eY, sum_eZ, sum_eTot, sum_eT = 0;

  TH1F *soma_energia_X = new TH1F("histograma1", "Soma das energias das particulas em estado final em X (s/ restricao)", 100, 0, 0);
  TH1F *soma_energia_Y = new TH1F("histograma2", "Soma das energias das particulas em estado final em Y (s/ restricao)",  100, 0, 0);
  TH1F *soma_energia_Z = new TH1F("histograma3", "Soma das energias das particulas em estado final em Z (s/ restricao)",  100, 0, 0);
  TH1F *energias_transversais = new TH1F("histograma4", "Energia total transversal das particulas finais", 100, 0, 0);
  TH1F *total_energy = new TH1F("histograma5", "Energia total obtida a partir das componentes (s/restricao)", 100, 0, 0);

    TFile *file = TFile::Open(fileName, "READ");

    if (!file || !file->IsOpen())
    {
        std::cerr << "Error: Could not open the file." << std::endl;
        return;
    }

    TTree *tree = dynamic_cast<TTree *>(file->Get("partTree"));

    if (!tree)
    {
        std::cerr << "Error: Could not find the TTree in the file." << std::endl;
        file->Close();
        return;
    }

    TClonesArray *fparticles = new TClonesArray("MyPart");
    tree->SetBranchAddress("particles", &fparticles);

    Long64_t ne = tree->GetEntries();

    // Correspondente ao loop de eventos
    for (Long64_t ni = 0; ni < ne; ++ni)
    {
        tree->GetEntry(ni);

        sum_eX = 0;
        sum_eY = 0;
        sum_eZ = 0;
        sum_eTot = 0;
        sum_eT = 0;

        // Correspondente ao loop de particulas
        for (Int_t nj = 0; nj < fparticles->GetEntries(); nj++)
        {
            MyPart *part = static_cast<MyPart *>(fparticles->At(nj));

            px = part->fPx;
            py = part->fPy;
            pz = part->fPz;
            mass = part->fMass;
            eta = part->fEta;
            eTot = part->fE;

            eT = eTot / TMath::CosH(eta);
            eX = TMath::Sqrt( TMath::Power( px, 2 ) + TMath::Power( mass, 2 ) );
            eY = TMath::Sqrt( TMath::Power( py, 2 ) + TMath::Power( mass, 2 ) );
            eZ = TMath::Sqrt( TMath::Power( pz, 2 ) + TMath::Power( mass, 2 ) );

            sum_eX += eX;
            sum_eY += eY;
            sum_eZ += eZ;
            sum_eTot += eTot;
            sum_eT += eT;
        }
        soma_energia_X->Fill(sum_eX);
        soma_energia_Y->Fill(sum_eY);
        soma_energia_Z->Fill(sum_eZ);
        total_energy->Fill(sum_eTot);
        energias_transversais->Fill(sum_eT);

        fparticles->Clear();
    }

    TCanvas *c1 = new TCanvas("c1", "Histogramas e distribuicao", 2500, 2500);
    c1->Divide(2, 3);

    c1->cd(1);
    total_energy->SetTitle("Soma das energias das particulas de estado final");
    total_energy->GetXaxis()->SetTitle("Soma por evento [GeV]");
    total_energy->GetYaxis()->SetTitle("Frequencia");
    total_energy->Draw();

    c1->cd(2);
    soma_energia_X->SetTitle("Soma das energias das particulas em X");
    soma_energia_X->GetXaxis()->SetTitle("Soma por evento [GeV]");
    soma_energia_X->GetYaxis()->SetTitle("Frequencia");
    soma_energia_X->Draw();

    c1->cd(3);
    soma_energia_Y->SetTitle("Soma das energias das particulas em Y");
    soma_energia_Y->GetXaxis()->SetTitle("Soma por evento [GeV]");
    soma_energia_Y->GetYaxis()->SetTitle("Frequencia");
    soma_energia_Y->Draw();

    c1->cd(4);
    soma_energia_Z->SetTitle("Soma das energias das particulas em Z");
    soma_energia_Z->GetXaxis()->SetTitle("Soma por evento [GeV]");
    soma_energia_Z->GetYaxis()->SetTitle("Frequencia");
    soma_energia_Z->Draw();

    c1->cd(5);
    energias_transversais->SetTitle("Soma das energias transversais das particulas em estado final");
    energias_transversais->GetXaxis()->SetTitle("Soma por evento [GeV]");
    energias_transversais->GetYaxis()->SetTitle("Frequencia");
    energias_transversais->Draw();

    file->Close();
}