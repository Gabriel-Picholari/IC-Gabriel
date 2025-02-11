#include "TSystem.h"
#include "TClonesArray.h"
#include "TDatabasePDG.h"
#include "TParticle.h"
#include "TPythia8.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"

void soma_de_energias(Int_t nev = 10000, Int_t ndeb = 1 /* Listagem */ )
{
  gSystem->Load("libEG");
  gSystem->Load("libEgPythia8");

  Double_t light_speed = 3e8;
  Double_t mass, eta;
  Double_t px, py, pz;
  Double_t energy_x, energy_y, energy_z, energy_tot, energy_t;
  Double_t sum_Ex, sum_Ey, sum_Ez, energy_sum_testing, energy_sum_t;

  TH1F *soma_energia_X = new TH1F("histograma1", "Soma das energias das particulas em estado final em X (s/ restricao)", 100, 0, 0);
  TH1F *soma_energia_Y = new TH1F("histograma2", "Soma das energias das particulas em estado final em Y (s/ restricao)",  100, 0, 0);
  TH1F *soma_energia_Z = new TH1F("histograma3", "Soma das energias das particulas em estado final em Z (s/ restricao)",  100, 0, 0);
  TH1F *energias_transversais = new TH1F("histograma4", "Energia total transversal das particulas finais", 100, 0, 0);
  TH1F *total_energy = new TH1F("histograma5", "Energia total obtida a partir das componentes (s/restricao)", 100, 0, 0);

  TPythia8 *pythia8 = new TPythia8();
  pythia8->ReadString("HardQCD:all = on");
  pythia8->ReadString("Random:setSeed = on");
  pythia8->ReadString("Random:seed = 42");

  pythia8->Initialize(2212 /* Proton */, 2212 /* Proton */, 14000 /* TeV */); /* 14000 TeV = 14000000 GeV */

  TClonesArray *particles = new TClonesArray("TParticle", 10000);

  for (Int_t iev = 0; iev < nev; iev++)
  {
    sum_Ex = 0;
    sum_Ey = 0;
    sum_Ez = 0;
    energy_sum_testing = 0;
    energy_sum_t = 0;
    pythia8->GenerateEvent();

    /* if (iev < ndeb) pythia8->EventListing(); */

    pythia8->ImportParticles(particles, "All");
    Int_t np = particles->GetEntriesFast();

    for (Int_t ip = 0; ip < np; ip++)
    {
      TParticle *part = (TParticle*) particles->At(ip); /* Criacao do ponteiro */

      Int_t status = part->GetStatusCode();
      Double_t pdg = part->GetPdgCode();

      if (status > 0) /* Particula de estado final */
      {
        px = part->Px();
        py = part->Py();
        pz = part->Pz();
        mass = part->GetMass();
        eta = part->Eta();

        energy_tot = part->Energy();
        energy_t = energy_tot / TMath::CosH(eta);

        energy_x = TMath::Sqrt( TMath::Power( px, 2 ) + TMath::Power( mass, 2 ) );
        energy_y = TMath::Sqrt( TMath::Power( py, 2 ) + TMath::Power( mass, 2 ) );
        energy_z = TMath::Sqrt( TMath::Power( pz, 2 ) + TMath::Power( mass, 2 ) );

        sum_Ex += energy_x;
        sum_Ey += energy_y;
        sum_Ez += energy_z;
        energy_sum_testing += energy_tot;
        energy_sum_t += energy_t;

      }
    }

    soma_energia_X->Fill(sum_Ex);
    soma_energia_Y->Fill(sum_Ey);
    soma_energia_Z->Fill(sum_Ez);
    total_energy->Fill(energy_sum_testing);
    energias_transversais->Fill(energy_sum_t);

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
}