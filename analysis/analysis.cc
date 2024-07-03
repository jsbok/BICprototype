#include "RootInterface.h"
#include "DRsimInterface.h"
#include "functions.h"

#include "TROOT.h"
#include "TStyle.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TPaveStats.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "TGraph.h"

#include <iostream>
#include <string>

int main(int argc, char* argv[]) {
  TString filename = argv[1];
  float low = std::stof(argv[2]);
  float high = std::stof(argv[3]);

  gStyle->SetOptFit(1);

  if (!filename.EndsWith(".root"))
    filename = filename + ".root";
  RootInterface<DRsimInterface::DRsimEventData>* drInterface = new RootInterface<DRsimInterface::DRsimEventData>(std::string(filename), true);
  drInterface->set("DRsim","DRsimEventData");
  //drInterface->GetChain("DRsim");

  TH1F* tEdep = new TH1F("totEdep",";MeV;Evt",100,low*1000.,high*1000.);
  tEdep->Sumw2(); tEdep->SetLineColor(kRed); tEdep->SetLineWidth(2);
  TH1F* tHit_S = new TH1F("Hit_S",";# of p.e.;Evt",200,12000*low*2,12000*high*2);
  tHit_S->Sumw2(); tHit_S->SetLineColor(kRed); tHit_S->SetLineWidth(2);
  TH1F* tHit_S_left = new TH1F("Hit_S_left",";# of p.e.;Evt",200,12000*low,12000*high);
  tHit_S_left->Sumw2(); tHit_S_left->SetLineColor(kRed); tHit_S_left->SetLineWidth(2);
  TH1F* tHit_S_right = new TH1F("Hit_S_right",";# of p.e.;Evt",200,12000*low,12000*high);
  tHit_S_right->Sumw2(); tHit_S_right->SetLineColor(kRed); tHit_S_right->SetLineWidth(2);
  TH1F* tP_leak = new TH1F("Pleak",";MeV;Evt",100,0.,1000.*high);
  tP_leak->Sumw2(); tP_leak->SetLineWidth(2);
  TH1F* tP_leak_nu = new TH1F("Pleak_nu",";MeV;Evt",100,0.,1000.*high);
  tP_leak_nu->Sumw2(); tP_leak_nu->SetLineWidth(2);

  TH1F* tT_S = new TH1F("time_S",";ns;p.e.",700,0.,70.);
  tT_S->Sumw2(); tT_S->SetLineColor(kRed); tT_S->SetLineWidth(2);
  TH1F* tT_S_left = new TH1F("time_S_left",";ns;p.e.",700,0.,70.);
  tT_S_left->Sumw2(); tT_S_left->SetLineColor(kRed); tT_S_left->SetLineWidth(2);
  TH1F* tT_S_right = new TH1F("time_S_right",";ns;p.e.",700,0.,70.);
  tT_S_right->Sumw2(); tT_S_right->SetLineColor(kRed); tT_S_right->SetLineWidth(2);
  TH1F* tWav_S = new TH1F("wavlen_S",";nm;p.e.",120,300.,900.);
  tWav_S->Sumw2(); tWav_S->SetLineColor(kRed); tWav_S->SetLineWidth(2);
  TH1F* tWav_S_left = new TH1F("wavlen_S_left",";nm;p.e.",120,300.,900.);
  tWav_S_left->Sumw2(); tWav_S_left->SetLineColor(kRed); tWav_S_left->SetLineWidth(2);
  TH1F* tWav_S_right = new TH1F("wavlen_S_right",";nm;p.e.",120,300.,900.);
  tWav_S_right->Sumw2(); tWav_S_right->SetLineColor(kRed); tWav_S_right->SetLineWidth(2);
  TH1F* tNhit_S = new TH1F("nHits_S",";p.e.;n",200,0.,200.);
  tNhit_S->Sumw2(); tNhit_S->SetLineColor(kRed); tNhit_S->SetLineWidth(2);

  TH2D* t2DhitS = new TH2D("2D Hit S", "", 420, -0.5, 419.5, 420, -0.5, 419.5); t2DhitS->Sumw2(); t2DhitS->SetStats(0);

  unsigned int entries = drInterface->entries();
  while (drInterface->numEvt() < entries) {
    if (drInterface->numEvt() % 100 == 0) printf("Analyzing %dth event ...\n", drInterface->numEvt());

    DRsimInterface::DRsimEventData drEvt;
    drInterface->read(drEvt);

    float Edep = 0.;
    for (auto edepItr = drEvt.Edeps.begin(); edepItr != drEvt.Edeps.end(); ++edepItr) {
      auto edep = *edepItr;
      Edep += edep.Edep;
    }
    tEdep->Fill(Edep);

    float Pleak = 0.;
    float Eleak_nu = 0.;
    for (auto leak : drEvt.leaks) {
      TLorentzVector leak4vec;
      leak4vec.SetPxPyPzE(leak.px,leak.py,leak.pz,leak.E);
      if ( std::abs(leak.pdgId)==12 || std::abs(leak.pdgId)==14 || std::abs(leak.pdgId)==16 ) {
        Eleak_nu += leak4vec.P();
      } else {
        Pleak += leak4vec.P();
      }
    }
    tP_leak->Fill(Pleak);
    tP_leak_nu->Fill(Eleak_nu);

    int nHitS = 0; int isLeft = 0;
    int nHitS_left = 0; int nHitS_right = 0;
    for (auto tower = drEvt.towers.begin(); tower != drEvt.towers.end(); ++tower) {
      int moduleNum = tower->ModuleNum;
      for (auto sipm = tower->SiPMs.begin(); sipm != tower->SiPMs.end(); ++sipm) {
        isLeft = sipm->isleft;
        int plateNum = sipm->x; int fiberNum = sipm->y; 

        tNhit_S->Fill(sipm->count);
        nHitS += sipm->count;
        
        if (isLeft==0) {nHitS_left += sipm->count;}
        else {nHitS_right += sipm->count;}

        t2DhitS->Fill(60*(moduleNum%7)+fiberNum, 60*(moduleNum/7)+plateNum, sipm->count);
        for (const auto timepair : sipm->timeStruct) {
          tT_S->Fill(timepair.first.first+0.05,timepair.second);
          if (isLeft==0) {tT_S_left->Fill(timepair.first.first+0.05,timepair.second);}
          else {tT_S_right->Fill(timepair.first.first+0.05,timepair.second);}
        }
        for (const auto wavpair : sipm->wavlenSpectrum) {
          tWav_S->Fill(wavpair.first.first,wavpair.second);
          if (isLeft==0) {tWav_S_left->Fill(wavpair.first.first,wavpair.second);}
          else {tWav_S_right->Fill(wavpair.first.first,wavpair.second);}
        }
      }
    }

    tHit_S->Fill(nHitS);
    tHit_S_left->Fill(nHitS_left);
    tHit_S_right->Fill(nHitS_right);
  } // event loop
  //drInterface->close();

  TCanvas* c = new TCanvas("c","");

  tEdep->Draw("Hist"); c->SaveAs(filename+"_Edep.png");

  c->SetLogy(1);
  tP_leak->Draw("Hist"); c->SaveAs(filename+"_Pleak.png");
  tP_leak_nu->Draw("Hist"); c->SaveAs(filename+"_Pleak_nu.png");
  c->SetLogy(0);

  tHit_S->Draw("Hist"); c->SaveAs(filename+"_nHitpEventS.pdf");
  tHit_S_left->Draw("Hist"); c->SaveAs(filename+"_nHitpEventS_left.pdf");
  tHit_S_right->Draw("Hist"); c->SaveAs(filename+"_nHitpEventS_right.pdf");

  t2DhitS->Draw("COLZ"); c->SaveAs(filename+"_n2DHitS.png");

  tT_S->Draw("Hist"); c->SaveAs(filename+"_tS.png");
  tT_S_left->Draw("Hist"); c->SaveAs(filename+"_tS_left.png");
  tT_S_right->Draw("Hist"); c->SaveAs(filename+"_tS_right.png");
  tWav_S->Draw("Hist"); c->SaveAs(filename+"_wavS.png");
  tWav_S_left->Draw("Hist"); c->SaveAs(filename+"_wavS_left.png");
  tWav_S_right->Draw("Hist"); c->SaveAs(filename+"_wavS_right.png");
  tNhit_S->Draw("Hist"); c->SaveAs(filename+"_nhitS.png");
}
