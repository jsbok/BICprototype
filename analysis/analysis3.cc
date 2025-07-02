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
  TString filename2 = "/u/user/changhui/work/BICprototype/rootfiles/3by5_e-_3_GeV_final/root/R3by5_e-_3_GeV_final_10M" ;
  TString filename1 = "/u/user/changhui/work/BICprototype/rootfiles/3by5_e-_3_GeV_z4_final/root/R3by5_e-_3_GeV_z4_final_10M" ;
  TString filename3 = "/u/user/changhui/work/BICprototype/rootfiles/3by5_e-_3_GeV_z8_final/root/R3by5_e-_3_GeV_z8_final_10M" ;
  TString filename4 = "/u/user/changhui/work/BICprototype/rootfiles/3by5_e-_3_GeV_z12_final/root/R3by5_e-_3_GeV_z12_final_10M" ;
  TString filename5 = "/u/user/changhui/work/BICprototype/rootfiles/3by5_e-_1_GeV_final/root/R4by4_e-_1_GeV_final_10M" ;
  TString filename6 = "/u/user/changhui/work/BICprototype/rootfiles/4by4_e-_500_MeV_final/root/R4by4_e-_500_MeV_final_10M" ;


  float inE = std::stof(argv[1]);
  float low = std::stof(argv[2]);
  float high = std::stof(argv[3]);

  gStyle->SetOptFit(1);

  RootInterface<DRsimInterface::DRsimEventData>* drInterface = new RootInterface<DRsimInterface::DRsimEventData>(std::string(filename1 + ".root"), true);
  drInterface->set("DRsim","DRsimEventData");
  //drInterface->GetChain("DRsim");

  TH1F* tEdep = new TH1F("totEdep",";MeV;Evt",100,low*1000.,high*1000.);
  tEdep->Sumw2(); tEdep->SetLineColor(kRed); tEdep->SetLineWidth(2);
  TH1F* tot_E = new TH1F("tot_E",";MeV;Evt",100,high*400.,high*600.);
  tot_E->Sumw2(); tot_E->SetLineColor(kRed); tot_E->SetLineWidth(2);
  TH1D* ratio_E = new TH1D("ratio_E",";evtN;ratio",100.,0,100.);
  ratio_E->Sumw2(); ratio_E->SetLineColor(kRed);

  TH1F* tEdep2 = new TH1F("totEdep2",";MeV;Evt",100,low*1000.,high*1000.);
  tEdep2->Sumw2(); tEdep2->SetLineColor(kBlue); tEdep2->SetLineWidth(2);
  TH1F* tot_E2 = new TH1F("tot_E2",";MeV;Evt",100,high*400.,high*600.);
  tot_E2->Sumw2(); tot_E2->SetLineColor(kBlue); tot_E2->SetLineWidth(2);
  TH1D* ratio_E2 = new TH1D("ratio_E",";evtN;ratio",100.,0,100.);
  ratio_E2->Sumw2(); ratio_E2->SetLineColor(kBlue);

  TH1F* tEdep3 = new TH1F("totEdep3",";MeV;Evt",100,low*1000.,high*1000.);
  tEdep3->Sumw2(); tEdep3->SetLineColor(kGreen); tEdep3->SetLineWidth(2);
  TH1F* tot_E3 = new TH1F("tot_E3",";MeV;Evt",100,high*400.,high*600.);
  tot_E3->Sumw2(); tot_E3->SetLineColor(kGreen); tot_E3->SetLineWidth(2);
  TH1D* ratio_E3 = new TH1D("ratio_E3",";evtN;ratio",100.,0,100.);
  ratio_E3->Sumw2(); ratio_E3->SetLineColor(kGreen);

  TH1F* tHit_S = new TH1F("Hit_S",";# of p.e.;Evt",200,12000*low*2,12000*high*2);
  tHit_S->Sumw2(); tHit_S->SetLineColor(kRed); tHit_S->SetLineWidth(2);
  TH1F* tHit_S2 = new TH1F("Hit_S2",";# of p.e.;Evt",200,12000*low*2,12000*high*2);
  tHit_S2->Sumw2(); tHit_S2->SetLineColor(kBlue); tHit_S2->SetLineWidth(2);
  TH1F* tHit_S3 = new TH1F("Hit_S3",";# of p.e.;Evt",200,12000*low*2,12000*high*2);
  tHit_S3->Sumw2(); tHit_S3->SetLineColor(kGreen); tHit_S3->SetLineWidth(2);


  TH1F* tHit_S_left = new TH1F("Hit_S_left",";# of p.e.;Evt",200,12000*low,12000*high);
  tHit_S_left->Sumw2(); tHit_S_left->SetLineColor(kRed); tHit_S_left->SetLineWidth(2);
  TH1F* tHit_S_left2 = new TH1F("Hit_S_left2",";# of p.e.;Evt",200,12000*low,12000*high);
  tHit_S_left2->Sumw2(); tHit_S_left2->SetLineColor(kBlue); tHit_S_left2->SetLineWidth(2);
  TH1F* tHit_S_left3 = new TH1F("Hit_S_left3",";# of p.e.;Evt",200,12000*low,12000*high);
  tHit_S_left3->Sumw2(); tHit_S_left3->SetLineColor(kGreen); tHit_S_left3->SetLineWidth(2);
  TH1F* tHit_S_right = new TH1F("Hit_S_right",";# of p.e.;Evt",200,12000*low,12000*high);
  tHit_S_right->Sumw2(); tHit_S_right->SetLineColor(kRed); tHit_S_right->SetLineWidth(2);
 TH1F* tHit_S_right2 = new TH1F("Hit_S_right2",";# of p.e.;Evt",200,12000*low,12000*high);
  tHit_S_right2->Sumw2(); tHit_S_right2->SetLineColor(kBlue); tHit_S_right2->SetLineWidth(2);
 TH1F* tHit_S_right3 = new TH1F("Hit_S_right3",";# of p.e.;Evt",200,12000*low,12000*high);
  tHit_S_right3->Sumw2(); tHit_S_right3->SetLineColor(kGreen); tHit_S_right3->SetLineWidth(2);

 TH1F* THits = new TH1F("THits",";ModuleN; # of p.e",20,0,20);
 THits->Sumw2(); THits->SetLineColor(kRed); THits->SetLineWidth(2);
 TH1F* THits2 = new TH1F("THits2",";ModuleN; # of p.e.",20,0,20);
 THits2->Sumw2(); THits2->SetLineColor(kBlue); THits2->SetLineWidth(2);
 TH1F* THits3 = new TH1F("THits3",";ModuleN; # of p.e.",20,0,20);
 THits3->Sumw2(); THits3->SetLineColor(kGreen); THits3->SetLineWidth(2);

  TH1F* tP_leak = new TH1F("Pleak",";MeV;Evt",100,0.,1000.*high);
  tP_leak->Sumw2(); tP_leak->SetLineWidth(2); tP_leak->SetLineColor(kRed);
  TH1F* tP_leak_nu = new TH1F("Pleak_nu",";MeV;Evt",100,0.,100.*high);
  tP_leak_nu->Sumw2(); tP_leak_nu->SetLineWidth(2); tP_leak_nu->SetLineColor(kRed);
  TH1F* tP_leak2 = new TH1F("Pleak2",";MeV;Evt",100,0.,1000.*high);
  tP_leak2->Sumw2(); tP_leak2->SetLineWidth(2); tP_leak2->SetLineColor(kBlue);
  TH1F* tP_leak_nu2 = new TH1F("Pleak_nu2",";MeV;Evt",100,0.,100.*high);
  tP_leak_nu2->Sumw2(); tP_leak_nu2->SetLineWidth(2); tP_leak_nu2->SetLineColor(kBlue);
  TH1F* tP_leak3 = new TH1F("Pleak3",";MeV;Evt",100,0.,1000.*high);
  tP_leak3->Sumw2(); tP_leak3->SetLineWidth(2); tP_leak3->SetLineColor(kGreen);
  TH1F* tP_leak_nu3 = new TH1F("Pleak_nu3",";MeV;Evt",100,0.,100.*high);
  tP_leak_nu3->Sumw2(); tP_leak_nu3->SetLineWidth(2);tP_leak_nu3->SetLineColor(kGreen);

  TH1F* leak_ratio = new TH1F("leak_ratio","; z0; leakage (%)",20,0.,20.);
  leak_ratio->Sumw2(); leak_ratio->SetLineWidth(2); leak_ratio->SetLineColor(kRed);
  TH1F* leak_ratio2 = new TH1F("leak_ratio","; z0; leakage (%)",20,0.,20.);
  leak_ratio2->Sumw2(); leak_ratio2->SetLineWidth(2); leak_ratio2->SetLineColor(kBlue);
  TH1F* leak_ratio3 = new TH1F("leak_ratio","; z0;leakage (%)",20,0.,20.);
  leak_ratio3->Sumw2(); leak_ratio3->SetLineWidth(2); leak_ratio3->SetLineColor(kGreen);

/*
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
*/  
TH1F* tNhit_S = new TH1F("nHits_S",";p.e.;n",200,0.,200.);
  tNhit_S->Sumw2(); tNhit_S->SetLineColor(kRed); tNhit_S->SetLineWidth(2);
TH1F* tNhit_S2 = new TH1F("nHits_S2",";p.e.;n",200,0.,200.);
  tNhit_S2->Sumw2(); tNhit_S2->SetLineColor(kBlue); tNhit_S2->SetLineWidth(2);
TH1F* tNhit_S3 = new TH1F("nHits_S3",";p.e.;n",200,0.,200.);
  tNhit_S3->Sumw2(); tNhit_S3->SetLineColor(kGreen); tNhit_S3->SetLineWidth(2);



  TH2D* t2DhitS = new TH2D("2D Hit S", "", 100, -0.5, 99.5, 100, -0.5, 99.5); t2DhitS->Sumw2(); t2DhitS->SetStats(0);
  TH2D* t2DhitS2 = new TH2D("2D Hit S2", "", 100, -0.5, 99.5, 100, -0.5, 99.5); t2DhitS2->Sumw2(); t2DhitS2->SetStats(0);
  TH2D* t2DhitS3 = new TH2D("2D Hit S3", "", 100, -0.5, 99.5, 100, -0.5, 99.5); t2DhitS3->Sumw2(); t2DhitS3->SetStats(0);

  unsigned int entries = drInterface->entries();
  while (drInterface->numEvt() < entries) {
    if (drInterface->numEvt() % 1000 == 0) printf("1st analyzing %dth event ...\n", drInterface->numEvt());

    DRsimInterface::DRsimEventData drEvt;
    drInterface->read(drEvt);

    float Edep = 0.; float totE = 0.; float ratE = 0;
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
    totE = Pleak + Eleak_nu + Edep;
    tot_E -> Fill(totE) ;
    ratE = totE / inE ;
    ratio_E -> Fill(drInterface->numEvt(),ratE);



    int nHitS = 0; int isLeft = 0; int nHits = 0;
    int nHitS_left = 0; int nHitS_right = 0;
    for (auto tower = drEvt.towers.begin(); tower != drEvt.towers.end(); ++tower) {
      int moduleNum = tower->ModuleNum;
      for (auto sipm = tower->SiPMs.begin(); sipm != tower->SiPMs.end(); ++sipm) {
        isLeft = sipm->isleft;
        int plateNum = sipm->x; int fiberNum = sipm->y; 

//std::cout <<" fiberNum : " << fiberNum << std::endl ;
//std::cout <<" plateNum : " << plateNum << std::endl ;

        tNhit_S->Fill(sipm->count);
        nHitS += sipm->count;
        nHits += sipm->count;   
        if (isLeft==0) {nHitS_left += sipm->count;}
        else {nHitS_right += sipm->count;}

        t2DhitS->Fill(22*(moduleNum%4)+fiberNum, 24*(moduleNum/4)+plateNum, sipm->count);
     /*   for (const auto timepair : sipm->timeStruct) {
          tT_S->Fill(timepair.first.first+0.05,timepair.second);
          if (isLeft==0) {tT_S_left->Fill(timepair.first.first+0.05,timepair.second);}
          else {tT_S_right->Fill(timepair.first.first+0.05,timepair.second);}
        }
        for (const auto wavpair : sipm->wavlenSpectrum) {
          tWav_S->Fill(wavpair.first.first,wavpair.second);
          if (isLeft==0) {tWav_S_left->Fill(wavpair.first.first,wavpair.second);}
          else {tWav_S_right->Fill(wavpair.first.first,wavpair.second);}
        }*/
THits -> Fill(moduleNum,nHits);
 
     }//SiPM loop
//THits -> Fill(moduleNum,nHits);
    }//tower loop

    tHit_S->Fill(nHitS);
    tHit_S_left->Fill(nHitS_left);
    tHit_S_right->Fill(nHitS_right);

  } // event loop
  //drInterface->close();

// second file
  RootInterface<DRsimInterface::DRsimEventData>* drInterface2 = new RootInterface<DRsimInterface::DRsimEventData>(std::string(filename2 + ".root"), true);
  drInterface2->set("DRsim","DRsimEventData"); 

  unsigned int entries2 = drInterface2->entries();
  while (drInterface2->numEvt() < entries) {
    if (drInterface2->numEvt() % 1000 == 0) printf("2nd analyzing %dth event ...\n", drInterface2->numEvt());

    DRsimInterface::DRsimEventData drEvt;
    drInterface2->read(drEvt);

    float Edep2 = 0.; float totE2 = 0.; float ratE2 = 0;
    for (auto edepItr2 = drEvt.Edeps.begin(); edepItr2 != drEvt.Edeps.end(); ++edepItr2) {
      auto edep2 = *edepItr2;
      Edep2 += edep2.Edep;
    }
    tEdep2->Fill(Edep2);

    float Pleak2 = 0.;
    float Eleak_nu2 = 0.;
    for (auto leak : drEvt.leaks) {
      TLorentzVector leak4vec;
      leak4vec.SetPxPyPzE(leak.px,leak.py,leak.pz,leak.E);
      if ( std::abs(leak.pdgId)==12 || std::abs(leak.pdgId)==14 || std::abs(leak.pdgId)==16 ) {
        Eleak_nu2 += leak4vec.P();
      } else {
        Pleak2 += leak4vec.P();
      }
    }
    tP_leak2->Fill(Pleak2);
    tP_leak_nu2->Fill(Eleak_nu2);
    totE2 = Pleak2 + Eleak_nu2 + Edep2;
    tot_E2 -> Fill(totE2) ;
    ratE2 = totE2 / inE ;
    ratio_E2 -> Fill(drInterface2->numEvt(),ratE2);

  int nHitS2 = 0; int isLeft2 = 0; int nHits2 = 0;
    int nHitS_left2 = 0; int nHitS_right2 = 0;
    for (auto tower2 = drEvt.towers.begin(); tower2 != drEvt.towers.end(); ++tower2) {
      int moduleNum2 = tower2->ModuleNum;
      for (auto sipm2 = tower2->SiPMs.begin(); sipm2 != tower2->SiPMs.end(); ++sipm2) {
        isLeft2 = sipm2->isleft;
        int plateNum2 = sipm2->x; int fiberNum2 = sipm2->y;
        tNhit_S2->Fill(sipm2->count);
        nHitS2 += sipm2->count;
	nHits2 += sipm2->count;
        if (isLeft2==0) {nHitS_left2 += sipm2->count;}
        else {nHitS_right2 += sipm2->count;}

        t2DhitS2->Fill(22*(moduleNum2%4)+fiberNum2, 24*(moduleNum2/4)+plateNum2, sipm2->count);
     THits2 -> Fill(moduleNum2,nHits2);   
    }
 //   THits2 -> Fill(moduleNum2,nHits2);
    }

    tHit_S2->Fill(nHitS2);
    tHit_S_left2->Fill(nHitS_left2);
    tHit_S_right2->Fill(nHitS_right2);
} // 2nd event root

  RootInterface<DRsimInterface::DRsimEventData>* drInterface3 = new RootInterface<DRsimInterface::DRsimEventData>(std::string(filename3 + ".root"), true);
  drInterface3->set("DRsim","DRsimEventData");

  unsigned int entries3 = drInterface3->entries();
  while (drInterface3->numEvt() < entries) {
    if (drInterface3->numEvt() % 1000 == 0) printf("3rd analyzing %dth event ...\n", drInterface3->numEvt());

    DRsimInterface::DRsimEventData drEvt;
    drInterface3->read(drEvt);

    float Edep3 = 0.; float totE3 = 0.; float ratE3 = 0;
    for (auto edepItr3 = drEvt.Edeps.begin(); edepItr3 != drEvt.Edeps.end(); ++edepItr3) {
      auto edep3 = *edepItr3;
      Edep3 += edep3.Edep;
    }
    tEdep3->Fill(Edep3);

    float Pleak3 = 0.;
    float Eleak_nu3 = 0.;
    for (auto leak : drEvt.leaks) {
      TLorentzVector leak4vec;
      leak4vec.SetPxPyPzE(leak.px,leak.py,leak.pz,leak.E);
      if ( std::abs(leak.pdgId)==12 || std::abs(leak.pdgId)==14 || std::abs(leak.pdgId)==16 ) {
        Eleak_nu3 += leak4vec.P();
      } else {
        Pleak3 += leak4vec.P();
      }
    }
    tP_leak3->Fill(Pleak3);
    tP_leak_nu3->Fill(Eleak_nu3);
    totE3 = Pleak3 + Eleak_nu3 + Edep3;
    tot_E3 -> Fill(totE3) ;
    ratE3 = totE3 / inE ;
    ratio_E3 -> Fill(drInterface3->numEvt(),ratE3);

  int nHitS3 = 0; int isLeft3 = 0; int nHits3 = 0;
    int nHitS_left3 = 0; int nHitS_right3 = 0;
    for (auto tower3 = drEvt.towers.begin(); tower3 != drEvt.towers.end(); ++tower3) {
      int moduleNum3 = tower3->ModuleNum;
      for (auto sipm3 = tower3->SiPMs.begin(); sipm3 != tower3->SiPMs.end(); ++sipm3) {
        isLeft3 = sipm3->isleft;
        int plateNum3 = sipm3->x; int fiberNum3 = sipm3->y;
        tNhit_S3->Fill(sipm3->count);
        nHitS3 += sipm3->count;
        nHits3 += sipm3->count;
        if (isLeft3==0) {nHitS_left3 += sipm3->count;}
        else {nHitS_right3 += sipm3->count;}

        t2DhitS3->Fill(22*(moduleNum3%4)+fiberNum3, 24*(moduleNum3/4)+plateNum3, sipm3->count);

   THits3 -> Fill(moduleNum3,nHits3);
}
//    THits3 -> Fill(moduleNum3,nHits3);
    }
    tHit_S3->Fill(nHitS3);
    tHit_S_left3->Fill(nHitS_left3);
    tHit_S_right3->Fill(nHitS_right3);

} // 3rd event root

  RootInterface<DRsimInterface::DRsimEventData>* drInterface4 = new RootInterface<DRsimInterface::DRsimEventData>(std::string(filename4 + ".root"), true);
  drInterface4->set("DRsim","DRsimEventData");

  unsigned int entries4 = drInterface4->entries();
  while (drInterface4->numEvt() < entries) {
    if (drInterface4->numEvt() % 1000 == 0) printf("4rd analyzing %dth event ...\n", drInterface4->numEvt());

    DRsimInterface::DRsimEventData drEvt;
    drInterface4->read(drEvt);

    float Edep4 = 0.; float totE4 = 0.; float ratE4 = 0;
    for (auto edepItr4 = drEvt.Edeps.begin(); edepItr4 != drEvt.Edeps.end(); ++edepItr4) {
      auto edep4 = *edepItr4;
      Edep4 += edep4.Edep;
    }
    tEdep4->Fill(Edep4);

    float Pleak4 = 0.;
    float Eleak_nu4 = 0.;
    for (auto leak : drEvt.leaks) {
      TLorentzVector leak4vec;
      leak4vec.SetPxPyPzE(leak.px,leak.py,leak.pz,leak.E);
      if ( std::abs(leak.pdgId)==12 || std::abs(leak.pdgId)==14 || std::abs(leak.pdgId)==16 ) {
        Eleak_nu4 += leak4vec.P();
      } else {
        Pleak4 += leak4vec.P();
      }
    }
    tP_leak4->Fill(Pleak4);
    tP_leak_nu4->Fill(Eleak_nu4);
    totE4 = Pleak3 + Eleak_nu4 + Edep4;
    tot_E4 -> Fill(totE4) ;
    ratE4 = totE4 / inE ;
    ratio_E4 -> Fill(drInterface4->numEvt(),ratE4);

  int nHitS4 = 0; int isLeft4 = 0; int nHits4 = 0;
    int nHitS_left4 = 0; int nHitS_right4 = 0;
    for (auto tower4 = drEvt.towers.begin(); tower4 != drEvt.towers.end(); ++tower4) {
      int moduleNum4 = tower4->ModuleNum;
      for (auto sipm4 = tower4->SiPMs.begin(); sipm4 != tower4->SiPMs.end(); ++sipm4) {
        isLeft4 = sipm4->isleft;
        int plateNum4 = sipm4->x; int fiberNum4 = sipm4->y;
        tNhit_S3->Fill(sipm4->count);
        nHitS4 += sipm4->count;
        nHits4 += sipm4->count;
        if (isLeft4==0) {nHitS_left4 += sipm4->count;}
        else {nHitS_right4 += sipm4->count;}

//        t2DhitS4->Fill(22*(moduleNum3%4)+fiberNum3, 24*(moduleNum3/4)+plateNum3, sipm3->count);

   THits4 -> Fill(moduleNum4,nHits4);
}
   }
    tHit_S4->Fill(nHitS4);
    tHit_S_left4->Fill(nHitS_left4);
    tHit_S_right4->Fill(nHitS_right4);

} // 4rd event root


  Float_t l1 =  tP_leak->GetMean() + tP_leak_nu->GetMean() ;
  Float_t l2 =  tP_leak2->GetMean() + tP_leak_nu2->GetMean();
  Float_t l3 =  tP_leak3->GetMean() + tP_leak_nu3->GetMean();
  leak_ratio -> Fill(8,l1/50);
  leak_ratio2 -> Fill(1,l2/50);
  leak_ratio3 -> Fill(15,l3/50);


  TCanvas* c = new TCanvas("c","");

  tEdep->Draw("Hist");  tEdep2->Draw("Hist same"); tEdep3->Draw("Hist same"); c->SaveAs(filename1+"compare_Edep.png");   
  c->SetLogy(1);
  tP_leak->Draw("Hist");tP_leak2->Draw("Hist same");tP_leak3->Draw("Hist same");  c->SaveAs(filename1+"compare_Pleak.png");
  tP_leak_nu->Draw("Hist"); tP_leak_nu2->Draw("Hist same"); tP_leak_nu3->Draw("Hist same");  c->SaveAs(filename1+"compare_Pleak_nu.png");
  c->SetLogy(0);
  tot_E3->Draw("Hist"); tot_E2->Draw("Hist same"); tot_E->Draw("Hist same"); c->SaveAs(filename1+"compare_totE.png");
  tot_E->SetMarkerColor(2); tot_E2->SetMarkerColor(3); tot_E3->SetMarkerColor(4);
  THits->Draw("Hist"); THits2->Draw("Hist same"); THits3->Draw("Hist same");  c->SaveAs(filename1+"compare_towerE.png");
  ratio_E -> SetMarkerStyle(20); ratio_E2 -> SetMarkerStyle(20); ratio_E3 -> SetMarkerStyle(20);
  ratio_E -> SetMarkerSize(0.3); ratio_E2 -> SetMarkerSize(0.3); ratio_E3 -> SetMarkerSize(0.3);
  ratio_E -> SetMarkerColor(kRed); ratio_E2 -> SetMarkerColor(kBlue);ratio_E3 -> SetMarkerColor(kGreen);
  ratio_E -> GetYaxis() -> SetRangeUser(0.7,1.3);
  ratio_E-> Draw("hist p"); ratio_E2-> Draw("hist p same"); ratio_E3-> Draw("hist p same");
  TLine *tl = new TLine (0,1,100,1);
  tl -> Draw(); c->SaveAs(filename1+"compare_ratioE.png");  
  
  tHit_S->Draw("Hist"); tHit_S2->Draw("Hist same");tHit_S3->Draw("Hist same");c->SaveAs(filename1+"compare_nHitpEventS.pdf");
  tHit_S_left->Draw("Hist"); tHit_S_left2->Draw("Hist same"); tHit_S_left3->Draw("Hist same"); c->SaveAs(filename1+"compare_nHitpEventS_left.pdf");
  tHit_S_right->Draw("Hist");tHit_S_right2->Draw("Hist same");tHit_S_right3->Draw("Hist same"); c->SaveAs(filename1+"compare_nHitpEventS_right.pdf");

  t2DhitS->Draw("COLZ"); c->SaveAs(filename1+"_n2DHitS.png");
  t2DhitS2->Draw("COLZ"); c->SaveAs(filename2+"_n2DHitS.png");
  t2DhitS3->Draw("COLZ"); c->SaveAs(filename3+"_n2DHitS.png");

 leak_ratio3 -> SetStats(0);
 leak_ratio3 -> Draw("hist");  leak_ratio2 -> Draw("hist same");  leak_ratio -> Draw("hist same");  c->SaveAs(filename1+"compare_LeakageP.png");

std::cout << "z1  :  "   << leak_ratio2->GetBinContent(2)   << "     z8  :  "     <<  leak_ratio->GetBinContent(9)   << "    z15  =  "    << leak_ratio3->GetBinContent(16)  << std::endl;

/*
  tT_S->Draw("Hist"); c->SaveAs(filename+"_tS.png");
  tT_S_left->Draw("Hist"); c->SaveAs(filename+"_tS_left.png");
  tT_S_right->Draw("Hist"); c->SaveAs(filename+"_tS_right.png");
  tWav_S->Draw("Hist"); c->SaveAs(filename+"_wavS.png");
  tWav_S_left->Draw("Hist"); c->SaveAs(filename+"_wavS_left.png");
  tWav_S_right->Draw("Hist"); c->SaveAs(filename+"_wavS_right.png");
  tNhit_S->Draw("Hist"); c->SaveAs(filename+"_nhitS.png");*/
}
