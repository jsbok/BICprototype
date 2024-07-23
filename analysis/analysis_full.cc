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
  TString filename1 = "/u/user/changhui/work/BICprototype/rootfiles/3by5_e-_500_MeV_final/root/R3by5_e-_500_MeV_final_10M" ;
  TString filename2 = "/u/user/changhui/work/BICprototype/rootfiles/3by5_e-_4_GeV_final/root/R4by4_e-_4_GeV_final_10M" ;
  TString filename3 = "/u/user/changhui/work/BICprototype/rootfiles/2by8_e-_3_GeV_final/root/R2by8_e-_3_GeV_final_10M" ;
  TString filename4 = "/u/user/changhui/work/BICprototype/rootfiles/3by5_e-_2_GeV_final/root/R3by5_e-_2_GeV_final_10M" ;
  TString filename5 = "/u/user/changhui/work/BICprototype/rootfiles/3by5_e-_1_GeV_final/root/R4by4_e-_1_GeV_final_10M" ;
  TString filename6 = "/u/user/changhui/work/BICprototype/rootfiles/4by4_e-_500_MeV_final/root/R4by4_e-_500_MeV_final_10M" ;
  int file1=1; int file2=0; int file3=0; int file4=1; int file5=1; int file6=1; // File loop on/off

  float inE = std::stof(argv[1]);
  float low = std::stof(argv[2]);
  float high = std::stof(argv[3]);

  int e1 = inE*5;int e2 = inE*4; int e3 = inE *3 ; int e4 = inE*2 ; int e5 = inE; int e6 = inE/2; // Beam E of each file
  int nLayers = 22;  int nFibers = 24;  int nRows = 2;  int nColumns = 8; // module configuration


  gStyle->SetOptFit(1);

  RootInterface<DRsimInterface::DRsimEventData>* drInterface = new RootInterface<DRsimInterface::DRsimEventData>(std::string(filename1 + ".root"), true);
  drInterface->set("DRsim","DRsimEventData");
  //drInterface->GetChain("DRsim");

  TH1F* tEdep = new TH1F("totEdep",";MeV;Evt",80,low*1000.,high*800.);
  tEdep->Sumw2(); tEdep->SetLineColor(kRed); tEdep->SetLineWidth(2);
  TH1F* tEdep4 = new TH1F("totEdep4",";MeV;Evt",80,low*1000.,high*800.);
  tEdep4->Sumw2(); tEdep4->SetLineColor(1); tEdep4->SetLineWidth(2);
  TH1F* tEdep5 = new TH1F("totEdep5",";MeV;Evt",80,low*1000.,high*800.);
  tEdep5->Sumw2(); tEdep5->SetLineColor(5); tEdep5->SetLineWidth(2);

  TH1F* tot_E = new TH1F("tot_E",";MeV;Evt",100,high*400.,high*600.);
  tot_E->Sumw2(); tot_E->SetLineColor(kRed); tot_E->SetLineWidth(2);
  TH1D* ratio_E = new TH1D("ratio_E",";evtN;ratio",100.,0,100.);
  ratio_E->Sumw2(); ratio_E->SetLineColor(kRed);

  TH1F* tEdep2 = new TH1F("totEdep2",";MeV;Evt",80,low*1000.,high*800.);
  tEdep2->Sumw2(); tEdep2->SetLineColor(kBlue); tEdep2->SetLineWidth(2);
  TH1F* tot_E2 = new TH1F("tot_E2",";MeV;Evt",100,high*400.,high*600.);
  tot_E2->Sumw2(); tot_E2->SetLineColor(kBlue); tot_E2->SetLineWidth(2);
  TH1D* ratio_E2 = new TH1D("ratio_E",";evtN;ratio",100.,0,100.);
  ratio_E2->Sumw2(); ratio_E2->SetLineColor(kBlue);

  TH1F* tEdep3 = new TH1F("totEdep3",";MeV;Evt",80,low*1000.,high*800.);
  tEdep3->Sumw2(); tEdep3->SetLineColor(kGreen); tEdep3->SetLineWidth(2);
  TH1F* tot_E3 = new TH1F("tot_E3",";MeV;Evt",100,high*400.,high*600.);
  tot_E3->Sumw2(); tot_E3->SetLineColor(kGreen); tot_E3->SetLineWidth(2);
  TH1D* ratio_E3 = new TH1D("ratio_E3",";evtN;ratio",100.,0,100.);
  ratio_E3->Sumw2(); ratio_E3->SetLineColor(kGreen);

  TH1F* tHit_S = new TH1F("Hit_S",";# of p.e.;Evt",150,12000*low*2,9000*high*2);
  tHit_S->Sumw2(); tHit_S->SetLineColor(2); tHit_S->SetLineWidth(2);
  TH1F* tHit_S2 = new TH1F("Hit_S2",";# of p.e.;Evt",150,12000*low*2,9000*high*2);
  tHit_S2->Sumw2(); tHit_S2->SetLineColor(kBlue); tHit_S2->SetLineWidth(2);
  TH1F* tHit_S3 = new TH1F("Hit_S3",";# of p.e.;Evt",150,12000*low*2,9000*high*2);
  tHit_S3->Sumw2(); tHit_S3->SetLineColor(kGreen); tHit_S3->SetLineWidth(2);
  TH1F* tHit_S4 = new TH1F("Hit_S4",";# of p.e.;Evt",150,12000*low*2,9000*high*2);
  tHit_S4->Sumw2(); tHit_S4->SetLineColor(1); tHit_S4->SetLineWidth(2);
  TH1F* tHit_S5 = new TH1F("Hit_S5",";# of p.e.;Evt",150,12000*low*2,9000*high*2);
  tHit_S5->Sumw2(); tHit_S5->SetLineColor(5); tHit_S5->SetLineWidth(2);


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
 TH1F* THits4 = new TH1F("THits4",";ModuleN; # of p.e.",20,0,20);
 THits4->Sumw2(); THits4->SetLineColor(1); THits4->SetLineWidth(2);
 TH1F* THits5 = new TH1F("THits5",";ModuleN; # of p.e.",20,0,20);
 THits5->Sumw2(); THits5->SetLineColor(5); THits5->SetLineWidth(2);

  TH1F* tP_leak = new TH1F("Pleak",";MeV;Evt",80,0.,800.*high);
  tP_leak->Sumw2(); tP_leak->SetLineWidth(2); tP_leak->SetLineColor(kRed);
  TH1F* tP_leak2 = new TH1F("Pleak2",";MeV;Evt",80,0.,800.*high);
  tP_leak2->Sumw2(); tP_leak2->SetLineWidth(2); tP_leak2->SetLineColor(kBlue);
  TH1F* tP_leak3 = new TH1F("Pleak3",";MeV;Evt",80,0.,800.*high);
  tP_leak3->Sumw2(); tP_leak3->SetLineWidth(2); tP_leak3->SetLineColor(kGreen);
  TH1F* tP_leak4 = new TH1F("Pleak4",";MeV;Evt",80,0.,800.*high);
  tP_leak4->Sumw2(); tP_leak4->SetLineWidth(2); tP_leak4->SetLineColor(1);
  TH1F* tP_leak5 = new TH1F("Pleak5",";MeV;Evt",80,0.,1000.*high);
  tP_leak5->Sumw2(); tP_leak5->SetLineWidth(2); tP_leak5->SetLineColor(5);
/*
  TH1F* leak_ratio = new TH1F("leak_ratio","; z0; leakage (%)",20,0.,20.);
  leak_ratio->Sumw2(); leak_ratio->SetLineWidth(2); leak_ratio->SetLineColor(kRed);
  TH1F* leak_ratio2 = new TH1F("leak_ratio","; z0; leakage (%)",20,0.,20.);
  leak_ratio2->Sumw2(); leak_ratio2->SetLineWidth(2); leak_ratio2->SetLineColor(kBlue);
  TH1F* leak_ratio3 = new TH1F("leak_ratio","; z0;leakage (%)",20,0.,20.);
  leak_ratio3->Sumw2(); leak_ratio3->SetLineWidth(2); leak_ratio3->SetLineColor(kGreen);
*/
  TH1F* Edep_ratio = new TH1F("Edep_ratio","; Beam E (MeV); Edep (%)",60,low*1000.,high*600.);
  Edep_ratio->Sumw2(); Edep_ratio->SetLineWidth(1); Edep_ratio->SetLineColor(2);
  TH1F* Edep_ratio2 = new TH1F("Edep_ratio2","; Beam E (MeV); Edep (%)",60,low*1000.,high*600.);
  Edep_ratio2->Sumw2(); Edep_ratio2->SetLineWidth(1); Edep_ratio2->SetLineColor(3);
  TH1F* Edep_ratio3 = new TH1F("Edep_ratio3","; Beam E (MeV); Edep (%)",100,low*1000.,high*1000.);
  Edep_ratio3->Sumw2(); Edep_ratio3->SetLineWidth(1); Edep_ratio3->SetLineColor(3);
  TH1F* Edep_ratio4 = new TH1F("Edep_ratio4","; Beam E (MeV); Edep (%)",100,low*1000.,high*1000.);
  Edep_ratio4->Sumw2(); Edep_ratio4->SetLineWidth(1); Edep_ratio4->SetLineColor(1);
  TH1F* Edep_ratio5 = new TH1F("Edep_ratio5","; Beam E (MeV); Edep (%)",60,low*1000.,high*600.);
  Edep_ratio5->Sumw2(); Edep_ratio5->SetLineWidth(1); Edep_ratio5->SetLineColor(5);
    
  TH1F* tNhit_S = new TH1F("nHits_S",";p.e.;n",200,0.,200.);
  tNhit_S->Sumw2(); tNhit_S->SetLineColor(kRed); tNhit_S->SetLineWidth(2);
  TH1F* tNhit_S2 = new TH1F("nHits_S2",";p.e.;n",200,0.,200.);
  tNhit_S2->Sumw2(); tNhit_S2->SetLineColor(kBlue); tNhit_S2->SetLineWidth(2);
  TH1F* tNhit_S3 = new TH1F("nHits_S3",";p.e.;n",200,0.,200.);
  tNhit_S3->Sumw2(); tNhit_S3->SetLineColor(kGreen); tNhit_S3->SetLineWidth(2);
 
  float rat_E = 0;  float rat_E2 = 0; float rat_E3 = 0; float rat_E4 = 0; float rat_E5 = 0;
  double xBins = nFibers * nColumns;  double xLower = 0;  double xUpper = nFibers * nColumns;  double yBins = nLayers * nRows;
  double yLower = 0;  double yUpper = nLayers * nRows;


 TH2D* t2DhitS = new TH2D("2D Hit S1", "", xBins, xLower, xUpper, yBins, yLower, yUpper); t2DhitS->Sumw2(); t2DhitS->SetStats(0);
 TH2D* t2DhitS2 = new TH2D("2D Hit S2", "", xBins, xLower, xUpper, yBins, yLower, yUpper); t2DhitS2->Sumw2(); t2DhitS2->SetStats(0);
 TH2D* t2DhitS3 = new TH2D("2D Hit S3", "", xBins, xLower, xUpper, yBins, yLower, yUpper); t2DhitS3->Sumw2(); t2DhitS3->SetStats(0);
 TH2D* t2DhitS4 = new TH2D("2D Hit S4", "", xBins, xLower, xUpper, yBins, yLower, yUpper); t2DhitS4->Sumw2(); t2DhitS4->SetStats(0);
 TH2D* t2DhitS5 = new TH2D("2D Hit S5", "", xBins, xLower, xUpper, yBins, yLower, yUpper); t2DhitS5->Sumw2(); t2DhitS5->SetStats(0);


if (file1==1) {
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

    float Eleak = 0.;
    float Pleak = 0.;
    float Eleak_nu = 0.;
    for (auto leak : drEvt.leaks) {
      TLorentzVector leak4vec;
      leak4vec.SetPxPyPzE(leak.px,leak.py,leak.pz,leak.E);
/*     
   if ( drInterface->numEvt() == 1207 || drInterface->numEvt() == 3530 || drInterface->numEvt() == 5799 || drInterface->numEvt() == 7490 || drInterface->numEvt() == 9749) {
        std::cout << " P : "<< leak4vec.P() << ", kinE " << leak.kE  << ", E :" << leak.E  << ", x :" << leak.vx << ", y : " << leak.vy << ", z : " << leak.vz << ", pid : " << leak.pdgId << ", Tid : " << leak.Tid << ", time : " << leak.vt << std::endl;
  }*/                                // leak particle informations
//        Pleak += leak4vec.P();     
          Eleak += leak.kE;
    }
// if (Pleak > 5100) {std::cout << "evtN :  " << drInterface->numEvt()  <<std::endl; }
       
     tP_leak->Fill(Pleak);
//     totE = Pleak + Edep;
//     tot_E -> Fill(totE) ;
//    ratE = totE / e ;
//    ratio_E -> Fill(drInterface->numEvt(),ratE);

    int nHitS = 0; int isLeft = 0; int nHits;
    int nHitS_left = 0; int nHitS_right = 0;
    for (auto tower = drEvt.towers.begin(); tower != drEvt.towers.end(); ++tower) {
      int moduleNum = tower->ModuleNum; nHits=0;
      for (auto sipm = tower->SiPMs.begin(); sipm != tower->SiPMs.end(); ++sipm) {
        isLeft = sipm->isleft;
        int plateNum = sipm->y; int fiberNum = 24 - sipm->x; 

//std::cout <<" fiberNum : " << fiberNum << std::endl ;
//std::cout <<" plateNum : " << plateNum << std::endl ;

//        tNhit_S->Fill(sipm->count);
        nHitS += sipm->count;
        nHits += sipm->count;   
//        if (isLeft==0) {nHitS_left += sipm->count;}
//        else {nHitS_right += sipm->count;}

        t2DhitS->Fill(nFibers*(moduleNum%nColumns)+fiberNum, nLayers*(moduleNum/nColumns)+plateNum, sipm->count);
 
     }//SiPM loop
THits -> Fill(moduleNum,nHits);
    }//tower loop

    tHit_S->Fill(nHitS);
//    tHit_S_left->Fill(nHitS_left);
//    tHit_S_right->Fill(nHitS_right);
  } // event loop
 
    rat_E = tEdep->GetMean() / e1 ;
    Edep_ratio -> SetBinContent(50,rat_E*100);
    float a = Edep_ratio->GetBinContent(50);
    std::cout << "BinCon = " << a << std::endl;
    Edep_ratio -> SetBinError(50, tEdep->GetStdDev()*100/e1);
    std::cout << "BinErr = " << tEdep->GetStdDev()*100/e1  << std::endl;

 }
if (file2==1) {
// second file
  RootInterface<DRsimInterface::DRsimEventData>* drInterface2 = new RootInterface<DRsimInterface::DRsimEventData>(std::string(filename2 + ".root"), true);
  drInterface2->set("DRsim","DRsimEventData"); 

  unsigned int entries2 = drInterface2->entries();
  while (drInterface2->numEvt() < entries2) {
    if (drInterface2->numEvt() % 1000 == 0) printf("2nd analyzing %dth event ...\n", drInterface2->numEvt());

    DRsimInterface::DRsimEventData drEvt;
    drInterface2->read(drEvt);

    float Edep2 = 0.; float totE2 = 0.; float ratE2 = 0; float rat_E2 = 0;
    for (auto edepItr2 = drEvt.Edeps.begin(); edepItr2 != drEvt.Edeps.end(); ++edepItr2) {
      auto edep2 = *edepItr2;
      Edep2 += edep2.Edep;
    }
    tEdep2->Fill(Edep2);

    float Pleak2 = 0.;
    float Eleak2 = 0.;
    for (auto leak : drEvt.leaks) {
      TLorentzVector leak4vec;
      leak4vec.SetPxPyPzE(leak.px,leak.py,leak.pz,leak.E);
      
//        Pleak2 += leak4vec.P();
        Eleak2 += leak.kE;
    }
    tP_leak2->Fill(Eleak2);
//    tP_leak_nu2->Fill(Eleak_nu2);
//    totE2 = Pleak2 + Eleak_nu2 + Edep2;
//    tot_E2 -> Fill(totE2) ;
//    ratE2 = totE2 / e2 ;
 //   ratio_E2 -> Fill(drInterface2->numEvt(),ratE2);

  int nHitS2 = 0; int isLeft2 = 0; int nHits2 = 0;
    int nHitS_left2 = 0; int nHitS_right2 = 0;
    for (auto tower2 = drEvt.towers.begin(); tower2 != drEvt.towers.end(); ++tower2) {
      int moduleNum2 = tower2->ModuleNum;
      for (auto sipm2 = tower2->SiPMs.begin(); sipm2 != tower2->SiPMs.end(); ++sipm2) {
        isLeft2 = sipm2->isleft;
        int plateNum2 = sipm2->y; int fiberNum2 = 24 - sipm2->x;
//        tNhit_S2->Fill(sipm2->count);
        nHitS2 += sipm2->count;
	nHits2 += sipm2->count;
//        if (isLeft2==0) {nHitS_left2 += sipm2->count;}
//        else {nHitS_right2 += sipm2->count;}

     t2DhitS2->Fill(nFibers*(moduleNum2%nColumns)+fiberNum2, nLayers*(moduleNum2/nColumns)+plateNum2, sipm2->count);
     THits2 -> Fill(moduleNum2,nHits2);   
    }
    }

    tHit_S2->Fill(nHitS2);
//    tHit_S_left2->Fill(nHitS_left2);
//    tHit_S_right2->Fill(nHitS_right2);
} // 2nd event root

    rat_E2 = tEdep2->GetMean() / e2 ;
    Edep_ratio -> SetBinContent(40,rat_E2*100);
    Edep_ratio -> SetBinError(40, tEdep2->GetStdDev()*100/e2);
    std::cout << "4GeV BinCon = " << rat_E2*100 << std::endl;
    std::cout << "4GeV BinErr = " << tEdep2->GetStdDev()*100/e2  << std::endl;

}
if (file3==1) {
  RootInterface<DRsimInterface::DRsimEventData>* drInterface3 = new RootInterface<DRsimInterface::DRsimEventData>(std::string(filename3 + ".root"), true);
  drInterface3->set("DRsim","DRsimEventData");

  unsigned int entries3 = drInterface3->entries();
  while (drInterface3->numEvt() < entries3) {
    if (drInterface3->numEvt() % 1000 == 0) printf("3rd analyzing %dth event ...\n", drInterface3->numEvt());

    DRsimInterface::DRsimEventData drEvt;
    drInterface3->read(drEvt);

    float Edep3 = 0.; float totE3 = 0.; float ratE3 = 0; float rat_E3 = 0; 
    for (auto edepItr3 = drEvt.Edeps.begin(); edepItr3 != drEvt.Edeps.end(); ++edepItr3) {
      auto edep3 = *edepItr3;
      Edep3 += edep3.Edep;
    }
    tEdep3->Fill(Edep3);

    float Pleak3 = 0.;
    float Eleak3 = 0.;
    for (auto leak : drEvt.leaks) {
      TLorentzVector leak4vec;
      leak4vec.SetPxPyPzE(leak.px,leak.py,leak.pz,leak.E);
        Eleak3 += leak.kE;
 //       Pleak3 += leak4vec.P();
/*if ( drInterface3->numEvt() == 327 || drInterface->numEvt() == 1445 || drInterface->numEvt() == 2726 || drInterface->numEvt() == 6043 || drInterface->numEvt() == 8115) {
        std::cout << " P : "<< leak4vec.P() << ", kinE " << leak.kE  << ", E :" << leak.E  << ", x :" << leak.vx << ", y : " << leak.vy << ", z : " << leak.vz << ", pid : " << leak.pdgId << ", Tid : " << leak.Tid << ", time : " << leak.vt << std::endl;
  }      
if (drInterface3->numEvt() == 9111) { std::cout << " P : "<< leak4vec.P() << ", kinE " << leak.kE  << ", E :" << leak.E  << ", x :" << leak.vx << ", y : " << leak.vy << ", z : " << leak.vz << ", pid : " << leak.pdgId << ", xxxxxxxxxx" << std::endl;
  }*/
}
    
//  if (Pleak3 > 3100) {std::cout << "evtN :  " << drInterface3->numEvt()  <<std::endl; }
    tP_leak3->Fill(Eleak3);
//    totE3 = Eleak3 + Edep3;
//    tot_E3 -> Fill(totE3) ;
//    ratE3 = totE3 / e3 ;
//    ratio_E3 -> Fill(drInterface3->numEvt(),ratE3);

  int nHitS3 = 0; int isLeft3 = 0; int nHits3;
    int nHitS_left3 = 0; int nHitS_right3 = 0;
    for (auto tower3 = drEvt.towers.begin(); tower3 != drEvt.towers.end(); ++tower3) {
      int moduleNum3 = tower3->ModuleNum; nHits3=0;
      for (auto sipm3 = tower3->SiPMs.begin(); sipm3 != tower3->SiPMs.end(); ++sipm3) {
        isLeft3 = sipm3->isleft;
         int plateNum = sipm3->y; int fiberNum = 24 - sipm3->x;
//        tNhit_S3->Fill(sipm3->count);
        nHitS3 += sipm3->count;
        nHits3 += sipm3->count;
//        if (isLeft3==0) {nHitS_left3 += sipm3->count;}
//        else {nHitS_right3 += sipm3->count;}
 t2DhitS3->Fill(nFibers*(moduleNum3%nColumns)+fiberNum, nLayers*(moduleNum3/nColumns)+plateNum, sipm3->count);
}
   THits3 -> Fill(moduleNum3,nHits3);
    }
    tHit_S3->Fill(nHitS3);
//    tHit_S_left3->Fill(nHitS_left3);
//    tHit_S_right3->Fill(nHitS_right3);

} // 3rd event root
}
if (file4==1) {
  RootInterface<DRsimInterface::DRsimEventData>* drInterface4 = new RootInterface<DRsimInterface::DRsimEventData>(std::string(filename4 + ".root"), true);
  drInterface4->set("DRsim","DRsimEventData");

  unsigned int entries4 = drInterface4->entries();
  while (drInterface4->numEvt() < entries4) {
    if (drInterface4->numEvt() % 1000 == 0) printf("4th analyzing %dth event ...\n", drInterface4->numEvt());

    DRsimInterface::DRsimEventData drEvt;
    drInterface4->read(drEvt);

    float Edep4 = 0.; float totE4 = 0.; float ratE4 = 0; float rat_E4 = 0;
    for (auto edepItr4 = drEvt.Edeps.begin(); edepItr4 != drEvt.Edeps.end(); ++edepItr4) {
      auto edep4 = *edepItr4;
      Edep4 += edep4.Edep;
    }
    tEdep4->Fill(Edep4);

   float Pleak4 = 0.;    float Eleak4 = 0.;
    for (auto leak : drEvt.leaks) {
      TLorentzVector leak4vec;
      leak4vec.SetPxPyPzE(leak.px,leak.py,leak.pz,leak.E);
     
//        Pleak4 += leak4vec.P();
        Eleak4 += leak.kE;
    }
    tP_leak4->Fill(Eleak4);
    
//    totE4 = Pleak4 + Edep4;
//    tot_E4 -> Fill(totE4) ;
//    ratE3 = totE4 / e4 ;
//    ratio_E4 -> Fill(drInterface4->numEvt(),ratE4);

  int nHitS4 = 0; int isLeft4 = 0; int nHits4 = 0;
    int nHitS_left4 = 0; int nHitS_right4 = 0;
    for (auto tower4 = drEvt.towers.begin(); tower4 != drEvt.towers.end(); ++tower4) {
      int moduleNum4 = tower4->ModuleNum;
      for (auto sipm4 = tower4->SiPMs.begin(); sipm4 != tower4->SiPMs.end(); ++sipm4) {
        isLeft4 = sipm4->isleft;
        int plateNum4 = sipm4->y; int fiberNum4 = 24 - sipm4->x;
//        tNhit_S4->Fill(sipm4->count);
        nHitS4 += sipm4->count;
        nHits4 += sipm4->count;
//        if (isLeft4==0) {nHitS_left4 += sipm4->count;}
//        else {nHitS_right4 += sipm4->count;}

        t2DhitS4->Fill(22*(moduleNum4%5)+fiberNum4, 24*(moduleNum4/5)+plateNum4, sipm4->count);

   THits4 -> Fill(moduleNum4,nHits4);
}
    }
    tHit_S4->Fill(nHitS4);
//    tHit_S_left3->Fill(nHitS_left3);
//    tHit_S_right3->Fill(nHitS_right3);

} // 4th event root
}
if (file5==1) {
    RootInterface<DRsimInterface::DRsimEventData>* drInterface5 = new RootInterface<DRsimInterface::DRsimEventData>(std::string(filename5 + ".root"), true);
  drInterface5->set("DRsim","DRsimEventData");

  unsigned int entries5 = drInterface5->entries();
  while (drInterface5->numEvt() < entries5) {
    if (drInterface5->numEvt() % 1000 == 0) printf("5th analyzing %dth event ...\n", drInterface5->numEvt());

    DRsimInterface::DRsimEventData drEvt;
    drInterface5->read(drEvt);

    float Edep5 = 0.; float totE5 = 0.; float ratE5 = 0; float rat_E5 = 0; float PE = 0;
    for (auto edepItr5 = drEvt.Edeps.begin(); edepItr5 != drEvt.Edeps.end(); ++edepItr5) {
      auto edep5 = *edepItr5;
      Edep5 += edep5.Edep;
    }
    tEdep5->Fill(Edep5);

    float Pleak5 = 0.;    float Eleak5 = 0.;
    for (auto leak : drEvt.leaks) {
      TLorentzVector leak4vec;
      leak4vec.SetPxPyPzE(leak.px,leak.py,leak.pz,leak.E);
      
//        Pleak5 += leak4vec.P();
        Eleak5 += leak.kE;

    }
    tP_leak5->Fill(Eleak5); 
//    totE5 = Pleak5 + Edep5;
 //   tot_E5 -> Fill(totE5) ;
//    ratE5 = totE5 / inE ;
//    ratio_E5 -> Fill(drInterface5->numEvt(),ratE5);

  int nHitS5 = 0; int isLeft5 = 0; int nHits5 = 0;
    int nHitS_left5 = 0; int nHitS_right5 = 0;
    for (auto tower5 = drEvt.towers.begin(); tower5 != drEvt.towers.end(); ++tower5) {
      int moduleNum5 = tower5->ModuleNum;
      for (auto sipm5 = tower5->SiPMs.begin(); sipm5 != tower5->SiPMs.end(); ++sipm5) {
        isLeft5 = sipm5->isleft;
        int plateNum5 = sipm5->y; int fiberNum5 = 24- sipm5->x;
//        tNhit_S5->Fill(sipm5->count);
        nHitS5 += sipm5->count;
        nHits5 += sipm5->count;
//        if (isLeft5==0) {nHitS_left5 += sipm5->count;}
//        else {nHitS_right5 += sipm5->count;}

        t2DhitS5->Fill(24*(moduleNum5%4)+fiberNum5, 22*(moduleNum5/4)+plateNum5, sipm5->count);

   THits5 -> Fill(moduleNum5,nHits5);
}
    }
    tHit_S5->Fill(nHitS5);
//    tHit_S_left5->Fill(nHitS_left5);
//    tHit_S_right5->Fill(nHitS_right5);
} // 5th event root
}
//  Float_t l1 =  tP_leak->GetMean() + tP_leak_nu->GetMean() ;
//  Float_t l2 =  tP_leak2->GetMean() + tP_leak_nu2->GetMean();
//  Float_t l3 =  tP_leak3->GetMean() + tP_leak_nu3->GetMean();
//  leak_ratio -> Fill(8,l1/50);
//  leak_ratio2 -> Fill(1,l2/50);
//  leak_ratio3 -> Fill(15,l3/50);

   rat_E3 = tEdep3->GetMean() / e3 ;
   Edep_ratio -> SetBinContent(30,rat_E3*100);
   Edep_ratio -> SetBinError(30, tEdep3->GetStdDev()*100/e3);
   rat_E4 = tEdep4->GetMean() / e3 ;
   Edep_ratio4 -> SetBinContent(30,rat_E4*100);
   Edep_ratio4 -> SetBinError(30, tEdep4->GetStdDev()*100/e3);
//   rat_E5 = tEdep5->GetMean() / inE ;
//   Edep_ratio5 -> SetBinContent(10,rat_E5*100);
//   Edep_ratio5 -> SetBinError(10, tEdep5->GetStdDev()*100/inE);

    std::cout << "3GeV BinCon = " << rat_E3*100 << std::endl;
    std::cout << "3GeV BinErr = " << tEdep3->GetStdDev()*100/e3  << std::endl;

//    std::cout << "2GeV BinCon = " << rat_E4*100 << std::endl;
//    Edep_ratio -> SetBinError(10, tEdep5->GetStdDev()*100/inE);
//    std::cout << "1GeV BinCon = " << rat_E5*100 << std::endl;

   Edep_ratio2 -> SetBinContent(50, 86.87);
   Edep_ratio2 -> SetBinError(50, 5.40);
   Edep_ratio2 -> SetBinContent(40,87.56);
   Edep_ratio2 -> SetBinError(40,5.29);   
   Edep_ratio2 -> SetBinContent(30,88.54);
   Edep_ratio2 -> SetBinError(30, 4.91);   
/*   Edep_ratio2 -> SetBinContent(20,83.7442);
   Edep_ratio2 -> SetBinError(20, 9.63355); 
   Edep_ratio2 -> SetBinContent(10,87.1821);
   Edep_ratio2 -> SetBinError(10, 7.83431);
 */

  TCanvas* c = new TCanvas("c","");
  tEdep3->Draw("Hist");//  tEdep2->Draw("Hist same");
tEdep->Draw("Hist same");  tEdep4->Draw("Hist same"); tEdep5->Draw("Hist same"); c->SaveAs(filename1+"compare_Edep_e.png");   
  c->SetLogy(1);
  tP_leak->Draw("Hist"); tP_leak3->Draw("Hist same");tP_leak4->Draw("Hist same"); //tP_leak4->Draw("Hist same"); tP_leak3->Draw("Hist same"); 
c->SaveAs(filename1+"compare_Pleak_e.png");
/*//  tP_leak_nu->Draw("Hist"); tP_leak_nu2->Draw("Hist same"); tP_leak_nu3->Draw("Hist same");  c->SaveAs(filename1+"compare_Pleak_nu_e.png");
//  c->SetLogy(0);
*///  tot_E3->Draw("Hist");/* tot_E2->Draw("Hist same"); tot_E->Draw("Hist same");*/ c->SaveAs(filename1+"compare_totE_e.png");
/*   tot_E->SetMarkerColor(2); tot_E2->SetMarkerColor(3);*/// tot_E3->SetMarkerColor(4);

  c->SetLogy(0);
  THits3->Draw("Hist"); THits->Draw("Hist same"); THits4->Draw("Hist same");//  THits4->Draw("Hist same"); THits5->Draw("Hist same"); 
  THits3->SetStats(0);
 c->SaveAs(filename1+"compare_towerE_e.png");
/*  ratio_E -> SetMarkerStyle(20); ratio_E2 -> SetMarkerStyle(20); ratio_E3 -> SetMarkerStyle(20);
  ratio_E -> SetMarkerSize(0.3); ratio_E2 -> SetMarkerSize(0.3); ratio_E3 -> SetMarkerSize(0.3);
  ratio_E -> SetMarkerColor(kRed); ratio_E2 -> SetMarkerColor(kBlue); ratio_E3 -> SetMarkerColor(kGreen);
  ratio_E -> GetYaxis() -> SetRangeUser(0.7,1.3);
  ratio_E-> Draw("hist p"); ratio_E2-> Draw("hist p same"); ratio_E3-> Draw("hist p same");
  TLine *tl = new TLine (0,1,100,1);
  tl -> Draw(); c->SaveAs(filename1+"compare_ratioE_e.png");  
  */
   tHit_S3->Draw("Hist"); tHit_S2->Draw("Hist same");tHit_S5->Draw("Hist same"); tHit_S4->Draw("Hist same"); tHit_S->Draw("Hist same"); c->SaveAs(filename1 +"compare_nHitpEventS_e.png");
/*  tHit_S_left->Draw("Hist"); tHit_S_left2->Draw("Hist same"); tHit_S_left3->Draw("Hist same"); c->SaveAs(filename1+"compare_nHitpEventS_left_e.pdf");
  tHit_S_right->Draw("Hist");tHit_S_right2->Draw("Hist same");tHit_S_right3->Draw("Hist same"); c->SaveAs(filename1+"compare_nHitpEventS_right_e.pdf");
*/
 
  t2DhitS->Draw("COLZ"); c->SaveAs(filename1+"_n2DHitS.png");
  t2DhitS2->Draw("COLZ"); c->SaveAs(filename2+"_n2DHitS.png");
  t2DhitS3->Draw("COLZ"); c->SaveAs(filename3+"_n2DHitS.png");  
  t2DhitS4->Draw("COLZ"); c->SaveAs(filename4+"_n2DHitS.png");
  t2DhitS5->Draw("COLZ"); c->SaveAs(filename5+"_n2DHitS.png");
/*
 leak_ratio3 -> SetStats(0);
 leak_ratio3 -> Draw("hist");  leak_ratio2 -> Draw("hist same");  leak_ratio -> Draw("hist same");  c->SaveAs(filename1+"compare_LeakageP_e.png");

std::cout << "z1  :  "   << leak_ratio2->GetBinContent(2)   << "     z8  :  "     <<  leak_ratio->GetBinContent(9)   << "    z15  =  "    << leak_ratio3->GetBinContent(16)  << std::endl;
*/
Edep_ratio->Draw("p E1"); //Edep_ratio2->Draw("p E1 same ");
Edep_ratio2->Draw("p E1 same"); Edep_ratio4->Draw("p E1 same ");/* Edep_ratio5->Draw("p E1 same"); */
gStyle->SetErrorX(0); Edep_ratio -> SetStats(0);  
Edep_ratio->SetMarkerSize(1.3); Edep_ratio->SetMarkerStyle(20); Edep_ratio->GetYaxis()->SetRangeUser(60,100);Edep_ratio->SetMarkerColor(2);
Edep_ratio2->SetMarkerSize(1.3); Edep_ratio2->SetMarkerStyle(20); Edep_ratio2->SetMarkerColor(3); //Edep_ratio3->SetMarkerSize(1.3); Edep_ratio3->SetMarkerStyle(20); Edep_ratio3->SetMarkerColor(4); Edep_ratio3->SetLineColor(4);
Edep_ratio4->SetMarkerSize(1.3); Edep_ratio4->SetMarkerStyle(20); Edep_ratio4->SetMarkerColor(1);
//Edep_ratio5->SetMarkerSize(1.3); Edep_ratio5->SetMarkerStyle(20);
c->SaveAs(filename1+"compare_Edep_ratio_e.png");

  for ( int Bin = 0 ; Bin < 16 ; ++Bin) {
 std::cout << "Nhits " << Bin << " : "  << THits3->GetBinContent(Bin) <<std::endl;
} // Bin counts


/*
  tT_S->Draw("Hist"); c->SaveAs(filename+"_tS.png");
  tT_S_left->Draw("Hist"); c->SaveAs(filename+"_tS_left.png");
  tT_S_right->Draw("Hist"); c->SaveAs(filename+"_tS_right.png"); // -> time 
  tWav_S->Draw("Hist"); c->SaveAs(filename+"_wavS.png");
  tWav_S_left->Draw("Hist"); c->SaveAs(filename+"_wavS_left.png");
  tWav_S_right->Draw("Hist"); c->SaveAs(filename+"_wavS_right.png");  // -> wave
  tNhit_S->Draw("Hist"); c->SaveAs(filename+"_nhitS.png");*/  // -> Hits per time
}
