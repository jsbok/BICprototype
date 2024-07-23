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
  TString filename1 = "/u/user/changhui/work/BICprototype/rootfiles/4by4_e-_5_GeV_final/root/R4by4_e-_5_GeV_final_10M" ;
  TString filename2 = "/u/user/changhui/work/BICprototype/rootfiles/4by4_e-_4_GeV_final/root/R4by4_e-_4_GeV_final_10M" ;
  TString filename3 = "/u/user/changhui/work/BICprototype/rootfiles/4by4_e-_3_GeV_final/root/R4by4_e-_3_GeV_final_10M" ;
  TString filename4 = "/u/user/changhui/work/BICprototype/rootfiles/4by4_e-_2_GeV_final/root/R4by4_e-_2_GeV_final_10M" ;
  TString filename5 = "/u/user/changhui/work/BICprototype/rootfiles/4by4_e-_1_GeV_final/root/R4by4_e-_1_GeV_final_10M" ;
  TString filename6 = "/u/user/changhui/work/BICprototype/rootfiles/4by4_e-_500_MeV_final/root/R4by4_e-_500_MeV_final_10M" ;
  int file1=0; int file2=0; int file3=1; int file4=0; int file5=0; int file6=0; // File loop on/off

  float inE = std::stof(argv[1]);
  float low = std::stof(argv[2]);
  float high = std::stof(argv[3]);

  int e1 = inE*5;int e2 = inE*4; int e3 = inE*3 ; int e4 = inE*2 ; int e5 = inE; int e6 = inE/2; // Beam E of each file
  int nLayers = 22;  int nFibers = 24;  int nRows = 3;  int nColumns = 5; // module configuration


  gStyle->SetOptFit(1);

  RootInterface<DRsimInterface::DRsimEventData>* drInterface = new RootInterface<DRsimInterface::DRsimEventData>(std::string(filename1 + ".root"), true);
  drInterface->set("DRsim","DRsimEventData");

  TH1F* tEdep = new TH1F("totEdep",";MeV;Evt",80,low*1000.,high*800.);
  tEdep->Sumw2(); tEdep->SetLineColor(kRed); tEdep->SetLineWidth(2);
  TH1F* tEdep4 = new TH1F("totEdep4",";MeV;Evt",80,low*1000.,high*800.);
  tEdep4->Sumw2(); tEdep4->SetLineColor(1); tEdep4->SetLineWidth(2);
  TH1F* tEdep5 = new TH1F("totEdep5",";MeV;Evt",80,low*1000.,high*800.);
  tEdep5->Sumw2(); tEdep5->SetLineColor(5); tEdep5->SetLineWidth(2);
  TH1F* tEdep2 = new TH1F("totEdep2",";MeV;Evt",80,low*1000.,high*800.);
  tEdep2->Sumw2(); tEdep2->SetLineColor(kBlue); tEdep2->SetLineWidth(2);
  TH1F* tEdep3 = new TH1F("totEdep3",";MeV;Evt",80,low*1000.,high*800.);
  tEdep3->Sumw2(); tEdep3->SetLineColor(kGreen); tEdep3->SetLineWidth(2);
  TH1F* tEdep6 = new TH1F("totEdep6",";MeV;Evt",80,low*1000.,high*800.);
  tEdep6->Sumw2(); tEdep6->SetLineColor(6); tEdep6->SetLineWidth(2);
 
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
  TH1F* tHit_S6 = new TH1F("Hit_S6",";# of p.e.;Evt",150,12000*low*2,9000*high*2);
  tHit_S6->Sumw2(); tHit_S6->SetLineColor(6); tHit_S6->SetLineWidth(2);

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
  TH1F* THits6 = new TH1F("THits6",";ModuleN; # of p.e.",20,0,20);
  THits6->Sumw2(); THits5->SetLineColor(6); THits5->SetLineWidth(2);
 
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
  TH1F* tP_leak6 = new TH1F("Pleak6",";MeV;Evt",80,0.,1000.*high);
  tP_leak6->Sumw2(); tP_leak6->SetLineWidth(2); tP_leak6->SetLineColor(6);

  TH1F* Edep_ratio = new TH1F("Edep_ratio","; Beam E (MeV); Edep (%)",60,low*1000.,high*600.);
  Edep_ratio->Sumw2(); Edep_ratio->SetLineWidth(1); Edep_ratio->SetLineColor(2);
  TH1F* Edep_ratio2 = new TH1F("Edep_ratio2","; Beam E (MeV); Edep (%)",60,low*1000.,high*600.);
  Edep_ratio2->Sumw2(); Edep_ratio2->SetLineWidth(1); Edep_ratio2->SetLineColor(3);
  TH1F* Edep_ratio3 = new TH1F("Edep_ratio3","; Beam E (MeV); Edep (%)",60,low*1000.,high*600.);
  Edep_ratio3->Sumw2(); Edep_ratio3->SetLineWidth(1); Edep_ratio3->SetLineColor(1);
 
  float rat_E = 0;  float rat_E2 = 0; float rat_E3 = 0; float rat_E4 = 0; float rat_E5 = 0; float rat_E6=0;
  double xBins = nFibers * nColumns;  double xLower = 0;  double xUpper = nFibers * nColumns;  double yBins = nLayers * nRows;
  double yLower = 0;  double yUpper = nLayers * nRows;

 TH2D* t2DhitS = new TH2D("2D Hit S1", "", xBins, xLower, xUpper, yBins, yLower, yUpper); t2DhitS->Sumw2(); t2DhitS->SetStats(0);
 TH2D* t2DhitS2 = new TH2D("2D Hit S2", "", xBins, xLower, xUpper, yBins, yLower, yUpper); t2DhitS2->Sumw2(); t2DhitS2->SetStats(0);
 TH2D* t2DhitS3 = new TH2D("2D Hit S3", "", xBins, xLower, xUpper, yBins, yLower, yUpper); t2DhitS3->Sumw2(); t2DhitS3->SetStats(0);
 TH2D* t2DhitS4 = new TH2D("2D Hit S4", "", xBins, xLower, xUpper, yBins, yLower, yUpper); t2DhitS4->Sumw2(); t2DhitS4->SetStats(0);
 TH2D* t2DhitS5 = new TH2D("2D Hit S5", "", xBins, xLower, xUpper, yBins, yLower, yUpper); t2DhitS5->Sumw2(); t2DhitS5->SetStats(0);
 TH2D* t2DhitS6 = new TH2D("2D Hit S6", "", xBins, xLower, xUpper, yBins, yLower, yUpper); t2DhitS6->Sumw2(); t2DhitS6->SetStats(0);

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

          Eleak += leak.kE;
    }
      
     tP_leak->Fill(Pleak);

    int nHitS = 0; int nHits;
    for (auto tower = drEvt.towers.begin(); tower != drEvt.towers.end(); ++tower) {
      int moduleNum = tower->ModuleNum; nHits=0;
      for (auto sipm = tower->SiPMs.begin(); sipm != tower->SiPMs.end(); ++sipm) {   
        int plateNum = sipm->y; int fiberNum = 24 - sipm->x; 
        nHitS += sipm->count;
        nHits += sipm->count;   

        t2DhitS->Fill(nFibers*(moduleNum%nColumns)+fiberNum, nLayers*(moduleNum/nColumns)+plateNum, sipm->count);
 
     }//SiPM loop
THits -> Fill(moduleNum,nHits);
    }//tower loop
    tHit_S->Fill(nHitS);
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
      
        Eleak2 += leak.kE;
    }
    tP_leak2->Fill(Eleak2);

  int nHitS2 = 0; int nHits2;
    for (auto tower2 = drEvt.towers.begin(); tower2 != drEvt.towers.end(); ++tower2) {
      int moduleNum2 = tower2->ModuleNum; nHits2=0;
      for (auto sipm2 = tower2->SiPMs.begin(); sipm2 != tower2->SiPMs.end(); ++sipm2) {
        int plateNum2 = sipm2->y; int fiberNum2 = 24 - sipm2->x;

        nHitS2 += sipm2->count;
	nHits2 += sipm2->count;

     t2DhitS2->Fill(nFibers*(moduleNum2%nColumns)+fiberNum2, nLayers*(moduleNum2/nColumns)+plateNum2, sipm2->count);
      
    }//fiber loop
     THits2 -> Fill(moduleNum2,nHits2);
    }//tower loop

    tHit_S2->Fill(nHitS2);
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
}    
    tP_leak3->Fill(Eleak3);

  int nHitS3 = 0; int nHits3;;
    for (auto tower3 = drEvt.towers.begin(); tower3 != drEvt.towers.end(); ++tower3) {
      int moduleNum3 = tower3->ModuleNum; nHits3=0;
      for (auto sipm3 = tower3->SiPMs.begin(); sipm3 != tower3->SiPMs.end(); ++sipm3) {
         int plateNum = sipm3->y; int fiberNum = 24 - sipm3->x;

        nHitS3 += sipm3->count;
        nHits3 += sipm3->count;

 t2DhitS3->Fill(nFibers*(moduleNum3%nColumns)+fiberNum, nLayers*(moduleNum3/nColumns)+plateNum, sipm3->count);
} // fiber loop
   THits3 -> Fill(moduleNum3,nHits3);
    }// tower loop
    tHit_S3->Fill(nHitS3);

} // 3rd event loop
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
     
        Eleak4 += leak.kE;
    }
    tP_leak4->Fill(Eleak4);
    
  int nHitS4 = 0; int nHits4;

    for (auto tower4 = drEvt.towers.begin(); tower4 != drEvt.towers.end(); ++tower4) {
      int moduleNum4 = tower4->ModuleNum; nHits4=0;
      for (auto sipm4 = tower4->SiPMs.begin(); sipm4 != tower4->SiPMs.end(); ++sipm4) {
        int plateNum4 = sipm4->y; int fiberNum4 = 24 - sipm4->x;

        nHitS4 += sipm4->count;
        nHits4 += sipm4->count;

         t2DhitS4->Fill(nFibers*(moduleNum4%nColumns)+fiberNum4, nLayers*(moduleNum4/nColumns)+plateNum4, sipm4->count);

} // fiber loop
   THits4 -> Fill(moduleNum4,nHits4);
    }// tower loop
    tHit_S4->Fill(nHitS4);

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
        Eleak5 += leak.kE;
    }
    tP_leak5->Fill(Eleak5); 

  int nHitS5 = 0; int nHits5;
    int nHitS_left5 = 0; int nHitS_right5 = 0; nHits5=0;
    for (auto tower5 = drEvt.towers.begin(); tower5 != drEvt.towers.end(); ++tower5) {
      int moduleNum5 = tower5->ModuleNum;
      for (auto sipm5 = tower5->SiPMs.begin(); sipm5 != tower5->SiPMs.end(); ++sipm5) {
        int plateNum5 = sipm5->y; int fiberNum5 = 24- sipm5->x;

        nHitS5 += sipm5->count;
        nHits5 += sipm5->count;

         t2DhitS5->Fill(nFibers*(moduleNum5%nColumns)+fiberNum5, nLayers*(moduleNum5/nColumns)+plateNum5, sipm5->count);
}// fiber loop
   THits5 -> Fill(moduleNum5,nHits5);
    }// tower loop
    tHit_S5->Fill(nHitS5);
} // 5th event root
} 

if (file6==1) {
    RootInterface<DRsimInterface::DRsimEventData>* drInterface6 = new RootInterface<DRsimInterface::DRsimEventData>(std::string(filename6 + ".root"), true);
  drInterface6->set("DRsim","DRsimEventData");

  unsigned int entries6 = drInterface6->entries();
  while (drInterface6->numEvt() < entries6) {
    if (drInterface6->numEvt() % 1000 == 0) printf("6th analyzing %dth event ...\n", drInterface6->numEvt());

    DRsimInterface::DRsimEventData drEvt;
    drInterface6->read(drEvt);

    float Edep6 = 0.; float totE6 = 0.; float ratE6 = 0; float rat_E6 = 0; float PE = 0;
    for (auto edepItr6 = drEvt.Edeps.begin(); edepItr6 != drEvt.Edeps.end(); ++edepItr6) {
      auto edep6 = *edepItr6;
      Edep6 += edep6.Edep;
    }
    tEdep6->Fill(Edep6);

    float Pleak6 = 0.;    float Eleak6 = 0.;
    for (auto leak : drEvt.leaks) {
      TLorentzVector leak4vec;
      leak4vec.SetPxPyPzE(leak.px,leak.py,leak.pz,leak.E);

        Eleak6 += leak.kE;
    }
    tP_leak6->Fill(Eleak6);

  int nHitS6 = 0; int nHits6;
    int nHitS_left6 = 0; int nHitS_right6 = 0; nHits6=0;
    for (auto tower6 = drEvt.towers.begin(); tower6 != drEvt.towers.end(); ++tower6) {
      int moduleNum6 = tower6->ModuleNum;
      for (auto sipm6 = tower6->SiPMs.begin(); sipm6 != tower6->SiPMs.end(); ++sipm6) {
        int plateNum6 = sipm6->y; int fiberNum6 = 24- sipm6->x;

        nHitS6 += sipm6->count;
        nHits6 += sipm6->count;

         t2DhitS6->Fill(nFibers*(moduleNum6%nColumns)+fiberNum6, nLayers*(moduleNum6/nColumns)+plateNum6, sipm6->count);
}// fiber loop
   THits6 -> Fill(moduleNum6,nHits6);
    }// tower loop
    tHit_S6->Fill(nHitS6);
} // 5th event root
}

   rat_E3 = tEdep3->GetMean() / e3 ;
   Edep_ratio -> SetBinContent(30,rat_E3*100);
   Edep_ratio -> SetBinError(30, tEdep3->GetStdDev()*100/e3);
   rat_E4 = tEdep4->GetMean() / e4 ;
   Edep_ratio -> SetBinContent(20,rat_E4*100);
   Edep_ratio -> SetBinError(20, tEdep4->GetStdDev()*100/e4);
   rat_E5 = tEdep5->GetMean() / e5 ;
   Edep_ratio -> SetBinContent(10,rat_E5*100);
   Edep_ratio -> SetBinError(10, tEdep5->GetStdDev()*100/e5);
   rat_E6 = tEdep6->GetMean() / e6 ;
   Edep_ratio -> SetBinContent(5,rat_E6*100);
   Edep_ratio -> SetBinError(5, tEdep6->GetStdDev()*100/e6);

    std::cout << "3GeV BinCon = " << rat_E3*100 << std::endl;
    std::cout << "3GeV BinErr = " << tEdep3->GetStdDev()*100/e3  << std::endl;
    std::cout << "2GeV BinCon = " << rat_E4*100 << std::endl;
    std::cout << "2GeV BinErr = " << tEdep4->GetStdDev()*100/e4  << std::endl;
    std::cout << "1GeV BinCon = " << rat_E5*100 << std::endl;
    std::cout << "1GeV BinErr = " << tEdep5->GetStdDev()*100/e5  << std::endl;
    std::cout << "0.5GeV BinCon = " << rat_E6*100 << std::endl;
    std::cout << "0.5GeV BinCon = " << tEdep6->GetStdDev()*100/e6 << std::endl;

   Edep_ratio2 -> SetBinContent(50, 86.87);
   Edep_ratio2 -> SetBinError(50, 5.40);
   Edep_ratio2 -> SetBinContent(40, 87.56);
   Edep_ratio2 -> SetBinError(40, 5.29);
   Edep_ratio2 -> SetBinContent(30, 88.54);
   Edep_ratio2 -> SetBinError(30, 4.91);
   Edep_ratio2 -> SetBinContent(20,89.69);
   Edep_ratio2 -> SetBinError(20, 4.64);
   Edep_ratio2 -> SetBinContent(10, 91.30);
   Edep_ratio2 -> SetBinError(10, 4.34);
   Edep_ratio2 -> SetBinContent(5,92.49);
   Edep_ratio2 -> SetBinError(5, 4.42);

   Edep_ratio3 -> SetBinContent(30,92.24);
   Edep_ratio3 -> SetBinError(30, 1.70);

  TCanvas* c = new TCanvas("c","");
  tEdep->Draw("Hist");  tEdep2->Draw("Hist same"); tEdep3->Draw("Hist same");  tEdep4->Draw("Hist same"); tEdep5->Draw("Hist same"); c->SaveAs(filename1+"compare_Edep_e.png");   

  c->SetLogy(1);
  tP_leak->Draw("Hist"); tP_leak2->Draw("Hist same");tP_leak3->Draw("Hist same"); tP_leak4->Draw("Hist same"); tP_leak5->Draw("Hist same"); 
c->SaveAs(filename1+"compare_Pleak_e.png");

  c->SetLogy(0);
  THits->Draw("Hist"); THits2->Draw("Hist same"); THits3->Draw("Hist same");  THits4->Draw("Hist same"); THits5->Draw("Hist same"); 
  THits->SetStats(0);
 c->SaveAs(filename1+"compare_towerE_e.png");

   tHit_S3->Draw("Hist"); tHit_S2->Draw("Hist same");tHit_S5->Draw("Hist same"); tHit_S4->Draw("Hist same"); tHit_S->Draw("Hist same"); c->SaveAs(filename1 +"compare_nHitpEventS_e.png");
 
  t2DhitS->Draw("COLZ"); c->SaveAs(filename1+"_n2DHitS.png");
  t2DhitS2->Draw("COLZ"); c->SaveAs(filename2+"_n2DHitS.png");
  t2DhitS3->Draw("COLZ"); c->SaveAs(filename3+"_n2DHitS.png");  
  t2DhitS4->Draw("COLZ"); c->SaveAs(filename4+"_n2DHitS.png");
  t2DhitS5->Draw("COLZ"); c->SaveAs(filename5+"_n2DHitS.png");

Edep_ratio->Draw("p E1"); Edep_ratio2->Draw("p E1 same ");Edep_ratio3->Draw("p E1 same"); 
gStyle->SetErrorX(0); Edep_ratio -> SetStats(0);  
Edep_ratio->SetMarkerSize(1.3); Edep_ratio->SetMarkerStyle(20); Edep_ratio->GetYaxis()->SetRangeUser(60,100); Edep_ratio->SetMarkerColor(2);
Edep_ratio2->SetMarkerSize(1.3); Edep_ratio2->SetMarkerStyle(20); Edep_ratio2->SetMarkerColor(3); 
Edep_ratio3->SetMarkerSize(1.3); Edep_ratio3->SetMarkerStyle(20); Edep_ratio3->SetMarkerColor(1);
c->SaveAs(filename1+"compare_Edep_ratio_e.png");

/*
  for ( int Bin = 0 ; Bin < 16 ; ++Bin) {
 std::cout << "Nhits " << Bin << " : "  << THits3->GetBinContent(Bin) <<std::endl;
}*/// Get bin counts
}

