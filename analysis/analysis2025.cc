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
  TString filename1 = "/u/user/changhui/BICprototype2025/rootfiles/3by8_e-_3_GeV_2025_slow/root/R3by8_e-_3_GeV_2025_slow_10M";
//  TString filename1 = "/u/user/changhui/BICprototype2025/rootfiles/3by5_e-_500_MeV_final_old/root/R3by5_e-_500_MeV_final_old_10M";
  TString filename3 = "/u/user/changhui/work/BICtest/rootfiles/3by5_e-_5_GeV_QE/root/R3by5_e-_5_GeV_QE_10M" ;
  TString filename4 = "/u/user/changhui/work/BICprototype/rootfiles/3by5_e-_2_GeV_final/root/R3by5_e-_2_GeV_final_10M" ;
  TString filename2 = "/u/user/changhui/BICprototype2025/rootfiles/3by8_e-_3_GeV_2025_fast/root/R3by8_e-_3_GeV_2025_fast_10M";
  TString filename5 = "/u/user/changhui/work/BICprototype/rootfiles/3by5_e-_4_GeV_final/root/R3by5_e-_4_GeV_final_10M" ;
  TString filename6 = "/u/user/changhui/work/BICprototype/rootfiles/3by5_e-_5_GeV_final/root/R3by5_e-_5_GeV_final_10M" ;

//  TString filename = "/u/user/changhui/work/BICtest/build/analysis";

  int file1=1; int file2=1; int file3=0; int file4=0; int file5=0; int file6=0; // File loop on/off

  float inE = std::stof(argv[1]);
  float low = std::stof(argv[2]);
  float high = std::stof(argv[3]);

  int e1 = inE; int e2 = inE; int e3 = inE*2 ; int e4 = inE*3 ; int e5 = inE*4; int e6 = inE*5; // Beam E of each file
  int nLayers = 22;  int nFibers = 24;  int nRows = 3;  int nColumns = 5; // module configuration

  gStyle->SetOptFit(1);

  RootInterface<DRsimInterface::DRsimEventData>* drInterface = new RootInterface<DRsimInterface::DRsimEventData>(std::string(filename1 + ".root"), true);
  drInterface->set("DRsim","DRsimEventData");

  TH1F* tEdep = new TH1F("totEdep",";MeV;Evt",100,low*80.,high*130.);
  tEdep->Sumw2(); tEdep->SetLineColor(2); tEdep->SetLineWidth(2);
  TH1F* tEdep4 = new TH1F("totEdep4",";MeV;Evt",80,low*1000.,high*200.);
  tEdep4->Sumw2(); tEdep4->SetLineColor(1); tEdep4->SetLineWidth(2);
  TH1F* tEdep5 = new TH1F("totEdep5",";MeV;Evt",80,low*1000.,high*800.);
  tEdep5->Sumw2(); tEdep5->SetLineColor(5); tEdep5->SetLineWidth(2);
  TH1F* tEdep2 = new TH1F("totEdep2",";MeV;Evt",100,low*80.,high*130.);
  tEdep2->Sumw2(); tEdep2->SetLineColor(kBlue); tEdep2->SetLineWidth(2);
  TH1F* tEdep3 = new TH1F("totEdep3",";MeV;Evt",80,low*1000.,high*800.);
  tEdep3->Sumw2(); tEdep3->SetLineColor(kGreen); tEdep3->SetLineWidth(2);
  TH1F* tEdep6 = new TH1F("totEdep6",";MeV;Evt",80,low*1000.,high*800.);
  tEdep6->Sumw2(); tEdep6->SetLineColor(6); tEdep6->SetLineWidth(2);
 

  TH1F* tHit_S = new TH1F("Hit_S",";# of p.e.;Evt",150,5000*low,5000*high);
  tHit_S->Sumw2(); tHit_S->SetLineColor(2); tHit_S->SetLineWidth(2);
  TH1F* tHit_S2 = new TH1F("Hit_S2",";# of p.e.;Evt",150,5000*low,5000*high);
  tHit_S2->Sumw2(); tHit_S2->SetLineColor(kBlue); tHit_S2->SetLineWidth(2);
  TH1F* tHit_S3 = new TH1F("Hit_S3",";# of p.e.;Evt",150,12000*low*2,100*high*2);
  tHit_S3->Sumw2(); tHit_S3->SetLineColor(kGreen); tHit_S3->SetLineWidth(2);
  TH1F* tHit_S4 = new TH1F("Hit_S4",";# of p.e.;Evt",150,12000*low*2,600*high*2);
  tHit_S4->Sumw2(); tHit_S4->SetLineColor(1); tHit_S4->SetLineWidth(2);
  TH1F* tHit_S5 = new TH1F("Hit_S5",";# of p.e.;Evt",150,12000*low*2,800*high*2);
  tHit_S5->Sumw2(); tHit_S5->SetLineColor(5); tHit_S5->SetLineWidth(2);
  TH1F* tHit_S6 = new TH1F("Hit_S6",";# of p.e.;Evt",150,12000*low*2,1000*high*2);
  tHit_S6->Sumw2(); tHit_S6->SetLineColor(6); tHit_S6->SetLineWidth(2);

  TH1F* Edep_M[41];
for (int i =1; i < 41; i++){
  Edep_M[i] = new TH1F(Form("Edep_M%d",i),";Energy(MeV);Evt",inE*0.12,0,inE*1.2);
  Edep_M[i]->Sumw2(); Edep_M[i]->SetLineColor(4); Edep_M[i]->SetLineWidth(2);}

  TH1F* Nhits_M_L[41]; TH1F* Nhits_M_R[41];
for (int j =1; j < 41; j++){
  Nhits_M_L[j] = new TH1F(Form("Nhits_L_M%d",j),";# of p.e.;Evt",0.2*inE,0,0.2*inE);
  Nhits_M_L[j]->Sumw2(); Nhits_M_L[j]->SetLineColor(3); Nhits_M_L[j]->SetLineWidth(2);
  Nhits_M_R[j] = new TH1F(Form("Nhits_R_M%d",j),";# of p.e.;Evt",0.2*inE,0,0.2*inE);
  Nhits_M_R[j]->Sumw2(); Nhits_M_R[j]->SetLineColor(3); Nhits_M_R[j]->SetLineWidth(2);
}

  TH1F* THits = new TH1F("THits",";ModuleN; # of p.e (%)",32,0,32);
  THits->Sumw2(); THits->SetLineColor(kRed); THits->SetLineWidth(2);
  TH1F* THits2 = new TH1F("THits2",";ModuleN; # of p.e.",32,0,32);
  THits2->Sumw2(); THits2->SetLineColor(kBlue); THits2->SetLineWidth(2);
  TH1F* THits3 = new TH1F("Nhit per Module",";ModuleN; # of p.e.",20,0,20);
  THits3->Sumw2(); THits3->SetLineColor(kGreen); THits3->SetLineWidth(2);
  TH1F* THits4 = new TH1F("THits4",";ModuleN; # of p.e.",20,0,20);
  THits4->Sumw2(); THits4->SetLineColor(1); THits4->SetLineWidth(2);
  TH1F* THits5 = new TH1F("THits5",";ModuleN; # of p.e.",20,0,20);
  THits5->Sumw2(); THits5->SetLineColor(5); THits5->SetLineWidth(2);
  TH1F* THits6 = new TH1F("THits6",";ModuleN; # of p.e.",20,0,20);
  THits6->Sumw2(); THits6->SetLineColor(6); THits6->SetLineWidth(2);
 
  TH1F* THitss = new TH1F("THits",";ModuleN; Edep(MeV)",32,0,32);
  THitss->Sumw2(); THitss->SetLineColor(kRed); THitss->SetLineWidth(2);
  TH1F* THitss2 = new TH1F("THits2",";ModuleN; # of p.e.",32,0,32);
  THitss2->Sumw2(); THitss2->SetLineColor(kBlue); THitss2->SetLineWidth(2);
  TH1F* THitss3 = new TH1F("Edep per Module",";ModuleN;Edep(MeV)/evt",20,0,20);
  THitss3->Sumw2(); THitss3->SetLineColor(kGreen); THitss3->SetLineWidth(2);
  TH1F* THitss4 = new TH1F("THits4",";ModuleN; # of p.e.",20,0,20);
  THitss4->Sumw2(); THitss4->SetLineColor(1); THitss4->SetLineWidth(2);
  TH1F* THitss5 = new TH1F("THits5",";ModuleN; # of p.e.",20,0,20);
  THitss5->Sumw2(); THitss5->SetLineColor(5); THitss5->SetLineWidth(2);
  TH1F* THitss6 = new TH1F("THits6",";ModuleN; # of p.e.",20,0,20);
  THitss6->Sumw2(); THitss6->SetLineColor(6); THitss6->SetLineWidth(2);

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
  Edep_ratio3->Sumw2(); Edep_ratio3->SetLineWidth(1); Edep_ratio3->SetLineColor(4);
  TH1F* Edep_ratio4 = new TH1F("Edep_ratio4","; Beam E (MeV); Edep (%)",60,low*1000.,high*600.);
  Edep_ratio4->Sumw2(); Edep_ratio4->SetLineWidth(1); Edep_ratio4->SetLineColor(28);
  TH1F* Edep_ratio5 = new TH1F("Edep_ratio5","; Beam E (MeV); Edep (%)",60,low*1000.,high*600.);
  Edep_ratio5->Sumw2(); Edep_ratio5->SetLineWidth(1); Edep_ratio5->SetLineColor(1);
  TH1F* Edep_ratio6 = new TH1F("Edep_ratio6","; Beam E (MeV); Edep (%)",60,low*1000.,high*600.);
  Edep_ratio6->Sumw2(); Edep_ratio6->SetLineWidth(1); Edep_ratio6->SetLineColor(6);
  TH1F* Edep_ratio7 = new TH1F("Edep_ratio7","; Beam E (MeV); Edep (%)",60,low*1000.,high*600.);
  Edep_ratio7->Sumw2(); Edep_ratio7->SetLineWidth(1); Edep_ratio7->SetLineColor(7);

  float rat_E = 0;  float rat_E2 = 0; float rat_E3 = 0; float rat_E4 = 0; float rat_E5 = 0; float rat_E6=0;
  double xBins = nFibers * nColumns;  double xLower = 0;  double xUpper = nFibers * nColumns;  double yBins = nLayers * nRows;
  double yLower = 0;  double yUpper = nLayers * nRows;

 TH2D* t2DhitS = new TH2D("2D Hit S1", "", xBins, xLower, xUpper, yBins, yLower, yUpper); t2DhitS->Sumw2(); t2DhitS->SetStats(0);
 TH2D* t2DhitS2 = new TH2D("2D Hit S2", "", xBins, xLower, xUpper, yBins, yLower, yUpper); t2DhitS2->Sumw2(); t2DhitS2->SetStats(0);
 TH2D* t2DhitS3 = new TH2D("2D Hit S3", "", xBins, xLower, xUpper, yBins, yLower, yUpper); t2DhitS3->Sumw2(); t2DhitS3->SetStats(0);
 TH2D* t2DhitS4 = new TH2D("2D Hit S4", "", xBins, xLower, xUpper, yBins, yLower, yUpper); t2DhitS4->Sumw2(); t2DhitS4->SetStats(0);
 TH2D* t2DhitS5 = new TH2D("2D Hit S5", "", xBins, xLower, xUpper, yBins, yLower, yUpper); t2DhitS5->Sumw2(); t2DhitS5->SetStats(0);
 TH2D* t2DhitS6 = new TH2D("2D Hit S6", "", xBins, xLower, xUpper, yBins, yLower, yUpper); t2DhitS6->Sumw2(); t2DhitS6->SetStats(0);

float Edep_Num[40] = {0};float Edep_Num2[40] = {0};float Edep_Num3[40] = {0};float Edep_Num4[40] = {0}; float Edep_Numm[40] = {0};

if (file1==1) {

 TFile *file = new TFile("/u/user/changhui/BICprototype2025/build/analysis/3x5_1_GeV_old.root", "RECREATE");
 TTree *tree = new TTree("Ttree", "Edep and Nhits ");
// tree->SetBranchAddress("edep.Edep", &edep);
// tree->SetBranchAddress("edep.Module", &moduleN);

    Float_t edep1[32]={0};  // 24개의모듈 데이터를 저장
    Int_t moduleN_E[32]={0} ; Int_t moduleN_H[32]={0} ; Float_t Nhits_L[32]={0}; Float_t Nhits_R[32]={0}; 
    tree->Branch("edep1", edep1, "edep1[32]/F");
    tree->Branch("moduleN_E",moduleN_E, "moduleN_E[32]/I");
    tree->Branch("moduleN_H",moduleN_H, "moduleN_H[32]/I");
    tree->Branch("Nhits_L", Nhits_L, "Nhits_L[32]/F");
    tree->Branch("Nhits_R", Nhits_R, "Nhits_R[32]/F");

    float totHit = 0; float nHits[32]= {0};
  unsigned int entries = 10000; //drInterface->entries();
  while (drInterface->numEvt() < entries) {
    if (drInterface->numEvt() % 1000 == 0) printf("1st analyzing %dth event ...\n", drInterface->numEvt());

    DRsimInterface::DRsimEventData drEvt;
    drInterface->read(drEvt);
//    float edep1[25]={0};  // 24개의모듈 데이터를 저장
//    float Nhits_L[25]={0}; float Nhits_R[25]={0};  int moduleN_E[25]={0} ; int moduleN_H[25]={0};

    float Edep = 0.; float totE = 0.; float ratE = 0;  float Edep_Num[32] = {0}; //float Edep_Numm[24] = {0};

 float Nhits_Numm_L[32]={0}; float Nhits_Numm_R[32]={0};
    for (auto edepItr = drEvt.Edeps.begin();  edepItr != drEvt.Edeps.end(); ++edepItr) {
      auto edep = *edepItr;

    for (int i=0; i<32; i++){
      if (edep.ModuleNum == i){
        Edep_Num[i] += edep.Edep;
        Edep_Numm[i] += edep.Edep;
        
	} 
 //     else {Edep_Num[i] = 0;}
    }
      Edep += edep.Edep;
    }
    tEdep->Fill(Edep); 
float Eleak = 0;
    for (auto leak : drEvt.leaks) {
            Eleak += leak.kE;
}
tP_leak->Fill(Eleak);
for (int j=1; j<33; j++){
    Edep_M[j]->Fill(Edep_Num[j-1]) ;
    edep1[j-1] = Edep_Num[j-1];     
    moduleN_E[j-1] = j ; 
//std::cout<<"Edep_Numm : "<<Edep_Num[j-1]<< "edep : " << edep1[j-1] << "moduleN : "<< moduleN_E[j-1] <<  j <<std::endl;
} 
//tree->Fill();
    int nHitS = 0; int isLeft=0;
    for (auto tower = drEvt.towers.begin(); tower != drEvt.towers.end(); ++tower) {
      int moduleNum = tower->ModuleNum;// nHits=0; 
      for (auto sipm = tower->SiPMs.begin(); sipm != tower->SiPMs.end(); ++sipm) {   
        isLeft = sipm->isleft;
        int plateNum = sipm->y; int fiberNum = 24 - sipm->x; 
        nHitS += sipm->count;
        nHits[moduleNum] += sipm->count;   

  for (int i=0;i<33;i++) {
//     if (i == 7 || i == 15 || i == 23 || i == 31) continue;
   if (moduleNum == i) {//moduleN_H=i;
      if (isLeft==0) {Nhits_Numm_L[i] += sipm->count;// Nhits_L[i]+=sipm->count; 
}
      else {Nhits_Numm_R[i] += sipm->count;// Nhits_R[i] += sipm->count;
}
    }
    }
        t2DhitS->Fill(nFibers*(moduleNum%nColumns)+fiberNum, nLayers*(moduleNum/nColumns)+plateNum, sipm->count); 
     }//SiPM loop
//  THits -> Fill(moduleNum,nHits/entries);
    }//tower loop
//      std::cout << "nHits  : " << nHitS << std::endl;

tHit_S->Fill(nHitS/*/36.306*//*/112.78*/);  totHit += nHitS ; 
//std::cout << "nHit  : " <<  nHits[1] << std::endl;
/*
  for (int j=1; j<33; j++){
 //   if (j == 8 || j == 16 || j == 24 || j == 15 || j == 14|| j == 23|| j == 22|| j == 32|| j == 31|| j == 30) continue;
//    if (j == 8 || j == 16 || j == 24 || j == 32) continue;
    nHitS += Nhits_Numm_L[j-1]  ;
    nHitS += Nhits_Numm_R[j-1] ;
    Nhits_M_L[j]->Fill(Nhits_Numm_L[j-1]/140.82);
    Nhits_M_R[j]->Fill(Nhits_Numm_R[j-1]/140.82);
    Nhits_L[j-1]=Nhits_Numm_L[j-1];
    Nhits_R[j-1]=Nhits_Numm_R[j-1];
    moduleN_H[j-1]=j;
//std::cout<<"Edep_Numm : "<<Edep_Num[j-1]<< "edep : " << edep1[j] << "moduleN : "<< moduleN_E[j] <<  j <<std::endl;

  } tree->Fill();*/
//      tHit_S->Fill(nHitS/140.82);
  } // event loop
//    THits->Scale(1/totHit);
  for (int z=0;z<32;z++) {  THits -> Fill(z,nHits[z]*100/totHit); }


    rat_E = tEdep->GetMean() / e1 ; std::cout << "ratE  : " <<  rat_E << std::endl;
    Edep_ratio -> SetBinContent(30,rat_E*100);
    float a = Edep_ratio->GetBinContent(30);
    std::cout << "BinCon = " << a << std::endl;
    Edep_ratio -> SetBinError(30, tEdep->GetStdDev()*100/e1);
    std::cout << "BinErr = " << tEdep->GetStdDev()*100/e1  << std::endl;
/*
  for (int jjj=1; jjj<25; jjj++){

std::cout<<"Edep_Numm : "<<Edep_Num[jjj-1]<< "edep : " << edep1[jjj] << "moduleN : "<< moduleN_E[jjj]  <<std::endl;
}
*/
//tree->Write();
//file->Close();
 }

if (file2==1) {
// second file
  RootInterface<DRsimInterface::DRsimEventData>* drInterface2 = new RootInterface<DRsimInterface::DRsimEventData>(std::string(filename2 + ".root"), true);
  drInterface2->set("DRsim","DRsimEventData"); 

  int totHit2=0; float nHits2[32]={0};
  unsigned int entries2 = 10000;// drInterface2->entries();
  while (drInterface2->numEvt() < entries2) {
    if (drInterface2->numEvt() % 1000 == 0) printf("2nd analyzing %dth event ...\n", drInterface2->numEvt());

    DRsimInterface::DRsimEventData drEvt;
    drInterface2->read(drEvt);

    float Edep2 = 0.; float totE2 = 0.; float ratE2 = 0; float rat_E2 = 0;
    for (auto edepItr2 = drEvt.Edeps.begin(); edepItr2 != drEvt.Edeps.end(); ++edepItr2) {
      auto edep2 = *edepItr2;

for (int i2=0; i2<32; i2++){
      if (edep2.ModuleNum == i2){
        Edep_Num2[i2] += edep2.Edep;
      }
    } 
   Edep2 += edep2.Edep;
    }
    tEdep2->Fill(Edep2);


  int nHitS2 = 0; 
    for (auto tower2 = drEvt.towers.begin(); tower2 != drEvt.towers.end(); ++tower2) {
      int moduleNum2 = tower2->ModuleNum; 
      for (auto sipm2 = tower2->SiPMs.begin(); sipm2 != tower2->SiPMs.end(); ++sipm2) {
        int plateNum2 = sipm2->y; int fiberNum2 = 24 - sipm2->x;

        nHitS2 += sipm2->count;
	nHits2[moduleNum2] += sipm2->count;
	

     t2DhitS2->Fill(nFibers*(moduleNum2%nColumns)+fiberNum2, nLayers*(moduleNum2/nColumns)+plateNum2, sipm2->count);
      
    }//fiber loop
//     THits2 -> Fill(moduleNum2,nHits2/entries2);
    }//tower loop
//std::cout << "Nhits :" << nHitS2 << std::endl;
    tHit_S2->Fill(nHitS2); totHit2 += nHitS2 ;
} // 2nd event root
    for (int zz=0;zz<32;zz++) {  THits2 -> Fill(zz,nHits2[zz]*100/totHit2); }
  
    rat_E2 = tEdep2->GetMean() / e2 ;
    Edep_ratio2-> SetBinContent(30,rat_E2*100);
    Edep_ratio2 -> SetBinError(30, tEdep2->GetStdDev()*100/e2);
//    std::cout << "4GeV BinCon = " << rat_E2*100 << std::endl;
//    std::cout << "4GeV BinErr = " << tEdep2->GetStdDev()*100/e2  << std::endl;

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
 
    for (int i3=0; i3<16; i3++){
      if (edep3.ModuleNum == i3){
        Edep_Num3[i3] += edep3.Edep/entries3;
      }
    }


     Edep3 += edep3.Edep;
    }
    tEdep3->Fill(Edep3);

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

for (int i4=0; i4<16; i4++){
      if (edep4.ModuleNum == i4){
        Edep_Num4[i4] += edep4.Edep;
      }
    }
      Edep4 += edep4.Edep;
    }
    tEdep4->Fill(Edep4);
   
  int nHitS4 = 0; int nHits4;

    for (auto tower4 = drEvt.towers.begin(); tower4 != drEvt.towers.end(); ++tower4) {
      int moduleNum4 = tower4->ModuleNum; nHits4=0;
      for (auto sipm4 = tower4->SiPMs.begin(); sipm4 != tower4->SiPMs.end(); ++sipm4) {
        int plateNum4 = sipm4->y; int fiberNum4 = 24 - sipm4->x;

        nHitS4 += sipm4->count;
        nHits4 += sipm4->count/entries4;

         t2DhitS4->Fill(nFibers*(moduleNum4%nColumns)+fiberNum4, nLayers*(moduleNum4/nColumns)+plateNum4, sipm4->count);

} // fiber loop
   THits4 -> Fill(moduleNum4,nHits4);
    }// tower loop
    tHit_S4->Fill(nHitS4/18.9);

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
    tHit_S5->Fill(nHitS5/18.9);
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
    tHit_S6->Fill(nHitS6/18.9);
} // 5th event root
}
 for (int ii=1; ii<33; ii++){
     THitss->SetBinContent(ii,Edep_Numm[ii-1]/10000);
    }
// for (int ii2=1; ii2<17; ii2++){
///     THitss2->SetBinContent(ii2,Edep_Num2[ii2-1]/10000);
//    }
// for (int ii3=1; ii3<17; ii3++){
//     THitss3->SetBinContent(ii3,Edep_Num3[ii3-1]/10000);
//    }
// for (int ii4=1; ii4<17; ii4++){
//     THitss4->SetBinContent(ii4,Edep_Num4[ii4-1]/10000);
//    }
 

  std::cout << "Module total1 = " << THitss->Integral() << std::endl;

   rat_E3 = tEdep3->GetMean() / e3 ;
   Edep_ratio3 -> SetBinContent(30,rat_E3*100);
   Edep_ratio3 -> SetBinError(30, tEdep3->GetStdDev()*100/e3);
   rat_E4 = tEdep4->GetMean() / e4 ;
   Edep_ratio4 -> SetBinContent(30,rat_E4*100);
   Edep_ratio4 -> SetBinError(30, tEdep4->GetStdDev()*100/e4);
   rat_E5 = tEdep5->GetMean() / e5 ;
   Edep_ratio -> SetBinContent(10,rat_E5*100);
   Edep_ratio -> SetBinError(10, tEdep5->GetStdDev()*100/e5);
   rat_E6 = tEdep6->GetMean() / e6 ;
   Edep_ratio7 -> SetBinContent(30,rat_E6*100);
   Edep_ratio7 -> SetBinError(30, tEdep6->GetStdDev()*100/e6);

    std::cout << "3GeV BinCon = " << rat_E3*100 << std::endl;
    std::cout << "3GeV BinErr = " << tEdep3->GetStdDev()*100/e3  << std::endl;
    std::cout << "2GeV BinCon = " << rat_E4*100 << std::endl;
    std::cout << "2GeV BinErr = " << tEdep4->GetStdDev()*100/e4  << std::endl;
    std::cout << "1GeV BinCon = " << rat_E5*100 << std::endl;
    std::cout << "1GeV BinErr = " << tEdep5->GetStdDev()*100/e5  << std::endl;
    std::cout << "0.5GeV BinCon = " << rat_E6*100 << std::endl;
    std::cout << "0.5GeV BinCon = " << tEdep6->GetStdDev()*100/e6 << std::endl;
/*
   Edep_ratio2 -> SetBinContent(50, 86.87);
   Edep_ratio2 -> SetBinError(50, 5.40);
   Edep_ratio2 -> SetBinContent(40, 87.56);
   Edep_ratio2 -> SetBinError(40, 5.29);*/
/*   Edep_ratio5 -> SetBinContent(30, 88.54);
   Edep_ratio5 -> SetBinError(30, 4.91);
   Edep_ratio2 -> SetBinContent(20,89.69);
   Edep_ratio2 -> SetBinError(20, 4.64);
   Edep_ratio2 -> SetBinContent(10, 91.30);
   Edep_ratio2 -> SetBinError(10, 4.34);
   Edep_ratio2 -> SetBinContent(5,92.49);
   Edep_ratio2 -> SetBinError(5, 4.42);

   Edep_ratio6 -> SetBinContent(30,92.24);
   Edep_ratio6 -> SetBinError(30, 1.70);
*/

    TF1 *gaussFit_E = new TF1("gaussFit_E", "gaus", 0, 5000);

  TCanvas* c = new TCanvas("c","");
  tEdep->Fit(gaussFit_E);  tEdep->SetStats(1);
  tEdep->Draw("Hist"); 
//  gaussFit_E->Draw("same");
     gStyle->SetOptFit(1);

 tEdep2->Draw("p same");  tEdep2->SetMarkerStyle(24); tEdep2->SetMarkerColor(kBlue);// tEdep3->Draw("Hist same");  tEdep4->Draw("Hist same"); tEdep5->Draw("Hist same"); 
 tEdep->SetStats(1);
c->SaveAs(filename1+"compare_EdepG_e.png");   

    TF1 *gaussFit_S = new TF1("gaussFit_S", "gaus", 0, 200000);
/*
for (int j=1; j<33; j++){ c->cd(j);
//    Edep_M[j]->Fit(gaussFit_E);
    Edep_M[j]->Draw("Hist");   
//    gaussFit_E->Draw("same");  
//    gStyle->SetOptFit(1);
 c->SaveAs(Form(filename1+"EdepG_M%d.png",j));  }
for (int jj=1; jj<33; jj++){  c->cd(jj);
//    Nhits_M[jj]->Fit(gaussFit_S);
    Nhits_M_L[jj]->Draw("Hist"); 
    Nhits_M_R[jj]->Draw("Hist");
//    gaussFit_S->Draw("same");
//    gStyle->SetOptFit(1);
 c->SaveAs(Form(filename1+"NhitsG_M%d.png",jj)); }
*/
    /*
  c->SetLogy(1);
  tP_leak->Draw("Hist"); tP_leak2->Draw("Hist same");tP_leak3->Draw("Hist same"); tP_leak4->Draw("Hist same"); tP_leak5->Draw("Hist same"); 
c->SaveAs(filename1+"compare_Pleak_e.png");
*/
  c->SetLogy(0);
  THits->Draw("Hist"); THits2->Draw("Hist p same"); THits2->SetMarkerStyle(24); THits2->SetMarkerColor(kBlue);
  /*THits3->Draw("Hist same");  THits4->Draw("Hist same"); THits5->Draw("Hist same"); */
  THits->SetStats(0);
 c->SaveAs(filename1+"compare_towerE_e.png");
  
  THitss->Draw("Hist");
   

/* THitss2->Draw("Hist same"); THitss3->Draw("Hist same");  THitss4->Draw("Hist same"); THitss5->Draw("Hist same");*/
//  THitss->SetStats(0);
//

  float Tot = THitss->Integral();
 for (int o=1; o<33; o++){ float val = THitss->GetBinContent(o); 
 std::cout << Form("Mod %d :  ",o) << val *100/Tot << "%,   " << val  << "MeV" <<  std::endl;   }
// c->SaveAs(filename1+"compare_tower_Edep.png");
 std::cout << "Total Edep, (%) : " << Tot <<  " MeV, "  << Tot*100/inE << "(%)" << std::endl;
  tHit_S->Fit(gaussFit_S);/* tHit_S2->Fit(gaussFit_S);tHit_S3->Fit(gaussFit_S);tHit_S4->Fit(gaussFit_S);tHit_S5->Fit(gaussFit_S);tHit_S6->Fit(gaussFit_S); */
   tHit_S->SetStats(1);
 //  tHit_S->GetYaxis()->SetRangeUser(0,100);
 
 tHit_S->SetMarkerColor(kRed);  tHit_S->SetMarkerStyle(7);  tHit_S->Draw("Hist"); tHit_S2->SetMarkerColor(kBlue);tHit_S3->SetMarkerColor(1);  tHit_S2->SetMarkerStyle(7); tHit_S3->SetMarkerStyle(7);
//    gaussFit_S->Draw("same");
    gStyle->SetOptFit(1);

  tHit_S2->Draw("Hist p same" ); tHit_S2->SetMarkerStyle(24); //tHit_S2->SetMarkerStyle(kBlue);//gaussFit_S->Draw("same");
//  tHit_S3->Draw("Hist p same");//gaussFit_S->Draw("same"); tHit_S4->Draw("Hist");gaussFit_S->Draw("same"); tHit_S5->Draw("Hist");gaussFit_S->Draw("same");tHit_S6->Draw("Hist");gaussFit_S->Draw("same");

c->SaveAs(filename1 +"compare_nHitpEventSG_cal.png");

    t2DhitS->Draw("COLZ"); c->SaveAs(filename1+"_n2DHitS.png");
//  t2DhitS2->Draw("COLZ"); c->SaveAs(filename2+"_n2DHitS.png");
//  t2DhitS3->Draw("COLZ"); c->SaveAs(filename3+"_n2DHitS.png");  
//  t2DhitS4->Draw("COLZ"); c->SaveAs(filename4+"_n2DHitS.png");
//  t2DhitS5->Draw("COLZ"); c->SaveAs(filename5+"_n2DHitS.png");

Edep_ratio->Draw("p E1"); Edep_ratio2->Draw("p E1 same ");/*Edep_ratio3->Draw("p E1 same"); Edep_ratio4->Draw("p E1 same"); Edep_ratio5->Draw("p E1 same ");Edep_ratio6->Draw("p E1 same"); Edep_ratio7->Draw("p E1 same");*/
gStyle->SetErrorX(0); Edep_ratio -> SetStats(0);  
Edep_ratio->SetMarkerSize(1.3); Edep_ratio->SetMarkerStyle(20); Edep_ratio->GetYaxis()->SetRangeUser(70,100); Edep_ratio->SetMarkerColor(2);
Edep_ratio2->SetMarkerSize(1.3); Edep_ratio2->SetMarkerStyle(20); Edep_ratio2->SetMarkerColor(3); 
Edep_ratio3->SetMarkerSize(1.3); Edep_ratio3->SetMarkerStyle(20); Edep_ratio3->SetMarkerColor(4);
Edep_ratio4->SetMarkerSize(1.3); Edep_ratio4->SetMarkerStyle(20); Edep_ratio4->SetMarkerColor(28);
Edep_ratio5->SetMarkerSize(1.3); Edep_ratio5->SetMarkerStyle(20); Edep_ratio5->SetMarkerColor(1);
Edep_ratio6->SetMarkerSize(1.3); Edep_ratio6->SetMarkerStyle(20); Edep_ratio6->SetMarkerColor(6);
Edep_ratio7->SetMarkerSize(1.3); Edep_ratio7->SetMarkerStyle(20); Edep_ratio7->SetMarkerColor(7);

/* 
   TFile *outputFile1 = new TFile("/u/user/changhui/work/BICprototype2/build/analysis/resol.root", "RECREATE");    
tHit_S->Write();tHit_S2->Write();tHit_S3->Write();tHit_S4->Write();tHit_S5->Write();tHit_S6->Write();
outputFile1->Close();   
*/
/*
   TFile *outputFile2 = new TFile("/u/user/changhui/work/BICKEK/build/analysis/3x5_1GeV_old.root", "RECREATE");
for (int jjj=1; jjj<33; jjj++){ 
    Edep_M[jjj]->Write();
}
for (int jjjj=1; jjjj<33; jjjj++){
    Nhits_M_L[jjjj]->Write(); Nhits_M_R[jjjj]->Write();
}
outputFile2->Close();
*/
/*
 TFile *file = new TFile("/u/user/changhui/work/BICprototype2/build/analysis/3x8_3GeV_tree.root", "RECREATE");
 TTree *tree = new TTree("hist_tree", "Tree with multiple histograms");
//    TH1F *hist = nullptr;     TH1F *hist2 = nullptr;     TH1F *hist3 = nullptr;
    float entries[360];
    float Edep[360];  // 100개의 bin 데이터를 저장
    int moduleN = 0;
    tree->Branch("Edep", Edep, "Edep[360]/F");
    tree->Branch("moduleN",&moduleN);
    tree->Branch("entries",entries, "entries[360]/F"); 
    for (int jjj = 1; jjj < 25; ++jjj) {
             entries[360]=0; Edep[360]=0;
       for (int bin = 1; bin < 361 ; ++bin) {
            entries[bin - 1] = Edep_M[jjj]->GetBinContent(bin);}
       for (int bin2 = 1; bin2 < 361 ; ++bin2) {
            Edep[bin2 - 1] = Edep_M[jjj]->GetBinCenter(bin2);}
            moduleN=jjj;
     tree -> Fill();
    }
    tree->Write();
    file->Close();
*/


c->SaveAs(filename1+"compare_Edep_ratio_e.png");
/*
  TCanvas* c1 = new TCanvas("c1","");
c1->SetCanvasSize(1200,480);
c1->Divide(5,3);

for (int j=1; j<16; j++){ c1->cd(j);// c1->Clear();// c1->SetMargin(0.01, 0.01, 0.01, 0.01);
        Edep_M[j]->Draw("Hist");
	gStyle -> SetOptFit(0); 
	TPad*pad=(TPad*)gPad;
        pad->SetTopMargin(0.01);
        pad->SetBottomMargin(0.12);
        pad->SetLeftMargin(0.15);
        pad->SetRightMargin(0.01);
        
        Edep_M[j]->GetXaxis()->SetTitleSize(0.06);
        Edep_M[j]->GetYaxis()->SetTitleSize(0.05);

//gPad->Update();

    TPaveStats *stats = (TPaveStats*)Edep_M[j]->GetListOfFunctions()->FindObject("stats");
    if (stats) { 
	double leftMargin = pad->GetLeftMargin();
        double rightMargin = pad->GetRightMargin();

        stats->SetTextSize(0.08);
stats->SetX1NDC(0.65);stats->SetX2NDC(1-rightMargin);stats->SetY1NDC(0.55);stats->SetY2NDC(0.99);
 }}

c1->SaveAs(filename1+"Edep_Allmod.png");
*/
/*
for (int jj=1; jj<16; jj++){  c1->cd(jj);
//    Nhits_M[jj]->Draw("Hist");
//    gStyle->SetOptFit(0);

TPad*pad1=(TPad*)gPad;
     pad1->SetTopMargin(0.01);
        pad1->SetBottomMargin(0.1);
        pad1->SetLeftMargin(0.15);
        pad1->SetRightMargin(0.01);

    TPaveStats *stats1 = (TPaveStats*)Nhits_M[jj]->GetListOfFunctions()->FindObject("stats");
    if (stats1) {
        stats1->SetTextSize(0.08);
stats1->SetX1NDC(0.55);stats1->SetX2NDC(0.99);stats1->SetY1NDC(0.55);stats1->SetY2NDC(0.99);

} }c1->SaveAs(filename1+"Nhits_Allmod.png");
*/

/*
  for ( int Bin = 0 ; Bin < 16 ; ++Bin) {
 std::cout << "Nhits " << Bin << " : "  << THits3->GetBinContent(Bin) <<std::endl;
}*/// Get bin counts
}

