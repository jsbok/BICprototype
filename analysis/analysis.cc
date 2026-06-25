#include "RootInterface.h"
#include "koBICInterface.h"
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
  TString filename1 = "/u/user/changhui/koBIC2026/BICprototype/rootfiles/3SF1B_e-_1_GeV_2026/root/R3SF1B_e-_1_GeV_2026_500";
//  TString filename1 = "/u/user/changhui/koBIC2025/BICprototype2025/rootfiles/4by8_e-_7_GeV_2025/root/R4by8_epi_7_GeV_2025_20M";
  TString filename2 = "/u/user/changhui/BICprototype2025/rootfiles/3by8_e-_1_GeV2024_slow/root/R3by8_e-_1_GeV2024_slow_10M";

  int file1=1; int file2=0; // File loop on/off

  float inE = std::stof(argv[1]);
  float low = std::stof(argv[2]);
  float high = std::stof(argv[3]);

  int e1 = inE; int e2 = inE; // Beam E of each file
  int nLayers = 22;  int nFibers = 24;  int nRows = 3;  int nColumns = 5; // module configuration

  gStyle->SetOptFit(1);

  RootInterface<koBICInterface::koBICEventData>* drInterface = new RootInterface<koBICInterface::koBICEventData>(std::string(filename1 + ".root"), true);
  drInterface->set("koBIC","koBICEventData");

  TH1F* tEdepC = new TH1F("S.F.",";Sampling Fraction; Evt",100,0.,0.2);
  tEdepC->Sumw2(); tEdepC->SetLineColor(2); tEdepC->SetLineWidth(2);


  TH1F* tEdep = new TH1F("totEdep",";MeV;Evt",100,low*80.,high*200.);
  tEdep->Sumw2(); tEdep->SetLineColor(2); tEdep->SetLineWidth(2);
  TH1F* tEdep2 = new TH1F("totEdep2",";MeV;Evt",100,low*80.,high*130.);
  tEdep2->Sumw2(); tEdep2->SetLineColor(kBlue); tEdep2->SetLineWidth(2);

  TH1F* tHit_S = new TH1F("Hit_S","; reconstructed E (MeV);Evt",150, 0,60000);
  tHit_S->Sumw2(); tHit_S->SetLineColor(2); tHit_S->SetLineWidth(2);
  TH1F* tHit_S2 = new TH1F("Hit_S2","; reconstructed E (MeV);Evt",300,0,21000);
  tHit_S2->Sumw2(); tHit_S2->SetLineColor(kBlue); tHit_S2->SetLineWidth(2);

  TH1F* Edep_M[51];
for (int i =1; i < 51; i++){
  Edep_M[i] = new TH1F(Form("Edep_M%d",i),";Energy(MeV); Evt",100,0,0.2);
  Edep_M[i]->Sumw2(); Edep_M[i]->SetLineColor(4); Edep_M[i]->SetLineWidth(2);}

  TH1F* Nhits_M_L[51]; TH1F* Nhits_M_R[51];
for (int j =1; j < 51; j++){
  Nhits_M_L[j] = new TH1F(Form("Nhits_L_M%d",j),"; Reconstructed E (MeV); Evt",0.3*inE,0,0.3*inE);
  Nhits_M_L[j]->Sumw2(); Nhits_M_L[j]->SetLineColor(3); Nhits_M_L[j]->SetLineWidth(2);
  Nhits_M_R[j] = new TH1F(Form("Nhits_R_M%d",j),"; Reconstructed E (MeV); Evt",0.3*inE,0,0.3*inE);
  Nhits_M_R[j]->Sumw2(); Nhits_M_R[j]->SetLineColor(3); Nhits_M_R[j]->SetLineWidth(2);
}

  TH1F* THits = new TH1F("THits",";ModuleN; # of p.e (%)",50,0,50);
  THits->Sumw2(); THits->SetLineColor(kRed); THits->SetLineWidth(2);
  TH1F* THits2 = new TH1F("THits2",";ModuleN; # of p.e.",50,0,50);
  THits2->Sumw2(); THits2->SetLineColor(kBlue); THits2->SetLineWidth(2);
 
  TH1F* THitss = new TH1F("THits",";ModuleN; Edep(MeV)",50,0,50);
  THitss->Sumw2(); THitss->SetLineColor(kRed); THitss->SetLineWidth(2);
  TH1F* THitss2 = new TH1F("THits2",";ModuleN; # of p.e.",50,0,50);
  THitss2->Sumw2(); THitss2->SetLineColor(kBlue); THitss2->SetLineWidth(2);

  TH1F* tP_leak = new TH1F("Pleak",";MeV;Evt",80,0.,800.*high);
  tP_leak->Sumw2(); tP_leak->SetLineWidth(2); tP_leak->SetLineColor(kRed);
  TH1F* tP_leak2 = new TH1F("Pleak2",";MeV;Evt",80,0.,800.*high);
  tP_leak2->Sumw2(); tP_leak2->SetLineWidth(2); tP_leak2->SetLineColor(kBlue);

  TH1F* Edep_ratio = new TH1F("Edep_ratio","; Beam E (MeV); Edep (%)",60,low*1000.,high*600.);
  Edep_ratio->Sumw2(); Edep_ratio->SetLineWidth(1); Edep_ratio->SetLineColor(2);
  TH1F* Edep_ratio2 = new TH1F("Edep_ratio2","; Beam E (MeV); Edep (%)",60,low*1000.,high*600.);
  Edep_ratio2->Sumw2(); Edep_ratio2->SetLineWidth(1); Edep_ratio2->SetLineColor(3);

  float rat_E = 0;  float rat_E2 = 0;
  double xBins = nFibers * nColumns;  double xLower = 0;  double xUpper = nFibers * nColumns;  double yBins = nLayers * nRows;
  double yLower = 0;  double yUpper = nLayers * nRows;

 TH2D* t2DhitS = new TH2D("2D Hit S1", "", xBins, xLower, xUpper, yBins, yLower, yUpper); t2DhitS->Sumw2(); t2DhitS->SetStats(0);
 TH2D* t2DhitS2 = new TH2D("2D Hit S2", "", xBins, xLower, xUpper, yBins, yLower, yUpper); t2DhitS2->Sumw2(); t2DhitS2->SetStats(0);

float Edep_Num[50] = {0};float Edep_Num2[50] = {0}; float Edep_Numm[50] = {0};

if (file1==1) {

 TFile *file = new TFile("/u/user/changhui/koBIC2026/BICprototype/build/analysis/100MeV_tree.root", "RECREATE");
 TTree *tree = new TTree("Ttree", "Edep and Nhits ");
// tree->SetBranchAddress("edep.Edep", &edep);
// tree->SetBranchAddress("edep.Module", &moduleN);

    Float_t edep1[50]={0};  // 24개의모듈 데이터를 저장
    Int_t moduleN_E[50]={0} ; Int_t moduleN_H[50]={0} ; Float_t Nhits_L[50]={0}; Float_t Nhits_R[50]={0}; Float_t Core_E[50]={0};
    tree->Branch("edep1", edep1, "edep1[24]/F");
    tree->Branch("moduleN_E",moduleN_E, "moduleN_E[24]/I");
    tree->Branch("moduleN_H",moduleN_H, "moduleN_H[24]/I");
    tree->Branch("Core_E", Core_E, "Core_E[24]/F");
    tree->Branch("Nhits_L", Nhits_L, "Nhits_L[24]/F");
    tree->Branch("Nhits_R", Nhits_R, "Nhits_R[24]/F");

    float totHit = 0; float nHits[50]= {0};
  unsigned int entries = drInterface->entries();
  while (drInterface->numEvt() < entries) {
    if (drInterface->numEvt() % 10 == 0) printf("1st analyzing %dth event ...\n", drInterface->numEvt());

    koBICInterface::koBICEventData drEvt;
    drInterface->read(drEvt);
//    float edep1[25]={0};  // 24개의모듈 데이터를 저장
//    float Nhits_L[25]={0}; float Nhits_R[25]={0};  int moduleN_E[25]={0} ; int moduleN_H[25]={0};

    float EdepCore = 0.;float Edep = 0.; float totE = 0.; float ratE = 0;  float Edep_Num[50] = {0}; float EdepC_Num[50] = {0};

 float Nhits_Numm_L[50]={0}; float Nhits_Numm_R[50]={0};
    for (auto edepItr = drEvt.Edeps.begin();  edepItr != drEvt.Edeps.end(); ++edepItr) {
      auto edep = *edepItr;

    for (int i=0; i<50; i++){
      if (edep.ModuleNum == i){
        Edep_Num[i] += edep.Edep;
        Edep_Numm[i] += edep.Edep;
        
	} 
 //     else {Edep_Num[i] = 0;}
    }
      Edep += edep.Edep;
//      EdepCore += edep.EdepCore;
    }

    tEdep->Fill(Edep);
float Eleak = 0;
    for (auto leak : drEvt.leaks) {


    for (int i=0; i<50; i++){
      if (leak.ModuleNum == i){
        EdepC_Num[i] += leak.EdepCore;
       }
       }
	    Eleak += leak.kE;
            EdepCore += leak.EdepCore;
}
    tEdepC->Fill(EdepCore/Edep);
tP_leak->Fill(Eleak);
for (int j=1; j<51; j++){
    Edep_M[j]->Fill(EdepC_Num[j-1]/Edep_Num[j-1]);
    edep1[j-1] = Edep_Num[j-1]; 
//    edepC[j-1] = EdepC_Num[j-1];    
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

  for (int i=0;i<50;i++) {
   if (moduleNum == i) {//moduleN_H=i;
      if (isLeft==0) {Nhits_Numm_L[i] += sipm->count; // Nhits_L[i]+=sipm->count; 
}
      else {Nhits_Numm_R[i] += sipm->count;// Nhits_R[i] += sipm->count;
}
    }
    }
        t2DhitS->Fill(nFibers*(moduleNum%nColumns)+fiberNum, nLayers*(moduleNum/nColumns)+plateNum, sipm->count); 
     }//SiPM loop
    }//tower loop

tHit_S->Fill(nHitS);  totHit += nHitS;

  for (int j=1; j<51; j++){
//    nHitS += Nhits_Numm_L[j-1]  ;
//    nHitS += Nhits_Numm_R[j-1] ;
    Nhits_M_L[j]->Fill(Nhits_Numm_L[j-1]/36.27);
    Nhits_M_R[j]->Fill(Nhits_Numm_R[j-1]/36.27);
    Core_E[j-1]=EdepC_Num[j-1];
    Nhits_L[j-1]=Nhits_Numm_L[j-1];
    Nhits_R[j-1]=Nhits_Numm_R[j-1];
    moduleN_H[j-1]=j;
//std::cout<<"Edep_Numm : "<<Edep_Num[j-1]<< "edep : " << edep1[j] << "moduleN : "<< moduleN_E[j] <<  j <<std::endl;
  } tree->Fill();
  } // event loop
    for (int z=0;z<50;z++) {  THits -> Fill(z,nHits[z]*100/totHit); }

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
tree->Write();
file->Close();
 }

if (file2==1) {
// second file
  RootInterface<koBICInterface::koBICEventData>* drInterface2 = new RootInterface<koBICInterface::koBICEventData>(std::string(filename2 + ".root"), true);
  drInterface2->set("koBIC","koBICEventData"); 

  int totHit2=0; float nHits2[24]={0};
  unsigned int entries2 = 10000;// drInterface2->entries();
  while (drInterface2->numEvt() < entries2) {
    if (drInterface2->numEvt() % 1000 == 0) printf("2nd analyzing %dth event ...\n", drInterface2->numEvt());

    koBICInterface::koBICEventData drEvt;
    drInterface2->read(drEvt);

    float Edep2 = 0.; float totE2 = 0.; float ratE2 = 0; float rat_E2 = 0;
    for (auto edepItr2 = drEvt.Edeps.begin(); edepItr2 != drEvt.Edeps.end(); ++edepItr2) {
      auto edep2 = *edepItr2;

for (int i2=0; i2<24; i2++){
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
    for (int zz=0;zz<24;zz++) {  THits2 -> Fill(zz,nHits2[zz]*100/totHit2); }
  
    rat_E2 = tEdep2->GetMean() / e2 ;
    Edep_ratio2-> SetBinContent(30,rat_E2*100);
    Edep_ratio2 -> SetBinError(30, tEdep2->GetStdDev()*100/e2);
//    std::cout << "4GeV BinCon = " << rat_E2*100 << std::endl;
//    std::cout << "4GeV BinErr = " << tEdep2->GetStdDev()*100/e2  << std::endl;

}
 for (int ii=1; ii<51; ii++){
     THitss->SetBinContent(ii,Edep_Numm[ii-1]/10000);
    }

  std::cout << "Module total1 = " << THitss->Integral() << std::endl;

    TF1 *gaussFit_E = new TF1("gaussFit_E", "gaus", 0, 5000);

  TCanvas* c = new TCanvas("c","");
  tEdep->Fit(gaussFit_E);  tEdep->SetStats(1);
  tEdep->Draw("Hist"); 
//  gaussFit_E->Draw("same");
     gStyle->SetOptFit(1);

 tEdep2->Draw("p same");  tEdep2->SetMarkerStyle(24); tEdep2->SetMarkerColor(kBlue); 
 tEdep->SetStats(1);
c->SaveAs(filename1+"compare_EdepG_e.png");   

    TF1 *gaussFit_S = new TF1("gaussFit_S", "gaus", 0, 200000);
/*
for (int j=1; j<25; j++){ c->cd(j);
//    Edep_M[j]->Fit(gaussFit_E);
    Edep_M[j]->Draw("Hist");   
//    gaussFit_E->Draw("same");  
//    gStyle->SetOptFit(1);
 c->SaveAs(Form(filename1+"EdepG_M%d.png",j));  }
for (int jj=1; jj<25; jj++){  c->cd(jj);
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
  
  THits->SetStats(0);
 c->SaveAs(filename1+"compare_towerE_e.png");
  
  THitss->Draw("Hist");
   

/* THitss2->Draw("Hist same"); */
//  THitss->SetStats(0);
//

  float Tot = THitss->Integral();
 for (int o=1; o<51; o++){ float val = THitss->GetBinContent(o); 
 std::cout << Form("Mod %d :  ",o) << val *100/Tot << "%,   " << val  << "MeV" <<  std::endl;   }
// c->SaveAs(filename1+"compare_tower_Edep.png");
 std::cout << "Total Edep, (%) : " << Tot <<  " MeV, "  << Tot*100/inE << "(%)" << std::endl;
  tHit_S->Fit(gaussFit_S);/* tHit_S2->Fit(gaussFit_S); */
   tHit_S->SetStats(1);
 //  tHit_S->GetYaxis()->SetRangeUser(0,100);
 
 tHit_S->SetMarkerColor(kBlue);  tHit_S->SetMarkerStyle(7);  tHit_S->Draw("Hist"); tHit_S2->SetMarkerColor(kBlue); tHit_S2->SetMarkerStyle(7);
    gaussFit_S->Draw("same");
    gStyle->SetOptFit(1);
//  tHit_S2->Draw("Hist p same" ); tHit_S2->SetMarkerStyle(24); 
c->SaveAs(filename1 +"compare_nHitpEventSG_cal.png");

    t2DhitS->Draw("COLZ"); c->SaveAs(filename1+"_n2DHitS.png");
//  t2DhitS2->Draw("COLZ"); c->SaveAs(filename2+"_n2DHitS.png");

Edep_ratio->Draw("p E1"); Edep_ratio2->Draw("p E1 same ");
gStyle->SetErrorX(0); Edep_ratio -> SetStats(0);  
Edep_ratio->SetMarkerSize(1.3); Edep_ratio->SetMarkerStyle(20); Edep_ratio->GetYaxis()->SetRangeUser(70,100); Edep_ratio->SetMarkerColor(2);
Edep_ratio2->SetMarkerSize(1.3); Edep_ratio2->SetMarkerStyle(20); Edep_ratio2->SetMarkerColor(3); 
/* 
   TFile *outputFile1 = new TFile("/u/user/changhui/work/BICprototype2/build/analysis/resol.root", "RECREATE");    
tHit_S->Write();tHit_S2->Write();tHit_S3->Write();tHit_S4->Write();tHit_S5->Write();tHit_S6->Write();
outputFile1->Close();   
*/

   TFile *outputFile2 = new TFile("/u/user/changhui/koBIC2026/BICprototype/build/analysis/100MeV_hist.root", "RECREATE");
tEdep->Write(); tEdepC->Write(); tHit_S->Write();
for (int jjj=1; jjj<51; jjj++){ 
    Edep_M[jjj]->Write();
}
for (int jjjj=1; jjjj<51; jjjj++){
    Nhits_M_L[jjjj]->Write(); Nhits_M_R[jjjj]->Write();
}
outputFile2->Close();

c->SaveAs(filename1+"compare_Edep_ratio_e.png");

  TCanvas* c1 = new TCanvas("c1","");
c1->SetCanvasSize(2400,1200);
c1->Divide(8,4);

for (int j=1; j<51; j++){   //  if (j<1 || j>50) continue;
	c1->cd(j); gPad->Clear();// c1->SetMargin(0.01, 0.01, 0.01, 0.01);
        Edep_M[j]->Draw("Hist");
        Edep_M[j]->GetYaxis()->SetMaxDigits(3);
	gStyle -> SetOptFit(0); 

	TPad*pad=(TPad*)gPad;
        pad->SetTopMargin(0.03);
        pad->SetBottomMargin(0.1);
        pad->SetLeftMargin(0.1);
        pad->SetRightMargin(0.01);
       
    TPaveStats *stats = (TPaveStats*)Edep_M[j]->GetListOfFunctions()->FindObject("stats");
    if (stats) { 
        stats->SetTextSize(0.06);
stats->SetX1NDC(0.55);stats->SetX2NDC(0.99);stats->SetY1NDC(0.55);stats->SetY2NDC(0.99);
 }}
c1->SaveAs(filename1+"Edep_Allmod.png");
/*
TCanvas* c2 = new TCanvas("c1","");
c2->SetCanvasSize(2400,1200);
c2->Divide(8,4);

for (int jj=1; jj<51; jj++){  
    c1->cd(jj); gPad->Clear();
    Nhits_M_R[jj]->Draw("Hist");
    Nhits_M_R[jj]->GetYaxis()->SetMaxDigits(3);
    gStyle->SetOptFit(0);

	TPad*pad1=(TPad*)gPad;
        pad1->SetTopMargin(0.03);
        pad1->SetBottomMargin(0.1);
        pad1->SetLeftMargin(0.1);
        pad1->SetRightMargin(0.01);

    TPaveStats *stats1 = (TPaveStats*)Nhits_M_R[jj]->GetListOfFunctions()->FindObject("stats");
    if (stats1) {
        stats1->SetTextSize(0.06);
stats1->SetX1NDC(0.55);stats1->SetX2NDC(0.99);stats1->SetY1NDC(0.55);stats1->SetY2NDC(0.99);

} }c1->SaveAs(filename1+"Nhits_Allmod_R.png");

for (int jj=1; jj<51; jj++){
    c1->cd(jj); gPad->Clear();
    Nhits_M_L[jj]->Draw("Hist");
    Nhits_M_L[jj]->GetYaxis()->SetMaxDigits(3);
    gStyle->SetOptFit(0);

        TPad*pad1=(TPad*)gPad;
        pad1->SetTopMargin(0.03);
        pad1->SetBottomMargin(0.1);
        pad1->SetLeftMargin(0.1);
        pad1->SetRightMargin(0.01);

    TPaveStats *stats1 = (TPaveStats*)Nhits_M_L[jj]->GetListOfFunctions()->FindObject("stats");
    if (stats1) {
        stats1->SetTextSize(0.06);
stats1->SetX1NDC(0.55);stats1->SetX2NDC(0.99);stats1->SetY1NDC(0.55);stats1->SetY2NDC(0.99);

 }c1->SaveAs(filename1+"Nhits_Allmod_L.png");


*/
/*
  for ( int Bin = 0 ; Bin < 16 ; ++Bin) {

}*/// Get bin counts
}

