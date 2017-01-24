void ComparePlot(){
  const float pi = 3.141592;
  const Int_t nfiles = 4;
  const Int_t nPtBins = 5;
  Int_t g,h,i,j;

  // create a new Root file
  TFile *top = new TFile("GeantPythiaStar.root","recreate");

  double ptBinLo[nPtBins] = { 3, 5, 9, 13, 17 };
  double ptBinHi[nPtBins] = { 4, 8, 12, 16, 24 };
  TString ptBinString[nPtBins] = { "0.5-1.0", "1.0-2.0", "2.0-3.0", "3.0-4.0", "4.0-6.0" };
  TString PtBinName[nPtBins] = {"assoc05to10", "assoc10to20", "assoc20to30", "assoc30to40", "assoc40to60"};
  TString minSubPt[nfiles] = {"7",  "8",  "9", "10"};
  TString minLeadPt[nfiles] = {"14", "16", "18", "20"};
  TString LeadSub[nfiles] = {"lead14_sub7", "lead16_sub8", "lead18_sub9", "lead20_sub10", };
  TString CTitle[nPtBins] = {"0.5 < Pt_{assoc} < 1.0","1.0 < Pt_{assoc} < 2.0","2.0 < Pt_{assoc} < 3.0","3.0 < Pt_{assoc} < 4.0","4.0 < Pt_{assoc} < 6.0"};
  int color[nfiles] = { 2, 51, 4, 8 };
  TDirectory *ptbin[nPtBins];

  //  HISTOGRAMS
  TH1D *GeantCorr[nPtBins][nfiles];
  TH1D *PythiaCorr[nPtBins][nfiles];
  TH1D *StarCorr[nPtBins][nfiles];
  TH1D* GeantPtPlot[nPtBins][nfiles];
  TH1D* PythiaPtPlot[nPtBins][nfiles];
  TH1D* StarPtPlot[nPtBins][nfiles];
  TF1 *GeantFit[nfiles];
  TF1 *PythiaFit[nfiles];
  TF1 *StarFit[nfiles];

  

  for (h=0; h<nPtBins; h++) {
    TString directory = PtBinName[h];
    ptbin[h] = top->mkdir(directory);   // create a new subdirectory for each Pt bin range
    ptbin[h]->cd();
    
    for (j=0; j<nfiles; j++){
      TString importName = "geant_";
      importName += LeadSub[j];
      importName += ".root";
      TFile* GeantDijetFILE = new TFile( importName, "READ");
      TH3D* GeantLeadCorr = (TH3D*) GeantDijetFILE->Get("ppleadjetcorr");
      TH3D* GeantSubCorr = (TH3D*) GeantDijetFILE->Get("ppsubjetcorr");
      TH3D* GeantLeadPt = (TH3D*) GeantDijetFILE->Get("ppleadjetpt");
      TH2D* GeantEvents = (TH2D*) GeantDijetFILE->Get("binvzdist");
      TH2D* GeantLead = (TH2D*) GeantLeadCorr->Project3D("ZY");
      TH2D* GeantSub = (TH2D*) GeantSubCorr->Project3D("ZY");
      
      TString nameSet = "GeantLeadCorr";
      GeantLead->SetName( nameSet );
      nameSet = "GeantSubCorr";
      GeantSub->SetName( nameSet );
      nameSet = "GeantLeadCorr_";
      nameSet += LeadSub[j];
      GeantLead->SetName( nameSet );

      ptbin[h]->cd();
      GeantCorr[h][j] = (TH1D*) GeantLead->ProjectionX( nameSet , ptBinLo[j], ptBinHi[j] );
      GeantCorr[h][j]->Scale(1/GeantCorr[i][j]->GetBinWidth(1));
      GeantCorr[h][j]->Scale( 1/double(GeantEvents->Integral()) );

      // FIT
      TString fiteq = "[0]+[1]*exp(-0.5*((x-[2])/[3])**2)+[4]*exp(-0.5*((x-[5])/[6])**2)";
      double phiMin = -pi + pi/2.0;
      double phiMax = pi + pi/2.0;
      nameSet = "GeantFit_";
      nameSet += LeadSub[j];
      
      GeantFit[j] = new TF1( nameSet, fiteq, phiMin, phiMax);
      GeantFit[j]->FixParameter(2, 0);
      GeantFit[j]->FixParameter(5, pi);
      GeantFit[j]->SetParameter(3, 0.2);
      GeantFit[j]->SetParameter(6, 0.2);
      GeantFit[j]->SetLineColor(color[i]);
      GeantFit[j]->SetLineWidth(2);
      GeantCorr[h][j]->Fit(GeantFit[j]);
      
      gROOT->SetEditHistograms();
      GeantCorr[h][j]->SetLineColor(color[i]);
      GeantCorr[h][j]->SetLineWidth(2);
      //GeantCorr[h][j]->SetMaximum(0.6);
      GeantCorr[h][j]->SetMinimum(0);
      gStyle->SetOptStat(0);

      // WRITE
      GeantCorr[h][j]->Write();
    }
    
    GeantLead->Write();
    GeantSub->Write();
    GeantEvents->Write();
    GeantLeadPt->Write();
    
    
  }

  delete top;
}

     

 
