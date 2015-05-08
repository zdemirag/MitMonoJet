//--------------------------------------------------------------------------------------------------
// Make reweighting histos given an observable, a target and a source
//
// Authors: L. Di Matteo                                                                  (Feb 2015)
//----------------

void getReweightingHistos() {

  TFile *OutFile = new TFile("weights_new.root","RECREATE");  
  TH1F* histoSource = new TH1F("histoSource","histoSource",20,0,1);
  histoSource->Sumw2();
  TH1F* histoTarget = (TH1F*) histoSource->Clone("histoTarget");
  
  std::vector<TString> inSourceFiles;
  inSourceFiles.push_back("/scratch5/dimatteo/cms/hist/boostedv-v12/merged-pro/boostedv-v12_GJets_HT-40To100_8TeV-madgraph+Summer12_DR53X-PU_S10_START53_V19-v1+AODSIM_noskim_flatntuple.root");
  inSourceFiles.push_back("/scratch5/dimatteo/cms/hist/boostedv-v12/merged-pro/boostedv-v12_GJets_HT-100To200_8TeV-madgraph+Summer12_DR53X-PU_S10_START53_V19-v1+AODSIM_noskim_flatntuple.root");
  inSourceFiles.push_back("/scratch5/dimatteo/cms/hist/boostedv-v12/merged-pro/boostedv-v12_GJets_HT-200To400_8TeV-madgraph+Summer12_DR53X-PU_S10_START53_V7A-v1+AODSIM_noskim_flatntuple.root");
  inSourceFiles.push_back("/scratch5/dimatteo/cms/hist/boostedv-v12/merged-pro/boostedv-v12_GJets_HT-400ToInf_8TeV-madgraph+Summer12_DR53X-PU_S10_START53_V7A-v1+AODSIM_noskim_flatntuple.root");
  std::vector<float> inSourceXsec;
  inSourceXsec.push_back(20930);
  inSourceXsec.push_back(5212);
  inSourceXsec.push_back(960);
  inSourceXsec.push_back(107);

  TString inTargetFile = "/scratch5/dimatteo/cms/hist/boostedv-v12/merged-pro/boostedv-v12_s12-zjets-ptz100-v7a_noskim_flatntuple.root";
  //TString inTargetFile = "/scratch5/dimatteo/cms/hist/boostedv-v12/merged-pro/boostedv-v12_s12-ww-v7a_noskim_flatntuple.root";
  
  for (unsigned int isample = 0; isample < inSourceFiles.size(); isample++) {
    TFile *MyFile = new TFile(inSourceFiles[isample],"READ");
    float nGenEvts = ((TH1D*)MyFile->FindObjectAny("hDAllEvents"))->GetEntries();
    TTree *MyTree = (TTree*)MyFile->FindObjectAny("DMSTree");
    //TString cut = Form("%f*(fjet1QGtag>-0.5)",inSourceXsec[isample]/nGenEvts);
    TString cut = Form("%f*(fjet1QGtag>-0.5 && fjet1Tau2/fjet1Tau1 < 0.5 && 60 < fjet1MassPruned && fjet1MassPruned < 110)",inSourceXsec[isample]/nGenEvts);
    OutFile->cd();
    MyTree->Draw("fjet1QGtag >>+ histoSource",cut,"goff");
  }//end loop on Gjets samples
  histoSource->Scale(1./histoSource->GetSumOfWeights());

  TFile *MyFile = new TFile(inTargetFile,"READ");
  float nGenEvts = ((TH1D*)MyFile->FindObjectAny("hDAllEvents"))->GetEntries();
  TTree *MyTree = (TTree*)MyFile->FindObjectAny("DMSTree");
  OutFile->cd();
  MyTree->Draw("fjet1QGtag >>+ histoTarget","fjet1QGtag>-0.5 && fjet1Tau2/fjet1Tau1 < 0.5 && 60 < fjet1MassPruned && fjet1MassPruned < 110","goff");
  histoTarget->Scale(1./histoTarget->GetSumOfWeights());
  
  TH1F* histoWeight = (TH1F*) histoTarget->Clone("histoWeight");
  histoWeight->Divide(histoSource);
  
  TF1 *fitfunc = new TF1("fitfunc","pol3",0,1); 
  histoWeight->Fit("fitfunc");
  
  TCanvas* can = new TCanvas();
  histoWeight->Draw();  

  TCanvas* can1 = new TCanvas();
  histoTarget->Draw();  

  TCanvas* can2 = new TCanvas();
  histoSource->Draw();  
  
  histoWeight->Write();
  
  
  return;
}
