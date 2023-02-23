#include "TriggerHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "TriggerCard.h"
#include "JDrawer.h"
#include "../src/TriggerHistograms.h"

/*
 * Macro for comparing weighted and un-weighted Monte Carlo simulations with data
 */
void triggerEfficiencyPlotter(){
  
  // File containing the base trigger jet spectrum and selected trigger jet spectra on top of that
  TString fileName = "veryCoolData_processed.root";
  
  // Open the file and check that it exists
  TFile* inputFile = TFile::Open(fileName);
  
  if(inputFile == NULL){
    cout << "Error! The file " << fileName.Data() << " does not exist!" << endl;
    cout << "Maybe you forgot the data/ folder path?" << endl;
    cout << "Will not execute the code" << endl;
    return;
  }
  
  // Check if we are using PbPb or pp data
  TriggerCard *systemCard = new TriggerCard(inputFile);
  TString collisionSystem = systemCard->GetDataType();
  bool isPbPbData = collisionSystem.Contains("PbPb");
  
  int nCentralityBins = systemCard->GetNCentralityBins();
  if(!isPbPbData) nCentralityBins = 1;
  
  // ====================================================
  //                    Configuration
  // ====================================================
  
  // Figure saving
  const bool saveFigures = false;  // Save figures
  const char* saveComment = "";   // Comment given for this specific file
  const char* figureFormat = "pdf"; // Format given for the figures
  
  // Create and setup a new histogram managers to project and handle the histograms
  TriggerHistogramManager *histograms;
  histograms = new TriggerHistogramManager(inputFile);
  histograms->SetLoadJetHistograms(true);
  histograms->LoadProcessedHistograms();
  
  // Jet pT histograms to calculate trigger turn on curve
  TH1D* hJetPt[TriggerHistograms::knTriggerTypes+1][nCentralityBins];   // Jet pT distributions
  TH1D* hJetTurnOn[TriggerHistograms::knTriggerTypes][nCentralityBins]; // Trigger efficiency curve calculated from jet pT distributions
  
  // Initialize all histograms to NULL
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iTrigger = 0; iTrigger < TriggerHistograms::knTriggerTypes; iTrigger++){
      hJetPt[iTrigger][iCentrality] = NULL;
      hJetTurnOn[iTrigger][iCentrality] = NULL;
    } // Data type loop
    hJetPt[TriggerHistograms::knTriggerTypes][iCentrality] = NULL;
  }
  
  // Read the histograms from the file
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iTrigger = 0; iTrigger <= TriggerHistograms::knTriggerTypes; iTrigger++){
      hJetPt[iTrigger][iCentrality] = histograms->GetHistogramJetPt(iCentrality, TriggerHistograms::kReconstructed, iTrigger);
    }
  }
  
  // Calculate the ratio to determine trigger efficiency
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iTrigger = 0; iTrigger < TriggerHistograms::knTriggerTypes; iTrigger++){
      hJetTurnOn[iTrigger][iCentrality] = (TH1D*) hJetPt[iTrigger][iCentrality]->Clone(Form("turnOn%d%d", iCentrality, iTrigger));
      hJetTurnOn[iTrigger][iCentrality]->Divide(hJetPt[TriggerHistograms::knTriggerTypes][iCentrality]);
    }
  }
  
  
  // ===============================================
  //        Draw the trigger turn-on curves
  // ===============================================
  
  JDrawer *drawer = new JDrawer();
  TLegend *legend;
  TLine *oneLine = new TLine(0,1,500,1);
  oneLine->SetLineStyle(2);
  oneLine->SetLineColor(kBlack);
  TLine *oneTwentyLine = new TLine(120,0,120,1.3);
  oneTwentyLine->SetLineStyle(2);
  oneTwentyLine->SetLineColor(kRed);
  TString centralityString;
  TString compactCentralityString;
  TString triggerString;
  
  TriggerHistograms *triggerProvider = new TriggerHistograms();
  
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    
    if(isPbPbData){
      centralityString = Form("C: %.0f-%.0f", systemCard->GetLowBinBorderCentrality(iCentrality), systemCard->GetHighBinBorderCentrality(iCentrality));
      compactCentralityString = Form("_C%.0f-%.0f", systemCard->GetLowBinBorderCentrality(iCentrality), systemCard->GetHighBinBorderCentrality(iCentrality));
    } else {
      centralityString = "pp";
      compactCentralityString = "_pp";
    }
    
    for(int iTrigger = 0; iTrigger < TriggerHistograms::knTriggerTypes; iTrigger++){
      
      // Create a legend for the figure
      legend = new TLegend(0.2,0.73,0.5,0.88);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      
      // Draw the turn-on curve
      hJetTurnOn[iTrigger][iCentrality]->GetYaxis()->SetRangeUser(0,2);
      drawer->DrawHistogram(hJetTurnOn[iTrigger][iCentrality], "Jet p_{T} (GeV)", "Trigger efficiency", " ");
      
      // Add information to legend
      legend->AddEntry((TObject*)0, centralityString.Data(), "");
      triggerString = Form("%s / %s", triggerProvider->GetTriggerName(iTrigger).Data(), triggerProvider->GetTriggerName(systemCard->GetBaseTrigger()).Data());
      legend->AddEntry((TObject*)0, triggerString.Data(), "");
      
      // Draw the legend
      legend->Draw();
      
      // Draw a line to one
      oneLine->Draw("same");
      oneTwentyLine->Draw("same");
      
      // Save the figures to a file
      if(saveFigures){
        gPad->GetCanvas()->SaveAs(Form("figures/jetTriggerTurnOn%s%s%s.%s", triggerProvider->GetTriggerName(iTrigger).Data(), saveComment, compactCentralityString.Data(), figureFormat));
      }
      
    } // Trigger selection loop
  } // Centrality loop
  
  delete triggerProvider;
}
