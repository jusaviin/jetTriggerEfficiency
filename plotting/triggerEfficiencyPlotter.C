#include "TriggerHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "TriggerCard.h"
#include "JDrawer.h"
#include "../src/TriggerHistograms.h"

/*
 * Macro for comparing weighted and un-weighted Monte Carlo simulations with data
 */
void triggerEfficiencyPlotter(){
  
  // File containing the base trigger jet spectrum and selected trigger jet spectra on top of that
  TString fileName = "data/triggerAnalysis_akFlowJets_includeLeading_eta1v6_baseCalo60_processed_2023-03-01.root";
  
  // triggerAnalysis_akFlowJets_cutBadPhi_eta1v6_baseCalo60_processed_2023-03-01.root
  // triggerAnalysis_akFlowJets_includeLeading_eta1v6_baseCalo60_processed_2023-03-01.root
  // triggerAnalysis_akFlowJets_eta1v6_baseCalo40_processed_2023-02-22.root
  
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
  
  // Select which triggers to compare
  vector<int> triggerComparisonIndex{TriggerHistograms::kCalo100, TriggerHistograms::kCalo80};
  
  // Select which histograms to draw
  const bool drawSingleTurnOns = false;
  const bool drawComparisons = true;
  
  // Select which jet distribution to use for drawing
  bool drawJets[TriggerHistogramManager::knJetTypes];
  drawJets[TriggerHistogramManager::kInclusiveJet] = true;
  drawJets[TriggerHistogramManager::kLeadingJet] = false;
  
  // Figure saving
  const bool saveFigures = true;  // Save figures
  const char* saveComment = "_base60Comparison";   // Comment given for this specific file
  const char* figureFormat = "pdf"; // Format given for the figures
  
  // Create and setup a new histogram managers to project and handle the histograms
  TriggerHistogramManager *histograms;
  histograms = new TriggerHistogramManager(inputFile);
  histograms->SetLoadJetHistograms(true);
  histograms->LoadProcessedHistograms();
  
  // Jet pT histograms to calculate trigger turn on curve
  TH1D* hJetPt[TriggerHistogramManager::knJetTypes][TriggerHistograms::knTriggerTypes+1][nCentralityBins];   // Jet pT distributions
  TH1D* hJetTurnOn[TriggerHistogramManager::knJetTypes][TriggerHistograms::knTriggerTypes][nCentralityBins]; // Trigger efficiency curve calculated from jet pT distributions
  
  // Initialize all histograms to NULL
  for(int iJetType = 0; iJetType < TriggerHistogramManager::knJetTypes; iJetType++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iTrigger = 0; iTrigger < TriggerHistograms::knTriggerTypes; iTrigger++){
        hJetPt[iJetType][iTrigger][iCentrality] = NULL;
        hJetTurnOn[iJetType][iTrigger][iCentrality] = NULL;
      } // Data type loop
      hJetPt[iJetType][TriggerHistograms::knTriggerTypes][iCentrality] = NULL;
    } // Centrality loop
  } // Jet type loop
  
  // Read the histograms from the file
  for(int iJetType = 0; iJetType < TriggerHistogramManager::knJetTypes; iJetType++){
    if(!drawJets[iJetType]) continue;
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iTrigger = 0; iTrigger <= TriggerHistograms::knTriggerTypes; iTrigger++){
        hJetPt[iJetType][iTrigger][iCentrality] = histograms->GetHistogramJetPt(iJetType, iCentrality, TriggerHistograms::kReconstructed, iTrigger);
      } // Trigger selection loop
    } // Centrality loop
  } // Jet type loop
  
  // Calculate the ratio to determine trigger efficiency
  for(int iJetType = 0; iJetType < TriggerHistogramManager::knJetTypes; iJetType++){
    if(!drawJets[iJetType]) continue;
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iTrigger = 0; iTrigger < TriggerHistograms::knTriggerTypes; iTrigger++){
        hJetTurnOn[iJetType][iTrigger][iCentrality] = (TH1D*) hJetPt[iJetType][iTrigger][iCentrality]->Clone(Form("turnOn%d%d%d", iJetType, iCentrality, iTrigger));
        hJetTurnOn[iJetType][iTrigger][iCentrality]->Divide(hJetPt[iJetType][TriggerHistograms::knTriggerTypes][iCentrality]);
      } // Trigger selection loop
    } // Centrality loop
  } // Jet type loop
  
  
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
  
  // Selection of colors and styles
  int color[] = {kBlack,kRed,kBlue,kGreen+4};
  int markerStyle[] = {kOpenSquare, kOpenCircle, kOpenDiamond, kOpenCross};
  
  TriggerHistograms *triggerProvider = new TriggerHistograms();
  
  // Draw turnons for all included triggers
  if(drawSingleTurnOns){
    for(int iJetType = 0; iJetType < TriggerHistogramManager::knJetTypes; iJetType++){
      if(!drawJets[iJetType]) continue;
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        
        if(isPbPbData){
          centralityString = Form("Centrality: %.0f-%.0f", systemCard->GetLowBinBorderCentrality(iCentrality), systemCard->GetHighBinBorderCentrality(iCentrality));
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
          hJetTurnOn[iJetType][iTrigger][iCentrality]->GetYaxis()->SetRangeUser(0,2);
          drawer->DrawHistogram(hJetTurnOn[iJetType][iTrigger][iCentrality], "Jet p_{T} (GeV)", "Trigger efficiency", " ");
          
          // Add information to legend
          legend->AddEntry((TObject*)0, centralityString.Data(), "");
          triggerString = Form("%s / %s", triggerProvider->GetTriggerName(iTrigger).Data(), triggerProvider->GetTriggerName(systemCard->GetBaseTrigger()).Data());
          legend->AddEntry((TObject*)0, triggerString.Data(), "");
          
          // Draw the legend
          legend->Draw();
          
          // Draw a lines to one and 120 GeV
          oneLine->Draw("same");
          oneTwentyLine->Draw("same");
          
          // Save the figures to a file
          if(saveFigures){
            gPad->GetCanvas()->SaveAs(Form("figures/%sTriggerTurnOn%s%s%s.%s", histograms->GetJetHistogramName(iJetType), triggerProvider->GetTriggerName(iTrigger).Data(), saveComment, compactCentralityString.Data(), figureFormat));
          }
          
        } // Trigger selection loop
      } // Centrality loop
    } // Jet type loop
  } // Drawing single distributions
  
  // Draw turn-ons for selected triggers to the same plot
  if(drawComparisons){
    for(int iJetType = 0; iJetType < TriggerHistogramManager::knJetTypes; iJetType++){
      if(!drawJets[iJetType]) continue;
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        
        if(isPbPbData){
          centralityString = Form("Centrality: %.0f-%.0f", systemCard->GetLowBinBorderCentrality(iCentrality), systemCard->GetHighBinBorderCentrality(iCentrality));
          compactCentralityString = Form("_C%.0f-%.0f", systemCard->GetLowBinBorderCentrality(iCentrality), systemCard->GetHighBinBorderCentrality(iCentrality));
        } else {
          centralityString = "pp";
          compactCentralityString = "_pp";
        }
        
        // Draw the turn-on curve
        hJetTurnOn[iJetType][triggerComparisonIndex.at(0)][iCentrality]->GetYaxis()->SetRangeUser(0,2);
        hJetTurnOn[iJetType][triggerComparisonIndex.at(0)][iCentrality]->SetMarkerStyle(markerStyle[0]);
        hJetTurnOn[iJetType][triggerComparisonIndex.at(0)][iCentrality]->SetMarkerColor(color[0]);
        hJetTurnOn[iJetType][triggerComparisonIndex.at(0)][iCentrality]->SetLineColor(color[0]);
        drawer->DrawHistogram(hJetTurnOn[iJetType][triggerComparisonIndex.at(0)][iCentrality], "Jet p_{T} (GeV)", "Trigger efficiency", " ");
        
        // Add other histograms to same figure
        for(int iComparison = 1; iComparison < triggerComparisonIndex.size(); iComparison++){
          hJetTurnOn[iJetType][triggerComparisonIndex.at(iComparison)][iCentrality]->SetMarkerStyle(markerStyle[iComparison]);
          hJetTurnOn[iJetType][triggerComparisonIndex.at(iComparison)][iCentrality]->SetMarkerColor(color[iComparison]);
          hJetTurnOn[iJetType][triggerComparisonIndex.at(iComparison)][iCentrality]->SetLineColor(color[iComparison]);
          hJetTurnOn[iJetType][triggerComparisonIndex.at(iComparison)][iCentrality]->Draw("same");
        }
        
        // Create a legend for the figure
        legend = new TLegend(0.16,0.84-0.06*triggerComparisonIndex.size(),0.43,0.9);
        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
        
        // Add information to legend
        legend->AddEntry((TObject*)0, centralityString.Data(), "");
        for(int iComparison = 0; iComparison < triggerComparisonIndex.size(); iComparison++){
          triggerString = Form("%s / %s", triggerProvider->GetTriggerName(triggerComparisonIndex.at(iComparison)).Data(), triggerProvider->GetTriggerName(systemCard->GetBaseTrigger()).Data());
          legend->AddEntry(hJetTurnOn[iJetType][triggerComparisonIndex.at(iComparison)][iCentrality], triggerString.Data(), "p");
        }
        
        // Draw the legend
        legend->Draw();
        
        // Draw a lines to one and 120 GeV
        oneLine->Draw("same");
        oneTwentyLine->Draw("same");
        
        // Save the figures to a file
        if(saveFigures){
          gPad->GetCanvas()->SaveAs(Form("figures/%sTriggerTurnOnComparison%s%s.%s", histograms->GetJetHistogramName(iJetType), saveComment, compactCentralityString.Data(), figureFormat));
        }
        
      } // Centrality loop
    } // Jet type loop
  } // Drawing comparison plots
  
  delete triggerProvider;
}
