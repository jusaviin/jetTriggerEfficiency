#include "TriggerCard.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "TriggerHistogramManager.h"

#include <bitset>

/*
 * Macro for projecting the histograms needed in the energy-energy correlator analysis from the THnSparses
 *
 *  Arguments:
 *   TString inputFileName = File from which the histograms are read
 *   const char* outputFileName = If we are producing output file, name of the output file
 */
void projectTriggerHistograms(TString inputFileName = "veryCoolData.root", const char* outputFileName = "veryCoolData_processed.root"){

  // Print the file name to console
  cout << "Projecting histograms from " << inputFileName.Data() << endl;
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  // If we write a file, define the output name and write mode
  const char* fileWriteMode = "UPDATE";
  
  // ====================================================
  //  Binning configuration for the projected histograms
  // ====================================================
  
  // Option to read all the binning information from TriggerCard used to create the file
  const bool readCentralityBinsFromFile = false;
  
  // If not reading the bins from the file, manually define new bin borders
  const int nCentralityBins = 4;
  double centralityBinBorders[nCentralityBins+1] = {0,10,30,50,90};   // Bin borders for centrality
  
  // Projected bin range
  int firstDrawnCentralityBin = 0;
  int lastDrawnCentralityBin = nCentralityBins-1;
  
  // ==================================================================
  // ===================== Configuration ready ========================
  // ==================================================================
  
  // Open the input file
  TFile *inputFile = TFile::Open(inputFileName);
  
  if(inputFile == NULL){
    cout << "Error! The file " << inputFileName.Data() << " does not exist!" << endl;
    cout << "Maybe you forgot the data/ folder path?" << endl;
    cout << "Will not execute the code" << endl;
    return;
  }
  
  // Load the card from the file and read the collision system
  TriggerCard *card = new TriggerCard(inputFile);
  TString collisionSystem = card->GetDataType();
  
  // Remove centrality selection from pp data
  if(collisionSystem.Contains("pp")){
    lastDrawnCentralityBin = 0;
    centralityBinBorders[0] = -0.5;
  }
  
  // If we change the binning, save the new binning to the card
  if(!readCentralityBinsFromFile) card->AddVector(TriggerCard::kCentralityBinEdges,nCentralityBins+1,centralityBinBorders);
  
  // Add information about the used input files to the card
  card->AddFileName(TriggerCard::kInputFileName,inputFileName);
  
  // The git hash here will be replaced by the latest commit hash by projectHistogramsInSteps.sh script
  const char* gitHash = "GITHASHHERE";
  card->AddProjectionGitHash(gitHash);
  
  // ============================== //
  //     TriggerHistogramManager    //
  // ============================== //
    
  // Create and setup a new histogram manager to project and handle the histograms
  TriggerHistogramManager *histograms = new TriggerHistogramManager(inputFile,card);
  
  // Set which histograms to project from the input file
  histograms->SetLoadEventInformation(true);
  histograms->SetLoadJetHistograms(true);
  histograms->SetLoad2DHistograms(true);

  // Set the binning information
  histograms->SetCentralityBins(readCentralityBinsFromFile,nCentralityBins,centralityBinBorders,true);
  if(!readCentralityBinsFromFile) histograms->SetCentralityBinRange(firstDrawnCentralityBin,lastDrawnCentralityBin);
  
  // Project the one dimensional histograms from the THnSparses
  histograms->LoadHistograms();
  
  // Save the histograms to an output file
  histograms->Write(outputFileName,fileWriteMode);
  
}
