/*
 * Implementation of TriggerHistogramManager
 */

// Root includes
#include <THnSparse.h>
#include <TPad.h>

// Own includes
#include "TriggerHistogramManager.h"

/*
 * Default constructor
 */
TriggerHistogramManager::TriggerHistogramManager() :
  fInputFile(NULL),
  fCard(NULL),
  fSystemAndEnergy(""),
  fCompactSystemAndEnergy(""),
  fLoadEventInformation(false),
  fLoadJets(false),
  fLoad2DHistograms(false),
  fFirstLoadedCentralityBin(0),
  fLastLoadedCentralityBin(1),
  fnCentralityBins(kMaxCentralityBins)
{
  
  // Default binning for centrality
  for(int iCentrality = 0; iCentrality < kMaxCentralityBins + 1; iCentrality++){
    fCentralityBinIndices[iCentrality] = iCentrality+1;
    fCentralityBinBorders[iCentrality] = 0;
  }
  
  // Initialize all the other histograms to null
  fhVertexZ = NULL;            // Vertex z position
  fhVertexZWeighted = NULL;    // Weighted vertex z-position (only meaningfull for MC)
  fhEvents = NULL;             // Number of events surviving different event cuts
  fhCentrality = NULL;         // Centrality of all events
  fhCentralityWeighted = NULL; // Weighted centrality distribution in all events (only meaningful for MC)
  fhPtHat = NULL;              // pT hat for MC events (only meaningful for MC)
  fhPtHatWeighted = NULL;      // Weighted pT hat distribution (only meaningful for MC)
  
  // Centrality loop
  for(int iCentrality = 0; iCentrality < kMaxCentralityBins; iCentrality++){
    for(int iDataLevel = 0; iDataLevel < TriggerHistograms::knDataLevels; iDataLevel++){
      for(int iTrigger = 0; iTrigger <= TriggerHistograms::knTriggerTypes; iTrigger++){
        for(int iJetType = 0; iJetType < knJetTypes; iJetType++){
          
          // Jet histograms
          fhJetPt[iJetType][iCentrality][iDataLevel][iTrigger] = NULL;      // Jet pT histograms
          fhJetPhi[iJetType][iCentrality][iDataLevel][iTrigger] = NULL;     // Jet phi histograms
          fhJetEta[iJetType][iCentrality][iDataLevel][iTrigger] = NULL;     // Jet eta histograms
          fhJetEtaPhi[iJetType][iCentrality][iDataLevel][iTrigger] = NULL;  // 2D eta-phi histogram for jets
          
        } // Jet type loop
      } // Trigger selection loop
    } // Data level loop
  } // Centrality loop
}

/*
 * Constructor
 */
TriggerHistogramManager::TriggerHistogramManager(TFile *inputFile) :
  TriggerHistogramManager()
{
  fInputFile = inputFile;
  
  // Read card from inputfile
  fCard = new TriggerCard(inputFile);
  
  // Initialize values using the information in card
  InitializeFromCard();
  
}

/*
 * Constructor
 */
TriggerHistogramManager::TriggerHistogramManager(TFile *inputFile, TriggerCard *card) :
  TriggerHistogramManager()
{
  fInputFile = inputFile;
  
  // Initialize values using the information in card
  fCard = card;
  InitializeFromCard();
  
}

/*
 * Initialize several member variables from TriggerCard
 */
void TriggerHistogramManager::InitializeFromCard(){
  
  // Read the collision system from the card
  TString collisionSystem = fCard->GetDataType();
  
  // Make a string for collision system based on information on the card
  fSystemAndEnergy = Form("%s 5.02 TeV",collisionSystem.Data());
  fCompactSystemAndEnergy = fSystemAndEnergy;
  fCompactSystemAndEnergy.ReplaceAll(" ","");
  fCompactSystemAndEnergy.ReplaceAll(".","v");
  
  // Read bins for centrality, track pT, jet pT for energy-energy correlators, and track pT for energy-energy correlators from the card
  fnCentralityBins = fCard->GetNCentralityBins();
  
  // Centrality binning
  for(int iCentrality = 0; iCentrality <= fnCentralityBins; iCentrality++){
    fCentralityBinBorders[iCentrality] = fCard->GetLowBinBorderCentrality(iCentrality);
  }
  fLastLoadedCentralityBin = fnCentralityBins-1;
  
  // Remove centrality selection from pp data and local testing
  if(collisionSystem.Contains("pp") || collisionSystem.Contains("localTest")){
    fLastLoadedCentralityBin = 0;
    fCentralityBinBorders[0] = -0.5;
  }
  
}

/*
 * Copy constructor
 */
TriggerHistogramManager::TriggerHistogramManager(const TriggerHistogramManager& in) :
  fInputFile(in.fInputFile),
  fCard(in.fCard),
  fSystemAndEnergy(in.fSystemAndEnergy),
  fCompactSystemAndEnergy(in.fCompactSystemAndEnergy),
  fLoadEventInformation(in.fLoadEventInformation),
  fLoadJets(in.fLoadJets),
  fLoad2DHistograms(in.fLoad2DHistograms),
  fFirstLoadedCentralityBin(in.fFirstLoadedCentralityBin),
  fLastLoadedCentralityBin(in.fLastLoadedCentralityBin),
  fhVertexZ(in.fhVertexZ),
  fhVertexZWeighted(in.fhVertexZWeighted),
  fhEvents(in.fhEvents),
  fhCentrality(in.fhCentrality),
  fhCentralityWeighted(in.fhCentralityWeighted)
{
  // Copy constructor
  
  // Copy binning for centrality
  for(int iCentrality = 0; iCentrality < kMaxCentralityBins + 1; iCentrality++){
    fCentralityBinIndices[iCentrality] = in.fCentralityBinIndices[iCentrality];
    fCentralityBinBorders[iCentrality] = in.fCentralityBinBorders[iCentrality];
  }
  
  // Centrality loop
  for(int iCentrality = 0; iCentrality < kMaxCentralityBins; iCentrality++){
    for(int iDataLevel = 0; iDataLevel < TriggerHistograms::knDataLevels; iDataLevel++){
      for(int iTrigger = 0; iTrigger <= TriggerHistograms::knTriggerTypes; iTrigger++){
        for(int iJetType = 0; iJetType < knJetTypes; iJetType++){
          
          // Jet histograms
          fhJetPt[iJetType][iCentrality][iDataLevel][iTrigger] = in.fhJetPt[iJetType][iCentrality][iDataLevel][iTrigger];         // Jet pT histograms
          fhJetPhi[iJetType][iCentrality][iDataLevel][iTrigger] = in.fhJetPhi[iJetType][iCentrality][iDataLevel][iTrigger];       // Jet phi histograms
          fhJetEta[iJetType][iCentrality][iDataLevel][iTrigger] = in.fhJetEta[iJetType][iCentrality][iDataLevel][iTrigger];       // Jet eta histograms
          fhJetEtaPhi[iJetType][iCentrality][iDataLevel][iTrigger] = in.fhJetEtaPhi[iJetType][iCentrality][iDataLevel][iTrigger]; // 2D eta-phi histogram for jets
          
        } // Jet type loop
      } // Trigger selection loop
    } // Data level loop
  } // Centrality loop
}

/*
 * Destructor
 */
TriggerHistogramManager::~TriggerHistogramManager(){
  delete fCard;
}

/*
 * Load all the selected histograms from the inputfile
 */
void TriggerHistogramManager::LoadHistograms(){
  
  // Always load the number of events histogram
  fhEvents = (TH1D*) fInputFile->Get("nEvents");                           // Number of events surviving different event cuts
  
  // Load the event information histograms
  if(fLoadEventInformation){
    fhVertexZ = (TH1D*) fInputFile->Get("vertexZ");                        // Vertex z position
    fhVertexZWeighted = (TH1D*) fInputFile->Get("vertexZweighted");        // MC weighted vertex z position
    fhCentrality = (TH1D*) fInputFile->Get("centrality");                  // Centrality in all events
    fhCentralityWeighted = (TH1D*) fInputFile->Get("centralityWeighted");  // MC weighted centrality in all events
    fhPtHat = (TH1D*) fInputFile->Get("pthat");                            // pT hat for MC events
    fhPtHatWeighted = (TH1D*) fInputFile->Get("pthatWeighted");            // Weighted pT hat for MC events
  }
  
  // Load jet histograms
  LoadJetHistograms();
  
}

/*
 * Loader for jet histograms
 *
 * THnSparse for jets:
 *
 *   Histogram name: inclusiveJet
 *
 *     Axis index               Content of axis
 * --------------------------------------------------------
 *       Axis 0                     Jet pT
 *       Axis 1                     Jet phi
 *       Axis 2                     Jet eta
 *       Axis 3                    Centrality
 *       Axis 4     Data level (Reconstructed / Generator level)
 *       Axis 5                Trigger selection
 */
void TriggerHistogramManager::LoadJetHistograms(){
  
  // Define helper variables
  int duplicateRemoverCentrality = -1;
  int lowerCentralityBin = 0;
  int higherCentralityBin = 0;
  
  // Define arrays to help find the histograms
  int axisIndices[3] = {0};
  int lowLimits[3] = {0};
  int highLimits[3] = {0};
  
  int nAxes = 3;           // Number of constraining axes for this iteration
  
  
  for(int iCentralityBin = fFirstLoadedCentralityBin; iCentralityBin <= fLastLoadedCentralityBin; iCentralityBin++){
    
    // Select the bin indices
    lowerCentralityBin = fCentralityBinIndices[iCentralityBin];
    higherCentralityBin = fCentralityBinIndices[iCentralityBin+1]+duplicateRemoverCentrality;
    
    
    axisIndices[0] = 3; lowLimits[0] = lowerCentralityBin; highLimits[0] = higherCentralityBin;  // Centrality
    
    for(int iDataLevel = 0; iDataLevel < TriggerHistograms::knDataLevels; iDataLevel++){
      
      // Select the correct data level
      axisIndices[1] = 4; lowLimits[1] = iDataLevel+1; highLimits[1] = iDataLevel+1;
      
      for(int iTrigger = 0; iTrigger <= TriggerHistograms::knTriggerTypes; iTrigger++){
        
        // Apply trigger selection on top of the base trigger
        axisIndices[2] = 5; lowLimits[2] = iTrigger+1; highLimits[2] = iTrigger+1;
        
        for(int iJetType = 0; iJetType < knJetTypes; iJetType++){
          
          // Always load jet pT histograms
          fhJetPt[iJetType][iCentralityBin][iDataLevel][iTrigger] = FindHistogram(fInputFile,fJetHistogramName[iJetType],0,nAxes,axisIndices,lowLimits,highLimits);
          
          
          if(!fLoadJets) continue;  // Only load the remaining jet histograms if selected
          
          fhJetPhi[iJetType][iCentralityBin][iDataLevel][iTrigger] = FindHistogram(fInputFile,fJetHistogramName[iJetType],1,nAxes,axisIndices,lowLimits,highLimits);
          fhJetEta[iJetType][iCentralityBin][iDataLevel][iTrigger] = FindHistogram(fInputFile,fJetHistogramName[iJetType],2,nAxes,axisIndices,lowLimits,highLimits);
          if(fLoad2DHistograms) fhJetEtaPhi[iJetType][iCentralityBin][iDataLevel][iTrigger] = FindHistogram2D(fInputFile,fJetHistogramName[iJetType],1,2,nAxes,axisIndices,lowLimits,highLimits);
          
        } // Jet type loop
      } // Trigger selection loop
      
    } // Data level loop
    
  } // Loop over centrality bins
}

/*
 * Extract a 2D histogram from a given centrality bin from THnSparseD
 *
 *  Arguments:
 *   TFile *inputFile = Inputfile containing the THnSparse to be read
 *   const char *name = Name of the THnSparse that is read
 *   int xAxis = Index for the axis in THnSparse that is projected to x-axis for TH2D
 *   int yAxis = Index for the axis in THnSparse that is projected to y-axis for TH2D
 *   int nAxes = Number of axes that are restained for the projection
 *   int *axisNumber = Indices for the axes in THnSparse that are used as a restriction for the projection
 *   int *lowBinIndex = Indices of the lowest considered bins in the restriction axis
 *   int *highBinIndex = Indices of the highest considered bins in the restriction axis
 *   const bool normalizeToBinWidth = Flag for normalizing the projected histogram to the bin width
 */
TH2D* TriggerHistogramManager::FindHistogram2D(TFile *inputFile, const char *name, int xAxis, int yAxis, int nAxes, int *axisNumber, int *lowBinIndex, int *highBinIndex, const bool normalizeToBinWidth){
  
  // Read the histogram with the given name from the file
  THnSparseD *histogramArray = (THnSparseD*) inputFile->Get(name);
  
  // If cannot find histogram, inform that it could not be found and return null
  if(histogramArray == nullptr){
    cout << "Could not find " << name << ". Skipping loading this histogram." << endl;
    return NULL;
  }
  
  // Apply the restrictions in the set of axes
  for(int i = 0; i < nAxes; i++) histogramArray->GetAxis(axisNumber[i])->SetRange(lowBinIndex[i],highBinIndex[i]);
  
  // Create a unique name for eeach histogram that is read from the file
  TString newName = histogramArray->GetName();
  for(int iBinIndex = 0; iBinIndex < nAxes; iBinIndex++){
    newName.Append(Form("_%d=%d-%d",axisNumber[iBinIndex],lowBinIndex[iBinIndex],highBinIndex[iBinIndex]));
  }
  
  // Project out the histogram and give it the created unique name
  TH2D *projectedHistogram = (TH2D*) histogramArray->Projection(yAxis,xAxis);
  projectedHistogram->SetName(newName.Data());
  
  // Apply bin width normalization to the projected histogram
  if(normalizeToBinWidth) projectedHistogram->Scale(1.0,"width");
  
  // Return the projected histogram
  return projectedHistogram;
}

/*
 * Extract a 2D histogram from a given centrality bin from THnSparseD
 *
 *  Arguments:
 *   TFile *inputFile = Inputfile containing the THnSparse to be read
 *   const char *name = Name of the THnSparse that is read
 *   int xAxis = Index for the axis in THnSparse that is projected to x-axis for TH2D
 *   int yAxis = Index for the axis in THnSparse that is projected to y-axis for TH2D
 *   int restrictionAxis = Index for the axis in THnSparse that is used as a restriction for the projection
 *   int lowBinIndex = Index of the lowest considered bin in the restriction axis
 *   int highBinIndex = Index of the highest considered bin in the restriction axis
 *   int restrictionAxis2 = Index for the axis in THnSparse that is used as a second restriction for the projection
 *   int lowBinIndex2 = Index of the lowest considered bin in the second restriction axis
 *   int highBinIndex2 = Index of the highest considered bin in the second restriction axis
 *   const bool normalizeToBinWidth = Flag for normalizing the projected histogram to the bin width
 */
TH2D* TriggerHistogramManager::FindHistogram2D(TFile *inputFile, const char *name, int xAxis, int yAxis, int restrictionAxis, int lowBinIndex, int highBinIndex, int restrictionAxis2, int lowBinIndex2, int highBinIndex2, const bool normalizeToBinWidth){
  int restrictionAxes[2] = {restrictionAxis,restrictionAxis2};
  int lowBinIndices[2] = {lowBinIndex,lowBinIndex2};
  int highBinIndices[2] = {highBinIndex,highBinIndex2};
  int nAxes = 2;
  if(highBinIndex2 == 0 && lowBinIndex2 == 0) nAxes = 1;
  return FindHistogram2D(inputFile,name,xAxis,yAxis,nAxes,restrictionAxes,lowBinIndices,highBinIndices,normalizeToBinWidth);
}

/*
 * Extract a histogram with given restrictions on other axes in THnSparse
 *
 *  Arguments:
 *   TFile *inputFile = Inputfile containing the THnSparse to be read
 *   const char *name = Name of the THnSparse that is read
 *   int xAxis = Index for the axis in THnSparse that is projected to x-axis for TH1D
 *   int nAxes = Number of axes that are restained for the projection
 *   int *axisNumber = Indices for the axes in THnSparse that are used as a restriction for the projection
 *   int *lowBinIndex = Indices of the lowest considered bins in the restriction axis
 *   int *highBinIndex = Indices of the highest considered bins in the restriction axis
 *   const bool normalizeToBinWidth = Flag for normalizing the projected histogram to the bin width
 */
TH1D* TriggerHistogramManager::FindHistogram(TFile *inputFile, const char *name, int xAxis, int nAxes, int *axisNumber, int *lowBinIndex, int *highBinIndex, const bool normalizeToBinWidth){
  
  // Read the histogram with the given name from the file
  THnSparseD *histogramArray = (THnSparseD*) inputFile->Get(name);
  
  // If cannot find histogram, inform that it could not be found and return null
  if(histogramArray == nullptr){
    cout << "Could not find " << name << ". Skipping loading this histogram." << endl;
    return NULL;
  }
  
  // Apply the restrictions in the set of axes
  for(int i = 0; i < nAxes; i++) histogramArray->GetAxis(axisNumber[i])->SetRange(lowBinIndex[i],highBinIndex[i]);
  
  // Create a unique name for each histogram that is read from the file
  TString newName = histogramArray->GetName();
  for(int iBinIndex = 0; iBinIndex < nAxes; iBinIndex++){
    newName.Append(Form("_%d=%d-%d",axisNumber[iBinIndex],lowBinIndex[iBinIndex],highBinIndex[iBinIndex]));
  }
  
  // Project out the histogram and give it the created unique name
  TH1D *projectedHistogram = NULL;
  
  // Check that we are not trying to project a non-existing axis
  if(xAxis < histogramArray->GetNdimensions()){
    projectedHistogram = (TH1D*) histogramArray->Projection(xAxis);
    projectedHistogram->SetName(newName.Data());
  
    // Apply bin width normalization to the projected histogram
    if(normalizeToBinWidth) projectedHistogram->Scale(1.0,"width");
  }
  
  // Return the projected histogram
  return projectedHistogram;
}

/*
 * Extract a histogram from a given centrality bin from THnSparseD
 *
 *  Arguments:
 *   TFile *inputFile = Inputfile containing the THnSparse to be read
 *   const char *name = Name of the THnSparse that is read
 *   int xAxis = Index for the axis in THnSparse that is projected to TH1D
 *   int restrictionAxis = Index for the axis in THnSparse that is used as a restriction for the projection
 *   int lowBinIndex = Index of the lowest considered bin in the restriction axis
 *   int highBinIndex = Index of the highest considered bin in the restriction axis
 *   int restrictionAxis2 = Index for the axis in THnSparse that is used as a second restriction for the projection
 *   int lowBinIndex2 = Index of the lowest considered bin in the second restriction axis
 *   int highBinIndex2 = Index of the highest considered bin in the second restriction axis
 *   const bool normalizeToBinWidth = Flag for normalizing the projected histogram to the bin width
 */
TH1D* TriggerHistogramManager::FindHistogram(TFile *inputFile, const char *name, int xAxis, int restrictionAxis, int lowBinIndex, int highBinIndex, int restrictionAxis2, int lowBinIndex2, int highBinIndex2, const bool normalizeToBinWidth){
  int restrictionAxes[2] = {restrictionAxis,restrictionAxis2};
  int lowBinIndices[2] = {lowBinIndex,lowBinIndex2};
  int highBinIndices[2] = {highBinIndex,highBinIndex2};
  int nAxes = 2;
  if(highBinIndex2 == 0 && lowBinIndex2 == 0) nAxes = 1;
  return FindHistogram(inputFile,name,xAxis,nAxes,restrictionAxes,lowBinIndices,highBinIndices,normalizeToBinWidth);
}

/*
 * Write all the loaded histograms into a file
 *
 *  const char* fileName = Name of the file to which the histograms are written
 *  const char* fileOption = Option given to the file when it is loaded
 */
void TriggerHistogramManager::Write(const char* fileName, const char* fileOption){
  
  // Create the output file
  TFile *outputFile = new TFile(fileName,fileOption);
  
  // Helper variable for renaming the saved histograms
  TString histogramNamer;
  
  // Write the event information histograms to the output file
  if(fLoadEventInformation){
    fhEvents->Write("",TObject::kOverwrite);             // Number of events surviving different event cuts
    fhVertexZ->Write("",TObject::kOverwrite);            // Vertex z position
    fhVertexZWeighted->Write("",TObject::kOverwrite);    // MC weighted vertex z position
    fhCentrality->Write("",TObject::kOverwrite);         // Centrality in all events
    fhCentralityWeighted->Write("",TObject::kOverwrite); // MC weighted centrality in all events
    fhPtHat->Write("",TObject::kOverwrite);              // pT hat for MC events (only meaningful for MC)
    fhPtHatWeighted->Write("",TObject::kOverwrite);      // Weighted pT hat distribution (only meaningful for MC)
  }
 
  // Write the jet histograms to the output file
  WriteJetHistograms();
  
  // Write the card to the output file if it is not already written
  if(!gDirectory->GetDirectory("JCard")) fCard->Write(outputFile);
  
  // Close the file after everything is written
  outputFile->Close();
  
  // Delete the outputFile object
  delete outputFile;
  
}

/*
 * Write the jet histograms to the file that is currently open
 */
void TriggerHistogramManager::WriteJetHistograms(){
  
  // Helper variable for histogram naming
  TString histogramNamer;
  TString triggerNamer;
  TriggerHistograms *triggerProvider = new TriggerHistograms();
  
  // Write the jet histograms to the output file
  if(!fLoadJets) return;  // Only write the jet histograms if they are loaded
  
  for(int iJetType = 0; iJetType < knJetTypes; iJetType++){
    for(int iDataLevel = 0; iDataLevel < TriggerHistograms::knDataLevels; iDataLevel++){
      
      // Create a directory for the histograms if it does not already exist
      histogramNamer = Form("%s%s", fJetHistogramName[iJetType], fDataLevelName[iDataLevel]);
      if(!gDirectory->GetDirectory(histogramNamer)) gDirectory->mkdir(histogramNamer);
      gDirectory->cd(histogramNamer);
      
      for(int iCentralityBin = fFirstLoadedCentralityBin; iCentralityBin <= fLastLoadedCentralityBin; iCentralityBin++){
        for(int iTrigger = 0; iTrigger <= TriggerHistograms::knTriggerTypes; iTrigger++){
          
          // Check that the histograms are actually there before trying to save them.
          if(fhJetPt[iJetType][iCentralityBin][iDataLevel][iTrigger] == NULL) {
            cout << "Could not find histograms of type " << fJetHistogramName[iJetType] << " to write. Will skip writing these." << endl;
            continue;
          }
          
          triggerNamer = "_" + triggerProvider->GetTriggerName(iTrigger);
          if(iTrigger == TriggerHistograms::knTriggerTypes) triggerNamer = "";
          
          // Jet pT
          histogramNamer = Form("%sPt%s_C%d",fJetHistogramName[iJetType], triggerNamer.Data(), iCentralityBin);
          if(fhJetPt[iJetType][iCentralityBin][iDataLevel][iTrigger]) fhJetPt[iJetType][iCentralityBin][iDataLevel][iTrigger]->Write(histogramNamer.Data(), TObject::kOverwrite);
          
          // Jet phi
          histogramNamer = Form("%sPhi%s_C%d",fJetHistogramName[iJetType], triggerNamer.Data(), iCentralityBin);
          if(fhJetPhi[iJetType][iCentralityBin][iDataLevel][iTrigger]) fhJetPhi[iJetType][iCentralityBin][iDataLevel][iTrigger]->Write(histogramNamer.Data(), TObject::kOverwrite);
          
          // Jet eta
          histogramNamer = Form("%sEta%s_C%d",fJetHistogramName[iJetType], triggerNamer.Data(), iCentralityBin);
          if(fhJetEta[iJetType][iCentralityBin][iDataLevel][iTrigger]) fhJetEta[iJetType][iCentralityBin][iDataLevel][iTrigger]->Write(histogramNamer.Data(), TObject::kOverwrite);
          
          // Jet eta-phi
          histogramNamer = Form("%sEtaPhi%s_C%d",fJetHistogramName[iJetType], triggerNamer.Data(), iCentralityBin);
          if(fLoad2DHistograms && fhJetEtaPhi[iJetType][iCentralityBin][iDataLevel][iTrigger]) fhJetEtaPhi[iJetType][iCentralityBin][iDataLevel][iTrigger]->Write(histogramNamer.Data(), TObject::kOverwrite);
          
        } // Trigger selection loop
      } // Loop over centrality bins
      
      // Return back to main directory
      gDirectory->cd("../");
      
    } // Data level loop
  } // Jet type loop
   
  delete triggerProvider;
  
}

/*
 * Load the selected histograms from a file containing readily processed histograms
 */
void TriggerHistogramManager::LoadProcessedHistograms(){
  
  // Helper variable for finding names of loaded histograms
  TString histogramNamer;
  TString triggerNamer;
  TriggerHistograms *triggerProvider = new TriggerHistograms();
  
  // Always load the number of events histogram
  fhEvents = (TH1D*) fInputFile->Get("nEvents");                           // Number of events surviving different event cuts
  
  // Load the event information histograms
  if(fLoadEventInformation){
    fhVertexZ = (TH1D*) fInputFile->Get("vertexZ");                        // Vertex z position
    fhVertexZWeighted = (TH1D*) fInputFile->Get("vertexZweighted");        // MC weighted vertex z position
    fhCentrality = (TH1D*) fInputFile->Get("centrality");                  // Centrality in all events
    fhCentralityWeighted = (TH1D*) fInputFile->Get("centralityWeighted");  // MC weighted centrality in all events
    fhPtHat = (TH1D*) fInputFile->Get("pthat");                            // pT hat for MC events
    fhPtHatWeighted = (TH1D*) fInputFile->Get("pthatWeighted");            // Weighted pT hat for MC events
  }
  
  // Load the jet histograms from the input file
  for(int iJetType = 0; iJetType < knJetTypes; iJetType++){
    for(int iDataLevel = 0; iDataLevel < TriggerHistograms::knDataLevels; iDataLevel++){
      for(int iCentralityBin = fFirstLoadedCentralityBin; iCentralityBin <= fLastLoadedCentralityBin; iCentralityBin++){
        
        for(int iTrigger = 0; iTrigger <= TriggerHistograms::knTriggerTypes; iTrigger++){
          
          triggerNamer = "_" + triggerProvider->GetTriggerName(iTrigger);
          if(iTrigger == TriggerHistograms::knTriggerTypes) triggerNamer = "";
          
          // Always load jet pT histograms
          histogramNamer = Form("%s%s/%sPt%s_C%d", fJetHistogramName[iJetType], fDataLevelName[iDataLevel], fJetHistogramName[iJetType], triggerNamer.Data(), iCentralityBin);
          fhJetPt[iJetType][iCentralityBin][iDataLevel][iTrigger] = (TH1D*) fInputFile->Get(histogramNamer.Data());
          
          if(!fLoadJets) continue;  // Only load the loaded the selected histograms
          
          // Jet phi
          histogramNamer = Form("%s%s/%sPhi%s_C%d", fJetHistogramName[iJetType], fDataLevelName[iDataLevel], fJetHistogramName[iJetType], triggerNamer.Data(), iCentralityBin);
          fhJetPhi[iJetType][iCentralityBin][iDataLevel][iTrigger] = (TH1D*) fInputFile->Get(histogramNamer.Data());
          
          // Jet eta
          histogramNamer = Form("%s%s/%sEta%s_C%d", fJetHistogramName[iJetType], fDataLevelName[iDataLevel], fJetHistogramName[iJetType], triggerNamer.Data(), iCentralityBin);
          fhJetEta[iJetType][iCentralityBin][iDataLevel][iTrigger] = (TH1D*) fInputFile->Get(histogramNamer.Data());
          
          // Jet eta-phi
          histogramNamer = Form("%s%s/%sEtaPhi%s_C%d", fJetHistogramName[iJetType], fDataLevelName[iDataLevel], fJetHistogramName[iJetType], triggerNamer.Data(), iCentralityBin);
          if(fLoad2DHistograms) fhJetEtaPhi[iJetType][iCentralityBin][iDataLevel][iTrigger] = (TH2D*) fInputFile->Get(histogramNamer.Data());
          
        } // Jet trigger selection loop
        
      } // Loop over centrality bins
    } // Data level loop
  } // Jet type loop
  
  delete triggerProvider;
}

/*
 * Read the bin indices for given bin borders
 *
 *  Arguments:
 *   const char* histogramName = Name of the histogram from which the bin indices are searched
 *   const int nBins = Number of bins for the indices
 *   int *binIndices = Array of integers to be filled with bin index information read from the file
 *   const double *binBorders = Array for bin borders that are searched from the file
 *   const int iAxis = Index of the axis used for reading bin indices
 */
void TriggerHistogramManager::SetBinIndices(const char* histogramName, const int nBins, int *binIndices, const double *binBorders, const int iAxis){
  TH1D* hBinner = FindHistogram(fInputFile,histogramName,iAxis,0,0,0);
  for(int iBin = 0; iBin < nBins+1; iBin++){
    binIndices[iBin] = hBinner->GetXaxis()->FindBin(binBorders[iBin]);
  }
}

/*
 * Read the bin indices for given bin borders
 *
 *  Arguments:
 *   const char* histogramName = Name of the histogram from which the bin indices are searched
 *   const int nBins = Number of bins for the indices
 *   double *copyBinBorders = Array to which a copy of bin borders is made
 *   int *binIndices = Array of integers to be filled with bin index information read from the file
 *   const double *binBorders = Array for bin borders that are searched from the file
 *   const int iAxis = Index of the axis used for reading bin indices
 *   const bool setIndices = true: Set both bin indices and bin borders, false: Set only bin borders
 */
void TriggerHistogramManager::SetBinBordersAndIndices(const char* histogramName, const int nBins, double *copyBinBorders, int *binIndices, const double *binBorders, const int iAxis, const bool setIndices){
  TH1D* hBinner;
  if(setIndices) hBinner = FindHistogram(fInputFile,histogramName,iAxis,0,0,0);
  for(int iBin = 0; iBin < nBins+1; iBin++){
    copyBinBorders[iBin] = binBorders[iBin];
    if(setIndices) binIndices[iBin] = hBinner->GetXaxis()->FindBin(binBorders[iBin]);
  }
}

/*
 * Set up generic bin borders and indices according to provided bin borders
 *
 *  const bool readBinsFromFile = True: Disregard given bins ans use the ones in fCard. False: Use given bins
 *  const char* histogramName = Name of the histogram from which the bin indices are searched
 *  const int iAxis = Axis from which the set indices can be found
 *  int nSetBins = Number of bins that is set
 *  double* setBinBorders = Bin borders that are set
 *  int* setBinIndices = Bin indices that are set
 *  const int nBins = New number of bins that is given
 *  const double *binBorders = New bin borders that are given
 *  const char* errorMessage = Type of the set bins to be printed in possible error message
 *  const int maxBins = Maximum number of allowed bins of this type
 *  const bool setIndices = Set the bin indices in THnSparse
 */
void TriggerHistogramManager::SetGenericBins(const bool readBinsFromFile, const char* histogramName, const int iAxis, int nSetBins, double* setBinBorders, int* setBinIndices, const int nBins, const double *binBorders, const char* errorMessage, const int maxBins, const bool setIndices){
  
  // If bins are read from file, do not use the given bin borders
  if(readBinsFromFile){
    if(setIndices) SetBinIndices(histogramName, nSetBins, setBinIndices, setBinBorders, iAxis);
  } else { // If the given bin borders are use, update the number of bins and bin borders according to input
    if(nBins <= maxBins){
      nSetBins = nBins;
      SetBinBordersAndIndices(histogramName, nSetBins, setBinBorders, setBinIndices, binBorders, iAxis, setIndices);
    } else {
      cout << "Error! Too many " << errorMessage << " bins given. Maximum number is " << maxBins << ". Will not set bins." << endl;
    }
  }
}

/*
 * Set up centrality bin borders and indices according to provided bin borders
 *
 *  const bool readBinsFromFile = True: Disregard given bins ans use the ones in fCard. False: Use given bins
 *  const int nBins = Number of given centrality bins
 *  const double *binBorders = New bin borders for centrality
 *  const bool setIndices = Set the bin indices in THnSparse
 */
void TriggerHistogramManager::SetCentralityBins(const bool readBinsFromFile, const int nBins, const double *binBorders, const bool setIndices){
  
  SetGenericBins(readBinsFromFile, fJetHistogramName[0], 3, fnCentralityBins, fCentralityBinBorders, fCentralityBinIndices, nBins, binBorders, "centrality", kMaxCentralityBins, setIndices);
  
}

// Setter for loading event information
void TriggerHistogramManager::SetLoadEventInformation(const bool loadOrNot){
  fLoadEventInformation = loadOrNot;
}

// Setter for loading jet histograms
void TriggerHistogramManager::SetLoadJetHistograms(const bool loadOrNot){
  fLoadJets = loadOrNot;
}

 // Setter for loading two-dimensional histograms
void TriggerHistogramManager::SetLoad2DHistograms(const bool loadOrNot){
  fLoad2DHistograms = loadOrNot;
}

// Setter for loaded centrality bins
void TriggerHistogramManager::SetCentralityBinRange(const int first, const int last){
  fFirstLoadedCentralityBin = first;
  fLastLoadedCentralityBin = last;
  
  // Sanity check for centrality bins
  BinSanityCheck(fnCentralityBins,fFirstLoadedCentralityBin,fLastLoadedCentralityBin);
}

// Sanity check for set bins
void TriggerHistogramManager::BinSanityCheck(const int nBins, int& first, int& last){
  if(first < 0) first = 0;
  if(last < first) last = first;
  if(last > nBins-1) last = nBins-1;
}

// Sanity check for input bin index
int TriggerHistogramManager::BinIndexCheck(const int nBins, const int binIndex) const{
  if(binIndex < 0) return 0;
  if(binIndex > nBins-1) return nBins-1;
  return binIndex;
}

// Getter for the number of centrality bins
int TriggerHistogramManager::GetNCentralityBins() const{
  return fnCentralityBins;
}

// Getter for the jet histogram name
const char* TriggerHistogramManager::GetJetHistogramName(const int iJetType) const{
  return fJetHistogramName[iJetType];
}

// Getter for name suitable for x-axis in a given jet histogram
const char* TriggerHistogramManager::GetJetAxisName(const int iJetType) const{
  return fJetAxisName[iJetType];
}

// Getter for collision system
TString TriggerHistogramManager::GetSystem() const{
  return fCard->GetDataType();
}

// Getter for i:th centrality bin border
double TriggerHistogramManager::GetCentralityBinBorder(const int iCentrality) const{
  return fCentralityBinBorders[iCentrality];
}

// Getters for event information histograms

// Getter for z-vertex histogram
TH1D* TriggerHistogramManager::GetHistogramVertexZ() const{
  return fhVertexZ;
}

// Getter for z-vertex histogram
TH1D* TriggerHistogramManager::GetHistogramVertexZWeighted() const{
  return fhVertexZWeighted;
}

// Getter for histogram for number of events surviving different event cuts
TH1D* TriggerHistogramManager::GetHistogramEvents() const{
  return fhEvents;
}

// Getter for centrality histogram in all events
TH1D* TriggerHistogramManager::GetHistogramCentrality() const{
  return fhCentrality;
}

// Getter for weighted centrality histogram in all events
TH1D* TriggerHistogramManager::GetHistogramCentralityWeighted() const{
  return fhCentralityWeighted;
}

// Getters for jet histograms

// Getter for jet pT histograms
TH1D* TriggerHistogramManager::GetHistogramJetPt(const int iJetType, int iCentrality, const int iDataLevel, const int iTrigger) const{
  if(fCard->GetDataType().Contains("pp",TString::kIgnoreCase)) iCentrality = 0;  // No centrality selection for pp
  return fhJetPt[iJetType][iCentrality][iDataLevel][iTrigger];
}

// Getter for jet phi histograms
TH1D* TriggerHistogramManager::GetHistogramJetPhi(const int iJetType, int iCentrality, const int iDataLevel, const int iTrigger) const{
  if(fCard->GetDataType().Contains("pp",TString::kIgnoreCase)) iCentrality = 0;  // No centrality selection for pp
  return fhJetPhi[iJetType][iCentrality][iDataLevel][iTrigger];
}

// Getter for jet eta histograms
TH1D* TriggerHistogramManager::GetHistogramJetEta(const int iJetType, int iCentrality, const int iDataLevel, const int iTrigger) const{
  if(fCard->GetDataType().Contains("pp",TString::kIgnoreCase)) iCentrality = 0;  // No centrality selection for pp
  return fhJetEta[iJetType][iCentrality][iDataLevel][iTrigger];
}

// Getter for 2D eta-phi histogram for jets
TH2D* TriggerHistogramManager::GetHistogramJetEtaPhi(const int iJetType, int iCentrality, const int iDataLevel, const int iTrigger) const{
  if(fCard->GetDataType().Contains("pp",TString::kIgnoreCase)) iCentrality = 0;  // No centrality selection for pp
  return fhJetEtaPhi[iJetType][iCentrality][iDataLevel][iTrigger];
}

// Get the first loaded centrality bin
int TriggerHistogramManager::GetFirstCentralityBin() const{
  return fFirstLoadedCentralityBin;
}

// Get the last loaded centrality bin
int TriggerHistogramManager::GetLastCentralityBin() const{
  return fLastLoadedCentralityBin;
}

// Getter for the number of events passing the cuts
int TriggerHistogramManager::GetNEvents() const{
  return fhEvents->GetBinContent(fhEvents->FindBin(TriggerHistograms::kVzCut));
}

// Getter for the JCard
TriggerCard* TriggerHistogramManager::GetCard() const{
  return fCard;
}

// Getter for integral over inclusive jet pT. Include the overflow bin in the integral.
double TriggerHistogramManager::GetJetPtIntegral(const int iJetType, const int iCentrality, const int iDataLevel, const int iTrigger) const{
  return fhJetPt[iJetType][iCentrality][iDataLevel][iTrigger]->Integral(1,fhJetPt[iJetType][iCentrality][iDataLevel][iTrigger]->GetNbinsX()+1,"width");
}

/*
 * Getter for integral over inclusive jet pT over specified range
 *
 *  const int iCentrality = Centrality bin
 *  const double minPt = Lower pT range for integral calculation
 *  const double maxPt = Higher pT range for integral calculation
 */
double TriggerHistogramManager::GetJetPtIntegral(const int iJetType, const int iCentrality, const int iDataLevel, const int iTrigger, const double minPt, const double maxPt) const{
  return fhJetPt[iJetType][iCentrality][iDataLevel][iTrigger]->Integral(fhJetPt[iJetType][iCentrality][iDataLevel][iTrigger]->FindBin(minPt+0.001), fhJetPt[iJetType][iCentrality][iDataLevel][iTrigger]->FindBin(maxPt-0.001), "width");
}
