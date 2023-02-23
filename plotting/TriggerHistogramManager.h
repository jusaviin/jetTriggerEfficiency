#ifndef TRIGGERHISTOGRAMMANAGER_H
#define TRIGGERHISTOGRAMMANAGER_H

// Root includes
#include <TFile.h>
#include <TDirectory.h>
#include <TH1.h>
#include <TH2.h>
#include <TString.h>
#include <TLegend.h>
#include <TStyle.h>

// Own includes
#include "TriggerCard.h"
#include "../src/TriggerHistograms.h"

/*
 * Class for managing the histograms produced in the trigger analysis
 */
class TriggerHistogramManager {
 
public:
  
  // Dimensions for histogram arrays
  static const int kMaxCentralityBins = 5;       // Maximum allowed number of centrality bins
  
private:
  
  // Naming for jet histograms
  const char* fJetHistogramName = "inclusiveJet";
  const char* fJetAxisName = "Jet";
  
  // Naming for data levels
  const char* fDataLevelName[TriggerHistograms::knDataLevels] = {"", "GeneratorLevel"};
  
public:
  
  TriggerHistogramManager();                                    // Default constructor
  TriggerHistogramManager(TFile *inputFile);                    // Constructor
  TriggerHistogramManager(TFile *inputFile, TriggerCard *card);   // Constructor with card
  TriggerHistogramManager(const TriggerHistogramManager& in);     // Copy constructor
  ~TriggerHistogramManager();                                   // Destructor
  
  void LoadHistograms();          // Load the histograms from the inputfile
  void Write(const char* fileName, const char* fileOption);          // Write all the loaded histograms into a file
  void LoadProcessedHistograms(); // Load processed histograms from the inputfile
  
  // Setters for binning information
  void SetCentralityBins(const bool readBinsFromFile, const int nBins, const double *binBorders, bool setIndices = true); // Set up centrality bin indices according to provided bin borders
  
  // Setters for event information and dijets
  void SetLoadEventInformation(const bool loadOrNot); // Setter for loading event information
  
  // Setters for jets
  void SetLoadJetHistograms(const bool loadOrNot);        // Setter for loading all jet histograms
  
  // Setter for loading additional histograms
  void SetLoad2DHistograms(const bool loadOrNot);           // Setter for loading two-dimensional histograms
  
  // Setters for ranges for different bins
  void SetCentralityBinRange(const int first, const int last); // Setter for centrality bin range
  
  // Getters for number of bins in histograms
  int GetNCentralityBins() const;  // Getter for the number of centrality bins
  double GetCentralityBinBorder(const int iCentrality) const;  // Getter for i:th centrality bin border
  
  // Getters for histogram and axis naming
  const char* GetJetHistogramName() const; // Getter for the jet histogram name
  const char* GetJetAxisName() const;      // Getter for name suitable for x-axis in a jet histogram
  
  TString GetSystem() const;  // Getter for collision system
  
  // Getters for event information histograms
  TH1D* GetHistogramVertexZ() const;            // Getter for z-vertex histogram
  TH1D* GetHistogramVertexZWeighted() const;    // Getter for weighted z-vertex histogram
  TH1D* GetHistogramEvents() const;             // Getter for histogram for number of events surviving different event cuts
  TH1D* GetHistogramCentrality() const;         // Getter for centrality histogram in all events
  TH1D* GetHistogramCentralityWeighted() const; // Getter for weighted centrality histogram in all events
  
  // Getters for jet histograms
  TH1D* GetHistogramJetPt(int iCentrality, const int iDataLevel, const int iTrigger) const;     // Jet pT histograms
  TH1D* GetHistogramJetPhi(int iCentrality, const int iDataLevel, const int iTrigger) const;    // Jet phi histograms
  TH1D* GetHistogramJetEta(int iCentrality, const int iDataLevel, const int iTrigger) const;    // Jet eta histograms
  TH2D* GetHistogramJetEtaPhi(int iCentrality, const int iDataLevel, const int iTrigger) const; // 2D eta-phi histogram for jets
  
  // Getters for the loaded centrality bins
  int GetFirstCentralityBin() const;  // Get the first loaded centrality bin
  int GetLastCentralityBin() const;   // Get the last loaded centrality bin
  
  // Getters for normalization information
  int GetNEvents() const;                      // Getter for the number of events passing the cuts
  double GetJetPtIntegral(const int iCentrality, const int iDataLevel, const int iTrigger) const; // Getter for integral over inclusive jet pT in a given centrality
  double GetJetPtIntegral(int iCentrality, const int iDataLevel, const int iTrigger, const double minPt, const double maxPt) const; // Getter for integral over inclusive jet pT in a given pT range within a given centrality bin
  
  // Getter for the card
  TriggerCard* GetCard() const;  // Getter for the JCard
  
private:
  
  // Data members
  TFile *fInputFile;                  // File from which the histograms are read
  TriggerCard *fCard;                   // Card inside the data file for binning, cut collision system etc. information
  TString fSystemAndEnergy;           // Collision system (pp,PbPb,pp MC,PbPb MC,localTest) and energy
  TString fCompactSystemAndEnergy;    // Same a before but without white spaces and dots
  
  // ==============================================
  // ======== Flags for histograms to load ========
  // ==============================================
  
  bool fLoadEventInformation;                              // Load the event information histograms
  bool fLoadJets;                                          // Load the jet histograms
  bool fLoad2DHistograms;                                  // Load also two-dimensional (eta,phi) histograms
  
  // ==============================================
  // ======== Ranges of histograms to load ========
  // ==============================================
  
  int fFirstLoadedCentralityBin;  // First centrality bin that is loaded
  int fLastLoadedCentralityBin;   // Last centrality bin that is loaded
  
  // =============================================
  // ============ Binning information ============
  // =============================================
  int fCentralityBinIndices[kMaxCentralityBins+1];    // Indices for centrality bins in centrality binned histograms
  double fCentralityBinBorders[kMaxCentralityBins+1]; // Centrality bin borders, from which bin indices are obtained
  int fnCentralityBins;                               // Number of centrality bins in the JCard of the data file
  
  // =============================================
  // ===== Histograms for the dijet analysis =====
  // =============================================
  
  // Event information histograms
  TH1D *fhVertexZ;            // Vertex z position
  TH1D *fhVertexZWeighted;    // Weighted vertex z-position (only meaningfull for MC)
  TH1D *fhEvents;             // Number of events surviving different event cuts
  TH1D *fhCentrality;         // Centrality of all events
  TH1D *fhCentralityWeighted; // Weighted centrality distribution in all events (only meaningful for MC)
  TH1D *fhPtHat;              // pT hat for MC events (only meaningful for MC)
  TH1D *fhPtHatWeighted;      // Weighted pT hat distribution (only meaningful for MC)
  
  // Histograms for jets
  TH1D *fhJetPt[kMaxCentralityBins][TriggerHistograms::knDataLevels][TriggerHistograms::knTriggerTypes+1];      // Jet pT histograms
  TH1D *fhJetPhi[kMaxCentralityBins][TriggerHistograms::knDataLevels][TriggerHistograms::knTriggerTypes+1];     // Jet phi histograms
  TH1D *fhJetEta[kMaxCentralityBins][TriggerHistograms::knDataLevels][TriggerHistograms::knTriggerTypes+1];     // Jet eta histograms
  TH2D *fhJetEtaPhi[kMaxCentralityBins][TriggerHistograms::knDataLevels][TriggerHistograms::knTriggerTypes+1];  // 2D eta-phi histogram for jets
  
  // Private methods
  void InitializeFromCard(); // Initialize several member variables from TriggerCard
  
  // Binning related methods
  void SetBinIndices(const char* histogramName, const int nBins, int *binIndices, const double *binBorders, const int iAxis); // Read the bin indices for given bin borders
  void SetBinBordersAndIndices(const char* histogramName, const int nBins, double *copyBinBorders, int *binIndices, const double *binBorders, const int iAxis, const bool setIndices); // Read the bin indices for given bin borders
  
  // Finders for histograms with different amount of restrictions
  TH2D* FindHistogram2D(TFile *inputFile, const char *name, int xAxis, int yAxis, int nAxes, int *axisNumber, int *lowBinIndex, int *highBinIndex, const bool normalizeToBinWidth = true); // Extract a 2D histogram using given axis restrictions from THnSparseD
  TH2D* FindHistogram2D(TFile *inputFile, const char *name, int xAxis, int yAxis, int restrictionAxis, int lowBinIndex, int highBinIndex, int restrictionAxis2 = 0, int lowBinIndex2 = 0, int highBinIndex2 = 0, const bool normalizeToBinWidth = true); // Extract a 2D histogram using given axis restrictions from THnSparseD
  TH1D* FindHistogram(TFile *inputFile, const char *name, int xAxis, int nAxes, int *axisNumber, int *lowBinIndex, int *highBinIndex, const bool normalizeToBinWidth = true); // Extract a histogram using given axis restrictions from THnSparseD
  TH1D* FindHistogram(TFile *inputFile, const char *name, int xAxis, int restrictionAxis, int lowBinIndex, int highBinIndex, int restrictionAxis2 = 0, int lowBinIndex2 = 0, int highBinIndex2 = 0, const bool normalizeToBinWidth = true); // Extract a histogram using given axis restrictions from THnSparseD
  
  // Loaders for different groups of histograms
  void LoadJetHistograms(); // Loader for jet histograms
  
  // Generic setter for bin indice and borders
  void SetGenericBins(const bool readBinsFromFile, const char* histogramName, const int iAxis, int nSetBins, double* setBinBorders, int* setBinIndices, const int nBins, const double *binBorders, const char* errorMessage, const int maxBins, const bool setIndices); // Generic bin setter
  
  // Methods for binning
  void BinSanityCheck(const int nBins, int& first, int& last); // Sanity check for given binning
  int BinIndexCheck(const int nBins, const int binIndex) const; // Check that given index is in defined range
  
  // Methods for histogram writing
  void WriteJetHistograms();                        // Write the jet histograms to the file that is currently open
  
};

#endif
