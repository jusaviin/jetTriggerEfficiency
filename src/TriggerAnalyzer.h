// Class for the main analysis algorithms for leading-subleading jet analysis

#ifndef TRIGGERANALYZER_H
#define TRIGGERANALYZER_H

// C++ includes
#include <vector>
#include <bitset>
#include <assert.h>   // Standard c++ debugging tool. Terminates the program if expression given evaluates to 0.
#include <tuple>      // For returning several arguments in a transparent manner
#include <fstream>
#include <string>

// Root includes
#include <TString.h>
#include <TRandom3.h>
#include <TMath.h>

// Own includes
#include "ConfigurationCard.h"
#include "TriggerHistograms.h"
#include "ForestReader.h"

class TriggerAnalyzer{
  
public:
  
  // Constructors and destructor
  TriggerAnalyzer(); // Default constructor
  TriggerAnalyzer(std::vector<TString> fileNameVector, ConfigurationCard *newCard); // Custom constructor
  TriggerAnalyzer(const TriggerAnalyzer& in); // Copy constructor
  virtual ~TriggerAnalyzer(); // Destructor
  TriggerAnalyzer& operator=(const TriggerAnalyzer& obj); // Equal sign operator
  
  // Methods
  void RunAnalysis();                     // Run the dijet analysis
  TriggerHistograms* GetHistograms() const;   // Getter for histograms
  
private:
  
  // Private methods
  void ReadConfigurationFromCard(); // Read all the configuration from the input card
  
  Bool_t PassEventCuts(ForestReader *eventReader); // Check if the event passes the event cuts
  Double_t GetVzWeight(const Double_t vz) const;  // Get the proper vz weighting depending on analyzed system
  Double_t GetCentralityWeight(const Int_t hiBin) const; // Get the proper centrality weighting depending on analyzed system
  Double_t GetJetPtWeight(const Double_t jetPt) const; // Get the proper jet pT weighting for 2017 and 2018 MC
  
  // Private data members
  ForestReader *fJetReader;                 // Reader for jets in the event
  std::vector<TString> fFileNames;          // Vector for all the files to loop over
  ConfigurationCard *fCard;                 // Configuration card for the analysis
  TriggerHistograms *fHistograms;           // Filled histograms
  TF1 *fVzWeightFunction;                   // Weighting function for vz. Needed for MC.
  TF1 *fCentralityWeightFunctionCentral;    // Weighting function for central centrality classes. Needed for MC.
  TF1 *fCentralityWeightFunctionPeripheral; // Weighting function for peripheral centrality classes. Needed for MC.
  TF1 *fPtWeightFunction;                   // Weighting function for jet pT. Needed for MC.
  
  // Analyzed data and forest types
  Int_t fDataType;                   // Analyzed data type
  Int_t fJetType;                    // Type of jets used for analysis. 0 = Calo jets, 1 = PF jets
  Int_t fBaseTrigger;                // Trigger index used as base trigger for efficiency study
  Int_t fDebugLevel;                 // Amount of debug messages printed to console
  
  // Weights for filling the MC histograms
  Double_t fVzWeight;                // Weight for vz in MC
  Double_t fCentralityWeight;        // Weight for centrality in MC
  Double_t fPtHatWeight;             // Weight for pT hat in MC
  Double_t fTotalEventWeight;        // Combined weight factor for MC
  
  // Jet and track selection cuts
  Int_t fJetAxis;                      // Used jet axis type. 0 = Anti-kT jet axis, 1 = Axis from leading PF candidate
  Double_t fVzCut;                     // Cut for vertez z-position in an event
  Double_t fMinimumPtHat;              // Minimum accepted pT hat value
  Double_t fMaximumPtHat;              // Maximum accepted pT hat value
  Double_t fJetEtaCut;                 // Eta cut around midrapidity
  Double_t fJetMinimumPtCut;           // Minimum pT cut for jets
  Double_t fJetMaximumPtCut;           // Maximum pT accepted for jets (and tracks)
  Bool_t fCutBadPhiRegion;             // Cut the phi region with bad tracker performance from the analysis
  Double_t fMinimumMaxTrackPtFraction; // Cut for jets consisting only from soft particles
  Double_t fMaximumMaxTrackPtFraction; // Cut for jets consisting only from one high pT

};

#endif
