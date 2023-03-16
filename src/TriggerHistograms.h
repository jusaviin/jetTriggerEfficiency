// Class for histograms needed in the trigger analysis

#ifndef TRIGGERHISTOGRAMS_H
#define TRIGGERHISTOGRAMS_H

// Root includes
#include <TH1.h>
#include <TH2.h>
#include <THnSparse.h>

// Own includes
#include "ConfigurationCard.h"

class TriggerHistograms{
  
public:
  
  // Enumeration for event types to event histogram and track cuts for track cut histogram
  enum enumEventTypes {kAll, kPrimaryVertex, kHfCoincidence, kClusterCompatibility, kBeamScraping, kCaloJet, kVzCut, knEventTypes};
  enum enumTriggerSelection {kCalo40, kCalo60, kCalo80, kCalo100, kPF60, kPF80, kPF100, knTriggerTypes};
  enum enumDataLevel {kReconstructed, kGeneratorLevel, knDataLevels};
  
  // Constructors and destructor
  TriggerHistograms(); // Default constructor
  TriggerHistograms(ConfigurationCard *newCard); // Custom constructor
  TriggerHistograms(const TriggerHistograms& in); // Copy constructor
  virtual ~TriggerHistograms(); // Destructor
  TriggerHistograms& operator=(const TriggerHistograms& obj); // Equal sign operator
  
  // Methods
  void CreateHistograms();                      // Create all histograms
  void Write() const;                           // Write the histograms to a file that is opened somewhere else
  void Write(TString outputFileName) const;     // Write the histograms to a file
  void SetCard(ConfigurationCard *newCard);     // Set a new configuration card for the histogram class
  TString GetTriggerName(Int_t iTrigger) const; // Getter for the trigger name
  
  // Histograms defined public to allow easier access to them. Should not be abused
  // Notation in comments: l = leading jet, s = subleading jet, inc - inclusive jet, uc = uncorrected, ptw = pT weighted
  TH1F *fhVertexZ;                 // Vertex z-position
  TH1F *fhVertexZWeighted;         // Weighted vertex z-position (only meaningfull for MC)
  TH1F *fhEvents;                  // Number of events. For binning see enumEventTypes.
  TH1F *fhCentrality;              // Centrality information. -0.5 for pp or PYTHIA.
  TH1F *fhCentralityWeighted;      // Weighted centrality distribution (only meaningful for MC)
  TH1F *fhPtHat;                   // pT hat for MC events (only meaningful for MC)
  TH1F *fhPtHatWeighted;           // Weighted pT hat distribution
  TH1F *fhGenJetDeltaR;            // Generator level jet deltaR between e-scheme and WTA
  THnSparseF *fhInclusiveJet;      // Inclusive jet information. Axes: [jet pT][jet phi][jet eta][cent][reco/gen][trigger]
  THnSparseF *fhLeadingJet;        // Leading jet information. Axes: [jet pT][jet phi][jet eta][cent][reco/gen][trigger]
  
private:
  
  ConfigurationCard *fCard;    // Card for binning info
  const TString kEventTypeStrings[knEventTypes] = {"All", "PrimVertex", "HfCoin2Th4", "ClustCompt", "BeamScrape", "CaloJet", "v_{z} cut"}; // Strings corresponding to event types
  const TString kTriggerStrings[knTriggerTypes] = {"CaloJet40", "CaloJet60", "CaloJet80", "CaloJet100", "PFJet60", "PFJet80", "PFJet100"};
  
};

#endif
