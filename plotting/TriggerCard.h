#ifndef TRIGGERCARD_H
#define TRIGGERCARD_H

// C++ includes
#include <iostream>

// Root includes
#include <TFile.h>
#include <TDirectory.h>
#include <TString.h>
#include <TObjString.h>
#include <TVectorT.h>

/*
 * TriggerCard class
 *
 * This class reads the ConfigurationCard from the input root file and decodes the
 * necessary information for the analysis.
 */
class TriggerCard {
  
public:
 
  // Indices for card entries
  enum enumCardEntries{
    kDataType,                  // Data type in the data file (pp, PbPb, pp MC, PbPb MC)
    kBaseTrigger,               // Base jet trigger with respect to which the other jet triggers are evaluated
    kJetType,                   // 0 = Calorimeter jets, 1 = PF jets
    kJetAxis,                   // 0 = Anti-kt axis, 1 = Leading PF candidate axis, 2 = WTA axis
    kJetEtaCut,                 // Eta cut for jets
    kMinPtCut,                  // Minimum allowed pT for the inclusive jets
    kMaxPtCut,                  // Maximum allowed pT for the inclusive jets
    kCutBadPhiRegion,           // Cut the phi region with bad tracked performance from the analysis
    kMinMaxTrackPtFraction,     // Minimum fraction of jet pT taken by the highest pT track in jet
    kMaxMaxTrackPtFraction,     // Maximum fraction of jet pT taken by the highest pT track in jet
    kZVertexCut,                // Maximum accepted vz in the event
    kLowPtHatCut,               // Minimum accepted pT hat
    kHighPtHatCut,              // Maximum accepted pT hat
    kCentralityBinEdges,        // Centrality bin edges
    kPtHatBinEdges,             // pT hat bin edges
    knEntries};                 // Number of entries in the card
  
  // Enumeration for input files used in postprocessing
  enum enumFileNames{kInputFileName,knFileNames};
  
private:
  
  // Names for each entry read from the configuration card
  const char *fCardEntryNames[knEntries] = {"DataType","BaseTrigger","JetType","JetAxis","JetEtaCut","MinJetPtCut","MaxJetPtCut","CutBadPhi","MinMaxTrackPtFraction","MaxMaxTrackPtFraction","ZVertexCut","LowPtHatCut","HighPtHatCut","CentralityBinEdges","PtHatBinEdges"};
  const char *fFileNameType[knFileNames] = {"input"};
  const char *fFileNameSaveName[knFileNames] = {"InputFile"};
  
  TFile *fInputFile;         // Input file from which all the data is read
  TString fCardDirectory;    // Path to the ConfigurationCard directory
  int fDataType;             // Total number of centrality bins in the analysis
  TString fDataTypeString;   // Total number of eta gaps in the analysis
  TString fAlternativeDataTypeString; // Alternative data type string
  
  void FindDataTypeString(); // Construct a data type string based on information on the card
  void ReadVectors();        // Read the vectors from the file
  
  // Strings for git hash
  TObjString *fGitHash;
  TObjString *fProjectionGitHash;
  
  // Vectors for all the lines inside the card
  TVectorT<float> *fCardEntries[knEntries];   // Array of all the vectors in the card
  TObjString *fFileNames[knFileNames];        // Array for filenames used in postprocessing
  
  // Private methods
  int GetNBins(const int index) const;                            // Get the number of bins for internal index
  double GetLowBinBorder(const int index, const int iBin) const;  // Get the low border of i:th bin from internal index
  double GetHighBinBorder(const int index, const int iBin) const; // Get the high border of i:th bin from internal index
  int GetBinIndex(const int index, const double value) const;     // Get the bin index in the i:th bin from internal index based on given value
   
public:
  
  TriggerCard(TFile *inFile); // Contructor with input file
  ~TriggerCard();             // Destructor
  
  TString GetDataType() const;             // Getter for data type string
  TString GetAlternativeDataType() const;  // Getter for alternative data type string
  void Write(TDirectory *file);            // Write the contents of the card to a file
  void Print() const;                      // Print the contents of the card to the console
  
  int GetNCentralityBins() const; // Get the number of centrality bins
  double GetLowBinBorderCentrality(const int iBin) const;  // Get the low border of i:th centrality bin
  double GetHighBinBorderCentrality(const int iBin) const; // Get the high border of i:th centrality bin
  int GetBinIndexCentrality(const double value) const;     // Get the bin index for a given centrality value
  int GetJetType() const;          // Get the jet type index
  double GetJetPtCut() const;      // Get the minimum jet pT cut
  int GetBaseTrigger() const;      // Get the index of the base trigger
  
  void AddOneDimensionalVector(int entryIndex, float entryContent); // Add one dimensional vector to the card
  void AddVector(int entryIndex, int dimension, double *contents); // Add a vector to the card
  void AddFileName(int entryIndex, TString fileName); // Add a file name to the card
  void AddProjectionGitHash(const char* gitHash); // Add a git hash used to project the histograms to the file
  
};

#endif
