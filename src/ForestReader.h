// Reader for jet trees from CMS data
//
//===========================================================
// ForestReader.h
//
// Author: Jussi Viinikainen
//===========================================================

#ifndef FORESTREADER_H
#define FORESTREADER_H

// C++ includes
#include <iostream>
#include <assert.h>
#include <vector>

// Root includes
#include <TString.h>
#include <TTree.h>
#include <TChain.h>
#include <TBranch.h>
#include <TFile.h>

// Own includes
#include "TriggerHistograms.h"

using namespace std;

class ForestReader{
  
private:
  static const Int_t fnMaxJet = 250;        // Maximum number of jets in an event
  
public:
  
  // Possible data types to be read with the reader class
  enum enumDataTypes{kPp, kPbPb, kPpMC, kPbPbMC, knDataTypes};
  
  // Constructors and destructors
  ForestReader();                                          // Default constructor
  ForestReader(Int_t dataType, Int_t jetType, Int_t jetAxis, Int_t baseTrigger); // Custom constructor
  ForestReader(const ForestReader& in);                    // Copy constructor
  ~ForestReader();                                 // Destructor
  ForestReader& operator=(const ForestReader& obj);        // Equal sign operator
  
  // Methods
  void GetEvent(Int_t nEvent);                 // Get the nth event in tree
  Int_t GetNEvents() const;                        // Get the number of events
  void ReadForestFromFile(TFile *inputFile);   // Read the forest from a file
  void ReadForestFromFileList(std::vector<TString> fileList);   // Read the forest from a file list
  void BurnForest();                           // Burn the forest
  
  // Getters for leaves in heavy ion tree
  Float_t GetVz() const;              // Getter for vertex z position
  Float_t GetCentrality() const;      // Getter for centrality
  Int_t GetHiBin() const;             // Getter for CMS hiBin
  Float_t GetPtHat() const;           // Getter for pT hat
  Float_t GetEventWeight() const;     // Getter for event weight in MC
  
  // Getters for leaves in jet tree
  Int_t GetNJets() const;                     // Getter for number of jets
  Float_t GetJetPt(Int_t iJet) const;         // Getter for jet pT
  Float_t GetJetPhi(Int_t iJet) const;        // Getter for jet phi
  Float_t GetJetEta(Int_t iJet) const;        // Getter for jet eta
  Float_t GetJetRawPt(Int_t iJet) const;      // Getter for jet raw pT
  Float_t GetJetMaxTrackPt(Int_t iJet) const; // Getter for maximum track pT inside a jet
  
  Int_t GetNGeneratorJets() const;                   // Getter for number of generator level jets
  Float_t GetGeneratorJetPt(Int_t iJet) const;       // Getter for generator level jet pT
  Float_t GetGeneratorJetPhi(Int_t iJet) const;      // Getter for generator level jet phi
  Float_t GetGeneratorJetEta(Int_t iJet) const;      // Getter for generator level jet eta
  Float_t GetGeneratorJetWTAPhi(Int_t iJet) const;   // Getter for generator level jet WTA phi
  Float_t GetGeneratorJetWTAEta(Int_t iJet) const;   // Getter for generator level jet WTA eta
  
  // Getters for leaves in HLT tree
  Int_t GetBaseJetFilterBit() const;                    // Getter for base jet filter bit
  Int_t GetJetFilterBit(Int_t iTrigger) const;          // Getter for the selected jet filter bit
  Double_t GetJetTriggerPrescale(Int_t iTrigger) const; // Getter for the prescale value of the chosen trigger
  
  // Getters for leaves in skim tree
  Int_t GetPrimaryVertexFilterBit() const;           // Getter for primary vertex filter bit
  Int_t GetBeamScrapingFilterBit() const;            // Getter got beam scraping filter bit
  Int_t GetHfCoincidenceFilterBit() const;           // Getter for hadronic forward coincidence filter bit
  Int_t GetClusterCompatibilityFilterBit() const;    // Getter for cluster compatibility filter bit
  
  // Setter for data type
  void SetDataType(Int_t dataType); // Setter for data type
  
private:
  
  // Methods
  void Initialize();      // Connect the branches to the tree
    
  Int_t fDataType;        // Type of data read with the tree. 0 = pp, 1 = PbPb, 2 = ppMC, 3 = PbPbMC, 4 = LocalTest
  Int_t fJetType;         // Choose the type of jets usedfor analysis. 0 = Calo jets, 1 = PF jets
  Int_t fJetAxis;         // Jet axis used for the jets. 0 = Anti-kT, 1 = Leading particle flow candidate, 2 = WTA
  Int_t fBaseTrigger;     // Trigger index that is used as a base trigger with respect to which other triggers are compared
  Bool_t fIsMiniAOD;      // Flag for type of the forest True = MiniAOD forest, False = AOD forest
  
  // Trees in the forest
  TTree *fHeavyIonTree;    // Tree for heavy ion event information
  TTree *fJetTree;         // Tree for jet information
  TTree *fHltTree;         // Tree for high level trigger information
  TTree *fSkimTree;        // Tree for event selection information
  
  // Branches for heavy ion tree
  TBranch *fHiVzBranch;                   // Branch for vertex z-position
  TBranch *fHiBinBranch;                  // Branch for centrality
  TBranch *fPtHatBranch;                  // Branch for pT hat
  TBranch *fEventWeightBranch;            // Branch for event weight
  
  // Branches for jet tree
  TBranch *fnJetsBranch;         // Branch for number of jets
  TBranch *fJetPtBranch;         // Branch for jet pT
  TBranch *fJetPhiBranch;        // Branch for jet phi
  TBranch *fJetEtaBranch;        // Branch for jet eta
  TBranch *fJetRawPtBranch;      // Branch for raw jet pT
  TBranch *fJetMaxTrackPtBranch; // Maximum pT for a track inside a jet
  
  TBranch *fnGenJetsBranch;      // Branch for number of generator level jets
  TBranch *fGenJetPtBranch;      // Branch for generator level jet pT
  TBranch *fGenJetPhiBranch;     // Branch for generator level jet phi
  TBranch *fGenJetEtaBranch;     // Branch for generator level jet eta
  TBranch *fGenJetWTAPhiBranch;     // Branch for generator level jet phi
  TBranch *fGenJetWTAEtaBranch;     // Branch for generator level jet eta
  
  // Branches for HLT tree
  TBranch *fJetFilterBranch[TriggerHistograms::knTriggerTypes];  // Branches for all jet trigger bits that are studied
  TBranch *fJetFilterPrescaleNumeratorBranch[TriggerHistograms::knTriggerTypes];    // Branches for jet trigger bit prescale numerators
  TBranch *fJetFilterPrescaleDenominatorBranch[TriggerHistograms::knTriggerTypes];  // Branches for jet trigger bit prescale denominators
  
  // Branches for skim tree
  TBranch *fPrimaryVertexBranch;           // Branch for primary vertex filter bit
  TBranch *fBeamScrapingBranch;            // Branch for beam scraping filter bit
  TBranch *fHfCoincidenceBranch;           // Branch for energy recorded in at least 3 HF calorimeter towers
  TBranch *fClusterCompatibilityBranch;    // Branch for cluster compatibility
    
  // Leaves for heavy ion tree
  Float_t fVertexZ;    // Vertex z-position
  Int_t fHiBin;        // HiBin = Centrality percentile * 2
  Float_t fPtHat;      // pT hat
  
  // Leaves for jet tree
  Int_t fnJets;          // number of jets in an event
  Int_t fnGenJets;       // Number of generator level jets in an event
  Float_t fEventWeight;  // jet weight in the MC tree
  
  Float_t fJetPtArray[fnMaxJet] = {0};         // pT:s of all the jets in an event
  Float_t fJetPhiArray[fnMaxJet] = {0};        // phis of all the jets in an event
  Float_t fJetEtaArray[fnMaxJet] = {0};        // etas of all the jets in an event
  Float_t fJetRawPtArray[fnMaxJet] = {0};      // raw jet pT for all the jets in an event
  Float_t fJetMaxTrackPtArray[fnMaxJet] = {0}; // maximum track pT inside a jet for all the jets in an event
  
  Float_t fGenJetPtArray[fnMaxJet] = {0};      // pT:s of the generator level jets in an event
  Float_t fGenJetPhiArray[fnMaxJet] = {0};     // phis of the generator level jets in an event
  Float_t fGenJetEtaArray[fnMaxJet] = {0};     // etas of the generator level jets in an event
  Float_t fGenJetWTAPhiArray[fnMaxJet] = {0};  // WTA phis of the generator level jets in an event
  Float_t fGenJetWTAEtaArray[fnMaxJet] = {0};  // WTA etas of the generator level jets in an event
  
  // Leaves for the HLT tree
  Int_t fJetFilterBit[TriggerHistograms::knTriggerTypes];  // Filter bits for all jet triggers
  Int_t fJetPrescaleNumerator[TriggerHistograms::knTriggerTypes];   // Prescale numerators for jet triggers
  Int_t fJetPrescaleDenominator[TriggerHistograms::knTriggerTypes]; // Prescale denominators for jet triggers
  
  // Leaves for the skim tree
  Int_t fPrimaryVertexFilterBit;           // Filter bit for primary vertex
  Int_t fBeamScrapingFilterBit;            // Filter bit for beam scraping
  Int_t fHfCoincidenceFilterBit;           // Filter bit for energy recorded in at least 3 HF calorimeter towers
  Int_t fClusterCompatibilityFilterBit;    // Filter bit for cluster compatibility
  
  
};

#endif
