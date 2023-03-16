// Histograms needed in the trigger analysis

// C++ includes
#include <assert.h>

// Root includes
#include <TFile.h>
#include <TMath.h>

// Own includes
#include "TriggerHistograms.h"

/*
 * Default constructor
 */
TriggerHistograms::TriggerHistograms() :
  fhVertexZ(0),
  fhVertexZWeighted(0),
  fhEvents(0),
  fhCentrality(0),
  fhCentralityWeighted(0),
  fhPtHat(0),
  fhPtHatWeighted(0),
  fhGenJetDeltaR(0),
  fhInclusiveJet(0),
  fhLeadingJet(0),
  fCard(0)
{
  // Default constructor
  
}

/*
 * Custom constructor
 */
TriggerHistograms::TriggerHistograms(ConfigurationCard *newCard) :
  fhVertexZ(0),
  fhVertexZWeighted(0),
  fhEvents(0),
  fhCentrality(0),
  fhCentralityWeighted(0),
  fhPtHat(0),
  fhPtHatWeighted(0),
  fhGenJetDeltaR(0),
  fhInclusiveJet(0),
  fhLeadingJet(0),
  fCard(newCard)
{
  // Custom constructor

}

/*
 * Copy constructor
 */
TriggerHistograms::TriggerHistograms(const TriggerHistograms& in) :
  fhVertexZ(in.fhVertexZ),
  fhVertexZWeighted(in.fhVertexZWeighted),
  fhEvents(in.fhEvents),
  fhCentrality(in.fhCentrality),
  fhCentralityWeighted(in.fhCentralityWeighted),
  fhPtHat(in.fhPtHat),
  fhPtHatWeighted(in.fhPtHatWeighted),
  fhGenJetDeltaR(in.fhGenJetDeltaR),
  fhInclusiveJet(in.fhInclusiveJet),
  fhLeadingJet(in.fhLeadingJet),
  fCard(in.fCard)
{
  // Copy constructor
  
}

/*
 * Assingment operator
 */
TriggerHistograms& TriggerHistograms::operator=(const TriggerHistograms& in){
  // Assingment operator
  
  if (&in==this) return *this;
  
  fhVertexZ = in.fhVertexZ;
  fhVertexZWeighted = in.fhVertexZWeighted;
  fhEvents = in.fhEvents;
  fhCentrality = in.fhCentrality;
  fhCentralityWeighted = in.fhCentralityWeighted;
  fhPtHat = in.fhPtHat;
  fhPtHatWeighted = in.fhPtHatWeighted;
  fhGenJetDeltaR = in.fhGenJetDeltaR;
  fhInclusiveJet = in.fhInclusiveJet;
  fhLeadingJet = in.fhLeadingJet;
  fCard = in.fCard;
  
  return *this;
}

/*
 * Destructor
 */
TriggerHistograms::~TriggerHistograms(){
  // destructor
  delete fhVertexZ;
  delete fhVertexZWeighted;
  delete fhEvents;
  delete fhCentrality;
  delete fhCentralityWeighted;
  delete fhPtHat;
  delete fhPtHatWeighted;
  delete fhGenJetDeltaR;
  delete fhInclusiveJet;
  delete fhLeadingJet;
}

/*
 * Set the configuration card used for the histogram class
 */
void TriggerHistograms::SetCard(ConfigurationCard *newCard){
  fCard = newCard;
}

/*
 * Getter for the trigger name for given index
 */
TString TriggerHistograms::GetTriggerName(Int_t iTrigger) const{
  if(iTrigger < 0 || iTrigger >= knTriggerTypes) return "IndexOutOfBounds";
  return kTriggerStrings[iTrigger];
}

/*
 * Create the necessary histograms
 */
void TriggerHistograms::CreateHistograms(){
  
  // ======== Common binning information for histograms =========
  
  // Centrality
  const Double_t minCentrality = -0.75;   // Minimum centrality bin, is negative since hiBin is -1 for pp
  const Double_t maxCentrality = 100.25;  // Maximum centrality bin
  const Int_t nCentralityBins = 202;      // Number of centrality bins
  
  // Jet pT
  const Double_t minPtJet = 0;     // Minimum jet pT
  const Double_t maxPtJet = 500;   // Maximum jet pT
  const Int_t nPtBinsJet = 100;    // Number of jet pT bins
  
  // Phi
  const Double_t minPhi = -TMath::Pi();  // Minimum phi
  const Double_t maxPhi = TMath::Pi();   // Maximum phi
  const Int_t nPhiBins = 64;             // Number of phi bins
  
  // Eta
  const Double_t minEta = -2.5;    // Minimum eta (current eta cut for tracks = 2.4)
  const Double_t maxEta = 2.5;     // Maximum eta (current eta cut for tracks = 2.4)
  const Int_t nEtaBins = 50;       // Number of eta bins
  
  // Vertex z-position
  const Double_t minVz = -20;   // Minimum vz
  const Double_t maxVz = 20;    // Maximum vz
  const Int_t nVzBins = 80;     // Number of vz bins
  
  // pT hat
  const Double_t minPtHat = 0;     // Minimum pT hat
  const Double_t maxPtHat = 460;   // Maximum pT hat
  const Int_t nFinePtHatBins = 230; // Number of fine pT hat bins
  
  // Reconstructed/generator level
  const Double_t minDataLevel = -0.5;
  const Double_t maxDataLevel = knDataLevels-0.5;
  const Int_t nDataLevelBins = knDataLevels;
  
  // Trigger selection for the jet histogram
  const Double_t minTriggerSelection = -0.5;
  const Double_t maxTriggerSelection = knTriggerTypes+0.5;
  const Int_t nTriggerSelectionBins = knTriggerTypes+1;
  
  // DeltaR for E-scheme axis vs. WTA axis
  const Double_t minDeltaR = 0;
  const Double_t maxDeltaR = 4;
  const Int_t nDeltaRbins = 80;
  
  // Centrality bins for THnSparses (We run into memory issues, if have all the bins)
  const Int_t nWideCentralityBins = fCard->GetNBin("CentralityBinEdges");
  Double_t wideCentralityBins[nWideCentralityBins+1];
  for(Int_t iCentrality = 0; iCentrality < nWideCentralityBins+1; iCentrality++){
    wideCentralityBins[iCentrality] = fCard->Get("CentralityBinEdges",iCentrality);
  }
  
  // Bins for the pT hat histogram
  const Int_t nPtHatBins = fCard->GetNBin("PtHatBinEdges");
  Double_t ptHatBins[nPtHatBins+1];
  for(Int_t iPtHat = 0; iPtHat < nPtHatBins+1; iPtHat++){
    ptHatBins[iPtHat] = fCard->Get("PtHatBinEdges",iPtHat);
  }
  
  const Int_t nAxesJet = 6;
  Int_t nBinsJet[nAxesJet];
  Double_t lowBinBorderJet[nAxesJet];
  Double_t highBinBorderJet[nAxesJet];
  
  // ======== Plain TH1 histograms ========
  
  fhVertexZ = new TH1F("vertexZ","vertexZ",nVzBins,minVz,maxVz); fhVertexZ->Sumw2();
  fhVertexZWeighted = new TH1F("vertexZweighted","vertexZweighted",nVzBins,minVz,maxVz); fhVertexZWeighted->Sumw2();
  fhEvents = new TH1F("nEvents","nEvents",knEventTypes,-0.5,knEventTypes-0.5); fhEvents->Sumw2();
  fhCentrality = new TH1F("centrality","centrality",nCentralityBins,minCentrality,maxCentrality); fhCentrality->Sumw2();
  fhCentralityWeighted = new TH1F("centralityWeighted","centralityWeighted",nCentralityBins,minCentrality,maxCentrality); fhCentralityWeighted->Sumw2();
  fhPtHat = new TH1F("pthat","pthat",nPtHatBins,ptHatBins); fhPtHat->Sumw2();
  fhPtHatWeighted = new TH1F("pthatWeighted","pthatWeighted",nFinePtHatBins,minPtHat,maxPtHat); fhPtHatWeighted->Sumw2();
  fhGenJetDeltaR = new TH1F("genJetDeltaR","genJetDeltaR",nDeltaRbins,minDeltaR,maxDeltaR); fhGenJetDeltaR->Sumw2();
  
  // For the event histogram, label each bin corresponding to an event cut
  for(Int_t i = 0; i < knEventTypes; i++){
    fhEvents->GetXaxis()->SetBinLabel(i+1,kEventTypeStrings[i]);
  }
  
  // If we are using PF jets, change the axis label for that
  if(fCard->Get("JetType") == 1) fhEvents->GetXaxis()->SetBinLabel(kCaloJet+1,"PFJet");
  
  // ======== THnSparse for all jets ========
  
  // Axis 0 for the jet histogram: jet pT
  nBinsJet[0] = nPtBinsJet;         // nBins for any jet pT
  lowBinBorderJet[0] = minPtJet;    // low bin border for any jet pT
  highBinBorderJet[0] = maxPtJet;   // high bin border for any jet pT
  
  // Axis 1 for the jet histogram: jet phi
  nBinsJet[1] = nPhiBins;        // nBins for any jet phi
  lowBinBorderJet[1] = minPhi;   // low bin border for any jet phi
  highBinBorderJet[1] = maxPhi;  // high bin border for any jet phi
  
  // Axis 2 for the jet histogram: jet eta
  nBinsJet[2] = nEtaBins;        // nBins for any jet eta
  lowBinBorderJet[2] = minEta;   // low bin border for any jet eta
  highBinBorderJet[2] = maxEta;  // high bin border for any jet eta
  
  // Axis 3 for the jet histogram: centrality
  nBinsJet[3] = nWideCentralityBins;   // nBins for wide centrality bins
  lowBinBorderJet[3] = minCentrality;  // low bin border for centrality
  highBinBorderJet[3] = maxCentrality; // high bin border for centrality
  
  // Axis 4 for the jet histogram: data level (reconstructed / generator level)
  nBinsJet[4] = nDataLevelBins;       // nBins different data levels
  lowBinBorderJet[4] = minDataLevel;  // low bin border for data levels
  highBinBorderJet[4] = maxDataLevel; // high bin border for data levels
  
  // Axis 5 for the jet histogram: trigger selection
  nBinsJet[5] = nTriggerSelectionBins;       // nBins different data levels
  lowBinBorderJet[5] = minTriggerSelection;  // low bin border for data levels
  highBinBorderJet[5] = maxTriggerSelection; // high bin border for data levels
  
  // Create the histogram for all jets using the above binning information
  fhInclusiveJet = new THnSparseF("inclusiveJet","inclusiveJet",nAxesJet,nBinsJet,lowBinBorderJet,highBinBorderJet); fhInclusiveJet->Sumw2();
  fhLeadingJet = new THnSparseF("leadingJet","leadingJet",nAxesJet,nBinsJet,lowBinBorderJet,highBinBorderJet); fhLeadingJet->Sumw2();

  // Set custom centrality bins for histograms
  fhInclusiveJet->SetBinEdges(3,wideCentralityBins);
  fhLeadingJet->SetBinEdges(3,wideCentralityBins);
}

/*
 * Write the histograms to file
 */
void TriggerHistograms::Write() const{
  
  // Write the histograms to file
  fhVertexZ->Write();
  fhVertexZWeighted->Write();
  fhEvents->Write();
  fhCentrality->Write();
  fhCentralityWeighted->Write();
  fhPtHat->Write();
  fhPtHatWeighted->Write();
  fhGenJetDeltaR->Write();
  fhInclusiveJet->Write();
  fhLeadingJet->Write();
  
}

/*
 * Write the histograms to a given file
 */
void TriggerHistograms::Write(TString outputFileName) const{
  
  // Define the output file
  TFile *outputFile = new TFile(outputFileName, "RECREATE");
  
  // Write the histograms to file
  Write();
  
  // Close the file after everything is written
  outputFile->Close();
  
  // Delete the outputFile object
  delete outputFile;
}


