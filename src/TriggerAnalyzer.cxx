// Class for the main analysis algorithms for the energy-energy correlator analysis

// Root includes
#include <TFile.h>
#include <TMath.h>

// Own includes
#include "TriggerAnalyzer.h"

using namespace std;

/*
 * Default constructor
 */
TriggerAnalyzer::TriggerAnalyzer() :
  fFileNames(0),
  fCard(0),
  fHistograms(0),
  fVzWeightFunction(0),
  fCentralityWeightFunctionCentral(0),
  fCentralityWeightFunctionPeripheral(0),
  fPtWeightFunction(0),
  fDataType(-1),
  fJetType(0),
  fBaseTrigger(1),
  fDebugLevel(0),
  fVzWeight(1),
  fCentralityWeight(1),
  fPtHatWeight(1),
  fTotalEventWeight(1),
  fJetAxis(0),
  fVzCut(0),
  fMinimumPtHat(0),
  fMaximumPtHat(0),
  fJetEtaCut(0),
  fJetMinimumPtCut(0),
  fJetMaximumPtCut(0),
  fCutBadPhiRegion(false),
  fMinimumMaxTrackPtFraction(0),
  fMaximumMaxTrackPtFraction(0)
{
  // Default constructor
  fHistograms = new TriggerHistograms();
  fHistograms->CreateHistograms();
  
  // Initialize readers to null
  fJetReader = NULL;
  
}

/*
 * Custom constructor
 */
TriggerAnalyzer::TriggerAnalyzer(std::vector<TString> fileNameVector, ConfigurationCard *newCard) :
  fFileNames(fileNameVector),
  fCard(newCard),
  fHistograms(0),
  fVzWeight(1),
  fCentralityWeight(1),
  fPtHatWeight(1),
  fTotalEventWeight(1)
{
  // Custom constructor
  fHistograms = new TriggerHistograms(fCard);
  fHistograms->CreateHistograms();
  
  // Initialize readers to null
  fJetReader = NULL;
  
  // Configurure the analyzer from input card
  ReadConfigurationFromCard();
  
  // pT weight function for Pythia to match 2017 MC and data pT spectra. Derived from all jets above 120 GeV
  fPtWeightFunction = new TF1("fPtWeightFunction","pol3",0,500);
  fPtWeightFunction->SetParameters(0.79572,0.0021861,-6.35407e-06,6.66435e-09); // From JECv6
  
  // Find the correct folder for track correction tables based on data type
  if(fDataType == ForestReader::kPp || fDataType == ForestReader::kPpMC){
    
    // Weight function for 2017 MC
    fVzWeightFunction = new TF1("fvz","pol6",-15,15);  // Weight function for 2017 MC
    fVzWeightFunction->SetParameters(0.973805, 0.00339418, 0.000757544, -1.37331e-06, -2.82953e-07, -3.06778e-10, 3.48615e-09);
    fCentralityWeightFunctionCentral = NULL;
    fCentralityWeightFunctionPeripheral = NULL;
    
  } else if (fDataType == ForestReader::kPbPb || fDataType == ForestReader::kPbPbMC){
    
    // The vz weight function is rederived from the miniAOD dataset.
    // Macro used for derivation: deriveMonteCarloWeights.C, Git hash: f95771aa3242a7a9ed385c1ff364500481088eec
    // Input files: eecAnalysis_akFlowJets_wtaAxis_cutBadPhi_miniAODtesting_processed_2023-01-30.root
    //              PbPbMC2018_RecoGen_eecAnalysis_akFlowJets_mAOD_4pC_wtaAxis_jetTrig_cutBadPhi_processed_2023-02-10.root
    fVzWeightFunction = new TF1("fvz","pol6",-15,15);
    fVzWeightFunction->SetParameters(1.0082, -0.0190011, 0.000779051, -2.15118e-05, -6.70894e-06, 1.47181e-07, 6.65274e-09);
    
    // The centrality weight function is rederived for the miniAOD dataset.
    // Macro used for derivation: deriveMonteCarloWeights.C, Git hash: f95771aa3242a7a9ed385c1ff364500481088eec
    // Input files: eecAnalysis_akFlowJets_wtaAxis_cutBadPhi_miniAODtesting_processed_2023-01-30.root
    //              PbPbMC2018_RecoGen_eecAnalysis_akFlowJets_mAOD_4pC_wtaAxis_jetTrig_cutBadPhi_processed_2023-02-10.root
    fCentralityWeightFunctionCentral = new TF1("fCentralWeight","pol6",0,30);
    fCentralityWeightFunctionCentral->SetParameters(4.44918,-0.0544424, -0.0248668,0.00254486,-0.000117819,2.65985e-06,-2.35606e-08);
    fCentralityWeightFunctionPeripheral = new TF1("fPeripheralWeight","pol6",30,90);
    fCentralityWeightFunctionPeripheral->SetParameters(3.41938,-0.0643178, -0.00186948,7.67356e-05,-1.06981e-06,7.04102e-09,-1.84554e-11);
    
    
  } else {
    fVzWeightFunction = NULL;
    fCentralityWeightFunctionCentral = NULL;
    fCentralityWeightFunctionPeripheral = NULL;
  }
  
}

/*
 * Copy constructor
 */
TriggerAnalyzer::TriggerAnalyzer(const TriggerAnalyzer& in) :
  fJetReader(in.fJetReader),
  fFileNames(in.fFileNames),
  fCard(in.fCard),
  fHistograms(in.fHistograms),
  fVzWeightFunction(in.fVzWeightFunction),
  fCentralityWeightFunctionCentral(in.fCentralityWeightFunctionCentral),
  fCentralityWeightFunctionPeripheral(in.fCentralityWeightFunctionPeripheral),
  fPtWeightFunction(in.fPtWeightFunction),
  fDataType(in.fDataType),
  fJetType(in.fJetType),
  fBaseTrigger(in.fBaseTrigger),
  fDebugLevel(in.fDebugLevel),
  fVzWeight(in.fVzWeight),
  fCentralityWeight(in.fCentralityWeight),
  fPtHatWeight(in.fPtHatWeight),
  fTotalEventWeight(in.fTotalEventWeight),
  fJetAxis(in.fJetAxis),
  fVzCut(in.fVzCut),
  fMinimumPtHat(in.fMinimumPtHat),
  fMaximumPtHat(in.fMaximumPtHat),
  fJetEtaCut(in.fJetEtaCut),
  fJetMinimumPtCut(in.fJetMinimumPtCut),
  fJetMaximumPtCut(in.fJetMaximumPtCut),
  fCutBadPhiRegion(in.fCutBadPhiRegion),
  fMinimumMaxTrackPtFraction(in.fMinimumMaxTrackPtFraction),
  fMaximumMaxTrackPtFraction(in.fMaximumMaxTrackPtFraction)
{
  // Copy constructor
  
}

/*
 * Assingment operator
 */
TriggerAnalyzer& TriggerAnalyzer::operator=(const TriggerAnalyzer& in){
  // Assingment operator
  
  if (&in==this) return *this;
  
  fJetReader = in.fJetReader;
  fFileNames = in.fFileNames;
  fCard = in.fCard;
  fHistograms = in.fHistograms;
  fVzWeightFunction = in.fVzWeightFunction;
  fCentralityWeightFunctionCentral = in.fCentralityWeightFunctionCentral;
  fCentralityWeightFunctionPeripheral = in.fCentralityWeightFunctionPeripheral;
  fPtWeightFunction = in.fPtWeightFunction;
  fDataType = in.fDataType;
  fJetType = in.fJetType;
  fBaseTrigger = in.fBaseTrigger;
  fDebugLevel = in.fDebugLevel;
  fVzWeight = in.fVzWeight;
  fCentralityWeight = in.fCentralityWeight;
  fPtHatWeight = in.fPtHatWeight;
  fTotalEventWeight = in.fTotalEventWeight;
  fJetAxis = in.fJetAxis;
  fVzCut = in.fVzCut;
  fMinimumPtHat = in.fMinimumPtHat;
  fMaximumPtHat = in.fMaximumPtHat;
  fJetEtaCut = in.fJetEtaCut;
  fJetMinimumPtCut = in.fJetMinimumPtCut;
  fJetMaximumPtCut = in.fJetMaximumPtCut;
  fCutBadPhiRegion = in.fCutBadPhiRegion;
  fMinimumMaxTrackPtFraction = in.fMinimumMaxTrackPtFraction;
  fMaximumMaxTrackPtFraction = in.fMaximumMaxTrackPtFraction;
  
  return *this;
}

/*
 * Destructor
 */
TriggerAnalyzer::~TriggerAnalyzer(){
  // destructor
  delete fHistograms;
  if(fVzWeightFunction) delete fVzWeightFunction;
  if(fCentralityWeightFunctionCentral) delete fCentralityWeightFunctionCentral;
  if(fCentralityWeightFunctionPeripheral) delete fCentralityWeightFunctionPeripheral;
  if(fPtWeightFunction) delete fPtWeightFunction;
  if(fJetReader) delete fJetReader;
}

/*
 * Read all the configuration from the input card
 */
void TriggerAnalyzer::ReadConfigurationFromCard(){
  
  //****************************************
  //     Analyzed data type and trigger
  //****************************************
  fDataType = fCard->Get("DataType");
  
  //****************************************
  //         Event selection cuts
  //****************************************
  
  fVzCut = fCard->Get("ZVertexCut");          // Event cut vor the z-position of the primary vertex
  fMinimumPtHat = fCard->Get("LowPtHatCut");  // Minimum accepted pT hat value
  fMaximumPtHat = fCard->Get("HighPtHatCut"); // Maximum accepted pT hat value
  fBaseTrigger = fCard->Get("BaseTrigger");   // Base trigger for trigger efficiency analysis
  
  //****************************************
  //          Jet selection cuts
  //****************************************
  
  fJetEtaCut = fCard->Get("JetEtaCut");           // Eta cut around midrapidity
  fJetMinimumPtCut = fCard->Get("MinJetPtCut");   // Minimum pT cut for jets
  fJetMaximumPtCut = fCard->Get("MaxJetPtCut");   // Maximum pT accepted for jets (and tracks)
  fMinimumMaxTrackPtFraction = fCard->Get("MinMaxTrackPtFraction");  // Cut for jets consisting only from soft particles
  fMaximumMaxTrackPtFraction = fCard->Get("MaxMaxTrackPtFraction");  // Cut for jets consisting only from one high pT particle
  fCutBadPhiRegion = (fCard->Get("CutBadPhi") == 1);   // Flag for cutting the phi region with bad tracking efficiency from the analysis
  

  //****************************************
  //            Jet selection
  //****************************************
  fJetType = fCard->Get("JetType");              // Select the type of analyzed jets (Calo, CSPF, PuPF, FlowPF)
  fJetAxis = fCard->Get("JetAxis");              // Select between escheme and WTA axes

  
  //************************************************
  //              Debug messages
  //************************************************
  fDebugLevel = fCard->Get("DebugLevel");
}

/*
 * Main analysis loop
 */
void TriggerAnalyzer::RunAnalysis(){
  
  //************************************************
  //  Define variables needed in the analysis loop
  //************************************************
  
  // Input files and forest readers for analysis
  TFile *inputFile;
  
  // Event variables
  Int_t nEvents = 0;                // Number of events
  Double_t vz = 0;                  // Vertex z-position
  Double_t centrality = 0;          // Event centrality
  Int_t hiBin = 0;                  // CMS hiBin (centrality * 2)
  Double_t ptHat = 0;               // pT hat for MC events
  
  // Variables for jets
  Int_t nJets = 0;                  // Number of jets in an event
  Double_t jetPt = 0;               // pT of the i:th jet in the event
  Double_t jetPhi = 0;              // phi of the i:th jet in the event
  Double_t jetEta = 0;              // eta of the i:th jet in the event
  Double_t jetWTAPhi = 0;           // WTA phi of the i:th jet in the event
  Double_t jetWTAEta = 0;           // WTA eta of the i:th jet in the event
  Double_t jetPtWeight = 1;         // Weighting for jet pT
  Double_t leadingJetPt = 0;        // Leading jet pT
  Double_t leadingJetEta = 0;       // Leading jet eta
  Double_t leadingJetPhi = 0;       // Leading jet phi
  Double_t deltaRaxes = 0;          // DeltaR between WTA and Escheme axes
  
  // File name helper variables
  TString currentFile;
  
  // Trigger selection array for the events
  Bool_t triggerSelection[TriggerHistograms::knTriggerTypes];   // Table of triggers that fires this event
  Double_t triggerPrescale[TriggerHistograms::knTriggerTypes];  // Prescale value for the chosen trigger
  
  // Fillers for THnSparses
  const Int_t nFillJet = 6;
  Double_t fillerJet[nFillJet];

  
  //************************************************
  //      Define forest reader for data files
  //************************************************
  
  fJetReader = new ForestReader(fDataType, fJetType, fJetAxis, fBaseTrigger);
  
  
  //************************************************
  //       Main analysis loop over all files
  //************************************************
  
  // Loop over files
  Int_t nFiles = fFileNames.size();
  for(Int_t iFile = 0; iFile < nFiles; iFile++) {
    
    //************************************************
    //              Find and open files
    //************************************************
    
    // Find the filename and open the input file
    currentFile = fFileNames.at(iFile);
    inputFile = TFile::Open(currentFile);
    
    // Check that the file exists
    if(!inputFile){
      cout << "Error! Could not find the file: " << currentFile.Data() << endl;
      assert(0);
    }

    // Check that the file is open
    if(!inputFile->IsOpen()){
      cout << "Error! Could not open the file: " << currentFile.Data() << endl;
      assert(0);
    }
    
    // Check that the file is not zombie
    if(inputFile->IsZombie()){
      cout << "Error! The following file is a zombie: " << currentFile.Data() << endl;
      assert(0);
    }
    

    // Print the used files
    if(fDebugLevel > 0) cout << "Reading from file: " << currentFile.Data() << endl;

    
    //************************************************
    //            Read forest from file
    //************************************************
    
    // If file is good, read the forest from the file
    fJetReader->ReadForestFromFile(inputFile);  // There might be a memory leak in handling the forest...
    nEvents = fJetReader->GetNEvents();

    //************************************************
    //         Main event loop for each file
    //************************************************
    
    for(Int_t iEvent = 0; iEvent < nEvents; iEvent++){ // nEvents
      
      //************************************************
      //         Read basic event information
      //************************************************
      
      // Print to console how the analysis is progressing
      if(fDebugLevel > 1 && iEvent % 1000 == 0) cout << "Analyzing event " << iEvent << endl;
      
      // Read the event to memory
      fJetReader->GetEvent(iEvent);

      // Get vz, centrality and pT hat information
      vz = fJetReader->GetVz();
      centrality = fJetReader->GetCentrality();
      hiBin = fJetReader->GetHiBin();
      ptHat = fJetReader->GetPtHat();
      
      // We need to apply pT hat cuts before getting pT hat weight. There might be rare events above the upper
      // limit from which the weights are calculated, which could cause the code to crash.
      if(ptHat < fMinimumPtHat || ptHat >= fMaximumPtHat) continue;
      
      // Get the weighting for the event
      fVzWeight = GetVzWeight(vz);
      fCentralityWeight = GetCentralityWeight(hiBin);
      fPtHatWeight = fJetReader->GetEventWeight();
      fTotalEventWeight = fVzWeight*fCentralityWeight*fPtHatWeight;
      
      // Fill event counter histogram
      fHistograms->fhEvents->Fill(TriggerHistograms::kAll);          // All the events looped over
      
      //  ============================================
      //  ===== Apply all the event quality cuts =====
      //  ============================================
      
      if(!PassEventCuts(fJetReader)) continue;
      
      // Fill the event information histograms for the events that pass the event cuts
      fHistograms->fhVertexZ->Fill(vz);                            // z vertex distribution from all events
      fHistograms->fhVertexZWeighted->Fill(vz,fVzWeight);          // z-vertex distribution weighted with the weight function
      fHistograms->fhCentrality->Fill(centrality);                 // Centrality filled from all events
      fHistograms->fhCentralityWeighted->Fill(centrality,fCentralityWeight); // Centrality weighted with the centrality weighting function
      fHistograms->fhPtHat->Fill(ptHat);                           // pT hat histogram
      fHistograms->fhPtHatWeighted->Fill(ptHat,fPtHatWeight);      // pT het histogram weighted with corresponding cross section and event number
      
      // Determine the trigger selection for all included triggers
      for(Int_t iTrigger = 0; iTrigger < TriggerHistograms::knTriggerTypes; iTrigger++){
        triggerSelection[iTrigger] = fJetReader->GetJetFilterBit(iTrigger);
        triggerPrescale[iTrigger] = fJetReader->GetJetTriggerPrescale(iTrigger);
      }
      
      // Set the prescale for the base trigger branch to one, as it is meaningless after the selection
      triggerPrescale[fBaseTrigger] = 1;
      
      // ======================================
      // ===== Event quality cuts applied =====
      // ======================================
      
      //***********************************************************************
      //    Loop over all jets and fill histograms for different triggers
      //***********************************************************************
      
      // Jet loop
      nJets = fJetReader->GetNJets();
      leadingJetPt = 0; leadingJetEta = 0; leadingJetPhi = 0;
      for(Int_t jetIndex = 0; jetIndex < nJets; jetIndex++) {
        
        jetPt = fJetReader->GetJetPt(jetIndex);
        jetPhi = fJetReader->GetJetPhi(jetIndex);
        jetEta = fJetReader->GetJetEta(jetIndex);
        
        //  ========================================
        //  ======== Apply jet quality cuts ========
        //  ========================================
        
        if(TMath::Abs(jetEta) >= fJetEtaCut) continue; // Cut for jet eta
        if(fCutBadPhiRegion && (jetPhi > -0.1 && jetPhi < 1.2)) continue; // Cut the area of large inefficiency in tracker
        
        if(fMinimumMaxTrackPtFraction >= fJetReader->GetJetMaxTrackPt(jetIndex)/fJetReader->GetJetRawPt(jetIndex)) {
          continue; // Cut for jets with only very low pT particles
        }
        if(fMaximumMaxTrackPtFraction <= fJetReader->GetJetMaxTrackPt(jetIndex)/fJetReader->GetJetRawPt(jetIndex)) {
          continue; // Cut for jets where all the pT is taken by one track
        }
        
        
        //  ========================================
        //  ======= Jet quality cuts applied =======
        //  ========================================
        
        // After the jet pT can been corrected, apply analysis jet pT cuts
        if(jetPt < fJetMinimumPtCut) continue;
        if(jetPt > fJetMaximumPtCut) continue;
        
        //************************************************
        //         Fill histograms for all jets
        //************************************************
        
        // Remember the leading jet
        if(jetPt > leadingJetPt) {
          leadingJetPt = jetPt; leadingJetEta = jetEta; leadingJetPhi = jetPhi;
        }
        
        // Find the pT weight for the jet
        jetPtWeight = GetJetPtWeight(jetPt);
        
        // Fill the axes in correct order
        fillerJet[0] = jetPt;          // Axis 0 = jet pT
        fillerJet[1] = jetPhi;         // Axis 1 = jet phi
        fillerJet[2] = jetEta;         // Axis 2 = jet eta
        fillerJet[3] = centrality;     // Axis 3 = centrality
        fillerJet[4] = TriggerHistograms::kReconstructed;  // Axis 4 = Reconstruction flag
        fillerJet[5] = TriggerHistograms::knTriggerTypes;  // Axis 5 = Trigger selection
        
        fHistograms->fhInclusiveJet->Fill(fillerJet,fTotalEventWeight*jetPtWeight); // Fill the data point to histogram
        
        // Apply a trigger selection on top of the base trigger to evaluate the trigger efficiency of the other triggers
        for(Int_t iTrigger = 0; iTrigger < TriggerHistograms::knTriggerTypes; iTrigger++){
          if(triggerSelection[iTrigger]){
            fillerJet[5] = iTrigger;
            
            fHistograms->fhInclusiveJet->Fill(fillerJet,fTotalEventWeight*jetPtWeight*triggerPrescale[iTrigger]); // Fill the data point to histogram
          }
        }
        
      } // End of jet loop
      
      // =============================== //
      // Fill the leading jet histograms //
      // =============================== //
      
      // Find the pT weight for the jet
      jetPtWeight = GetJetPtWeight(leadingJetPt);
      
      // Fill the axes in correct order
      fillerJet[0] = leadingJetPt;          // Axis 0 = leading jet pT
      fillerJet[1] = leadingJetPhi;         // Axis 1 = leading jet phi
      fillerJet[2] = leadingJetEta;         // Axis 2 = leading jet eta
      fillerJet[3] = centrality;            // Axis 3 = centrality
      fillerJet[4] = TriggerHistograms::kReconstructed;  // Axis 4 = Reconstruction flag
      fillerJet[5] = TriggerHistograms::knTriggerTypes;  // Axis 5 = Trigger selection
      
      fHistograms->fhLeadingJet->Fill(fillerJet,fTotalEventWeight*jetPtWeight); // Fill the data point to histogram
      
      // Apply a trigger selection on top of the base trigger to evaluate the trigger efficiency of the other triggers
      for(Int_t iTrigger = 0; iTrigger < TriggerHistograms::knTriggerTypes; iTrigger++){
        if(triggerSelection[iTrigger]){
          fillerJet[5] = iTrigger;
          
          fHistograms->fhLeadingJet->Fill(fillerJet,fTotalEventWeight*jetPtWeight*triggerPrescale[iTrigger]); // Fill the data point to histogram
        }
      }
      
      // For MC, do another jet loop using generator level jets
      if(fDataType == ForestReader::kPpMC || fDataType == ForestReader::kPbPbMC){
        
        // Generator level jet loop
        nJets = fJetReader->GetNGeneratorJets();
        leadingJetPt = 0; leadingJetEta = 0; leadingJetPhi = 0;
        for(Int_t jetIndex = 0; jetIndex < nJets; jetIndex++) {
          
          jetPt = fJetReader->GetGeneratorJetPt(jetIndex);
          jetPhi = fJetReader->GetGeneratorJetPhi(jetIndex);
          jetEta = fJetReader->GetGeneratorJetEta(jetIndex);
          jetWTAPhi = fJetReader->GetGeneratorJetWTAPhi(jetIndex);
          jetWTAEta = fJetReader->GetGeneratorJetWTAEta(jetIndex);
          
          //  ==========================================
          //  ======== Apply jet kinematic cuts ========
          //  ==========================================
          
          if(TMath::Abs(jetEta) >= fJetEtaCut) continue; // Cut for jet eta
          if(jetPt < fJetMinimumPtCut) continue;
          if(jetPt > fJetMaximumPtCut) continue;
          
          //************************************************
          //     Fill histograms for generator level jets
          //************************************************
          
          // Remember the leading jet
          if(jetPt > leadingJetPt) {
            leadingJetPt = jetPt; leadingJetEta = jetEta; leadingJetPhi = jetPhi;
          }
          
          // Find the pT weight for the jet
          jetPtWeight = GetJetPtWeight(jetPt);
          
          // Fill the axes in correct order
          fillerJet[0] = jetPt;          // Axis 0 = generator level jet pT
          fillerJet[1] = jetPhi;         // Axis 1 = generator level jet phi
          fillerJet[2] = jetEta;         // Axis 2 = generator level jet eta
          fillerJet[3] = centrality;     // Axis 3 = centrality
          fillerJet[4] = TriggerHistograms::kGeneratorLevel;   // Axis 4 = Generator level flag
          fillerJet[5] = TriggerHistograms::knTriggerTypes;    // Axis 5 = Trigger selection
          
          fHistograms->fhInclusiveJet->Fill(fillerJet,fTotalEventWeight*jetPtWeight); // Fill the data point to histogram
          
          // Fill the deltaR between WTA and E-schame axes
          deltaRaxes = GetDeltaR(jetEta, jetPhi, jetWTAEta, jetWTAPhi);
          fHistograms->fhGenJetDeltaR->Fill(deltaRaxes, fTotalEventWeight*jetPtWeight);
          
          // Apply a trigger selection on top of the base trigger to evaluate the trigger efficiency of the other triggers
          for(Int_t iTrigger = 0; iTrigger < TriggerHistograms::knTriggerTypes; iTrigger++){
            if(triggerSelection[iTrigger]){
              fillerJet[5] = iTrigger;
              
              fHistograms->fhInclusiveJet->Fill(fillerJet,fTotalEventWeight*jetPtWeight*triggerPrescale[iTrigger]); // Fill the data point to histogram
            }
          }
          
        } // End of jet loop
        
        // =============================================== //
        // Fill the leading generator level jet histograms //
        // =============================================== //
        
        // Find the pT weight for the jet
        jetPtWeight = GetJetPtWeight(leadingJetPt);
        
        // Fill the axes in correct order
        fillerJet[0] = leadingJetPt;          // Axis 0 = leading generator level jet pT
        fillerJet[1] = leadingJetPhi;         // Axis 1 = leading generator level jet phi
        fillerJet[2] = leadingJetEta;         // Axis 2 = leading generator level jet eta
        fillerJet[3] = centrality;            // Axis 3 = centrality
        fillerJet[4] = TriggerHistograms::kGeneratorLevel; // Axis 4 = Generator level flag
        fillerJet[5] = TriggerHistograms::knTriggerTypes;  // Axis 5 = Trigger selection
        
        fHistograms->fhLeadingJet->Fill(fillerJet,fTotalEventWeight*jetPtWeight); // Fill the data point to histogram
        
        // Apply a trigger selection on top of the base trigger to evaluate the trigger efficiency of the other triggers
        for(Int_t iTrigger = 0; iTrigger < TriggerHistograms::knTriggerTypes; iTrigger++){
          if(triggerSelection[iTrigger]){
            fillerJet[5] = iTrigger;
            
            fHistograms->fhLeadingJet->Fill(fillerJet,fTotalEventWeight*jetPtWeight*triggerPrescale[iTrigger]); // Fill the data point to histogram
          }
        }
        
      } // MC if
      
      
    } // Event loop
    
    //************************************************
    //      Cleanup at the end of the file loop
    //************************************************
    
    // Close the input files after the event has been read
    inputFile->Close();
    
  } // File loop
  
}


/*
 * Get the proper vz weighting depending on analyzed system
 *
 *  Arguments:
 *   const Double_t vz = Vertex z position for the event
 *
 *   return: Multiplicative correction factor for vz
 */
Double_t TriggerAnalyzer::GetVzWeight(const Double_t vz) const{
  if(fDataType == ForestReader::kPp || fDataType == ForestReader::kPbPb) return 1;  // No correction for real data
  if(fDataType == ForestReader::kPbPbMC || fDataType == ForestReader::kPpMC) return fVzWeightFunction->Eval(vz); // Weight for 2018 MC
  return -1; // Return crazy value for unknown data types, so user will not miss it
}

/*
 * Get the proper centrality weighting depending on analyzed system
 *
 *  Arguments:
 *   const Int_t hiBin = CMS hiBin
 *
 *   return: Multiplicative correction factor for the given CMS hiBin
 */
Double_t TriggerAnalyzer::GetCentralityWeight(const Int_t hiBin) const{
  if(fDataType != ForestReader::kPbPbMC) return 1;
  
  // No weighting for the most peripheral centrality bins. Different weight function for central and peripheral.
  if(hiBin < 60) return fCentralityWeightFunctionCentral->Eval(hiBin/2.0);
  return (hiBin < 194) ? fCentralityWeightFunctionPeripheral->Eval(hiBin/2.0) : 1;
}

/*
 * Get the proper jet pT weighting depending on analyzed system
 *
 *  Arguments:
 *   const Double_t jetPt = Jet pT for the weighted jet
 *
 *   return: Multiplicative correction factor for the jet pT
 */
Double_t TriggerAnalyzer::GetJetPtWeight(const Double_t jetPt) const{
  if(fDataType == ForestReader::kPbPb || fDataType == ForestReader::kPp) return 1.0;  // No weight for data
  
  return fPtWeightFunction->Eval(jetPt);
}

/*
 * Check if the event passes all the track cuts
 *
 *  Arguments:
 *   ForestReader *eventReader = ForestReader containing the event information checked for event cuts
 *
 *   return = True if all event cuts are passed, false otherwise
 */
Bool_t TriggerAnalyzer::PassEventCuts(ForestReader *eventReader){

  // Primary vertex has at least two tracks, is within 25 cm in z-rirection and within 2 cm in xy-direction. Only applied for data.
  if(eventReader->GetPrimaryVertexFilterBit() == 0) return false;
  fHistograms->fhEvents->Fill(TriggerHistograms::kPrimaryVertex);
  
  // Have at least two HF towers on each side of the detector with an energy deposit of 4 GeV. Only applied for PbPb data.
  if(eventReader->GetHfCoincidenceFilterBit() == 0) return false;
  fHistograms->fhEvents->Fill(TriggerHistograms::kHfCoincidence);
  
  // Calculated from pixel clusters. Ensures that measured and predicted primary vertices are compatible. Only applied for PbPb data.
  if(eventReader->GetClusterCompatibilityFilterBit() == 0) return false;
  fHistograms->fhEvents->Fill(TriggerHistograms::kClusterCompatibility);
  
  // Cut for beam scraping. Only applied for pp data.
  if(eventReader->GetBeamScrapingFilterBit() == 0) return false;
  fHistograms->fhEvents->Fill(TriggerHistograms::kBeamScraping);
  
  // Jet trigger requirement.
  if(eventReader->GetBaseJetFilterBit() == 0) return false;
  fHistograms->fhEvents->Fill(TriggerHistograms::kCaloJet);
  
  // Cut for vertex z-position
  if(TMath::Abs(eventReader->GetVz()) > fVzCut) return false;
  fHistograms->fhEvents->Fill(TriggerHistograms::kVzCut);
  
  return true;
  
}

/*
 * Getter for trigger histograms
 */
TriggerHistograms* TriggerAnalyzer::GetHistograms() const{
  return fHistograms;
}

/*
 * Get deltaR between two objects
 *
 *  Arguments:
 *   const Double_t eta1 = Eta of the first object
 *   const Double_t phi1 = Phi of the first object
 *   const Double_t eta2 = Eta of the second object
 *   const Double_t phi2 = Phi of the second object
 *
 *  return: DeltaR between the two objects
 */
Double_t TriggerAnalyzer::GetDeltaR(const Double_t eta1, const Double_t phi1, const Double_t eta2, const Double_t phi2) const{

  Double_t deltaEta = eta1 - eta2;
  Double_t deltaPhi = phi1 - phi2;
  
  // Transform deltaPhi to interval [-pi,pi]
  while(deltaPhi > TMath::Pi()){deltaPhi += -2*TMath::Pi();}
  while(deltaPhi < -TMath::Pi()){deltaPhi += 2*TMath::Pi();}
  
  // Return the distance between the objects
  return TMath::Sqrt(deltaPhi*deltaPhi + deltaEta*deltaEta);
  
}
