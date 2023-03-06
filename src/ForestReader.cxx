// Implementation for ForestReader

// Own includes
#include "ForestReader.h"

/*
 * Default constructor
 */
ForestReader::ForestReader() :
  fDataType(0),
  fJetType(0),
  fJetAxis(0),
  fBaseTrigger(0),
  fIsMiniAOD(false),
  fHeavyIonTree(0),
  fJetTree(0),
  fHltTree(0),
  fSkimTree(0),
  fHiVzBranch(0),
  fHiBinBranch(0),
  fPtHatBranch(0),
  fEventWeightBranch(0),
  fnJetsBranch(0),
  fJetPtBranch(0),
  fJetPhiBranch(0),
  fJetEtaBranch(0),
  fJetRawPtBranch(0),
  fJetMaxTrackPtBranch(0),
  fnGenJetsBranch(0),
  fGenJetPtBranch(0),
  fGenJetPhiBranch(0),
  fGenJetEtaBranch(0),
  fPrimaryVertexBranch(0),
  fBeamScrapingBranch(0),
  fHfCoincidenceBranch(0),
  fClusterCompatibilityBranch(0),
  fVertexZ(-100),
  fHiBin(-1),
  fPtHat(0),
  fnJets(0),
  fnGenJets(0),
  fEventWeight(1),
  fJetPtArray(),
  fJetPhiArray(),
  fJetEtaArray(),
  fJetRawPtArray(),
  fPrimaryVertexFilterBit(0),
  fBeamScrapingFilterBit(0),
  fHfCoincidenceFilterBit(0),
  fClusterCompatibilityFilterBit(0)

{
  // Default constructor
  
  // Initialize fJetMaxTrackPtArray to -1
  for(Int_t i = 0; i < fnMaxJet; i++){
    fJetMaxTrackPtArray[i] = -1;
  }
  
  // Initialize jet filter branches to NULL and bits to 1
  for(Int_t iTrigger = 0; iTrigger < TriggerHistograms::knTriggerTypes; iTrigger++){
    fJetFilterBranch[iTrigger] = NULL;
    fJetFilterPrescaleNumeratorBranch[iTrigger] = NULL;
    fJetFilterPrescaleDenominatorBranch[iTrigger] = NULL;
    fJetFilterBit[iTrigger] = 1;
    fJetPrescaleNumerator[iTrigger] = 1;
    fJetPrescaleDenominator[iTrigger] = 1;
  }
  
}

/*
 * Custom constructor
 *
 *  Arguments:
 *   Int_t dataType: 0 = pp, 1 = PbPb, 2 = pp MC, 3 = PbPb MC
 *   Int_t jetType: 0 = Calo jets, 1 = CSPF jets, 2 = PuPF jets, 3 = Flow jets
 *   Int_t jetAxis: 0 = Anti-kT axis, 1 = WTA axis
 *   Int_t baseTrigger: 0 = CaloJet40, 1 = CaloJet60, 2 = CaloJet80, 3 = CaloJet100, 4 = PFJet60, 5 = PFJet80, 6 = PFJet100
 */
ForestReader::ForestReader(Int_t dataType, Int_t jetType, Int_t jetAxis, Int_t baseTrigger) :
  fDataType(0),
  fJetType(jetType),
  fJetAxis(jetAxis),
  fBaseTrigger(baseTrigger),
  fIsMiniAOD(false),
  fHeavyIonTree(0),
  fJetTree(0),
  fHltTree(0),
  fSkimTree(0),
  fHiVzBranch(0),
  fHiBinBranch(0),
  fPtHatBranch(0),
  fEventWeightBranch(0),
  fnJetsBranch(0),
  fJetPtBranch(0),
  fJetPhiBranch(0),
  fJetEtaBranch(0),
  fJetRawPtBranch(0),
  fJetMaxTrackPtBranch(0),
  fnGenJetsBranch(0),
  fGenJetPtBranch(0),
  fGenJetPhiBranch(0),
  fGenJetEtaBranch(0),
  fPrimaryVertexBranch(0),
  fBeamScrapingBranch(0),
  fHfCoincidenceBranch(0),
  fClusterCompatibilityBranch(0),
  fVertexZ(-100),
  fHiBin(-1),
  fPtHat(0),
  fnJets(0),
  fnGenJets(0),
  fEventWeight(1),
  fJetPtArray(),
  fJetPhiArray(),
  fJetEtaArray(),
  fJetRawPtArray(),
  fPrimaryVertexFilterBit(0),
  fBeamScrapingFilterBit(0),
  fHfCoincidenceFilterBit(0),
  fClusterCompatibilityFilterBit(0)
{
  // Custom constructor
  
  SetDataType(dataType);
  
  // Initialize fJetMaxTrackPtArray to -1
  for(int i = 0; i < fnMaxJet; i++){
    fJetMaxTrackPtArray[i] = -1;
  }
 
  // Initialize jet filter branches to NULL and bits to 1
  for(Int_t iTrigger = 0; iTrigger < TriggerHistograms::knTriggerTypes; iTrigger++){
    fJetFilterBranch[iTrigger] = NULL;
    fJetFilterPrescaleNumeratorBranch[iTrigger] = NULL;
    fJetFilterPrescaleDenominatorBranch[iTrigger] = NULL;
    fJetFilterBit[iTrigger] = 1;
    fJetPrescaleNumerator[iTrigger] = 1;
    fJetPrescaleDenominator[iTrigger] = 1;
  }
}

/*
 * Copy constructor
 */
ForestReader::ForestReader(const ForestReader& in) :
  fDataType(in.fDataType),
  fJetType(in.fJetType),
  fJetAxis(in.fJetAxis),
  fBaseTrigger(in.fBaseTrigger),
  fIsMiniAOD(in.fIsMiniAOD),
  fHeavyIonTree(in.fHeavyIonTree),
  fJetTree(in.fJetTree),
  fHltTree(in.fHltTree),
  fSkimTree(in.fSkimTree),
  fHiVzBranch(in.fHiVzBranch),
  fHiBinBranch(in.fHiBinBranch),
  fPtHatBranch(in.fPtHatBranch),
  fEventWeightBranch(in.fEventWeightBranch),
  fnJetsBranch(in.fnJetsBranch),
  fJetPtBranch(in.fJetPtBranch),
  fJetPhiBranch(in.fJetPhiBranch),
  fJetEtaBranch(in.fJetEtaBranch),
  fJetRawPtBranch(in.fJetRawPtBranch),
  fJetMaxTrackPtBranch(in.fJetMaxTrackPtBranch),
  fnGenJetsBranch(in.fnGenJetsBranch),
  fGenJetPtBranch(in.fGenJetPtBranch),
  fGenJetPhiBranch(in.fGenJetPhiBranch),
  fGenJetEtaBranch(in.fGenJetEtaBranch),
  fPrimaryVertexBranch(in.fPrimaryVertexBranch),
  fBeamScrapingBranch(in.fBeamScrapingBranch),
  fHfCoincidenceBranch(in.fHfCoincidenceBranch),
  fClusterCompatibilityBranch(in.fClusterCompatibilityBranch),
  fVertexZ(in.fVertexZ),
  fHiBin(in.fHiBin),
  fPtHat(in.fPtHat),
  fnJets(in.fnJets),
  fnGenJets(in.fnGenJets),
  fEventWeight(in.fEventWeight),
  fPrimaryVertexFilterBit(in.fPrimaryVertexFilterBit),
  fBeamScrapingFilterBit(in.fBeamScrapingFilterBit),
  fHfCoincidenceFilterBit(in.fHfCoincidenceFilterBit),
  fClusterCompatibilityFilterBit(in.fClusterCompatibilityFilterBit)
{
  // Copy constructor
  
  for(Int_t i = 0; i < fnMaxJet; i++){
    fJetPtArray[i] = in.fJetPtArray[i];
    fJetPhiArray[i] = in.fJetPhiArray[i];
    fJetEtaArray[i] = in.fJetEtaArray[i];
    fJetRawPtArray[i] = in.fJetRawPtArray[i];
    fJetMaxTrackPtArray[i] = in.fJetMaxTrackPtArray[i];
  }
  
  // Copy trigger branches and bits
  for(Int_t iTrigger = 0; iTrigger < TriggerHistograms::knTriggerTypes; iTrigger++){
    fJetFilterBranch[iTrigger] = in.fJetFilterBranch[iTrigger];
    fJetFilterPrescaleNumeratorBranch[iTrigger] = in.fJetFilterPrescaleNumeratorBranch[iTrigger];
    fJetFilterPrescaleDenominatorBranch[iTrigger] = in.fJetFilterPrescaleDenominatorBranch[iTrigger];
    fJetFilterBit[iTrigger] = in.fJetFilterBit[iTrigger];
    fJetPrescaleNumerator[iTrigger] = in.fJetPrescaleNumerator[iTrigger];
    fJetPrescaleDenominator[iTrigger] = in.fJetPrescaleDenominator[iTrigger];
  }
}

/*
 * Assignment operator
 */
ForestReader& ForestReader::operator=(const ForestReader& in){
  // Assignment operator
  
  if (&in==this) return *this;
  
  fDataType = in.fDataType;
  fJetType = in.fJetType;
  fJetAxis = in.fJetAxis;
  fBaseTrigger = in.fBaseTrigger;
  fIsMiniAOD = in.fIsMiniAOD;
  fHeavyIonTree = in.fHeavyIonTree;
  fJetTree = in.fJetTree;
  fHltTree = in.fHltTree;
  fSkimTree = in.fSkimTree;
  fHiVzBranch = in.fHiVzBranch;
  fHiBinBranch = in.fHiBinBranch;
  fPtHatBranch = in.fPtHatBranch;
  fEventWeightBranch = in.fEventWeightBranch;
  fnJetsBranch = in.fnJetsBranch;
  fJetPtBranch = in.fJetPtBranch;
  fJetPhiBranch = in.fJetPhiBranch;
  fJetEtaBranch = in.fJetEtaBranch;
  fJetRawPtBranch = in.fJetRawPtBranch;
  fJetMaxTrackPtBranch = in.fJetMaxTrackPtBranch;
  fnGenJetsBranch = in.fnGenJetsBranch;
  fGenJetPtBranch = in.fGenJetPtBranch;
  fGenJetPhiBranch = in.fGenJetPhiBranch;
  fGenJetEtaBranch = in.fGenJetEtaBranch;
  fPrimaryVertexBranch = in.fPrimaryVertexBranch;
  fBeamScrapingBranch = in.fBeamScrapingBranch;
  fHfCoincidenceBranch = in.fHfCoincidenceBranch;
  fClusterCompatibilityBranch = in.fClusterCompatibilityBranch;
  fVertexZ = in.fVertexZ;
  fHiBin = in.fHiBin;
  fPtHat = in.fPtHat;
  fnJets = in.fnJets;
  fnGenJets = in.fnGenJets;
  fEventWeight = in.fEventWeight;
  fPrimaryVertexFilterBit = in.fPrimaryVertexFilterBit;
  fBeamScrapingFilterBit = in.fBeamScrapingFilterBit;
  fHfCoincidenceFilterBit = in.fHfCoincidenceFilterBit;
  fClusterCompatibilityFilterBit = in.fClusterCompatibilityFilterBit;
  
  for(Int_t i = 0; i < fnMaxJet; i++){
    fJetPtArray[i] = in.fJetPtArray[i];
    fJetPhiArray[i] = in.fJetPhiArray[i];
    fJetEtaArray[i] = in.fJetEtaArray[i];
    fJetRawPtArray[i] = in.fJetRawPtArray[i];
    fJetMaxTrackPtArray[i] = in.fJetMaxTrackPtArray[i];
  }
  
  // Copy trigger branches and bits
  for(Int_t iTrigger = 0; iTrigger < TriggerHistograms::knTriggerTypes; iTrigger++){
    fJetFilterBranch[iTrigger] = in.fJetFilterBranch[iTrigger];
    fJetFilterPrescaleNumeratorBranch[iTrigger] = in.fJetFilterPrescaleNumeratorBranch[iTrigger];
    fJetFilterPrescaleDenominatorBranch[iTrigger] = in.fJetFilterPrescaleDenominatorBranch[iTrigger];
    fJetFilterBit[iTrigger] = in.fJetFilterBit[iTrigger];
    fJetPrescaleNumerator[iTrigger] = in.fJetPrescaleNumerator[iTrigger];
    fJetPrescaleDenominator[iTrigger] = in.fJetPrescaleDenominator[iTrigger];
  }
  
  return *this;
}

/*
 * Destructor
 */
ForestReader::~ForestReader(){
  // destructor
}

/*
 * Initialization, meaning that the branches are connected to the tree
 */
void ForestReader::Initialize(){
  
  // Connect the branches of the heavy ion tree
  fHeavyIonTree->SetBranchStatus("*",0);
  fHeavyIonTree->SetBranchStatus("vz",1);
  fHeavyIonTree->SetBranchAddress("vz",&fVertexZ,&fHiVzBranch);
  fHeavyIonTree->SetBranchStatus("hiBin",1);
  fHeavyIonTree->SetBranchAddress("hiBin",&fHiBin,&fHiBinBranch);
  if(fDataType == kPpMC || fDataType == kPbPbMC){
    fHeavyIonTree->SetBranchStatus("pthat",1);
    fHeavyIonTree->SetBranchAddress("pthat",&fPtHat,&fPtHatBranch); // pT hat only for MC
    fHeavyIonTree->SetBranchStatus("weight",1);
    fHeavyIonTree->SetBranchAddress("weight",&fEventWeight,&fEventWeightBranch); // event weight only for MC
  } else {
    fPtHat = 0; // We do not have pT hat information for real data
    fEventWeight = 1;
  }
  
  // Connect the branches to the jet tree
  const char *jetAxis[2] = {"jt", "WTA"};
  const char *genJetAxis[2] = {"", "WTA"};
  char branchName[30];
  
  fJetTree->SetBranchStatus("*",0);
  fJetTree->SetBranchStatus("jtpt",1);
  fJetTree->SetBranchAddress("jtpt",&fJetPtArray,&fJetPtBranch);
  
  // If specified, select WTA axis for jet phi
  sprintf(branchName,"%sphi",jetAxis[fJetAxis]);
  fJetTree->SetBranchStatus(branchName,1);
  fJetTree->SetBranchAddress(branchName,&fJetPhiArray,&fJetPhiBranch);
  
  // If specified, select WTA axis for jet eta
  sprintf(branchName,"%seta",jetAxis[fJetAxis]);
  fJetTree->SetBranchStatus(branchName,1);
  fJetTree->SetBranchAddress(branchName,&fJetEtaArray,&fJetEtaBranch);
  
  fJetTree->SetBranchStatus("nref",1);
  fJetTree->SetBranchAddress("nref",&fnJets,&fnJetsBranch);
  fJetTree->SetBranchStatus("rawpt",1);
  fJetTree->SetBranchAddress("rawpt",&fJetRawPtArray,&fJetRawPtBranch);
  fJetTree->SetBranchStatus("trackMax",1);
  fJetTree->SetBranchAddress("trackMax",&fJetMaxTrackPtArray,&fJetMaxTrackPtBranch);
  
  // If we are looking at Monte Carlo, connect the reference pT and parton arrays
  if(fDataType > kPbPb){
    fJetTree->SetBranchStatus("genpt",1);
    fJetTree->SetBranchAddress("genpt",&fGenJetPtArray,&fGenJetPtBranch);
    
    // If specified, select WTA axis for jet phi
    sprintf(branchName,"%sgenphi",genJetAxis[fJetAxis]);
    fJetTree->SetBranchStatus(branchName,1);
    fJetTree->SetBranchAddress(branchName,&fGenJetPhiArray,&fGenJetPhiBranch);
    
    // If specified, select WTA axis for jet eta
    sprintf(branchName,"%sgeneta",genJetAxis[fJetAxis]);
    fJetTree->SetBranchStatus(branchName,1);
    fJetTree->SetBranchAddress(branchName,&fGenJetEtaArray,&fGenJetEtaBranch);
    
    fJetTree->SetBranchStatus("ngen",1);
    fJetTree->SetBranchAddress("ngen",&fnGenJets,&fnGenJetsBranch);
  }
  
  // Event selection summary
  //
  //         tree                      branch                         What it is
  //  hltanalysis/HltTree   HLT_HIPuAK4CaloJet100_Eta5p1_v1      Event selection for PbPb
  //  hltanalysis/HltTree      HLT_AK4CaloJet80_Eta5p1_v1         Event selection for pp
  // skimanalysis/HltTree         pprimaryVertexFilter           Event selection for PbPb
  // skimanalysis/HltTree    HBHENoiseFilterResultRun2Loose   Event selection for pp and PbPb
  // skimanalysis/HltTree         pPAprimaryVertexFilter          Event selection for pp
  // skimanalysis/HltTree           pBeamScrapingFilter           Event selection for pp
  
  // Connect the branches to the HLT tree
  fHltTree->SetBranchStatus("*",0);
  TriggerHistograms *triggerProvider = new TriggerHistograms();
  
  if(fDataType == kPp || fDataType == kPpMC){ // pp data or MC
    
    for(Int_t iTrigger = TriggerHistograms::kCalo40; iTrigger <= TriggerHistograms::kPF100; iTrigger++){
      fHltTree->SetBranchStatus(Form("HLT_HIAK4%s_v1", triggerProvider->GetTriggerName(iTrigger).Data()),1);
      fHltTree->SetBranchAddress(Form("HLT_HIAK4%s_v1", triggerProvider->GetTriggerName(iTrigger).Data()), &fJetFilterBit[iTrigger], &fJetFilterBranch[iTrigger]);
      fHltTree->SetBranchStatus(Form("HLT_HIAK4%s_v1_Prescl", triggerProvider->GetTriggerName(iTrigger).Data()),1);
      fHltTree->SetBranchAddress(Form("HLT_HIAK4%s_v1_Prescl", triggerProvider->GetTriggerName(iTrigger).Data()), &fJetPrescaleNumerator[iTrigger], &fJetFilterPrescaleNumeratorBranch[iTrigger]);
      
      
      // Only integer prescales for AOD
      fJetPrescaleDenominator[iTrigger] = 1;
    }

  } else { // PbPb data or MC
    
    for(Int_t iTrigger = TriggerHistograms::kCalo40; iTrigger <= TriggerHistograms::kCalo100; iTrigger++){
      fHltTree->SetBranchStatus(Form("HLT_HIPuAK4%sEta5p1_v1", triggerProvider->GetTriggerName(iTrigger).Data()),1);
      fHltTree->SetBranchAddress(Form("HLT_HIPuAK4%sEta5p1_v1", triggerProvider->GetTriggerName(iTrigger).Data()), &fJetFilterBit[iTrigger], &fJetFilterBranch[iTrigger]);
      fHltTree->SetBranchStatus(Form("HLT_HIPuAK4%sEta5p1_v1_PrescaleNumerator", triggerProvider->GetTriggerName(iTrigger).Data()),1);
      fHltTree->SetBranchAddress(Form("HLT_HIPuAK4%sEta5p1_v1_PrescaleNumerator", triggerProvider->GetTriggerName(iTrigger).Data()), &fJetPrescaleNumerator[iTrigger], &fJetFilterPrescaleNumeratorBranch[iTrigger]);
      fHltTree->SetBranchStatus(Form("HLT_HIPuAK4%sEta5p1_v1_PrescaleDenominator", triggerProvider->GetTriggerName(iTrigger).Data()),1);
      fHltTree->SetBranchAddress(Form("HLT_HIPuAK4%sEta5p1_v1_PrescaleDenominator", triggerProvider->GetTriggerName(iTrigger).Data()), &fJetPrescaleDenominator[iTrigger], &fJetFilterPrescaleDenominatorBranch[iTrigger]);
    }
    for(Int_t iTrigger = TriggerHistograms::kPF60; iTrigger <= TriggerHistograms::kPF100; iTrigger++){
      fHltTree->SetBranchStatus(Form("HLT_HICsAK4%sEta1p5_v1", triggerProvider->GetTriggerName(iTrigger).Data()),1);
      fHltTree->SetBranchAddress(Form("HLT_HICsAK4%sEta1p5_v1", triggerProvider->GetTriggerName(iTrigger).Data()), &fJetFilterBit[iTrigger], &fJetFilterBranch[iTrigger]);
      fHltTree->SetBranchStatus(Form("HLT_HICsAK4%sEta1p5_v1_PrescaleNumerator", triggerProvider->GetTriggerName(iTrigger).Data()),1);
      fHltTree->SetBranchAddress(Form("HLT_HICsAK4%sEta1p5_v1_PrescaleNumerator", triggerProvider->GetTriggerName(iTrigger).Data()), &fJetPrescaleNumerator[iTrigger], &fJetFilterPrescaleNumeratorBranch[iTrigger]);
      fHltTree->SetBranchStatus(Form("HLT_HICsAK4%sEta1p5_v1_PrescaleDenominator", triggerProvider->GetTriggerName(iTrigger).Data()),1);
      fHltTree->SetBranchAddress(Form("HLT_HICsAK4%sEta1p5_v1_PrescaleDenominator", triggerProvider->GetTriggerName(iTrigger).Data()), &fJetPrescaleDenominator[iTrigger], &fJetFilterPrescaleDenominatorBranch[iTrigger]);
    }
    
  }
  
  delete triggerProvider;
  
  // Connect the branches to the skim tree (different for pp and PbPb data and Monte Carlo)
  fSkimTree->SetBranchStatus("*",0);
  // pprimaryVertexFilter && phfCoincFilter2Th4 && pclusterCompatibilityFilter
  if(fDataType == kPp || fDataType == kPpMC){ // pp data or MC
    fSkimTree->SetBranchStatus("pPAprimaryVertexFilter",1);
    fSkimTree->SetBranchAddress("pPAprimaryVertexFilter",&fPrimaryVertexFilterBit,&fPrimaryVertexBranch);
    fSkimTree->SetBranchStatus("pBeamScrapingFilter",1);
    fSkimTree->SetBranchAddress("pBeamScrapingFilter",&fBeamScrapingFilterBit,&fBeamScrapingBranch);
    fHfCoincidenceFilterBit = 1; // No HF energy coincidence requirement for pp
    fClusterCompatibilityFilterBit = 1; // No cluster compatibility requirement for pp
  } else { // PbPb data or MC
    
    // Primary vertex has at least two tracks, is within 25 cm in z-rirection and within 2 cm in xy-direction
    fSkimTree->SetBranchStatus("pprimaryVertexFilter",1);
    fSkimTree->SetBranchAddress("pprimaryVertexFilter",&fPrimaryVertexFilterBit,&fPrimaryVertexBranch);
    
    // Cut on noise on HCAL
    if(fIsMiniAOD){
      // Have at least two HF towers on each side of the detector with an energy deposit of 4 GeV
      fSkimTree->SetBranchStatus("pphfCoincFilter2Th4",1);
      fSkimTree->SetBranchAddress("pphfCoincFilter2Th4", &fHfCoincidenceFilterBit, &fHfCoincidenceBranch);
      
    } else {
      
      // Have at least two HF towers on each side of the detector with an energy deposit of 4 GeV
      fSkimTree->SetBranchStatus("phfCoincFilter2Th4",1);
      fSkimTree->SetBranchAddress("phfCoincFilter2Th4", &fHfCoincidenceFilterBit, &fHfCoincidenceBranch);
    }
    
    // Calculated from pixel clusters. Ensures that measured and predicted primary vertices are compatible
    fSkimTree->SetBranchStatus("pclusterCompatibilityFilter",1);
    fSkimTree->SetBranchAddress("pclusterCompatibilityFilter",&fClusterCompatibilityFilterBit,&fClusterCompatibilityBranch);
    
    fBeamScrapingFilterBit = 1;  // No beam scraping filter for PbPb
  }
  
}


/*
 * Setter for fDataType
 */
void ForestReader::SetDataType(Int_t dataType){
  
  //Sanity check for given data type
  if(dataType < 0 || dataType > knDataTypes-1){
    cout << "ERROR: Data type input " << dataType << " is invalid in ForestReader.cxx!" << endl;
    cout << "Please give integer between 0 and " << knDataTypes-1 << "." << endl;
    cout << "Setting data type to 0 (pp)." << endl;
    fDataType = 0;
  } else {
    
    // If the sanity check passes, set the given data type
    fDataType = dataType;
  }
}

/*
 * Connect a new tree to the reader
 */
void ForestReader::ReadForestFromFile(TFile *inputFile){
  
  // When reading a forest, we need to check if it is AOD or MiniAOD forest as there are some differences
  // The HiForest tree is renamed to HiForestInfo in MiniAODs, so we can determine the forest type from this.
  TTree* miniAODcheck = (TTree*)inputFile->Get("HiForestInfo/HiForest");
  fIsMiniAOD = !(miniAODcheck == NULL);
  
  // Helper variable for finding the correct tree
  const char *treeName[4] = {"none","none","none","none"};
  
  // Connect a trees from the file to the reader
  fHeavyIonTree = (TTree*)inputFile->Get("hiEvtAnalyzer/HiTree");
  fHltTree = (TTree*)inputFile->Get("hltanalysis/HltTree");
  fSkimTree = (TTree*)inputFile->Get("skimanalysis/HltTree");
  
  // The jet tree has different name in different datasets
  if(fDataType == kPp || fDataType == kPpMC){
    treeName[0] = "ak4CaloJetAnalyzer/t"; // Tree for calo jets
    treeName[1] = "ak4PFJetAnalyzer/t";   // Tree for PF jets
  } else if (fDataType == kPbPb || fDataType == kPbPbMC){
    treeName[0] = "akPu4CaloJetAnalyzer/t";     // Tree for calo jets
    treeName[1] = "akCs4PFJetAnalyzer/t";       // Tree for csPF jets
    treeName[2] = "akPu4PFJetAnalyzer/t";       // Tree for puPF jets
    treeName[3] = "akFlowPuCs4PFJetAnalyzer/t"; // Tree for flow subtracted csPF jets
  }
  
  fJetTree = (TTree*)inputFile->Get(treeName[fJetType]);
  
  Initialize();
}

/*
 * Connect a new tree to the reader
 */
void ForestReader::ReadForestFromFileList(std::vector<TString> fileList){
  TFile *inputFile = TFile::Open(fileList.at(0));
  ReadForestFromFile(inputFile);
}

/*
 * Burn the current forest.
 */
void ForestReader::BurnForest(){
  fHeavyIonTree->Delete();
  fHltTree->Delete();
  fSkimTree->Delete();
  fJetTree->Delete();
}

/*
 * Load an event to memory
 */
void ForestReader::GetEvent(Int_t nEvent){
  fHeavyIonTree->GetEntry(nEvent);
  fJetTree->GetEntry(nEvent);
  fHltTree->GetEntry(nEvent);
  fSkimTree->GetEntry(nEvent);
}

// Getter for number of events in the tree
Int_t ForestReader::GetNEvents() const{
  return fJetPtBranch->GetEntries();
}

// Getter for number of jets in an event
Int_t ForestReader::GetNJets() const{
  return fnJets;
}

// Getter for number of jets in an event
Int_t ForestReader::GetNGeneratorJets() const{
  return fnGenJets;
}

// Getter for jet pT
Float_t ForestReader::GetJetPt(Int_t iJet) const{
  return fJetPtArray[iJet];
}

// Getter for jet phi
Float_t ForestReader::GetJetPhi(Int_t iJet) const{
  return fJetPhiArray[iJet];
}

// Getter for jet eta
Float_t ForestReader::GetJetEta(Int_t iJet) const{
  return fJetEtaArray[iJet];
}

// Getter for jet raw pT
Float_t ForestReader::GetJetRawPt(Int_t iJet) const{
  return fJetRawPtArray[iJet];
}

// Getter for maximum track pT inside a jet
Float_t ForestReader::GetJetMaxTrackPt(Int_t iJet) const{
  return fJetMaxTrackPtArray[iJet];
}

// Getter for generator level jet pT
Float_t ForestReader::GetGeneratorJetPt(Int_t iJet) const{
  return fGenJetPtArray[iJet];
}

// Getter for generator level jet phi
Float_t ForestReader::GetGeneratorJetPhi(Int_t iJet) const{
  return fGenJetPhiArray[iJet];
}

// Getter for generator level jet eta
Float_t ForestReader::GetGeneratorJetEta(Int_t iJet) const{
  return fGenJetEtaArray[iJet];
}

// Getter for vertex z position
Float_t ForestReader::GetVz() const{
  return fVertexZ;
}

// Getter for centrality. CMS has integer centrality bins from 0 to 200, thus division by 2.
Float_t ForestReader::GetCentrality() const{
  return fHiBin/2.0;
}

// Getter for hiBin. Return 1 for negative values (for easier handling of tracking efficiency correction)
Int_t ForestReader::GetHiBin() const{
  if(fHiBin < 0) return 1;
  return fHiBin;
}

// Getter for pT hat
Float_t ForestReader::GetPtHat() const{
  return fPtHat;
}

// Getter for pT hat
Float_t ForestReader::GetEventWeight() const{
  return fEventWeight;
}

// Getter for calorimeter jet filter bit. Always 1 for MC (set in the initializer).
Int_t ForestReader::GetBaseJetFilterBit() const{
  return GetJetFilterBit(fBaseTrigger);
}

// Getter for the selected jet filter bit
Int_t ForestReader::GetJetFilterBit(Int_t iTrigger) const{
  if(iTrigger < 0 || iTrigger >= TriggerHistograms::knTriggerTypes) return -1;
  return fJetFilterBit[iTrigger];
}

// Getter for the prescale value of the chosen trigger
Double_t ForestReader::GetJetTriggerPrescale(Int_t iTrigger) const{
  if(iTrigger < 0 || iTrigger >= TriggerHistograms::knTriggerTypes) return -1;
  return (fJetPrescaleNumerator[iTrigger]*1.0)/fJetPrescaleDenominator[iTrigger];
}

// Getter for primary vertex filter bit. Always 1 for MC (set in the initializer).
Int_t ForestReader::GetPrimaryVertexFilterBit() const{
  return fPrimaryVertexFilterBit;
}

// Getter for beam scraping filter bit. Always 1 for MC and PbPb (set in the initializer).
Int_t ForestReader::GetBeamScrapingFilterBit() const{
  return fBeamScrapingFilterBit;
}

// Getter for HF energy coincidence filter bit. Always 1 for MC and pp (set in the initializer).
Int_t ForestReader::GetHfCoincidenceFilterBit() const{
  return fHfCoincidenceFilterBit;
}

// Getter for cluster compatibility filter bit. Always 1 for MC and pp (set in the initializer).
Int_t ForestReader::GetClusterCompatibilityFilterBit() const{
  return fClusterCompatibilityFilterBit;
}
