const Double_t Mmuon = 0.1056583745;
const Double_t Mpion = 0.13957018;
const Double_t Mkaon = 0.493677;
const Double_t Mproton = 0.938272081;
const Double_t Melectron = 0.0005109989461;

const Double_t PI = TMath::Pi();

// 0 HLT_HIL1DoubleMuOpen_OS_Centrality_40_100 (HIDoubleMuon, HIDoubleMuonPsiPeri) -> Events with at least two opposite electric-charge L1 muons with loose trigger quality (requires at least two muon stations with measurements) within a centrality percentile range between 40-100%
// 1 HLT_HIL1DoubleMuOpen_Centrality_50_100 (HIDoubleMuon, HIDoubleMuonPsiPeri) -> Events with at least two L1 muons with loose trigger quality (requires at least two muon stations with measurements) within a centrality percentile range between 50-100%
// 4 HLT_HIUPC_DoubleMu0_NotMBHF2AND (HIForward) -> Events with at least two L1 muons with medium trigger quality  (requires at least three muon stations with measurements) with low energy in the HF calorimeters (NotMBHF2AND means !(HF+ > 15 ADC Counts  AND  HF- > 19 ADC Counts))
// 5 HLT_HIL1MuOpen_Centrality_80_100 (HIForward) -> Events with at least one L1 muon with loose trigger quality  (requires at least two muon stations with measurements) within a centrality percentile range between 80-100%
// 7 HLT_HIUPC_SingleMuOpen_NotMBHF2AND -> Events with at least one L1 muon with loose trigger quality  (requires at least two muon stations with measurements) within low energy in the HF calorimeters (NotMBHF2AND means !(HF+ > 15 ADC Counts  AND  HF- > 19 ADC Counts))
const Int_t nTrigs = 9;

//const Int_t   trigIdx = 4;
//const TString trigName = "DoubleMuUPC";

const Int_t   trigIdx = 7;
const TString trigName = "SingleMuUPC";

const Int_t zbTrigIdx = 9;
const Int_t emptyBXTrigIdx = 10;

const UInt_t   mRunNbCut = 326776; // get this number from Quan (L = 1520.302020023633 ub^{-1} for runID >= 326776)
const Double_t mCMSLum = 1520.302020023633; // ub^{-1} for HLT_HIUPC_SingleMuOpen_NotMBHF2AND_v with runID >= 326776 recorded


const Int_t  mCentralityCut = 180;
const Int_t  mNtrkofflineCut = 2;

const Int_t nDirs = 2; // 0 - Plus; 1 - Minus;
const TString  mDir[nDirs] = {"Plus", "Minus"};
const Double_t mHFsumETCut[nDirs] = {12, 12};
const Double_t mZdcFitLow[nDirs]  = {4.2e3, 6.e3};
const Double_t mZdcFitHi[nDirs]   = {25.e3, 37.5e3};

const Int_t nNeus = 3;
const Int_t nPts = nNeus*(nNeus+1)/2;

const Double_t mNeuZDCLow[nDirs][nNeus] = {
    {0, 4.2e3, 10.e3},
    {0, 6.0e3, 16.e3}
};
const Double_t mNeuZDCHi[nDirs][nNeus] = {
    {4.2e3, 10.e3, 5.e5},
    {6.0e3, 16.e3, 6.e5}
};

// expected neutron peak poisition
// Plus(*e3):   7.33, 14.66, 21.99, 29.32, 36.65, 43.98, 51.31, 58.64,  65.97,  73.3
// Minus(*e3): 11.42, 22.84, 34.26, 45.68, 57.10, 68.52, 79.94, 91.36, 102.78, 114.2
const Int_t  nRegions = 8;
const Double_t mNeuZDC[nDirs][nRegions+1] = {
    {0, 4.2e3, 10.7e3, 17.0e3, 24.9e3, 32.3e3, 46.9e3, 61.6e3,  80.e3}, 
    {0, 6.e3,  17.3e3, 27.0e3, 38.8e3, 50.2e3, 73.1e3, 95.9e3, 120.e3}
};

const Double_t mTwoNeuZDCLow[nDirs] = {11e3, 17.5e3};
const Double_t mTwoNeuZDCHi[nDirs]  = {17e3, 27.0e3};
const Double_t mThreeNeuZDCLow[nDirs] = {21e3, 32e3};
const Double_t mThreeNeuZDCHi[nDirs]  = {27e3, 40e3};

const Int_t    nMBins = 3;
const Double_t massLow[nMBins] = { 8, 20, 30};
const Double_t massHi[nMBins] = { 20, 30, 60};

//track quality cuts
const Double_t mPtCut = 3.5;
const Double_t mEtaCut = 2.4;

const Double_t mVtxProbCut = 1.e-6;
const Double_t mPairYCut = 2.4;

const Double_t mAlphaCut = 6.e-3;
