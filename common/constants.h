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
const Int_t nTrigs = 8;

//const Int_t   trigIdx = 4;
//const TString trigName = "DoubleMuUPC";

const Int_t   trigIdx = 7;
const TString trigName = "SingleMuUPC";

const UInt_t mRunNbCut = 326776; // get this number from Quan (L = 1570.6796 ub^{-1} for runID >= 326776)
const Double_t mCMSLum = 1570.6796; // ub^{-1} for HLT_HIUPC_SingleMuOpen_NotMBHF2AND_v with runID >= 326776
const Int_t  mCentralityCut = 180;
const Int_t  mNtrkofflineCut = 2;

const Int_t nDirs = 2; // 0 - Plus; 1 - Minus;
const TString  mDir[nDirs] = {"Plus", "Minus"};
const Double_t mHFsumETCut[nDirs] = {12, 12};
const Double_t mZdcFitLow[nDirs]  = {4.2e3, 6.e3};
const Double_t mZdcFitHi[nDirs]   = {25.e3, 37.5e3};

//const Int_t nNeus = 3;
const Int_t nNeus = 2;
const Int_t nNeuMults = nNeus*(nNeus+1)/2;

//const Double_t mNeuZDCLow[nDirs][nNeus] = {
//    {0, 4.2e3, 10.e3},
//    {0, 6.0e3, 16.e3}
//};
//const Double_t mNeuZDCHi[nDirs][nNeus] = {
//    {4.2e3, 10.e3, 1.e6},
//    {6.0e3, 16.e3, 1.e6}
//};

const Double_t mNeuZDCLow[nDirs][nNeus] = {
    {0, 4.2e3},
    {0, 6.0e3}
};
const Double_t mNeuZDCHi[nDirs][nNeus] = {
    {4.2e3, 1.e6},
    {6.0e3, 1.e6}
};

const Double_t mTwoNeuZDCLow[nDirs] = {11e3, 17.5e3};
const Double_t mTwoNeuZDCHi[nDirs]  = {17e3, 27.0e3};

//const Int_t    nRapBins = 16;
//const Double_t mRapLow[nRapBins] = {-2.4, -2.3, -2.2, -2.1, -2.0, -1.9, -1.8, -1.7, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3};
//const Double_t mRapHi[nRapBins]  = {-2.3, -2.2, -2.1, -2.0, -1.9, -1.8, -1.7, -1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4};

const Int_t    nRapBins = 12;
const Double_t mRapLow[nRapBins] = {-2.4, -2.2, -2.1, -2.0, -1.9, -1.8, 1.6, 1.8, 1.9, 2.0, 2.1, 2.2};
const Double_t mRapHi[nRapBins]  = {-2.2, -2.1, -2.0, -1.9, -1.8, -1.6, 1.8, 1.9, 2.0, 2.1, 2.2, 2.4};

//const Int_t    nDiffRapBins = 6;
//const Double_t mDiffRapLow[nRapBins] = {-2.4, -2.1, -1.8, 1.6, 1.8, 2.1};
//const Double_t mDiffRapHi[nRapBins]  = {-2.1, -1.8, -1.6, 1.8, 2.1, 2.4};

const Int_t    nDiffRapBins = 4;
const Double_t mDiffRapLow[nRapBins] = {-2.4, -2.0, 1.6, 2.0};
const Double_t mDiffRapHi[nRapBins]  = {-2.0, -1.6, 2.0, 2.4};

const Double_t mMassLow4MuonAccStudy = 2;
const Double_t mMassHi4MuonAccStudy  = 5;

const Double_t mLowMassBandLow = 2.75;
const Double_t mLowMassBandHi  = 2.9;
//const Double_t mJpsiMassLow = 2.85;
//const Double_t mJpsiMassHi  = 3.35;
const Double_t mJpsiMassLow = 2.95;
const Double_t mJpsiMassHi  = 3.25;
const Double_t mHiMassBandLow = 3.3;
const Double_t mHiMassBandHi  = 3.45;

const Double_t mPairYCut = 2.4;
const Double_t mAlphaCut = 6.e-3;

const Int_t    nSmearScan = 200;
const Double_t mInitPar0  = 0.01;
const Double_t mSmearStep  = 5.e-5;
const Double_t mPar0 = 0.0145;
