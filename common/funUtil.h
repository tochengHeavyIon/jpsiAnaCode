#include "constants.h"

// *** Initialization ***
TH3D* hPosMu3DMthEff;
TH3D* hNegMu3DMthEff;
TH3D* hPosMu3DTrigEff;
TH3D* hNegMu3DTrigEff;
// **********************

//// Details of matrixIdx can be found in the analysis note
//const Int_t matrixIdx[nNeus][nNeus] = {
//    {0, 1, 3},
//    {2, 5, 6},
//    {4, 7, 8}
//};

//const Int_t nSces = 9;
//Double_t    corrMatrix[nSces][nSces];

TF1 *funPtMeanShift;
TF1 *funRawPtRes;
TF1 *funTunedPtRes;

// selType represents the dissociative event selection criteria in zerobias sample
//Bool_t init(TString selType="Default"){
Bool_t init(){
    TFile *fRawPtRes = TFile::Open("/Users/syang/work/run2/upcJpsi/ptSmear/rawPtRes/ptRes/rawPtRes.root");
    funPtMeanShift = (TF1 *)fRawPtRes->Get("funPtMeanShift");
    funRawPtRes = (TF1 *)fRawPtRes->Get("funRawPtRes");

    funTunedPtRes = new TF1("funTunedPtRes", "sqrt([0]*[0]/x/x+[1]*[1])", 0, 5);
    funTunedPtRes->SetParameters(mPar0, funRawPtRes->GetParameter(1));

    //TFile *fEff = TFile::Open("/Users/syang/work/run2/upcDimuon/simulation/effPlots_GammaGamma/3DMthEffAndTrigEff.root");
    ////TFile *fEff = TFile::Open("/Users/syang/work/run2/upcDimuon/simulation/effPlots_GammaGamma_woVtxAndNtkHPSel/3DMthEffAndTrigEff.root");
    //if(!fEff->IsOpen()){
    //    cout<<"Failed to open 3-D single efficiencies !"<<endl;
    //    return kFALSE;
    //}
    //else{
    //    hPosMu3DMthEff  = (TH3D *)fEff->Get("hPosMu3DMthEff");
    //    hNegMu3DMthEff  = (TH3D *)fEff->Get("hNegMu3DMthEff");
    //    hPosMu3DTrigEff = (TH3D *)fEff->Get("hPosMu3DTrigEff");
    //    hNegMu3DTrigEff = (TH3D *)fEff->Get("hNegMu3DTrigEff");

    //    //hPosMu3DMthEff  = (TH3D *)fEff->Get("hPosMu3DMthEff_RebPhi");
    //    //hNegMu3DMthEff  = (TH3D *)fEff->Get("hNegMu3DMthEff_RebPhi");
    //    //hPosMu3DTrigEff = (TH3D *)fEff->Get("hPosMu3DTrigEff_RebPhi");
    //    //hNegMu3DTrigEff = (TH3D *)fEff->Get("hNegMu3DTrigEff_RebPhi");
    //}

    //memset(corrMatrix, 0, sizeof(corrMatrix));

    //ifstream inData(Form("/Users/syang/work/run2/upcDimuon/zerobias/corrMatrix_%s.dat", selType.Data()));
    //if(!inData.is_open()){
    //    cout<<"Failed to open correction matrix data !"<<endl;
    //}

    //Int_t icount=0;
    //Double_t ele;
    //while(inData>>ele){
    //    Int_t iraw = icount/nSces;
    //    Int_t icol = icount%nSces;

    //    corrMatrix[iraw][icol] = ele;

    //    icount++;
    //}

    //for(Int_t i=0; i<nSces; i++){
    //    for(Int_t j=0; j<nSces; j++){
    //        cout<< std::right << setw(15) << corrMatrix[i][j];
    //    }
    //    cout<<endl;
    //}

    return kTRUE;
}

Double_t trigAcc(Double_t *x, Double_t *par){
    par = NULL;

    //const Int_t nPts = 10;
    //Double_t mEta[nPts] = {-2.4, -2.1, -1.6, -1.4, -1.0, 1.0, 1.4, 1.6, 2.1, 2.4};
    //Double_t mPt[nPts]  = { 1.2,  1.2,  2.1,  2.1,  3.3, 3.3, 2.1, 2.1, 1.2, 1.2};

    const Int_t nPts = 14;
    Double_t mEta[nPts] = {-2.4, -2.1, -1.65, -1.45, -1.1, -0.3, -0.3,  0.3, 0.3, 1.1, 1.45, 1.65, 2.1, 2.4};
    Double_t mPt[nPts]  = { 1.2,  1.2,  2.15,  2.15,  3.3,  3.3, 3.45, 3.45, 3.3, 3.3, 2.15, 2.15, 1.2, 1.2};

    Int_t iseg = -1;
    for(Int_t i=0; i<nPts-1; i++){
        if(x[0]>=mEta[i] && x[0]<mEta[i+1]){
            iseg = i;
            break;
        }
    }

    if(x[0]==mEta[nPts-1]) iseg = nPts - 2;

    if(iseg<0) return 999999.;

    Double_t mSlope = (mPt[iseg+1] - mPt[iseg]) / (mEta[iseg+1] - mEta[iseg]);
    Double_t mPtTh  = mSlope * (x[0] - mEta[iseg]) + mPt[iseg];

    return mPtTh;
}
TF1 *fTrigAcc = new TF1("fTrigAcc", trigAcc, -2.5, 2.5, 0);

Double_t trkAcc(Double_t *x, Double_t *par){
    par = NULL;

    const Int_t nPts = 10;
    Double_t mEta[nPts] = {-2.4, -1.7,  -1.3, -1.3, -1.0, 1.0, 1.3,  1.3, 1.7, 2.4};
    Double_t mPt[nPts]  = { 1.0,  1.0,  1.53,  2.1,  3.3, 3.3, 2.1, 1.53, 1.0, 1.0};

    Int_t iseg = -1;
    for(Int_t i=0; i<nPts-1; i++){
        if(x[0]>=mEta[i] && x[0]<mEta[i+1]){
            iseg = i;
            break;
        }
    }

    if(x[0]==mEta[nPts-1]) iseg = nPts - 2;

    if(iseg<0) return 999999.;

    Double_t mSlope = (mPt[iseg+1] - mPt[iseg]) / (mEta[iseg+1] - mEta[iseg]);
    Double_t mPtTh  = mSlope * (x[0] - mEta[iseg]) + mPt[iseg];

    return mPtTh;
}
TF1 *fTrkAcc = new TF1("fTrkAcc", trkAcc, -2.5, 2.5, 0);

Double_t oldTrkAcc(Double_t *x, Double_t *par){ // slightly tighter
    par = NULL;

    const Int_t nPts = 10;
    Double_t mEta[nPts] = {-2.4, -1.7, -1.4, -1.4, -1.0, 1.0, 1.4, 1.4, 1.7, 2.4};
    Double_t mPt[nPts]  = { 1.0,  1.0,  1.4,  2.1,  3.3, 3.3, 2.1, 1.4, 1.0, 1.0};

    Int_t iseg = -1;
    for(Int_t i=0; i<nPts-1; i++){
        if(x[0]>=mEta[i] && x[0]<mEta[i+1]){
            iseg = i;
            break;
        }
    }

    if(x[0]==mEta[nPts-1]) iseg = nPts - 2;

    if(iseg<0) return 999999.;

    Double_t mSlope = (mPt[iseg+1] - mPt[iseg]) / (mEta[iseg+1] - mEta[iseg]);
    Double_t mPtTh  = mSlope * (x[0] - mEta[iseg]) + mPt[iseg];

    return mPtTh;
}
TF1 *fOldTrkAcc = new TF1("fOldTrkAcc", oldTrkAcc, -2.5, 2.5, 0);
