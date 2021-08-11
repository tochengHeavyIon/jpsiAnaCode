#include "constants.h"

TF1 *funPtMeanShift;
TF1 *funRawPtRes;
TF1 *funTunedPtRes;

Bool_t init(){
    std::string cwd = gSystem->ExpandPathName(gSystem->pwd());
    auto comDir = cwd.substr(0,cwd.find("jpsiAnaCode"))+"jpsiAnaCode/common/";

    TFile *fRawPtRes = TFile::Open(Form("%s/rawPtRes.root", comDir.c_str()));
    funPtMeanShift = (TF1 *)fRawPtRes->Get("funPtMeanShift");
    funRawPtRes = (TF1 *)fRawPtRes->Get("funRawPtRes");

    funTunedPtRes = new TF1("funTunedPtRes", "sqrt([0]*[0]/x/x+[1]*[1])", 0, 5);
    funTunedPtRes->SetParameters(mPar0, funRawPtRes->GetParameter(1));

    return kTRUE;
}

Double_t trigAcc(Double_t *x, Double_t *par){
    par = NULL;

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
