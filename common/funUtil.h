#include "/Users/syang/work/run2/upcDimuon/common/constants.h"

// *** Initialization ***
TH1D* hZDC[nDirs];
TF1*  multiGaus[nDirs];
TF1*  singleGaus[nDirs][nNeus-1];

TH3D* hPosMu3DMthEff;
TH3D* hNegMu3DMthEff;
TH3D* hPosMu3DTrigEff;
TH3D* hNegMu3DTrigEff;

TH2D* hLumvsLSvsRun_UPC;
TH2D* hLumvsLSvsRun_ZB;
// **********************

// Details of matrixIdx can be found in the analysis note
const Int_t matrixIdx[nNeus][nNeus] = {
    {0, 1, 3},
    {2, 5, 6},
    {4, 7, 8}
};

const Int_t nSces = 9;
Double_t    corrMatrix[nSces][nSces];

Bool_t init(TString hfVetoType="Default"){
    TFile *fNeu = TFile::Open("/Users/syang/work/run2/upcDimuon/neuDecouple/neuDecouple/neutronDecouple.root");
    if(!fNeu->IsOpen()){
        cout<<"Failed to open neutronDecouple.root !"<<endl;
        return kFALSE;
    }
    else{
        for(Int_t idir = 0; idir < nDirs; idir++){
            hZDC[idir]       = (TH1D *)fNeu->Get(Form("hZDC%s", mDir[idir].Data()));
            multiGaus[idir]  = (TF1 *)fNeu->Get(Form("multiGaus_%s", mDir[idir].Data())); 

            for(Int_t ineu = 0; ineu < nNeus - 1; ineu++){
                singleGaus[idir][ineu]  = (TF1 *)fNeu->Get(Form("singleGaus_%s%d", mDir[idir].Data(), ineu));
            }
        }
    }

    TFile *fEff = TFile::Open("/Users/syang/work/run2/upcDimuon/simulation/effPlots_GammaGamma/3DMthEffAndTrigEff.root");
    //TFile *fEff = TFile::Open("/Users/syang/work/run2/upcDimuon/simulation/effPlots_GammaGamma_woVtxAndNtkHPSel/3DMthEffAndTrigEff.root");
    if(!fEff->IsOpen()){
        cout<<"Failed to open 3-D single efficiencies !"<<endl;
        return kFALSE;
    }
    else{
        hPosMu3DMthEff  = (TH3D *)fEff->Get("hPosMu3DMthEff");
        hNegMu3DMthEff  = (TH3D *)fEff->Get("hNegMu3DMthEff");
        hPosMu3DTrigEff = (TH3D *)fEff->Get("hPosMu3DTrigEff");
        hNegMu3DTrigEff = (TH3D *)fEff->Get("hNegMu3DTrigEff");

        //hPosMu3DMthEff  = (TH3D *)fEff->Get("hPosMu3DMthEff_RebPhi");
        //hNegMu3DMthEff  = (TH3D *)fEff->Get("hNegMu3DMthEff_RebPhi");
        //hPosMu3DTrigEff = (TH3D *)fEff->Get("hPosMu3DTrigEff_RebPhi");
        //hNegMu3DTrigEff = (TH3D *)fEff->Get("hNegMu3DTrigEff_RebPhi");
    }

    // For generating pileup matrix in Zero Bias data
    TFile *fLumi_UPC = TFile::Open("/Users/syang/work/run2/upcDimuon/anaUPCLumi/InstantLumvsLSNbvsRunNb_singleMuUPC.root");
    if(!fLumi_UPC->IsOpen()){
        cout<<"Failed to open UPC instantaneous luminosity file !"<<endl;
        return kFALSE;
    }
    else{
        hLumvsLSvsRun_UPC = (TH2D *)fLumi_UPC->Get("hLumvsLSvsRun_UPC");
    }

    TFile *fLumi_ZB = TFile::Open("/Users/syang/work/run2/upcDimuon/anaZBLumi/InstantLumvsLSNbvsRunNb_ZeroBias.root");
    if(!fLumi_ZB->IsOpen()){
        cout<<"Failed to open Zero Bias instantaneous luminosity file !"<<endl;
        return kFALSE;
    }
    else{
        hLumvsLSvsRun_ZB = (TH2D *)fLumi_ZB->Get("hLumvsLSvsRun_ZB");
    }

    memset(corrMatrix, 0, sizeof(corrMatrix));

    ifstream inData(Form("/Users/syang/work/run2/upcDimuon/zerobias/corrMatrix_%s.dat", hfVetoType.Data()));
    if(!inData.is_open()){
        cout<<"Failed to open correction matrix data !"<<endl;
        return kFALSE;
    }

    Int_t icount=0;
    Double_t ele;
    while(inData>>ele){
        Int_t iraw = icount/nSces;
        Int_t icol = icount%nSces;

        corrMatrix[iraw][icol] = ele;

        icount++;
    }

    //for(Int_t i=0; i<nSces; i++){
    //    for(Int_t j=0; j<nSces; j++){
    //        cout<< std::right << setw(15) << corrMatrix[i][j];
    //    }
    //    cout<<endl;
    //}

    return kTRUE;
}

Double_t grabNeutronProb(Int_t dirIdx, Double_t zdc, Int_t neuNum)
{
    Int_t    binIdx     = hZDC[dirIdx]->GetXaxis()->FindBin(zdc);
    Double_t binCenter  = hZDC[dirIdx]->GetXaxis()->GetBinCenter(binIdx);
    Double_t binContent = hZDC[dirIdx]->GetBinContent(binIdx);

    Double_t prob = -1;
    if(neuNum == 0){
        if(zdc<=mZdcFitLow[dirIdx]){
            prob = (binContent - multiGaus[dirIdx]->Eval(binCenter)) / binContent;
        }
        else{
            prob = 0;
        }
    }
    else if(neuNum == 1 || neuNum == 2){
        if(zdc < mZdcFitLow[dirIdx] || zdc > mZdcFitHi[dirIdx]){
            prob = singleGaus[dirIdx][neuNum-1]->Eval(binCenter) / binContent;
        }
        else{ // mZdcFitLow <= zdc <= mZdcFitHi using fit as denominator
            prob = singleGaus[dirIdx][neuNum-1]->Eval(zdc) / multiGaus[dirIdx]->Eval(zdc);
        }
    }
    else if(neuNum == 3){
        if(zdc < mZdcFitLow[dirIdx]){
            prob = singleGaus[dirIdx][neuNum-1]->Eval(binCenter) / binContent;
        }
        else if(zdc > mZdcFitHi[dirIdx]){
            prob = (binContent - singleGaus[dirIdx][0]->Eval(binCenter) - singleGaus[dirIdx][1]->Eval(binCenter)) / binContent;
        }
        else{ // mZdcFitLow <= zdc <= mZdcFitHi using fit as denominator
            prob = (multiGaus[dirIdx]->Eval(zdc) - singleGaus[dirIdx][0]->Eval(zdc) - singleGaus[dirIdx][1]->Eval(zdc)) / multiGaus[dirIdx]->Eval(zdc);
        }
    }
    else{
        cout<<"The neutron number should be 0, 1, 2, 3"<<endl;
    }

    return prob;
}
