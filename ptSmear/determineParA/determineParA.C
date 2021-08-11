#include "../../common/headers.h"
#include "../../common/function.C"
#include "../../common/funUtil.h"

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooFormulaVar.h"
#include "RooHistPdf.h"
#include "RooGenericPdf.h"
#include "RooAddPdf.h"
#include "RooPolynomial.h"
#include "RooChi2Var.h"
#include "RooMinimizer.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooPlot.h"
#include "RooFitResult.h"

using namespace RooFit;
Bool_t  mUseMinos = kFALSE;

const Bool_t  mStorePDF = kFALSE;

const Double_t mTinyNum = 1.e-6;
const Double_t mOffSet = 0.1;

Int_t    mTextFont  = 42;
Double_t mTextSize  = 0.045;
Int_t    mTextColor = 1;

Double_t mMarkerStyle = 20;
Double_t mMarkerSize  = 0.8;

Double_t mTitleSize    = 0.06;
Double_t mXTitleOffset = 0.95;
Double_t mYTitleOffset = 0.95;
Double_t mLabelSize    = 0.05;
Double_t mTickLength   = 0.02;
Int_t    mXNdivisions  = 210;
Int_t    mYNdivisions  = 208;

Int_t    mLineWidth = 1;
Int_t    cohJpsiColor   = kBlue,        cohJpsiStyle = 1;
Int_t    incohJpsiColor = kViolet-1,    incohJpsiStyle = 1;
Int_t    dissoJpsiColor = kMagenta-4,   dissoJpsiStyle = 1;
Int_t    feeddownJpsiColor = kAzure+10, feeddownJpsiStyle = 1;
Int_t    qedColor = kGreen+3, qedStyle = 1;

const Double_t mJpsi = 3.0969, mPsi = 3.686097;

Double_t freject(Double_t *x, Double_t *par){
    if(x[0]>2.9 && x[0]<3.3){
        TF1::RejectPoint();
        return 0;
    }
    return par[0] + par[1]*x[0] + par[2]*pow(x[0],2) + par[3]*pow(x[0],3);
}

void determineParA(Double_t massLow=2.5, Double_t massHi=3.5, Double_t mMaxPt = 0.15)
{
    gStyle->SetOptFit(1111);

    TFile *f = TFile::Open("../../anaData/jpsiHistos/rawSig.root");
    TFile *fTemp = TFile::Open("smearedTemp/smearedTemp.root");

    TString dirName  = "ptSmearCheck";

    system(Form("mkdir -p %s", dirName.Data()));
    system(Form("rm -rf %s/*", dirName.Data()));

    TCanvas* c1 = new TCanvas("c1", "c1", 0, 0, 800, 600);
    setPad(0.12, 0.08, 0.07, 0.13);

    TPDF *ps = new TPDF(Form("%s/DataVsJpsiTemp.pdf", dirName.Data()), 111);
    ps->Off();

    Int_t nColumns = 2;
    Int_t nRaws = 2;
    Int_t nPads = nColumns * nRaws;
    TCanvas* c2 = new TCanvas("c2", "c2", 1200, 900);
    c2->Divide(nColumns, nRaws);
    for(Int_t ipad=0; ipad<nPads; ipad++){
        c2->cd(ipad+1);
        setPad(0.12, 0.08, 0.07, 0.13);
    }

    TH3D *hMvsPtvsRap = (TH3D *)f->Get("hMvsPtvsRap");

    TH1D *hMass = nullptr;
    TH1D *hMass_Rap[nRapBins];

    TH1D *hMass_CohJpsi[nSmearScan], *hMass_CohPsi2S[nSmearScan], *hMass_QED[nSmearScan];
    TH1D *hMass_CohJpsi_Rap[nRapBins][nSmearScan], *hMass_CohPsi2S_Rap[nRapBins][nSmearScan], *hMass_QED_Rap[nRapBins][nSmearScan];
    for(Int_t irap=0; irap<nRapBins; irap++){
        Int_t rapBinLow  = hMvsPtvsRap->GetXaxis()->FindBin(mRapLow[irap] + mTinyNum);
        Int_t rapBinHi   = hMvsPtvsRap->GetXaxis()->FindBin(mRapHi[irap] - mTinyNum);
        Int_t ptBinLow   = 1;
        Int_t ptBinHi    = hMvsPtvsRap->GetYaxis()->FindBin(mMaxPt - mTinyNum);

        hMass_Rap[irap] = (TH1D *)hMvsPtvsRap->ProjectionZ(Form("hMass_RapBin%d", irap), rapBinLow, rapBinHi, ptBinLow, ptBinHi);

        if(irap==0){
            hMass = (TH1D *)hMass_Rap[irap]->Clone("hMass");
            hMass->SetTitle(Form("%1.1f < |y| < %1.1f", mRapLow[nRapBins/2], mRapHi[nRapBins-1]));
        }
        else{
            hMass->Add(hMass_Rap[irap]);
        }

        for(Int_t iscan=0; iscan<nSmearScan; iscan++){
            if(irap==0)
            {
                hMass_CohJpsi[iscan]  = (TH1D *)fTemp->Get(Form("hMass_CohJpsi_Scan%d", iscan));
                hMass_CohPsi2S[iscan] = (TH1D *)fTemp->Get(Form("hMass_CohPsi2S_Scan%d", iscan));
                hMass_QED[iscan]      = (TH1D *)fTemp->Get(Form("hMass_QED_Scan%d", iscan));
            }

            hMass_CohJpsi_Rap[irap][iscan]  = (TH1D *)fTemp->Get(Form("hMass_CohJpsi_Rap%d_Scan%d", irap, iscan));
            hMass_CohPsi2S_Rap[irap][iscan] = (TH1D *)fTemp->Get(Form("hMass_CohPsi2S_Rap%d_Scan%d", irap, iscan));
            hMass_QED_Rap[irap][iscan]      = (TH1D *)fTemp->Get(Form("hMass_QED_Rap%d_Scan%d", irap, iscan));
        }
    }

    //hMass->Rebin(2);
    //for(Int_t irap=0; irap<nRapBins; irap++){
    //    hMass_Rap[irap]->Rebin(4);
    //}

    TF1 *fQED_Rej = new TF1("fQED_Rej", freject, 0, 5, 4);
    TF1 *fQED = new TF1("fQED", "[0] + [1]*x + [2]*x*x +[3]*x*x*x", 0, 5);
    setFun(fQED, kAzure+10, 2, 3);

    hMass->GetXaxis()->SetRangeUser(massLow, massHi);
    hMass->Fit(fQED_Rej, "QR0", "", massLow, massHi);
    fQED->SetParameters(fQED_Rej->GetParameters());

    //for(int i=0; i<4; i++){
    //    cout<<fQED_Rej->GetParameter(i)<<endl;
    //}

    RooRealVar mMass("mMass", "M_{#mu#mu} (GeV)", massLow, massHi);
    RooRealVar cbAlpha("cbAlpha", "cbAlpha", 2, 0, 10);
    RooRealVar cbN("cbN", "cbN", 5, 0, 20);
    RooRealVar jpsiMu("jpsiMu", "jpsiMu", 3.096, 3.0, 3.2);
    RooRealVar gausN("gausN", "gausN", 5, 0, 20);
    RooRealVar jpsiSigma("jpsiSigma", "jpsiSigma", 0.04, 0, 0.12);
    RooRealVar sigmaRatio("sigmaRatio", "sigmaRatio", 1.7, 1, 5);
    RooRealVar nJpsi("nJpsi", "nJpsi", 3.e4, 0, 8.e4);

    RooRealVar mP0("mP0", "mP0", -3.e4, -1.e5, 1.e5);
    RooRealVar mP1("mP1", "mP1",   3e4, -1.e5, 1.e5);
    RooRealVar mP2("mP2", "mP2", -8000, -1.e5, 1.e5);
    RooRealVar mP3("mP3", "mP3",   800, -1.e4, 1.e4);
    RooRealVar nQED("nQED", "nQED", 3.e4, 0, 8.e4);

    RooGenericPdf *jpsiPdf = new RooGenericPdf("jpsiPdf", "jpsiPdf", "ROOT::Math::crystalball_function(mMass,cbAlpha,cbN,jpsiSigma*sigmaRatio,jpsiMu) + gausN*TMath::Gaus(mMass, jpsiMu, jpsiSigma)", RooArgSet(mMass, cbAlpha, cbN, jpsiSigma, sigmaRatio, jpsiMu, gausN));

    RooGenericPdf *qedPdf = new RooGenericPdf("qedPdf", "qedPdf", "mP0 + mP1*mMass + mP2*mMass*mMass + mP3*mMass*mMass*mMass", RooArgSet(mP0, mP1, mP2, mP3, mMass));

    RooAddPdf totPdf("totPdf", "totPdf", RooArgList(*jpsiPdf, *qedPdf), RooArgList(nJpsi, nQED)); 

    RooDataHist dataMass("dataMass", "dataMass", mMass, hMass); 

    //totPdf.fitTo(dataMass,Range(massLow, massHi),Extended(kTRUE),SumW2Error(kTRUE),Hesse(kTRUE),Minos(kTRUE),Save()); 
    totPdf.fitTo(dataMass,Range(massLow, massHi),Extended(kTRUE),SumW2Error(kTRUE),Hesse(kTRUE),Minos(kFALSE),Save()); 

    Int_t nFrameMBins = (massHi - massLow) / hMass->GetBinWidth(1);

    c1->cd();
    RooPlot *frameMass = mMass.frame(Range(massLow, massHi), Title(""), Bins(nFrameMBins));
    dataMass.plotOn(frameMass, MarkerStyle(20), MarkerSize(1), MarkerColor(1), LineColor(1), LineWidth(2), DrawOption("pz"));
    totPdf.plotOn(frameMass, LineColor(2), LineStyle(1), LineWidth(2));
    totPdf.plotOn(frameMass, Components(RooArgSet(*jpsiPdf)), LineColor(kBlue), LineStyle(5), LineWidth(2));
    totPdf.plotOn(frameMass, Components(RooArgSet(*qedPdf)), LineColor(kMagenta), LineStyle(2), LineWidth(2));

    //cout<<endl;
    //cout<<"******** Print frame ********"<<endl;
    //frameMass->Print();
    //cout<<"******** End ********"<<endl;
    //cout<<endl;

    Double_t chi2ndf = frameMass->chiSquare("totPdf_Norm[mMass]_Range[fit_nll_totPdf_dataMass]_NormRange[fit_nll_totPdf_dataMass]", "h_dataMass", 12);

    frameMass->GetYaxis()->SetRangeUser(0, 1.2*hMass->GetMaximum());
    frameMass->Draw() ;
    fQED->Draw("same");
    drawLatex(0.2, 0.84, Form("%1.1f < |y| < %1.1f", mRapLow[nRapBins/2], mRapHi[nRapBins-1]), mTextFont, 0.05, mTextColor);
    drawLatex(0.2, 0.72, Form("#chi^{2}/ndf = %1.1f", chi2ndf), mTextFont, mTextSize, mTextColor);
    drawLatex(0.2, 0.64, Form("N_{J/#psi} = %d #pm %d", TMath::Nint(nJpsi.getVal()), TMath::Nint(nJpsi.getError())), mTextFont, mTextSize, mTextColor);
    drawLatex(0.2, 0.56, Form("N_{QED} = %d #pm %d", TMath::Nint(nQED.getVal()), TMath::Nint(nQED.getError())), mTextFont, mTextSize, mTextColor);

    c1->SaveAs(Form("%s/massSpec_NoTemp.pdf", dirName.Data()));

    //nQED.setConstant(kTRUE);
    mP0.setConstant(kTRUE);
    mP1.setConstant(kTRUE);
    mP2.setConstant(kTRUE);
    mP3.setConstant(kTRUE);

    TH1D *hChi2VsScanIdx = new TH1D("hChi2VsScanIdx", "hChi2VsScanIdx; scanIdx; #chi^{2}", nSmearScan, -0.5, nSmearScan-0.5);

    Int_t nFreePars = 2;
    for(Int_t iscan=0; iscan<nSmearScan; iscan++){
        RooDataHist hCohJpsiRooHist("hCohJpsiRooHist", "hCohJpsiRooHist", mMass, hMass_CohJpsi[iscan]);
        RooHistPdf  cohJpsiPdf("cohJpsiPdf", "cohJpsiPdf", mMass, hCohJpsiRooHist, 2);

        //hMass_QED[iscan]->RebinX(5);
        //RooDataHist hQEDMassRooHist("hQEDMassRooHist", "hQEDMassRooHist", mMass, hMass_QED[iscan]);
        //RooHistPdf  qedPdf("qedPdf", "qedPdf", mMass, hQEDMassRooHist, 2); // RebinX and interpolation order to make the QED pdf smooth 

        RooAddPdf totMassPdf("totMassPdf", "totMassPdf", RooArgList(cohJpsiPdf, *qedPdf), RooArgList(nJpsi, nQED)); 
        //RooAddPdf totMassPdf("totMassPdf", "totMassPdf", RooArgList(cohJpsiPdf, qedPdf), RooArgList(nJpsi, nQED)); 

        RooDataHist dataMass("dataMass", "dataMass", mMass, hMass); 
        if(mUseMinos) totMassPdf.fitTo(dataMass,Range(massLow, massHi),Extended(kTRUE),SumW2Error(kTRUE),Hesse(kTRUE),Minos(kTRUE),Save()); 
        else          totMassPdf.fitTo(dataMass,Range(massLow, massHi),Extended(kTRUE),SumW2Error(kTRUE),Hesse(kTRUE),Minos(kFALSE),Save());

        c2->cd(iscan%nPads+1);
        RooPlot *frameMass = mMass.frame(Range(massLow, massHi), Title(""), Bins(nFrameMBins));
        dataMass.plotOn(frameMass, MarkerStyle(20), MarkerSize(1), MarkerColor(1), LineColor(1), LineWidth(2), DrawOption("pz"));
        totMassPdf.plotOn(frameMass, LineColor(2), LineStyle(1), LineWidth(2));
        totMassPdf.plotOn(frameMass, Components(RooArgSet(cohJpsiPdf)), LineColor(kBlue), LineStyle(5), LineWidth(2));
        totMassPdf.plotOn(frameMass, Components(RooArgSet(*qedPdf)), LineColor(kMagenta), LineStyle(2), LineWidth(2));
        //totMassPdf.plotOn(frameMass, Components(RooArgSet(qedPdf)), LineColor(kMagenta), LineStyle(2), LineWidth(2));

        //cout<<endl;
        //cout<<"******** Print frame ********"<<endl;
        //frameMass->Print();
        //cout<<"******** End ********"<<endl;
        //cout<<endl;

        Int_t   ndf = hMass->GetXaxis()->FindBin(massHi-mTinyNum) - hMass->GetXaxis()->FindBin(massLow+mTinyNum) + 1 - nFreePars;
        Double_t chi2ndf = frameMass->chiSquare("totMassPdf_Norm[mMass]_Range[fit_nll_totMassPdf_dataMass]_NormRange[fit_nll_totMassPdf_dataMass]", "h_dataMass", nFreePars);
        Double_t chi2 = chi2ndf * ndf;

        hChi2VsScanIdx->SetBinContent(iscan+1, chi2);
        hChi2VsScanIdx->SetBinError(iscan+1, 0);

        Double_t par0 = mInitPar0 + iscan*mSmearStep;

        frameMass->GetYaxis()->SetRangeUser(0, 1.2*hMass->GetMaximum());
        frameMass->DrawClone() ;
        drawLatex(0.28, 0.95, Form("scanIdx = %d (par0 = %5.5f)", iscan, par0), mTextFont, 0.06, 4);
        drawLatex(0.2, 0.84, Form("%1.1f < |y| < %1.1f", mRapLow[nRapBins/2], mRapHi[nRapBins-1]), mTextFont, 0.05, mTextColor);
        drawLatex(0.2, 0.72, Form("#chi^{2}/ndf = %1.2f/%d (%1.2f)", chi2, ndf, chi2ndf), mTextFont, mTextSize, mTextColor);
        drawLatex(0.2, 0.64, Form("N_{J/#psi} = %d #pm %d", TMath::Nint(nJpsi.getVal()), TMath::Nint(nJpsi.getError())), mTextFont, mTextSize, mTextColor);
        drawLatex(0.2, 0.56, Form("N_{QED} = %d #pm %d", TMath::Nint(nQED.getVal()), TMath::Nint(nQED.getError())), mTextFont, mTextSize, mTextColor);

        if(iscan%nPads == nPads-1) pdfAction(c2, ps);
    }
    if(nSmearScan%nPads != 0) pdfAction(c2, ps, kTRUE);
    else{
        ps->On();
        ps->Close();
    }

    TFile *fOut = new TFile(Form("%s/chi2VsScanIdx.root", dirName.Data()), "recreate");
    fOut->cd();
    hChi2VsScanIdx->Write();
    fOut->Close();

    cout << "End of program !" << endl;
}
