#include "../common/headers.h"
#include "../common/function.C"
#include "../common/funUtil.h"

const Int_t    mFont = 42;
const Double_t mTinyNum = 1.e-6;

void plotAcc(Bool_t  mStorePDF = kTRUE)
{
    gStyle->SetOptFit(1111);
    gStyle->SetLabelOffset(0.005,"z");

    TFile *f = TFile::Open("jpsiHistos/rawSig.root");

    TString dirName  = "accStudy";
    system(Form("mkdir -p %s", dirName.Data()));
    system(Form("rm -rf %s/*", dirName.Data()));

    TCanvas* c1 = new TCanvas("c1", "c1", 0, 0, 800, 600);
    setPad(0.12, 0.12, 0.07, 0.13);

    Int_t nColumns = 2;
    Int_t nRaws = 2;
    Int_t nPads = nColumns * nRaws;
    TCanvas* c2 = new TCanvas("c2", "c2", 1200, 900);
    c2->Divide(nColumns, nRaws);
    for(Int_t ipad=0; ipad<nPads; ipad++){
        c2->cd(ipad+1);
        setPad(0.12, 0.12, 0.08, 0.13);
    }

    TH3D *hMuPtvsEtavsRap          = (TH3D *)f->Get("hMuPtvsEtavsRap");
    TH3D *hTrigMuPtvsEtavsRap      = (TH3D *)f->Get("hTrigMuPtvsEtavsRap");
    TH3D *hNegMuPtvsPosMuPtvsRap   = (TH3D *)f->Get("hNegMuPtvsPosMuPtvsRap");
    TH3D *hNegMuEtavsPosMuEtavsRap = (TH3D *)f->Get("hNegMuEtavsPosMuEtavsRap");
    TH3D *hDeltaPtvsDeltaEtavsRap  = (TH3D *)f->Get("hDeltaPtvsDeltaEtavsRap");

    fTrigAcc->SetNpx(500);
    setFun(fTrigAcc, kViolet, 2, 1);

    fTrkAcc->SetNpx(500);
    setFun(fTrkAcc, 2, 2, 2);

    fOldTrkAcc->SetNpx(500);
    setFun(fOldTrkAcc, 1, 2, 2);

    TLegend *leg = new TLegend(0.32, 0.15, 0.54, 0.25);
    leg->SetFillStyle(0);
    leg->SetFillColor(10);
    leg->SetTextSize(0.05);
    leg->SetTextFont(mFont);
    leg->AddEntry(fTrigAcc, "Trigger Muon's Acc.", "l");
    leg->AddEntry(fTrkAcc, "Soft Muon's Acc.", "l");

    c1->cd();
    gPad->SetLogz(1);
    TH2D *hMuPtvsEta = (TH2D *)hMuPtvsEtavsRap->Project3D("hMuPtvsEta_zy");
    hMuPtvsEta->GetXaxis()->SetTitleOffset(0.8);
    hMuPtvsEta->GetYaxis()->SetTitleOffset(0.85);
    hMuPtvsEta->GetYaxis()->SetRangeUser(0.8, 4);
    hMuPtvsEta->Draw("colz");
    fTrigAcc->Draw("same");
    fTrkAcc->Draw("same");
    //fOldTrkAcc->Draw("same");
    drawLatex(0.39, 0.95, "Daughter #mu", mFont, 0.07, 1);
    drawLatex(0.38, 0.32, "2 < M_{#mu#mu} < 5 GeV", mFont, 0.05, 1);
    drawLatex(0.38, 0.38, "1.6 < |y_{#mu#mu}| < 2.4", mFont, 0.05, 1);
    drawLatex(0.72, 0.82, "Data", mFont, 0.07, 1);
    leg->Draw("same");
    if(mStorePDF) c1->SaveAs(Form("%s/JpsiAnaAcc_DaughterMu.pdf", dirName.Data()));
    else          c1->SaveAs(Form("%s/JpsiAnaAcc_DaughterMu.png", dirName.Data()));

    c1->cd();
    gPad->SetLogz(1);
    TH2D *hTrigMuPtvsEta = (TH2D *)hTrigMuPtvsEtavsRap->Project3D("hTrigMuPtvsEta_zy");
    hTrigMuPtvsEta->GetXaxis()->SetTitleOffset(0.8);
    hTrigMuPtvsEta->GetYaxis()->SetTitleOffset(0.85);
    hTrigMuPtvsEta->GetYaxis()->SetRangeUser(0.8, 4);
    hTrigMuPtvsEta->Draw("colz");
    drawLatex(0.40, 0.95, "Trigger #mu", mFont, 0.07, 4);
    drawLatex(0.38, 0.32, "2 < M_{#mu#mu} < 5 GeV", mFont, 0.05, 1);
    drawLatex(0.38, 0.38, "1.6 < |y_{#mu#mu}| < 2.4", mFont, 0.05, 1);
    drawLatex(0.72, 0.82, "Data", mFont, 0.07, 1);
    fTrigAcc->Draw("same");
    fTrkAcc->Draw("same");
    //fOldTrkAcc->Draw("same");
    leg->Draw("same");
    if(mStorePDF) c1->SaveAs(Form("%s/JpsiAnaAcc_TriggerMu.pdf", dirName.Data()));
    else          c1->SaveAs(Form("%s/JpsiAnaAcc_TriggerMu.png", dirName.Data()));

    leg->SetX1NDC(0.15);
    leg->SetX2NDC(0.40);
    leg->SetY1NDC(0.76);
    leg->SetY2NDC(0.88);
    leg->SetTextSize(0.06);

    TH2D *hMuPtvsEta_Rap[nRapBins];
    TH2D *hTrigMuPtvsEta_Rap[nRapBins];
    TH2D *hNegMuPtvsPosMuPt_Rap[nRapBins];
    TH2D *hNegMuEtavsPosMuEta_Rap[nRapBins];
    TH2D *hDeltaPtvsDeltaEta_Rap[nRapBins];
    for(Int_t irap=0; irap<nRapBins; irap++){
        Int_t rapBinLow  = hMuPtvsEtavsRap->GetXaxis()->FindBin(mRapLow[irap] + mTinyNum);
        Int_t rapBinHi   = hMuPtvsEtavsRap->GetXaxis()->FindBin(mRapHi[irap] - mTinyNum);

        c2->cd(1);
        gPad->SetLogz(1);
        //hMuPtvsEtavsRap->GetXaxis()->SetRange(rapBinLow, rapBinHi);
        hMuPtvsEtavsRap->GetXaxis()->SetRangeUser(mRapLow[irap], mRapHi[irap]);
        hMuPtvsEta_Rap[irap] = (TH2D *)hMuPtvsEtavsRap->Project3D(Form("hMuPtvsEta_RapBin%d_zy", irap));
        hMuPtvsEta_Rap[irap]->GetYaxis()->SetTitleOffset(0.85);
        hMuPtvsEta_Rap[irap]->GetYaxis()->SetRangeUser(0.8, 4);
        hMuPtvsEta_Rap[irap]->Draw("colz");
        drawLatex(0.39, 0.95, "Daughter #mu", mFont, 0.07, 1);
        drawLatex(0.35, 0.28, Form("%1.1f < y_{#mu#mu} < %1.1f", mRapLow[irap], mRapHi[irap]), mFont, 0.07, 1);
        drawLatex(0.34, 0.20, "2 < M_{#mu#mu} < 5 GeV", mFont, 0.07, 1);
        drawLatex(0.72, 0.82, "Data", mFont, 0.08, 1);
        fTrigAcc->Draw("same");
        fTrkAcc->Draw("same");
        leg->Draw("same");

        c2->cd(2);
        gPad->SetLogz(1);
        hTrigMuPtvsEtavsRap->GetXaxis()->SetRangeUser(mRapLow[irap], mRapHi[irap]);
        hTrigMuPtvsEta_Rap[irap] = (TH2D *)hTrigMuPtvsEtavsRap->Project3D(Form("hTrigMuPtvsEta_RapBin%d_zy", irap));
        hTrigMuPtvsEta_Rap[irap]->SetTitle(Form("%1.1f < y < %1.1f", mRapLow[irap], mRapHi[irap]));
        hTrigMuPtvsEta_Rap[irap]->GetYaxis()->SetTitleOffset(0.85);
        hTrigMuPtvsEta_Rap[irap]->GetYaxis()->SetRangeUser(0.8, 4);
        hTrigMuPtvsEta_Rap[irap]->Draw("colz");
        drawLatex(0.40, 0.95, "Trigger #mu", mFont, 0.07, 4);
        fTrigAcc->Draw("same");
        fTrkAcc->Draw("same");

        c2->cd(3);
        gPad->SetLogz(1);
        hNegMuPtvsPosMuPtvsRap->GetXaxis()->SetRangeUser(mRapLow[irap], mRapHi[irap]);
        hNegMuPtvsPosMuPt_Rap[irap] = (TH2D *)hNegMuPtvsPosMuPtvsRap->Project3D(Form("hNegMuPtvsPosMuPt_RapBin%d_zy", irap));
        hNegMuPtvsPosMuPt_Rap[irap]->SetTitle(Form("%1.1f < y < %1.1f", mRapLow[irap], mRapHi[irap]));
        hNegMuPtvsPosMuPt_Rap[irap]->GetYaxis()->SetTitleOffset(0.85);
        hNegMuPtvsPosMuPt_Rap[irap]->GetXaxis()->SetRangeUser(0, 4.0);
        hNegMuPtvsPosMuPt_Rap[irap]->GetYaxis()->SetRangeUser(0, 4.0);
        hNegMuPtvsPosMuPt_Rap[irap]->Draw("colz");
        drawLatex(0.39, 0.95, "Daughter #mu", mFont, 0.07, 1);

        c2->cd(4);
        gPad->SetLogz(1);
        hNegMuEtavsPosMuEtavsRap->GetXaxis()->SetRangeUser(mRapLow[irap], mRapHi[irap]);
        hNegMuEtavsPosMuEta_Rap[irap] = (TH2D *)hNegMuEtavsPosMuEtavsRap->Project3D(Form("hNegMuEtavsPosMuEta_RapBin%d_zy", irap));
        hNegMuEtavsPosMuEta_Rap[irap]->SetTitle(Form("%1.1f < y < %1.1f", mRapLow[irap], mRapHi[irap]));
        hNegMuEtavsPosMuEta_Rap[irap]->GetYaxis()->SetTitleOffset(0.85);
        if(mRapLow[irap]>0){
            hNegMuEtavsPosMuEta_Rap[irap]->GetXaxis()->SetRangeUser(1, 2.5);
            hNegMuEtavsPosMuEta_Rap[irap]->GetYaxis()->SetRangeUser(1, 2.5);
        }
        else{
            hNegMuEtavsPosMuEta_Rap[irap]->GetXaxis()->SetRangeUser(-2.5, -1);
            hNegMuEtavsPosMuEta_Rap[irap]->GetYaxis()->SetRangeUser(-2.5, -1);
        }
        hNegMuEtavsPosMuEta_Rap[irap]->Draw("colz");
        drawLatex(0.39, 0.95, "Daughter #mu", mFont, 0.07, 1);

        //c2->cd(4);
        //gPad->SetLogz(1);
        //hDeltaPtvsDeltaEtavsRap->GetXaxis()->SetRangeUser(mRapLow[irap], mRapHi[irap]);
        //hDeltaPtvsDeltaEta_Rap[irap] = (TH2D *)hDeltaPtvsDeltaEtavsRap->Project3D(Form("hDeltaPtvsDeltaEta_RapBin%d_zy", irap));
        //hDeltaPtvsDeltaEta_Rap[irap]->SetTitle(Form("%1.1f < y < %1.1f", mRapLow[irap], mRapHi[irap]));
        //hDeltaPtvsDeltaEta_Rap[irap]->GetYaxis()->SetTitleOffset(0.85);
        //hDeltaPtvsDeltaEta_Rap[irap]->GetXaxis()->SetRangeUser(-1.5, 1.5);
        //hDeltaPtvsDeltaEta_Rap[irap]->GetYaxis()->SetRangeUser(-0.5, 0.5);
        //hDeltaPtvsDeltaEta_Rap[irap]->Draw("colz");
        //drawLatex(0.36, 0.95, Form("%1.1f < y < %1.1f", mRapLow[irap], mRapHi[irap]), mFont, 0.07, 1);

        if(mStorePDF) c2->SaveAs(Form("%s/JpsiAnaAcc_RapBin%d.pdf", dirName.Data(), irap));
        else          c2->SaveAs(Form("%s/JpsiAnaAcc_RapBin%d.png", dirName.Data(), irap));
    }

    cout << "End of program !" << endl;
}
