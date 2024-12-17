#include "/Users/syang/Tools/Macro/headers.h"
#include "/Users/syang/Tools/Macro/function.C"
#include "/Users/syang/work/run2/upcDimuon/common/funUtil.h"

const Double_t mTinyNum = 1.e-6;

Int_t    mGenMarker = 20;
Int_t    mGenColor = 1;
Int_t    mGenWidth = 2;
Double_t mGenSize = 0.8;

void genEff(TString filename="GammaGamma")
{
    gStyle->SetOptFit(1111);

    TFile* f = TFile::Open(Form("dimuonHistos.%s.root", filename.Data()));

    TH3D *hPosMuPhivsEtavsPt_Gen    = (TH3D *)f->Get("hPosMuPhivsEtavsPt_Gen");
    TH3D *hMthPosMuPhivsEtavsPt_Gen = (TH3D *)f->Get("hMthPosMuPhivsEtavsPt_Gen");
    TH3D *hMthPosMuPhivsEtavsPt     = (TH3D *)f->Get("hMthPosMuPhivsEtavsPt");
    TH3D *hTrigPosMuPhivsEtavsPt    = (TH3D *)f->Get("hTrigPosMuPhivsEtavsPt");
    TH3D *hNegMuPhivsEtavsPt_Gen    = (TH3D *)f->Get("hNegMuPhivsEtavsPt_Gen");
    TH3D *hMthNegMuPhivsEtavsPt_Gen = (TH3D *)f->Get("hMthNegMuPhivsEtavsPt_Gen");
    TH3D *hMthNegMuPhivsEtavsPt     = (TH3D *)f->Get("hMthNegMuPhivsEtavsPt");
    TH3D *hTrigNegMuPhivsEtavsPt    = (TH3D *)f->Get("hTrigNegMuPhivsEtavsPt");

    TString dir = Form("effPlots_%s", filename.Data());
    system(Form("mkdir -p %s", dir.Data()));
    system(Form("rm -rf %s/*", dir.Data()));

    TH3D *hPosMu3DMthEff = (TH3D *)hMthPosMuPhivsEtavsPt_Gen->Clone("hPosMu3DMthEff");
    hPosMu3DMthEff->SetNameTitle("hPosMu3DMthEff", "p_{T}, #eta, #phi");
    hPosMu3DMthEff->Reset();
    hPosMu3DMthEff->Divide(hMthPosMuPhivsEtavsPt_Gen, hPosMuPhivsEtavsPt_Gen, 1, 1, "B");

    TH3D *hNegMu3DMthEff = (TH3D *)hMthNegMuPhivsEtavsPt_Gen->Clone("hNegMu3DMthEff");
    hNegMu3DMthEff->SetNameTitle("hNegMu3DMthEff", "p_{T}, #eta, #phi");
    hNegMu3DMthEff->Reset();
    hNegMu3DMthEff->Divide(hMthNegMuPhivsEtavsPt_Gen, hNegMuPhivsEtavsPt_Gen, 1, 1, "B");

    TH3D *hPosMu3DTrigEff = (TH3D *)hTrigPosMuPhivsEtavsPt->Clone("hPosMu3DTrigEff");
    hPosMu3DTrigEff->SetNameTitle("hPosMu3DTrigEff", "p_{T}, #eta, #phi");
    hPosMu3DTrigEff->Reset();
    hPosMu3DTrigEff->Divide(hTrigPosMuPhivsEtavsPt, hMthPosMuPhivsEtavsPt, 1, 1, "B");

    TH3D *hNegMu3DTrigEff = (TH3D *)hTrigNegMuPhivsEtavsPt->Clone("hNegMu3DTrigEff");
    hNegMu3DTrigEff->SetNameTitle("hNegMu3DTrigEff", "p_{T}, #eta, #phi");
    hNegMu3DTrigEff->Reset();
    hNegMu3DTrigEff->Divide(hTrigNegMuPhivsEtavsPt, hMthNegMuPhivsEtavsPt, 1, 1, "B");

    // to address Ota's comments
    Int_t nRebPhi = 60;
    TH3D *hPosMuPhivsEtavsPt_Gen_RebPhi    = (TH3D *)hPosMuPhivsEtavsPt_Gen->Clone("hPosMuPhivsEtavsPt_Gen_RebPhi");
    TH3D *hMthPosMuPhivsEtavsPt_Gen_RebPhi = (TH3D *)hMthPosMuPhivsEtavsPt_Gen->Clone("hMthPosMuPhivsEtavsPt_Gen_RebPhi");
    TH3D *hMthPosMuPhivsEtavsPt_RebPhi     = (TH3D *)hMthPosMuPhivsEtavsPt->Clone("hMthPosMuPhivsEtavsPt_RebPhi");
    TH3D *hTrigPosMuPhivsEtavsPt_RebPhi    = (TH3D *)hTrigPosMuPhivsEtavsPt->Clone("hTrigPosMuPhivsEtavsPt_RebPhi");
    TH3D *hNegMuPhivsEtavsPt_Gen_RebPhi    = (TH3D *)hNegMuPhivsEtavsPt_Gen->Clone("hNegMuPhivsEtavsPt_Gen_RebPhi");
    TH3D *hMthNegMuPhivsEtavsPt_Gen_RebPhi = (TH3D *)hMthNegMuPhivsEtavsPt_Gen->Clone("hMthNegMuPhivsEtavsPt_Gen_RebPhi");
    TH3D *hMthNegMuPhivsEtavsPt_RebPhi     = (TH3D *)hMthNegMuPhivsEtavsPt->Clone("hMthNegMuPhivsEtavsPt_RebPhi");
    TH3D *hTrigNegMuPhivsEtavsPt_RebPhi    = (TH3D *)hTrigNegMuPhivsEtavsPt->Clone("hTrigNegMuPhivsEtavsPt_RebPhi");
    hPosMuPhivsEtavsPt_Gen_RebPhi->RebinZ(nRebPhi);
    hMthPosMuPhivsEtavsPt_Gen_RebPhi->RebinZ(nRebPhi);
    hMthPosMuPhivsEtavsPt_RebPhi->RebinZ(nRebPhi);
    hTrigPosMuPhivsEtavsPt_RebPhi->RebinZ(nRebPhi);
    hNegMuPhivsEtavsPt_Gen_RebPhi->RebinZ(nRebPhi);
    hMthNegMuPhivsEtavsPt_Gen_RebPhi->RebinZ(nRebPhi);
    hMthNegMuPhivsEtavsPt_RebPhi->RebinZ(nRebPhi);
    hTrigNegMuPhivsEtavsPt_RebPhi->RebinZ(nRebPhi);

    TH3D *hPosMu3DMthEff_RebPhi = (TH3D *)hMthPosMuPhivsEtavsPt_Gen_RebPhi->Clone("hPosMu3DMthEff_RebPhi");
    hPosMu3DMthEff_RebPhi->SetNameTitle("hPosMu3DMthEff_RebPhi", "p_{T}, #eta, #phi");
    hPosMu3DMthEff_RebPhi->Reset();
    hPosMu3DMthEff_RebPhi->Divide(hMthPosMuPhivsEtavsPt_Gen_RebPhi, hPosMuPhivsEtavsPt_Gen_RebPhi, 1, 1, "B");

    TH3D *hNegMu3DMthEff_RebPhi = (TH3D *)hMthNegMuPhivsEtavsPt_Gen_RebPhi->Clone("hNegMu3DMthEff_RebPhi");
    hNegMu3DMthEff_RebPhi->SetNameTitle("hNegMu3DMthEff_RebPhi", "p_{T}, #eta, #phi");
    hNegMu3DMthEff_RebPhi->Reset();
    hNegMu3DMthEff_RebPhi->Divide(hMthNegMuPhivsEtavsPt_Gen_RebPhi, hNegMuPhivsEtavsPt_Gen_RebPhi, 1, 1, "B");

    TH3D *hPosMu3DTrigEff_RebPhi = (TH3D *)hTrigPosMuPhivsEtavsPt_RebPhi->Clone("hPosMu3DTrigEff_RebPhi");
    hPosMu3DTrigEff_RebPhi->SetNameTitle("hPosMu3DTrigEff_RebPhi", "p_{T}, #eta, #phi");
    hPosMu3DTrigEff_RebPhi->Reset();
    hPosMu3DTrigEff_RebPhi->Divide(hTrigPosMuPhivsEtavsPt_RebPhi, hMthPosMuPhivsEtavsPt_RebPhi, 1, 1, "B");

    TH3D *hNegMu3DTrigEff_RebPhi = (TH3D *)hTrigNegMuPhivsEtavsPt_RebPhi->Clone("hNegMu3DTrigEff_RebPhi");
    hNegMu3DTrigEff_RebPhi->SetNameTitle("hNegMu3DTrigEff_RebPhi", "p_{T}, #eta, #phi");
    hNegMu3DTrigEff_RebPhi->Reset();
    hNegMu3DTrigEff_RebPhi->Divide(hTrigNegMuPhivsEtavsPt_RebPhi, hMthNegMuPhivsEtavsPt_RebPhi, 1, 1, "B");

    TFile *fOut = new TFile(Form("%s/3DMthEffAndTrigEff.root", dir.Data()), "recreate");
    fOut->cd();

    hPosMu3DMthEff->Write();
    hNegMu3DMthEff->Write();
    hPosMu3DTrigEff->Write();
    hNegMu3DTrigEff->Write();

    hPosMu3DMthEff_RebPhi->Write();
    hNegMu3DMthEff_RebPhi->Write();
    hPosMu3DTrigEff_RebPhi->Write();
    hNegMu3DTrigEff_RebPhi->Write();

    fOut->Close();

    const Int_t mPtBinLow  = hPosMuPhivsEtavsPt_Gen->GetXaxis()->FindBin(3.5 + mTinyNum);
    const Int_t mPtBinHi   = hPosMuPhivsEtavsPt_Gen->GetNbinsX();
    const Int_t mEtaBinLow = hPosMuPhivsEtavsPt_Gen->GetYaxis()->FindBin(-2.4 + mTinyNum);
    const Int_t mEtaBinHi  = hPosMuPhivsEtavsPt_Gen->GetYaxis()->FindBin(2.4 - mTinyNum);

    TH2D *hPosMuEtavsPt_Gen = (TH2D *)hPosMuPhivsEtavsPt_Gen->Project3D("hPosMuEtavsPt_Gen_yx");
    TH2D *hPosMuPhivsPt_Gen = (TH2D *)hPosMuPhivsEtavsPt_Gen->Project3D("hPosMuPhivsPt_Gen_zx");
    TH2D *hMthPosMuEtavsPt_Gen = (TH2D *)hMthPosMuPhivsEtavsPt_Gen->Project3D("hMthPosMuEtavsPt_Gen_yx");
    TH2D *hMthPosMuPhivsPt_Gen = (TH2D *)hMthPosMuPhivsEtavsPt_Gen->Project3D("hMthPosMuPhivsPt_Gen_zx");

    hPosMuPhivsEtavsPt_Gen->GetXaxis()->SetRange(mPtBinLow, mPtBinHi);
    hMthPosMuPhivsEtavsPt_Gen->GetXaxis()->SetRange(mPtBinLow, mPtBinHi);
    TH2D *hPosMuPhivsEta_Gen    = (TH2D *)hPosMuPhivsEtavsPt_Gen->Project3D("hPosMuPhivsEta_Gen_zy");
    TH2D *hMthPosMuPhivsEta_Gen = (TH2D *)hMthPosMuPhivsEtavsPt_Gen->Project3D("hMthPosMuPhivsEta_Gen_zy");
    TH1D *hPosMuPhi_Gen         = (TH1D *)hPosMuPhivsEtavsPt_Gen->ProjectionZ("hPosMuPhi_Gen", mPtBinLow, mPtBinHi, mEtaBinLow, mEtaBinHi);
    TH1D *hMthPosMuPhi_Gen      = (TH1D *)hMthPosMuPhivsEtavsPt_Gen->ProjectionZ("hMthPosMuPhi_Gen", mPtBinLow, mPtBinHi, mEtaBinLow, mEtaBinHi);

    TH2D *hMthPosMuEtavsPt  = (TH2D *)hMthPosMuPhivsEtavsPt->Project3D("hMthPosMuEtavsPt_yx");
    TH2D *hMthPosMuPhivsPt  = (TH2D *)hMthPosMuPhivsEtavsPt->Project3D("hMthPosMuPhivsPt_zx");
    TH2D *hTrigPosMuEtavsPt = (TH2D *)hTrigPosMuPhivsEtavsPt->Project3D("hTrigPosMuEtavsPt_yx");
    TH2D *hTrigPosMuPhivsPt = (TH2D *)hTrigPosMuPhivsEtavsPt->Project3D("hTrigPosMuPhivsPt_zx");

    hMthPosMuPhivsEtavsPt->GetXaxis()->SetRange(mPtBinLow, mPtBinHi);
    hTrigPosMuPhivsEtavsPt->GetXaxis()->SetRange(mPtBinLow, mPtBinHi);
    TH2D *hMthPosMuPhivsEta  = (TH2D *)hMthPosMuPhivsEtavsPt->Project3D("hMthPosMuPhivsEta_zy");
    TH2D *hTrigPosMuPhivsEta = (TH2D *)hTrigPosMuPhivsEtavsPt->Project3D("hTrigPosMuPhivsEta_zy");
    TH1D *hMthPosMuPhi       = (TH1D *)hMthPosMuPhivsEtavsPt->ProjectionZ("hMthPosMuPhi", mPtBinLow, mPtBinHi, mEtaBinLow, mEtaBinHi);
    TH1D *hTrigPosMuPhi      = (TH1D *)hTrigPosMuPhivsEtavsPt->ProjectionZ("hTrigPosMuPhi", mPtBinLow, mPtBinHi, mEtaBinLow, mEtaBinHi);

    TH2D *hNegMuEtavsPt_Gen = (TH2D *)hNegMuPhivsEtavsPt_Gen->Project3D("hNegMuEtavsPt_Gen_yx");
    TH2D *hNegMuPhivsPt_Gen = (TH2D *)hNegMuPhivsEtavsPt_Gen->Project3D("hNegMuPhivsPt_Gen_zx");
    TH2D *hMthNegMuEtavsPt_Gen = (TH2D *)hMthNegMuPhivsEtavsPt_Gen->Project3D("hMthNegMuEtavsPt_Gen_yx");
    TH2D *hMthNegMuPhivsPt_Gen = (TH2D *)hMthNegMuPhivsEtavsPt_Gen->Project3D("hMthNegMuPhivsPt_Gen_zx");

    hNegMuPhivsEtavsPt_Gen->GetXaxis()->SetRange(mPtBinLow, mPtBinHi);
    hMthNegMuPhivsEtavsPt_Gen->GetXaxis()->SetRange(mPtBinLow, mPtBinHi);
    TH2D *hNegMuPhivsEta_Gen    = (TH2D *)hNegMuPhivsEtavsPt_Gen->Project3D("hNegMuPhivsEta_Gen_zy");
    TH2D *hMthNegMuPhivsEta_Gen = (TH2D *)hMthNegMuPhivsEtavsPt_Gen->Project3D("hMthNegMuPhivsEta_Gen_zy");
    TH1D *hNegMuPhi_Gen         = (TH1D *)hNegMuPhivsEtavsPt_Gen->ProjectionZ("hNegMuPhi_Gen", mPtBinLow, mPtBinHi, mEtaBinLow, mEtaBinHi);
    TH1D *hMthNegMuPhi_Gen      = (TH1D *)hMthNegMuPhivsEtavsPt_Gen->ProjectionZ("hMthNegMuPhi_Gen", mPtBinLow, mPtBinHi, mEtaBinLow, mEtaBinHi);

    TH2D *hMthNegMuEtavsPt  = (TH2D *)hMthNegMuPhivsEtavsPt->Project3D("hMthNegMuEtavsPt_yx");
    TH2D *hMthNegMuPhivsPt  = (TH2D *)hMthNegMuPhivsEtavsPt->Project3D("hMthNegMuPhivsPt_zx");
    TH2D *hTrigNegMuEtavsPt = (TH2D *)hTrigNegMuPhivsEtavsPt->Project3D("hTrigNegMuEtavsPt_yx");
    TH2D *hTrigNegMuPhivsPt = (TH2D *)hTrigNegMuPhivsEtavsPt->Project3D("hTrigNegMuPhivsPt_zx");

    hMthNegMuPhivsEtavsPt->GetXaxis()->SetRange(mPtBinLow, mPtBinHi);
    hTrigNegMuPhivsEtavsPt->GetXaxis()->SetRange(mPtBinLow, mPtBinHi);
    TH2D *hMthNegMuPhivsEta  = (TH2D *)hMthNegMuPhivsEtavsPt->Project3D("hMthNegMuPhivsEta_zy");
    TH2D *hTrigNegMuPhivsEta = (TH2D *)hTrigNegMuPhivsEtavsPt->Project3D("hTrigNegMuPhivsEta_zy");
    TH1D *hMthNegMuPhi       = (TH1D *)hMthNegMuPhivsEtavsPt->ProjectionZ("hMthNegMuPhi", mPtBinLow, mPtBinHi, mEtaBinLow, mEtaBinHi);
    TH1D *hTrigNegMuPhi      = (TH1D *)hTrigNegMuPhivsEtavsPt->ProjectionZ("hTrigNegMuPhi", mPtBinLow, mPtBinHi, mEtaBinLow, mEtaBinHi);

    TH2D *hPosMu2DMthEff_EtavsPt = (TH2D *)hMthPosMuEtavsPt_Gen->Clone("hPosMu2DMthEff_EtavsPt");
    hPosMu2DMthEff_EtavsPt->Reset();
    hPosMu2DMthEff_EtavsPt->Divide(hMthPosMuEtavsPt_Gen, hPosMuEtavsPt_Gen, 1, 1, "B");

    TH2D *hPosMu2DMthEff_PhivsPt = (TH2D *)hMthPosMuPhivsPt_Gen->Clone("hPosMu2DMthEff_PhivsPt");
    hPosMu2DMthEff_PhivsPt->Reset();
    hPosMu2DMthEff_PhivsPt->Divide(hMthPosMuPhivsPt_Gen, hPosMuPhivsPt_Gen, 1, 1, "B");

    TH2D *hPosMu2DMthEff_PhivsEta = (TH2D *)hMthPosMuPhivsEta_Gen->Clone("hPosMu2DMthEff_PhivsEta");
    hPosMu2DMthEff_PhivsEta->Reset();
    hPosMu2DMthEff_PhivsEta->Divide(hMthPosMuPhivsEta_Gen, hPosMuPhivsEta_Gen, 1, 1, "B");

    TH1D *hPosMuMthEffvsPhi = (TH1D *)hMthPosMuPhi_Gen->Clone("hPosMuMthEffvsPhi");
    hPosMuMthEffvsPhi->Reset();
    hPosMuMthEffvsPhi->Divide(hMthPosMuPhi_Gen, hPosMuPhi_Gen, 1, 1, "B");

    TH2D *hNegMu2DMthEff_EtavsPt = (TH2D *)hMthNegMuEtavsPt_Gen->Clone("hNegMu2DMthEff_EtavsPt");
    hNegMu2DMthEff_EtavsPt->Reset();
    hNegMu2DMthEff_EtavsPt->Divide(hMthNegMuEtavsPt_Gen, hNegMuEtavsPt_Gen, 1, 1, "B");

    TH2D *hNegMu2DMthEff_PhivsPt = (TH2D *)hMthNegMuPhivsPt_Gen->Clone("hNegMu2DMthEff_PhivsPt");
    hNegMu2DMthEff_PhivsPt->Reset();
    hNegMu2DMthEff_PhivsPt->Divide(hMthNegMuPhivsPt_Gen, hNegMuPhivsPt_Gen, 1, 1, "B");

    TH2D *hNegMu2DMthEff_PhivsEta = (TH2D *)hMthNegMuPhivsEta_Gen->Clone("hNegMu2DMthEff_PhivsEta");
    hNegMu2DMthEff_PhivsEta->Reset();
    hNegMu2DMthEff_PhivsEta->Divide(hMthNegMuPhivsEta_Gen, hNegMuPhivsEta_Gen, 1, 1, "B");

    TH1D *hNegMuMthEffvsPhi = (TH1D *)hMthNegMuPhi_Gen->Clone("hNegMuMthEffvsPhi");
    hNegMuMthEffvsPhi->Reset();
    hNegMuMthEffvsPhi->Divide(hMthNegMuPhi_Gen, hNegMuPhi_Gen, 1, 1, "B");

    TH1D *hMuMthEffRatiovsPhi = (TH1D *)hPosMuMthEffvsPhi->Clone("hMuMthEffRatiovsPhi");
    hMuMthEffRatiovsPhi->Reset();
    hMuMthEffRatiovsPhi->Divide(hPosMuMthEffvsPhi, hNegMuMthEffvsPhi, 1, 1);

    TCanvas* c1 = new TCanvas("c1", "c1", 0, 0, 800, 600);
    c1->Divide(2, 2);

    c1->cd(1);
    gPad->SetLogy(0);
    gPad->SetLogz(0);
    hPosMu2DMthEff_EtavsPt->GetYaxis()->SetRangeUser(-2.4+mTinyNum, 2.4-mTinyNum);
    hPosMu2DMthEff_EtavsPt->Draw("colz");
    drawLatex(0.28, 0.945, "Reconstruted Efficiency (#mu^{+})", 22, 0.06, 4);
    drawLatex(0.78, 0.18, "(a)", 22, 0.08, 1);
    c1->cd(2);
    gPad->SetLogy(0);
    gPad->SetLogz(0);
    hPosMu2DMthEff_PhivsPt->Draw("colz");
    drawLatex(0.28, 0.945, "Reconstruted Efficiency (#mu^{+})", 22, 0.06, 4);
    drawLatex(0.78, 0.18, "(b)", 22, 0.08, 1);
    c1->cd(3);
    gPad->SetLogy(0);
    gPad->SetLogz(0);
    hNegMu2DMthEff_EtavsPt->GetYaxis()->SetRangeUser(-2.4+mTinyNum, 2.4-mTinyNum);
    hNegMu2DMthEff_EtavsPt->Draw("colz");
    drawLatex(0.28, 0.945, "Reconstruted Efficiency (#mu^{-})", 22, 0.06, 4);
    drawLatex(0.78, 0.18, "(c)", 22, 0.08, 1);
    c1->cd(4);
    gPad->SetLogy(0);
    gPad->SetLogz(0);
    hNegMu2DMthEff_PhivsPt->Draw("colz");
    drawLatex(0.28, 0.945, "Reconstruted Efficiency (#mu^{-})", 22, 0.06, 4);
    drawLatex(0.78, 0.18, "(d)", 22, 0.08, 1);
    c1->SaveAs(Form("%s/2DMthEff.pdf", dir.Data()));
    c1->SaveAs(Form("%s/2DMthEff.png", dir.Data()));

    TLegend* leg = new TLegend(0.18, 0.18, 0.36, 0.4);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.1);

    c1->cd(1);
    gPad->SetLogy(0);
    gPad->SetLogz(0);
    hPosMu2DMthEff_PhivsEta->GetXaxis()->SetRangeUser(-2.4+mTinyNum, 2.4-mTinyNum);
    hPosMu2DMthEff_PhivsEta->Draw("colz");
    drawLatex(0.28, 0.945, "Reconstruted Efficiency (#mu^{+})", 22, 0.06, 4);
    drawLatex(0.78, 0.18, "(a)", 22, 0.08, 1);
    c1->cd(2);
    gPad->SetLogy(0);
    gPad->SetLogz(0);
    hNegMu2DMthEff_PhivsEta->GetXaxis()->SetRangeUser(-2.4+mTinyNum, 2.4-mTinyNum);
    hNegMu2DMthEff_PhivsEta->Draw("colz");
    drawLatex(0.28, 0.945, "Reconstruted Efficiency (#mu^{-})", 22, 0.06, 4);
    drawLatex(0.78, 0.18, "(b)", 22, 0.08, 1);
    c1->cd(3);
    gPad->SetLogy(0);
    gPad->SetLogz(0);
    setHisto(hPosMuMthEffvsPhi, 20, 0.5, 1, 1, 2);
    setHisto(hNegMuMthEffvsPhi, 24, 0.5, 4, 4, 2);
    hPosMuMthEffvsPhi->GetYaxis()->SetTitle("Reconstruted Efficiency");
    hPosMuMthEffvsPhi->GetYaxis()->SetRangeUser(0.6, 1.05);
    hPosMuMthEffvsPhi->Draw("hist");
    hNegMuMthEffvsPhi->Draw("histsame");
    leg->AddEntry(hPosMuMthEffvsPhi, "#mu^{+}", "l");
    leg->AddEntry(hNegMuMthEffvsPhi, "#mu^{-}", "l");
    leg->Draw("same");
    drawLatex(0.52, 0.2, "p_{T}^{#mu} > 3.5 GeV", 22, 0.09, 1);
    c1->cd(4);
    gPad->SetLogy(0);
    gPad->SetLogz(0);
    setHisto(hMuMthEffRatiovsPhi, 20, 0.8, 1, 1, 2);
    hMuMthEffRatiovsPhi->GetYaxis()->SetRangeUser(0.9, 1.1);
    hMuMthEffRatiovsPhi->GetYaxis()->SetTitle("#mu^{+}/#mu^{-} Efficiency Ratio");
    hMuMthEffRatiovsPhi->Fit("pol0", "", "", -TMath::Pi(), TMath::Pi());
    hMuMthEffRatiovsPhi->Draw("psame");
    c1->SaveAs(Form("%s/2DMthEff_PhivsEta.pdf", dir.Data()));
    c1->SaveAs(Form("%s/2DMthEff_PhivsEta.png", dir.Data()));

    TH2D *hPosMu2DTrigEff_EtavsPt = (TH2D *)hTrigPosMuEtavsPt->Clone("hPosMu2DTrigEff_EtavsPt");
    hPosMu2DTrigEff_EtavsPt->Reset();
    hPosMu2DTrigEff_EtavsPt->Divide(hTrigPosMuEtavsPt, hMthPosMuEtavsPt, 1, 1, "B");

    TH2D *hPosMu2DTrigEff_PhivsPt = (TH2D *)hTrigPosMuPhivsPt->Clone("hPosMu2DTrigEff_PhivsPt");
    hPosMu2DTrigEff_PhivsPt->Reset();
    hPosMu2DTrigEff_PhivsPt->Divide(hTrigPosMuPhivsPt, hMthPosMuPhivsPt, 1, 1, "B");

    TH2D *hPosMu2DTrigEff_PhivsEta = (TH2D *)hTrigPosMuPhivsEta->Clone("hPosMu2DTrigEff_PhivsEta");
    hPosMu2DTrigEff_PhivsEta->Reset();
    hPosMu2DTrigEff_PhivsEta->Divide(hTrigPosMuPhivsEta, hMthPosMuPhivsEta, 1, 1, "B");

    TH1D *hPosMuTrigEffvsPhi = (TH1D *)hTrigPosMuPhi->Clone("hPosMuTrigEffvsPhi");
    hPosMuTrigEffvsPhi->Reset();
    hPosMuTrigEffvsPhi->Divide(hTrigPosMuPhi, hMthPosMuPhi, 1, 1, "B");

    TH2D *hNegMu2DTrigEff_EtavsPt = (TH2D *)hTrigNegMuEtavsPt->Clone("hNegMu2DTrigEff_EtavsPt");
    hNegMu2DTrigEff_EtavsPt->Reset();
    hNegMu2DTrigEff_EtavsPt->Divide(hTrigNegMuEtavsPt, hMthNegMuEtavsPt, 1, 1, "B");

    TH2D *hNegMu2DTrigEff_PhivsPt = (TH2D *)hTrigNegMuPhivsPt->Clone("hNegMu2DTrigEff_PhivsPt");
    hNegMu2DTrigEff_PhivsPt->Reset();
    hNegMu2DTrigEff_PhivsPt->Divide(hTrigNegMuPhivsPt, hMthNegMuPhivsPt, 1, 1, "B");

    TH2D *hNegMu2DTrigEff_PhivsEta = (TH2D *)hTrigNegMuPhivsEta->Clone("hNegMu2DTrigEff_PhivsEta");
    hNegMu2DTrigEff_PhivsEta->Reset();
    hNegMu2DTrigEff_PhivsEta->Divide(hTrigNegMuPhivsEta, hMthNegMuPhivsEta, 1, 1, "B");

    TH1D *hNegMuTrigEffvsPhi = (TH1D *)hTrigNegMuPhi->Clone("hNegMuTrigEffvsPhi");
    hNegMuTrigEffvsPhi->Reset();
    hNegMuTrigEffvsPhi->Divide(hTrigNegMuPhi, hMthNegMuPhi, 1, 1, "B");

    TH1D *hMuTrigEffRatiovsPhi = (TH1D *)hPosMuTrigEffvsPhi->Clone("hMuTrigEffRatiovsPhi");
    hMuTrigEffRatiovsPhi->Reset();
    hMuTrigEffRatiovsPhi->Divide(hPosMuTrigEffvsPhi, hNegMuTrigEffvsPhi, 1, 1);

    c1->cd(1);
    gPad->SetLogy(0);
    gPad->SetLogz(0);
    hPosMu2DTrigEff_EtavsPt->GetYaxis()->SetRangeUser(-2.4+mTinyNum, 2.4-mTinyNum);
    hPosMu2DTrigEff_EtavsPt->Draw("colz");
    drawLatex(0.32, 0.945, "Trigger Efficiency (#mu^{+})", 22, 0.06, 4);
    drawLatex(0.78, 0.18, "(a)", 22, 0.08, 1);
    c1->cd(2);
    gPad->SetLogy(0);
    gPad->SetLogz(0);
    hPosMu2DTrigEff_PhivsPt->Draw("colz");
    drawLatex(0.32, 0.945, "Trigger Efficiency (#mu^{+})", 22, 0.06, 4);
    drawLatex(0.78, 0.18, "(b)", 22, 0.08, 1);
    c1->cd(3);
    gPad->SetLogy(0);
    gPad->SetLogz(0);
    hNegMu2DTrigEff_EtavsPt->GetYaxis()->SetRangeUser(-2.4+mTinyNum, 2.4-mTinyNum);
    hNegMu2DTrigEff_EtavsPt->Draw("colz");
    drawLatex(0.32, 0.945, "Trigger Efficiency (#mu^{-})", 22, 0.06, 4);
    drawLatex(0.78, 0.18, "(c)", 22, 0.08, 1);
    c1->cd(4);
    gPad->SetLogy(0);
    gPad->SetLogz(0);
    hNegMu2DTrigEff_PhivsPt->Draw("colz");
    drawLatex(0.32, 0.945, "Trigger Efficiency (#mu^{-})", 22, 0.06, 4);
    drawLatex(0.78, 0.18, "(d)", 22, 0.08, 1);
    c1->SaveAs(Form("%s/2DTrigEff.pdf", dir.Data()));
    c1->SaveAs(Form("%s/2DTrigEff.png", dir.Data()));

    c1->cd(1);
    gPad->SetLogy(0);
    gPad->SetLogz(0);
    hPosMu2DTrigEff_PhivsEta->GetXaxis()->SetRangeUser(-2.4+mTinyNum, 2.4-mTinyNum);
    hPosMu2DTrigEff_PhivsEta->Draw("colz");
    drawLatex(0.32, 0.945, "Trigger Efficiency (#mu^{+})", 22, 0.06, 4);
    drawLatex(0.78, 0.18, "(a)", 22, 0.08, 1);
    c1->cd(2);
    gPad->SetLogy(0);
    gPad->SetLogz(0);
    hNegMu2DTrigEff_PhivsEta->GetXaxis()->SetRangeUser(-2.4+mTinyNum, 2.4-mTinyNum);
    hNegMu2DTrigEff_PhivsEta->Draw("colz");
    drawLatex(0.32, 0.945, "Trigger Efficiency (#mu^{-})", 22, 0.06, 4);
    drawLatex(0.78, 0.18, "(b)", 22, 0.08, 1);
    c1->cd(3);
    gPad->SetLogy(0);
    gPad->SetLogz(0);
    setHisto(hPosMuTrigEffvsPhi, 20, 0.5, 1, 1, 2);
    setHisto(hNegMuTrigEffvsPhi, 24, 0.5, 4, 4, 2);
    hPosMuTrigEffvsPhi->GetYaxis()->SetTitle("Trigger Efficiency");
    hPosMuTrigEffvsPhi->GetYaxis()->SetRangeUser(0.6, 1.05);
    hPosMuTrigEffvsPhi->Draw("hist");
    hNegMuTrigEffvsPhi->Draw("histsame");
    leg->Draw("same");
    drawLatex(0.52, 0.2, "p_{T}^{#mu} > 3.5 GeV", 22, 0.09, 1);
    c1->cd(4);
    gPad->SetLogy(0);
    gPad->SetLogz(0);
    setHisto(hMuTrigEffRatiovsPhi, 20, 0.8, 1, 1, 2);
    hMuTrigEffRatiovsPhi->GetYaxis()->SetRangeUser(0.9, 1.1);
    hMuTrigEffRatiovsPhi->GetYaxis()->SetTitle("#mu^{+}/#mu^{-} Efficiency Ratio");
    hMuTrigEffRatiovsPhi->Fit("pol0", "", "", -TMath::Pi(), TMath::Pi());
    hMuTrigEffRatiovsPhi->Draw("psame");
    c1->SaveAs(Form("%s/2DTrigEff_PhivsEta.pdf", dir.Data()));
    c1->SaveAs(Form("%s/2DTrigEff_PhivsEta.png", dir.Data()));

    TCanvas* c2 = new TCanvas("c2", "c2", 1000, 600);
    TPDF *ps = new TPDF(Form("%s/3DMthEffAndTrigEff.pdf", dir.Data()),111);
    ps->Off();

    Int_t nColumns = 5;
    Int_t nRaws    = 4;
    Int_t nPads = nColumns * nRaws;
    c2->Divide(nColumns, nRaws);

    hPosMuPhivsEtavsPt_Gen->GetXaxis()->UnZoom();
    hMthPosMuPhivsEtavsPt_Gen->GetXaxis()->UnZoom();
    hMthPosMuPhivsEtavsPt->GetXaxis()->UnZoom();
    hTrigPosMuPhivsEtavsPt->GetXaxis()->UnZoom();
    hNegMuPhivsEtavsPt_Gen->GetXaxis()->UnZoom();
    hMthNegMuPhivsEtavsPt_Gen->GetXaxis()->UnZoom();
    hMthNegMuPhivsEtavsPt->GetXaxis()->UnZoom();
    hTrigNegMuPhivsEtavsPt->GetXaxis()->UnZoom();

    const Int_t nEtaBins = mEtaBinHi - mEtaBinLow + 1;
    const Int_t nPhiBins = hPosMuPhivsEtavsPt_Gen->GetNbinsZ();

    TH1D *hPosMuPt_Gen[nEtaBins][nPhiBins];
    TH1D *hMthPosMuPt_Gen[nEtaBins][nPhiBins];
    TH1D *hPosMuMthEffvsPt[nEtaBins][nPhiBins];
    TH1D *hMthPosMuPt[nEtaBins][nPhiBins];
    TH1D *hTrigPosMuPt[nEtaBins][nPhiBins];
    TH1D *hPosMuTrigEffvsPt[nEtaBins][nPhiBins];
    TH1D *hNegMuPt_Gen[nEtaBins][nPhiBins];
    TH1D *hMthNegMuPt_Gen[nEtaBins][nPhiBins];
    TH1D *hNegMuMthEffvsPt[nEtaBins][nPhiBins];
    TH1D *hMthNegMuPt[nEtaBins][nPhiBins];
    TH1D *hTrigNegMuPt[nEtaBins][nPhiBins];
    TH1D *hNegMuTrigEffvsPt[nEtaBins][nPhiBins];
    for(Int_t ieta=0; ieta<nEtaBins; ieta++){
        for(Int_t iphi=0; iphi<nPhiBins; iphi++){
            hPosMuPt_Gen[ieta][iphi]    = (TH1D *)hPosMuPhivsEtavsPt_Gen->ProjectionX(Form("hPosMuPt_Gen_EtaBin%d_PhiBin%d", ieta, iphi), ieta+mEtaBinLow, ieta+mEtaBinLow, iphi+1, iphi+1);
            hMthPosMuPt_Gen[ieta][iphi] = (TH1D *)hMthPosMuPhivsEtavsPt_Gen->ProjectionX(Form("hMthPosMuPt_Gen_EtaBin%d_PhiBin%d", ieta, iphi), ieta+mEtaBinLow, ieta+mEtaBinLow, iphi+1, iphi+1);
            hPosMuMthEffvsPt[ieta][iphi] = (TH1D *)hMthPosMuPt_Gen[ieta][iphi]->Clone(Form("hPosMuMthEffvsPt_EtaBin%d_PhiBin%d", ieta, iphi));
            hPosMuMthEffvsPt[ieta][iphi]->GetYaxis()->SetTitle("Reconstructed Efficiency");
            hPosMuMthEffvsPt[ieta][iphi]->Reset();
            hPosMuMthEffvsPt[ieta][iphi]->Divide(hMthPosMuPt_Gen[ieta][iphi], hPosMuPt_Gen[ieta][iphi], 1, 1, "B");
            //hPosMuMthEffvsPt[ieta][iphi]->Divide(hMthPosMuPt_Gen[ieta][iphi], hPosMuPt_Gen[ieta][iphi], 1, 1);

            hMthPosMuPt[ieta][iphi]  = (TH1D *)hMthPosMuPhivsEtavsPt->ProjectionX(Form("hMthPosMuPt_EtaBin%d_PhiBin%d", ieta, iphi), ieta+mEtaBinLow, ieta+mEtaBinLow, iphi+1, iphi+1);
            hTrigPosMuPt[ieta][iphi] = (TH1D *)hTrigPosMuPhivsEtavsPt->ProjectionX(Form("hTrigPosMuPt_EtaBin%d_PhiBin%d", ieta, iphi), ieta+mEtaBinLow, ieta+mEtaBinLow, iphi+1, iphi+1);
            hPosMuTrigEffvsPt[ieta][iphi] = (TH1D *)hTrigPosMuPt[ieta][iphi]->Clone(Form("hPosMuTrigEffvsPt_EtaBin%d_PhiBin%d", ieta, iphi));
            hPosMuTrigEffvsPt[ieta][iphi]->GetYaxis()->SetTitle("Trigger Efficiency");
            hPosMuTrigEffvsPt[ieta][iphi]->Reset();
            hPosMuTrigEffvsPt[ieta][iphi]->Divide(hTrigPosMuPt[ieta][iphi], hMthPosMuPt[ieta][iphi], 1, 1, "B");

            hNegMuPt_Gen[ieta][iphi]    = (TH1D *)hNegMuPhivsEtavsPt_Gen->ProjectionX(Form("hNegMuPt_Gen_EtaBin%d_PhiBin%d", ieta, iphi), ieta+mEtaBinLow, ieta+mEtaBinLow, iphi+1, iphi+1);
            hMthNegMuPt_Gen[ieta][iphi] = (TH1D *)hMthNegMuPhivsEtavsPt_Gen->ProjectionX(Form("hMthNegMuPt_Gen_EtaBin%d_PhiBin%d", ieta, iphi), ieta+mEtaBinLow, ieta+mEtaBinLow, iphi+1, iphi+1);
            hNegMuMthEffvsPt[ieta][iphi] = (TH1D *)hMthNegMuPt_Gen[ieta][iphi]->Clone(Form("hNegMuMthEffvsPt_EtaBin%d_PhiBin%d", ieta, iphi));
            hNegMuMthEffvsPt[ieta][iphi]->GetYaxis()->SetTitle("Reconstructed Efficiency");
            hNegMuMthEffvsPt[ieta][iphi]->Reset();
            hNegMuMthEffvsPt[ieta][iphi]->Divide(hMthNegMuPt_Gen[ieta][iphi], hNegMuPt_Gen[ieta][iphi], 1, 1, "B");

            hMthNegMuPt[ieta][iphi]  = (TH1D *)hMthNegMuPhivsEtavsPt->ProjectionX(Form("hMthNegMuPt_EtaBin%d_PhiBin%d", ieta, iphi), ieta+mEtaBinLow, ieta+mEtaBinLow, iphi+1, iphi+1);
            hTrigNegMuPt[ieta][iphi] = (TH1D *)hTrigNegMuPhivsEtavsPt->ProjectionX(Form("hTrigNegMuPt_EtaBin%d_PhiBin%d", ieta, iphi), ieta+mEtaBinLow, ieta+mEtaBinLow, iphi+1, iphi+1);
            hNegMuTrigEffvsPt[ieta][iphi] = (TH1D *)hTrigNegMuPt[ieta][iphi]->Clone(Form("hNegMuTrigEffvsPt_EtaBin%d_PhiBin%d", ieta, iphi));
            hNegMuTrigEffvsPt[ieta][iphi]->GetYaxis()->SetTitle("Trigger Efficiency");
            hNegMuTrigEffvsPt[ieta][iphi]->Reset();
            hNegMuTrigEffvsPt[ieta][iphi]->Divide(hTrigNegMuPt[ieta][iphi], hMthNegMuPt[ieta][iphi], 1, 1, "B");
        }
    }

    setLegend(leg, 0.6, 0.2, 0.9, 0.42, 0.10);

    for(Int_t ieta=0; ieta<nEtaBins; ieta++){
        Double_t mEtaLow = hPosMuPhivsEtavsPt_Gen->GetYaxis()->GetBinLowEdge(ieta+mEtaBinLow);
        Double_t mEtaHi  = hPosMuPhivsEtavsPt_Gen->GetYaxis()->GetBinUpEdge(ieta+mEtaBinLow);

        for(Int_t iphi=0; iphi<nPhiBins; iphi++){
            c2->cd(iphi%nPads + 1);

            Double_t mPhiLow = hPosMuPhivsEtavsPt_Gen->GetZaxis()->GetBinLowEdge(iphi+1);
            Double_t mPhiHi  = hPosMuPhivsEtavsPt_Gen->GetZaxis()->GetBinUpEdge(iphi+1);

            setHisto(hPosMuMthEffvsPt[ieta][iphi], mGenMarker, mGenSize, mGenColor, mGenColor, mGenWidth);
            setHisto(hNegMuMthEffvsPt[ieta][iphi], mGenMarker+4, mGenSize, mGenColor+3, mGenColor+3, mGenWidth);

            hPosMuMthEffvsPt[ieta][iphi]->GetYaxis()->SetRangeUser(0, 1.1);
            hPosMuMthEffvsPt[ieta][iphi]->Draw("p");
            hNegMuMthEffvsPt[ieta][iphi]->Draw("psame");

            if(ieta+iphi==0){
                leg->AddEntry(hPosMuMthEffvsPt[ieta][iphi], "#mu^{+}", "p");
                leg->AddEntry(hNegMuMthEffvsPt[ieta][iphi], "#mu^{-}", "p");
            }

            if(iphi%nPads==0) leg->Draw("same");

            drawLatex(0.2, 0.95, Form("%1.1f<#eta<%1.1f,  %1.3f<#phi<%1.3f", mEtaLow, mEtaHi, mPhiLow, mPhiHi), 22, 0.07, 2);

            if(iphi%nPads == nPads-1) pdfAction(c2, ps);
        }
        if(nPhiBins%nPads != 0) pdfAction(c2, ps);
    }

    for(Int_t ieta=0; ieta<nEtaBins; ieta++){
        Double_t mEtaLow = hMthPosMuPhivsEtavsPt->GetYaxis()->GetBinLowEdge(ieta+mEtaBinLow);
        Double_t mEtaHi  = hMthPosMuPhivsEtavsPt->GetYaxis()->GetBinUpEdge(ieta+mEtaBinLow);

        for(Int_t iphi=0; iphi<nPhiBins; iphi++){
            c2->cd(iphi%nPads + 1);

            Double_t mPhiLow = hMthPosMuPhivsEtavsPt->GetZaxis()->GetBinLowEdge(iphi+1);
            Double_t mPhiHi  = hMthPosMuPhivsEtavsPt->GetZaxis()->GetBinUpEdge(iphi+1);

            setHisto(hPosMuTrigEffvsPt[ieta][iphi], mGenMarker, mGenSize, mGenColor, mGenColor, mGenWidth);
            setHisto(hNegMuTrigEffvsPt[ieta][iphi], mGenMarker+4, mGenSize, mGenColor+3, mGenColor+3, mGenWidth);

            hPosMuTrigEffvsPt[ieta][iphi]->GetYaxis()->SetRangeUser(0, 1.1);
            hPosMuTrigEffvsPt[ieta][iphi]->Draw("p");
            hNegMuTrigEffvsPt[ieta][iphi]->Draw("psame");

            if(ieta+iphi==0){
                leg->Clear();
                leg->AddEntry(hPosMuTrigEffvsPt[ieta][iphi], "#mu^{+}", "p");
                leg->AddEntry(hNegMuTrigEffvsPt[ieta][iphi], "#mu^{-}", "p");
            }

            if(iphi%nPads==0) leg->Draw("same");

            drawLatex(0.2, 0.95, Form("%1.1f<#eta<%1.1f,  %1.3f<#phi<%1.3f", mEtaLow, mEtaHi, mPhiLow, mPhiHi), 22, 0.07, 2);

            if(iphi%nPads == nPads-1) pdfAction(c2, ps);
        }
        if(nPhiBins%nPads != 0) pdfAction(c2, ps);
    }

    ps->On();
    ps->Close();

    cout << "End of program !" << endl;
}
