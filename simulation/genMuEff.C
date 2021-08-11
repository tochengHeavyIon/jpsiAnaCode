#include "../common/headers.h"
#include "../common/function.C"
#include "../common/funUtil.h"

const Bool_t  mStorePDF = kTRUE;

const Double_t mTinyNum = 1.e-6;

Int_t    mGenMarker = 20;
Int_t    mGenColor = 1;
Int_t    mGenWidth = 2;
Double_t mGenSize = 0.8;

void genMuEff(TString filename="LowMassGammaGamma")
{
    gStyle->SetOptFit(1111);

    TFile* f = TFile::Open(Form("mcHistos/dimuonHistos.%s.root", filename.Data()));

    TH3D *hPosMuPhivsEtavsPt_Gen    = (TH3D *)f->Get("hPosMuPhivsEtavsPt_Gen");
    TH3D *hMthPosMuPhivsEtavsPt_Gen = (TH3D *)f->Get("hMthPosMuPhivsEtavsPt_Gen");
    TH3D *hMthPosMuPhivsEtavsPt     = (TH3D *)f->Get("hMthPosMuPhivsEtavsPt");
    TH3D *hTrigPosMuPhivsEtavsPt    = (TH3D *)f->Get("hTrigPosMuPhivsEtavsPt");
    TH3D *hNegMuPhivsEtavsPt_Gen    = (TH3D *)f->Get("hNegMuPhivsEtavsPt_Gen");
    TH3D *hMthNegMuPhivsEtavsPt_Gen = (TH3D *)f->Get("hMthNegMuPhivsEtavsPt_Gen");
    TH3D *hMthNegMuPhivsEtavsPt     = (TH3D *)f->Get("hMthNegMuPhivsEtavsPt");
    TH3D *hTrigNegMuPhivsEtavsPt    = (TH3D *)f->Get("hTrigNegMuPhivsEtavsPt");

    TH3D *hMthPosMuPhivsEtavsPtInpair  = (TH3D *)f->Get("hMthPosMuPhivsEtavsPtInpair");
    TH3D *hTrigPosMuPhivsEtavsPtInpair = (TH3D *)f->Get("hTrigPosMuPhivsEtavsPtInpair");
    TH3D *hMthNegMuPhivsEtavsPtInpair  = (TH3D *)f->Get("hMthNegMuPhivsEtavsPtInpair");
    TH3D *hTrigNegMuPhivsEtavsPtInpair = (TH3D *)f->Get("hTrigNegMuPhivsEtavsPtInpair");

    TString dir = Form("muEffPlots/%s", filename.Data());
    system(Form("mkdir -p %s", dir.Data()));
    system(Form("rm -rf %s/*", dir.Data()));

    TH2D *hPosMuEtavsPt_Gen = (TH2D *)hPosMuPhivsEtavsPt_Gen->Project3D("hPosMuEtavsPt_Gen_yx");
    TH2D *hPosMuPhivsPt_Gen = (TH2D *)hPosMuPhivsEtavsPt_Gen->Project3D("hPosMuPhivsPt_Gen_zx");
    TH2D *hMthPosMuEtavsPt_Gen = (TH2D *)hMthPosMuPhivsEtavsPt_Gen->Project3D("hMthPosMuEtavsPt_Gen_yx");
    TH2D *hMthPosMuPhivsPt_Gen = (TH2D *)hMthPosMuPhivsEtavsPt_Gen->Project3D("hMthPosMuPhivsPt_Gen_zx");

    TH2D *hPosMuPhivsEta_Gen    = (TH2D *)hPosMuPhivsEtavsPt_Gen->Project3D("hPosMuPhivsEta_Gen_zy");
    TH2D *hMthPosMuPhivsEta_Gen = (TH2D *)hMthPosMuPhivsEtavsPt_Gen->Project3D("hMthPosMuPhivsEta_Gen_zy");
    TH1D *hPosMuPhi_Gen         = (TH1D *)hPosMuPhivsEtavsPt_Gen->ProjectionZ("hPosMuPhi_Gen");
    TH1D *hMthPosMuPhi_Gen      = (TH1D *)hMthPosMuPhivsEtavsPt_Gen->ProjectionZ("hMthPosMuPhi_Gen");

    TH2D *hMthPosMuEtavsPt  = (TH2D *)hMthPosMuPhivsEtavsPt->Project3D("hMthPosMuEtavsPt_yx");
    TH2D *hMthPosMuPhivsPt  = (TH2D *)hMthPosMuPhivsEtavsPt->Project3D("hMthPosMuPhivsPt_zx");
    TH2D *hTrigPosMuEtavsPt = (TH2D *)hTrigPosMuPhivsEtavsPt->Project3D("hTrigPosMuEtavsPt_yx");
    TH2D *hTrigPosMuPhivsPt = (TH2D *)hTrigPosMuPhivsEtavsPt->Project3D("hTrigPosMuPhivsPt_zx");

    TH2D *hMthPosMuPhivsEta  = (TH2D *)hMthPosMuPhivsEtavsPt->Project3D("hMthPosMuPhivsEta_zy");
    TH2D *hTrigPosMuPhivsEta = (TH2D *)hTrigPosMuPhivsEtavsPt->Project3D("hTrigPosMuPhivsEta_zy");
    TH1D *hMthPosMuPhi       = (TH1D *)hMthPosMuPhivsEtavsPt->ProjectionZ("hMthPosMuPhi");
    TH1D *hTrigPosMuPhi      = (TH1D *)hTrigPosMuPhivsEtavsPt->ProjectionZ("hTrigPosMuPhi");

    TH2D *hNegMuEtavsPt_Gen = (TH2D *)hNegMuPhivsEtavsPt_Gen->Project3D("hNegMuEtavsPt_Gen_yx");
    TH2D *hNegMuPhivsPt_Gen = (TH2D *)hNegMuPhivsEtavsPt_Gen->Project3D("hNegMuPhivsPt_Gen_zx");
    TH2D *hMthNegMuEtavsPt_Gen = (TH2D *)hMthNegMuPhivsEtavsPt_Gen->Project3D("hMthNegMuEtavsPt_Gen_yx");
    TH2D *hMthNegMuPhivsPt_Gen = (TH2D *)hMthNegMuPhivsEtavsPt_Gen->Project3D("hMthNegMuPhivsPt_Gen_zx");

    TH2D *hNegMuPhivsEta_Gen    = (TH2D *)hNegMuPhivsEtavsPt_Gen->Project3D("hNegMuPhivsEta_Gen_zy");
    TH2D *hMthNegMuPhivsEta_Gen = (TH2D *)hMthNegMuPhivsEtavsPt_Gen->Project3D("hMthNegMuPhivsEta_Gen_zy");
    TH1D *hNegMuPhi_Gen         = (TH1D *)hNegMuPhivsEtavsPt_Gen->ProjectionZ("hNegMuPhi_Gen");
    TH1D *hMthNegMuPhi_Gen      = (TH1D *)hMthNegMuPhivsEtavsPt_Gen->ProjectionZ("hMthNegMuPhi_Gen");

    TH2D *hMthNegMuEtavsPt  = (TH2D *)hMthNegMuPhivsEtavsPt->Project3D("hMthNegMuEtavsPt_yx");
    TH2D *hMthNegMuPhivsPt  = (TH2D *)hMthNegMuPhivsEtavsPt->Project3D("hMthNegMuPhivsPt_zx");
    TH2D *hTrigNegMuEtavsPt = (TH2D *)hTrigNegMuPhivsEtavsPt->Project3D("hTrigNegMuEtavsPt_yx");
    TH2D *hTrigNegMuPhivsPt = (TH2D *)hTrigNegMuPhivsEtavsPt->Project3D("hTrigNegMuPhivsPt_zx");

    TH2D *hMthNegMuPhivsEta  = (TH2D *)hMthNegMuPhivsEtavsPt->Project3D("hMthNegMuPhivsEta_zy");
    TH2D *hTrigNegMuPhivsEta = (TH2D *)hTrigNegMuPhivsEtavsPt->Project3D("hTrigNegMuPhivsEta_zy");
    TH1D *hMthNegMuPhi       = (TH1D *)hMthNegMuPhivsEtavsPt->ProjectionZ("hMthNegMuPhi");
    TH1D *hTrigNegMuPhi      = (TH1D *)hTrigNegMuPhivsEtavsPt->ProjectionZ("hTrigNegMuPhi");

    Int_t mEtaBinLow = 1;
    Int_t mEtaBinMid = hMthPosMuPhivsEtavsPt->GetYaxis()->FindBin(0.);
    Int_t mEtaBinHi  = hMthPosMuPhivsEtavsPt->GetNbinsY();
    TH1D *hMthPosMuPhi_NegEta   = (TH1D *)hMthPosMuPhivsEtavsPt->ProjectionZ("hMthPosMuPhi_NegEta", 0, -1, mEtaBinLow, mEtaBinMid);
    TH1D *hMthPosMuPhi_PosEta   = (TH1D *)hMthPosMuPhivsEtavsPt->ProjectionZ("hMthPosMuPhi_PosEta", 0, -1, mEtaBinMid, mEtaBinHi);
    TH1D *hTrigPosMuPhi_NegEta  = (TH1D *)hTrigPosMuPhivsEtavsPt->ProjectionZ("hTrigPosMuPhi_NegEta", 0, -1, mEtaBinLow, mEtaBinMid);
    TH1D *hTrigPosMuPhi_PosEta  = (TH1D *)hTrigPosMuPhivsEtavsPt->ProjectionZ("hTrigPosMuPhi_PosEta", 0, -1, mEtaBinMid, mEtaBinHi);
    TH1D *hMthNegMuPhi_NegEta   = (TH1D *)hMthNegMuPhivsEtavsPt->ProjectionZ("hMthNegMuPhi_NegEta", 0, -1, mEtaBinLow, mEtaBinMid);
    TH1D *hMthNegMuPhi_PosEta   = (TH1D *)hMthNegMuPhivsEtavsPt->ProjectionZ("hMthNegMuPhi_PosEta", 0, -1, mEtaBinMid, mEtaBinHi);
    TH1D *hTrigNegMuPhi_NegEta  = (TH1D *)hTrigNegMuPhivsEtavsPt->ProjectionZ("hTrigNegMuPhi_NegEta", 0, -1, mEtaBinLow, mEtaBinMid);
    TH1D *hTrigNegMuPhi_PosEta  = (TH1D *)hTrigNegMuPhivsEtavsPt->ProjectionZ("hTrigNegMuPhi_PosEta", 0, -1, mEtaBinMid, mEtaBinHi);

    TH2D *hPosMu2DMthEff_EtavsPt = (TH2D *)hMthPosMuEtavsPt_Gen->Clone("hPosMu2DMthEff_EtavsPt");
    hPosMu2DMthEff_EtavsPt->Divide(hMthPosMuEtavsPt_Gen, hPosMuEtavsPt_Gen, 1, 1, "B");

    TH2D *hPosMu2DMthEff_PhivsPt = (TH2D *)hMthPosMuPhivsPt_Gen->Clone("hPosMu2DMthEff_PhivsPt");
    hPosMu2DMthEff_PhivsPt->Divide(hMthPosMuPhivsPt_Gen, hPosMuPhivsPt_Gen, 1, 1, "B");

    TH2D *hPosMu2DMthEff_PhivsEta = (TH2D *)hMthPosMuPhivsEta_Gen->Clone("hPosMu2DMthEff_PhivsEta");
    hPosMu2DMthEff_PhivsEta->Divide(hMthPosMuPhivsEta_Gen, hPosMuPhivsEta_Gen, 1, 1, "B");

    TH1D *hPosMuMthEffvsPhi = (TH1D *)hMthPosMuPhi_Gen->Clone("hPosMuMthEffvsPhi");
    hPosMuMthEffvsPhi->Divide(hMthPosMuPhi_Gen, hPosMuPhi_Gen, 1, 1, "B");

    TH2D *hNegMu2DMthEff_EtavsPt = (TH2D *)hMthNegMuEtavsPt_Gen->Clone("hNegMu2DMthEff_EtavsPt");
    hNegMu2DMthEff_EtavsPt->Divide(hMthNegMuEtavsPt_Gen, hNegMuEtavsPt_Gen, 1, 1, "B");

    TH2D *hNegMu2DMthEff_PhivsPt = (TH2D *)hMthNegMuPhivsPt_Gen->Clone("hNegMu2DMthEff_PhivsPt");
    hNegMu2DMthEff_PhivsPt->Divide(hMthNegMuPhivsPt_Gen, hNegMuPhivsPt_Gen, 1, 1, "B");

    TH2D *hNegMu2DMthEff_PhivsEta = (TH2D *)hMthNegMuPhivsEta_Gen->Clone("hNegMu2DMthEff_PhivsEta");
    hNegMu2DMthEff_PhivsEta->Divide(hMthNegMuPhivsEta_Gen, hNegMuPhivsEta_Gen, 1, 1, "B");

    TH1D *hNegMuMthEffvsPhi = (TH1D *)hMthNegMuPhi_Gen->Clone("hNegMuMthEffvsPhi");
    hNegMuMthEffvsPhi->Divide(hMthNegMuPhi_Gen, hNegMuPhi_Gen, 1, 1, "B");

    TH1D *hMuMthEffRatiovsPhi = (TH1D *)hPosMuMthEffvsPhi->Clone("hMuMthEffRatiovsPhi");
    hMuMthEffRatiovsPhi->Divide(hPosMuMthEffvsPhi, hNegMuMthEffvsPhi, 1, 1);

    TCanvas* c1 = new TCanvas("c1", "c1", 0, 0, 800, 600);
    c1->Divide(2, 2);

    Double_t xPos = 0.15;
    Double_t yPos = 0.84;

    c1->cd(1);
    gPad->SetLogy(0);
    gPad->SetLogz(0);
    hPosMu2DMthEff_EtavsPt->GetYaxis()->SetRangeUser(-2.4+mTinyNum, 2.4-mTinyNum);
    hPosMu2DMthEff_EtavsPt->Draw("colz");
    drawLatex(0.28, 0.945, "Reconstruted Efficiency (#mu^{+})", 22, 0.06, 4);
    drawLatex(xPos, yPos, "(a)", 22, 0.08, 1);
    c1->cd(2);
    gPad->SetLogy(0);
    gPad->SetLogz(0);
    hPosMu2DMthEff_PhivsPt->Draw("colz");
    drawLatex(0.28, 0.945, "Reconstruted Efficiency (#mu^{+})", 22, 0.06, 4);
    drawLatex(xPos, yPos, "(b)", 22, 0.08, 1);
    c1->cd(3);
    gPad->SetLogy(0);
    gPad->SetLogz(0);
    hNegMu2DMthEff_EtavsPt->GetYaxis()->SetRangeUser(-2.4+mTinyNum, 2.4-mTinyNum);
    hNegMu2DMthEff_EtavsPt->Draw("colz");
    drawLatex(0.28, 0.945, "Reconstruted Efficiency (#mu^{-})", 22, 0.06, 4);
    drawLatex(xPos, yPos, "(c)", 22, 0.08, 1);
    c1->cd(4);
    gPad->SetLogy(0);
    gPad->SetLogz(0);
    hNegMu2DMthEff_PhivsPt->Draw("colz");
    drawLatex(0.28, 0.945, "Reconstruted Efficiency (#mu^{-})", 22, 0.06, 4);
    drawLatex(xPos, yPos, "(d)", 22, 0.08, 1);
    c1->SaveAs(Form("%s/2DMthEff.png", dir.Data()));
    if(mStorePDF) c1->SaveAs(Form("%s/2DMthEff.pdf", dir.Data()));

    TLegend* leg = new TLegend(0.18, 0.18, 0.36, 0.4);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.1);

    c1->cd(1);
    gPad->SetLogy(0);
    gPad->SetLogz(0);
    hPosMu2DMthEff_PhivsEta->GetXaxis()->SetRangeUser(-2.4+mTinyNum, 2.4-mTinyNum);
    hPosMu2DMthEff_PhivsEta->Draw("colz");
    drawLatex(0.28, 0.945, "Reconstruted Efficiency (#mu^{+})", 22, 0.06, 4);
    drawLatex(0.48, 0.16, "(a)", 22, 0.08, 1);
    c1->cd(2);
    gPad->SetLogy(0);
    gPad->SetLogz(0);
    hNegMu2DMthEff_PhivsEta->GetXaxis()->SetRangeUser(-2.4+mTinyNum, 2.4-mTinyNum);
    hNegMu2DMthEff_PhivsEta->Draw("colz");
    drawLatex(0.28, 0.945, "Reconstruted Efficiency (#mu^{-})", 22, 0.06, 4);
    drawLatex(0.48, 0.16, "(b)", 22, 0.08, 1);
    c1->cd(3);
    gPad->SetLogy(0);
    gPad->SetLogz(0);
    setHisto(hPosMuMthEffvsPhi, 20, 0.5, 1, 1, 2);
    setHisto(hNegMuMthEffvsPhi, 24, 0.5, 4, 4, 2);
    hPosMuMthEffvsPhi->GetYaxis()->SetTitle("Reconstruted Efficiency");
    if(filename.EqualTo("LowMassGammaGamma")) hPosMuMthEffvsPhi->GetYaxis()->SetRangeUser(0.45, 0.65);
    else hPosMuMthEffvsPhi->GetYaxis()->SetRangeUser(0.6, 0.75);
    hPosMuMthEffvsPhi->Draw("histe");
    hNegMuMthEffvsPhi->Draw("histesame");
    leg->AddEntry(hPosMuMthEffvsPhi, "#mu^{+}", "l");
    leg->AddEntry(hNegMuMthEffvsPhi, "#mu^{-}", "l");
    leg->Draw("same");
    drawLatex(0.48, 0.16, "(c)", 22, 0.08, 1);
    c1->cd(4);
    gPad->SetLogy(0);
    gPad->SetLogz(0);
    setHisto(hMuMthEffRatiovsPhi, 20, 0.8, 1, 1, 2);
    hMuMthEffRatiovsPhi->GetYaxis()->SetRangeUser(0.8, 1.2);
    hMuMthEffRatiovsPhi->GetYaxis()->SetTitle("#mu^{+}/#mu^{-} Efficiency Ratio");
    hMuMthEffRatiovsPhi->Fit("pol0", "", "", -TMath::Pi(), TMath::Pi());
    hMuMthEffRatiovsPhi->Draw("psame");
    drawLatex(0.48, 0.16, "(d)", 22, 0.08, 1);
    c1->SaveAs(Form("%s/2DMthEff_PhivsEta.png", dir.Data()));
    if(mStorePDF) c1->SaveAs(Form("%s/2DMthEff_PhivsEta.pdf", dir.Data()));

    TH2D *hPosMu2DTrigEff_EtavsPt = (TH2D *)hTrigPosMuEtavsPt->Clone("hPosMu2DTrigEff_EtavsPt");
    hPosMu2DTrigEff_EtavsPt->Divide(hTrigPosMuEtavsPt, hMthPosMuEtavsPt, 1, 1, "B");

    TH2D *hPosMu2DTrigEff_PhivsPt = (TH2D *)hTrigPosMuPhivsPt->Clone("hPosMu2DTrigEff_PhivsPt");
    hPosMu2DTrigEff_PhivsPt->Divide(hTrigPosMuPhivsPt, hMthPosMuPhivsPt, 1, 1, "B");

    TH2D *hPosMu2DTrigEff_PhivsEta = (TH2D *)hTrigPosMuPhivsEta->Clone("hPosMu2DTrigEff_PhivsEta");
    hPosMu2DTrigEff_PhivsEta->Divide(hTrigPosMuPhivsEta, hMthPosMuPhivsEta, 1, 1, "B");

    TH1D *hPosMuTrigEffvsPhi = (TH1D *)hTrigPosMuPhi->Clone("hPosMuTrigEffvsPhi");
    hPosMuTrigEffvsPhi->Divide(hTrigPosMuPhi, hMthPosMuPhi, 1, 1, "B");

    TH1D *hPosMuTrigEffvsPhi_PosEta = (TH1D *)hTrigPosMuPhi_PosEta->Clone("hPosMuTrigEffvsPhi_PosEta");
    hPosMuTrigEffvsPhi_PosEta->Divide(hTrigPosMuPhi_PosEta, hMthPosMuPhi_PosEta, 1, 1, "B");

    TH1D *hPosMuTrigEffvsPhi_NegEta = (TH1D *)hTrigPosMuPhi_NegEta->Clone("hPosMuTrigEffvsPhi_NegEta");
    hPosMuTrigEffvsPhi_NegEta->Divide(hTrigPosMuPhi_NegEta, hMthPosMuPhi_NegEta, 1, 1, "B");

    TH2D *hNegMu2DTrigEff_EtavsPt = (TH2D *)hTrigNegMuEtavsPt->Clone("hNegMu2DTrigEff_EtavsPt");
    hNegMu2DTrigEff_EtavsPt->Divide(hTrigNegMuEtavsPt, hMthNegMuEtavsPt, 1, 1, "B");

    TH2D *hNegMu2DTrigEff_PhivsPt = (TH2D *)hTrigNegMuPhivsPt->Clone("hNegMu2DTrigEff_PhivsPt");
    hNegMu2DTrigEff_PhivsPt->Divide(hTrigNegMuPhivsPt, hMthNegMuPhivsPt, 1, 1, "B");

    TH2D *hNegMu2DTrigEff_PhivsEta = (TH2D *)hTrigNegMuPhivsEta->Clone("hNegMu2DTrigEff_PhivsEta");
    hNegMu2DTrigEff_PhivsEta->Divide(hTrigNegMuPhivsEta, hMthNegMuPhivsEta, 1, 1, "B");

    TH1D *hNegMuTrigEffvsPhi = (TH1D *)hTrigNegMuPhi->Clone("hNegMuTrigEffvsPhi");
    hNegMuTrigEffvsPhi->Divide(hTrigNegMuPhi, hMthNegMuPhi, 1, 1, "B");

    TH1D *hNegMuTrigEffvsPhi_PosEta = (TH1D *)hTrigNegMuPhi_PosEta->Clone("hNegMuTrigEffvsPhi_PosEta");
    hNegMuTrigEffvsPhi_PosEta->Divide(hTrigNegMuPhi_PosEta, hMthNegMuPhi_PosEta, 1, 1, "B");

    TH1D *hNegMuTrigEffvsPhi_NegEta = (TH1D *)hTrigNegMuPhi_NegEta->Clone("hNegMuTrigEffvsPhi_NegEta");
    hNegMuTrigEffvsPhi_NegEta->Divide(hTrigNegMuPhi_NegEta, hMthNegMuPhi_NegEta, 1, 1, "B");

    TH1D *hMuTrigEffRatiovsPhi_PosEta = (TH1D *)hPosMuTrigEffvsPhi_PosEta->Clone("hMuTrigEffRatiovsPhi_PosEta");
    hMuTrigEffRatiovsPhi_PosEta->Divide(hPosMuTrigEffvsPhi_PosEta, hNegMuTrigEffvsPhi_PosEta, 1, 1);

    TH1D *hMuTrigEffRatiovsPhi_NegEta = (TH1D *)hPosMuTrigEffvsPhi_NegEta->Clone("hMuTrigEffRatiovsPhi_NegEta");
    hMuTrigEffRatiovsPhi_NegEta->Divide(hPosMuTrigEffvsPhi_NegEta, hNegMuTrigEffvsPhi_NegEta, 1, 1);

    TH1D *hMuTrigEffRatiovsPhi = (TH1D *)hPosMuTrigEffvsPhi->Clone("hMuTrigEffRatiovsPhi");
    hMuTrigEffRatiovsPhi->Divide(hPosMuTrigEffvsPhi, hNegMuTrigEffvsPhi, 1, 1);

    c1->cd(1);
    gPad->SetLogy(0);
    gPad->SetLogz(0);
    hPosMu2DTrigEff_EtavsPt->GetYaxis()->SetRangeUser(-2.4+mTinyNum, 2.4-mTinyNum);
    hPosMu2DTrigEff_EtavsPt->Draw("colz");
    drawLatex(0.32, 0.945, "Trigger Efficiency (#mu^{+})", 22, 0.06, 4);
    drawLatex(xPos, yPos, "(a)", 22, 0.08, 1);
    c1->cd(2);
    gPad->SetLogy(0);
    gPad->SetLogz(0);
    hPosMu2DTrigEff_PhivsPt->Draw("colz");
    drawLatex(0.32, 0.945, "Trigger Efficiency (#mu^{+})", 22, 0.06, 4);
    drawLatex(xPos, yPos, "(b)", 22, 0.08, 1);
    c1->cd(3);
    gPad->SetLogy(0);
    gPad->SetLogz(0);
    hNegMu2DTrigEff_EtavsPt->GetYaxis()->SetRangeUser(-2.4+mTinyNum, 2.4-mTinyNum);
    hNegMu2DTrigEff_EtavsPt->Draw("colz");
    drawLatex(0.32, 0.945, "Trigger Efficiency (#mu^{-})", 22, 0.06, 4);
    drawLatex(xPos, yPos, "(c)", 22, 0.08, 1);
    c1->cd(4);
    gPad->SetLogy(0);
    gPad->SetLogz(0);
    hNegMu2DTrigEff_PhivsPt->Draw("colz");
    drawLatex(0.32, 0.945, "Trigger Efficiency (#mu^{-})", 22, 0.06, 4);
    drawLatex(xPos, yPos, "(d)", 22, 0.08, 1);
    c1->SaveAs(Form("%s/2DTrigEff.png", dir.Data()));
    if(mStorePDF) c1->SaveAs(Form("%s/2DTrigEff.pdf", dir.Data()));

    c1->cd(1);
    gPad->SetLogy(0);
    gPad->SetLogz(0);
    hPosMu2DTrigEff_PhivsEta->GetXaxis()->SetRangeUser(-2.4+mTinyNum, 2.4-mTinyNum);
    hPosMu2DTrigEff_PhivsEta->Draw("colz");
    drawLatex(0.32, 0.945, "Trigger Efficiency (#mu^{+})", 22, 0.06, 4);
    drawLatex(0.46, 0.16, "(a)", 22, 0.08, 1);
    c1->cd(2);
    gPad->SetLogy(0);
    gPad->SetLogz(0);
    hNegMu2DTrigEff_PhivsEta->GetXaxis()->SetRangeUser(-2.4+mTinyNum, 2.4-mTinyNum);
    hNegMu2DTrigEff_PhivsEta->Draw("colz");
    drawLatex(0.32, 0.945, "Trigger Efficiency (#mu^{-})", 22, 0.06, 4);
    drawLatex(0.46, 0.16, "(b)", 22, 0.08, 1);
    c1->cd(3);
    gPad->SetLogy(0);
    gPad->SetLogz(0);
    setHisto(hPosMuTrigEffvsPhi, 20, 0.5, 1, 1, 2);
    setHisto(hNegMuTrigEffvsPhi, 24, 0.5, 4, 4, 2);
    hPosMuTrigEffvsPhi->GetYaxis()->SetTitle("Trigger Efficiency");
    hPosMuTrigEffvsPhi->GetYaxis()->SetRangeUser(0.1, 0.35);
    hPosMuTrigEffvsPhi->Draw("histe");
    hNegMuTrigEffvsPhi->Draw("histesame");
    leg->Draw("same");
    drawLatex(0.46, 0.16, "(c)", 22, 0.08, 1);
    c1->cd(4);
    gPad->SetLogy(0);
    gPad->SetLogz(0);
    setHisto(hMuTrigEffRatiovsPhi, 20, 0.8, 1, 1, 2);
    hMuTrigEffRatiovsPhi->GetYaxis()->SetRangeUser(0.6, 1.4);
    hMuTrigEffRatiovsPhi->GetYaxis()->SetTitle("#mu^{+}/#mu^{-} Efficiency Ratio");
    hMuTrigEffRatiovsPhi->Fit("pol0", "", "", -TMath::Pi(), TMath::Pi());
    hMuTrigEffRatiovsPhi->Draw("psame");
    drawLatex(0.46, 0.16, "(d)", 22, 0.08, 1);
    c1->SaveAs(Form("%s/2DTrigEff_PhivsEta.png", dir.Data()));
    if(mStorePDF) c1->SaveAs(Form("%s/2DTrigEff_PhivsEta.pdf", dir.Data()));

    c1->cd(1);
    gPad->SetLogy(0);
    gPad->SetLogz(0);
    setHisto(hPosMuTrigEffvsPhi_PosEta, 20, 0.5, 1, 1, 2);
    setHisto(hNegMuTrigEffvsPhi_PosEta, 24, 0.5, 4, 4, 2);
    hPosMuTrigEffvsPhi_PosEta->GetYaxis()->SetTitle("Trigger Efficiency");
    hPosMuTrigEffvsPhi_PosEta->GetYaxis()->SetRangeUser(0.1, 0.35);
    hPosMuTrigEffvsPhi_PosEta->Draw("histe");
    hNegMuTrigEffvsPhi_PosEta->Draw("histesame");
    leg->Draw("same");
    drawLatex(xPos, yPos, "(a)", 22, 0.08, 1);
    drawLatex(0.45, 0.95, "|#eta| > 0", 22, 0.08, 4);
    c1->cd(2);
    gPad->SetLogy(0);
    gPad->SetLogz(0);
    setHisto(hMuTrigEffRatiovsPhi_PosEta, 20, 0.8, 1, 1, 2);
    hMuTrigEffRatiovsPhi_PosEta->GetYaxis()->SetRangeUser(0.6, 1.4);
    hMuTrigEffRatiovsPhi_PosEta->GetYaxis()->SetTitle("#mu^{+}/#mu^{-} Efficiency Ratio");
    hMuTrigEffRatiovsPhi_PosEta->Fit("pol0", "", "", -TMath::Pi(), TMath::Pi());
    hMuTrigEffRatiovsPhi_PosEta->Draw("psame");
    drawLatex(xPos, yPos, "(b)", 22, 0.08, 1);
    drawLatex(0.45, 0.95, "|#eta| > 0", 22, 0.08, 4);
    c1->cd(3);
    gPad->SetLogy(0);
    gPad->SetLogz(0);
    setHisto(hPosMuTrigEffvsPhi_NegEta, 20, 0.5, 1, 1, 2);
    setHisto(hNegMuTrigEffvsPhi_NegEta, 24, 0.5, 4, 4, 2);
    hPosMuTrigEffvsPhi_NegEta->GetYaxis()->SetTitle("Trigger Efficiency");
    hPosMuTrigEffvsPhi_NegEta->GetYaxis()->SetRangeUser(0.1, 0.35);
    hPosMuTrigEffvsPhi_NegEta->Draw("histe");
    hNegMuTrigEffvsPhi_NegEta->Draw("histesame");
    leg->Draw("same");
    drawLatex(xPos, yPos, "(c)", 22, 0.08, 1);
    drawLatex(0.45, 0.95, "|#eta| < 0", 22, 0.08, 4);
    c1->cd(4);
    gPad->SetLogy(0);
    gPad->SetLogz(0);
    setHisto(hMuTrigEffRatiovsPhi_NegEta, 20, 0.8, 1, 1, 2);
    hMuTrigEffRatiovsPhi_NegEta->GetYaxis()->SetRangeUser(0.6, 1.4);
    hMuTrigEffRatiovsPhi_NegEta->GetYaxis()->SetTitle("#mu^{+}/#mu^{-} Efficiency Ratio");
    hMuTrigEffRatiovsPhi_NegEta->Fit("pol0", "", "", -TMath::Pi(), TMath::Pi());
    hMuTrigEffRatiovsPhi_NegEta->Draw("psame");
    drawLatex(xPos, yPos, "(d)", 22, 0.08, 1);
    drawLatex(0.45, 0.95, "|#eta| < 0", 22, 0.08, 4);
    c1->SaveAs(Form("%s/TrigEff_PosEtavsNegEta.png", dir.Data()));
    if(mStorePDF) c1->SaveAs(Form("%s/TrigEff_PosEtavsNegEta.pdf", dir.Data()));

    TCanvas* c2 = new TCanvas("c2", "c2", 1000, 750);
    TPDF *ps = new TPDF(Form("%s/2DMthEffAndTrigEff.pdf", dir.Data()),111);
    ps->Off();

    Int_t nColumns = 3;
    Int_t nRaws    = 3;
    Int_t nPads = nColumns * nRaws;
    c2->Divide(nColumns, nRaws);

    mEtaBinLow = hPosMuPhivsEtavsPt_Gen->GetYaxis()->FindBin(-2.4 + mTinyNum);
    mEtaBinHi  = hPosMuPhivsEtavsPt_Gen->GetYaxis()->FindBin(2.4 - mTinyNum);
    const Int_t nEtaBins = mEtaBinHi - mEtaBinLow + 1;
    const Int_t nPhiBins = hPosMuPhivsEtavsPt_Gen->GetNbinsZ();

    TH1D *hPosMuPt_Gen[nEtaBins];
    TH1D *hMthPosMuPt_Gen[nEtaBins];
    TH1D *hPosMuMthEffvsPt[nEtaBins];
    TH1D *hMthPosMuPt[nEtaBins];
    TH1D *hTrigPosMuPt[nEtaBins];
    TH1D *hPosMuTrigEffvsPt[nEtaBins];
    TH1D *hMthPosMuPtInpair[nEtaBins];
    TH1D *hTrigPosMuPtInpair[nEtaBins];
    TH1D *hPosMuTrigEffvsPtInpair[nEtaBins];
    TH1D *hNegMuPt_Gen[nEtaBins];
    TH1D *hMthNegMuPt_Gen[nEtaBins];
    TH1D *hNegMuMthEffvsPt[nEtaBins];
    TH1D *hMthNegMuPt[nEtaBins];
    TH1D *hTrigNegMuPt[nEtaBins];
    TH1D *hNegMuTrigEffvsPt[nEtaBins];
    TH1D *hMthNegMuPtInpair[nEtaBins];
    TH1D *hTrigNegMuPtInpair[nEtaBins];
    TH1D *hNegMuTrigEffvsPtInpair[nEtaBins];
    for(Int_t ieta=0; ieta<nEtaBins; ieta++){
        hPosMuPt_Gen[ieta]    = (TH1D *)hPosMuPhivsEtavsPt_Gen->ProjectionX(Form("hPosMuPt_Gen_EtaBin%d", ieta), ieta+mEtaBinLow, ieta+mEtaBinLow, 1, nPhiBins);
        hMthPosMuPt_Gen[ieta] = (TH1D *)hMthPosMuPhivsEtavsPt_Gen->ProjectionX(Form("hMthPosMuPt_Gen_EtaBin%d", ieta), ieta+mEtaBinLow, ieta+mEtaBinLow, 1, nPhiBins);
        hPosMuMthEffvsPt[ieta] = (TH1D *)hMthPosMuPt_Gen[ieta]->Clone(Form("hPosMuMthEffvsPt_EtaBin%d", ieta));
        hPosMuMthEffvsPt[ieta]->GetYaxis()->SetTitle("Reconstructed Efficiency");
        hPosMuMthEffvsPt[ieta]->Divide(hMthPosMuPt_Gen[ieta], hPosMuPt_Gen[ieta], 1, 1, "B");

        hMthPosMuPt[ieta]  = (TH1D *)hMthPosMuPhivsEtavsPt->ProjectionX(Form("hMthPosMuPt_EtaBin%d", ieta), ieta+mEtaBinLow, ieta+mEtaBinLow, 1, nPhiBins);
        hTrigPosMuPt[ieta] = (TH1D *)hTrigPosMuPhivsEtavsPt->ProjectionX(Form("hTrigPosMuPt_EtaBin%d", ieta), ieta+mEtaBinLow, ieta+mEtaBinLow, 1, nPhiBins);
        hPosMuTrigEffvsPt[ieta] = (TH1D *)hTrigPosMuPt[ieta]->Clone(Form("hPosMuTrigEffvsPt_EtaBin%d", ieta));
        hPosMuTrigEffvsPt[ieta]->GetYaxis()->SetTitle("Trigger Efficiency");
        hPosMuTrigEffvsPt[ieta]->Divide(hTrigPosMuPt[ieta], hMthPosMuPt[ieta], 1, 1, "B");

        hMthPosMuPtInpair[ieta]  = (TH1D *)hMthPosMuPhivsEtavsPtInpair->ProjectionX(Form("hMthPosMuPtInpair_EtaBin%d", ieta), ieta+mEtaBinLow, ieta+mEtaBinLow, 1, nPhiBins);
        hTrigPosMuPtInpair[ieta] = (TH1D *)hTrigPosMuPhivsEtavsPtInpair->ProjectionX(Form("hTrigPosMuPtInpair_EtaBin%d", ieta), ieta+mEtaBinLow, ieta+mEtaBinLow, 1, nPhiBins);
        hPosMuTrigEffvsPtInpair[ieta] = (TH1D *)hTrigPosMuPtInpair[ieta]->Clone(Form("hPosMuTrigEffvsPtInpair_EtaBin%d", ieta));
        hPosMuTrigEffvsPtInpair[ieta]->GetYaxis()->SetTitle("Trigger Efficiency");
        hPosMuTrigEffvsPtInpair[ieta]->Divide(hTrigPosMuPtInpair[ieta], hMthPosMuPtInpair[ieta], 1, 1, "B");

        hNegMuPt_Gen[ieta]    = (TH1D *)hNegMuPhivsEtavsPt_Gen->ProjectionX(Form("hNegMuPt_Gen_EtaBin%d", ieta), ieta+mEtaBinLow, ieta+mEtaBinLow, 1, nPhiBins);
        hMthNegMuPt_Gen[ieta] = (TH1D *)hMthNegMuPhivsEtavsPt_Gen->ProjectionX(Form("hMthNegMuPt_Gen_EtaBin%d", ieta), ieta+mEtaBinLow, ieta+mEtaBinLow, 1, nPhiBins);
        hNegMuMthEffvsPt[ieta] = (TH1D *)hMthNegMuPt_Gen[ieta]->Clone(Form("hNegMuMthEffvsPt_EtaBin%d", ieta));
        hNegMuMthEffvsPt[ieta]->GetYaxis()->SetTitle("Reconstructed Efficiency");
        hNegMuMthEffvsPt[ieta]->Divide(hMthNegMuPt_Gen[ieta], hNegMuPt_Gen[ieta], 1, 1, "B");

        hMthNegMuPt[ieta]  = (TH1D *)hMthNegMuPhivsEtavsPt->ProjectionX(Form("hMthNegMuPt_EtaBin%d", ieta), ieta+mEtaBinLow, ieta+mEtaBinLow, 1, nPhiBins);
        hTrigNegMuPt[ieta] = (TH1D *)hTrigNegMuPhivsEtavsPt->ProjectionX(Form("hTrigNegMuPt_EtaBin%d", ieta), ieta+mEtaBinLow, ieta+mEtaBinLow, 1, nPhiBins);
        hNegMuTrigEffvsPt[ieta] = (TH1D *)hTrigNegMuPt[ieta]->Clone(Form("hNegMuTrigEffvsPt_EtaBin%d", ieta));
        hNegMuTrigEffvsPt[ieta]->GetYaxis()->SetTitle("Trigger Efficiency");
        hNegMuTrigEffvsPt[ieta]->Divide(hTrigNegMuPt[ieta], hMthNegMuPt[ieta], 1, 1, "B");

        hMthNegMuPtInpair[ieta]  = (TH1D *)hMthNegMuPhivsEtavsPtInpair->ProjectionX(Form("hMthNegMuPtInpair_EtaBin%d", ieta), ieta+mEtaBinLow, ieta+mEtaBinLow, 1, nPhiBins);
        hTrigNegMuPtInpair[ieta] = (TH1D *)hTrigNegMuPhivsEtavsPtInpair->ProjectionX(Form("hTrigNegMuPtInpair_EtaBin%d", ieta), ieta+mEtaBinLow, ieta+mEtaBinLow, 1, nPhiBins);
        hNegMuTrigEffvsPtInpair[ieta] = (TH1D *)hTrigNegMuPtInpair[ieta]->Clone(Form("hNegMuTrigEffvsPtInpair_EtaBin%d", ieta));
        hNegMuTrigEffvsPtInpair[ieta]->GetYaxis()->SetTitle("Trigger Efficiency");
        hNegMuTrigEffvsPtInpair[ieta]->Divide(hTrigNegMuPtInpair[ieta], hMthNegMuPtInpair[ieta], 1, 1, "B");
    }

    setLegend(leg, 0.6, 0.2, 0.9, 0.42, 0.10);

    Int_t ipad = 0;
    for(Int_t ieta=0; ieta<nEtaBins; ieta++){
        Double_t mEtaLow = hPosMuPhivsEtavsPt_Gen->GetYaxis()->GetBinLowEdge(ieta+mEtaBinLow);
        Double_t mEtaHi  = hPosMuPhivsEtavsPt_Gen->GetYaxis()->GetBinUpEdge(ieta+mEtaBinLow);

        Double_t etaCenter = hPosMuPhivsEtavsPt_Gen->GetYaxis()->GetBinCenter(ieta+mEtaBinLow);
        if(TMath::Abs(etaCenter)<1.3) continue;

        c2->cd(ipad%nPads + 1);
        setHisto(hPosMuMthEffvsPt[ieta], mGenMarker, mGenSize, mGenColor, mGenColor, mGenWidth);
        setHisto(hNegMuMthEffvsPt[ieta], mGenMarker+4, mGenSize, mGenColor+3, mGenColor+3, mGenWidth);
        hPosMuMthEffvsPt[ieta]->GetXaxis()->SetRangeUser(1, 4);
        hPosMuMthEffvsPt[ieta]->GetYaxis()->SetRangeUser(0, 1.1);
        hPosMuMthEffvsPt[ieta]->Draw("p");
        hNegMuMthEffvsPt[ieta]->Draw("psame");

        if(ipad==0){
            leg->AddEntry(hPosMuMthEffvsPt[ieta], "#mu^{+}", "p");
            leg->AddEntry(hNegMuMthEffvsPt[ieta], "#mu^{-}", "p");
        }
        if(ipad%nPads==0) leg->Draw("same");

        drawLatex(0.34, 0.95, Form("%1.2f < #eta < %1.2f", mEtaLow, mEtaHi), 22, 0.07, 2);

        if(ipad%nPads == nPads-1) pdfAction(c2, ps);

        ipad++;
    }
    if(ipad%nPads != 0) pdfAction(c2, ps);

    ipad = 0;
    for(Int_t ieta=0; ieta<nEtaBins; ieta++){
        Double_t mEtaLow = hTrigPosMuPhivsEtavsPt->GetYaxis()->GetBinLowEdge(ieta+mEtaBinLow);
        Double_t mEtaHi  = hTrigPosMuPhivsEtavsPt->GetYaxis()->GetBinUpEdge(ieta+mEtaBinLow);

        Double_t etaCenter = hPosMuPhivsEtavsPt_Gen->GetYaxis()->GetBinCenter(ieta+mEtaBinLow);
        if(TMath::Abs(etaCenter)<1.3) continue;

        c2->cd(ipad%nPads + 1);
        setHisto(hPosMuTrigEffvsPt[ieta], mGenMarker, mGenSize, mGenColor, mGenColor, mGenWidth);
        setHisto(hNegMuTrigEffvsPt[ieta], mGenMarker+4, mGenSize, mGenColor+3, mGenColor+3, mGenWidth);
        setHisto(hPosMuTrigEffvsPtInpair[ieta], mGenMarker, 0.1, mGenColor, mGenColor, 1);
        setHisto(hNegMuTrigEffvsPtInpair[ieta], mGenMarker+4, 0.1, mGenColor+3, mGenColor+3, 1);
        hPosMuTrigEffvsPtInpair[ieta]->GetXaxis()->SetRangeUser(1, 4);
        hPosMuTrigEffvsPtInpair[ieta]->GetYaxis()->SetRangeUser(0, 1.1);
        hPosMuTrigEffvsPtInpair[ieta]->Draw("hist");
        hNegMuTrigEffvsPtInpair[ieta]->Draw("histsamesame");
        hPosMuTrigEffvsPt[ieta]->Draw("psame");
        hNegMuTrigEffvsPt[ieta]->Draw("psame");

        //if(ipad==0){
        //    leg->AddEntry(hPosMuTrigEffvsPt[ieta], "#mu^{+}", "p");
        //    leg->AddEntry(hNegMuTrigEffvsPt[ieta], "#mu^{-}", "p");
        //}
        if(ipad%nPads==0) leg->Draw("same");

        drawLatex(0.34, 0.95, Form("%1.2f < #eta < %1.2f", mEtaLow, mEtaHi), 22, 0.07, 2);

        if(ipad%nPads == nPads-1) pdfAction(c2, ps);

        ipad++;
    }
    if(ipad%nPads != 0) pdfAction(c2, ps);

    ps->On();
    ps->Close();

    cout << "End of program !" << endl;
}
