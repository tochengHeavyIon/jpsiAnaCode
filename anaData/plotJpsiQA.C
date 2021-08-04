#include "../common/headers.h"
#include "../common/function.C"
#include "../common/funUtil.h"

const Bool_t   mStorePDF = kFALSE;
const Int_t    mFont = 42;
const Double_t mTinyNum = 1.e-6;

void plotJpsiQA(Bool_t incHadron = kFALSE, TString hfVetoType="Default")
{
    gStyle->SetOptFit(1111);

    if(!hfVetoType.EqualTo("Default") && !hfVetoType.EqualTo("Tight") && !hfVetoType.EqualTo("Loose") && !hfVetoType.EqualTo("removeHF")){
        cout<<"Please input the correct hfVetoType string: 'Default' OR 'Tight' OR 'Loose' OR 'removeHF'"<<endl;
        return;
    }

    if(!init()) {
        cout<<"Initialization failed !"<<endl;
        return;
    }

    TString histoDir = "jpsiHistos";
    TString plotDir  = "jpsiQAPlots";

    TString fileName = Form("%s/rawSig", histoDir.Data());
    TString dirName  = Form("%s/rawSig", plotDir.Data());

     if(incHadron){
        fileName += ".incHadron";
        dirName  += "_incHadron";
    }

    if(hfVetoType.EqualTo("Tight")){
        fileName += ".tightHF";
        dirName  += "_tightHF";
    }
    else if(hfVetoType.EqualTo("Loose")){
        fileName += ".looseHF";
        dirName  += "_looseHF";
    }
    else if(hfVetoType.EqualTo("removeHF")){
        fileName += ".removeHF";
        dirName  += "_removeHF";
    }

    cout<<"fileName: "<<fileName<<endl;
    cout<<"dirName: "<<dirName<<endl;

    TFile *f = TFile::Open(Form("%s.root", fileName.Data()));

    system(Form("mkdir -p %s", dirName.Data()));
    system(Form("rm -rf %s/*", dirName.Data()));

    TCanvas* c1 = new TCanvas("c1", "c1", 0, 0, 800, 600);
    setPad(0.12, 0.08, 0.07, 0.13);

    Int_t nColumns = 2;
    Int_t nRaws = 2;
    Int_t nPads = nColumns * nRaws;
    TCanvas* c2 = new TCanvas("c2", "c2", 1200, 900);
    c2->Divide(nColumns, nRaws);
    for(Int_t ipad=0; ipad<nPads; ipad++){
        c2->cd(ipad+1);
        setPad(0.12, 0.08, 0.07, 0.13);
    }

    TCanvas* c3 = new TCanvas("c3", "c3", 1200, 600);
    c3->Divide(3, 2);

    TCanvas* c4 = new TCanvas("c4", "c4", 1200, 900);
    c4->Divide(nNeus, nNeus);

    TCanvas* c5 = new TCanvas("c5", "c5", 1600, 900);
    c5->Divide(4, 3);

    TLegend* leg = new TLegend(0.15, 0.72, 0.4, 0.9);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.05);

    TH1D* hnEvts          = (TH1D *)f->Get("hnEvts");

    TH3D* hVzvsVyvsVx     = (TH3D*)f->Get("hVzvsVyvsVx");
    TH3D* hVzvsVyvsVx_Sel = (TH3D*)f->Get("hVzvsVyvsVx_Sel"); // Vertex selection

    TH3D* hHFMinusvsHFPlusvsCen     = (TH3D*)f->Get("hHFMinusvsHFPlusvsCen");
    TH3D* hHFMinusvsHFPlusvsCen_Sel = (TH3D*)f->Get("hHFMinusvsHFPlusvsCen_Sel");

    TH3D* hNtrkHPvsNtrkofflinevsCen     = (TH3D *)f->Get("hNtrkHPvsNtrkofflinevsCen");
    TH3D* hNtrkHPvsNtrkofflinevsCen_Sel = (TH3D *)f->Get("hNtrkHPvsNtrkofflinevsCen_Sel");
    TH1D* hNtrkHP_2SoftMuons            = (TH1D *)f->Get("hNtrkHP_2SoftMuons");

    TH2D* hZDCMinusvsZDCPlus_LS            = (TH2D *)f->Get("hZDCMinusvsZDCPlus_LS");
    TH2D* hZDCMinusvsZDCPlus_Sel_LS        = (TH2D *)f->Get("hZDCMinusvsZDCPlus_Sel_LS");
    TH2D* hZDCMinusvsZDCPlus_Only2MuTrk_LS = (TH2D *)f->Get("hZDCMinusvsZDCPlus_Only2MuTrk_LS");

    TH3D *hPosRap_Mass_2_5GeV_MuPhivsEtavsPt = (TH3D *)f->Get("hPosRap_Mass_2_5GeV_MuPhivsEtavsPt");
    TH3D *hNegRap_Mass_2_5GeV_MuPhivsEtavsPt = (TH3D *)f->Get("hNegRap_Mass_2_5GeV_MuPhivsEtavsPt");

    TH3D *hMvsPtvsRap     = (TH3D *)f->Get("hMvsPtvsRap");
    TH3D *hMvsAsyPhivsRap = (TH3D *)f->Get("hMvsAsyPhivsRap");
    TH3D *hMvsPtvsRap_WS = (TH3D *)f->Get("hMvsPtvsRap_WS");

    TH3D *hMvsPtvsRap_NeuDir[nNeus][nNeus];
    TH3D *hMvsAsyPhivsRap_NeuDir[nNeus][nNeus];
    for (Int_t ip = 0; ip < nNeus; ip++) {
        for (Int_t im = 0; im < nNeus; im++) {
            hMvsPtvsRap_NeuDir[ip][im]     = (TH3D*)f->Get(Form("hMvsPtvsRap_NeuDir%dp%dm", ip, im));
            hMvsAsyPhivsRap_NeuDir[ip][im] = (TH3D*)f->Get(Form("hMvsAsyPhivsRap_NeuDir%dp%dm", ip, im));
        }
    }

    Int_t neuIdx=0;
    TString neuPlus, neuMinus;

    const Int_t nPts = nNeus*nNeus;
    TString neuMultDirName[nPts];
    for(Int_t ip=0; ip<nNeus; ip++){
        for(Int_t im=0; im<nNeus; im++){
            if(ip == nNeus-1) neuPlus = "X";
            else              neuPlus = Form("%d", ip);
            if(im == nNeus-1) neuMinus = "X";
            else              neuMinus = Form("%d", im);

            neuMultDirName[neuIdx] = Form("Plus:%sn, Minus:%sn", neuPlus.Data(), neuMinus.Data());

            neuIdx++;
        }
    }

    TString neuMultName[nNeuMults];
    neuIdx = 0;
    for(Int_t ip=0; ip<nNeus; ip++){
        for(Int_t im=ip; im<nNeus; im++){
            if(ip == nNeus-1) neuPlus = "X";
            else              neuPlus = Form("%d", ip);
            if(im == nNeus-1) neuMinus = "X";
            else              neuMinus = Form("%d", im);

            neuMultName[neuIdx] = Form("%sn%sn", neuPlus.Data(), neuMinus.Data());

            neuIdx++;
        }
    }

    c1->cd();
    hnEvts->SetLineWidth(2);
    hnEvts->GetXaxis()->SetLabelSize(0.045);
    hnEvts->GetYaxis()->SetRangeUser(0, 0.9e6);
    hnEvts->SetMarkerColor(2);
    hnEvts->SetMarkerSize(1.6);
    hnEvts->Draw("histtext");
    c1->SaveAs(Form("%s/nEvts.png", dirName.Data()));

    TH2D* hVyvsVx = (TH2D*)hVzvsVyvsVx->Project3D("hVyvsVx_yx");
    TH1D* hVz = (TH1D*)hVzvsVyvsVx->ProjectionZ("hVz");
    TH2D* hVyvsVx_Sel = (TH2D*)hVzvsVyvsVx_Sel->Project3D("hVyvsVx_Sel_yx");
    TH1D* hVz_Sel = (TH1D*)hVzvsVyvsVx_Sel->ProjectionZ("hVz_Sel");
    c2->cd(1);
    gPad->SetLogz(1);
    hVyvsVx->GetXaxis()->SetRangeUser(-0.5, 0.5);
    hVyvsVx->GetYaxis()->SetRangeUser(-0.5, 0.5);
    hVyvsVx->Draw("col");
    drawLatex(0.36, 0.95, "w/o Vertex Selection", mFont, 0.06, 1);
    c2->cd(2);
    hVz->GetXaxis()->SetRangeUser(-30, 30);
    hVz->Draw("hist");
    drawLatex(0.36, 0.95, "w/o Vertex Selection", mFont, 0.06, 1);
    c2->cd(3);
    gPad->SetLogz(1);
    hVyvsVx_Sel->GetXaxis()->SetRangeUser(-0.5, 0.5);
    hVyvsVx_Sel->GetYaxis()->SetRangeUser(-0.5, 0.5);
    hVyvsVx_Sel->Draw("col");
    drawLatex(0.36, 0.95, "w/ Vertex Selection", mFont, 0.06, 2);
    c2->cd(4);
    hVz_Sel->GetXaxis()->SetRangeUser(-30, 30);
    hVz_Sel->Draw("hist");
    drawLatex(0.36, 0.95, "w/ Vertex Selection", mFont, 0.06, 2);
    drawLatex(0.16, 0.2, Form("Event Rejection: %1.1f%%", 100 * (1 - hVz_Sel->GetEntries() * 1. / hVz->GetEntries())), mFont, 0.07, 4);
    c2->SaveAs(Form("%s/vertex.png", dirName.Data()));

    clearPad(c2, nPads);

    TH1D *hCen       = (TH1D *)hHFMinusvsHFPlusvsCen->ProjectionX("hCen");
    TH1D *hCen_Sel   = (TH1D *)hHFMinusvsHFPlusvsCen_Sel->ProjectionX("hCen_Sel");
    TH1D* hCen_Final = (TH1D *)f->Get("hCen_Final");

    c1->cd();
    gPad->SetLogy(1);
    setHisto(hCen, 20, 0.1, 1, 1, 2);
    setHisto(hCen_Sel, 20, 0.1, 2, 2, 2);
    setHisto(hCen_Final, 20, 0.1, 4, 4, 2);
    hCen->SetAxisRange(120,200,"X");
    hCen->GetYaxis()->SetTitle("Entries");
    hCen->Draw("histe");
    hCen_Sel->Draw("histesame");
    hCen_Final->Draw("histesame");
    leg->AddEntry(hCen, "!BeamHalo & validVtx", "l");
    leg->AddEntry(hCen_Sel, "+ max HF Eenergy < 7.6 (7.3) GeV", "l");
    leg->AddEntry(hCen_Final, "+ N_{trk}^{HP} == 2 & N_{#mu}^{soft} == 2", "l");
    leg->Draw("same");
    c1->SaveAs(Form("%s/centrality.png", dirName.Data()));

    TH2D *hNtrkHPvsNtrkoffline     = (TH2D *)hNtrkHPvsNtrkofflinevsCen->Project3D("hNtrkHPvsNtrkoffline_zy");
    TH2D *hNtrkHPvsNtrkoffline_Sel = (TH2D *)hNtrkHPvsNtrkofflinevsCen_Sel->Project3D("hNtrkHPvsNtrkoffline_Sel_zy");
    TH1D *hNtrkHP                  = (TH1D *)hNtrkHPvsNtrkofflinevsCen->ProjectionZ("hNtrkHP");
    TH1D *hNtrkHP_Sel              = (TH1D *)hNtrkHPvsNtrkofflinevsCen_Sel->ProjectionZ("hNtrkHP_Sel");

    setLegend(leg, 0.24, 0.72, 0.6, 0.90, 0.055);

    c2->cd(1);
    gPad->SetLogz(1);
    hNtrkHPvsNtrkoffline->GetXaxis()->SetRangeUser(0, 100);
    hNtrkHPvsNtrkoffline->GetYaxis()->SetRangeUser(0, 100);
    hNtrkHPvsNtrkoffline->Draw("col");
    drawLatex(0.36, 0.95, "w/ vertex selection", mFont, 0.06, 1);
    c2->cd(2);
    gPad->SetLogz(1);
    hNtrkHPvsNtrkoffline_Sel->GetXaxis()->SetRangeUser(0, 30);
    hNtrkHPvsNtrkoffline_Sel->GetYaxis()->SetRangeUser(0, 30);
    hNtrkHPvsNtrkoffline_Sel->Draw("col");
    drawLatex(0.26, 0.95, "w/ vertex & HFVeto selection", mFont, 0.06, 1);
    c2->cd(3);
    gPad->SetLogy(1);
    setHisto(hNtrkHP, 20, 0.1, 1, 1, 2);
    setHisto(hNtrkHP_Sel, 20, 0.1, 2, 2, 2);
    setHisto(hNtrkHP_2SoftMuons, 20, 0.1, 4, 4, 2);
    hNtrkHP->GetXaxis()->SetRangeUser(0, 50);
    hNtrkHP->GetYaxis()->SetTitle("Entries");
    hNtrkHP->SetMinimum(0.5);
    hNtrkHP->Draw("histe");
    hNtrkHP_Sel->Draw("histesame");
    hNtrkHP_2SoftMuons->Draw("histesame");
    leg->AddEntry(hNtrkHP, "!BeamHalo & validVtx", "l");
    leg->AddEntry(hNtrkHP_Sel, "+ max HF Eenergy < 7.6 (7.3) GeV", "l");
    leg->AddEntry(hNtrkHP_2SoftMuons, "+ N_{#mu}^{soft} == 2", "l");
    leg->Draw("same");
    c2->SaveAs(Form("%s/NtrkHPvsNtrkoffline.png", dirName.Data()));

    c1->cd();
    gPad->SetLogy(1);
    setHisto(hNtrkHP, 20, 0.1, 1, 1, 2);
    setHisto(hNtrkHP_Sel, 20, 0.1, 2, 2, 2);
    setHisto(hNtrkHP_2SoftMuons, 20, 0.1, 4, 4, 2);
    hNtrkHP->GetXaxis()->SetRangeUser(0, 50);
    hNtrkHP->GetYaxis()->SetTitle("Entries");
    hNtrkHP->SetMinimum(0.5);
    hNtrkHP->Draw("histe");
    hNtrkHP_Sel->Draw("histesame");
    hNtrkHP_2SoftMuons->Draw("histesame");
    leg->Draw("same");
    c1->SaveAs(Form("%s/NtrkHP.png", dirName.Data()));

    TH1D *hZDCPlus_LS  = (TH1D *)hZDCMinusvsZDCPlus_LS->ProjectionX("hZDCPlus_LS");
    TH1D *hZDCMinus_LS = (TH1D *)hZDCMinusvsZDCPlus_LS->ProjectionY("hZDCMinus_LS");
    TH1D *hZDCPlus_Sel_LS  = (TH1D *)hZDCMinusvsZDCPlus_Sel_LS->ProjectionX("hZDCPlus_Sel_LS");
    TH1D *hZDCMinus_Sel_LS = (TH1D *)hZDCMinusvsZDCPlus_Sel_LS->ProjectionY("hZDCMinus_Sel_LS");
    TH1D *hZDCPlus_Only2MuTrk_LS  = (TH1D *)hZDCMinusvsZDCPlus_Only2MuTrk_LS->ProjectionX("hZDCPlus_Only2MuTrk_LS");
    TH1D *hZDCMinus_Only2MuTrk_LS = (TH1D *)hZDCMinusvsZDCPlus_Only2MuTrk_LS->ProjectionY("hZDCMinus_Only2MuTrk_LS");

    c1->cd();
    gPad->SetLogy(1);
    setHisto(hZDCPlus_LS, 20, 0.1, 1, 1, 2);
    setHisto(hZDCPlus_Sel_LS, 20, 0.1, 2, 2, 2);
    setHisto(hZDCPlus_Only2MuTrk_LS, 20, 0.1, 4, 4, 2);
    hZDCPlus_LS->SetMinimum(1);
    hZDCPlus_LS->GetXaxis()->SetRangeUser(0, 8e5);
    hZDCPlus_LS->Draw("hist");
    hZDCPlus_Sel_LS->Draw("histsame");
    hZDCPlus_Only2MuTrk_LS->Draw("histsame");
    setLegend(leg, 0.36, 0.68, 0.84, 0.88, 0.06);
    leg->AddEntry(hZDCPlus_LS, "w/ vtx selection", "l");
    leg->AddEntry(hZDCPlus_Sel_LS, "+ HF Rejection", "l");
    leg->AddEntry(hZDCPlus_Only2MuTrk_LS, "+ N_{trk}^{HP} == 2 && N_{#mu} == 2", "l");
    leg->Draw("same");
    c1->SaveAs(Form("%s/ZDCPlus_Compare.png", dirName.Data()));

    c1->cd();
    gPad->SetLogy(1);
    setHisto(hZDCMinus_LS, 20, 0.1, 1, 1, 2);
    setHisto(hZDCMinus_Sel_LS, 20, 0.1, 2, 2, 2);
    setHisto(hZDCMinus_Only2MuTrk_LS, 20, 0.1, 4, 4, 2);
    hZDCMinus_LS->SetMinimum(1);
    hZDCMinus_LS->GetXaxis()->SetRangeUser(0, 1e6);
    hZDCMinus_LS->GetYaxis()->SetTitle("Entries");
    hZDCMinus_LS->Draw("hist");
    hZDCMinus_Sel_LS->Draw("histsame");
    hZDCMinus_Only2MuTrk_LS->Draw("histsame");
    leg->Draw("same");
    c1->SaveAs(Form("%s/ZDCMinus_Compare.png", dirName.Data()));

    TH2D *hMvsPt_CS = (TH2D*)hMvsPtvsRap->Project3D("hMvsPt_CS_zy");
    TH1D *hMass_CS  = (TH1D*)hMvsPtvsRap->ProjectionZ("hMass_CS", 0, -1, 0, -1);

    TH2D *hMvsPt_WS = (TH2D*)hMvsPtvsRap_WS->Project3D("hMvsPt_WS_zy");
    TH1D *hMass_WS  = (TH1D*)hMvsPtvsRap_WS->ProjectionZ("hMass_WS", 0, -1, 0, -1);

    Int_t ptBinLow = 1, ptBinHi = hMvsPtvsRap->GetNbinsY();
    Int_t massBinLow = 1, massBinHi = hMvsPtvsRap->GetNbinsZ();
    TH1D *hCS_Rap  = (TH1D *)hMvsPtvsRap->ProjectionX("hCS_Rap", ptBinLow, ptBinHi, massBinLow, massBinHi);

    c2->cd(1);
    gPad->SetLogy(0);
    gPad->SetLogz(1);
    hMvsPt_CS->Draw("col");
    drawLatex(0.4, 0.95, "Correct-Sign", mFont, 0.06, 1);
    c2->cd(2);
    gPad->SetLogy(0);
    gPad->SetLogz(1);
    //hMvsPt_WS->SetFillColor(2);
    hMvsPt_WS->Draw("box");
    drawLatex(0.4, 0.95, "Wrong-Sign", mFont, 0.06, 1);
    c2->cd(3);
    gPad->SetLogy(1);
    setLegend(leg, 0.5, 0.76, 0.84, 0.88, 0.06);
    leg->AddEntry(hMass_CS, "Correct-Sign", "p");
    leg->AddEntry(hMass_WS, "Wrong-Sign", "l");
    setHisto(hMass_CS, 20, 0.5, 1, 1, 2);
    setHisto(hMass_WS, 24, 0.5, 2, 2, 2);
    hMass_WS->SetLineStyle(2);
    hMass_CS->RebinX(2);
    hMass_WS->RebinX(2);
    hMass_CS->GetYaxis()->SetRangeUser(0.5, 5e4);
    hMass_CS->GetYaxis()->SetTitle(Form("Entries / (%1.2f GeV/c^{2})", hMass_CS->GetXaxis()->GetBinWidth(1)));
    hMass_CS->Draw("p");
    hMass_WS->Draw("histsame");
    leg->Draw("same");
    c2->cd(4);
    gPad->SetLogy(1);
    setHisto(hCS_Rap, 20, 0.5, 1, 1, 2);
    hCS_Rap->GetYaxis()->SetRangeUser(0.5, 5e4);
    hCS_Rap->Draw("p");
    drawLatex(0.4, 0.95, "Correct-Sign", mFont, 0.06, 1);
    drawLatex(0.38, 0.85, "2 < M_{#mu#mu} < 5 GeV", mFont, 0.06, 1);
    c2->SaveAs(Form("%s/MvsPt.png", dirName.Data()));

    TH1D *hRawMass_Rap[nRapBins];
    TH1D *hRawPt_Rap[nRapBins];
    for(Int_t irap=0; irap<nRapBins; irap++){
        Int_t ptBinLow   = 1;
        Int_t ptBinHi    = hMvsPtvsRap->GetNbinsY();
        Int_t massBinLow = hMvsPtvsRap->GetZaxis()->FindBin(mJpsiMassLow + mTinyNum);
        Int_t massBinHi  = hMvsPtvsRap->GetZaxis()->FindBin(mJpsiMassHi - mTinyNum);
        Int_t rapBinLow  = hMvsPtvsRap->GetXaxis()->FindBin(mRapLow[irap] + mTinyNum);
        Int_t rapBinHi   = hMvsPtvsRap->GetXaxis()->FindBin(mRapHi[irap] - mTinyNum);
        hRawPt_Rap[irap]   = (TH1D *)hMvsPtvsRap->ProjectionY(Form("hRawPt_RapBin%d", irap), rapBinLow, rapBinHi, massBinLow, massBinHi);
        hRawMass_Rap[irap] = (TH1D *)hMvsPtvsRap->ProjectionZ(Form("hRawMass_RapBin%d", irap), rapBinLow, rapBinHi, ptBinLow, ptBinHi);

        c5->cd(irap+1);
        gPad->SetLogy(1);
        setHisto(hRawMass_Rap[irap], 20, 0.6, 1, 1, 1);
        hRawMass_Rap[irap]->RebinX(2);
        //hRawMass_Rap[irap]->SetMinimum(0);
        hRawMass_Rap[irap]->GetYaxis()->SetRangeUser(0.5, hRawMass_Rap[irap]->GetMaximum()*5);
        hRawMass_Rap[irap]->GetXaxis()->SetTitleOffset(0.9);
        hRawMass_Rap[irap]->GetYaxis()->SetTitle("Entries");
        hRawMass_Rap[irap]->Draw("p");
        drawLatex(0.39, 0.95, Form("%1.1f < y < %1.1f", mRapLow[irap], mRapHi[irap]), mFont, 0.06, 4);
    }
    c5->SaveAs(Form("%s/rawMassvsRap.png", dirName.Data()));
    if(mStorePDF) c5->SaveAs(Form("%s/rawMassvsRap.pdf", dirName.Data()));

    for(Int_t irap=0; irap<nRapBins; irap++){
        c5->cd(irap+1);
        gPad->SetLogy(1);
        setHisto(hRawPt_Rap[irap], 20, 0.6, 1, 1, 1);
        hRawPt_Rap[irap]->RebinX(2);
        hRawPt_Rap[irap]->GetXaxis()->SetRangeUser(0, 3);
        hRawPt_Rap[irap]->GetYaxis()->SetRangeUser(0.5, hRawPt_Rap[irap]->GetMaximum()*5);
        hRawPt_Rap[irap]->GetXaxis()->SetTitleOffset(0.9);
        hRawPt_Rap[irap]->GetYaxis()->SetTitle("Entries");
        hRawPt_Rap[irap]->Draw("p");
        drawLatex(0.36, 0.95, Form("%1.1f < y < %1.1f", mRapLow[irap], mRapHi[irap]), mFont, 0.07, 4);
        if(irap==0) drawLatex(0.42, 0.82, Form("%1.2f < M_{#mu#mu} < %1.2f GeV", mJpsiMassLow, mJpsiMassHi), mFont, 0.065, 1);
    }
    c5->SaveAs(Form("%s/rawPtvsRap.png", dirName.Data()));
    if(mStorePDF) c5->SaveAs(Form("%s/rawPtvsRap.pdf", dirName.Data()));

    TH1D *hMass_NeuDir[nNeus][nNeus][nDiffRapBins];
    TH1D *hPt_NeuDir[nNeus][nNeus][nDiffRapBins];
    TH1D *hLowMassBandPt_NeuDir[nNeus][nNeus][nDiffRapBins];
    TH1D *hJpsiPt_NeuDir[nNeus][nNeus][nDiffRapBins];
    TH1D *hHiMassBandPt_NeuDir[nNeus][nNeus][nDiffRapBins];

    TH1D *hMass_Neu[nNeus][nNeus][nDiffRapBins];
    TH1D *hPt_Neu[nNeus][nNeus][nDiffRapBins];
    for(Int_t ip=0; ip<nNeus; ip++){
        for(Int_t im=0; im<nNeus; im++){
            for(Int_t irap=0; irap<nDiffRapBins; irap++){
                Int_t rapBinLow  = hMvsPtvsRap_NeuDir[ip][im]->GetXaxis()->FindBin(mDiffRapLow[irap] + mTinyNum);
                Int_t rapBinHi   = hMvsPtvsRap_NeuDir[ip][im]->GetXaxis()->FindBin(mDiffRapHi[irap] - mTinyNum);
                Int_t ptBinLow   = 1;
                Int_t ptBinHi    = hMvsPtvsRap_NeuDir[ip][im]->GetNbinsY();
                Int_t massBinLow = 1;
                Int_t massBinHi  = hMvsPtvsRap_NeuDir[ip][im]->GetNbinsZ();
                hMass_NeuDir[ip][im][irap] = (TH1D *)hMvsPtvsRap_NeuDir[ip][im]->ProjectionZ(Form("hMass_NeuDir%dp%dm_RapBin%d", ip, im, irap), rapBinLow, rapBinHi, ptBinLow, ptBinHi);
                hPt_NeuDir[ip][im][irap]   = (TH1D *)hMvsPtvsRap_NeuDir[ip][im]->ProjectionY(Form("hPt_NeuDir%dp%dm_RapBin%d", ip, im, irap), rapBinLow, rapBinHi, massBinLow, massBinHi);

                massBinLow = hMvsPtvsRap_NeuDir[ip][im]->GetZaxis()->FindBin(mLowMassBandLow + mTinyNum);
                massBinHi  = hMvsPtvsRap_NeuDir[ip][im]->GetZaxis()->FindBin(mLowMassBandHi - mTinyNum);
                hLowMassBandPt_NeuDir[ip][im][irap] = (TH1D *)hMvsPtvsRap_NeuDir[ip][im]->ProjectionY(Form("hLowMassBandPt_NeuDir%dp%dm_RapBin%d", ip, im, irap), rapBinLow, rapBinHi, massBinLow, massBinHi);
                hLowMassBandPt_NeuDir[ip][im][irap]->SetTitle(Form("%1.1f < y < %1.1f", mDiffRapLow[irap], mDiffRapHi[irap]));

                massBinLow = hMvsPtvsRap_NeuDir[ip][im]->GetZaxis()->FindBin(mJpsiMassLow + mTinyNum);
                massBinHi  = hMvsPtvsRap_NeuDir[ip][im]->GetZaxis()->FindBin(mJpsiMassHi - mTinyNum);
                hJpsiPt_NeuDir[ip][im][irap] = (TH1D *)hMvsPtvsRap_NeuDir[ip][im]->ProjectionY(Form("hJpsiPt_NeuDir%dp%dm_RapBin%d", ip, im, irap), rapBinLow, rapBinHi, massBinLow, massBinHi);
                hJpsiPt_NeuDir[ip][im][irap]->SetTitle(Form("%1.1f < y < %1.1f", mDiffRapLow[irap], mDiffRapHi[irap]));

                massBinLow = hMvsPtvsRap_NeuDir[ip][im]->GetZaxis()->FindBin(mHiMassBandLow + mTinyNum);
                massBinHi  = hMvsPtvsRap_NeuDir[ip][im]->GetZaxis()->FindBin(mHiMassBandHi - mTinyNum);
                hHiMassBandPt_NeuDir[ip][im][irap] = (TH1D *)hMvsPtvsRap_NeuDir[ip][im]->ProjectionY(Form("hHiMassBandPt_NeuDir%dp%dm_RapBin%d", ip, im, irap), rapBinLow, rapBinHi, massBinLow, massBinHi);
                hHiMassBandPt_NeuDir[ip][im][irap]->SetTitle(Form("%1.1f < y < %1.1f", mDiffRapLow[irap], mDiffRapHi[irap]));

                if(ip >= im){
                    hMass_Neu[im][ip][irap]   = (TH1D *)hMass_NeuDir[im][ip][irap]->Clone(Form("hMass_Neu%dn%dn_RapBin%d", im, ip, irap));
                    hPt_Neu[im][ip][irap] = (TH1D *)hPt_NeuDir[im][ip][irap]->Clone(Form("hPt_Neu%dn%dn_RapBin%d", im, ip, irap));

                    if(ip > im){
                        hMass_Neu[im][ip][irap]->Add(hMass_NeuDir[ip][im][irap]);
                        hPt_Neu[im][ip][irap]->Add(hPt_NeuDir[ip][im][irap]);
                    }
                }
            }
        }
    }

    Int_t ineu = 0;
    for(Int_t ip=0; ip<nNeus; ip++){
        for(Int_t im=0; im<nNeus; im++){

            for(Int_t irap=0; irap<nDiffRapBins; irap++){
                c2->cd(irap+1);
                gPad->SetLogy(1);
                hMass_NeuDir[ip][im][irap]->RebinX(2);
                hMass_NeuDir[ip][im][irap]->GetYaxis()->SetRangeUser(0.5, hMass_NeuDir[ip][im][irap]->GetMaximum()*5);
                hMass_NeuDir[ip][im][irap]->GetXaxis()->SetTitleOffset(0.9);
                hMass_NeuDir[ip][im][irap]->GetYaxis()->SetTitle("Entries");

                setHisto(hMass_NeuDir[ip][im][irap], 20, 0.6, 1, 1);
                hMass_NeuDir[ip][im][irap]->Draw("p");

                drawLatex(0.36, 0.95, Form("%1.1f < y < %1.1f", mDiffRapLow[irap], mDiffRapHi[irap]), mFont, 0.07, 4);
                if(irap==0) drawLatex(0.56, 0.84, Form("%s", neuMultDirName[ineu].Data()), mFont, 0.06, 1);
            }

            ineu++;

            c2->SaveAs(Form("%s/massvsRap_NeuDir%dp%dm.png", dirName.Data(), ip, im));
            if(mStorePDF) c2->SaveAs(Form("%s/massvsRap_NeuDir%dp%dm.pdf", dirName.Data(), ip, im));
        }
    }

    ineu = 0;
    for(Int_t ip=0; ip<nNeus; ip++){
        for(Int_t im=0; im<nNeus; im++){

            for(Int_t irap=0; irap<nDiffRapBins; irap++){
                c2->cd(irap+1);
                gPad->SetLogy(1);
                hPt_NeuDir[ip][im][irap]->RebinX(2);
                hPt_NeuDir[ip][im][irap]->GetYaxis()->SetRangeUser(0.5, hPt_NeuDir[ip][im][irap]->GetMaximum()*5);
                hPt_NeuDir[ip][im][irap]->GetXaxis()->SetTitleOffset(0.9);
                hPt_NeuDir[ip][im][irap]->GetYaxis()->SetTitle("Entries");

                setHisto(hPt_NeuDir[ip][im][irap], 20, 0.6, 1, 1);
                hPt_NeuDir[ip][im][irap]->GetXaxis()->SetRangeUser(0, 3);
                hPt_NeuDir[ip][im][irap]->Draw("p");

                drawLatex(0.36, 0.95, Form("%1.1f < y < %1.1f", mDiffRapLow[irap], mDiffRapHi[irap]), mFont, 0.07, 4);
                if(irap==0){
                    drawLatex(0.56, 0.84, Form("%s", neuMultDirName[ineu].Data()), mFont, 0.06, 1);
                    drawLatex(0.56, 0.76, "2 < M_{#mu#mu} < 5 GeV", mFont, 0.06, 1);
                }
            }

            ineu++;

            c2->SaveAs(Form("%s/ptvsRap_NeuDir%dp%dm.png", dirName.Data(), ip, im));
            if(mStorePDF) c2->SaveAs(Form("%s/ptvsRap_NeuDir%dp%dm.pdf", dirName.Data(), ip, im));
        }
    }

    ineu = 0;
    for(Int_t ip=0; ip<nNeus; ip++){
        for(Int_t im=0; im<nNeus; im++){

            for(Int_t irap=0; irap<nDiffRapBins; irap++){
                c2->cd(irap+1);
                gPad->SetLogy(1);
                hLowMassBandPt_NeuDir[ip][im][irap]->RebinX(2);
                hHiMassBandPt_NeuDir[ip][im][irap]->RebinX(2);
                hJpsiPt_NeuDir[ip][im][irap]->RebinX(2);

                //hLowMassBandPt_NeuDir[ip][im][irap]->Scale(hJpsiPt_NeuDir[ip][im][irap]->GetEntries()*1./ hLowMassBandPt_NeuDir[ip][im][irap]->GetEntries());
                //hHiMassBandPt_NeuDir[ip][im][irap]->Scale(hJpsiPt_NeuDir[ip][im][irap]->GetEntries()*1./ hHiMassBandPt_NeuDir[ip][im][irap]->GetEntries());

                hJpsiPt_NeuDir[ip][im][irap]->GetYaxis()->SetRangeUser(0.5, hJpsiPt_NeuDir[ip][im][irap]->GetMaximum()*5);
                hJpsiPt_NeuDir[ip][im][irap]->GetXaxis()->SetTitleOffset(0.9);
                hJpsiPt_NeuDir[ip][im][irap]->GetYaxis()->SetTitle("Entries");

                setHisto(hJpsiPt_NeuDir[ip][im][irap], 20, 0.6, 1, 1);
                setHisto(hLowMassBandPt_NeuDir[ip][im][irap], 24, 0.6, 2, 2);
                setHisto(hHiMassBandPt_NeuDir[ip][im][irap], 25, 0.6, 4, 4);
                hJpsiPt_NeuDir[ip][im][irap]->GetXaxis()->SetRangeUser(0, 3);
                hJpsiPt_NeuDir[ip][im][irap]->Draw("p");
                hLowMassBandPt_NeuDir[ip][im][irap]->Draw("psame");
                hHiMassBandPt_NeuDir[ip][im][irap]->Draw("psame");

                drawLatex(0.36, 0.95, Form("%1.1f < y < %1.1f", mDiffRapLow[irap], mDiffRapHi[irap]), mFont, 0.07, 4);
                if(irap==0){
                    drawLatex(0.56, 0.84, Form("%s", neuMultDirName[ineu].Data()), mFont, 0.06, 1);

                    setLegend(leg, 0.48, 0.64, 0.8, 0.82, 0.05);
                    leg->AddEntry(hJpsiPt_NeuDir[ip][im][irap], Form("%1.2f < M_{#mu#mu} < %1.2f GeV", mJpsiMassLow, mJpsiMassHi), "pl");
                    leg->AddEntry(hLowMassBandPt_NeuDir[ip][im][irap], Form("%1.2f < M_{#mu#mu} < %1.2f GeV", mLowMassBandLow, mLowMassBandHi), "pl");
                    leg->AddEntry(hHiMassBandPt_NeuDir[ip][im][irap], Form("%1.2f < M_{#mu#mu} < %1.2f GeV", mHiMassBandLow, mHiMassBandHi), "pl");
                    leg->Draw("same");
                }
            }

            ineu++;

            c2->SaveAs(Form("%s/jpsiPtvsRap_NeuDir%dp%dm.png", dirName.Data(), ip, im));
            if(mStorePDF) c2->SaveAs(Form("%s/jpsiPtvsRap_NeuDir%dp%dm.pdf", dirName.Data(), ip, im));
        }
    }

    ineu = 0;
    for(Int_t ip=0; ip<nNeus; ip++){
        for(Int_t im=ip; im<nNeus; im++){

            for(Int_t irap=0; irap<nDiffRapBins; irap++){
                c2->cd(irap+1);
                gPad->SetLogy(1);
                hMass_Neu[ip][im][irap]->RebinX(2);
                hMass_Neu[ip][im][irap]->GetYaxis()->SetRangeUser(0.5, hMass_Neu[ip][im][irap]->GetMaximum()*5);
                hMass_Neu[ip][im][irap]->GetXaxis()->SetTitleOffset(0.9);
                hMass_Neu[ip][im][irap]->GetYaxis()->SetTitle("Entries");

                setHisto(hMass_Neu[ip][im][irap], 20, 0.6, 1, 1);
                hMass_Neu[ip][im][irap]->Draw("p");

                drawLatex(0.36, 0.95, Form("%1.1f < y < %1.1f", mDiffRapLow[irap], mDiffRapHi[irap]), mFont, 0.07, 4);
                if(irap==0) drawLatex(0.14, 0.84, Form("%s", neuMultName[ineu].Data()), mFont, 0.07, 4);
            }

            ineu++;

            c2->SaveAs(Form("%s/massvsRap_Neu%dn%dn.png", dirName.Data(), ip, im));
            if(mStorePDF) c2->SaveAs(Form("%s/massvsRap_Neu%dn%dn.pdf", dirName.Data(), ip, im));
        }
    }

    ineu = 0;
    for(Int_t ip=0; ip<nNeus; ip++){
        for(Int_t im=ip; im<nNeus; im++){

            for(Int_t irap=0; irap<nDiffRapBins; irap++){
                c2->cd(irap+1);
                gPad->SetLogy(1);
                hPt_Neu[ip][im][irap]->RebinX(2);
                hPt_Neu[ip][im][irap]->GetYaxis()->SetRangeUser(0.5, hPt_Neu[ip][im][irap]->GetMaximum()*5);
                hPt_Neu[ip][im][irap]->GetXaxis()->SetTitleOffset(0.9);
                hPt_Neu[ip][im][irap]->GetYaxis()->SetTitle("Entries");

                setHisto(hPt_Neu[ip][im][irap], 20, 0.6, 1, 1);
                hPt_Neu[ip][im][irap]->GetXaxis()->SetRangeUser(0, 3);
                hPt_Neu[ip][im][irap]->Draw("p");

                drawLatex(0.36, 0.95, Form("%1.1f < y < %1.1f", mDiffRapLow[irap], mDiffRapHi[irap]), mFont, 0.07, 4);
                if(irap==0){
                    drawLatex(0.32, 0.84, Form("%s", neuMultName[ineu].Data()), mFont, 0.07, 4);
                    drawLatex(0.32, 0.76, "2 < M_{#mu#mu} < 5 GeV", mFont, 0.07, 4);
                }
            }

            ineu++;

            c2->SaveAs(Form("%s/ptvsRap_Neu%dn%dn.png", dirName.Data(), ip, im));
            if(mStorePDF) c2->SaveAs(Form("%s/ptvsRap_Neu%dn%dn.pdf", dirName.Data(), ip, im));
        }
    }

    TFile *fOut = new TFile(Form("%s/sideBand.root", dirName.Data()), "recreate");
    for(Int_t ip=0; ip<nNeus; ip++){
        for(Int_t im=0; im<nNeus; im++){
            for(Int_t irap=0; irap<nDiffRapBins; irap++){
                hJpsiPt_NeuDir[ip][im][irap]->GetXaxis()->UnZoom();
                hLowMassBandPt_NeuDir[ip][im][irap]->GetXaxis()->UnZoom();
                hHiMassBandPt_NeuDir[ip][im][irap]->GetXaxis()->UnZoom();

                hJpsiPt_NeuDir[ip][im][irap]->Write();
                hLowMassBandPt_NeuDir[ip][im][irap]->Write();
                hHiMassBandPt_NeuDir[ip][im][irap]->Write();
            }
        }
    }

    cout << "End of program !" << endl;
}
