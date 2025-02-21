#include "/afs/ihep.ac.cn/users/z/zhangyu1/bishe/jpsiAnaCode/common/headers.h"
#include "/afs/ihep.ac.cn/users/z/zhangyu1/bishe/jpsiAnaCode/common/function.C"
#include "/afs/ihep.ac.cn/users/z/zhangyu1/bishe/jpsiAnaCode/common/funUtil.h"

const Double_t mTinyOffNum = 1.e-6;
const Double_t mOffSet = 0.1;

Double_t neuPur[nDirs][nNeus];
Double_t neuContamLow[nDirs][nNeus]; // e.g. for neuContamLow[0][2]: In Plus direction, for nNeu = 2 case, the contamination cased by nNeu = 1
Double_t neuContamHi[nDirs][nNeus];  // e.g. for neuContamHi[0][2]:  In Plus direction, for nNeu = 2 case, the contamination cased by nNeu = 3

const Int_t mMaxColors = 4;

const Double_t asyPhiCoreLow = 0;
const Double_t asyPhiCoreHi  = 0.006;
const Double_t asyPtCoreLow = 0;
const Double_t asyPtCoreHi  = 0.03;
const Double_t asyPt2CoreLow = 0;
const Double_t asyPt2CoreHi  = 0.01;

Int_t    mMarker[nMBins] = { 20, 21, 33 };
Int_t    mColor[nMBins]  = { 1, 2, 4 };
Double_t mSize[nMBins]   = { 1.2, 1.2, 1.6 };

void estimatePurity();

void plotDiMuon(Bool_t effCorr = kTRUE, Bool_t applyTnPSF = kFALSE, Bool_t incHadron = kFALSE, TString hfVetoType="Default")
{
    gStyle->SetOptFit(1111);

    if(!hfVetoType.EqualTo("Default") && !hfVetoType.EqualTo("Tight") && !hfVetoType.EqualTo("removeHF")){
        cout<<"Please input the correct hfVetoType string: 'Default' OR 'Tight' OR 'removeHF'"<<endl;
        return;
    }

    if(!effCorr && applyTnPSF){
        cout<<"Cannot apply TnP scaling factor without implementing efficiency corrections !"<<endl;
        return;
    }

    if(!init(hfVetoType)) {
        cout<<"Initialization failed !"<<endl;
        return;
    }

    estimatePurity();
    //neuPur[0][1] = 0.9;
    //neuPur[1][1] = 0.9;
    //neuContamHi[0][1] = 0.1;
    //neuContamHi[1][1] = 0.1;
    for(Int_t idir = 0; idir < nDirs; idir++){
        for(Int_t ineu = 0; ineu < nNeus; ineu++){
            cout << Form("%s - %dn: ", mDir[idir].Data(), ineu) << "   purity:" << neuPur[idir][ineu] << "   neuContamLow:" << neuContamLow[idir][ineu] << "   neuContamHi:" << neuContamHi[idir][ineu] << endl;
        }
    }

    TString histoDir = "dimuonHistos";
    TString plotDir  = "dimuonPlots";

    TString fileName = Form("%s/", histoDir.Data());
    TString dirName  = Form("%s/", plotDir.Data());

    if(effCorr){
        fileName += "effCorr";
        dirName  += "effCorr";

        if(applyTnPSF){
            fileName += ".applyTnPSF";
            dirName  += "_applyTnPSF";
        }
    }
    else{
        fileName += "rawSig";
        dirName  += "rawSig";
    }

    if(incHadron){
        fileName += ".incHadron";
        dirName  += "_incHadron";
    }

    if(hfVetoType.EqualTo("Tight")){
        fileName += ".tightHF";
        dirName  += "_tightHF";
    }

    if(hfVetoType.EqualTo("removeHF")){
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

    TCanvas* c3 = new TCanvas("c3", "c3", 1200, 900);
    c3->Divide(nNeus, nNeus);

    TLegend* leg = new TLegend(0.15, 0.72, 0.4, 0.9);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.06);

    TH1D* hnEvts = (TH1D *)f->Get("hnEvts");
    TH2D* hCenvsTrig = (TH2D*)f->Get("hCenvsTrig");
    TH2D* hPreScalevsTrig = (TH2D*)f->Get("hPreScalevsTrig");

    TH3D* hVzvsVyvsVx = (TH3D*)f->Get("hVzvsVyvsVx");
    TH3D* hVzvsVyvsVx_Sel = (TH3D*)f->Get("hVzvsVyvsVx_Sel"); // Vertex selection

    TH3D* hHFMinusvsHFPlusvsCen = (TH3D*)f->Get("hHFMinusvsHFPlusvsCen");
    TH3D* hHFMinusvsHFPlusvsCen_Sel = (TH3D*)f->Get("hHFMinusvsHFPlusvsCen_Sel");

    TH3D* hNtrkHPvsNtrkofflinevsCen = (TH3D *)f->Get("hNtrkHPvsNtrkofflinevsCen");
    TH3D* hNtrkHPvsNtrkofflinevsCen_Sel = (TH3D *)f->Get("hNtrkHPvsNtrkofflinevsCen_Sel");
    TH1D* hNtrkHP_2SoftMuons = (TH1D *)f->Get("hNtrkHP_2SoftMuons");

    TH2D* hZDCMinusvsZDCPlus = (TH2D*)f->Get("hZDCMinusvsZDCPlus");
    TH2D* hZDCMinusvsZDCPlus_Sel = (TH2D*)f->Get("hZDCMinusvsZDCPlus_Sel");
    TH2D* hZDCMinusvsZDCPlus_LS = (TH2D*)f->Get("hZDCMinusvsZDCPlus_LS");
    TH2D* hZDCMinusvsZDCPlus_Sel_LS = (TH2D*)f->Get("hZDCMinusvsZDCPlus_Sel_LS");
    TH2D* hZDCMinusvsZDCPlus_Only2MuTrk_LS = (TH2D *)f->Get("hZDCMinusvsZDCPlus_Only2MuTrk_LS");
    TH2D* hZDCPlusvsNeuNum  = (TH2D *)f->Get("hZDCPlusvsNeuNum");
    TH2D* hZDCMinusvsNeuNum = (TH2D *)f->Get("hZDCMinusvsNeuNum");

    TH3D* hCS_MvsPtvsCen  = (TH3D*)f->Get("hCS_MvsPtvsCen");
    TH3D* hCS_MvsPtvsRap  = (TH3D*)f->Get("hCS_MvsPtvsRap");
    TH3D* hWS_MvsPtvsCen  = (TH3D*)f->Get("hWS_MvsPtvsCen");

    TH3D* hMvsPtvsRap     = (TH3D*)f->Get("hMvsPtvsRap");
    TH2D* hAsyPhivsM   = (TH2D*)f->Get("hAsyPhivsM");

    TH1D *hSigmavsCostheta = (TH1D *)f->Get("hSigmavsCostheta");
    TH1D *hSigmavsM_mumu = (TH1D *)f->Get("hSigmavsM_mumu");
    TH1D *hSigmavsk_max = (TH1D *)f->Get("hSigmavsk_max");
    TH1D *hSigmavsk_min = (TH1D *)f->Get("hSigmavsk_min");

    TH2D* hAsyPhivsM_Neu[nNeus][nNeus];
    for (Int_t ip = 0; ip < nNeus; ip++) {
        for (Int_t im = 0; im < nNeus; im++) {
            hAsyPhivsM_Neu[ip][im] = (TH2D*)f->Get(Form("hAsyPhivsM_Neu%dp%dm", ip, im));
        }
    }

    TH2D* hAsyPhivsM_Exact2Neu[nNeus][nDirs];
    for (Int_t ineu = 0; ineu < nNeus; ineu++) {
        for (Int_t idir = 0; idir < nDirs; idir++) {
            hAsyPhivsM_Exact2Neu[ineu][idir] = (TH2D*)f->Get(Form("hAsyPhivsM_Neu%dnExact2%s", ineu, mDir[idir].Data()));
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

    c1->cd();
    gPad->SetLogz(1);
    hCenvsTrig->Draw("colz");
    c1->SaveAs(Form("%s/CenvsTrig.png", dirName.Data()));

    c1->cd();
    gPad->SetLogz(1);
    hPreScalevsTrig->Draw("colz");
    c1->SaveAs(Form("%s/PreScalevsTrig.png", dirName.Data()));

    c1->cd();
    hSigmavsCostheta->GetXaxis()->SetLabelSize(0.045);
    hSigmavsCostheta->GetXaxis()->SetRangeUser(-1, 1);
    hSigmavsCostheta->GetYaxis()->SetTitle("#frac{d#sigma}{dcos#theta}");
    hSigmavsCostheta->SetMarkerStyle(20);
    hSigmavsCostheta->SetMarkerColor(kRed);
    hSigmavsCostheta->SetMarkerSize(1);
    hSigmavsCostheta->Draw("p");
    c1->SaveAs(Form("%s/SigmavsCostheta.pdf", dirName.Data()));

    c1->cd();
    hSigmavsM_mumu->GetXaxis()->SetLabelSize(0.045);
    hSigmavsM_mumu->GetXaxis()->SetRangeUser(8, 60);
    hSigmavsM_mumu->GetYaxis()->SetTitle("#frac{d#sigma}{d#it{M}_{#mu#mu}}");
    hSigmavsCostheta->SetMarkerStyle(20);
    hSigmavsCostheta->SetMarkerColor(kRed);
    hSigmavsCostheta->SetMarkerSize(1);
    hSigmavsM_mumu->Draw("p");
    c1->SaveAs(Form("%s/SigmavsM_mumu.pdf", dirName.Data()));

    c1->cd();
    gPad->SetLogx(1);
    gPad->SetLogy(1);
    hSigmavsk_max->GetXaxis()->SetTitle("#it{k}_{min,max}[GeV]"); 
    hSigmavsk_max->GetYaxis()->SetTitle("#frac{d#sigma}{d#it{k}_{min,max}} [#mub/GeV]");
    hSigmavsk_max->SetMarkerStyle(20);
    hSigmavsk_max->SetMarkerColor(kRed);
    hSigmavsk_max->SetMarkerSize(1);
    hSigmavsk_max->GetXaxis()->SetLabelSize(0.045);
    hSigmavsk_max->GetXaxis()->SetRangeUser(0.1, 1e4);
    hSigmavsk_max->GetYaxis()->SetRangeUser(1e-4, 1e2);
    hSigmavsk_max->Draw("p");

    hSigmavsk_min->SetMarkerStyle(21);
    hSigmavsk_min->SetMarkerColor(kBlue);
    hSigmavsk_min->SetMarkerSize(1);
    hSigmavsk_min->GetXaxis()->SetLabelSize(0.045);
    hSigmavsk_min->GetXaxis()->SetRangeUser(0.1, 1e4);
    hSigmavsk_min->GetYaxis()->SetRangeUser(1e-4, 1e2);
    hSigmavsk_min->Draw("samep");
    TLegend *leg2 = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg2->AddEntry(hSigmavsk_max, "#it{k}_{max}", "p");
    leg2->AddEntry(hSigmavsk_min, "#it{k}_{min}", "p");
    leg2->SetTextSize(0.04);
    leg2->Draw();
    c1->SaveAs(Form("%s/hSigmavsk.pdf", dirName.Data()));

    TH2D* hVyvsVx = (TH2D*)hVzvsVyvsVx->Project3D("hVyvsVx_yx");
    TH1D* hVz = (TH1D*)hVzvsVyvsVx->ProjectionZ("hVz");
    TH2D* hVyvsVx_Sel = (TH2D*)hVzvsVyvsVx_Sel->Project3D("hVyvsVx_Sel_yx");
    TH1D* hVz_Sel = (TH1D*)hVzvsVyvsVx_Sel->ProjectionZ("hVz_Sel");
    c2->cd(1);
    gPad->SetLogz(1);
    hVyvsVx->GetXaxis()->SetRangeUser(-0.5, 0.5);
    hVyvsVx->GetYaxis()->SetRangeUser(-0.5, 0.5);
    hVyvsVx->Draw("colz");
    drawLatex(0.36, 0.95, "w/o Vertex Selection", 22, 0.06, 1);
    c2->cd(2);
    hVz->GetXaxis()->SetRangeUser(-30, 30);
    hVz->Draw("hist");
    drawLatex(0.36, 0.95, "w/o Vertex Selection", 22, 0.06, 1);
    c2->cd(3);
    gPad->SetLogz(1);
    hVyvsVx_Sel->GetXaxis()->SetRangeUser(-0.5, 0.5);
    hVyvsVx_Sel->GetYaxis()->SetRangeUser(-0.5, 0.5);
    hVyvsVx_Sel->Draw("colz");
    drawLatex(0.36, 0.95, "w/ Vertex Selection", 22, 0.06, 2);
    c2->cd(4);
    hVz_Sel->GetXaxis()->SetRangeUser(-30, 30);
    hVz_Sel->Draw("hist");
    drawLatex(0.36, 0.95, "w/ Vertex Selection", 22, 0.06, 2);
    drawLatex(0.16, 0.2, Form("Event Rejection: %1.1f%%", 100 * (1 - hVz_Sel->GetEntries() * 1. / hVz->GetEntries())), 22, 0.07, 4);
    c2->SaveAs(Form("%s/vertex.png", dirName.Data()));

    TH2D *hHFMinusvsHFPlus = (TH2D*)hHFMinusvsHFPlusvsCen->Project3D("hHFMinusvsHFPlus_zy");
    TH2D *hHFMinusvsHFPlus_Sel = (TH2D*)hHFMinusvsHFPlusvsCen_Sel->Project3D("hHFMinusvsHFPlus_Sel_zy");

    TH1D *hCen       = (TH1D *)hHFMinusvsHFPlusvsCen->ProjectionX("hCen");
    TH1D *hCen_Sel   = (TH1D *)hHFMinusvsHFPlusvsCen_Sel->ProjectionX("hCen_Sel");
    TH1D* hCen_Final = (TH1D *)f->Get("hCen_Final");

    clearPad(c2, nPads);

    c2->cd(1);
    gPad->SetLogz(1);
    hHFMinusvsHFPlus->Draw("colz");
    drawLatex(0.36, 0.95, "w/ vertex selection", 22, 0.06, 1);
    c2->cd(2);
    gPad->SetLogz(1);
    hHFMinusvsHFPlus_Sel->GetXaxis()->SetRangeUser(0, 20);
    hHFMinusvsHFPlus_Sel->GetYaxis()->SetRangeUser(0, 20);
    hHFMinusvsHFPlus_Sel->Draw("colz");
    drawLatex(0.28, 0.95, "w/ vertex&maxHFE selection", 22, 0.06, 1);
    c2->cd(3);
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
    leg->AddEntry(hCen_Final, "+ N_{trk}^{HP} == 2 & N_{#mu} == 2", "l");
    leg->Draw("same");
    c2->SaveAs(Form("%s/HFsumET.png", dirName.Data()));

    TH2D *hNtrkHPvsNtrkoffline     = (TH2D *)hNtrkHPvsNtrkofflinevsCen->Project3D("hNtrkHPvsNtrkoffline_zy");
    TH2D *hNtrkHPvsNtrkoffline_Sel = (TH2D *)hNtrkHPvsNtrkofflinevsCen_Sel->Project3D("hNtrkHPvsNtrkoffline_Sel_zy");
    TH1D *hNtrkHP                  = (TH1D *)hNtrkHPvsNtrkofflinevsCen->ProjectionZ("hNtrkHP");
    TH1D *hNtrkHP_Sel              = (TH1D *)hNtrkHPvsNtrkofflinevsCen_Sel->ProjectionZ("hNtrkHP_Sel");

    setLegend(leg, 0.24, 0.72, 0.6, 0.90, 0.055);

    c2->cd(1);
    gPad->SetLogz(1);
    hNtrkHPvsNtrkoffline->GetXaxis()->SetRangeUser(0, 100);
    hNtrkHPvsNtrkoffline->GetYaxis()->SetRangeUser(0, 100);
    hNtrkHPvsNtrkoffline->Draw("colz");
    drawLatex(0.36, 0.95, "w/ vertex selection", 22, 0.06, 1);
    c2->cd(2);
    gPad->SetLogz(1);
    hNtrkHPvsNtrkoffline_Sel->GetXaxis()->SetRangeUser(0, 30);
    hNtrkHPvsNtrkoffline_Sel->GetYaxis()->SetRangeUser(0, 30);
    hNtrkHPvsNtrkoffline_Sel->Draw("colz");
    drawLatex(0.28, 0.95, "w/ vertex&maxHFE selection", 22, 0.06, 1);
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
    leg->AddEntry(hNtrkHP_2SoftMuons, "+ N_{#mu} == 2", "l");
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
    //hNtrkHP_2SoftMuons->Draw("histesame");
    leg->Clear();
    leg->AddEntry(hNtrkHP, "!BeamHalo & validVtx", "l");
    leg->AddEntry(hNtrkHP_Sel, "+ max HF Eenergy < 7.6 (7.3) GeV", "l");
    //leg->AddEntry(hNtrkHP_2SoftMuons, "+ N_{#mu} == 2", "l");
    leg->AddEntry((TObject*)0, "", "");
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

    c1->cd();
    gPad->SetLogy(0);
    gPad->SetLogz(1);
    hZDCPlusvsNeuNum->GetXaxis()->SetRangeUser(-0.5, nNeus-0.5);
    hZDCPlusvsNeuNum->Draw("colz");
    c1->SaveAs(Form("%s/ZDCPlusvsNeuNum.png", dirName.Data()));

    c1->cd();
    gPad->SetLogy(0);
    gPad->SetLogz(1);
    hZDCMinusvsNeuNum->GetXaxis()->SetRangeUser(-0.5, nNeus-0.5);
    hZDCMinusvsNeuNum->Draw("colz");
    c1->SaveAs(Form("%s/ZDCMinusvsNeuNum.png", dirName.Data()));

    TH2D* hCS_MvsPt = (TH2D*)hCS_MvsPtvsCen->Project3D("hCS_MvsPt_zy");
    TH2D* hWS_MvsPt = (TH2D*)hWS_MvsPtvsCen->Project3D("hWS_MvsPt_zy");
    TH1D* hCS_Mass = (TH1D*)hCS_MvsPtvsCen->ProjectionZ("hCS_Mass", 0, -1, 0, -1);
    TH1D* hWS_Mass = (TH1D*)hWS_MvsPtvsCen->ProjectionZ("hWS_Mass", 0, -1, 0, -1);

    c2->cd(1);
    gPad->SetLogy(0);
    gPad->SetLogz(1);
    hCS_MvsPt->GetYaxis()->SetRangeUser(5, 60);
    hCS_MvsPt->Draw("col");
    drawLatex(0.4, 0.95, "Correct-Sign", 22, 0.06, 1);
    c2->cd(2);
    gPad->SetLogy(0);
    gPad->SetLogz(1);
    hWS_MvsPt->GetYaxis()->SetRangeUser(5, 60);
    hWS_MvsPt->Draw("col");
    drawLatex(0.4, 0.95, "Wrong-Sign", 22, 0.06, 1);
    c2->cd(3);
    gPad->SetLogy(1);
    setLegend(leg, 0.5, 0.76, 0.84, 0.88, 0.06);
    leg->AddEntry(hCS_Mass, "Correct-Sign", "p");
    leg->AddEntry(hWS_Mass, "Wrong-Sign", "l");
    setHisto(hCS_Mass, 20, 0.5, 1, 1, 2);
    setHisto(hWS_Mass, 24, 0.5, 2, 2, 2);
    hWS_Mass->SetLineStyle(2);
    hCS_Mass->RebinX(2);
    hWS_Mass->RebinX(2);
    //hCS_Mass->GetXaxis()->SetRangeUser(6, 80);
    hCS_Mass->GetXaxis()->SetRangeUser(6, 60);
    hCS_Mass->GetYaxis()->SetRangeUser(0.5, 5e3);
    hCS_Mass->GetYaxis()->SetTitle(Form("Entries / (%1.1f GeV/c^{2})", hCS_Mass->GetXaxis()->GetBinWidth(1)));
    hCS_Mass->Draw("p");
    hWS_Mass->Draw("histsame");
    leg->Draw("same");
    c2->SaveAs(Form("%s/MvsPt.png", dirName.Data()));

    TH1D* hRawRap[nMBins];
    TH1D* hRap[nMBins];
    TH1D* hAsyPhi[nMBins];
    for (Int_t i = 0; i < nMBins; i++) {
        Int_t massBinLow = hMvsPtvsRap->GetZaxis()->FindBin(massLow[i] + mTinyOffNum);
        Int_t massBinHi  = hMvsPtvsRap->GetZaxis()->FindBin(massHi[i] - mTinyOffNum);
        hRawRap[i] = (TH1D*)hCS_MvsPtvsRap->ProjectionX(Form("hRawRap_massBin%d", i), 0, -1, massBinLow, massBinHi);
        hRap[i]    = (TH1D*)hMvsPtvsRap->ProjectionX(Form("hRap_massBin%d", i), 0, -1, massBinLow, massBinHi);

        massBinLow = hAsyPhivsM->GetXaxis()->FindBin(massLow[i] + mTinyOffNum);
        massBinHi  = hAsyPhivsM->GetXaxis()->FindBin(massHi[i] - mTinyOffNum);
        hAsyPhi[i] = (TH1D*)hAsyPhivsM->ProjectionY(Form("hAsyPhi_massBin%d", i), massBinLow, massBinHi);
        //cout << "AsyPhi OverFlow Fraction: " << hAsyPhi[i]->GetBinContent(hAsyPhi[i]->GetNbinsX() + 1) * 1. / hAsyPhi[i]->Integral(0, -1) << endl;
    }

    c1->cd();
    gPad->SetLogy(1);
    setLegend(leg, 0.32, 0.16, 0.64, 0.32, 0.05);
    for(Int_t i=0; i<nMBins; i++){
        setHisto(hRawRap[i], mMarker[i], mSize[i], mColor[i], mColor[i], 2);
        hRawRap[i]->GetYaxis()->SetRangeUser(5e-1, 1e4);
        if(i==0){
            hRawRap[i]->Draw("p");
            hRawRap[i]->GetYaxis()->SetTitle(Form("Entries / %1.1f", hRawRap[i]->GetXaxis()->GetBinWidth(1)));
        }
        else     hRawRap[i]->Draw("psame");
        leg->AddEntry(hRawRap[i], Form("%1.0f<M_{#mu#mu}<%1.0f GeV", massLow[i], massHi[i]), "p");
    }
    leg->Draw("same");
    c1->SaveAs(Form("%s/rawRapidityvsMass.png", dirName.Data()));

    c1->cd();
    gPad->SetLogy(1);
    for(Int_t i=0; i<nMBins; i++){
        setHisto(hRap[i], mMarker[i], mSize[i], mColor[i], mColor[i], 2);
        hRap[i]->GetYaxis()->SetRangeUser(5e-1, 1e4);
        if(i==0){
            hRap[i]->Draw("p");
            hRap[i]->GetYaxis()->SetTitle(Form("Entries / %1.1f", hRawRap[i]->GetXaxis()->GetBinWidth(1)));
        }
        else     hRap[i]->Draw("psame");
    }
    leg->Draw("same");
    c1->SaveAs(Form("%s/rapidityvsMass.png", dirName.Data()));

    c1->cd();
    gPad->SetLogy(1);
    setLegend(leg, 0.56, 0.72, 0.88, 0.88, 0.05);
    for (Int_t i = 0; i < nMBins; i++) {
        setHisto(hAsyPhi[i], mMarker[i], mSize[i], mColor[i], mColor[i], 2);
        hAsyPhi[i]->Scale(1. / hAsyPhi[i]->GetBinContent(1));
        hAsyPhi[i]->SetMinimum(5.e-4);
        hAsyPhi[i]->GetXaxis()->SetRangeUser(asyPhiCoreLow, asyPhiCoreHi);
        if (i == 0)
            hAsyPhi[i]->Draw("hist");
        else
            hAsyPhi[i]->Draw("histsame");

        leg->AddEntry(hAsyPhi[i], Form("%1.0f<M_{#mu#mu}<%1.0f GeV", massLow[i], massHi[i]), "l");
    }
    leg->DrawClone("same");
    c1->SaveAs(Form("%s/AsyPhivsMass.png", dirName.Data()));

    TH2D *hNorAsyPhivsM_Neu[nNeus][nNeus];
    TH2D *hAsyPhivsM_ContaCorr_Neu[nNeus][nNeus];  // correct for neutron contamination
    TH2D *hAsyPhivsM_PileupCorr_Neu[nNeus][nNeus]; // correct for dissociative pileup
    TH2D *hAsyPhivsM_ContaPileupCorr_Neu[nNeus][nNeus]; // correct for neutron contamination and dissociative pileup
    Double_t meaStats[nNeus][nNeus];
    for(Int_t ip=0; ip<nNeus; ip++){
        for(Int_t im=0; im<nNeus; im++){
            meaStats[ip][im] = hAsyPhivsM_Neu[ip][im]->Integral(0, -1, 0, -1);

            hNorAsyPhivsM_Neu[ip][im] = (TH2D *)hAsyPhivsM_Neu[ip][im]->Clone(Form("hNorAsyPhivsM_Neu%dp%dm", ip, im));
            hNorAsyPhivsM_Neu[ip][im]->Scale(1./meaStats[ip][im]);

            hAsyPhivsM_ContaCorr_Neu[ip][im] = (TH2D *)hAsyPhivsM_Neu[ip][im]->Clone(Form("hAsyPhivsM_ContaCorr_Neu%dp%dm", ip, im));
            hAsyPhivsM_PileupCorr_Neu[ip][im] = (TH2D *)hAsyPhivsM_Neu[ip][im]->Clone(Form("hAsyPhivsM_PileupCorr_Neu%dp%dm", ip, im));
            hAsyPhivsM_ContaPileupCorr_Neu[ip][im] = (TH2D *)hAsyPhivsM_Neu[ip][im]->Clone(Form("hAsyPhivsM_ContaPileupCorr_Neu%dp%dm", ip, im));
        }
    }

    // *** correct neutron contamination caused by ZDC resolution ***
    // normalized to unit
    for(Int_t ineu=0; ineu<nNeus; ineu++){
        for(Int_t idir=0; idir<nDirs; idir++){
            hAsyPhivsM_Exact2Neu[ineu][idir]->Scale(1./hAsyPhivsM_Exact2Neu[ineu][idir]->Integral(0, -1, 0, -1));
        }
    }

    for(Int_t ip=0; ip<nNeus; ip++){
        for(Int_t im=0; im<nNeus; im++){
            if(ip!=1 && im!=1) continue;

            hAsyPhivsM_ContaCorr_Neu[ip][im]->Reset();

            if((ip==1 && im!=1) || (im==1 && ip!=1)){
                Int_t idir=-1, ineu=-1;
                if(ip==1){
                    idir = 0; // neutrons in ZDCPlus side has contamination
                    ineu = im;
                }
                else if(im==1){
                    idir = 1; // neutrons in ZDCMinus side has contamination
                    ineu = ip;
                }

                // also fill underflow and overflow
                for(Int_t ibinX=0; ibinX<=hAsyPhivsM_ContaCorr_Neu[ip][im]->GetNbinsX()+1; ibinX++){
                    for(Int_t ibinY=0; ibinY<=hAsyPhivsM_ContaCorr_Neu[ip][im]->GetNbinsY()+1; ibinY++){
                        Double_t num1    = hNorAsyPhivsM_Neu[ip][im]->GetBinContent(ibinX, ibinY);
                        Double_t numErr1 = hNorAsyPhivsM_Neu[ip][im]->GetBinError(ibinX, ibinY);
                        Double_t num2    = hAsyPhivsM_Exact2Neu[ineu][idir]->GetBinContent(ibinX, ibinY);
                        Double_t numErr2 = hAsyPhivsM_Exact2Neu[ineu][idir]->GetBinError(ibinX, ibinY);

                        Double_t num    = num1 - neuContamHi[idir][1]*num2;
                        Double_t numErr = sqrt( pow(numErr1, 2) + pow(neuContamHi[idir][1]*numErr2, 2) );
                        Double_t den    = neuPur[idir][1];

                        hAsyPhivsM_ContaCorr_Neu[ip][im]->SetBinContent(ibinX, ibinY, num/den);
                        hAsyPhivsM_ContaCorr_Neu[ip][im]->SetBinError(ibinX, ibinY, numErr/den);
                    }
                }
            }

            if(ip==1 && im==1){
                // also fill underflow and overflow
                for(Int_t ibinX=0; ibinX<=hAsyPhivsM_ContaCorr_Neu[ip][im]->GetNbinsX()+1; ibinX++){
                    for(Int_t ibinY=0; ibinY<=hAsyPhivsM_ContaCorr_Neu[ip][im]->GetNbinsY()+1; ibinY++){
                        Double_t num1    = hNorAsyPhivsM_Neu[ip][im]->GetBinContent(ibinX, ibinY);
                        Double_t numErr1 = hNorAsyPhivsM_Neu[ip][im]->GetBinError(ibinX, ibinY);
                        Double_t num2    = hAsyPhivsM_Exact2Neu[im][0]->GetBinContent(ibinX, ibinY); // Plus side has exact two neutron
                        Double_t numErr2 = hAsyPhivsM_Exact2Neu[im][0]->GetBinError(ibinX, ibinY);   // Plus side has exact two neutron
                        Double_t num3    = hAsyPhivsM_Exact2Neu[ip][1]->GetBinContent(ibinX, ibinY); // Minus side has exact two neutron
                        Double_t numErr3 = hAsyPhivsM_Exact2Neu[ip][1]->GetBinError(ibinX, ibinY);   // Minus side has exact two neutron

                        Double_t num    = num1 - neuPur[1][1]*neuContamHi[0][1]*num2 - neuPur[0][1]*neuContamHi[1][1]*num3;
                        Double_t numErr = sqrt( pow(numErr1, 2) + pow(neuPur[1][1]*neuContamHi[0][1]*numErr2, 2) + pow(neuPur[0][1]*neuContamHi[1][1]*numErr3, 2) );
                        Double_t den    = neuPur[0][1]*neuPur[1][1];

                        hAsyPhivsM_ContaCorr_Neu[ip][im]->SetBinContent(ibinX, ibinY, num/den);
                        hAsyPhivsM_ContaCorr_Neu[ip][im]->SetBinError(ibinX, ibinY, numErr/den);
                    }
                }
            }

            //cout << Form("(ip, im) = (%d, %d): ", ip, im) << hAsyPhi_ContaCorr_Neu[ip][im]->Integral(0, -1) << endl;

            hAsyPhivsM_ContaCorr_Neu[ip][im]->Scale(meaStats[ip][im]); // scale back to measured statistics
        }
    }

    // *** correct pileup effect caused by neuclear dissociation without any activity in CMS tracker ***
    for(Int_t ip=0; ip<nNeus; ip++){
        for(Int_t im=0; im<nNeus; im++){
            cout<<Form("%dp%dm:", ip, im)<<endl;

            Int_t iRow = matrixIdx[ip][im];

            // also fill underflow and overflow
            for(Int_t ibinX=0; ibinX<=hAsyPhivsM_PileupCorr_Neu[ip][im]->GetNbinsX()+1; ibinX++){
                for(Int_t ibinY=0; ibinY<=hAsyPhivsM_PileupCorr_Neu[ip][im]->GetNbinsY()+1; ibinY++){
                    Double_t binCon=0;
                    Double_t binErr=0;
                    for(Int_t jp=0; jp<nNeus; jp++){
                        for(Int_t jm=0; jm<nNeus; jm++){
                            Int_t jCol = matrixIdx[jp][jm];
                            binCon += corrMatrix[iRow][jCol] * hAsyPhivsM_Neu[jp][jm]->GetBinContent(ibinX, ibinY);
                            binErr += pow(corrMatrix[iRow][jCol] * hAsyPhivsM_Neu[jp][jm]->GetBinError(ibinX, ibinY), 2);
                        }
                    }
                    binErr = sqrt(binErr);

                    hAsyPhivsM_PileupCorr_Neu[ip][im]->SetBinContent(ibinX, ibinY, binCon);
                    hAsyPhivsM_PileupCorr_Neu[ip][im]->SetBinError(ibinX, ibinY, binErr);
                }
            }

            cout<<"rawStat(withNeuConta): "<<hAsyPhivsM_Neu[ip][im]->Integral(0, -1, 0, -1)<<"        pileupCorrStat(withNeuConta): "<<hAsyPhivsM_PileupCorr_Neu[ip][im]->Integral(0, -1, 0, -1)<<endl;

            // also fill underflow and overflow
            for(Int_t ibinX=0; ibinX<=hAsyPhivsM_ContaPileupCorr_Neu[ip][im]->GetNbinsX()+1; ibinX++){
                for(Int_t ibinY=0; ibinY<=hAsyPhivsM_ContaPileupCorr_Neu[ip][im]->GetNbinsY()+1; ibinY++){
                    Double_t binCon=0;
                    Double_t binErr=0;
                    for(Int_t jp=0; jp<nNeus; jp++){
                        for(Int_t jm=0; jm<nNeus; jm++){
                            Int_t jCol = matrixIdx[jp][jm];
                            binCon += corrMatrix[iRow][jCol] * hAsyPhivsM_ContaCorr_Neu[jp][jm]->GetBinContent(ibinX, ibinY);
                            binErr += pow(corrMatrix[iRow][jCol] * hAsyPhivsM_ContaCorr_Neu[jp][jm]->GetBinError(ibinX, ibinY), 2);
                        }
                    }
                    binErr = sqrt(binErr);

                    hAsyPhivsM_ContaPileupCorr_Neu[ip][im]->SetBinContent(ibinX, ibinY, binCon);
                    hAsyPhivsM_ContaPileupCorr_Neu[ip][im]->SetBinError(ibinX, ibinY, binErr);
                }
            }

            cout<<"rawStat(corrNeuConta): "<<hAsyPhivsM_ContaCorr_Neu[ip][im]->Integral(0, -1, 0, -1)<<"        pileupCorrStat(corrNeuConta): "<<hAsyPhivsM_ContaPileupCorr_Neu[ip][im]->Integral(0, -1, 0, -1)<<endl;
        }
    }
    // *** done ***

    TH1D *hM_Neu[nNeus][nNeus];
    TH1D *hM_PileupCorr_Neu[nNeus][nNeus];
    TH1D *hM_ContaPileupCorr_Neu[nNeus][nNeus];
    TH1D *hM_Compare_Neu[nNeus][nNeus];
    TH1D *hM_PileupCorr_Compare_Neu[nNeus][nNeus];
    TH1D *hM_ContaPileupCorr_Compare_Neu[nNeus][nNeus];
    for(Int_t ip=0; ip<nNeus; ip++){
        for(Int_t im=0; im<nNeus; im++){
            hM_Neu[ip][im] = (TH1D *)hAsyPhivsM_Neu[ip][im]->ProjectionX(Form("hM_Neu%dp%dm", ip, im), 0, -1);
            hM_PileupCorr_Neu[ip][im] = (TH1D *)hAsyPhivsM_PileupCorr_Neu[ip][im]->ProjectionX(Form("hM_PileupCorr_Neu%dp%dm", ip, im), 0, -1);
            hM_ContaPileupCorr_Neu[ip][im] = (TH1D *)hAsyPhivsM_ContaPileupCorr_Neu[ip][im]->ProjectionX(Form("hM_ContaPileupCorr_Neu%dp%dm", ip, im), 0, -1);

            hM_PileupCorr_Compare_Neu[ip][im] = (TH1D *)hM_PileupCorr_Neu[ip][im]->Clone(Form("hM_PileupCorr_Compare_Neu%dp%dm", ip, im));
            hM_Compare_Neu[ip][im] = (TH1D *)hM_Neu[ip][im]->Clone(Form("hM_Compare_Neu%dp%dm", ip, im));
            hM_ContaPileupCorr_Compare_Neu[ip][im] = (TH1D *)hM_ContaPileupCorr_Neu[ip][im]->Clone(Form("hM_ContaPileupCorr_Compare_Neu%dp%dm", ip, im));
            hM_PileupCorr_Compare_Neu[ip][im]->RebinX(20);
            hM_Compare_Neu[ip][im]->RebinX(20);
            hM_ContaPileupCorr_Compare_Neu[ip][im]->RebinX(20);
        }
    }

    setLegend(leg, 0.15, 0.15, 0.35, 0.28, 0.065);
    leg->SetNColumns(1);

    TString neuPlus, neuMinus;
    Int_t   neuIdx = 0;
    for(Int_t ip=0; ip<nNeus; ip++){
        for(Int_t im=0; im<nNeus; im++){
            c3->cd(neuIdx+1);
            gPad->SetLogy(1);
            setHisto(hM_PileupCorr_Compare_Neu[ip][im], 20, 1, 1, 1, 2);
            setHisto(hM_Compare_Neu[ip][im], 24, 1, 2, 2, 2);

            Double_t ymin = 0.5;
            Double_t ymax = hM_PileupCorr_Compare_Neu[ip][im]->GetBinContent(hM_PileupCorr_Compare_Neu[ip][im]->GetXaxis()->FindBin(massLow[0]+mTinyOffNum))*5;
            hM_PileupCorr_Compare_Neu[ip][im]->GetXaxis()->SetRangeUser(massLow[0], massHi[nMBins-1]);
            hM_PileupCorr_Compare_Neu[ip][im]->GetYaxis()->SetRangeUser(ymin, ymax);
            hM_PileupCorr_Compare_Neu[ip][im]->GetYaxis()->SetTitle(Form("Entries / (%1.1f GeV)", hM_PileupCorr_Compare_Neu[ip][im]->GetXaxis()->GetBinWidth(1)));
            hM_PileupCorr_Compare_Neu[ip][im]->Draw("p");
            hM_Compare_Neu[im][ip]->Draw("psame");

            //if(ip == nNeus-1) neuPlus = Form("#geq%d", ip);
            if(ip == nNeus-1) neuPlus = "X";
            else              neuPlus = Form("%d", ip);
            //if(im == nNeus-1) neuMinus = Form("#geq%d", im);
            if(im == nNeus-1) neuMinus = "X";
            else              neuMinus = Form("%d", im);

            if(ip+im==0){
                leg->AddEntry(hM_PileupCorr_Compare_Neu[ip][im], "Corrected pileup", "p");
                leg->AddEntry(hM_Compare_Neu[ip][im], "With pileup", "p");
                leg->Draw("same");
            }

            drawLatex(0.46, 0.95, Form("%sp%sm", neuPlus.Data(), neuMinus.Data()), 22, 0.07, 4);

            neuIdx++;
        }
    }
    c3->SaveAs(Form("%s/massCompare_PileupCorr.png", dirName.Data()));
    c3->SaveAs(Form("%s/massCompare_PileupCorr.pdf", dirName.Data()));

    neuIdx=0;
    for(Int_t ip=0; ip<nNeus; ip++){
        for(Int_t im=0; im<nNeus; im++){
            c3->cd(neuIdx+1);
            gPad->SetLogy(1);
            setHisto(hM_PileupCorr_Compare_Neu[ip][im], 20, 1, 1, 1, 2);
            setHisto(hM_ContaPileupCorr_Compare_Neu[ip][im], 24, 1, 2, 2, 2);

            //Int_t binLow = hM_PileupCorr_Compare_Neu[ip][im]->GetXaxis()->FindBin(massLow[0] + mTinyOffNum);
            //Int_t binHi  = hM_PileupCorr_Compare_Neu[ip][im]->GetXaxis()->FindBin(massHi[nMBins-1] + mTinyOffNum);
            //cout<<ip<<"p"<<im<<"m:"<<endl;
            //cout<<"PileupCorr: "<<hM_PileupCorr_Compare_Neu[ip][im]->Integral(binLow, binHi)<<endl;
            //cout<<"ContaPileupCorr: "<<hM_ContaPileupCorr_Compare_Neu[ip][im]->Integral(binLow, binHi)<<endl;
            //cout<<endl;

            hM_PileupCorr_Compare_Neu[ip][im]->Draw("p");
            hM_ContaPileupCorr_Compare_Neu[im][ip]->Draw("psame");

            //if(ip == nNeus-1) neuPlus = Form("#geq%d", ip);
            if(ip == nNeus-1) neuPlus = "X";
            else              neuPlus = Form("%d", ip);
            //if(im == nNeus-1) neuMinus = Form("#geq%d", im);
            if(im == nNeus-1) neuMinus = "X";
            else              neuMinus = Form("%d", im);

            if(ip+im==0){
                leg->Clear();
                leg->AddEntry(hM_PileupCorr_Compare_Neu[ip][im], "Corrected pileup", "p");
                leg->AddEntry(hM_ContaPileupCorr_Compare_Neu[ip][im], "Corrected pileup&contamination", "p");
                leg->Draw("same");
            }

            drawLatex(0.46, 0.95, Form("%sp%sm", neuPlus.Data(), neuMinus.Data()), 22, 0.07, 4);

            neuIdx++;
        }
    }
    c3->SaveAs(Form("%s/massCompare_ContaPileupCorr.png", dirName.Data()));
    c3->SaveAs(Form("%s/massCompare_ContaPileupCorr.pdf", dirName.Data()));

    neuIdx = 0;
    for(Int_t ip=0; ip<nNeus; ip++){
        for(Int_t im=ip+1; im<nNeus; im++){
            c2->cd(neuIdx+1);
            gPad->SetLogy(1);
            setHisto(hM_PileupCorr_Compare_Neu[ip][im], 20, 1, 1, 1, 2);
            setHisto(hM_PileupCorr_Compare_Neu[im][ip], 24, 1, 2, 2, 2);
            hM_PileupCorr_Compare_Neu[ip][im]->GetXaxis()->SetRangeUser(massLow[0], massHi[nMBins-1]);
            hM_PileupCorr_Compare_Neu[ip][im]->Draw("p");
            hM_PileupCorr_Compare_Neu[im][ip]->Draw("psame");

            neuPlus = Form("%d", ip);
            if(im == nNeus-1) neuMinus = "X";
            else              neuMinus = Form("%d", im);

            leg->Clear();
            leg->AddEntry(hM_PileupCorr_Compare_Neu[ip][im], Form("%sp%sm", neuPlus.Data(), neuMinus.Data()), "p");
            leg->AddEntry(hM_PileupCorr_Compare_Neu[im][ip], Form("%sp%sm", neuMinus.Data(), neuPlus.Data()), "p");
            leg->DrawClone("same");

            neuIdx++;
        }
    }
    c2->SaveAs(Form("%s/massCompare_NeuAsy.png", dirName.Data()));
    c2->SaveAs(Form("%s/massCompare_NeuAsy.pdf", dirName.Data()));

    TGraphErrors* grMeanMvsNeuNum = new TGraphErrors(nPts);
    grMeanMvsNeuNum->SetName("grMeanMvsNeuNum");

    TGraphErrors* grMeanMvsNeuNum_PileupCorr = new TGraphErrors(nPts);
    grMeanMvsNeuNum_PileupCorr->SetName("grMeanMvsNeuNum_PileupCorr");

    TGraphErrors* grMeanMvsNeuNum_ContaPileupCorr = new TGraphErrors(nPts);
    grMeanMvsNeuNum_ContaPileupCorr->SetName("grMeanMvsNeuNum_ContaPileupCorr");

    TH1D *hM_Clone_Neu[nNeus][nNeus];
    TH1D *hM_PileupCorr_Clone_Neu[nNeus][nNeus];
    TH1D *hM_ContaPileupCorr_Clone_Neu[nNeus][nNeus];

    neuIdx = 0;
    for(Int_t ip=0; ip<nNeus; ip++){
        for(Int_t im=ip; im<nNeus; im++){
            hM_Neu[ip][im]->SetName(Form("hM_Neu%dn%dn", ip, im));
            hM_PileupCorr_Neu[ip][im]->SetName(Form("hM_PileupCorr_Neu%dn%dn", ip, im));
            hM_ContaPileupCorr_Neu[ip][im]->SetName(Form("hM_ContaPileupCorr_Neu%dn%dn", ip, im));

            if(im > ip){
                hM_PileupCorr_Neu[ip][im]->Add(hM_PileupCorr_Neu[im][ip]);
                hM_PileupCorr_Neu[im][ip]->Reset();

                hM_Neu[ip][im]->Add(hM_Neu[im][ip]);
                hM_Neu[im][ip]->Reset();

                hM_ContaPileupCorr_Neu[ip][im]->Add(hM_ContaPileupCorr_Neu[im][ip]);
                hM_ContaPileupCorr_Neu[im][ip]->Reset();
            }

            hM_Clone_Neu[ip][im] = (TH1D *)hM_Neu[ip][im]->Clone(Form("hM_Clone_Neu%dn%dn", ip, im));
            hM_PileupCorr_Clone_Neu[ip][im] = (TH1D *)hM_PileupCorr_Neu[ip][im]->Clone(Form("hM_PileupCorr_Clone_Neu%dn%dn", ip, im));
            hM_ContaPileupCorr_Clone_Neu[ip][im] = (TH1D *)hM_ContaPileupCorr_Neu[ip][im]->Clone(Form("hM_ContaPileupCorr_Clone_Neu%dn%dn", ip, im));

            hM_Clone_Neu[ip][im]->GetXaxis()->SetRangeUser(massLow[0], massHi[nMBins-1]);
            grMeanMvsNeuNum->SetPoint(neuIdx, neuIdx+mOffSet, hM_Clone_Neu[ip][im]->GetMean());
            grMeanMvsNeuNum->SetPointError(neuIdx, 0, hM_Clone_Neu[ip][im]->GetMeanError());

            hM_PileupCorr_Clone_Neu[ip][im]->GetXaxis()->SetRangeUser(massLow[0], massHi[nMBins-1]);
            grMeanMvsNeuNum_PileupCorr->SetPoint(neuIdx, neuIdx, hM_PileupCorr_Clone_Neu[ip][im]->GetMean());
            grMeanMvsNeuNum_PileupCorr->SetPointError(neuIdx, 0, hM_PileupCorr_Clone_Neu[ip][im]->GetMeanError());

            hM_ContaPileupCorr_Clone_Neu[ip][im]->GetXaxis()->SetRangeUser(massLow[0], massHi[nMBins-1]);
            grMeanMvsNeuNum_ContaPileupCorr->SetPoint(neuIdx, neuIdx+mOffSet, hM_ContaPileupCorr_Clone_Neu[ip][im]->GetMean());
            grMeanMvsNeuNum_ContaPileupCorr->SetPointError(neuIdx, 0, hM_ContaPileupCorr_Clone_Neu[ip][im]->GetMeanError());

            hM_Clone_Neu[ip][im]->RebinX(20);
            hM_PileupCorr_Clone_Neu[ip][im]->RebinX(20);
            hM_ContaPileupCorr_Clone_Neu[ip][im]->RebinX(20);

            neuIdx++;
        }
    }

    setLegend(leg, 0.42, 0.76, 0.88, 0.88, 0.06);
    leg->SetNColumns(3);

    c1->cd();
    neuIdx = 0;
    for(Int_t ip=0; ip<nNeus; ip++){
        for(Int_t im=ip; im<nNeus; im++){
            Int_t icolor = neuIdx+1;
            if(icolor==5)  icolor = 42;
            if(icolor==10) icolor = 28;

            hM_PileupCorr_Clone_Neu[ip][im]->GetXaxis()->SetRangeUser(massLow[0], massHi[nMBins-1]);
            Int_t binIdx = hM_PileupCorr_Clone_Neu[ip][im]->GetXaxis()->FindBin(massLow[0] + mTinyOffNum);
            hM_PileupCorr_Clone_Neu[ip][im]->Scale(1./hM_PileupCorr_Clone_Neu[ip][im]->GetBinContent(binIdx));

            setHisto(hM_PileupCorr_Clone_Neu[ip][im], 20, 0.8, icolor, icolor, 2);
            if(ip+im==0) hM_PileupCorr_Clone_Neu[ip][im]->Draw("p");
            else         hM_PileupCorr_Clone_Neu[ip][im]->Draw("psame");

            if(ip == nNeus-1) neuPlus = "X";
            else              neuPlus = Form("%d", ip);
            if(im == nNeus-1) neuMinus = "X";
            else              neuMinus = Form("%d", im);
            leg->AddEntry(hM_PileupCorr_Clone_Neu[ip][im], Form("%sn%sn", neuPlus.Data(), neuMinus.Data()), "p");

            neuIdx++;
        }
    }
    leg->Draw("same");
    c1->SaveAs(Form("%s/massvsNeuNum_PileupCorr.png", dirName.Data()));

    c1->cd();
    neuIdx = 0;
    for(Int_t ip=0; ip<nNeus; ip++){
        for(Int_t im=ip; im<nNeus; im++){
            Int_t icolor = neuIdx+1;
            if(icolor==5)  icolor = 42;
            if(icolor==10) icolor = 28;

            hM_Clone_Neu[ip][im]->GetXaxis()->SetRangeUser(massLow[0], massHi[nMBins-1]);
            Int_t binIdx = hM_Clone_Neu[ip][im]->GetXaxis()->FindBin(massLow[0] + mTinyOffNum);
            hM_Clone_Neu[ip][im]->Scale(1./hM_Clone_Neu[ip][im]->GetBinContent(binIdx));

            setHisto(hM_Clone_Neu[ip][im], 20, 0.8, icolor, icolor, 2);
            if(ip+im==0) hM_Clone_Neu[ip][im]->Draw("p");
            else         hM_Clone_Neu[ip][im]->Draw("psame");

            neuIdx++;
        }
    }
    leg->Draw("same");
    c1->SaveAs(Form("%s/massvsNeuNum.png", dirName.Data()));

    c1->cd();
    neuIdx = 0;
    for(Int_t ip=0; ip<nNeus; ip++){
        for(Int_t im=ip; im<nNeus; im++){
            Int_t icolor = neuIdx+1;
            if(icolor==5)  icolor = 42;
            if(icolor==10) icolor = 28;

            hM_ContaPileupCorr_Clone_Neu[ip][im]->GetXaxis()->SetRangeUser(massLow[0], massHi[nMBins-1]);
            Int_t binIdx = hM_ContaPileupCorr_Clone_Neu[ip][im]->GetXaxis()->FindBin(massLow[0] + mTinyOffNum);
            hM_ContaPileupCorr_Clone_Neu[ip][im]->Scale(1./hM_ContaPileupCorr_Clone_Neu[ip][im]->GetBinContent(binIdx));

            setHisto(hM_ContaPileupCorr_Clone_Neu[ip][im], 20, 0.8, icolor, icolor, 2);
            if(ip+im==0) hM_ContaPileupCorr_Clone_Neu[ip][im]->Draw("p");
            else         hM_ContaPileupCorr_Clone_Neu[ip][im]->Draw("psame");

            neuIdx++;
        }
    }
    leg->Draw("same");
    c1->SaveAs(Form("%s/massvsNeuNum_ContaPileup.png", dirName.Data()));

    TH2D *ddM = (TH2D *)histo("ddM", nPts, -0.5, nPts-0.5, 60, 12, 15, "", "#LTM_{#mu#mu}#GT (GeV/c^{2})");
    neuIdx = 0;
    for(Int_t ip=0; ip<nNeus; ip++){
        for(Int_t im=ip; im<nNeus; im++){
            if(ip == nNeus-1) neuPlus = "X";
            else              neuPlus = Form("%d", ip);
            if(im == nNeus-1) neuMinus = "X";
            else              neuPlus = Form("%d", ip);
            ddM->GetXaxis()->SetBinLabel(neuIdx+1, Form("%sn%sn", neuPlus.Data(), neuMinus.Data()));
            neuIdx++;
        }
    }
    ddM->GetXaxis()->LabelsOption("d");
    ddM->GetXaxis()->SetLabelSize(0.07);

    setLegend(leg, 0.15, 0.78, 0.35, 0.88, 0.05);
    leg->SetNColumns(1);

    c1->cd();
    gPad->SetLogy(0);
    ddM->Draw("c");
    setGraph(grMeanMvsNeuNum_PileupCorr, 20, 1.0, 1, 1, 2);
    setGraph(grMeanMvsNeuNum, 24, 1.0, 2, 2, 2);
    grMeanMvsNeuNum_PileupCorr->Draw("pzsame");
    grMeanMvsNeuNum->Draw("pzsame");
    leg->AddEntry(grMeanMvsNeuNum_PileupCorr, "Corrected pileup", "p");
    leg->AddEntry(grMeanMvsNeuNum, "With pileup", "p");
    leg->Draw("same");
    drawLatex(0.40, 0.955, Form("%1.0f<M_{#mu#mu}<%1.0f GeV", massLow[0], massHi[nMBins-1]), 22, 0.06, 1);
    c1->SaveAs(Form("%s/meanMassvsNeuNum_PileupCorr.pdf", dirName.Data()));
    c1->SaveAs(Form("%s/meanMassvsNeuNum_PileupCorr.png", dirName.Data()));

    c1->cd();
    gPad->SetLogy(0);
    ddM->Draw("c");
    setGraph(grMeanMvsNeuNum_PileupCorr, 20, 1.0, 1, 1, 2);
    setGraph(grMeanMvsNeuNum_ContaPileupCorr, 24, 1.0, 2, 2, 2);
    grMeanMvsNeuNum_PileupCorr->Draw("pzsame");
    grMeanMvsNeuNum_ContaPileupCorr->Draw("pzsame");
    leg->Clear();
    leg->AddEntry(grMeanMvsNeuNum_PileupCorr, "Corrected pileup", "p");
    leg->AddEntry(grMeanMvsNeuNum_ContaPileupCorr, "Corrected pileu&contamination", "p");
    leg->Draw("same");
    drawLatex(0.40, 0.955, Form("%1.0f<M_{#mu#mu}<%1.0f GeV", massLow[0], massHi[nMBins-1]), 22, 0.06, 1);
    c1->SaveAs(Form("%s/meanMassvsNeuNum_ContaPileupCorr.pdf", dirName.Data()));
    c1->SaveAs(Form("%s/meanMassvsNeuNum_ContaPileupCorr.png", dirName.Data()));

    TH1D *hAsyPhi_Neu[nNeus][nNeus];
    TH1D *hAsyPhi_PileupCorr_Neu[nNeus][nNeus];
    TH1D *hAsyPhi_ContaPileupCorr_Neu[nNeus][nNeus];
    TH1D *hAsyPhi_Compare_Neu[nNeus][nNeus];
    TH1D *hAsyPhi_PileupCorr_Compare_Neu[nNeus][nNeus];
    TH1D *hAsyPhi_ContaPileupCorr_Compare_Neu[nNeus][nNeus];
    for(Int_t ip=0; ip<nNeus; ip++){
        for(Int_t im=0; im<nNeus; im++){
            Int_t massBinLow = hAsyPhivsM_Neu[ip][im]->GetXaxis()->FindBin(massLow[0] + mTinyOffNum);
            Int_t massBinHi  = hAsyPhivsM_Neu[ip][im]->GetXaxis()->FindBin(massHi[nMBins-1] - mTinyOffNum);

            hAsyPhi_Neu[ip][im] = (TH1D *)hAsyPhivsM_Neu[ip][im]->ProjectionY(Form("hAsyPhi_Neu%dp%dm", ip, im), massBinLow, massBinHi);
            hAsyPhi_PileupCorr_Neu[ip][im] = (TH1D *)hAsyPhivsM_PileupCorr_Neu[ip][im]->ProjectionY(Form("hAsyPhi_PileupCorr_Neu%dp%dm", ip, im), massBinLow, massBinHi);
            hAsyPhi_ContaPileupCorr_Neu[ip][im] = (TH1D *)hAsyPhivsM_ContaPileupCorr_Neu[ip][im]->ProjectionY(Form("hAsyPhi_ContaPileupCorr_Neu%dp%dm", ip, im), massBinLow, massBinHi);

            hAsyPhi_Compare_Neu[ip][im] = (TH1D *)hAsyPhi_Neu[ip][im]->Clone(Form("hAsyPhi_Compare_Neu%dp%dm", ip, im));
            hAsyPhi_PileupCorr_Compare_Neu[ip][im] = (TH1D *)hAsyPhi_PileupCorr_Neu[ip][im]->Clone(Form("hAsyPhi_PileupCorr_Compare_Neu%dp%dm", ip, im));
            hAsyPhi_ContaPileupCorr_Compare_Neu[ip][im] = (TH1D *)hAsyPhi_ContaPileupCorr_Neu[ip][im]->Clone(Form("hAsyPhi_ContaPileupCorr_Compare_Neu%dp%dm", ip, im));
            hAsyPhi_PileupCorr_Compare_Neu[ip][im]->RebinX(5);
            hAsyPhi_Compare_Neu[ip][im]->RebinX(5);
            hAsyPhi_ContaPileupCorr_Compare_Neu[ip][im]->RebinX(5);
        }
    }

    setLegend(leg, 0.15, 0.15, 0.35, 0.28, 0.065);
    leg->SetNColumns(1);

    neuIdx = 0;
    for(Int_t ip=0; ip<nNeus; ip++){
        for(Int_t im=0; im<nNeus; im++){
            c3->cd(neuIdx+1);
            gPad->SetLogy(1);
            setHisto(hAsyPhi_PileupCorr_Compare_Neu[ip][im], 20, 1, 1, 1, 2);
            setHisto(hAsyPhi_Compare_Neu[ip][im], 24, 1, 2, 2, 2);
            hAsyPhi_PileupCorr_Compare_Neu[ip][im]->GetXaxis()->SetRangeUser(0, 0.02);
            hAsyPhi_PileupCorr_Compare_Neu[ip][im]->GetYaxis()->SetRangeUser(0.5, hAsyPhi_PileupCorr_Compare_Neu[ip][im]->GetMaximum()*5);
            hAsyPhi_PileupCorr_Compare_Neu[ip][im]->GetYaxis()->SetTitle(Form("Entries / (%1.4f)", hAsyPhi_PileupCorr_Compare_Neu[ip][im]->GetXaxis()->GetBinWidth(1)));
            hAsyPhi_PileupCorr_Compare_Neu[ip][im]->Draw("p");
            hAsyPhi_Compare_Neu[im][ip]->Draw("psame");

            if(ip == nNeus-1) neuPlus = "X";
            else              neuPlus = Form("%d", ip);
            if(im == nNeus-1) neuMinus = "X";
            else              neuMinus = Form("%d", im);

            if(ip+im==0){
                drawLatex(0.52, 0.8, Form("%1.0f<M_{#mu#mu}<%1.0f GeV", massLow[0], massHi[nMBins-1]), 22, 0.07, 1);

                leg->AddEntry(hAsyPhi_PileupCorr_Compare_Neu[ip][im], "Corrected pileup", "p");
                leg->AddEntry(hAsyPhi_Compare_Neu[ip][im], "With pileup", "p");
                leg->Draw("same");
            }

            drawLatex(0.46, 0.95, Form("%sp%sm", neuPlus.Data(), neuMinus.Data()), 22, 0.07, 4);

            neuIdx++;
        }
    }
    c3->SaveAs(Form("%s/asyPhiCompare_PileupCorr.png", dirName.Data()));
    c3->SaveAs(Form("%s/asyPhiCompare_PileupCorr.pdf", dirName.Data()));

    neuIdx = 0;
    for(Int_t ip=0; ip<nNeus; ip++){
        for(Int_t im=0; im<nNeus; im++){
            c3->cd(neuIdx+1);
            gPad->SetLogy(1);
            setHisto(hAsyPhi_PileupCorr_Compare_Neu[ip][im], 20, 1, 1, 1, 2);
            setHisto(hAsyPhi_ContaPileupCorr_Compare_Neu[ip][im], 24, 1, 2, 2, 2);

            hAsyPhi_PileupCorr_Compare_Neu[ip][im]->Draw("p");
            hAsyPhi_ContaPileupCorr_Compare_Neu[im][ip]->Draw("psame");

            if(ip == nNeus-1) neuPlus = "X";
            else              neuPlus = Form("%d", ip);
            if(im == nNeus-1) neuMinus = "X";
            else              neuMinus = Form("%d", im);

            if(ip+im==0){
                drawLatex(0.52, 0.8, Form("%1.0f<M_{#mu#mu}<%1.0f GeV", massLow[0], massHi[nMBins-1]), 22, 0.07, 1);

                leg->Clear();
                leg->AddEntry(hAsyPhi_PileupCorr_Compare_Neu[ip][im], "Corrected pileup", "p");
                leg->AddEntry(hAsyPhi_ContaPileupCorr_Compare_Neu[ip][im], "Corrected pileup&contamination", "p");
                leg->Draw("same");
            }

            drawLatex(0.46, 0.95, Form("%sp%sm", neuPlus.Data(), neuMinus.Data()), 22, 0.07, 4);

            neuIdx++;
        }
    }
    c3->SaveAs(Form("%s/asyPhiCompare_ContaPileupCorr.png", dirName.Data()));
    c3->SaveAs(Form("%s/asyPhiCompare_ContaPileupCorr.pdf", dirName.Data()));

    neuIdx = 0;
    for(Int_t ip=0; ip<nNeus; ip++){
        for(Int_t im=ip+1; im<nNeus; im++){
            c2->cd(neuIdx+1);
            gPad->SetLogy(1);
            setHisto(hAsyPhi_PileupCorr_Compare_Neu[ip][im], 20, 1, 1, 1, 2);
            setHisto(hAsyPhi_PileupCorr_Compare_Neu[im][ip], 24, 1, 2, 2, 2);
            hAsyPhi_PileupCorr_Compare_Neu[ip][im]->GetXaxis()->SetRangeUser(0, 20e-3);
            hAsyPhi_PileupCorr_Compare_Neu[ip][im]->Draw("p");
            hAsyPhi_PileupCorr_Compare_Neu[im][ip]->Draw("psame");

            neuPlus = Form("%d", ip);
            if(im == nNeus-1) neuMinus = "X";
            else              neuMinus = Form("%d", im);

            leg->Clear();
            leg->AddEntry(hAsyPhi_PileupCorr_Compare_Neu[ip][im], Form("%sp%sm", neuPlus.Data(), neuMinus.Data()), "p");
            leg->AddEntry(hAsyPhi_PileupCorr_Compare_Neu[im][ip], Form("%sp%sm", neuMinus.Data(), neuPlus.Data()), "p");
            leg->DrawClone("same");

            drawLatex(0.40, 0.955, Form("%1.0f<M_{#mu#mu}<%1.0f GeV", massLow[0], massHi[nMBins-1]), 22, 0.06, 1);
            neuIdx++;
        }
    }
    c2->SaveAs(Form("%s/asyPhiCompare_NeuAsy.png", dirName.Data()));
    c2->SaveAs(Form("%s/asyPhiCompare_NeuAsy.pdf", dirName.Data()));

    TGraphErrors* grMeanAsyPhivsNeuNum = new TGraphErrors(nPts); // asyPhi range [0, 6e-3]
    grMeanAsyPhivsNeuNum->SetName("grMeanAsyPhivsNeuNum");

    TGraphErrors* grMeanAsyPhivsNeuNum_PileupCorr = new TGraphErrors(nPts); // asyPhi range [0, 6e-3]
    grMeanAsyPhivsNeuNum_PileupCorr->SetName("grMeanAsyPhivsNeuNum_PileupCorr");

    TGraphErrors* grMeanAsyPhivsNeuNum_ContaPileupCorr = new TGraphErrors(nPts); // asyPhi range [0, 6e-3]
    grMeanAsyPhivsNeuNum_ContaPileupCorr->SetName("grMeanAsyPhivsNeuNum_ContaPileupCorr");

    TH1D *hAsyPhi_Clone_Neu[nNeus][nNeus];
    TH1D *hAsyPhi_PileupCorr_Clone_Neu[nNeus][nNeus];
    TH1D *hAsyPhi_ContaPileupCorr_Clone_Neu[nNeus][nNeus];

    neuIdx = 0;
    for(Int_t ip=0; ip<nNeus; ip++){
        for(Int_t im=ip; im<nNeus; im++){
            hAsyPhi_Neu[ip][im]->SetName(Form("hAsyPhi_Neu%dn%dn", ip, im));
            hAsyPhi_PileupCorr_Neu[ip][im]->SetName(Form("hAsyPhi_PileupCorr_Neu%dn%dn", ip, im));
            hAsyPhi_ContaPileupCorr_Neu[ip][im]->SetName(Form("hAsyPhi_ContaPileupCorr_Neu%dn%dn", ip, im));

            if(im > ip){
                //hAsyPhi_PileupCorr_Neu[ip][im]->GetXaxis()->SetRangeUser(asyPhiCoreLow, asyPhiCoreHi);
                //hAsyPhi_PileupCorr_Neu[im][ip]->GetXaxis()->SetRangeUser(asyPhiCoreLow, asyPhiCoreHi);
                //cout << Form("(%dp, %dm): ", ip, im) << hAsyPhi_PileupCorr_Neu[ip][im]->GetMean() << " +/- " << hAsyPhi_PileupCorr_Neu[ip][im]->GetMeanError() << endl;
                //cout << Form("(%dp, %dm): ", im, ip) << hAsyPhi_PileupCorr_Neu[im][ip]->GetMean() << " +/- " << hAsyPhi_PileupCorr_Neu[im][ip]->GetMeanError() << endl;
                //cout << endl;

                hAsyPhi_Neu[ip][im]->Add(hAsyPhi_Neu[im][ip]);
                hAsyPhi_Neu[im][ip]->Reset();

                hAsyPhi_PileupCorr_Neu[ip][im]->Add(hAsyPhi_PileupCorr_Neu[im][ip]);
                hAsyPhi_PileupCorr_Neu[im][ip]->Reset();

                hAsyPhi_ContaPileupCorr_Neu[ip][im]->Add(hAsyPhi_ContaPileupCorr_Neu[im][ip]);
                hAsyPhi_ContaPileupCorr_Neu[im][ip]->Reset();
            }

            hAsyPhi_Clone_Neu[ip][im] = (TH1D *)hAsyPhi_Neu[ip][im]->Clone(Form("hAsyPhi_Clone_Neu%dn%dn", ip, im));
            hAsyPhi_PileupCorr_Clone_Neu[ip][im] = (TH1D *)hAsyPhi_PileupCorr_Neu[ip][im]->Clone(Form("hAsyPhi_PileupCorr_Clone_Neu%dn%dn", ip, im));
            hAsyPhi_ContaPileupCorr_Clone_Neu[ip][im] = (TH1D *)hAsyPhi_ContaPileupCorr_Neu[ip][im]->Clone(Form("hAsyPhi_ContaPileupCorr_Clone_Neu%dn%dn", ip, im));

            hAsyPhi_Clone_Neu[ip][im]->GetXaxis()->SetRangeUser(asyPhiCoreLow, asyPhiCoreHi);
            grMeanAsyPhivsNeuNum->SetPoint(neuIdx, neuIdx+mOffSet, hAsyPhi_Clone_Neu[ip][im]->GetMean());
            grMeanAsyPhivsNeuNum->SetPointError(neuIdx, 0, hAsyPhi_Clone_Neu[ip][im]->GetMeanError());

            hAsyPhi_PileupCorr_Clone_Neu[ip][im]->GetXaxis()->SetRangeUser(asyPhiCoreLow, asyPhiCoreHi);
            grMeanAsyPhivsNeuNum_PileupCorr->SetPoint(neuIdx, neuIdx, hAsyPhi_PileupCorr_Clone_Neu[ip][im]->GetMean());
            grMeanAsyPhivsNeuNum_PileupCorr->SetPointError(neuIdx, 0, hAsyPhi_PileupCorr_Clone_Neu[ip][im]->GetMeanError());

            hAsyPhi_ContaPileupCorr_Clone_Neu[ip][im]->GetXaxis()->SetRangeUser(asyPhiCoreLow, asyPhiCoreHi);
            grMeanAsyPhivsNeuNum_ContaPileupCorr->SetPoint(neuIdx, neuIdx+mOffSet, hAsyPhi_ContaPileupCorr_Clone_Neu[ip][im]->GetMean());
            grMeanAsyPhivsNeuNum_ContaPileupCorr->SetPointError(neuIdx, 0, hAsyPhi_ContaPileupCorr_Clone_Neu[ip][im]->GetMeanError());

            hAsyPhi_Clone_Neu[ip][im]->RebinX(5);
            hAsyPhi_PileupCorr_Clone_Neu[ip][im]->RebinX(5);
            hAsyPhi_ContaPileupCorr_Clone_Neu[ip][im]->RebinX(5);

            neuIdx++;
        }
    }

    setLegend(leg, 0.42, 0.76, 0.88, 0.88, 0.06);
    leg->SetNColumns(3);

    c1->cd();
    gPad->SetLogy(1);
    neuIdx = 0;
    for(Int_t ip=0; ip<nNeus; ip++){
        for(Int_t im=ip; im<nNeus; im++){
            Int_t icolor = neuIdx+1;
            if(icolor==5)  icolor = 42;
            if(icolor==10) icolor = 28;

            hAsyPhi_PileupCorr_Clone_Neu[ip][im]->Scale(1./hAsyPhi_PileupCorr_Clone_Neu[ip][im]->GetBinContent(1));

            hAsyPhi_PileupCorr_Clone_Neu[ip][im]->GetXaxis()->SetRangeUser(0, 0.01);
            hAsyPhi_PileupCorr_Clone_Neu[ip][im]->SetMinimum(1.e-3);
            setHisto(hAsyPhi_PileupCorr_Clone_Neu[ip][im], 20, 0.8, icolor, icolor, 2);

            if(ip+im==0) hAsyPhi_PileupCorr_Clone_Neu[ip][im]->Draw("p");
            else         hAsyPhi_PileupCorr_Clone_Neu[ip][im]->Draw("psame");

            if(ip == nNeus-1) neuPlus = "X";
            else              neuPlus = Form("%d", ip);
            if(im == nNeus-1) neuMinus = "X";
            else              neuMinus = Form("%d", im);
            leg->AddEntry(hAsyPhi_PileupCorr_Clone_Neu[ip][im], Form("%sn%sn", neuPlus.Data(), neuMinus.Data()), "p");

            neuIdx++;
        }
    }
    leg->Draw("same");
    drawLatex(0.40, 0.955, Form("%1.0f<M_{#mu#mu}<%1.0f GeV", massLow[0], massHi[nMBins-1]), 22, 0.06, 1);
    c1->SaveAs(Form("%s/asyPhivsNeuNum_PileupCorr.png", dirName.Data()));

    c1->cd();
    gPad->SetLogy(1);
    neuIdx = 0;
    for(Int_t ip=0; ip<nNeus; ip++){
        for(Int_t im=ip; im<nNeus; im++){
            Int_t icolor = neuIdx+1;
            if(icolor==5)  icolor = 42;
            if(icolor==10) icolor = 28;

            hAsyPhi_Clone_Neu[ip][im]->Scale(1./hAsyPhi_Clone_Neu[ip][im]->GetBinContent(1));

            hAsyPhi_Clone_Neu[ip][im]->GetXaxis()->SetRangeUser(0, 0.01);
            hAsyPhi_Clone_Neu[ip][im]->SetMinimum(1.e-3);
            setHisto(hAsyPhi_Clone_Neu[ip][im], 20, 0.8, icolor, icolor, 2);

            if(ip+im==0) hAsyPhi_Clone_Neu[ip][im]->Draw("p");
            else         hAsyPhi_Clone_Neu[ip][im]->Draw("psame");

            neuIdx++;
        }
    }
    leg->Draw("same");
    drawLatex(0.40, 0.955, Form("%1.0f<M_{#mu#mu}<%1.0f GeV", massLow[0], massHi[nMBins-1]), 22, 0.06, 1);
    c1->SaveAs(Form("%s/asyPhivsNeuNum.png", dirName.Data()));

    c1->cd();
    gPad->SetLogy(1);
    neuIdx = 0;
    for(Int_t ip=0; ip<nNeus; ip++){
        for(Int_t im=ip; im<nNeus; im++){
            Int_t icolor = neuIdx+1;
            if(icolor==5)  icolor = 42;
            if(icolor==10) icolor = 28;

            hAsyPhi_ContaPileupCorr_Clone_Neu[ip][im]->Scale(1./hAsyPhi_ContaPileupCorr_Clone_Neu[ip][im]->GetBinContent(1));

            hAsyPhi_ContaPileupCorr_Clone_Neu[ip][im]->GetXaxis()->SetRangeUser(0, 0.01);
            hAsyPhi_ContaPileupCorr_Clone_Neu[ip][im]->SetMinimum(1.e-3);
            setHisto(hAsyPhi_ContaPileupCorr_Clone_Neu[ip][im], 20, 0.8, icolor, icolor, 2);

            if(ip+im==0) hAsyPhi_ContaPileupCorr_Clone_Neu[ip][im]->Draw("p");
            else         hAsyPhi_ContaPileupCorr_Clone_Neu[ip][im]->Draw("psame");

            neuIdx++;
        }
    }
    leg->Draw("same");
    drawLatex(0.40, 0.955, Form("%1.0f<M_{#mu#mu}<%1.0f GeV", massLow[0], massHi[nMBins-1]), 22, 0.06, 1);
    c1->SaveAs(Form("%s/asyPhivsNeuNum_ContaPileupCorr.png", dirName.Data()));

    TH2D *ddAsyPhi = (TH2D *)histo("ddAsyPhi", nPts, -0.5, nPts-0.5, 50, 1.1e-3, 1.6e-3, "", "#LT#alpha#GT");
    neuIdx = 0;
    for(Int_t ip=0; ip<nNeus; ip++){
        for(Int_t im=ip; im<nNeus; im++){
            if(ip == nNeus-1) neuPlus = "X";
            else              neuPlus = Form("%d", ip);
            if(im == nNeus-1) neuMinus = "X";
            else              neuMinus = Form("%d", im);
            ddAsyPhi->GetXaxis()->SetBinLabel(neuIdx+1, Form("%sn%sn", neuPlus.Data(), neuMinus.Data()));
            neuIdx++;
        }
    }
    ddAsyPhi->GetXaxis()->LabelsOption("d");
    ddAsyPhi->GetXaxis()->SetLabelSize(0.07);

    setLegend(leg, 0.15, 0.78, 0.35, 0.88, 0.05);
    leg->SetNColumns(1);

    c1->cd();
    gPad->SetLogy(0);
    ddAsyPhi->Draw("c");
    setGraph(grMeanAsyPhivsNeuNum_PileupCorr, 20, 1.2, 1, 1, 2);
    setGraph(grMeanAsyPhivsNeuNum, 24, 1.2, 1, 1, 2);
    grMeanAsyPhivsNeuNum_PileupCorr->Draw("pzsame");
    grMeanAsyPhivsNeuNum->Draw("pzsame");
    leg->AddEntry(grMeanAsyPhivsNeuNum_PileupCorr, "Corrected pileup", "p");
    leg->AddEntry(grMeanAsyPhivsNeuNum, "With pileup", "p");
    leg->Draw("same");
    drawLatex(0.40, 0.955, Form("%1.0f<M_{#mu#mu}<%1.0f GeV", massLow[0], massHi[nMBins-1]), 22, 0.06, 1);
    drawLatex(0.72, 0.18, Form("#alpha < %1.3f", asyPhiCoreHi), 22, 0.06, 1);
    c1->SaveAs(Form("%s/meanAsyPhivsNeuNum_PileupCorr.pdf", dirName.Data()));
    c1->SaveAs(Form("%s/meanAsyPhivsNeuNum_PileupCorr.png", dirName.Data()));

    c1->cd();
    gPad->SetLogy(0);
    ddAsyPhi->Draw("c");
    setGraph(grMeanAsyPhivsNeuNum_PileupCorr, 20, 1.2, 1, 1, 2);
    setGraph(grMeanAsyPhivsNeuNum_ContaPileupCorr, 24, 1.2, 1, 1, 2);
    grMeanAsyPhivsNeuNum_PileupCorr->Draw("pzsame");
    grMeanAsyPhivsNeuNum_ContaPileupCorr->Draw("pzsame");
    leg->Clear();
    leg->AddEntry(grMeanAsyPhivsNeuNum_PileupCorr, "Corrected pileup", "p");
    leg->AddEntry(grMeanAsyPhivsNeuNum_ContaPileupCorr, "Corrected pileup&contamination", "p");
    leg->Draw("same");
    drawLatex(0.40, 0.955, Form("%1.0f<M_{#mu#mu}<%1.0f GeV", massLow[0], massHi[nMBins-1]), 22, 0.06, 1);
    drawLatex(0.72, 0.18, Form("#alpha < %1.3f", asyPhiCoreHi), 22, 0.06, 1);
    c1->SaveAs(Form("%s/meanAsyPhivsNeuNum_ContaPileupCorr.pdf", dirName.Data()));
    c1->SaveAs(Form("%s/meanAsyPhivsNeuNum_ContaPileupCorr.png", dirName.Data()));

    TFile *fOut = new TFile(Form("%s/asyPhiAndMass_vs_NeuNum.root", dirName.Data()), "recreate");
    for(Int_t ip=0; ip<nNeus; ip++){
        for(Int_t im=ip; im<nNeus; im++){
            hAsyPhi_Neu[ip][im]->Write();
            hAsyPhi_PileupCorr_Neu[ip][im]->Write();
            hAsyPhi_ContaPileupCorr_Neu[ip][im]->Write();
            hM_Neu[ip][im]->Write();
            hM_PileupCorr_Neu[ip][im]->Write();
            hM_ContaPileupCorr_Neu[ip][im]->Write();
        }
    }

    cout << "End of program !" << endl;
}

void estimatePurity()
{
    memset(neuPur, 0, sizeof(neuPur));
    memset(neuContamLow, 0, sizeof(neuPur));
    memset(neuContamHi,  0, sizeof(neuPur));

    Int_t binLow, binHi;
    for(Int_t idir = 0; idir < nDirs; idir++){
        Double_t binWidth = hZDC[idir]->GetBinWidth(1);

        for(Int_t ineu = 0; ineu < nNeus; ineu++){
            binLow   = hZDC[idir]->FindBin(mNeuZDCLow[idir][ineu] + mTinyOffNum);
            binHi    = hZDC[idir]->FindBin(mNeuZDCHi[idir][ineu] - mTinyOffNum);
            
            //当ZDC判断中子数为0n时 没有更低的中子数对判断造成影响 and neuContamHi是 1n Gauss分布的尾巴在0n区间的部分
            if(ineu == 0){
                neuContamLow[idir][ineu] = 0;
                neuContamHi[idir][ineu]  = singleGaus[idir][ineu]->Integral(mNeuZDCLow[idir][ineu], mNeuZDCHi[idir][ineu]) / binWidth / hZDC[idir]->Integral(binLow, binHi); // the gaussian function index for "nNeu = 1" is 0, and so on
            }
            //当ZDC判断中子数为Xn时 没有更高的中子数对判断造成影响 and neuContamLow是 1n Gauss分布的尾巴在Xn区间的部分
            else if(ineu == nNeus-1){
                neuContamLow[idir][ineu] = singleGaus[idir][ineu-2]->Integral(mNeuZDCLow[idir][ineu], mNeuZDCHi[idir][ineu]) / binWidth / hZDC[idir]->Integral(binLow, binHi);
                neuContamHi[idir][ineu]  = 0;
            }
            //当ZDC判断中子数为1n时 0n部分的分布不会污染1n部分 and neuContamHi是 Xn Gauss分布的尾巴在1n区间的部分
            else{
                if(ineu == 1) neuContamLow[idir][ineu] = 0;
                else          neuContamLow[idir][ineu] = singleGaus[idir][ineu-2]->Integral(mNeuZDCLow[idir][ineu], mNeuZDCHi[idir][ineu]) / multiGaus[idir]->Integral(mNeuZDCLow[idir][ineu], mNeuZDCHi[idir][ineu]);

                neuContamHi[idir][ineu] = singleGaus[idir][ineu]->Integral(mNeuZDCLow[idir][ineu], mNeuZDCHi[idir][ineu]) / multiGaus[idir]->Integral(mNeuZDCLow[idir][ineu], mNeuZDCHi[idir][ineu]);
            }
            //这里的纯度定义只关注于 ZDC判断中子数区间内 别的部分对其判断造成的影响 - ZDC判断的准确性
            neuPur[idir][ineu] = 1 - neuContamLow[idir][ineu] - neuContamHi[idir][ineu];
        }
    }
}
