#include "../common/headers.h"
#include "../common/function.C"
#include "../common/funUtil.h"

const Bool_t  mStorePDF = kTRUE;

const Double_t mTinyNum = 1.e-6;

Int_t    mGenMarker = 20;
Int_t    mGenColor = 1;
Int_t    mGenWidth = 2;
Double_t mGenSize = 0.8;

void checkEffAndTemp_NueMult()
{
    gStyle->SetOptFit(1111);
    Int_t mFont = 42;

    const Int_t nSpecs = 4;
    TString specName[nSpecs]  = {"CohJpsi", "CohJpsi_0n0n", "CohJpsi_0nXn", "CohJpsi_XnXn", };
    TString specTitle[nSpecs] = {"Coh. J/#psi", "Coh. J/#psi (0n0n)", "Coh. J/#psi (0nXn)", "Coh. J/#psi (XnXn)"};

    TFile* f[nSpecs];
    TH3D *hMvsPtvsRap_Gen[nSpecs];
    TH3D *hMvsPtvsRap_woEvtSel[nSpecs];
    TH3D *hMvsPtvsRap[nSpecs];
    TH1D *hRap_Gen[nSpecs];
    TH1D *hRap_woEvtSel[nSpecs];
    TH1D *hRap[nSpecs];

    TH1D *hEffvsRap_woEvtSel[nSpecs];
    TH1D *hEffvsRap[nSpecs];
    TH1D *hEffRatiovsRap[nSpecs];
    TH1D *hEvtSelEffvsRap[nSpecs];

    TH1D *hMass[nSpecs];
    TH1D *hPt[nSpecs];
    TH1D *hMass_Rap[nSpecs][nDiffRapBins];
    TH1D *hPt_Rap[nSpecs][nDiffRapBins];
    for(Int_t i=0; i<nSpecs; i++){
        f[i] = TFile::Open(Form("mcHistos/dimuonHistos.%s.root", specName[i].Data()));
        hMvsPtvsRap_Gen[i]      = (TH3D *)f[i]->Get("hMvsPtvsRap_Gen");
        hMvsPtvsRap_woEvtSel[i] = (TH3D *)f[i]->Get("hMvsPtvsRap_woEvtSel");
        hMvsPtvsRap[i]          = (TH3D *)f[i]->Get("hMvsPtvsRap");

        Int_t massBinLow = 1;
        Int_t massBinHi  = hMvsPtvsRap_Gen[i]->GetNbinsZ();
        hRap_Gen[i]      = (TH1D *)hMvsPtvsRap_Gen[i]->ProjectionX(Form("hRap_Gen_%s", specName[i].Data()), 0, -1, massBinLow, massBinHi);
        hRap_woEvtSel[i] = (TH1D *)hMvsPtvsRap_woEvtSel[i]->ProjectionX(Form("hRap_woEvtSel_%s", specName[i].Data()), 0, -1, massBinLow, massBinHi);
        hRap[i]          = (TH1D *)hMvsPtvsRap[i]->ProjectionX(Form("hRap_%s", specName[i].Data()), 0, -1, massBinLow, massBinHi);

        hEffvsRap_woEvtSel[i] = (TH1D *)hRap_woEvtSel[i]->Clone(Form("hEffvsRap_woEvtSel_%s", specName[i].Data()));
        hEffvsRap_woEvtSel[i]->Divide(hRap_woEvtSel[i], hRap_Gen[i], 1, 1, "B");

        hEffvsRap[i] = (TH1D *)hRap[i]->Clone(Form("hRap_%s", specName[i].Data()));
        hEffvsRap[i]->Divide(hRap[i], hRap_Gen[i], 1, 1, "B");

        hEffRatiovsRap[i] = (TH1D *)hEffvsRap[i]->Clone(Form("hEffRatiovsRap_%s", specName[i].Data()));
        hEffRatiovsRap[i]->Divide(hEffvsRap[i], hEffvsRap[0], 1, 1, "B");

        hEvtSelEffvsRap[i] = (TH1D *)hRap[i]->Clone(Form("hEvtSelEffvsRap_%s", specName[i].Data())); 
        hEvtSelEffvsRap[i]->Divide(hRap[i], hRap_woEvtSel[i], 1, 1, "B");

        for(Int_t irap=0; irap<nDiffRapBins; irap++){
            Int_t rapBinLow  = hMvsPtvsRap[i]->GetXaxis()->FindBin(mDiffRapLow[irap] + mTinyNum);
            Int_t rapBinHi   = hMvsPtvsRap[i]->GetXaxis()->FindBin(mDiffRapHi[irap] - mTinyNum);
            hPt_Rap[i][irap]   = (TH1D *)hMvsPtvsRap[i]->ProjectionY(Form("hPt_%s_RapBin%d", specName[i].Data(), irap), rapBinLow, rapBinHi, massBinLow, massBinHi);
            hMass_Rap[i][irap] = (TH1D *)hMvsPtvsRap[i]->ProjectionZ(Form("hMass_%s_RapBin%d", specName[i].Data(), irap), rapBinLow, rapBinHi, 0, -1);

            if(irap==0) hMass[i] = (TH1D *)hMass_Rap[i][irap]->Clone(Form("hMass_%s", specName[i].Data()));
            else        hMass[i]->Add(hMass_Rap[i][irap]);

            if(irap==0) hPt[i] = (TH1D *)hPt_Rap[i][irap]->Clone(Form("hPt_%s", specName[i].Data()));
            else        hPt[i]->Add(hPt_Rap[i][irap]);
        }
    }

    TString dir = "effAndTemp_NeuMult";
    system(Form("mkdir -p %s", dir.Data()));
    system(Form("rm -rf %s/*", dir.Data()));

    TCanvas* c1 = new TCanvas("c1", "c1", 0, 0, 800, 600);

    TCanvas* c2 = new TCanvas("c2", "c2", 0, 0, 800, 300);
    c2->Divide(2, 1);

    TCanvas* c3 = new TCanvas("c3", "c3", 0, 0, 800, 600);
    c3->Divide(2, 2);

    Double_t xPos = 0.15;
    Double_t yPos = 0.84;

    TLegend* leg1 = new TLegend(0.36, 0.65, 0.6, 0.86);
    leg1->SetFillStyle(0);
    leg1->SetTextSize(0.06);

    TLegend* leg2 = new TLegend(0.36, 0.72, 0.66, 0.86);
    leg1->SetFillStyle(0);
    leg1->SetTextSize(0.06);

    for(Int_t i=0; i<nSpecs; i++){
        c2->cd(1);
        gPad->SetLogy(1);
        setHisto(hRap_Gen[i], 20, 0.8, 1, 1, 2);
        setHisto(hRap_woEvtSel[i], 24, 0.8, 2, 2, 2);
        setHisto(hRap[i], 25, 0.8, 4, 4, 2);
        hRap_Gen[i]->GetYaxis()->SetTitle("Entries");
        hRap_Gen[i]->SetMinimum(1);
        hRap_Gen[i]->Draw("p");
        hRap_woEvtSel[i]->Draw("psame");
        hRap[i]->Draw("psame");
        if(i==0){
            leg1->AddEntry(hRap_Gen[i], "GEN", "pl");
            leg1->AddEntry(hRap_woEvtSel[i], "#varepsilon_{reco}#times#varepsilon_{trig}", "pl");
            leg1->AddEntry(hRap[i], "#varepsilon_{reco}#times#varepsilon_{trig}#times#varepsilon_{evtSel}", "pl");
        }
        leg1->Draw("same");
        drawLatex(0.38, 0.95, Form("%s", specTitle[i].Data()), mFont, 0.06, 1);

        c2->cd(2);
        gPad->SetLogy(1);
        setHisto(hEffvsRap_woEvtSel[i], 24, 0.8, 2, 2, 2);
        setHisto(hEffvsRap[i], 25, 0.8, 4, 4, 2);
        hEffvsRap_woEvtSel[i]->GetYaxis()->SetTitle("Efficiency");
        hEffvsRap_woEvtSel[i]->SetMaximum(2);
        hEffvsRap_woEvtSel[i]->Draw("p");
        hEffvsRap[i]->Draw("psame");
        if(i==0){
            leg2->AddEntry(hEffvsRap_woEvtSel[i], "#varepsilon_{reco}#times#varepsilon_{trig}", "pl");
            leg2->AddEntry(hEffvsRap[i], "#varepsilon_{reco}#times#varepsilon_{trig}#times#varepsilon_{evtSel}", "pl");
        }
        leg2->Draw("same");
        drawLatex(0.38, 0.95, Form("%s", specTitle[i].Data()), mFont, 0.06, 1);

        if(mStorePDF) c2->SaveAs(Form("%s/EffvsRap_%s.pdf", dir.Data(), specName[i].Data()));
        c2->SaveAs(Form("%s/EffvsRap_%s.png", dir.Data(), specName[i].Data()));
    }

    Int_t    mMStyle[nSpecs] = {20, 24, 25, 27};
    Double_t mMSize[nSpecs] = {1.0, 1.0, 1.0, 1.5};
    Int_t    mColor[nSpecs]  = {1, 2, 4, 6};

    setLegend(leg1, 0.36, 0.62, 0.6, 0.86, 0.055);
    setLegend(leg2, 0.36, 0.68, 0.6, 0.86, 0.055);

    c2->cd(1);
    gPad->SetLogy(1);
    for(Int_t i=0; i<nSpecs; i++){
        setHisto(hEffvsRap[i], mMStyle[i], mMSize[i], mColor[i], mColor[i]);
        if(i==0){
            hEffvsRap[i]->GetXaxis()->SetRangeUser(-2.5, 2.5);
            hEffvsRap[i]->GetYaxis()->SetRangeUser(1.e-5, 1);
            hEffvsRap[i]->GetYaxis()->SetTitle("#varepsilon_{reco}#times#varepsilon_{trig}#times#varepsilon_{evtSel}");
            hEffvsRap[i]->Draw("p");
        }
        else{
            hEffvsRap[i]->Draw("psame");
        }
        leg1->AddEntry(hEffvsRap[i], specTitle[i].Data(), "pl");
    }
    leg1->Draw("same");
    c2->cd(2);
    gPad->SetLogy(0);
    for(Int_t i=1; i<nSpecs; i++){
        setHisto(hEffRatiovsRap[i], mMStyle[i], mMSize[i], mColor[i], mColor[i]);
        if(i==1){
            hEffRatiovsRap[i]->GetXaxis()->SetRangeUser(-2.5, 2.5);
            hEffRatiovsRap[i]->GetYaxis()->SetTitle("Eff. Ratio to Coh. J/#psi");
            hEffRatiovsRap[i]->GetYaxis()->SetRangeUser(0.8, 1.2);
            hEffRatiovsRap[i]->Draw("p");
            drawLine(-2.5, 1, 2.5, 1, 1, 2, 2);
            hEffRatiovsRap[i]->Draw("psame");
        }
        else{
            hEffRatiovsRap[i]->Draw("psame");
        }
    }
    if(mStorePDF) c2->SaveAs(Form("%s/EffComparisonvsRap.pdf", dir.Data()));
    c2->SaveAs(Form("%s/EffComparisonvsRap.png", dir.Data()));

    setLegend(leg1, 0.16, 0.16, 0.4, 0.4, 0.055);

    TH1D *hMassRatio_Rap[nSpecs][nDiffRapBins];
    TH1D *hPtRatio_Rap[nSpecs][nDiffRapBins];
    for(Int_t irap=0; irap<nDiffRapBins; irap++){
        gPad->SetLogy(1);
        for(Int_t i=0; i<nSpecs; i++){
            c3->cd(1);
            gPad->SetLogy(1);
            setHisto(hMass_Rap[i][irap], mMStyle[i], 0.5, mColor[i], mColor[i]);
            hMass_Rap[i][irap]->Scale(1./hMass_Rap[i][irap]->GetEntries());
            //hMass_Rap[i][irap]->GetYaxis()->SetRangeUser(0.5, hMass_Rap[i][irap]->GetMaximum()*2);
            hMass_Rap[i][irap]->GetXaxis()->SetTitleOffset(0.9);
            hMass_Rap[i][irap]->GetYaxis()->SetTitle("a.u.");
            if(i==0){
                hMass_Rap[i][irap]->GetXaxis()->SetRangeUser(2.7, 3.5);
                hMass_Rap[i][irap]->Draw("p");
            }
            else     hMass_Rap[i][irap]->Draw("psame");
            drawLatex(0.39, 0.95, Form("%1.1f < y < %1.1f", mDiffRapLow[irap], mDiffRapHi[irap]), mFont, 0.06, 1);

            hMassRatio_Rap[i][irap] = (TH1D *)hMass_Rap[i][irap]->Clone(Form("hMassRatio_%s_Rap%d", specName[i].Data(), irap));
            hMassRatio_Rap[i][irap]->Divide(hMass_Rap[i][irap], hMass_Rap[0][irap], 1, 1, "B");

            c3->cd(2);
            gPad->SetLogy(0);
            setHisto(hMassRatio_Rap[i][irap], mMStyle[i], 0.5, mColor[i], mColor[i]);
            hMassRatio_Rap[i][irap]->GetYaxis()->SetRangeUser(0.5, 1.5);
            hMassRatio_Rap[i][irap]->GetXaxis()->SetTitleOffset(0.9);
            hMassRatio_Rap[i][irap]->GetYaxis()->SetTitle("Temp. Ratio to inclu. J/#psi");
            if(i==1){
                hMassRatio_Rap[i][irap]->GetXaxis()->SetRangeUser(2.7, 3.5);
                hMassRatio_Rap[i][irap]->Draw("p");
                drawLine(2.7, 1, 3.5, 1, 1, 2, 2);
                hMassRatio_Rap[i][irap]->Draw("psame");
            }
            else if(i>1) hMassRatio_Rap[i][irap]->Draw("psame");
            drawLatex(0.39, 0.95, Form("%1.1f < y < %1.1f", mDiffRapLow[irap], mDiffRapHi[irap]), mFont, 0.06, 1);

            c3->cd(3);
            gPad->SetLogy(1);
            setHisto(hPt_Rap[i][irap], mMStyle[i], 0.5, mColor[i], mColor[i]);
            hPt_Rap[i][irap]->Scale(1./hPt_Rap[i][irap]->GetEntries());
            //hPt_Rap[i][irap]->GetYaxis()->SetRangeUser(0.5, hPt_Rap[i][irap]->GetMaximum()*2);
            hPt_Rap[i][irap]->GetXaxis()->SetTitleOffset(0.9);
            hPt_Rap[i][irap]->GetYaxis()->SetTitle("a.u.");
            if(i==0){
                hPt_Rap[i][irap]->GetXaxis()->SetRangeUser(0, 0.4);
                hPt_Rap[i][irap]->Draw("p");
            }
            else     hPt_Rap[i][irap]->Draw("psame");
            drawLatex(0.39, 0.95, Form("%1.1f < y < %1.1f", mDiffRapLow[irap], mDiffRapHi[irap]), mFont, 0.06, 1);
            if(irap==0){
                leg1->AddEntry(hPt_Rap[i][irap], specTitle[i].Data(), "pl");
            }
            leg1->Draw("same");

            hPtRatio_Rap[i][irap] = (TH1D *)hPt_Rap[i][irap]->Clone(Form("hPtRatio_%s_Rap%d", specName[i].Data(), irap));
            hPtRatio_Rap[i][irap]->Divide(hPt_Rap[i][irap], hPt_Rap[0][irap], 1, 1, "B");

            c3->cd(4);
            gPad->SetLogy(0);
            setHisto(hPtRatio_Rap[i][irap], mMStyle[i], 0.5, mColor[i], mColor[i]);
            hPtRatio_Rap[i][irap]->GetYaxis()->SetRangeUser(0.5, 1.5);
            hPtRatio_Rap[i][irap]->GetXaxis()->SetTitleOffset(0.9);
            hPtRatio_Rap[i][irap]->GetYaxis()->SetTitle("Temp. Ratio to inclu. J/#psi");
            if(i==1){
                hPtRatio_Rap[i][irap]->GetXaxis()->SetRangeUser(0, 0.4);
                hPtRatio_Rap[i][irap]->Draw("p");
                drawLine(0, 1, 0.4, 1, 1, 2, 2);
                hPtRatio_Rap[i][irap]->Draw("psame");
            }
            else if(i>1) hPtRatio_Rap[i][irap]->Draw("psame");
            drawLatex(0.39, 0.95, Form("%1.1f < y < %1.1f", mDiffRapLow[irap], mDiffRapHi[irap]), mFont, 0.06, 1);
        }

        if(mStorePDF) c3->SaveAs(Form("%s/TempComparison_RapBin%d.pdf", dir.Data(), irap));
        c3->SaveAs(Form("%s/TempComparison_RapBin%d.png", dir.Data(), irap));
    }

    cout << "End of program !" << endl;
}
