#include "../common/headers.h"
#include "../common/function.C"
#include "../common/constants.h"

const Double_t tinyOffset = 1.e-6;

TH1D *hZDC[nDirs];

const Int_t nGaus = 3;
TF1 *multiGaus[nDirs];
TF1 *singleGaus[nDirs][nGaus];

void plotZDC(Bool_t isPaper=kTRUE, Bool_t drawNeuRange = 1)
{
    gStyle->SetPalette(1);

    Int_t    mFont = 42;
    gStyle->SetTextFont(mFont);
    gStyle->SetLegendFont(mFont);	
    gStyle->SetLabelFont(mFont,"xyz");
    gStyle->SetTitleFont(mFont,"xyz");

    gStyle->SetOptFit(111);

    TFile* f = TFile::Open("../anaData/jpsiHistos/rawSig.root");

    TString dir = Form("zdcPlots");
    system(Form("mkdir -p %s", dir.Data()));
    system(Form("rm -rf %s/*", dir.Data()));

    TCanvas* c1 = new TCanvas("c1", "c1", 1200, 900);
    c1->Divide(2, 2);
    for(Int_t ipad=0; ipad<4; ipad++){
        c1->cd(ipad+1);
        if(ipad==0) setPad(0.12, 0.12, 0.08, 0.13);
        else        setPad(0.12, 0.08, 0.08, 0.13);
    }

    TH2D* hZDCMinusvsZDCPlus = (TH2D*)f->Get("hZDCMinusvsZDCPlus"); // In principle, the neutron peak should not be affected by UPC or hadronic collisions

    TString triplegaus = "[0]*TMath::Gaus(x,[1],[2],1)"
        "+ [3]*TMath::Gaus(x,[4],[5],1)"
        "+ [6]*TMath::Gaus(x,[7],[5]*[8],1)";

    for(Int_t idir=0; idir<nDirs; idir++){
        multiGaus[idir] = new TF1(Form("multiGaus_%s", mDir[idir].Data()), triplegaus.Data(), 0, 5.e4);
        multiGaus[idir]->SetParNames("N1", "#mu1", "#sigma1", "N2", "#mu2", "#sigma2", "N3", "#mu3", "scale w.r.t. #sigma2");
        multiGaus[idir]->SetParLimits(8, 1, 2);
        setFun(multiGaus[idir], kBlue-4, 2);
    }

    Double_t pars[nDirs][9] = {
        { 2e7, 7.5e3, 1.6e3, 6e6, 15e3, 2.3e3, 1e6, 22.5e3, 1.4 },
        { 2e7, 12e3, 3.0e3, 6e6, 25e3, 4.2e3, 1e6, 35e3, 1.4 }
    };

    const Int_t mGausColor[3] = {kRed-4, kGreen+2, kViolet-3};
    for(Int_t idir=0; idir<nDirs; idir++){
        for (Int_t i = 0; i < nGaus; i++) {
            singleGaus[idir][i] = new TF1(Form("singleGaus_%s%d", mDir[idir].Data(), i), "[0]*TMath::Gaus(x,[1],[2],1)", 0, 5e4);
            setFun(singleGaus[idir][i], mGausColor[i], 2, 2);
        }
    }

    hZDC[0] = (TH1D*)hZDCMinusvsZDCPlus->ProjectionX(Form("hZDC%s", mDir[0].Data()));
    hZDC[1] = (TH1D*)hZDCMinusvsZDCPlus->ProjectionY(Form("hZDC%s", mDir[1].Data()));

    const Double_t arrY = 3e3;
    const Double_t arrSize = 7e-3;
    const Int_t    arrColor = kBlue - 4;
    const Int_t    arrWidth = 2;

    Double_t mTextSize = 0.07;

    c1->cd(1);
    gPad->SetLogz(1);
    hZDCMinusvsZDCPlus->Rebin2D(2, 2);
    hZDCMinusvsZDCPlus->GetXaxis()->SetTitleFont(mFont);
    hZDCMinusvsZDCPlus->GetYaxis()->SetTitleFont(mFont);
    hZDCMinusvsZDCPlus->GetXaxis()->SetTitle("ZDC_{Plus} (a.u.)");
    hZDCMinusvsZDCPlus->GetYaxis()->SetTitle("ZDC_{Minus} (a.u.)");
    hZDCMinusvsZDCPlus->GetXaxis()->SetLabelSize(0.055);
    hZDCMinusvsZDCPlus->GetXaxis()->SetTitleSize(0.065);
    hZDCMinusvsZDCPlus->GetXaxis()->SetTitleOffset(0.9);
    hZDCMinusvsZDCPlus->GetYaxis()->SetLabelSize(0.055);
    hZDCMinusvsZDCPlus->GetYaxis()->SetTitleSize(0.065);
    hZDCMinusvsZDCPlus->GetYaxis()->SetTitleOffset(0.9);
    hZDCMinusvsZDCPlus->GetXaxis()->SetRangeUser(0, 3e4);
    hZDCMinusvsZDCPlus->GetYaxis()->SetRangeUser(0, 3e4);
    hZDCMinusvsZDCPlus->GetZaxis()->SetLabelSize(0.055);
    hZDCMinusvsZDCPlus->GetZaxis()->SetLabelOffset(0.002);
    hZDCMinusvsZDCPlus->Draw("colz");
    gPad->Update();
    TPaletteAxis *palette = (TPaletteAxis*)hZDCMinusvsZDCPlus->GetListOfFunctions()->FindObject("palette");
    palette->SetX1NDC(0.885);
    palette->SetX2NDC(0.935);
    palette->SetY1NDC(0.21);
    palette->SetY2NDC(0.92);
    drawLatex(0.56, 0.94, "PbPb 5.02 TeV", mFont, mTextSize, 1);
    drawLatex(0.20, 0.94, "CMS", 62, 0.085, 1);
    if(!isPaper) drawLatex(0.31, 0.94, "Preliminary", 52, mTextSize, 1);

    TGraphErrors *grZdcRatio[nDirs];
    TGraph       *grZdcPull[nDirs];
    TH1D         *hZDCPull[nDirs];

    TLegend *leg = new TLegend(0.35, 0.3, 0.75, 0.76);
    leg->SetFillColor(10);
    leg->SetFillStyle(0);
    leg->SetTextFont(mFont);
    leg->SetTextSize(0.1);

    gStyle->SetStatX(0.92);
    gStyle->SetStatY(0.92);
    gStyle->SetStatH(0.17);
    gStyle->SetStatW(0.21);

    for(Int_t idir=0; idir<nDirs; idir++){
        c1->cd(2+idir);
        gPad->SetLogy(1);
        setHisto(hZDC[idir], 20, 0.8, 1, 1, 1);

        if(idir==0) hZDC[idir]->GetXaxis()->SetRangeUser(0, 0.35e5);
        else        hZDC[idir]->GetXaxis()->SetRangeUser(0, 0.5e5);
        hZDC[idir]->GetXaxis()->SetLabelSize(0.055);
        hZDC[idir]->GetXaxis()->SetTitleSize(0.065);
        hZDC[idir]->GetXaxis()->SetTitleOffset(0.9);
        hZDC[idir]->GetXaxis()->SetTitleFont(mFont);
        hZDC[idir]->GetYaxis()->SetLabelSize(0.055);
        hZDC[idir]->GetYaxis()->SetTitleSize(0.065);
        hZDC[idir]->GetYaxis()->SetTitleOffset(1.0);
        hZDC[idir]->GetYaxis()->SetTitleFont(mFont);
        hZDC[idir]->GetYaxis()->SetRangeUser(9.999, 1e5);
        multiGaus[idir]->SetParameters(pars[idir]);
        multiGaus[idir]->SetRange(mZdcFitLow[idir], mZdcFitHi[idir]);
        hZDC[idir]->Fit(multiGaus[idir], "R");
        hZDC[idir]->Draw("p");
        if(drawNeuRange){
            for(Int_t ineu=0; ineu<nNeus-1; ineu++){
                drawArrow(mNeuZDCLow[idir][ineu], arrY, mNeuZDCHi[idir][ineu], arrY, arrSize, "<|>", arrColor, arrWidth, 1);
                drawLine(mNeuZDCHi[idir][ineu], 10, mNeuZDCHi[idir][ineu], 2e4, arrColor, arrWidth, 5);
            }
            if(idir==0){
                drawArrow(mNeuZDCLow[idir][nNeus-1], arrY, mNeuZDCLow[idir][nNeus-1]+8e3, arrY, arrSize, "|>", arrColor, arrWidth, 1);
            }
            else{
                drawArrow(mNeuZDCLow[idir][nNeus-1], arrY, mNeuZDCLow[idir][nNeus-1]+8e3*5/3.5, arrY, arrSize, "|>", arrColor, arrWidth, 1);
            }

            drawLatex(0.16, 0.68, "0n", mFont, mTextSize, arrColor);
            drawLatex(0.28, 0.68, "Xn", mFont, mTextSize, arrColor);
        }

        hZDC[idir]->Draw("psame");
        hZDC[idir]->GetXaxis()->SetTitle(Form("ZDC_{%s} (a.u.)", mDir[idir].Data()));

        drawLatex(0.60, 0.94, "PbPb 5.02 TeV", mFont, mTextSize, 1);
        drawLatex(0.12, 0.94, "CMS", 62, 0.085, 1);
        if(!isPaper) drawLatex(0.24, 0.94, "Preliminary", 52, mTextSize, 1);
        hZDC[idir]->GetYaxis()->SetTitle("Entries");
        multiGaus[idir]->DrawClone("same");
        if(idir==0) leg->AddEntry(multiGaus[idir], "Total fit", "l");
        for (Int_t i = 0; i < nGaus; i++) {
            if(i < nGaus-1) singleGaus[idir][i]->SetParameters(multiGaus[idir]->GetParameter(3*i), multiGaus[idir]->GetParameter(3*i + 1), multiGaus[idir]->GetParameter(3*i + 2));
            else            singleGaus[idir][i]->SetParameters(multiGaus[idir]->GetParameter(3*i), multiGaus[idir]->GetParameter(3*i + 1), multiGaus[idir]->GetParameter(3*i + 2)*multiGaus[idir]->GetParameter(3*i - 1));

            singleGaus[idir][i]->DrawClone("same");

            if(idir==0) leg->AddEntry(singleGaus[idir][i], Form("%dn", i+1), "l");
        }

        Int_t binLow = hZDC[idir]->GetXaxis()->FindBin(mZdcFitLow[idir] + tinyOffset);
        Int_t binHi  = hZDC[idir]->GetXaxis()->FindBin(mZdcFitHi[idir] - tinyOffset);
        grZdcRatio[idir] = new TGraphErrors(binHi - binLow + 1);
        grZdcPull[idir]  = new TGraph(binHi - binLow + 1);
        hZDCPull[idir]   = new TH1D(Form("hZDCPull%d",idir),Form("; ZDC %s Pull; Entries", mDir[idir].Data()), 50, -10, 10);
        for(Int_t binIdx=binLow; binIdx<=binHi; binIdx++){
            Double_t x = hZDC[idir]->GetBinCenter(binIdx);
            Double_t y = hZDC[idir]->GetBinContent(binIdx);
            Double_t yErr = hZDC[idir]->GetBinError(binIdx);
            Double_t zdcFit = multiGaus[idir]->Eval(x);

            grZdcRatio[idir]->SetPoint(binIdx-binLow, x, y/zdcFit);
            grZdcRatio[idir]->SetPointError(binIdx-binLow, 0, yErr/zdcFit);

            grZdcPull[idir]->SetPoint(binIdx-binLow, x, (y-zdcFit)/yErr);
            hZDCPull[idir]->Fill((y-zdcFit)/yErr);
        }
    }

    c1->cd(4);
    leg->Draw("same");
    if(isPaper) c1->SaveAs(Form("%s/ZDC_Energy.pdf", dir.Data()));
    else        c1->SaveAs(Form("%s/ZDC_Energy_Preliminary.pdf", dir.Data()));

    TH2D *ddRatio = (TH2D *)histo("ddRatio", 0, 40e3, -10, 10, "", "Data/Fit");
    ddRatio->GetXaxis()->SetTitleSize(0.06);
    ddRatio->GetYaxis()->SetTitleSize(0.06);

    for(Int_t idir=0; idir<nDirs; idir++){
        c1->cd(idir*2+1);
        gPad->SetLogy(0);
        ddRatio->GetXaxis()->SetTitle(Form("ZDC %s", mDir[idir].Data()));
        ddRatio->GetXaxis()->SetRangeUser(mZdcFitLow[idir], mZdcFitHi[idir]);
        ddRatio->GetYaxis()->SetTitle("Pull");
        ddRatio->GetYaxis()->SetRangeUser(-6, 6);
        ddRatio->DrawClone("c");
        drawLine(mZdcFitLow[idir], 0, mZdcFitHi[idir], 0, 2, 3, 2);
        setGraph(grZdcPull[idir], 24, 0.6, 1, 1, 2);
        grZdcPull[idir]->Draw("pzsame");

        c1->cd(idir*2+2);
        gPad->SetLogy(0);
        setHisto(hZDCPull[idir], 20, 1.2, 1, 1, 2);
        hZDCPull[idir]->GetXaxis()->SetRangeUser(-6, 8);
        hZDCPull[idir]->Fit("gaus", "", "", -4, 4);
        hZDCPull[idir]->Draw("pe");
    }
    c1->SaveAs(Form("%s/zdcPull.pdf", dir.Data()));
    c1->SaveAs(Form("%s/zdcPull.png", dir.Data()));

    TFile* fOut = new TFile("neutronDecouple.root", "recreate");

    fOut->cd();
    hZDC[0]->Write();
    hZDC[1]->Write(); 

    multiGaus[0]->Write();
    multiGaus[1]->Write();
    for(Int_t iGaus=0; iGaus<nGaus; iGaus++){
       singleGaus[0][iGaus]->Write();
       singleGaus[1][iGaus]->Write();
    }
    fOut->Close();

    cout << "End of program !" << endl;
}
