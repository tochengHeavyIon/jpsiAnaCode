#include "../../common/headers.h"
#include "../../common/function.C"

void plotPtRes(){
    Int_t mFont = 42;

    gStyle->SetOptFit(1111);

    TF1 *sg = new TF1("sg","gaus", -0.5, 0.5);
    sg->SetParNames("N", "#mu", "#sigma");
    setFun(sg, 2, 2, 1);

    TFile *fMc   = new TFile("../../simulation/mcHistos/dimuonHistos.LowMassGammaGamma.root");
    TH2D *hPtResVsGenPt = (TH2D *)fMc->Get("hPtResvsGenPt");

    //hPtResVsGenPt->Sumw2();
    hPtResVsGenPt->GetYaxis()->SetTitle("(p_{T}^{RECO}-p_{T}^{GEN})/p_{T}^{GEN}");
    hPtResVsGenPt->RebinX(2);

    TString dir = "ptRes";
    system(Form("mkdir -p %s", dir.Data()));
    system(Form("rm -rf %s/*", dir.Data()));

    TCanvas *c1 = new TCanvas("c1", "c1",0,0,800,600);
    TPDF *ps = new TPDF(Form("%s/ptResFitDetails.pdf", dir.Data()),111);
    ps->Off();

    Int_t nRaws = 5;
    Int_t nColumns = 5;
    Int_t nPads = nRaws*nColumns;
    TCanvas *c3 = new TCanvas("c3","c3",1200,900);
    c3->Divide(nColumns, nRaws);
    for(Int_t ipad=0; ipad<nPads; ipad++){
        c3->cd(ipad+1);
        setPad(0.12, 0.06, 0.09, 0.16);
    }

    const Int_t maxPts = 500;
    const Int_t statTh = 500;

    Double_t ptCenter[maxPts];
    Double_t ptErr[maxPts];
    Double_t mean[maxPts];
    Double_t meanErr[maxPts];
    Double_t sigma[maxPts];
    Double_t sigmaErr[maxPts];
    Double_t sigmaRatio[maxPts];
    Double_t sigmaRatioErr[maxPts];

    Int_t nPts = 0;
    TH1D *hPtRes[maxPts];
    for(Int_t i=0; i<hPtResVsGenPt->GetNbinsX(); i++){
        hPtRes[i] = (TH1D *)hPtResVsGenPt->ProjectionY(Form("hPtRes_ptBin%d", i), i+1, i+1);

        if(hPtRes[i]->GetEntries()<statTh) {
            //cout<<hPtResVsGenPt->GetXaxis()->GetBinCenter(i+1)<<"     "<<hPtRes[i]->GetEntries()<<endl;
            continue;
        }

        ptCenter[nPts] = hPtResVsGenPt->GetXaxis()->GetBinCenter(i+1);
        ptErr[nPts]    = hPtResVsGenPt->GetXaxis()->GetBinWidth(i+1)/2;

        if(ptCenter[nPts]>3.2)       hPtRes[i]->RebinX(5);
        else if(ptCenter[nPts]>2)    hPtRes[i]->RebinX(2);
        else if(ptCenter[nPts]<0.95) hPtRes[i]->RebinX(4);

        c3->cd(nPts%nPads+1);
        hPtRes[i]->GetXaxis()->SetTitleOffset(1.05);
        hPtRes[i]->GetYaxis()->SetTitle("Entries");
        setHisto(hPtRes[i],20,0.6,1,1);

        Double_t nSimgaLow = 2, nSigmaHi = 2;
        Double_t xmean = hPtRes[i]->GetMean();
        Double_t xrms  = hPtRes[i]->GetRMS();
        hPtRes[i]->SetAxisRange(-0.1,0.1,"X");
        hPtRes[i]->Fit(sg,"QR","",xmean-nSimgaLow*xrms,xmean+nSigmaHi*xrms);
        xmean = sg->GetParameter(1);
        xrms  = sg->GetParameter(2);
        hPtRes[i]->Fit(sg,"QR","",xmean-nSimgaLow*xrms,xmean+nSigmaHi*xrms);
        hPtRes[i]->Draw("pe");
        sg->SetRange(xmean-nSimgaLow*xrms,xmean+nSigmaHi*xrms);
        sg->DrawClone("same");

        mean[nPts]     = sg->GetParameter(1);
        meanErr[nPts]  = sg->GetParError(1);
        sigma[nPts]    = sg->GetParameter(2);
        sigmaErr[nPts] = sg->GetParError(2);

        drawLatex(0.36, 0.955, Form("%1.2f<p_{T}<%1.2f GeV", ptCenter[nPts]-ptErr[nPts], ptCenter[nPts]+ptErr[nPts]), mFont, 0.06, 1);

        if(nPts%nPads == nPads-1) pdfAction(c3, ps);

        nPts++;
    }
    if(nPts%nPads!=0) pdfAction(c3, ps, kTRUE);

    TGraphErrors *grMean   = new TGraphErrors(nPts,ptCenter,mean,0,meanErr);
    TGraphErrors *grPtRes = new TGraphErrors(nPts,ptCenter,sigma,0,sigmaErr);

    TH2D *dd = (TH2D *)histo("dd",0.5,4.5,-0.02,0.08,"muon p_{T}^{GEN} (GeV)","");
    dd->GetXaxis()->SetNdivisions(210);
    dd->GetYaxis()->SetNdivisions(210);
    dd->GetXaxis()->SetTitleOffset(1.0);
    dd->GetYaxis()->SetTitleOffset(1.0);

    TLegend *leg = new TLegend(0.3, 0.74, 0.5, 0.86);
    leg->SetFillStyle(0);
    leg->SetFillColor(10);
    leg->SetTextFont(mFont);
    leg->SetTextSize(0.06);

    TF1 *funPtMeanShift = new TF1("funPtMeanShift","[0]/pow(x, [1])+[2]",0,5);
    funPtMeanShift->SetParameter(1, 1.8);
    setFun(funPtMeanShift,4,2,1);

    TF1 *funRawPtRes = new TF1("funRawPtRes","sqrt([0]*[0]/x/x+[1]*[1])",0,5);
    setFun(funRawPtRes,2,2,1);

    c1->cd();
    setPad(0.12, 0.06, 0.08, 0.14);
    dd->GetYaxis()->SetRangeUser(-5e-3, 10e-3);
    dd->GetYaxis()->SetTitle("#mu of (p_{T}^{RECO}-p_{T}^{GEN})/p_{T}^{GEN}");
    dd->Draw("c");
    setGraph(grMean,20,0.8,1,1);
    grMean->Fit(funPtMeanShift,"R","",1,4);
    grMean->Draw("pzsame");
    c1->SaveAs(Form("%s/ptMeanShift.pdf",dir.Data()));
    c1->SaveAs(Form("%s/ptMeanShift.png",dir.Data()));

    c1->cd();
    setPad(0.12, 0.06, 0.08, 0.14);
    dd->GetYaxis()->SetRangeUser(10e-3, 30e-3);
    dd->GetYaxis()->SetTitle("#sigma of (p_{T}^{RECO}-p_{T}^{GEN})/p_{T}^{GEN}");
    dd->Draw("c");
    setGraph(grPtRes,20,0.8,1,1);
    grPtRes->Fit(funRawPtRes,"R","",1,4);
    grPtRes->Draw("pzsame");
    setLegend(leg,0.16,0.74,0.4,0.86,0.06);
    c1->SaveAs(Form("%s/rawPtRes.pdf",dir.Data()));
    c1->SaveAs(Form("%s/rawPtRes.png",dir.Data()));

    TFile *fOut = new TFile(Form("%s/rawPtRes.root",dir.Data()),"recreate");
    hPtResVsGenPt->SetTitle("muon (pT_RECO-pT_GEN)/pT_GEN distribution from simulation");
    hPtResVsGenPt->Write();
    funPtMeanShift->SetTitle("muon pT mean shift from simulation ");
    funPtMeanShift->Write();
    funRawPtRes->SetTitle("muon pT resolution from simulation");
    funRawPtRes->Write();
    fOut->Close();

    cout<<"End of program !"<<endl;
    return;
}
