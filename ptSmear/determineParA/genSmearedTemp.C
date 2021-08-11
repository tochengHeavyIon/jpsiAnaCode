#include "../../common/headers.h"
#include "../../common/function.C"
#include "../../common/funUtil.h"

const Double_t mTinyNum = 1.e-6;

void genSmearedTemp(Double_t mMaxPt = 0.15)
{
    TFile *fJpsiTemp  = TFile::Open("../makeSmear/smearedHistos/dimuonHistos.CohJpsi.root");
    TFile *fPsi2STemp = TFile::Open("../makeSmear/smearedHistos/dimuonHistos.CohPsi2S.root");
    TFile *fQEDTemp   = TFile::Open("../makeSmear/smearedHistos/dimuonHistos.LowMassGammaGamma.root");

    TString dirName  = "smearedTemp";
    system(Form("mkdir -p %s", dirName.Data()));
    system(Form("rm -rf %s/*", dirName.Data()));

    TH3D *hMvsPtvsRap_CohJpsi[nSmearScan];
    TH3D *hMvsPtvsRap_CohPsi2S[nSmearScan];
    TH3D *hMvsPtvsRap_QED[nSmearScan];

    TH1D *hMass_CohJpsi[nSmearScan], *hMass_CohPsi2S[nSmearScan], *hMass_QED[nSmearScan];
    TH1D *hMass_CohJpsi_Rap[nRapBins][nSmearScan], *hMass_CohPsi2S_Rap[nRapBins][nSmearScan], *hMass_QED_Rap[nRapBins][nSmearScan];
    for(Int_t iscan=0; iscan<nSmearScan; iscan++)
    {
        Double_t par0 = mInitPar0 + iscan*mSmearStep;

        hMvsPtvsRap_CohJpsi[iscan] = (TH3D *)fJpsiTemp->Get(Form("hMvsPtvsRap_Scan%d", iscan));
        hMvsPtvsRap_CohJpsi[iscan]->SetName(Form("hMvsPtvsRap_CohJpsi_Scan%d", iscan));

        hMvsPtvsRap_CohPsi2S[iscan] = (TH3D *)fPsi2STemp->Get(Form("hMvsPtvsRap_Scan%d", iscan));
        hMvsPtvsRap_CohPsi2S[iscan]->SetName(Form("hMvsPtvsRap_CohPsi2S_Scan%d", iscan));

        hMvsPtvsRap_QED[iscan] = (TH3D *)fQEDTemp->Get(Form("hMvsPtvsRap_Scan%d", iscan));
        hMvsPtvsRap_QED[iscan]->SetName(Form("hMvsPtvsRap_QED_Scan%d", iscan));

        for(Int_t irap=0; irap<nRapBins; irap++)
        {
            Int_t rapBinLow = hMvsPtvsRap_CohJpsi[iscan]->GetXaxis()->FindBin(mRapLow[irap] + mTinyNum);
            Int_t rapBinHi  = hMvsPtvsRap_CohJpsi[iscan]->GetXaxis()->FindBin(mRapHi[irap] - mTinyNum);
            Int_t ptBinLow  = 1;
            Int_t ptBinHi   = hMvsPtvsRap_CohJpsi[iscan]->GetYaxis()->FindBin(mMaxPt - mTinyNum);

            hMass_CohJpsi_Rap[irap][iscan]  = (TH1D *)hMvsPtvsRap_CohJpsi[iscan]->ProjectionZ(Form("hMass_CohJpsi_Rap%d_Scan%d", irap, iscan), rapBinLow, rapBinHi, ptBinLow, ptBinHi);
            hMass_CohJpsi_Rap[irap][iscan]->SetTitle(Form("par0 = %5.5f, %1.1f < y < %1.1f", par0, mRapLow[irap], mRapHi[irap]));

            hMass_CohPsi2S_Rap[irap][iscan] = (TH1D *)hMvsPtvsRap_CohPsi2S[iscan]->ProjectionZ(Form("hMass_CohPsi2S_Rap%d_Scan%d", irap, iscan), rapBinLow, rapBinHi, ptBinLow, ptBinHi);
            hMass_CohPsi2S_Rap[irap][iscan]->SetTitle(Form("par0 = %5.5f, %1.1f < y < %1.1f", par0, mRapLow[irap], mRapHi[irap]));

            hMass_QED_Rap[irap][iscan]      = (TH1D *)hMvsPtvsRap_QED[iscan]->ProjectionZ(Form("hMass_QED_Rap%d_Scan%d", irap, iscan), rapBinLow, rapBinHi, ptBinLow, ptBinHi);
            hMass_QED_Rap[irap][iscan]->SetTitle(Form("par0 = %5.5f, %1.1f < y < %1.1f", par0, mRapLow[irap], mRapHi[irap]));

            if(irap==0){
                hMass_CohJpsi[iscan]  = (TH1D *)hMass_CohJpsi_Rap[irap][iscan]->Clone(Form("hMass_CohJpsi_Scan%d", iscan));
                hMass_CohJpsi[iscan]->SetTitle(Form("par0 = %5.5f, %1.1f < y < %1.1f", par0, mRapLow[nRapBins/2], mRapHi[nRapBins-1]));

                hMass_CohPsi2S[iscan] = (TH1D *)hMass_CohPsi2S_Rap[irap][iscan]->Clone(Form("hMass_CohPsi2S_Scan%d", iscan));
                hMass_CohPsi2S[iscan]->SetTitle(Form("par0 = %5.5f, %1.1f < y < %1.1f", par0, mRapLow[nRapBins/2], mRapHi[nRapBins-1]));

                hMass_QED[iscan]      = (TH1D *)hMass_QED_Rap[irap][iscan]->Clone(Form("hMass_QED_Scan%d", iscan));
                hMass_QED[iscan]->SetTitle(Form("par0 = %5.5f, %1.1f < y < %1.1f", par0, mRapLow[nRapBins/2], mRapHi[nRapBins-1]));
            }
            else{
                hMass_CohJpsi[iscan]->Add(hMass_CohJpsi_Rap[irap][iscan]);
                hMass_CohPsi2S[iscan]->Add(hMass_CohPsi2S_Rap[irap][iscan]);
                hMass_QED[iscan]->Add(hMass_QED_Rap[irap][iscan]);
            }
        }
    }

    TFile *fOut = new TFile(Form("%s/smearedTemp.root", dirName.Data()), "recreate");
    fOut->cd();
    for(Int_t iscan=0; iscan<nSmearScan; iscan++){
        hMass_CohJpsi[iscan]->Write();
        hMass_CohPsi2S[iscan]->Write();
        hMass_QED[iscan]->Write();
        for(Int_t irap=0; irap<nRapBins; irap++){
            hMass_CohJpsi_Rap[irap][iscan]->Write();
            hMass_CohPsi2S_Rap[irap][iscan]->Write();
            hMass_QED_Rap[irap][iscan]->Write();
        }
    }

    cout << "End of program !" << endl;
}
