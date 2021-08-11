#include "../../common/headers.h"
#include "../../common/funUtil.h"
#include "../../common/VertexCompositeTree.h"

TH1D *hnEvts;
TH1D *hDeltaR;
TH1D *hDeltaPt;
TH3D *hMvsPtvsRap_Scan[nSmearScan];

Bool_t   goodMuPair(VertexCompositeTree& evtTree, const int icand);
Bool_t   goodMuPair(const TLorentzVector posFourMom, const Bool_t isPosTrig, const TLorentzVector negFourMom, const Bool_t isNegTrig);

void bookHistos();
void writeHistos(TString fileName = "test");

void smearMcEvt(TString fileName = "CohJpsi")
{
    TH1::SetDefaultSumw2(kTRUE);

    std::string inputFile;
    if(fileName.EqualTo("LowMassGammaGamma"))     inputFile = "../../rootfiles/VertexCompositeTree_STARLIGHT_LowMassGammaGammaToMuMu_GenFilter_DiMuMC_20200906.root";
    else if(fileName.EqualTo("CohJpsi"))          inputFile = "../../rootfiles/VertexCompositeTree_STARLIGHT_CohJpsiToMuMu_GenFilter_DiMuMC_20200906.root";
    else if(fileName.EqualTo("CohPsi2S"))         inputFile = "../../rootfiles/VertexCompositeTree_STARLIGHT_CohPsi2SToMuMu_GenFilter_DiMuMC_20200906.root";
    else{
        cout<<"The filename string is wrong! Should be 'CohJpsi', 'CohPsi2S', OR 'LowMassGammaGamma' "<<endl;
        return;
    }

    const auto& csTreeDir = "dimucontana_mc";

    // Extract the tree
    VertexCompositeTree csTree;
    if (!csTree.GetTree(inputFile, csTreeDir)) {
        cout << "Invalid Correct-Sign tree!" << endl;
        return;
    }

    if(!init()) {
        cout<<"Initialization failed !"<<endl;
        return;
    }

    TF1  *funNewPtRes[nSmearScan];
    for(Int_t i=0; i<nSmearScan; i++){
        Double_t par0 = mInitPar0 + i*mSmearStep;
        funNewPtRes[i] = new TF1(Form("funNewPtRes_Scan%d", i), "sqrt([0]*[0]/x/x+[1]*[1])", 0, 5);
        funNewPtRes[i]->SetParameters(par0, funRawPtRes->GetParameter(1));
    }

    bookHistos();

    for (Long64_t jentry = 1; jentry < csTree.GetEntries(); jentry++) {
        if (jentry % (csTree.GetEntries() / 10) == 0)
            cout << "begin " << jentry << "th entry...." << endl;

        // Get the entry
        if (csTree.GetEntry(jentry) < 0) {
            cout << "Invalid correct-sign entry!" << endl;
            return;
        }

        hnEvts->Fill(0.5);

        // Select Event - require this event has a valid vertex and is not beam-halo event
        //if(!csTree.evtSel()[2] || !csTree.evtSel()[3]) continue;
        Bool_t goodVtx = (csTree.evtSel()[2] && csTree.evtSel()[3]);
        if(goodVtx) hnEvts->Fill(1.5);

        //evtSel[4-15]
        //[4]=0: HFPlusMaxTower < 3 GeV;  [4]=1: HFPlusMaxTower > 3 GeV
        //[5]=0: HFMinusMaxTower < 3 GeV;  [5]=1: HFMinusMaxTower > 3 GeV
        //[6] is for Plus & [7] is for Minus; Threshold = 4 GeV
        //[8] is for Plus & [9] is for Minus; Threshold = 5 GeV
        //[10] is for Plus & [11] is for Minus; Threshold = 6 GeV
        //[12] is for Plus & [13] is for Minus; Threshold = 7 GeV
        //[14] is for Plus & [15] is for Minus; Threshold = 8 GeV
        //[16] is for Plus (Th = 7.6 GeV) & [17] is for Minus (Th = 7.3 GeV); 

        //if(csTree.evtSel()[16] || csTree.evtSel()[17]) continue;
        Bool_t goodHFVeto = (!csTree.evtSel()[16] && !csTree.evtSel()[17]);
        if(goodVtx & goodHFVeto) hnEvts->Fill(2.5);

        Int_t nTrkHP = csTree.NtrkHP();

        //if(nTrkHP != 2) continue;
        Bool_t passEvtSel = goodVtx && goodHFVeto && (nTrkHP==2);

        if(passEvtSel) hnEvts->Fill(3.5);
        else           continue;

        // Loop over the correct-sign candidates
        for (UInt_t icand = 0; icand < csTree.candSize_gen(); icand++) {
            Double_t posPt_gen  = csTree.chargeD1_gen()[icand] > 0 ? csTree.pTD1_gen()[icand] : csTree.pTD2_gen()[icand]; 
            Double_t posEta_gen = csTree.chargeD1_gen()[icand] > 0 ? csTree.EtaD1_gen()[icand] : csTree.EtaD2_gen()[icand]; 
            Double_t posPhi_gen = csTree.chargeD1_gen()[icand] > 0 ? csTree.PhiD1_gen()[icand] : csTree.PhiD2_gen()[icand]; 
            Double_t negPt_gen  = csTree.chargeD1_gen()[icand] < 0 ? csTree.pTD1_gen()[icand] : csTree.pTD2_gen()[icand]; 
            Double_t negEta_gen = csTree.chargeD1_gen()[icand] < 0 ? csTree.EtaD1_gen()[icand] : csTree.EtaD2_gen()[icand]; 
            Double_t negPhi_gen = csTree.chargeD1_gen()[icand] < 0 ? csTree.PhiD1_gen()[icand] : csTree.PhiD2_gen()[icand]; 

            TVector3 genPosMom, genNegMom; 
            genPosMom.SetPtEtaPhi(posPt_gen, posEta_gen, posPhi_gen);
            genNegMom.SetPtEtaPhi(negPt_gen, negEta_gen, negPhi_gen);

            Int_t recoIdx = csTree.RecIdx_gen()[icand];

            if(recoIdx<0 || csTree.candSize()<=0) continue;

            Double_t posPt     = csTree.chargeD1()[recoIdx] > 0 ? csTree.pTD1()[recoIdx] : csTree.pTD2()[recoIdx]; 
            Double_t posEta    = csTree.chargeD1()[recoIdx] > 0 ? csTree.EtaD1()[recoIdx] : csTree.EtaD2()[recoIdx]; 
            Double_t posPhi    = csTree.chargeD1()[recoIdx] > 0 ? csTree.PhiD1()[recoIdx] : csTree.PhiD2()[recoIdx]; 
            Double_t isPosTrig = csTree.chargeD1()[recoIdx] > 0 ? csTree.trigMuon1()[trigIdx][recoIdx] : csTree.trigMuon2()[trigIdx][recoIdx]; 
            Double_t negPt     = csTree.chargeD1()[recoIdx] < 0 ? csTree.pTD1()[recoIdx] : csTree.pTD2()[recoIdx]; 
            Double_t negEta    = csTree.chargeD1()[recoIdx] < 0 ? csTree.EtaD1()[recoIdx] : csTree.EtaD2()[recoIdx]; 
            Double_t negPhi    = csTree.chargeD1()[recoIdx] < 0 ? csTree.PhiD1()[recoIdx] : csTree.PhiD2()[recoIdx]; 
            Double_t isNegTrig = csTree.chargeD1()[recoIdx] < 0 ? csTree.trigMuon1()[trigIdx][recoIdx] : csTree.trigMuon2()[trigIdx][recoIdx]; 

            if(!csTree.softCand(recoIdx)) continue;

            TVector3 recoPosMom, recoNegMom; 
            recoPosMom.SetPtEtaPhi(posPt, posEta, posPhi);
            recoNegMom.SetPtEtaPhi(negPt, negEta, negPhi);

            Double_t posDeltaR = genPosMom.DeltaR(recoPosMom);
            Double_t negDeltaR = genNegMom.DeltaR(recoNegMom);

            hDeltaR->Fill(posDeltaR); 
            hDeltaR->Fill(negDeltaR);

            hDeltaPt->Fill((posPt-posPt_gen)/posPt_gen); 
            hDeltaPt->Fill((negPt-negPt_gen)/negPt_gen);

            //shift mean
            posPt -= funPtMeanShift->Eval(posPt_gen)*posPt_gen;
            negPt -= funPtMeanShift->Eval(negPt_gen)*negPt_gen;

            Double_t posPtNew[nSmearScan];
            Double_t negPtNew[nSmearScan];
            for(Int_t iscan=0; iscan<nSmearScan; iscan++){
                //smear width
                posPtNew[iscan] = posPt_gen + (posPt - posPt_gen)/funRawPtRes->Eval(posPt_gen)*funNewPtRes[iscan]->Eval(posPt_gen);
                negPtNew[iscan] = negPt_gen + (negPt - negPt_gen)/funRawPtRes->Eval(negPt_gen)*funNewPtRes[iscan]->Eval(negPt_gen);

                TLorentzVector posFourMom, negFourMom, pairFourMom;
                posFourMom.SetPtEtaPhiM(posPtNew[iscan], posEta, posPhi, Mmuon);
                negFourMom.SetPtEtaPhiM(negPtNew[iscan], negEta, negPhi, Mmuon);
                pairFourMom = posFourMom + negFourMom;

                Double_t pt   = pairFourMom.Pt();
                Double_t eta  = pairFourMom.Eta();
                Double_t phi  = pairFourMom.Phi();
                Double_t mass = pairFourMom.M();
                Double_t y    = pairFourMom.Rapidity();

                //if(!goodMuPair(csTree, recoIdx)) continue;
                if(!goodMuPair(posFourMom, isPosTrig, negFourMom, isNegTrig)) continue;

                hMvsPtvsRap_Scan[iscan]->Fill(y, pt, mass);
            }
        }
    }

    system("mkdir -p smearedHistos");
    writeHistos(Form("smearedHistos/dimuonHistos.%s", fileName.Data()));
}

void bookHistos()
{
    const Int_t    mHistRapBins = 50;
    const Double_t mHistRapLow = -2.5;
    const Double_t mHistRapHi = 2.5;
    const Int_t    mHistPtBins = 80;
    const Double_t mHistPtLow = 0;
    const Double_t mHistPtHi = 0.8;
    const Int_t    mHistMassBins = 300;
    const Double_t mHistMassLow = 2;
    const Double_t mHistMassHi = 5;

    // event level
    hnEvts = new TH1D("hnEvts", "hnEvts;", 5, 0, 5);
    hnEvts->GetXaxis()->SetBinLabel(1, "trigEvt");
    hnEvts->GetXaxis()->SetBinLabel(2, "validVtx & !Beam-halo");
    hnEvts->GetXaxis()->SetBinLabel(3, "HFMaxE <= 7.6(7.3) GeV");
    hnEvts->GetXaxis()->SetBinLabel(4, "N_{trk}^{HP} == 2");
    hnEvts->GetXaxis()->LabelsOption("d");
    hnEvts->GetXaxis()->SetLabelSize(0.06);

    hDeltaR              = new TH1D("hDeltaR", "hDeltaR; #DeltaR", 2000, 0, 0.2);
    hDeltaPt             = new TH1D("hDeltaPt", "hDeltaPt; (p_{T}^{Rc}-p_{T}^{Gen})/p_{T}^{Gen}", 5000, -0.5, 0.5);

    for(Int_t i=0; i<nSmearScan; i++){
        Double_t par0 = mInitPar0 + i*mSmearStep;

        hMvsPtvsRap_Scan[i] = new TH3D(Form("hMvsPtvsRap_Scan%d", i), Form("par0 = %5.5f; Rapidity; p_{T} (GeV); m_{#mu#mu} (GeV)", par0), mHistRapBins, mHistRapLow, mHistRapHi, mHistPtBins, mHistPtLow, mHistPtHi, mHistMassBins, mHistMassLow, mHistMassHi);
    }
}

void writeHistos(TString fileName)
{
    TFile* fOut = new TFile(Form("%s.root", fileName.Data()), "recreate");

    fOut->cd();

    hnEvts->Write();
    hDeltaR->Write();
    hDeltaPt->Write();

    for(Int_t i=0; i<nSmearScan; i++){
        hMvsPtvsRap_Scan[i]->Write();
    }

    fOut->Close();
}

Bool_t goodMuPair(VertexCompositeTree& evtTree, const int icand)
{
    Double_t mTrkPtTh1  = fTrkAcc->Eval(evtTree.EtaD1()[icand]);
    Double_t mTrkPtTh2  = fTrkAcc->Eval(evtTree.EtaD2()[icand]);
    Double_t mTrigPtTh1 = fTrigAcc->Eval(evtTree.EtaD1()[icand]);
    Double_t mTrigPtTh2 = fTrigAcc->Eval(evtTree.EtaD2()[icand]);

    if(evtTree.pTD1()[icand] < mTrkPtTh1 || evtTree.pTD2()[icand] < mTrkPtTh2) return kFALSE;

    Bool_t isTrigAcc1 = kFALSE, isTrigAcc2 = kFALSE;
    if(evtTree.trigMuon1()[trigIdx][icand] && evtTree.pTD1()[icand] >= mTrigPtTh1) isTrigAcc1 = kTRUE;
    if(evtTree.trigMuon2()[trigIdx][icand] && evtTree.pTD2()[icand] >= mTrigPtTh2) isTrigAcc2 = kTRUE;

    if(!isTrigAcc1 && !isTrigAcc2) return kFALSE;

    return kTRUE;
}

Bool_t goodMuPair(const TLorentzVector posFourMom, const Bool_t isPosTrig, const TLorentzVector negFourMom, const Bool_t isNegTrig)
{
    Double_t mPosTrkPtTh  = fTrkAcc->Eval(posFourMom.Eta());
    Double_t mNegTrkPtTh  = fTrkAcc->Eval(negFourMom.Eta());
    Double_t mPosTrigPtTh = fTrigAcc->Eval(posFourMom.Eta());
    Double_t mNegTrigPtTh = fTrigAcc->Eval(negFourMom.Eta());

    if(posFourMom.Pt() < mPosTrkPtTh || negFourMom.Pt() < mNegTrkPtTh) return kFALSE;

    Bool_t isPosTrigAcc = kFALSE, isNegTrigAcc = kFALSE;
    if(isPosTrig && posFourMom.Pt() >= mPosTrigPtTh) isPosTrigAcc = kTRUE;
    if(isNegTrig && negFourMom.Pt() >= mNegTrigPtTh) isNegTrigAcc = kTRUE;

    if(!isPosTrigAcc && !isNegTrigAcc) return kFALSE;

    return kTRUE;
}
