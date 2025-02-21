#include "/afs/ihep.ac.cn/users/z/zhangyu1/bishe/jpsiAnaCode/common/headers.h"
#include "/afs/ihep.ac.cn/users/z/zhangyu1/bishe/jpsiAnaCode/common/VertexCompositeTree.h"
#include "/afs/ihep.ac.cn/users/z/zhangyu1/bishe/jpsiAnaCode/common/funUtil.h"

TH1D *hnEvts;
TH3D *hVzvsVyvsVx;
TH2D *hHFMinusvsHFPlus;
TH3D *hNCanvsNtrkvsCen;
TH2D *hNtrkofflinevsNtrkHP;

TH3D *hPhivsEtavsPt_Gen;
TH2D *hNegPtvsPosPt_Gen;
TH2D *hNegEtavsPosEta_Gen;
TH2D *hNegPhivsPosPhi_Gen;
TH2D *hPtResvsGenPt;
TH2D *hEtaResvsGenEta;
TH2D *hPhiResvsGenPhi;

TH2D *hNRcMuvsNGenMu;
TH3D *hPhivsEtavsPt;

TH2D *hRcPairPtvsGenPairPt;
TH2D *hRcPairEtavsGenPairEta;
TH2D *hRcPairPhivsGenPairPhi;
TH2D *hMassResvsGenMass;
TH2D *hRapResvsGenRap;

TH3D *hRawMvsPtvsRap_Gen;
TH3D *hMvsPtvsRap_Gen;
TH2D *hPt2vsM_Gen;
TH2D *hAsyPhivsM_Gen;
TH2D *hAsyPtvsM_Gen;
TH2D *hDeltaPhivsM_Gen;

TH2D *hAsyPhivsPt_Gen;
TH1D *hModifiedPairPt_Neu[nNeus][nNeus];

TH3D *hMvsPtvsRap;
TH3D *hMvsAsyPhivsRap;
TH2D *hPt2vsM;
TH2D *hAsyPhivsM;
TH2D *hAsyPhivsM_EffCorr;
TH2D *hAsyPtvsM;
TH2D *hDeltaPhivsM;

TH3D *hPosMuPhivsEtavsPt;
TH3D *hNegMuPhivsEtavsPt;

TH2D *hTnpDenEtavsPt;
TH2D *hTnpPassEtavsPt;
TH2D *hTnpFailEtavsPt;

// To calculate 3D efficiency
TH3D *hPosMuPhivsEtavsPt_Gen;
TH3D *hMthPosMuPhivsEtavsPt_Gen;
TH3D *hMthPosMuPhivsEtavsPt;
TH3D *hTrigPosMuPhivsEtavsPt;

TH3D *hNegMuPhivsEtavsPt_Gen;
TH3D *hMthNegMuPhivsEtavsPt_Gen;
TH3D *hMthNegMuPhivsEtavsPt;
TH3D *hTrigNegMuPhivsEtavsPt;

TH3D *hMthPosMuPhivsEtavsPtInpair;
TH3D *hTrigPosMuPhivsEtavsPtInpair;
TH3D *hMthNegMuPhivsEtavsPtInpair;
TH3D *hTrigNegMuPhivsEtavsPtInpair;

TH1D *hEvtvsCostheta;
TH1D *hEvtvsM_mumu;
TH1D *hEvtvsk_max;
TH1D *hEvtvsk_min;

Bool_t   goodMcTrack(VertexCompositeTree& evtTree, const int icand);
Bool_t   goodRcTrack(VertexCompositeTree& evtTree, const int icand);
Bool_t   matchedTrack(Double_t genPt, Double_t genEta, Double_t genPhi, Double_t rcPt, Double_t rcEta, Double_t rcPhi);
Double_t shiftDeltaPhi(Double_t dPhi);
Double_t shiftToPi(Double_t dPhi);

double_t etaToY(double_t eta, double_t mass, double_t pt);
double_t yyToy_mumu(double_t y1, double_t y2);
double_t PhotonEnergy(double_t y_mu_mu, double_t SysM_mumu);
double_t m_mumu(double_t posPt, double_t posEta, double_t posPhi, double_t posy,double_t negPt, double_t negEta, double_t negPhi, double_t negy);

void bookHistos();
void writeHistos(TString fileName = "test");

void anaMcEvt(TString fileName = "GammaGamma")
{
    TH1::SetDefaultSumw2(kTRUE);

    std::string inputFile;
    if(fileName.EqualTo("GammaGamma")){
        inputFile = "/publicfs/cms/user/tocheng/HeavyIon/UPC/2018A/dimuana_mc.root";
    }
    else if(fileName.EqualTo("GammaGamma_XnXn")){
        inputFile = "../rootfiles/VertexCompositeTree_STARLIGHT_GGToMuMu_XnXn_woPtCut_DiMuMC_20191125.root";
    }
    else if(fileName.EqualTo("CohY1S")){
        inputFile = "../rootfiles/VertexCompositeTree_STARLIGHT_Ups1SToMuMu_Coherent_wIF_woPtCut_DiMuMC_20191125.root";
    }
    else{
        cout<<"fileName is invalid !"<<endl;
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

    TF1 *funSTARlight = new TF1("funSTARlight", "[0]*exp(-x/[1])", 0, 1);
    funSTARlight->SetParameters(1, 1.3489e-3);

    Double_t slopeData[nPts] = {0.00118441, 0.00131362, 0.00134486, 0.00140852, 0.00147138, 0.00155273};
    TF1 *funData[nPts];
    for(Int_t i=0; i<nPts; i++){
        funData[i] = new TF1(Form("funData_Sce%d", i), "[0]*exp(-x/[1])", 0, 0.1);
        funData[i]->SetParameters(1, slopeData[i]);
        funData[i]->SetNpx(10000);
    }

    bookHistos();

    TRandom3 *rnd = new TRandom3(0);

    vector<int> anaRuns;
    vector<int> badZDCRuns;
    anaRuns.clear();
    badZDCRuns.clear();

    for (Long64_t jentry = 1; jentry < csTree.GetEntries(); jentry++) {
        if (jentry % (csTree.GetEntries() / 10) == 0)
            cout << "begin " << jentry << "th entry...." << endl;

        // Get the entry
        if (csTree.GetEntry(jentry) < 0) {
            cout << "Invalid correct-sign entry!" << endl;
            return;
        }

        Int_t   cen   = csTree.centrality();
        Float_t vtxX  = csTree.bestvtxX();
        Float_t vtxY  = csTree.bestvtxY();
        Float_t vtxZ  = csTree.bestvtxZ();
        Float_t hfsumETPlus  = csTree.HFsumETPlus();
        Float_t hfsumETMinus = csTree.HFsumETMinus();
        Int_t   nTrkoffline  = csTree.Ntrkoffline();
        Int_t   nTrkHP       = csTree.NtrkHP();

        hnEvts->Fill(0.5);
        hVzvsVyvsVx->Fill(vtxX, vtxY, vtxZ);
        hHFMinusvsHFPlus->Fill(hfsumETPlus, hfsumETMinus);
        hNCanvsNtrkvsCen->Fill(cen, nTrkoffline, csTree.candSize());
        hNtrkofflinevsNtrkHP->Fill(nTrkHP, nTrkoffline);

        // Select Event - require this event has a valid vertex and is not beam-halo event
        if(!csTree.evtSel()[2] || !csTree.evtSel()[3]) continue;
        hnEvts->Fill(1.5);

        //evtSel[4-15]
        //[4]=0: HFPlusMaxTower < 3 GeV;  [4]=1: HFPlusMaxTower > 3 GeV
        //[5]=0: HFMinusMaxTower < 3 GeV;  [5]=1: HFMinusMaxTower > 3 GeV
        //[6] is for Plus & [7] is for Minus; Threshold = 4 GeV
        //[8] is for Plus & [9] is for Minus; Threshold = 5 GeV
        //[10] is for Plus & [11] is for Minus; Threshold = 6 GeV
        //[12] is for Plus & [13] is for Minus; Threshold = 7 GeV
        //[14] is for Plus & [15] is for Minus; Threshold = 8 GeV
        //[16] is for Plus (Th = 7.6 GeV) & [17] is for Minus (Th = 7.3 GeV); 

        if(csTree.evtSel()[16] || csTree.evtSel()[17]) continue;
        hnEvts->Fill(2.5);

        if(nTrkHP != 2) continue;
        hnEvts->Fill(3.5);

        // Loop over the correct-sign candidates
        Int_t nSoftMuon = 0;
        for (UInt_t icand = 0; icand < csTree.candSize_gen(); icand++) {
            Double_t posPt_gen  = csTree.chargeD1_gen()[icand] > 0 ? csTree.pTD1_gen()[icand] : csTree.pTD2_gen()[icand]; 
            Double_t posEta_gen = csTree.chargeD1_gen()[icand] > 0 ? csTree.EtaD1_gen()[icand] : csTree.EtaD2_gen()[icand]; 
            Double_t posPhi_gen = csTree.chargeD1_gen()[icand] > 0 ? csTree.PhiD1_gen()[icand] : csTree.PhiD2_gen()[icand]; 
            Double_t negPt_gen  = csTree.chargeD1_gen()[icand] < 0 ? csTree.pTD1_gen()[icand] : csTree.pTD2_gen()[icand]; 
            Double_t negEta_gen = csTree.chargeD1_gen()[icand] < 0 ? csTree.EtaD1_gen()[icand] : csTree.EtaD2_gen()[icand]; 
            Double_t negPhi_gen = csTree.chargeD1_gen()[icand] < 0 ? csTree.PhiD1_gen()[icand] : csTree.PhiD2_gen()[icand]; 

            hPhivsEtavsPt_Gen->Fill(posPt_gen, posEta_gen, posPhi_gen);
            hPhivsEtavsPt_Gen->Fill(negPt_gen, negEta_gen, negPhi_gen);

            hPosMuPhivsEtavsPt_Gen->Fill(posPt_gen, posEta_gen, posPhi_gen);
            hNegMuPhivsEtavsPt_Gen->Fill(negPt_gen, negEta_gen, negPhi_gen);

            Int_t nPosMth = 0;
            Int_t nNegMth = 0;
            for(UInt_t imu=0; imu<csTree.candSize_mu(); imu++){
                if(!csTree.softMuon_mu()[imu]) continue;

                if(icand==0) nSoftMuon++;

                Double_t muPt  = csTree.pT_mu()[imu];
                Double_t muEta = csTree.eta_mu()[imu];
                Double_t muPhi = csTree.phi_mu()[imu];
                Bool_t   isTrigMu = csTree.trigMuon_mu()[trigIdx][imu];

                if(matchedTrack(posPt_gen, posEta_gen, posPhi_gen, muPt, muEta, muPhi)){
                    hPtResvsGenPt->Fill(posPt_gen, (muPt-posPt_gen)/posPt_gen);
                    hEtaResvsGenEta->Fill(posEta_gen, muEta-posEta_gen);
                    if(posPt_gen>mPtCut && TMath::Abs(posEta_gen)<mEtaCut){
                        hPhiResvsGenPhi->Fill(posPhi_gen, muPhi-posPhi_gen);
                    }

                    hMthPosMuPhivsEtavsPt_Gen->Fill(posPt_gen, posEta_gen, posPhi_gen);
                    hMthPosMuPhivsEtavsPt->Fill(muPt, muEta, muPhi);

                    if(isTrigMu){
                        hTrigPosMuPhivsEtavsPt->Fill(muPt, muEta, muPhi);
                    }

                    nPosMth++;
                }

                if(matchedTrack(negPt_gen, negEta_gen, negPhi_gen, muPt, muEta, muPhi)){
                    hPtResvsGenPt->Fill(negPt_gen, (muPt-negPt_gen)/negPt_gen);
                    hEtaResvsGenEta->Fill(negEta_gen, muEta-negEta_gen);
                    if(negPt_gen>mPtCut && TMath::Abs(negEta_gen)<mEtaCut){
                        hPhiResvsGenPhi->Fill(negPhi_gen, muPhi-negPhi_gen);
                    }

                    hMthNegMuPhivsEtavsPt_Gen->Fill(negPt_gen, negEta_gen, negPhi_gen);
                    hMthNegMuPhivsEtavsPt->Fill(muPt, muEta, muPhi);

                    if(isTrigMu){
                        hTrigNegMuPhivsEtavsPt->Fill(muPt, muEta, muPhi);
                    }

                    nNegMth++;
                }
            }

            if(nPosMth>1){
                cout<<"More than 1 reconstructed tracks matched with positive generated track !"<<endl;
            }

            if(nNegMth>1){
                cout<<"More than 1 reconstructed tracks matched with negative generated track !"<<endl;
            }

            hNegPtvsPosPt_Gen->Fill(posPt_gen, negPt_gen);
            hNegEtavsPosEta_Gen->Fill(posEta_gen, negEta_gen);
            hNegPhivsPosPhi_Gen->Fill(posPhi_gen, negPhi_gen);

            TLorentzVector posFourMom_gen, negFourMom_gen, pairFourMom_gen;
            posFourMom_gen.SetPtEtaPhiM(posPt_gen, posEta_gen, posPhi_gen, Mmuon);
            negFourMom_gen.SetPtEtaPhiM(negPt_gen, negEta_gen, negPhi_gen, Mmuon);
            pairFourMom_gen = posFourMom_gen + negFourMom_gen;

            Double_t pt_gen   = pairFourMom_gen.Pt();
            Double_t eta_gen  = pairFourMom_gen.Eta();
            Double_t phi_gen  = pairFourMom_gen.Phi();
            Double_t mass_gen = pairFourMom_gen.M();
            Double_t y_gen    = pairFourMom_gen.Rapidity();

            Double_t asyPhi = 1 - TMath::Abs(shiftDeltaPhi(posFourMom_gen.DeltaPhi(negFourMom_gen))) / PI; //acoplanarity
            Double_t asyPt  = TMath::Abs((posFourMom_gen.Pt() - negFourMom_gen.Pt()) / (posFourMom_gen.Pt() + negFourMom_gen.Pt()));

	    double_t posy_gen   = etaToY(posEta_gen, mass_gen, posPt_gen);
            double_t negy_gen   = etaToY(negEta_gen, mass_gen, negPt_gen);
            
            double_t costheta = 0;
            if (posy_gen > 0) {
                costheta = tanh(0.5 * (posy_gen - negy_gen));
            } else {
                costheta = tanh(0.5 * (negy_gen - posy_gen));
            }

            double_t y_mumu = yyToy_mumu(posy_gen, negy_gen);

            double_t SysM_mumu = m_mumu(posPt_gen, posEta_gen, posPhi_gen, posy_gen, negPt_gen, negEta_gen, negPhi_gen, negy_gen);

            //from the final state muon pair, calculate the photon energy
            double_t k_max = PhotonEnergy(y_mumu, SysM_mumu);
            double_t k_min = PhotonEnergy(-y_mumu, SysM_mumu);

            TVector3 muMomDiff_gen = posFourMom_gen.Vect() - negFourMom_gen.Vect();
            TVector3 pairMom_gen = pairFourMom_gen.Vect();
            Double_t phiDiff_gen = shiftDeltaPhi(pairMom_gen.DeltaPhi(muMomDiff_gen));

            hRawMvsPtvsRap_Gen->Fill(y_gen, pt_gen, mass_gen);

            if(goodMcTrack(csTree, icand) && TMath::Abs(y_gen) <= mPairYCut){
                hMvsPtvsRap_Gen->Fill(y_gen, pt_gen, mass_gen);
                hPt2vsM_Gen->Fill(mass_gen, pt_gen * pt_gen);
                hAsyPhivsM_Gen->Fill(mass_gen, asyPhi);
                hAsyPtvsM_Gen->Fill(mass_gen, asyPt);
                hDeltaPhivsM_Gen->Fill(mass_gen, phiDiff_gen);

                if(mass_gen>=8 && mass_gen<=60){
                    hAsyPhivsPt_Gen->Fill(posPt_gen, asyPhi);

                    Int_t neuIdx = 0;

                    for(Int_t ip=0; ip<nNeus; ip++){
                        for(Int_t im=ip; im<nNeus; im++){
                            Double_t weight = funData[neuIdx]->Eval(asyPhi) / funSTARlight->Eval(asyPhi);
                            hModifiedPairPt_Neu[ip][im]->Fill(pt_gen, weight);

                            neuIdx++;
                        }
                    }
                }
            }

            Int_t recIdx = csTree.RecIdx_gen()[icand];

            if(recIdx<0 || csTree.candSize()<=0) continue;

            Double_t posPt     = csTree.chargeD1()[recIdx] > 0 ? csTree.pTD1()[recIdx] : csTree.pTD2()[recIdx]; 
            Double_t posEta    = csTree.chargeD1()[recIdx] > 0 ? csTree.EtaD1()[recIdx] : csTree.EtaD2()[recIdx]; 
            Double_t posPhi    = csTree.chargeD1()[recIdx] > 0 ? csTree.PhiD1()[recIdx] : csTree.PhiD2()[recIdx]; 
            Double_t posIsTrig = csTree.chargeD1()[recIdx] > 0 ? csTree.trigMuon1()[trigIdx][recIdx] : csTree.trigMuon2()[trigIdx][recIdx]; 
            Double_t posIsSoft = csTree.chargeD1()[recIdx] > 0 ? csTree.softMuon1()[recIdx] : csTree.softMuon2()[recIdx]; 
            Double_t negPt     = csTree.chargeD1()[recIdx] < 0 ? csTree.pTD1()[recIdx] : csTree.pTD2()[recIdx]; 
            Double_t negEta    = csTree.chargeD1()[recIdx] < 0 ? csTree.EtaD1()[recIdx] : csTree.EtaD2()[recIdx]; 
            Double_t negPhi    = csTree.chargeD1()[recIdx] < 0 ? csTree.PhiD1()[recIdx] : csTree.PhiD2()[recIdx]; 
            Double_t negIsTrig = csTree.chargeD1()[recIdx] < 0 ? csTree.trigMuon1()[trigIdx][recIdx] : csTree.trigMuon2()[trigIdx][recIdx]; 
            Double_t negIsSoft = csTree.chargeD1()[recIdx] < 0 ? csTree.softMuon1()[recIdx] : csTree.softMuon2()[recIdx]; 

            if (!csTree.softCand(recIdx)) continue;

            hMthPosMuPhivsEtavsPtInpair->Fill(posPt, posEta, posPhi);
            if(posIsTrig) hTrigPosMuPhivsEtavsPtInpair->Fill(posPt, posEta, posPhi);

            hMthNegMuPhivsEtavsPtInpair->Fill(negPt, negEta, negPhi);
            if(negIsTrig) hTrigNegMuPhivsEtavsPtInpair->Fill(negPt, negEta, negPhi);

            Double_t pt   = csTree.pT()[recIdx];
            Double_t eta  = csTree.eta()[recIdx];
            Double_t phi  = csTree.phi()[recIdx];
            Double_t mass = csTree.mass()[recIdx];
            Double_t y    = csTree.y()[recIdx];

            hRcPairPtvsGenPairPt->Fill(pt_gen, pt);
            hRcPairEtavsGenPairEta->Fill(phi_gen, phi);
            hRcPairPhivsGenPairPhi->Fill(eta_gen, eta);

            TLorentzVector posFourMom, negFourMom, pairFourMom;
            posFourMom.SetPtEtaPhiM(posPt, posEta, posPhi, Mmuon);
            negFourMom.SetPtEtaPhiM(negPt, negEta, negPhi, Mmuon);
            pairFourMom = posFourMom + negFourMom;

            //cout<<"pt_rc_read: "<<pt<<"    "<<"pt_rc_cal: "<<pairFourMom.Pt()<<endl;

            asyPhi = 1 - TMath::Abs(shiftDeltaPhi(posFourMom.DeltaPhi(negFourMom))) / PI; //acoplanarity
            asyPt  = TMath::Abs((posFourMom.Pt() - negFourMom.Pt()) / (posFourMom.Pt() + negFourMom.Pt()));

            //if (!csTree.softCand(recIdx))              continue;
            if (!goodRcTrack(csTree, recIdx))            continue;
            if (!csTree.trigCand(trigIdx, recIdx, true)) continue; // true: single muon UPC; false: dimuon UPC
            if (TMath::Abs(y) > mPairYCut)               continue;

            hMassResvsGenMass->Fill(mass_gen, (mass-mass_gen)/mass_gen);
            hRapResvsGenRap->Fill(y_gen, y-y_gen);

            // NOTE, in the pair selection, at least one triggered muon has already been deposited
            if(!posIsTrig){
                hTnpDenEtavsPt->Fill(posPt, TMath::Abs(posEta));
                hTnpFailEtavsPt->Fill(posPt, TMath::Abs(posEta));
            }
            if(!negIsTrig){
                hTnpDenEtavsPt->Fill(negPt, TMath::Abs(negEta));
                hTnpFailEtavsPt->Fill(negPt, TMath::Abs(negEta));
            }
            if(posIsTrig && negIsTrig){
                Double_t rndNum = rnd->Uniform(-1, 1);
                if(rndNum>0){ // positive muon is probe muon
                    hTnpDenEtavsPt->Fill(posPt, TMath::Abs(posEta));
                    hTnpPassEtavsPt->Fill(posPt, TMath::Abs(posEta));
                }
                else{ // negative muon is probe muon
                    hTnpDenEtavsPt->Fill(negPt, TMath::Abs(negEta));
                    hTnpPassEtavsPt->Fill(negPt, TMath::Abs(negEta));
                }
            }

            if(
                    mass>massLow[0] && mass<massHi[nMBins-1]
                    && asyPhi<mAlphaCut
              ){
                hPosMuPhivsEtavsPt->Fill(posPt, posEta, posPhi);
                hNegMuPhivsEtavsPt->Fill(negPt, negEta, negPhi);
            }

            hMvsPtvsRap->Fill(y, pt, mass);
            hMvsAsyPhivsRap->Fill(y, asyPhi, mass);
            hPt2vsM->Fill(mass, pt * pt);
            hAsyPhivsM->Fill(mass, asyPhi);
            hAsyPtvsM->Fill(mass, asyPt);

            TVector3 muMomDiff = posFourMom.Vect() - negFourMom.Vect();
            TVector3 pairMom = pairFourMom.Vect();

            Double_t phiDiff = shiftDeltaPhi(pairMom.DeltaPhi(muMomDiff));
            hDeltaPhivsM->Fill(mass, phiDiff);

            Int_t posPtBin  = hPosMu3DMthEff->GetXaxis()->FindBin(posPt);
            Int_t posEtaBin = hPosMu3DMthEff->GetYaxis()->FindBin(posEta);
            Int_t posPhiBin = hPosMu3DMthEff->GetZaxis()->FindBin(posPhi);
            Int_t negPtBin  = hNegMu3DMthEff->GetXaxis()->FindBin(negPt);
            Int_t negEtaBin = hNegMu3DMthEff->GetYaxis()->FindBin(negEta);
            Int_t negPhiBin = hNegMu3DMthEff->GetZaxis()->FindBin(negPhi);

            if(posPtBin > hPosMu3DMthEff->GetNbinsX()) posPtBin = hPosMu3DMthEff->GetNbinsX();
            if(negPtBin > hNegMu3DMthEff->GetNbinsX()) negPtBin = hNegMu3DMthEff->GetNbinsX();

            Double_t trkEff = hPosMu3DMthEff->GetBinContent(posPtBin, posEtaBin, posPhiBin) * hNegMu3DMthEff->GetBinContent(negPtBin, negEtaBin, negPhiBin);
            Double_t trigEff = 1 - (1 - hPosMu3DTrigEff->GetBinContent(posPtBin, posEtaBin, posPhiBin)) * (1 - hNegMu3DTrigEff->GetBinContent(negPtBin, negEtaBin, negPhiBin));
            Double_t totEff = trkEff * trigEff;
            hAsyPhivsM_EffCorr->Fill(mass, asyPhi, 1/totEff);

            hEvtvsCostheta->Fill(costheta, 1/totEff);
            hEvtvsM_mumu->Fill(SysM_mumu, 1/totEff);
            hEvtvsk_max->Fill(k_max, 1/totEff);
            hEvtvsk_min->Fill(k_min, 1/totEff);
        }

        hNRcMuvsNGenMu->Fill(csTree.candSize_gen()*2, nSoftMuon);
    }

    writeHistos(Form("dimuonHistos.%s", fileName.Data()));
}

void bookHistos()
{
    const Int_t    mHistRapBins = 60;
    const Double_t mHistRapLow = -3;
    const Double_t mHistRapHi = 3;
    const Int_t    mHistPtBins = 200;
    const Double_t mHistPtLow = 0;
    const Double_t mHistPtHi = 1;
    const Int_t    mHistMassBins = 2000;
    const Double_t mHistMassLow = 0;
    const Double_t mHistMassHi = 100;
    const Int_t    mHistPt2Bins = 100;
    const Double_t mHistPt2Low = 0;
    const Double_t mHistPt2Hi = 0.1;
    const Int_t    mHistAsyPhiBins = 500;
    const Double_t mHistAsyPhiLow = 0;
    const Double_t mHistAsyPhiHi = 0.05;
    const Int_t    mHistAsyPtBins = 500;
    const Double_t mHistAsyPtLow = 0;
    const Double_t mHistAsyPtHi = 0.5;

    // event level
    hnEvts = new TH1D("hnEvts", "hnEvts;", 5, 0, 5);
    hnEvts->GetXaxis()->SetBinLabel(1, "trigEvt");
    hnEvts->GetXaxis()->SetBinLabel(2, "validVtx & !Beam-halo");
    hnEvts->GetXaxis()->SetBinLabel(3, "HFMaxE <= 7.6(7.3) GeV");
    hnEvts->GetXaxis()->SetBinLabel(4, "N_{trk}^{HP} == 2");
    hnEvts->GetXaxis()->LabelsOption("d");
    hnEvts->GetXaxis()->SetLabelSize(0.06);

    hVzvsVyvsVx = new TH3D("hVzvsVyvsVx", "hVzvsVyvsVx; V_{x} (cm); V_{y} (cm); V_{z} (cm)", 100, -1, 1, 100, -1, 1, 300, -30, 30);
    hHFMinusvsHFPlus = new TH2D("hHFMinusvsHFPlus", "hHFMinusvsHFPlus; HFsumETPlus; HFsumETMinus", 200, 0, 20, 200, 0, 20);
    hNCanvsNtrkvsCen = new TH3D("hNCanvsNtrkvsCen", "hNCanvsNtrkvsCen; Centrality; N_{trk}^{offline}; # CS pair", 50, 150, 200, 10, 0, 10, 5, 0, 5);
    hNtrkofflinevsNtrkHP = new TH2D("hNtrkofflinevsNtrkHP", "hNtrkofflinevsNtrkHP; N_{trk}^{HP}; N_{trk}^{offline}", 10, 0, 10, 10, 0, 10);

    hPhivsEtavsPt_Gen   = new TH3D("hPhivsEtavsPt_Gen", "hPhivsEtavsPt_Gen; p_{T} (GeV/c); #eta; #phi", 500, 0, 50, 300, -3, 3, 180, -PI, PI);
    hNegPtvsPosPt_Gen    = new TH2D("hNegPtvsPosPt_Gen", "hNegPtvsPosPt_Gen; #mu^{+} p_{T} (GeV/c); #mu^{-} p_{T} (GeV/c);", 500, 0, 50, 500, 0, 50);
    hNegEtavsPosEta_Gen  = new TH2D("hNegEtavsPosEta_Gen", "hNegEtavsPosEta_Gen; #mu^{+} #eta; #mu^{-} #eta;", 300, -3, 3, 300, -3, 3);
    hNegPhivsPosPhi_Gen  = new TH2D("hNegPhivsPosPhi_Gen", "hNegPhivsPosPhi_Gen; #mu^{+} #phi; #mu^{-} #phi;", 180, -PI, PI, 180, -PI, PI);
    hPtResvsGenPt        = new TH2D("hPtResvsGenPt", "hPtResvsGenPt; p_{T}^{Gen} (GeV/c); (p_{T}^{Rc}-p_{T}^{Gen})/p_{T}^{Gen};", 500, 0, 50, 300, -0.15, 0.15);
    hEtaResvsGenEta      = new TH2D("hEtaResvsGenEta", "hEtaResvsGenEta; #eta^{Gen}; #eta^{Rc}-#eta^{Gen};", 300, -3, 3, 200, -0.01, 0.01);
    hPhiResvsGenPhi      = new TH2D("hPhiResvsGenPhi", "hPhiResvsGenPhi; #phi^{Gen}; #phi^{Rc}-#phi^{Gen};", 180, -PI, PI, 200, -0.01, 0.01);

    hNRcMuvsNGenMu          = new TH2D("hNRcMuvsNGenMu", "hNRcMuvsNGenMu; # of muon (GEN); # of muon (Soft)", 10, 0, 10, 10, 0, 10);
    hPhivsEtavsPt           = new TH3D("hPhivsEtavsPt", "hPhivsEtavsPt; p_{T} (GeV/c); #eta; #phi", 500, 0, 50, 300, -3, 3, 180, -PI, PI);

    hRcPairPtvsGenPairPt   = new TH2D("hRcPairPtvsGenPairPt", "hRcPairPtvsGenPairPt; pair p_{T}^{Gen} (GeV/c); pair p_{T}^{Rc} (GeV/c);", mHistPtBins, mHistPtLow, mHistPtHi, mHistPtBins, mHistPtLow, mHistPtHi);
    hRcPairEtavsGenPairEta = new TH2D("hRcPairEtavsGenPairEta", "hRcPairEtavsGenPairEta; pair #eta^{Gen}; pair #eta^{Rc};", 160, -4, 4, 160, -4, 4);
    hRcPairPhivsGenPairPhi = new TH2D("hRcPairPhivsGenPairPhi", "hRcPairPhivsGenPairPhi; pair #phi^{Gen}; pair #phi^{Rc};", 120, -PI, PI, 120, -PI, PI);
    hMassResvsGenMass      = new TH2D("hMassResvsGenMass", "hMassResvsGenMass; M_{#mu#mu}^{Gen} (GeV/c^{2}); (M_{#mu#mu}^{Rc} - M_{#mu#mu}^{Gen})/M_{#mu#mu}^{Gen};", mHistMassBins, mHistMassLow, mHistMassHi, 300, -0.15, 0.15);
    hRapResvsGenRap        = new TH2D("hRapResvsGenRap", "hRapResvsGenRap; pair y^{Gen}; y^{Rc} - y^{Gen};", 600, -3, 3, 1000, -0.1, 0.1);

    hRawMvsPtvsRap_Gen = new TH3D("hRawMvsPtvsRap_Gen", "hRawMvsPtvsRap_Gen; Rapidity; p_{T} (GeV/c); M_{#mu#mu} (GeV/c^{2})", mHistRapBins, mHistRapLow, mHistRapHi, mHistPtBins, mHistPtLow, mHistPtHi, mHistMassBins, mHistMassLow, mHistMassHi);
    hMvsPtvsRap_Gen    = new TH3D("hMvsPtvsRap_Gen", "hMvsPtvsRap_Gen; Rapidity; p_{T} (GeV/c); M_{#mu#mu} (GeV/c^{2})", mHistRapBins, mHistRapLow, mHistRapHi, mHistPtBins, mHistPtLow, mHistPtHi, mHistMassBins, mHistMassLow, mHistMassHi);
    hPt2vsM_Gen        = new TH2D("hPt2vsM_Gen", "hPt2vsM_Gen; M_{#mu#mu} (GeV/c^{2}); p_{T}^{2} ((GeV/c)^{2})", mHistMassBins, mHistMassLow, mHistMassHi, mHistPt2Bins, mHistPt2Low, mHistPt2Hi);
    hAsyPhivsM_Gen     = new TH2D("hAsyPhivsM_Gen", "hAsyPhivsM_Gen; M_{#mu#mu} (GeV/c^{2}); #alpha", mHistMassBins, mHistMassLow, mHistMassHi, mHistAsyPhiBins, mHistAsyPhiLow, mHistAsyPhiHi);
    hAsyPtvsM_Gen      = new TH2D("hAsyPtvsM_Gen", "hAsyPtvsM_Gen; M_{#mu#mu} (GeV/c^{2}); p_{T} Asymmetry", mHistMassBins, mHistMassLow, mHistMassHi, mHistAsyPtBins, mHistAsyPtLow, mHistAsyPtHi);
    hDeltaPhivsM_Gen   = new TH2D("hDeltaPhivsM_Gen", "hDeltaPhivsM_Gen; M_{#mu#mu} (GeV/c^{2}); #phi_{#mu^{+}+#mu^{-}} - #phi_{#mu^{+}-#mu^{-}}", mHistMassBins, mHistMassLow, mHistMassHi, 120, -PI, PI);

    hAsyPhivsPt_Gen     = new TH2D("hAsyPhivsPt_Gen", "hAsyPhivsPt_Gen; p_{T} (GeV/c); #alpha", 1000, 0, 100, mHistAsyPhiBins, mHistAsyPhiLow, mHistAsyPhiHi);

    hEvtvsCostheta = new TH1D("hEvtvsCostheta", "hEvtvsCostheta; cos #theta", 100, -1, 1);
    hEvtvsM_mumu = new TH1D("hEvtvsM_mumu", "hEvtvsM_mumu; M_{#mu#mu} (GeV/c^{2})", 52, 8, 60);
    hEvtvsk_max = new TH1D("hEvtvsk_max", "hEvtvsk_max; k_{max} (GeV)", 52, 8, 60);
    hEvtvsk_min = new TH1D("hEvtvsk_min", "hEvtvsk_min; k_{min} (GeV)", 100, 0, 100);

    for(Int_t ip=0; ip<nNeus; ip++){
        for(Int_t im=ip; im<nNeus; im++){
            hModifiedPairPt_Neu[ip][im] = new TH1D(Form("hModifiedPairPt_Neu%dn%dn", ip, im), "hModifiedPairPt; p_{T} (GeV/c)", 10000, mHistPtLow, mHistPtHi);
        }
    }

    hMvsPtvsRap = new TH3D("hMvsPtvsRap", "hMvsPtvsRap; Rapidity; p_{T} (GeV/c); M_{#mu#mu} (GeV/c^{2})", mHistRapBins, mHistRapLow, mHistRapHi, mHistPtBins, mHistPtLow, mHistPtHi, mHistMassBins, mHistMassLow, mHistMassHi);
    hMvsAsyPhivsRap = new TH3D("hMvsAsyPhivsRap", "hMvsAsyPhivsRap; Rapidity; #alpha; M_{#mu#mu} (GeV/c^{2})", mHistRapBins, mHistRapLow, mHistRapHi, mHistAsyPhiBins, mHistAsyPhiLow, mHistAsyPhiHi, mHistMassBins, mHistMassLow, mHistMassHi);
    hPt2vsM = new TH2D("hPt2vsM", "hPt2vsM; M_{#mu#mu} (GeV/c^{2}); p_{T}^{2} ((GeV/c)^{2})", mHistMassBins, mHistMassLow, mHistMassHi, mHistPt2Bins, mHistPt2Low, mHistPt2Hi);
    hAsyPhivsM = new TH2D("hAsyPhivsM", "hAsyPhivsM; M_{#mu#mu} (GeV/c^{2}); #alpha", mHistMassBins, mHistMassLow, mHistMassHi, mHistAsyPhiBins, mHistAsyPhiLow, mHistAsyPhiHi);
    hAsyPhivsM_EffCorr = new TH2D("hAsyPhivsM_EffCorr", "hAsyPhivsM_EffCorr; M_{#mu#mu} (GeV/c^{2}); #alpha", mHistMassBins, mHistMassLow, mHistMassHi, mHistAsyPhiBins, mHistAsyPhiLow, mHistAsyPhiHi);
    hAsyPtvsM = new TH2D("hAsyPtvsM", "hAsyPtvsM; M_{#mu#mu} (GeV/c^{2}); p_{T} Asymmetry", mHistMassBins, mHistMassLow, mHistMassHi, mHistAsyPtBins, mHistAsyPtLow, mHistAsyPtHi);
    hDeltaPhivsM = new TH2D("hDeltaPhivsM", "hDeltaPhivsM; M_{#mu#mu} (GeV/c^{2}); #phi_{#mu^{+}+#mu^{-}} - #phi_{#mu^{+}-#mu^{-}}", mHistMassBins, mHistMassLow, mHistMassHi, 120, -PI, PI);


    const Int_t nPtTailBins = 10;
    Double_t    PtTail[nPtTailBins+1] = {3.2, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 8.0, 10.0, 20.0};
    Double_t    mPtStep = 0.2;

    Int_t       mBinIdx = PtTail[0]/mPtStep;
    const Int_t nPtBins = mBinIdx + nPtTailBins;
    Double_t Pt[1000];
    for(Int_t ibin=0; ibin<mBinIdx; ibin++){
        Pt[ibin] = mPtStep*ibin;
    }
    for(Int_t ibin=0; ibin<=nPtTailBins; ibin++){
        Pt[mBinIdx+ibin] = PtTail[ibin];
    }

    const Int_t nEtaBins = 20;
    Double_t Eta[nEtaBins+1];
    Double_t mEtaLow = -3, mEtaHi = 3;
    Double_t mEtaStep = (mEtaHi - mEtaLow)/nEtaBins;
    for(Int_t ieta=0; ieta<=nEtaBins; ieta++){
        Eta[ieta] = mEtaLow + ieta*mEtaStep;
    }

    const Int_t nPhiBins = 60;
    Double_t Phi[nPhiBins+1];
    Double_t mPhiLow = -PI, mPhiHi = PI;
    Double_t mPhiStep = (mPhiHi - mPhiLow)/nPhiBins;
    for(Int_t iphi=0; iphi<=nPhiBins; iphi++){
        Phi[iphi] = mPhiLow + iphi*mPhiStep;
    }

    cout<<endl;
    cout<<"Pt Bin Boundary: "<<endl;
    for(Int_t ibin=0; ibin<nPtBins; ibin++){
        cout<<Pt[ibin]<<", ";
    }
    cout<<Pt[nPtBins]<<endl;
    cout<<endl;

    cout<<"Eta Bin Boundary: "<<endl;
    for(Int_t ibin=0; ibin<nEtaBins; ibin++){
        cout<<Eta[ibin]<<", ";
    }
    cout<<Eta[nEtaBins]<<endl;
    cout<<endl;

    cout<<"Phi Bin Boundary: "<<endl;
    for(Int_t ibin=0; ibin<nPhiBins; ibin++){
        cout<<Phi[ibin]<<", ";
    }
    cout<<Phi[nPhiBins]<<endl;
    cout<<endl;

    hPosMuPhivsEtavsPt_Gen    = new TH3D("hPosMuPhivsEtavsPt_Gen", "hPosMuPhivsEtavsPt_Gen; p_{T} (GeV/c); #eta; #phi", nPtBins, Pt, nEtaBins, Eta, nPhiBins, Phi);
    hMthPosMuPhivsEtavsPt_Gen = new TH3D("hMthPosMuPhivsEtavsPt_Gen", "hMthPosMuPhivsEtavsPt_Gen; p_{T} (GeV/c); #eta; #phi", nPtBins, Pt, nEtaBins, Eta, nPhiBins, Phi);
    hMthPosMuPhivsEtavsPt     = new TH3D("hMthPosMuPhivsEtavsPt", "hMthPosMuPhivsEtavsPt; p_{T} (GeV/c); #eta; #phi", nPtBins, Pt, nEtaBins, Eta, nPhiBins, Phi);
    hTrigPosMuPhivsEtavsPt    = new TH3D("hTrigPosMuPhivsEtavsPt", "hTrigPosMuPhivsEtavsPt; p_{T} (GeV/c); #eta; #phi", nPtBins, Pt, nEtaBins, Eta, nPhiBins, Phi);
    hNegMuPhivsEtavsPt_Gen    = new TH3D("hNegMuPhivsEtavsPt_Gen", "hNegMuPhivsEtavsPt_Gen; p_{T} (GeV/c); #eta; #phi", nPtBins, Pt, nEtaBins, Eta, nPhiBins, Phi);
    hMthNegMuPhivsEtavsPt_Gen = new TH3D("hMthNegMuPhivsEtavsPt_Gen", "hMthNegMuPhivsEtavsPt_Gen; p_{T} (GeV/c); #eta; #phi", nPtBins, Pt, nEtaBins, Eta, nPhiBins, Phi);
    hMthNegMuPhivsEtavsPt     = new TH3D("hMthNegMuPhivsEtavsPt", "hMthNegMuPhivsEtavsPt; p_{T} (GeV/c); #eta; #phi", nPtBins, Pt, nEtaBins, Eta, nPhiBins, Phi);
    hTrigNegMuPhivsEtavsPt    = new TH3D("hTrigNegMuPhivsEtavsPt", "hTrigNegMuPhivsEtavsPt; p_{T} (GeV/c); #eta; #phi", nPtBins, Pt, nEtaBins, Eta, nPhiBins, Phi);

    hMthPosMuPhivsEtavsPtInpair  = new TH3D("hMthPosMuPhivsEtavsPtInpair", "hMthPosMuPhivsEtavsPtInpair; p_{T} (GeV/c); #eta; #phi", nPtBins, Pt, nEtaBins, Eta, nPhiBins, Phi);
    hTrigPosMuPhivsEtavsPtInpair = new TH3D("hTrigPosMuPhivsEtavsPtInpair", "hTrigPosMuPhivsEtavsPtInpair; p_{T} (GeV/c); #eta; #phi", nPtBins, Pt, nEtaBins, Eta, nPhiBins, Phi);
    hMthNegMuPhivsEtavsPtInpair  = new TH3D("hMthNegMuPhivsEtavsPtInpair", "hMthNegMuPhivsEtavsPtInpair; p_{T} (GeV/c); #eta; #phi", nPtBins, Pt, nEtaBins, Eta, nPhiBins, Phi);
    hTrigNegMuPhivsEtavsPtInpair = new TH3D("hTrigNegMuPhivsEtavsPtInpair", "hTrigNegMuPhivsEtavsPtInpair; p_{T} (GeV/c); #eta; #phi", nPtBins, Pt, nEtaBins, Eta, nPhiBins, Phi);

    hPosMuPhivsEtavsPt = new TH3D("hPosMuPhivsEtavsPt", "hPosMuPhivsEtavsPt; p_{T} (GeV/c); #eta; #phi", 100, 0, 50, 24, -2.4, 2.4, 60, -PI, PI);
    hNegMuPhivsEtavsPt = new TH3D("hNegMuPhivsEtavsPt", "hNegMuPhivsEtavsPt; p_{T} (GeV/c); #eta; #phi", 100, 0, 50, 24, -2.4, 2.4, 60, -PI, PI);

    hTnpDenEtavsPt = new TH2D("hTnpDenEtavsPt", "hTnpDenEtavsPt; p_{T} (GeV/c); #eta", 300, 0, 30, 24, 0, 2.4);
    hTnpPassEtavsPt = new TH2D("hTnpPassEtavsPt", "hTnpPassEtavsPt; p_{T} (GeV/c); #eta", 300, 0, 30, 24, 0, 2.4);
    hTnpFailEtavsPt = new TH2D("hTnpFailEtavsPt", "hTnpFailEtavsPt; p_{T} (GeV/c); #eta", 300, 0, 30, 24, 0, 2.4);
}

void writeHistos(TString fileName)
{
    TFile* fOut = new TFile(Form("%s.root", fileName.Data()), "recreate");

    fOut->cd();

    hnEvts->Write();
    hVzvsVyvsVx->Write();
    hHFMinusvsHFPlus->Write();
    hNCanvsNtrkvsCen->Write();
    hNtrkofflinevsNtrkHP->Write();

    hPhivsEtavsPt_Gen->Write();
    hNegPtvsPosPt_Gen->Write();
    hNegEtavsPosEta_Gen->Write();
    hNegPhivsPosPhi_Gen->Write();
    hPtResvsGenPt->Write();
    hEtaResvsGenEta->Write();
    hPhiResvsGenPhi->Write();

    hNRcMuvsNGenMu->Write();
    hPhivsEtavsPt->Write();

    hRcPairPtvsGenPairPt->Write();
    hRcPairEtavsGenPairEta->Write();
    hRcPairPhivsGenPairPhi->Write();
    hMassResvsGenMass->Write();
    hRapResvsGenRap->Write();

    hRawMvsPtvsRap_Gen->Write();
    hMvsPtvsRap_Gen->Write();
    hPt2vsM_Gen->Write();
    hAsyPhivsM_Gen->Write();
    hAsyPtvsM_Gen->Write();
    hDeltaPhivsM_Gen->Write();

    hAsyPhivsPt_Gen->Write();
    for(Int_t ip=0; ip<nNeus; ip++){
        for(Int_t im=ip; im<nNeus; im++){
            hModifiedPairPt_Neu[ip][im]->Write();
        }
    }

    hMvsPtvsRap->Write();
    hMvsAsyPhivsRap->Write();
    hPt2vsM->Write();
    hAsyPhivsM->Write();
    hAsyPhivsM_EffCorr->Write();
    hAsyPtvsM->Write();
    hDeltaPhivsM->Write();

    hPosMuPhivsEtavsPt_Gen->Write();
    hMthPosMuPhivsEtavsPt_Gen->Write();
    hMthPosMuPhivsEtavsPt->Write();
    hTrigPosMuPhivsEtavsPt->Write();

    hNegMuPhivsEtavsPt_Gen->Write();
    hMthNegMuPhivsEtavsPt_Gen->Write();
    hMthNegMuPhivsEtavsPt->Write();
    hTrigNegMuPhivsEtavsPt->Write();

    hMthPosMuPhivsEtavsPtInpair->Write();
    hTrigPosMuPhivsEtavsPtInpair->Write();
    hMthNegMuPhivsEtavsPtInpair->Write();
    hTrigNegMuPhivsEtavsPtInpair->Write();

    hPosMuPhivsEtavsPt->Write();
    hNegMuPhivsEtavsPt->Write();

    hTnpDenEtavsPt->Write();
    hTnpPassEtavsPt->Write();
    hTnpFailEtavsPt->Write();
   
    hEvtvsCostheta->Write();
    hEvtvsM_mumu->Write();
    hEvtvsk_max->Write();
    hEvtvsk_min->Write();


    fOut->Close();
}

Bool_t goodMcTrack(VertexCompositeTree& evtTree, const int icand)
{
    if(evtTree.pTD1_gen()[icand] < mPtCut || evtTree.pTD2_gen()[icand] < mPtCut) return kFALSE;
    if(TMath::Abs(evtTree.EtaD1_gen()[icand]) > mEtaCut || TMath::Abs(evtTree.EtaD2_gen()[icand]) > mEtaCut) return kFALSE;

    return kTRUE;
}

Bool_t goodRcTrack(VertexCompositeTree& evtTree, const int icand)
{
    if(evtTree.pTD1()[icand] < mPtCut || evtTree.pTD2()[icand] < mPtCut) return kFALSE;
    if(TMath::Abs(evtTree.EtaD1()[icand]) > mEtaCut || TMath::Abs(evtTree.EtaD2()[icand]) > mEtaCut) return kFALSE;

    return kTRUE;
}

Bool_t matchedTrack(Double_t genPt, Double_t genEta, Double_t genPhi, Double_t rcPt, Double_t rcEta, Double_t rcPhi)
{
    if(TMath::Abs((rcPt-genPt)/genPt) > 0.15) return kFALSE;
    if(TMath::Abs(rcEta - genEta) > 0.01)     return kFALSE;
    if(TMath::Abs(rcPhi - genPhi) > 0.01)     return kFALSE;

    return kTRUE;
}

Double_t shiftDeltaPhi(Double_t dPhi)
{
    if (isnan(dPhi))
        return -999;

    while (dPhi < -PI)
        dPhi += 2 * PI;
    while (dPhi >= PI)
        dPhi -= 2 * PI;

    return dPhi;
}

Double_t shiftToPi(Double_t dPhi)
{
    if (isnan(dPhi))
        return -999;

    Double_t deltaPhi = shiftDeltaPhi(dPhi);
    if (deltaPhi < 0)
        deltaPhi += PI;

    return deltaPhi;
}

double_t etaToY(double_t eta, double_t mass, double_t pt) 
{
    double_t theta = 2 * atan(exp(-eta));
    double_t pz = pt / tan(theta); // 计算纵向动量
    double_t energy = sqrt(pt * pt + pz * pz + mass * mass);
    double_t rapidity = 0.5 * log((energy + pz) / (energy - pz));
    return rapidity;

}

double_t yyToy_mumu(double_t y1, double_t y2) 
{
    return 0.5 * log((exp(y1) + exp(y2)) / (exp(-y1) + exp(-y2)));
}

double_t PhotonEnergy(double_t y_mu_mu, double_t SysM_mumu) 
{
    return 0.5 * SysM_mumu * exp(y_mu_mu); // 计算 k
}

double_t m_mumu(double_t posPt, double_t posEta, double_t posPhi, double_t posy,double_t negPt, double_t negEta, double_t negPhi, double_t negy) 
{
    double_t mass = 0.10566;  // 缪子静止质量 (GeV/c^2)
    
    // 计算动量分量和总能量
    double_t posPz = posPt * sinh(posy);
    double_t negPz = negPt * sinh(negy);
    double_t posE = sqrt(posPt * posPt + posPz * posPz + mass * mass);
    double_t negE = sqrt(negPt * negPt + negPz * negPz + mass * mass);
    double_t posPx = posPt * cos(posPhi), posPy = posPt * sin(posPhi);
    double_t negPx = negPt * cos(negPhi), negPy = negPt * sin(negPhi);

    // 总能量与总动量
    double_t totalE = posE + negE;
    double_t totalPx = posPx + negPx, totalPy = posPy + negPy, totalPz = posPz + negPz;
    double_t totalP2 = totalPx * totalPx + totalPy * totalPy + totalPz * totalPz;

    return sqrt(totalE * totalE - totalP2); // 系统质量
}


