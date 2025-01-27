#include "/afs/ihep.ac.cn/users/z/zhangyu1/bishe/jpsiAnaCode/common/headers.h"
#include "/afs/ihep.ac.cn/users/z/zhangyu1/bishe/jpsiAnaCode/common/VertexCompositeTree.h"
#include "/afs/ihep.ac.cn/users/z/zhangyu1/bishe/jpsiAnaCode/common/funUtil.h"
#include "/afs/ihep.ac.cn/users/z/zhangyu1/bishe/jpsiAnaCode/common/tnp_weight_lowptPbPb.h"

TH1D *hnEvts;
TH2D *hCenvsTrig;
TH2D *hPreScalevsTrig;

TH3D *hVzvsVyvsVx;
TH3D *hVzvsVyvsVx_Sel;

TH1D *hRawCen;
TH1D *hCen_Final;

TH3D *hHFMinusvsHFPlusvsCen;
TH3D *hHFMinusvsHFPlusvsCen_Sel;

TH2D *hFlagvsBit;

TH3D *hNtrkHPvsNtrkofflinevsCen;
TH3D *hNtrkHPvsNtrkofflinevsCen_Sel;
TH1D *hNtrkHP_2SoftMuons;

TH2D *hRawZDCMinusvsZDCPlus;
TH2D *hZDCMinusvsZDCPlus;
TH2D *hZDCMinusvsZDCPlus_Sel;
TH2D *hZDCMinusvsZDCPlus_Only2MuTrk;

TH2D *hZDCMinusvsZDCPlus_LS;
TH2D *hZDCMinusvsZDCPlus_Sel_LS;
TH2D *hZDCMinusvsZDCPlus_Only2MuTrk_LS;

TH2D *hZDCvsNeuNum[nDirs];
TH2D *hZDCvsRegIdx[nDirs];

TH2D *hNeuNumMinusvsNeuNumPlus;

TH2D *hDeltaEtavsDetaPhi;

TH3D *hCS_MvsPtvsCen;
TH3D *hCS_MvsPtvsRap;
TH3D *hCS_MvsAsyPhivsRap;
TH3D *hWS_MvsPtvsCen;
TH3D *hWS_MvsPtvsRap;
TH3D *hMvsPtvsCen;
TH3D *hMvsPtvsRap;

TH3D *hPosMuPhivsEtavsPt;
TH3D *hNegMuPhivsEtavsPt;

TH2D *hTnpDenEtavsPt;
TH2D *hTnpPassEtavsPt;
TH2D *hTnpFailEtavsPt;

TH1D *hM_FineBin;
TH1D *hM_CoarseBin;

TH2D *hPt2vsM;
TH2D *hAsyPhivsM;
TH2D *hAsyPhivsRap;
TH2D *hAsyPtvsM;
TH2D *hDeltaPhivsM;
TH2D *hDeltaPhivsPt;

TH3D *hRapvsNMinusNeuvsNPlusNeu;

TH2D *hAsyPhivsM_Neu[nNeus][nNeus];
TH2D *hAsyPhivsM_Neu2[nNeus][nDirs];

TH2D *hAsyPhivsRap_Neu[nNeus][nNeus];
TH2D *hAsyPhivsRap_Neu2[nNeus][nDirs];

TH1D *hEvtvsCostheta;
TH1D *hEvtvsM_mumu;
TH1D *hEvtvsk_max;
TH1D *hEvtvsk_min;

Bool_t   goodTrack(VertexCompositeTree& evtTree, const int icand);
Int_t    grabNeutronNum(Int_t dirIdx, Double_t ranNum, Double_t prob[][nNeus]);
Int_t    grabNeutronNum(Int_t dirIdx, Double_t zdc);
Int_t    grabZDCRegionIdx(Int_t dirIdx, Double_t zdc);
Double_t shiftDeltaPhi(Double_t dPhi);
Double_t shiftToPi(Double_t dPhi);

void bookHistos();
void writeHistos(TString fileName = "test");

double_t etaToY(double_t eta, double_t mass, double_t pt);

double_t yyToy_mumu(double_t y1, double_t y2);

double_t PhotonEnergy(double_t y_mu_mu);

double_t m_mumu(double_t posPt, double_t posEta, double_t posPhi, double_t posy,double_t negPt, double_t negEta, double_t negPhi, double_t negy);


// *** load UPC TnP ***
const Int_t nEtaBins = 4;
Double_t    etaBoundary[nEtaBins+1] = {0, 1.2, 1.8, 2.1, 2.4};

TFile *fUpcTnp = 0x0;
TH1D  *hSFStat_MuIdTnp[nEtaBins];
TH1D  *hSFStat_TrigTnp[nEtaBins];

// Bool_t loadUpcTnPFile();
Int_t  grabEtaIdx(Double_t eta);
// ******

void anaEvt(Bool_t effCorr = kTRUE, Bool_t applyTnPSF = kFALSE, Bool_t incHadron = kFALSE, TString hfVetoType="Default")
{
    TH1::SetDefaultSumw2(kTRUE);


    if(!hfVetoType.EqualTo("Default") && !hfVetoType.EqualTo("Tight") && !hfVetoType.EqualTo("removeHF")){
        cout<<"Please input the correct hfVetoType string: 'Default' OR 'Tight' OR 'removeHF'"<<endl;
        return;
    }

    if(!effCorr && applyTnPSF){
        cout<<"Cannot apply TnP scaling factor without implementing efficiency corrections !"<<endl;
        return;
    }

    const auto& inputFile = "/publicfs/cms/user/tocheng/HeavyIon/UPC/2018A/dimuana_data.root";

    const auto& csTreeDir = "dimucontana";           // For MC use dimucontana_mc
    const auto& wsTreeDir = "dimucontana_wrongsign"; // For MC use dimucontana_wrongsign_mc

    // Extract the tree
    VertexCompositeTree csTree;
    VertexCompositeTree wsTree;

    if (!csTree.GetTree(inputFile, csTreeDir)) {
        cout << "Invalid Correct-Sign tree!" << endl;
        return;
    }
    if (!wsTree.GetTree(inputFile, wsTreeDir)) {
        cout << "Invalid Wrong-Sign tree!" << endl;
        return;
    }

    if(!init()) {
        cout<<"Initialization failed !"<<endl;
        return;
    }

    // if(!loadUpcTnPFile()){
    //     cout<<"Failed to load UPC TnP file !"<<endl;
    //     return;
    // }

    bookHistos();

    TRandom3 *rnd = new TRandom3(0);

    vector<int> anaRuns;
    anaRuns.clear();

    // Loop over csTree and wsTree ---> NOTE, these two trees have the exactly same event information
    for (Long64_t jentry = 1; jentry < csTree.GetEntries(); jentry++) {
        if (jentry % (csTree.GetEntries() / 10) == 0)
            cout << "begin " << jentry << "th entry...." << endl;

        // Get the entry
        if (csTree.GetEntry(jentry) < 0) {
            cout << "Invalid correct-sign entry!" << endl;
            return;
        }

        UInt_t  runNb = csTree.RunNb();
        Int_t   cen   = csTree.centrality();
        Float_t vtxX  = csTree.bestvtxX();
        Float_t vtxY  = csTree.bestvtxY();
        Float_t vtxZ  = csTree.bestvtxZ();
        Float_t zdcPlus  = csTree.ZDCPlus();
        Float_t zdcMinus = csTree.ZDCMinus();
        Float_t hfsumETPlus  = csTree.HFsumETPlus();
        Float_t hfsumETMinus = csTree.HFsumETMinus();
        Int_t   nTrkoffline  = csTree.Ntrkoffline();
        Int_t   nTrkHP       = csTree.NtrkHP();

        auto runIt = std::find(anaRuns.begin(), anaRuns.end(), runNb);
        if(runIt == anaRuns.end()) anaRuns.push_back(runNb);

        if(runNb < mRunNbCut) continue; // ZDC calibration only be valid since this run

        for(size_t itrig = 0; itrig < nTrigs; itrig++){
            if(csTree.trigHLT()[itrig]){
                hCenvsTrig->Fill(itrig, cen);
                hPreScalevsTrig->Fill(itrig, csTree.trigPrescale()[itrig]);
            }
        }

        // Select trigger
        if(!csTree.trigHLT()[trigIdx]) continue;
        hnEvts->Fill(0.5);

        hRawCen->Fill(cen);
        hRawZDCMinusvsZDCPlus->Fill(zdcPlus, zdcMinus);

        hVzvsVyvsVx->Fill(vtxX, vtxY, vtxZ);

        // Select Event - require this event has a valid vertex and is not beam-halo event
        if(!csTree.evtSel()[2] || !csTree.evtSel()[3]) continue;
        hnEvts->Fill(1.5);

        hVzvsVyvsVx_Sel->Fill(vtxX, vtxY, vtxZ);

        hZDCMinusvsZDCPlus->Fill(zdcPlus, zdcMinus);
        hZDCMinusvsZDCPlus_LS->Fill(zdcPlus, zdcMinus);
        hHFMinusvsHFPlusvsCen->Fill(cen, hfsumETPlus, hfsumETMinus);
        hNtrkHPvsNtrkofflinevsCen->Fill(cen, nTrkoffline, nTrkHP);

        Int_t nSoftMuonPair = 0;
        for (UInt_t icand = 0; icand < csTree.candSize(); icand++) {
            if (!csTree.softCand(icand)) continue;
            nSoftMuonPair++;
        }

        if(nSoftMuonPair==1 && nTrkHP==2){
            for(Int_t ibit=0; ibit<18; ibit++){
                if(csTree.evtSel()[ibit]) hFlagvsBit->Fill(ibit, 1);
                else                      hFlagvsBit->Fill(ibit, 0);
            }
        }

        //evtSel[4-15]
        //[4]=0: HFPlusMaxTower < 3 GeV;  [4]=1: HFPlusMaxTower > 3 GeV
        //[5]=0: HFMinusMaxTower < 3 GeV;  [5]=1: HFMinusMaxTower > 3 GeV
        //[6] is for Plus & [7] is for Minus; Threshold = 4 GeV
        //[8] is for Plus & [9] is for Minus; Threshold = 5 GeV
        //[10] is for Plus & [11] is for Minus; Threshold = 6 GeV
        //[12] is for Plus & [13] is for Minus; Threshold = 7 GeV
        //[14] is for Plus & [15] is for Minus; Threshold = 8 GeV
        //[16] is for Plus (Th = 7.3 GeV) & [17] is for Minus (Th = 7.6 GeV); 
        //
        //if(csTree.evtSel()[8] || csTree.evtSel()[9]) continue; // leading HF Energy < 5 GeV
        //if(csTree.evtSel()[12] || csTree.evtSel()[13]) continue; // leading HF Energy < 7 GeV

        Bool_t hfVeto = kTRUE;
        if(hfVetoType.EqualTo("Default"))  hfVeto = !csTree.evtSel()[16] && !csTree.evtSel()[17]; // leading HF Energy < 7.3 GeV (Plus), < 7.6 GeV (Minus)
        if(hfVetoType.EqualTo("Tight"))    hfVeto = !csTree.evtSel()[8] && !csTree.evtSel()[9];   // leading HF Energy < 5 GeV (Plus), < 5 GeV (Minus)
        if(hfVetoType.EqualTo("removeHF")) hfVeto = kTRUE; // remove HF veto
        if(!hfVeto) continue; 
        hnEvts->Fill(2.5);

        hZDCMinusvsZDCPlus_Sel->Fill(zdcPlus, zdcMinus);
        hZDCMinusvsZDCPlus_Sel_LS->Fill(zdcPlus, zdcMinus);
        hHFMinusvsHFPlusvsCen_Sel->Fill(cen, hfsumETPlus, hfsumETMinus);
        hNtrkHPvsNtrkofflinevsCen_Sel->Fill(cen, nTrkoffline, nTrkHP);

        if( nSoftMuonPair==1 ){
            hNtrkHP_2SoftMuons->Fill(nTrkHP);
            if(nTrkHP == 2){
                hZDCMinusvsZDCPlus_Only2MuTrk->Fill(zdcPlus, zdcMinus);
                hZDCMinusvsZDCPlus_Only2MuTrk_LS->Fill(zdcPlus, zdcMinus);
                hCen_Final->Fill(cen);
            }
        }

        if(!incHadron && nTrkHP!=2) continue;
        hnEvts->Fill(3.5);

        Int_t  neuIdx[nDirs]; // 0 - Plus; 1 - Minus;
        Int_t  regIdx[nDirs]; // 0 - Plus; 1 - Minus;
        Bool_t isTwoNueTrue[nDirs];
        memset(neuIdx, -1, sizeof(neuIdx));
        memset(regIdx, -1, sizeof(regIdx));
        memset(isTwoNueTrue, 0, sizeof(isTwoNueTrue));

        Double_t zdc = zdcPlus;
        for(Int_t idir = 0; idir < nDirs; idir++){
            if(idir == 1) zdc = zdcMinus;

            neuIdx[idir] = grabNeutronNum(idir, zdc);
            hZDCvsNeuNum[idir]->Fill(neuIdx[idir], zdc);

            regIdx[idir] = grabZDCRegionIdx(idir, zdc);
            hZDCvsRegIdx[idir]->Fill(regIdx[idir], zdc);

            if(zdc >= mTwoNeuZDCLow[idir] && zdc <= mTwoNeuZDCHi[idir]) isTwoNueTrue[idir] = 1;
        }

        hNeuNumMinusvsNeuNumPlus->Fill(neuIdx[0], neuIdx[1]);

        // Loop over the correct-sign candidates
        for (UInt_t icand = 0; icand < csTree.candSize(); icand++) {
            Float_t pt = csTree.pT()[icand];
            Float_t mass = csTree.mass()[icand];
            Float_t y = csTree.y()[icand];

            if (trigIdx == 4 && !csTree.trigCand(trigIdx, icand, false)) continue; // dimuon UPC
            if (trigIdx == 7 && !csTree.trigCand(trigIdx, icand, true))  continue; // single muon UPC

            if (!csTree.softCand(icand)) continue;
            //if (csTree.VtxProb()[icand] < mVtxProbCut) continue;
            if (!goodTrack(csTree, icand)) continue;
            if (TMath::Abs(y) > mPairYCut) continue;

            TVector3 muPlusMom, muMinusMom;
            if (csTree.chargeD1()[icand] > 0) {
                muPlusMom.SetPtEtaPhi(csTree.pTD1()[icand], csTree.EtaD1()[icand], csTree.PhiD1()[icand]);
                muMinusMom.SetPtEtaPhi(csTree.pTD2()[icand], csTree.EtaD2()[icand], csTree.PhiD2()[icand]);
            } else {
                muPlusMom.SetPtEtaPhi(csTree.pTD2()[icand], csTree.EtaD2()[icand], csTree.PhiD2()[icand]);
                muMinusMom.SetPtEtaPhi(csTree.pTD1()[icand], csTree.EtaD1()[icand], csTree.PhiD1()[icand]);
            }

            TVector3 reverseMuMinusMom = -muMinusMom;
            Double_t deltaPhi = shiftDeltaPhi(muPlusMom.Phi() - reverseMuMinusMom.Phi());
            Double_t deltaEta = muPlusMom.Eta() - reverseMuMinusMom.Eta();
            hDeltaEtavsDetaPhi->Fill(deltaPhi, deltaEta);

            Double_t asyPhi = 1 - TMath::Abs(shiftDeltaPhi(muPlusMom.DeltaPhi(muMinusMom))) / PI; //acoplanarity
            Double_t asyPt = TMath::Abs((muPlusMom.Pt() - muMinusMom.Pt()) / (muPlusMom.Pt() + muMinusMom.Pt()));

            hCS_MvsPtvsCen->Fill(cen, pt, mass);
            hCS_MvsPtvsRap->Fill(y, pt, mass);
            hCS_MvsAsyPhivsRap->Fill(y, asyPhi, mass);

            Double_t posPt     = csTree.chargeD1()[icand] > 0 ? csTree.pTD1()[icand]  : csTree.pTD2()[icand]; 
            Double_t posEta    = csTree.chargeD1()[icand] > 0 ? csTree.EtaD1()[icand] : csTree.EtaD2()[icand]; 
            Double_t posPhi    = csTree.chargeD1()[icand] > 0 ? csTree.PhiD1()[icand] : csTree.PhiD2()[icand]; 
            Bool_t   isPosTrig = csTree.chargeD1()[icand] > 0 ? csTree.trigMuon1()[trigIdx][icand] : csTree.trigMuon2()[trigIdx][icand];
            Double_t negPt     = csTree.chargeD1()[icand] < 0 ? csTree.pTD1()[icand]  : csTree.pTD2()[icand]; 
            Double_t negEta    = csTree.chargeD1()[icand] < 0 ? csTree.EtaD1()[icand] : csTree.EtaD2()[icand]; 
            Double_t negPhi    = csTree.chargeD1()[icand] < 0 ? csTree.PhiD1()[icand] : csTree.PhiD2()[icand]; 
            Bool_t   isNegTrig = csTree.chargeD1()[icand] < 0 ? csTree.trigMuon1()[trigIdx][icand] : csTree.trigMuon2()[trigIdx][icand];

            double_t posy      = etaToY(posEta, mass, posPt);
            double_t negy      = etaToY(negEta, mass, negPt);
            
            double_t costheta  = 0;
            if (posy > 0)  double costheta = tanh(0.5*(posy-negy));
            else double costheta = tanh(0.5*(negy-posy));

            double_t y_mumu = yyToy_mumu(posy, negy);

            //from the final state muon pair, calculate the photon energy
            double_t k_max = PhotonEnergy(y_mumu);
            double_t k_min = PhotonEnergy(-y_mumu);

            double_t SysM_mumu = m_mumu(posPt, posEta, posPhi, posy, negPt, negEta, negPhi, negy);
            
            if(
                    mass>massLow[0] && mass<massHi[nMBins-1]
                    && asyPhi<mAlphaCut
              ){
                hPosMuPhivsEtavsPt->Fill(posPt, posEta, posPhi);
                hNegMuPhivsEtavsPt->Fill(negPt, negEta, negPhi);

                // NOTE, in the pair selection, at least one triggered muon has already been deposited
                if(!isPosTrig){
                    hTnpDenEtavsPt->Fill(posPt, TMath::Abs(posEta));
                    hTnpFailEtavsPt->Fill(posPt, TMath::Abs(posEta));
                }
                if(!isNegTrig){
                    hTnpDenEtavsPt->Fill(negPt, TMath::Abs(negEta));
                    hTnpFailEtavsPt->Fill(negPt, TMath::Abs(negEta));
                }
                if(isPosTrig && isNegTrig){
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
            }

            Int_t posPtBin  = hPosMu3DMthEff->GetXaxis()->FindBin(posPt);
            Int_t posEtaBin = hPosMu3DMthEff->GetYaxis()->FindBin(posEta);
            Int_t posPhiBin = hPosMu3DMthEff->GetZaxis()->FindBin(posPhi);
            Int_t negPtBin  = hNegMu3DMthEff->GetXaxis()->FindBin(negPt);
            Int_t negEtaBin = hNegMu3DMthEff->GetYaxis()->FindBin(negEta);
            Int_t negPhiBin = hNegMu3DMthEff->GetZaxis()->FindBin(negPhi);

            if(posPtBin > hPosMu3DMthEff->GetNbinsX()) posPtBin = hPosMu3DMthEff->GetNbinsX();
            if(negPtBin > hNegMu3DMthEff->GetNbinsX()) negPtBin = hNegMu3DMthEff->GetNbinsX();

            Double_t totEff = 1;
            if(effCorr){
                Double_t trkEff = 1;
                Double_t trigEff = 1;

                trkEff  = hPosMu3DMthEff->GetBinContent(posPtBin, posEtaBin, posPhiBin) * hNegMu3DMthEff->GetBinContent(negPtBin, negEtaBin, negPhiBin);
                trigEff = 1 - (1 - hPosMu3DTrigEff->GetBinContent(posPtBin, posEtaBin, posPhiBin)) * (1 - hNegMu3DTrigEff->GetBinContent(negPtBin, negEtaBin, negPhiBin));
                

                totEff = trkEff * trigEff;
                //cout << "EffCorr: " << 1/totEff << endl;
            }

            hM_FineBin->Fill(mass, 1/totEff);
            hM_CoarseBin->Fill(mass, 1/totEff);

            hMvsPtvsCen->Fill(cen, pt, mass, 1/totEff);
            hMvsPtvsRap->Fill(y, pt, mass, 1/totEff);

            hPt2vsM->Fill(mass, pt*pt, 1/totEff);
            hAsyPhivsM->Fill(mass, asyPhi, 1/totEff);
            hAsyPtvsM->Fill(mass, asyPt, 1/totEff);

            hEvtvsCostheta->Fill(costheta, 1/totEff);
            hEvtvsM_mumu->Fill(SysM_mumu, 1/totEff);
            hEvtvsk_max->Fill(k_max, 1/totEff);
            hEvtvsk_min->Fill(k_min, 1/totEff);

            if(mass>massLow[0] && mass<massHi[nMBins-1]){
                hAsyPhivsRap->Fill(y, asyPhi, 1/totEff);
            }

            if(neuIdx[0] >= 0 && neuIdx[1] >= 0){
                hAsyPhivsM_Neu[neuIdx[0]][neuIdx[1]]->Fill(mass, asyPhi, 1/totEff);
                if(mass>massLow[0] && mass<massHi[nMBins-1]){
                    hAsyPhivsRap_Neu[neuIdx[0]][neuIdx[1]]->Fill(y, asyPhi, 1/totEff);
                }

                for(Int_t idir = 0; idir < nDirs; idir++){
                    if(isTwoNueTrue[idir]){
                        hAsyPhivsM_Neu2[neuIdx[1-idir]][idir]->Fill(mass, asyPhi, 1/totEff);;
                        if(mass>massLow[0] && mass<massHi[nMBins-1]){
                            hAsyPhivsRap_Neu2[neuIdx[1-idir]][idir]->Fill(y, asyPhi, 1/totEff);
                        }
                    }
                }

                hRapvsNMinusNeuvsNPlusNeu->Fill(neuIdx[0], neuIdx[1], y, 1/totEff);
            }

            //TVector3 muMomDiff;
            //if(muPlusMom.Mag()>muMinusMom.Mag()){
            //    muMomDiff = muPlusMom - muMinusMom;
            //}
            //else{
            //    muMomDiff = muMinusMom - muPlusMom;
            //}

            TVector3 muMomDiff = muPlusMom - muMinusMom;

            TVector3 pairMom;
            pairMom.SetPtEtaPhi(csTree.pT()[icand], csTree.eta()[icand], csTree.phi()[icand]);

            Double_t phiDiff = shiftDeltaPhi(pairMom.DeltaPhi(muMomDiff));
            //Double_t phiDiff = shiftToPi(pairMom.DeltaPhi(muMomDiff));

            if(asyPhi < mAlphaCut) hDeltaPhivsM->Fill(mass, phiDiff, 1/totEff);
            if(mass>massLow[0] && mass<massHi[nMBins-1]) hDeltaPhivsPt->Fill(pt, phiDiff, 1/totEff);
        }

        // Loop over the wrong-sign candidates
        if (wsTree.GetEntry(jentry) < 0) {
            cout << "Invalid wrong-sign entry!" << endl;
            return;
        }

        for (UInt_t icand = 0; icand < wsTree.candSize(); icand++) {
            Float_t pt = wsTree.pT()[icand];
            Float_t mass = wsTree.mass()[icand];
            Float_t y = wsTree.y()[icand];

            if (trigIdx == 4 && !wsTree.trigCand(trigIdx, icand, false)) continue; // dimuon UPC
            if (trigIdx == 7 && !wsTree.trigCand(trigIdx, icand, true))  continue; // single muon UPC

            if (!wsTree.softCand(icand)) continue;
            //if (wsTree.VtxProb()[icand] < mVtxProbCut) continue;
            if (!goodTrack(wsTree, icand)) continue;
            if (TMath::Abs(y) > mPairYCut) continue;

            hWS_MvsPtvsCen->Fill(cen, pt, mass);
            hWS_MvsPtvsRap->Fill(y, pt, mass);
        }
    }

    //cout<<endl;
    //cout<<"Analyzed runs: "<<endl;
    //sort(anaRuns.begin(), anaRuns.end());
    //for(auto anarun : anaRuns){
    //    cout << anarun << endl;
    //}

    TString dirName = "dimuonHistos";
    system(Form("mkdir -p %s", dirName.Data()));

    TString fileName = "";
    if(effCorr){
        fileName += "effCorr";

        if(applyTnPSF){
            fileName += ".applyTnPSF";
        }
    }
    else{
        fileName += "rawSig";
    }

    if(incHadron){
        fileName += ".incHadron";
    }

    if(hfVetoType.EqualTo("Tight")){
        fileName += ".tightHF";
    }

    if(hfVetoType.EqualTo("removeHF")){
        fileName += ".removeHF";
    }

    cout<<"fileName: "<<fileName<<endl;

    writeHistos(Form("%s/%s", dirName.Data(), fileName.Data()));
}

void bookHistos()
{
    const Int_t mHistRapBins = 4;
    const Double_t mHistRapLow = -2.4;
    const Double_t mHistRapHi = 2.4;
    const Int_t mHistCenBins = 20;
    const Double_t mHistCenLow = 0;
    const Double_t mHistCenHi = 200;
    const Int_t mHistPtBins = 200;
    const Double_t mHistPtLow = 0;
    const Double_t mHistPtHi = 1;
    const Int_t mHistMassBins = 2000;
    const Double_t mHistMassLow = 0;
    const Double_t mHistMassHi = 100;
    const Int_t mHistZdcBins = 500;
    const Double_t mHistZdcLow = 0;
    const Double_t mHistZdcHi = 5.e4;
    const Int_t mHistPt2Bins = 5000;
    const Double_t mHistPt2Low = 0;
    const Double_t mHistPt2Hi = 1.25;
    const Int_t mHistAsyPhiBins = 5000;
    const Double_t mHistAsyPhiLow = 0;
    const Double_t mHistAsyPhiHi = 0.5;
    const Int_t mHistAsyPtBins = 1000;
    const Double_t mHistAsyPtLow = 0;
    const Double_t mHistAsyPtHi = 1;

    // event level
    hnEvts = new TH1D("hnEvts", "hnEvts;", 5, 0, 5);
    hnEvts->GetXaxis()->SetBinLabel(1, "trigEvt");
    hnEvts->GetXaxis()->SetBinLabel(2, "!Beam-halo & GoodVtx");
    hnEvts->GetXaxis()->SetBinLabel(3, "HFMaxE <= 7.6(7.3) GeV");
    hnEvts->GetXaxis()->SetBinLabel(4, "N_{trk}^{HP} == 2");
    hnEvts->GetXaxis()->LabelsOption("d");
    hnEvts->GetXaxis()->SetLabelSize(0.055);

    hCenvsTrig = new TH2D("hCenvsTrig", "hCenvsTrig; trigIdx; Centrality", 8, -0.5, 7.5, 200, 0, 200);
    hPreScalevsTrig = new TH2D("hPreScalevsTrig", "hPreScalevsTrig; trigIdx; pre-Scale", 8, -0.5, 7.5, 20, 0, 20);

    hVzvsVyvsVx = new TH3D("hVzvsVyvsVx", "hVzvsVyvsVx; V_{x} (cm); V_{y} (cm); V_{z} (cm)", 100, -1, 1, 100, -1, 1, 300, -30, 30);
    hVzvsVyvsVx_Sel = new TH3D("hVzvsVyvsVx_Sel", "hVzvsVyvsVx_Sel; V_{x} (cm); V_{y} (cm); V_{z} (cm)", 100, -1, 1, 100, -1, 1, 300, -30, 30);

    hRawCen = new TH1D("hRawCen", "hRawCen; Centrality", 200, 0, 200);
    hCen_Final = new TH1D("hCen_Final", "hCen_Final; Centrality", 200, 0, 200);

    hRawZDCMinusvsZDCPlus = new TH2D("hRawZDCMinusvsZDCPlus", "hRawZDCMinusvsZDCPlus; ZDCPlus; ZDCMinus", 5000, 0, 5.e5, 5000, 0, 5.e5);
    hZDCMinusvsZDCPlus = new TH2D("hZDCMinusvsZDCPlus", "hZDCMinusvsZDCPlus; ZDCPlus; ZDCMinus", mHistZdcBins*4, mHistZdcLow, mHistZdcHi*4, mHistZdcBins*4, mHistZdcLow, mHistZdcHi*4);
    hZDCMinusvsZDCPlus_Sel = new TH2D("hZDCMinusvsZDCPlus_Sel", "hZDCMinusvsZDCPlus_Sel; ZDCPlus; ZDCMinus", mHistZdcBins*4, mHistZdcLow, mHistZdcHi*4, mHistZdcBins*4, mHistZdcLow, mHistZdcHi*4);
    hZDCMinusvsZDCPlus_Only2MuTrk = new TH2D("hZDCMinusvsZDCPlus_Only2MuTrk", "hZDCMinusvsZDCPlus_Only2MuTrk; ZDCPlus; ZDCMinus", mHistZdcBins*4, mHistZdcLow, mHistZdcHi*4, mHistZdcBins*4, mHistZdcLow, mHistZdcHi*4);

    hZDCMinusvsZDCPlus_LS = new TH2D("hZDCMinusvsZDCPlus_LS", "hZDCMinusvsZDCPlus_LS; ZDCPlus; ZDCMinus", 500, 0, 5.e6, 500, 0, 5.e6);
    hZDCMinusvsZDCPlus_Sel_LS = new TH2D("hZDCMinusvsZDCPlus_Sel_LS", "hZDCMinusvsZDCPlus_Sel_LS; ZDCPlus; ZDCMinus", 500, 0, 5.e6, 500, 0, 5.e6);
    hZDCMinusvsZDCPlus_Only2MuTrk_LS = new TH2D("hZDCMinusvsZDCPlus_Only2MuTrk_LS", "hZDCMinusvsZDCPlus_Only2MuTrk_LS; ZDCPlus; ZDCMinus", 500, 0, 5.e6, 500, 0, 5.e6);

    hHFMinusvsHFPlusvsCen = new TH3D("hHFMinusvsHFPlusvsCen", "hHFMinusvsHFPlusvsCen; Centrality; HFsumETPlus; HFsumETMinus", 200, 0, 200, 200, 0, 100, 200, 0, 100);
    hHFMinusvsHFPlusvsCen_Sel = new TH3D("hHFMinusvsHFPlusvsCen_Sel", "hHFMinusvsHFPlusvsCen_Sel; Centrality; HFsumETPlus; HFsumETMinus", 200, 0, 200, 200, 0, 100, 200, 0, 100);

    hFlagvsBit = new TH2D("hFlagvsBit", "hFlagvsBit; Bit; Flag", 20, -0.5, 20-0.5, 2, 0, 2);

    hNtrkHPvsNtrkofflinevsCen = new TH3D("hNtrkHPvsNtrkofflinevsCen", "hNtrkHPvsNtrkofflinevsCen; Centrality; N_{trk}^{offline}; N_{trk}^{HP}", 200, 0, 200, 100, 0, 100, 100, 0, 100);
    hNtrkHPvsNtrkofflinevsCen_Sel = new TH3D("hNtrkHPvsNtrkofflinevsCen_Sel", "hNtrkHPvsNtrkofflinevsCen_Sel; Centrality; N_{trk}^{offline}; N_{trk}^{HP}", 200, 0, 200, 100, 0, 100, 100, 0, 100);
    hNtrkHP_2SoftMuons = new TH1D("hNtrkHP_2SoftMuons", "hNtrkHP_2SoftMuons; N_{trk}^{HP}", 100, 0, 100);

    hEvtvsCostheta = new TH1D("hEvtvsCostheta", "hEvtvsCostheta; cos #theta", 100, -1, 1);
    hEvtvsM_mumu = new TH1D("hEvtvsM_mumu", "hEvtvsM_mumu; M_{#mu#mu} (GeV/c^{2})", 52, 8, 60);
    hEvtvsk_max = new TH1D("hEvtvsk_max", "hEvtvsk_max; k_{max} (GeV)", 52, 8, 60);
    hEvtvsk_min = new TH1D("hEvtvsk_min", "hEvtvsk_min; k_{min} (GeV)", 100, 0, 100);

    for(Int_t idir = 0; idir < nDirs; idir++){
        hZDCvsNeuNum[idir] = new TH2D(Form("hZDC%svsNeuNum", mDir[idir].Data()), Form("; # of Neutruon; ZDC%s", mDir[idir].Data()), 5, -0.5, 4.5, mHistZdcBins*4, mHistZdcLow*4, mHistZdcHi*4);

        hZDCvsRegIdx[idir] = new TH2D(Form("hZDC%svsRegIdx", mDir[idir].Data()), Form("; ZDC Region Index; ZDC%s", mDir[idir].Data()), 10, -0.5, 9.5, mHistZdcBins*4, mHistZdcLow*4, mHistZdcHi*4);
    }

    hNeuNumMinusvsNeuNumPlus = new TH2D("hNeuNumMinusvsNeuNumPlus", "hNeuNumMinusvsNeuNumPlus; NeuNumPlus; NeuNumMinus", 3, 0, 3, 3, 0, 3);

    hDeltaEtavsDetaPhi = new TH2D("hDeltaEtavsDetaPhi", "hDeltaEtavsDetaPhi; #Delta#phi; #Delta#eta", 180, -PI, PI, 500, -5, 5);

    hCS_MvsPtvsCen = new TH3D("hCS_MvsPtvsCen", "hCS_MvsPtvsCen; Centrality; p_{T} (GeV/c); M_{#mu#mu} (GeV/c^{2})", mHistCenBins, mHistCenLow, mHistCenHi, mHistPtBins, mHistPtLow, mHistPtHi, mHistMassBins, mHistMassLow, mHistMassHi);
    hCS_MvsPtvsRap = new TH3D("hCS_MvsPtvsRap", "hCS_MvsPtvsRap; Rapidity; p_{T} (GeV/c); M_{#mu#mu} (GeV/c^{2})", 24, -2.4, 2.4, mHistPtBins, mHistPtLow, mHistPtHi, mHistMassBins, mHistMassLow, mHistMassHi);
    hCS_MvsAsyPhivsRap = new TH3D("hCS_MvsAsyPhivsRap", "hCS_MvsAsyPhivsRap; Rapidity; #alpha; M_{#mu#mu} (GeV/c^{2})", 24, -2.4, 2.4, 500, 0, 0.05, mHistMassBins, mHistMassLow, mHistMassHi);
    hWS_MvsPtvsCen = new TH3D("hWS_MvsPtvsCen", "hWS_MvsPtvsCen; Centrality; p_{T} (GeV/c); M_{#mu#mu} (GeV/c^{2})", mHistCenBins, mHistCenLow, mHistCenHi, mHistPtBins, mHistPtLow, mHistPtHi, mHistMassBins, mHistMassLow, mHistMassHi);
    hWS_MvsPtvsRap = new TH3D("hWS_MvsPtvsRap", "hWS_MvsPtvsRap; Rapidity; p_{T} (GeV/c); M_{#mu#mu} (GeV/c^{2})", 24, -2.4, 2.4, mHistPtBins, mHistPtLow, mHistPtHi, mHistMassBins, mHistMassLow, mHistMassHi);
    hMvsPtvsCen = new TH3D("hMvsPtvsCen", "hMvsPtvsCen; Centrality; p_{T} (GeV/c); M_{#mu#mu} (GeV/c^{2})", mHistCenBins, mHistCenLow, mHistCenHi, mHistPtBins, mHistPtLow, mHistPtHi, mHistMassBins, mHistMassLow, mHistMassHi);
    hMvsPtvsRap = new TH3D("hMvsPtvsRap", "hMvsPtvsRap; Rapidity; p_{T} (GeV/c); M_{#mu#mu} (GeV/c^{2})", 24, -2.4, 2.4, mHistPtBins, mHistPtLow, mHistPtHi, mHistMassBins, mHistMassLow, mHistMassHi);

    hPosMuPhivsEtavsPt = new TH3D("hPosMuPhivsEtavsPt", "hPosMuPhivsEtavsPt; p_{T} (GeV/c); #eta; #phi", 100, 0, 50, 24, -2.4, 2.4, 60, -PI, PI);
    hNegMuPhivsEtavsPt = new TH3D("hNegMuPhivsEtavsPt", "hNegMuPhivsEtavsPt; p_{T} (GeV/c); #eta; #phi", 100, 0, 50, 24, -2.4, 2.4, 60, -PI, PI);

    hTnpDenEtavsPt = new TH2D("hTnpDenEtavsPt", "hTnpDenEtavsPt; p_{T} (GeV/c); #eta", 300, 0, 30, 24, 0, 2.4);
    hTnpPassEtavsPt = new TH2D("hTnpPassEtavsPt", "hTnpPassEtavsPt; p_{T} (GeV/c); #eta", 300, 0, 30, 24, 0, 2.4);
    hTnpFailEtavsPt = new TH2D("hTnpFailEtavsPt", "hTnpFailEtavsPt; p_{T} (GeV/c); #eta", 300, 0, 30, 24, 0, 2.4);

    hM_FineBin   = new TH1D("hM_FineBin", "hM_FineBin; M_{#mu#mu}", mHistMassBins, mHistMassLow, mHistMassHi);
    hM_CoarseBin = new TH1D("hM_CoarseBin", "hM_CoarseBin; M_{#mu#mu}", 100, mHistMassLow, mHistMassHi);

    hPt2vsM = new TH2D("hPt2vsM", "hPt2vsM; M_{#mu#mu} (GeV/c^{2}); p_{T}^{2} ((GeV/c)^{2})", mHistMassBins, mHistMassLow, mHistMassHi, mHistPt2Bins, mHistPt2Low, mHistPt2Hi);
    hAsyPhivsM = new TH2D("hAsyPhivsM", "hAsyPhivsM; M_{#mu#mu} (GeV/c^{2}); #alpha", mHistMassBins, mHistMassLow, mHistMassHi, mHistAsyPhiBins, mHistAsyPhiLow, mHistAsyPhiHi);
    hAsyPhivsRap = new TH2D("hAsyPhivsRap", "hAsyPhivsRap; Rapidity; #alpha", mHistRapBins, mHistRapLow, mHistRapHi, mHistAsyPhiBins, mHistAsyPhiLow, mHistAsyPhiHi);
    hAsyPtvsM = new TH2D("hAsyPtvsM", "hAsyPtvsM; M_{#mu#mu} (GeV/c^{2}); p_{T} Asymmetry", mHistMassBins, mHistMassLow, mHistMassHi, mHistAsyPtBins, mHistAsyPtLow, mHistAsyPtHi);
    hDeltaPhivsM = new TH2D("hDeltaPhivsM", "hDeltaPhivsM; M_{#mu#mu} (GeV/c^{2}); #phi_{#mu^{+}+#mu^{-}} - #phi_{#mu^{+}-#mu^{-}}", mHistMassBins, mHistMassLow, mHistMassHi, 120, -PI, PI);
    hDeltaPhivsPt = new TH2D("hDeltaPhivsPt", "hDeltaPhivsPt; p_{T} (GeV/c); #phi_{#mu^{+}+#mu^{-}} - #phi_{#mu^{+}-#mu^{-}}", mHistPtBins*5, mHistPtLow, mHistPtHi*5, 120, -PI, PI);

    hRapvsNMinusNeuvsNPlusNeu = new TH3D("hRapvsNMinusNeuvsNPlusNeu", "hRapvsNMinusNeuvsNPlusNeu; nNeu (Plus); nNeu (Minus); Rapidity", 5, 0, 5, 5, 0, 5, 48, -2.4, 2.4);

    for(Int_t ip=0; ip<nNeus; ip++){
        for(Int_t im=0; im<nNeus; im++){
            hAsyPhivsM_Neu[ip][im] = new TH2D(Form("hAsyPhivsM_Neu%dp%dm", ip, im), "; M_{#mu#mu} (GeV/c^{2}); #alpha", mHistMassBins, mHistMassLow, mHistMassHi, mHistAsyPhiBins, mHistAsyPhiLow, mHistAsyPhiHi);
            hAsyPhivsRap_Neu[ip][im] = new TH2D(Form("hAsyPhivsRap_Neu%dp%dm", ip, im), "; Rapidity; #alpha", mHistRapBins, mHistRapLow, mHistRapHi, mHistAsyPhiBins, mHistAsyPhiLow, mHistAsyPhiHi);
        }
    }

    for(Int_t ineu=0; ineu<nNeus; ineu++){
        for(Int_t idir=0; idir<nDirs; idir++){
            hAsyPhivsM_Neu2[ineu][idir] = new TH2D(Form("hAsyPhivsM_Neu%dnExact2%s", ineu, mDir[idir].Data()), "; M_{#mu#mu} (GeV/c^{2}); #alpha", mHistMassBins, mHistMassLow, mHistMassHi, mHistAsyPhiBins, mHistAsyPhiLow, mHistAsyPhiHi);
            hAsyPhivsRap_Neu2[ineu][idir] = new TH2D(Form("hAsyPhivsRap_Neu%dnExact2%s", ineu, mDir[idir].Data()), "; Rapidity; #alpha", mHistRapBins, mHistRapLow, mHistRapHi, mHistAsyPhiBins, mHistAsyPhiLow, mHistAsyPhiHi);
        }
    }
}

void writeHistos(TString fileName)
{
    TFile* fOut = new TFile(Form("%s.root", fileName.Data()), "recreate");

    fOut->cd();

    for(Int_t ineu=0; ineu<nNeus; ineu++){
        for(Int_t idir=0; idir<nDirs; idir++){
            hAsyPhivsM_Neu2[ineu][idir]->Write();
            hAsyPhivsRap_Neu2[ineu][idir]->Write();
        }
    }

    for(Int_t ip=0; ip<nNeus; ip++){
        for(Int_t im=0; im<nNeus; im++){
            hAsyPhivsM_Neu[ip][im]->Write();
            hAsyPhivsRap_Neu[ip][im]->Write();
        }
    }
    hRapvsNMinusNeuvsNPlusNeu->Write();

    hnEvts->Write();
    hCenvsTrig->Write();
    hPreScalevsTrig->Write();

    hVzvsVyvsVx->Write();
    hVzvsVyvsVx_Sel->Write();

    hRawCen->Write();
    hCen_Final->Write();

    hHFMinusvsHFPlusvsCen->Write();
    hHFMinusvsHFPlusvsCen_Sel->Write();

    hFlagvsBit->Write();

    hRawZDCMinusvsZDCPlus->Write();
    hZDCMinusvsZDCPlus->Write();
    hZDCMinusvsZDCPlus_Sel->Write();
    hZDCMinusvsZDCPlus_Only2MuTrk->Write();

    hZDCMinusvsZDCPlus_LS->Write();
    hZDCMinusvsZDCPlus_Sel_LS->Write();
    hZDCMinusvsZDCPlus_Only2MuTrk_LS->Write();

    hNtrkHPvsNtrkofflinevsCen->Write();
    hNtrkHPvsNtrkofflinevsCen_Sel->Write();
    hNtrkHP_2SoftMuons->Write();

    for(Int_t idir = 0; idir < nDirs; idir++){
        hZDCvsNeuNum[idir]->Write();
        hZDCvsRegIdx[idir]->Write();
    }

    hNeuNumMinusvsNeuNumPlus->Write();

    hDeltaEtavsDetaPhi->Write();

    hCS_MvsPtvsCen->Write();
    hCS_MvsPtvsRap->Write();
    hCS_MvsAsyPhivsRap->Write();
    hWS_MvsPtvsCen->Write();
    hWS_MvsPtvsRap->Write();

    hMvsPtvsCen->Write();
    hMvsPtvsRap->Write();

    hPosMuPhivsEtavsPt->Write();
    hNegMuPhivsEtavsPt->Write();

    hTnpDenEtavsPt->Write();
    hTnpPassEtavsPt->Write();
    hTnpFailEtavsPt->Write();

    hM_FineBin->Write();
    hM_CoarseBin->Write();

    hPt2vsM->Write();
    hAsyPhivsM->Write();
    hAsyPhivsRap->Write();
    hAsyPtvsM->Write();
    hDeltaPhivsM->Write();
    hDeltaPhivsPt->Write();

    hEvtvsCostheta->Write();
    hEvtvsM_mumu->Write();
    hEvtvsk_max->Write();
    hEvtvsk_min->Write();

    fOut->Close();
}

Bool_t goodTrack(VertexCompositeTree& evtTree, const int icand)
{
    if (evtTree.pTD1()[icand] < mPtCut || evtTree.pTD2()[icand] < mPtCut)
        return kFALSE;
    if (TMath::Abs(evtTree.EtaD1()[icand]) > mEtaCut || TMath::Abs(evtTree.EtaD2()[icand]) > mEtaCut)
        return kFALSE;

    return kTRUE;
}

Int_t grabNeutronNum(Int_t dirIdx, Double_t ranNum, Double_t prob[][nNeus])
{
    Double_t probTh = 0;
    for(Int_t ineu = 0; ineu < nNeus; ineu++){
        probTh += prob[dirIdx][ineu];

        if(ranNum < probTh) return ineu;
    }

    return -1;
}

Int_t grabNeutronNum(Int_t dirIdx, Double_t zdc)
{
    for(Int_t ineu = 0; ineu < nNeus; ineu++){
        if(zdc > mNeuZDCLow[dirIdx][ineu] && zdc <= mNeuZDCHi[dirIdx][ineu]) return ineu;
    }

    return -1;
}

Int_t grabZDCRegionIdx(Int_t dirIdx, Double_t zdc)
{
    for(Int_t ireg = 0; ireg < nRegions; ireg++){
        if(zdc > mNeuZDC[dirIdx][ireg] && zdc <= mNeuZDC[dirIdx][ireg+1]) return ireg;
    }

    return -1;
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

// Bool_t loadUpcTnPFile()
// {
//     Bool_t loadFlag = kTRUE;

//     //fUpcTnp = TFile::Open("/Users/syang/work/run2/upcDimuon/TnPStudy/plots_DatavsMC_woNtrkHPSel/upcTnP.root");
//     fUpcTnp = TFile::Open("/Users/syang/work/run2/upcDimuon/TnPStudy/plots_DatavsMC_wNtrkHPSel/upcTnP.root");
//     if(!fUpcTnp){
//         loadFlag = kFALSE;
//     }
//     else{
//         for(Int_t ieta=0; ieta<nEtaBins; ieta++){
//             hSFStat_MuIdTnp[ieta] = (TH1D *)fUpcTnp->Get(Form("hSFStat_MuIdTnp_EtaBin%d", ieta));
//             hSFStat_TrigTnp[ieta] = (TH1D *)fUpcTnp->Get(Form("hSFStat_TrigTnp_EtaBin%d", ieta));
//         }
//     }

//     return loadFlag;
// }

Int_t grabEtaIdx(Double_t eta)
{
    Int_t etaIdx = -1;
    for(Int_t ieta=0; ieta<nEtaBins; ieta++){
        if(TMath::Abs(eta)>=etaBoundary[ieta] && TMath::Abs(eta)<etaBoundary[ieta+1]){
            etaIdx = ieta;
            break;
        }
    }

    return etaIdx;
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

double_t PhotonEnergy(double_t y_mu_mu) 
{
    double_t m_mumu_const = 211.32 / 1000.0;  // 缪子对的静止质量 (GeV/c^2)
    return 0.5 * m_mumu_const * exp(y_mu_mu); // 计算 k
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

