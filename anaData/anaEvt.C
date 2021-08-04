#include "../common/headers.h"
#include "../common/funUtil.h"
#include "../common/VertexCompositeTree.h"

TH1D *hnEvts;
TH2D *hCenvsTrig;
TH2D *hPreScalevsTrig;

TH3D *hVzvsVyvsVx;
TH3D *hVzvsVyvsVx_Sel;

TH1D *hRawCen;
TH1D *hCen_Final;

TH3D *hHFMinusvsHFPlusvsCen;
TH3D *hHFMinusvsHFPlusvsCen_Sel;

TH3D *hNtrkHPvsNtrkofflinevsCen;
TH3D *hNtrkHPvsNtrkofflinevsCen_Sel;
TH1D *hNtrkHP_2SoftMuons;

TH2D *hZDCMinusvsZDCPlus;
TH2D *hZDCMinusvsZDCPlus_Sel;
TH2D *hZDCMinusvsZDCPlus_Only2MuTrk;

TH2D *hZDCMinusvsZDCPlus_LS;
TH2D *hZDCMinusvsZDCPlus_Sel_LS;
TH2D *hZDCMinusvsZDCPlus_Only2MuTrk_LS;

TH2D *hZDCvsNeuNum[nDirs];
TH2D *hNeuNumMinusvsNeuNumPlus;

// without muon acceptance selection
TH3D *hMuPtvsEtavsRap;
TH3D *hTrigMuPtvsEtavsRap;
TH3D *hNegMuPtvsPosMuPtvsRap;
TH3D *hNegMuEtavsPosMuEtavsRap;
TH3D *hDeltaPtvsDeltaEtavsRap;

// check cosmic ray contamination
TH2D *hDeltaEtavsDetaPhi4Cosmic;

TH3D *hMvsPtvsRap;
TH3D *hMvsAsyPhivsRap;
TH3D *hMvsPtvsRap_WS;
TH3D *hMvsPtvsRap_NeuDir[nNeus][nNeus];
TH3D *hMvsAsyPhivsRap_NeuDir[nNeus][nNeus];

Bool_t   goodMuPair(VertexCompositeTree& evtTree, const int icand);
Int_t    grabNeutronNum(Int_t dirIdx, Double_t zdc);
Double_t shiftDeltaPhi(Double_t dPhi);
Double_t shiftToPi(Double_t dPhi);

void bookHistos();
void writeHistos(TString fileName = "test");

void anaEvt(Bool_t incHadron = kFALSE, TString hfVetoType="Default")
{
    TH1::SetDefaultSumw2(kTRUE);

    if(!hfVetoType.EqualTo("Default") && !hfVetoType.EqualTo("Tight") && !hfVetoType.EqualTo("Loose") && !hfVetoType.EqualTo("removeHF")){
        cout<<"Please input the correct hfVetoType string: 'Default' OR 'Tight' OR 'Loose' OR 'removeHF'"<<endl;
        return;
    }

    const auto& inputFile = "../rootfiles/VertexCompositeTree_HIForward_HIRun2018_04Apr2019_DiMuMassMin0_20191206.root";

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

    bookHistos();

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
        hVzvsVyvsVx->Fill(vtxX, vtxY, vtxZ);

        // Select Event - require this event has a valid vertex and is not beam-halo event
        if(!csTree.evtSel()[2] || !csTree.evtSel()[3]) continue;
        hnEvts->Fill(1.5);

        hVzvsVyvsVx_Sel->Fill(vtxX, vtxY, vtxZ);

        hZDCMinusvsZDCPlus->Fill(zdcPlus, zdcMinus);
        hZDCMinusvsZDCPlus_LS->Fill(zdcPlus, zdcMinus);

        hHFMinusvsHFPlusvsCen->Fill(cen, hfsumETPlus, hfsumETMinus);
        hNtrkHPvsNtrkofflinevsCen->Fill(cen, nTrkoffline, nTrkHP);

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
        if(hfVetoType.EqualTo("Loose"))    hfVeto = !csTree.evtSel()[14] && !csTree.evtSel()[15];   // leading HF Energy < 8 GeV (Plus), < 8 GeV (Minus)
        if(hfVetoType.EqualTo("removeHF")) hfVeto = kTRUE; // remove HF veto
        if(!hfVeto) continue; 
        hnEvts->Fill(2.5);

        hZDCMinusvsZDCPlus_Sel->Fill(zdcPlus, zdcMinus);
        hZDCMinusvsZDCPlus_Sel_LS->Fill(zdcPlus, zdcMinus);
        hHFMinusvsHFPlusvsCen_Sel->Fill(cen, hfsumETPlus, hfsumETMinus);
        hNtrkHPvsNtrkofflinevsCen_Sel->Fill(cen, nTrkoffline, nTrkHP);

        Int_t nSoftMuonPair = 0;
        for (UInt_t icand = 0; icand < csTree.candSize(); icand++) {
            if (!csTree.softCand(icand)) continue;
            nSoftMuonPair++;
        }

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
        memset(neuIdx, -1, sizeof(neuIdx));

        Double_t zdc = zdcPlus;
        for(Int_t idir = 0; idir < nDirs; idir++){
            if(idir == 1) zdc = zdcMinus;

            neuIdx[idir] = grabNeutronNum(idir, zdc);
            hZDCvsNeuNum[idir]->Fill(neuIdx[idir], zdc);
        }

        hNeuNumMinusvsNeuNumPlus->Fill(neuIdx[0], neuIdx[1]);

        // Loop over the correct-sign candidates
        for (UInt_t icand = 0; icand < csTree.candSize(); icand++) {
            Float_t pt = csTree.pT()[icand];
            Float_t mass = csTree.mass()[icand];
            Float_t y = csTree.y()[icand];

            Double_t posPt     = csTree.chargeD1()[icand] > 0 ? csTree.pTD1()[icand]  : csTree.pTD2()[icand]; 
            Double_t posEta    = csTree.chargeD1()[icand] > 0 ? csTree.EtaD1()[icand] : csTree.EtaD2()[icand]; 
            Double_t posPhi    = csTree.chargeD1()[icand] > 0 ? csTree.PhiD1()[icand] : csTree.PhiD2()[icand]; 
            Bool_t   isPosTrig = csTree.chargeD1()[icand] > 0 ? csTree.trigMuon1()[trigIdx][icand] : csTree.trigMuon2()[trigIdx][icand];
            Double_t negPt     = csTree.chargeD1()[icand] < 0 ? csTree.pTD1()[icand]  : csTree.pTD2()[icand]; 
            Double_t negEta    = csTree.chargeD1()[icand] < 0 ? csTree.EtaD1()[icand] : csTree.EtaD2()[icand]; 
            Double_t negPhi    = csTree.chargeD1()[icand] < 0 ? csTree.PhiD1()[icand] : csTree.PhiD2()[icand]; 
            Bool_t   isNegTrig = csTree.chargeD1()[icand] < 0 ? csTree.trigMuon1()[trigIdx][icand] : csTree.trigMuon2()[trigIdx][icand];

            // soft muon selection
            if(!csTree.softCand(icand))   continue;

            // muon acceptance study in UPC Jpsi analysis
            if(mass>mMassLow4MuonAccStudy && mass<mMassHi4MuonAccStudy){
                hMuPtvsEtavsRap->Fill(y, posEta, posPt);
                hMuPtvsEtavsRap->Fill(y, negEta, negPt);

                if(isPosTrig) hTrigMuPtvsEtavsRap->Fill(y, posEta, posPt);
                if(isNegTrig) hTrigMuPtvsEtavsRap->Fill(y, negEta, negPt);

                Double_t deltaPt  = posPt - negPt;
                Double_t deltaEta = posEta - negEta;
                hNegMuPtvsPosMuPtvsRap->Fill(y, posPt, negPt);
                hNegMuEtavsPosMuEtavsRap->Fill(y, posEta, negEta);
                hDeltaPtvsDeltaEtavsRap->Fill(y, deltaEta, deltaPt);
            }

            // apply muon acceptance and trigger selections
            if(!goodMuPair(csTree, icand)) continue; //

            TVector3 muPlusMom, muMinusMom;
            if(csTree.chargeD1()[icand] > 0) {
                muPlusMom.SetPtEtaPhi(csTree.pTD1()[icand], csTree.EtaD1()[icand], csTree.PhiD1()[icand]);
                muMinusMom.SetPtEtaPhi(csTree.pTD2()[icand], csTree.EtaD2()[icand], csTree.PhiD2()[icand]);
            } 
            else{
                muPlusMom.SetPtEtaPhi(csTree.pTD2()[icand], csTree.EtaD2()[icand], csTree.PhiD2()[icand]);
                muMinusMom.SetPtEtaPhi(csTree.pTD1()[icand], csTree.EtaD1()[icand], csTree.PhiD1()[icand]);
            }

            TVector3 reversedMuMinusMom = -muMinusMom;
            Double_t deltaPhi4Cosmic = shiftDeltaPhi(muPlusMom.Phi() - reversedMuMinusMom.Phi());
            Double_t deltaEta4Cosmic = muPlusMom.Eta() - reversedMuMinusMom.Eta();
            hDeltaEtavsDetaPhi4Cosmic->Fill(deltaPhi4Cosmic, deltaEta4Cosmic);

            Double_t asyPhi = 1 - TMath::Abs(shiftDeltaPhi(muPlusMom.DeltaPhi(muMinusMom))) / PI; //acoplanarity

            hMvsPtvsRap->Fill(y, pt, mass);
            hMvsAsyPhivsRap->Fill(y, asyPhi, mass);

            if(neuIdx[0] >= 0 && neuIdx[1] >= 0){
                hMvsPtvsRap_NeuDir[neuIdx[0]][neuIdx[1]]->Fill(y, pt, mass);
                hMvsAsyPhivsRap_NeuDir[neuIdx[0]][neuIdx[1]]->Fill(y, asyPhi, mass);
            }
            else{
                cout<<"neutron index is wrong ---> "<<"plusIdx: "<<neuIdx[0]<<"  minusIdx:"<<neuIdx[1]<<endl;
            }
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

            if(!wsTree.softCand(icand))   continue;
            if(!goodMuPair(wsTree, icand)) continue;

            hMvsPtvsRap_WS->Fill(y, pt, mass);
        }
    }

    //cout<<endl;
    //cout<<"Analyzed runs: "<<endl;
    //sort(anaRuns.begin(), anaRuns.end());
    //for(auto anarun : anaRuns){
    //    cout << anarun << endl;
    //}

    TString dirName = "jpsiHistos";
    system(Form("mkdir -p %s", dirName.Data()));

    TString fileName = "rawSig";

    if(incHadron){
        fileName += ".incHadron";
    }

    if(hfVetoType.EqualTo("Tight")){
        fileName += ".tightHF";
    }
    else if(hfVetoType.EqualTo("Loose")){
        fileName += ".looseHF";
    }
    else if(hfVetoType.EqualTo("removeHF")){
        fileName += ".removeHF";
    }

    cout<<"fileName: "<<fileName<<endl;

    writeHistos(Form("%s/%s", dirName.Data(), fileName.Data()));
}

void bookHistos()
{
    const Int_t    mHistRapBins = 50;
    const Double_t mHistRapLow = -2.5;
    const Double_t mHistRapHi = 2.5;
    const Int_t    mHistCenBins = 10;
    const Double_t mHistCenLow = 0;
    const Double_t mHistCenHi = 200;
    const Int_t    mHistPtBins = 600;
    const Double_t mHistPtLow = 0;
    const Double_t mHistPtHi = 6;
    const Int_t    mHistMassBins = 300;
    const Double_t mHistMassLow = 2;
    const Double_t mHistMassHi = 5;
    const Int_t    mHistZdcBins = 500;
    const Double_t mHistZdcLow = 0;
    const Double_t mHistZdcHi = 5.e4;
    const Int_t    mHistAsyPhiBins = 3000;
    const Double_t mHistAsyPhiLow = 0;
    const Double_t mHistAsyPhiHi = 0.6;

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

    hZDCMinusvsZDCPlus = new TH2D("hZDCMinusvsZDCPlus", "hZDCMinusvsZDCPlus; ZDCPlus; ZDCMinus", mHistZdcBins*4, mHistZdcLow, mHistZdcHi*4, mHistZdcBins*4, mHistZdcLow, mHistZdcHi*4);
    hZDCMinusvsZDCPlus_Sel = new TH2D("hZDCMinusvsZDCPlus_Sel", "hZDCMinusvsZDCPlus_Sel; ZDCPlus; ZDCMinus", mHistZdcBins*4, mHistZdcLow, mHistZdcHi*4, mHistZdcBins*4, mHistZdcLow, mHistZdcHi*4);
    hZDCMinusvsZDCPlus_Only2MuTrk = new TH2D("hZDCMinusvsZDCPlus_Only2MuTrk", "hZDCMinusvsZDCPlus_Only2MuTrk; ZDCPlus; ZDCMinus", mHistZdcBins*4, mHistZdcLow, mHistZdcHi*4, mHistZdcBins*4, mHistZdcLow, mHistZdcHi*4);

    hZDCMinusvsZDCPlus_LS = new TH2D("hZDCMinusvsZDCPlus_LS", "hZDCMinusvsZDCPlus_LS; ZDCPlus; ZDCMinus", 1000, 0, 1.e6, 1000, 0, 1.e6);
    hZDCMinusvsZDCPlus_Sel_LS = new TH2D("hZDCMinusvsZDCPlus_Sel_LS", "hZDCMinusvsZDCPlus_Sel_LS; ZDCPlus; ZDCMinus", 1000, 0, 1.e6, 1000, 0, 1.e6);
    hZDCMinusvsZDCPlus_Only2MuTrk_LS = new TH2D("hZDCMinusvsZDCPlus_Only2MuTrk_LS", "hZDCMinusvsZDCPlus_Only2MuTrk_LS; ZDCPlus; ZDCMinus", 1000, 0, 1.e6, 1000, 0, 1.e6);

    hHFMinusvsHFPlusvsCen = new TH3D("hHFMinusvsHFPlusvsCen", "hHFMinusvsHFPlusvsCen; Centrality; HFsumETPlus; HFsumETMinus", 200, 0, 200, 200, 0, 100, 200, 0, 100);
    hHFMinusvsHFPlusvsCen_Sel = new TH3D("hHFMinusvsHFPlusvsCen_Sel", "hHFMinusvsHFPlusvsCen_Sel; Centrality; HFsumETPlus; HFsumETMinus", 200, 0, 200, 200, 0, 100, 200, 0, 100);

    hNtrkHPvsNtrkofflinevsCen = new TH3D("hNtrkHPvsNtrkofflinevsCen", "hNtrkHPvsNtrkofflinevsCen; Centrality; N_{trk}^{offline}; N_{trk}^{HP}", 200, 0, 200, 100, 0, 100, 100, 0, 100);
    hNtrkHPvsNtrkofflinevsCen_Sel = new TH3D("hNtrkHPvsNtrkofflinevsCen_Sel", "hNtrkHPvsNtrkofflinevsCen_Sel; Centrality; N_{trk}^{offline}; N_{trk}^{HP}", 200, 0, 200, 100, 0, 100, 100, 0, 100);
    hNtrkHP_2SoftMuons = new TH1D("hNtrkHP_2SoftMuons", "hNtrkHP_2SoftMuons; N_{trk}^{HP}", 100, 0, 100);

    for(Int_t idir = 0; idir < nDirs; idir++){
        hZDCvsNeuNum[idir] = new TH2D(Form("hZDC%svsNeuNum", mDir[idir].Data()), Form("; # of Neutruon; ZDC%s", mDir[idir].Data()), 5, -0.5, 4.5, mHistZdcBins*4, mHistZdcLow*4, mHistZdcHi*4);
    }

    hNeuNumMinusvsNeuNumPlus = new TH2D("hNeuNumMinusvsNeuNumPlus", "hNeuNumMinusvsNeuNumPlus; NeuNumPlus; NeuNumMinus", 3, 0, 3, 3, 0, 3);

    hMuPtvsEtavsRap          = new TH3D("hMuPtvsEtavsRap", "hMuPtvsEtavsRap; Rapidity; #eta; p_{T} (GeV)", mHistRapBins, mHistRapLow, mHistRapHi, 250, -2.5, 2.5, 200, 0, 10);
    hTrigMuPtvsEtavsRap      = new TH3D("hTrigMuPtvsEtavsRap", "hTrigMuPtvsEtavsRap; Rapidity; #eta; p_{T} (GeV)", mHistRapBins, mHistRapLow, mHistRapHi, 250, -2.5, 2.5, 200, 0, 10);
    hNegMuPtvsPosMuPtvsRap   = new TH3D("hNegMuPtvsPosMuPtvsRap", "hNegMuPtvsPosMuPtvsRap; Rapidity; #mu^{+} p_{T} (GeV); #mu^{-} p_{T} (GeV)", mHistRapBins, mHistRapLow, mHistRapHi, 200, 0, 10, 200, 0, 10);
    hNegMuEtavsPosMuEtavsRap = new TH3D("hNegMuEtavsPosMuEtavsRap", "hNegMuEtavsPosMuEtavsRap; Rapidity; #mu^{+} #eta; #mu^{-} #eta", mHistRapBins, mHistRapLow, mHistRapHi, 250, -2.5, 2.5, 250, -2.5, 2.5);
    hDeltaPtvsDeltaEtavsRap  = new TH3D("hDeltaPtvsDeltaEtavsRap", "hDeltaPtvsDeltaEtavsRap; Rapidity; #Delta#eta; #Deltap_{T} (GeV)", mHistRapBins, mHistRapLow, mHistRapHi, 250, -2.5, 2.5, 250, -2.5, 2.5);

    hDeltaEtavsDetaPhi4Cosmic = new TH2D("hDeltaEtavsDetaPhi4Cosmic", "hDeltaEtavsDetaPhi4Cosmic; #Delta#phi; #Delta#eta", 180, -PI, PI, 500, -5, 5);

    hMvsPtvsRap = new TH3D("hMvsPtvsRap", "hMvsPtvsRap; Rapidity; p_{T} (GeV); M_{#mu#mu} (GeV)", mHistRapBins, mHistRapLow, mHistRapHi, mHistPtBins, mHistPtLow, mHistPtHi, mHistMassBins, mHistMassLow, mHistMassHi);
    hMvsAsyPhivsRap = new TH3D("hMvsAsyPhivsRap", "hMvsAsyPhivsRap; Rapidity; #alpha; M_{#mu#mu} (GeV)", mHistRapBins, mHistRapLow, mHistRapHi, mHistAsyPhiBins, mHistAsyPhiLow, mHistAsyPhiHi, mHistMassBins, mHistMassLow, mHistMassHi);
    hMvsPtvsRap_WS = new TH3D("hMvsPtvsRap_WS", "hMvsPtvsRap_WS; Rapidity; p_{T} (GeV); M_{#mu#mu} (GeV)", mHistRapBins, mHistRapLow, mHistRapHi, mHistPtBins, mHistPtLow, mHistPtHi, mHistMassBins, mHistMassLow, mHistMassHi);

    for(Int_t ip=0; ip<nNeus; ip++){
        for(Int_t im=0; im<nNeus; im++){
            hMvsPtvsRap_NeuDir[ip][im]     = new TH3D(Form("hMvsPtvsRap_NeuDir%dp%dm", ip, im), "; Rapidity; p_{T} (GeV); M_{#mu#mu} (GeV)", mHistRapBins, mHistRapLow, mHistRapHi, mHistPtBins, mHistPtLow, mHistPtHi, mHistMassBins, mHistMassLow, mHistMassHi);
            hMvsAsyPhivsRap_NeuDir[ip][im] = new TH3D(Form("hMvsAsyPhivsRap_NeuDir%dp%dm", ip, im), "; Rapidity; #alpha; M_{#mu#mu} (GeV)", mHistRapBins, mHistRapLow, mHistRapHi, mHistAsyPhiBins, mHistAsyPhiLow, mHistAsyPhiHi, mHistMassBins, mHistMassLow, mHistMassHi);
        }
    }
}

void writeHistos(TString fileName)
{
    TFile* fOut = new TFile(Form("%s.root", fileName.Data()), "recreate");

    fOut->cd();

    hnEvts->Write();
    hCenvsTrig->Write();
    hPreScalevsTrig->Write();

    hVzvsVyvsVx->Write();
    hVzvsVyvsVx_Sel->Write();

    hRawCen->Write();
    hCen_Final->Write();

    hHFMinusvsHFPlusvsCen->Write();
    hHFMinusvsHFPlusvsCen_Sel->Write();

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
    }

    hNeuNumMinusvsNeuNumPlus->Write();

    hMuPtvsEtavsRap->Write();
    hTrigMuPtvsEtavsRap->Write();
    hNegMuPtvsPosMuPtvsRap->Write();
    hNegMuEtavsPosMuEtavsRap->Write();
    hDeltaPtvsDeltaEtavsRap->Write();

    hDeltaEtavsDetaPhi4Cosmic->Write();

    hMvsPtvsRap->Write();
    hMvsAsyPhivsRap->Write();
    hMvsPtvsRap_WS->Write();

    for(Int_t ip=0; ip<nNeus; ip++){
        for(Int_t im=0; im<nNeus; im++){
            hMvsPtvsRap_NeuDir[ip][im]->Write();
            hMvsAsyPhivsRap_NeuDir[ip][im]->Write();
        }
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

Int_t grabNeutronNum(Int_t dirIdx, Double_t zdc)
{
    for(Int_t ineu = 0; ineu < nNeus; ineu++){
        if(zdc > mNeuZDCLow[dirIdx][ineu] && zdc <= mNeuZDCHi[dirIdx][ineu]) return ineu;
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
