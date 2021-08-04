#!/bin/bash

##void anaEvt(Bool_t incHadron = kFALSE, TString hfVetoType="Default")
root -l -b -q 'anaEvt.C+(0, "Default")'
#root -l -b -q 'anaEvt.C+(0, "Tight")'
#root -l -b -q 'anaEvt.C+(0, "Loose")'
#root -l -b -q 'anaEvt.C+(0, "removeHF")'
#root -l -b -q 'anaEvt.C+(1, "Default")'

##void plotJpsiQA(Bool_t incHadron = kFALSE, TString hfVetoType="Default")
root -l -b -q 'plotJpsiQA.C+(0, "Default")'
#root -l -b -q 'plotJpsiQA.C+(0, "Tight")'
#root -l -b -q 'plotJpsiQA.C+(0, "Loose")'
#root -l -b -q 'plotJpsiQA.C+(0, "removeHF")'
#root -l -b -q 'plotJpsiQA.C+(1, "Default")'
