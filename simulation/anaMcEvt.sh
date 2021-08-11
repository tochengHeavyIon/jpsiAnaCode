#!/bin/bash

root -l -b -q 'anaMcEvt.C+("LowMassGammaGamma")'
root -l -b -q 'anaMcEvt.C+("CohJpsi")'
root -l -b -q 'anaMcEvt.C+("CohJpsi_0n0n")'
root -l -b -q 'anaMcEvt.C+("CohJpsi_0nXn")'
root -l -b -q 'anaMcEvt.C+("CohJpsi_XnXn")'
root -l -b -q 'anaMcEvt.C+("InCohJpsi")'
root -l -b -q 'anaMcEvt.C+("CohPsi2S")'
root -l -b -q 'anaMcEvt.C+("CohPsi2SFeeddown")'
root -l -b -q 'anaMcEvt.C+("InCohPsi2S")'
