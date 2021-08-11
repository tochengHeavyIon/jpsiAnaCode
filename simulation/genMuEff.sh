#!/bin/bash

root -l -b -q 'genMuEff.C+("LowMassGammaGamma")'
root -l -b -q 'genMuEff.C+("CohJpsi")'
root -l -b -q 'genMuEff.C+("CohJpsi_0n0n")'
root -l -b -q 'genMuEff.C+("CohJpsi_0nXn")'
root -l -b -q 'genMuEff.C+("CohJpsi_XnXn")'
root -l -b -q 'genMuEff.C+("InCohJpsi")'
root -l -b -q 'genMuEff.C+("CohPsi2SFeeddown")'
root -l -b -q 'genMuEff.C+("CohPsi2S")'
root -l -b -q 'genMuEff.C+("InCohPsi2S")'
