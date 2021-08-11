#!/bin/bash

date

nAnaEvts=10000000

root -l -b -q 'anaSTARlight.C+("CohJpsi", '${nAnaEvts}')'
root -l -b -q 'anaSTARlight.C+("CohJpsi_0n0n", '${nAnaEvts}')'
root -l -b -q 'anaSTARlight.C+("CohJpsi_0nXn", '${nAnaEvts}')'
root -l -b -q 'anaSTARlight.C+("CohJpsi_XnXn", '${nAnaEvts}')'
root -l -b -q 'anaSTARlight.C+("InCohJpsi", '${nAnaEvts}')'
root -l -b -q 'anaSTARlight.C+("CohPsi2S", '${nAnaEvts}')'
root -l -b -q 'anaSTARlight.C+("InCohPsi2S", '${nAnaEvts}')'
