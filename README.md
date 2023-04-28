#  MSSM_Htobb-combine

Run2 Legacy MSSM Htobb analysis tools for statistical combination

## Instructions for Run 2 combination:

## 1. Setting up the area

```bash

export CMSSW_GIT_REFERENCE=/nfs/dust/cms/user/${USER}/.cmsgit-cache
source /cvmfs/cms.cern.ch/cmsset_default.sh

export SCRAM_ARCH=slc7_amd64_gcc700

cmsrel CMSSW_10_2_13
cd CMSSW_10_2_13/src
cmsenv

git clone https://github.com/desy-cms/analysis-models.git Analysis/Models
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit

cd $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit
git fetch origin
git checkout v8.2.0

cd ../../

scram b -j 16 USER_CXXFLAGS="-Wno-misleading-indentation"

cd ../../

```


## 2. Clone 2016 analysis repository to get corresponding datacards and workspaces


```bash
git clone https://gitlab.cern.ch/cms-hcg/cadi/hig-16-018.git
```


## 3. Calculate limits

`CopyAll.sh` - copy and organize the datacards including year subindex (year = 2016, 2017, and 2018) to facilitate combination

`runLimits.sh` - as the name suggests macro to run the limit computation

`limits.txt` - filelist for limit computation files

`PlotLimits.C` - plotting macro 

