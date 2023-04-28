#!/bin/sh
echo "  "
echo " Running combine tool for all the mass points "
echo "  "

mass=mass
massf=massf

for i in 125 130 140 160 180 200 250 

do

  subrange=SR1

  echo "  "
  echo "  "
  
  eval "mass=${i}"
  

  echo "Computing limits for Ma = $mass GeV"
 
  cp hbb_mbb${mass}_SR1_mssm-13TeV_2017_SL.txt hbb_mbb${mass}_SR1_mssm-13TeV_Run2.txt

  combine -M AsymptoticLimits -t -1 hbb_mbb${mass}_SR1_mssm-13TeV_Run2.txt --setParameterRanges r=-1,1 --setParameters mask_SR=1,r=0 --cminDefaultMinimizerStrategy 0 -n Hbb_SR1 -m ${i} -v 5 

done


for i in 250 

do

  subrange=SR2

  echo "  "
  echo "  "
  
  eval "mass=${i}"
  

  echo "Computing limits for Ma = $mass GeV"
 
  cp hbb_mbb${mass}_SR2_mssm-13TeV_2017_SL.txt hbb_mbb${mass}_SR2_mssm-13TeV_Run2.txt

  combine -M AsymptoticLimits -t -1 hbb_mbb${mass}_SR2_mssm-13TeV_Run2.txt --setParameterRanges r=-1,1 --setParameters mask_SR=1,r=0 --cminDefaultMinimizerStrategy 0 -n Hbb_SR2 -m ${i} -v 5 

done


for i in 300 350 400 

do

  subrange=SR2

  echo "  "
  echo "  "
  
  eval "mass=${i}"
  

  echo "Computing limits for Ma = $mass GeV"
  
  if [ "$i" != 400 ]; then
  
  combineCards.py Run2016_FH=hbb_mbb${mass}_SR1_mssm-13TeV_2016_FH.txt Run2017_FH=hbb_mbb${mass}_SR1_mssm-13TeV_2017_FH.txt Run2018_FH=hbb_mbb${mass}_SR1_mssm-13TeV_2018_FH.txt Run2017_SL=hbb_mbb${mass}_SR2_mssm-13TeV_2017_SL.txt > hbb_mbb${mass}_SR2_mssm-13TeV_Run2.txt

  combine -M AsymptoticLimits -t -1 hbb_mbb${mass}_SR2_mssm-13TeV_Run2.txt --setParameterRanges r=-1,1 --setParameters mask_SR=1,r=0 --cminDefaultMinimizerStrategy 0 -n Hbb_SR2 -m ${i} -v 5 

  else
  
  combineCards.py Run2017_FH=hbb_mbb${mass}_SR1_mssm-13TeV_2017_FH.txt Run2018_FH=hbb_mbb${mass}_SR1_mssm-13TeV_2018_FH.txt Run2017_SL=hbb_mbb${mass}_SR2_mssm-13TeV_2017_SL.txt > hbb_mbb${mass}_SR2_mssm-13TeV_Run2.txt

  #combine -M AsymptoticLimits -t -1 hbb_mbb${mass}_SR2_mssm-13TeV_Run2.txt --setParameterRanges r=-1,1 --setParameters mask_SR=1,r=0 --cminDefaultMinimizerStrategy 0 -n Hbb_SR2 -m ${i} -v 5

  fi

done


for i in 400 

do

  subrange=SR3

  echo "  "
  echo "  "
  
  eval "mass=${i}"
  

  echo "Computing limits for Ma = $mass GeV"
  
  combineCards.py Run2016_FH=hbb_mbb${mass}_SR2_mssm-13TeV_2016_FH.txt Run2017_FH=hbb_mbb${mass}_SR2_mssm-13TeV_2017_FH.txt Run2018_FH=hbb_mbb${mass}_SR2_mssm-13TeV_2018_FH.txt Run2017_SL=hbb_mbb${mass}_SR3_mssm-13TeV_2017_SL.txt > hbb_mbb${mass}_SR3_mssm-13TeV_Run2.txt

  combine -M AsymptoticLimits -t -1 hbb_mbb${mass}_SR3_mssm-13TeV_Run2.txt --setParameterRanges r=-1,1 --setParameters mask_SR=1,r=0 --cminDefaultMinimizerStrategy 0 -n Hbb_SR3 -m ${i} -v 5 

done


for i in 450 500 600 700

do

  subrange=SR3

  echo "  "
  echo "  "
  
  eval "mass=${i}"
  

  echo "Computing limits for Ma = $mass GeV"
  
  
  if [ "$i" != 450 -a "$i" != 700 ]; then
  
    combineCards.py Run2016_FH=hbb_mbb${mass}_SR2_mssm-13TeV_2016_FH.txt Run2017_FH=hbb_mbb${mass}_SR2_mssm-13TeV_2017_FH.txt Run2018_FH=hbb_mbb${mass}_SR2_mssm-13TeV_2018_FH.txt Run2017_SL=hbb_mbb${mass}_SR3_mssm-13TeV_2017_SL.txt > hbb_mbb${mass}_SR3_mssm-13TeV_Run2.txt

    combine -M AsymptoticLimits -t -1 hbb_mbb${mass}_SR3_mssm-13TeV_Run2.txt --setParameterRanges r=-1,1 --setParameters mask_SR=1,r=0 --cminDefaultMinimizerStrategy 0 -n Hbb_SR3 -m ${i} -v 5 
    
  else 
  
    combineCards.py Run2017_FH=hbb_mbb${mass}_SR2_mssm-13TeV_2017_FH.txt Run2018_FH=hbb_mbb${mass}_SR2_mssm-13TeV_2018_FH.txt Run2017_SL=hbb_mbb${mass}_SR3_mssm-13TeV_2017_SL.txt > hbb_mbb${mass}_SR3_mssm-13TeV_Run2.txt

    #combine -M AsymptoticLimits -t -1 hbb_mbb${mass}_SR3_mssm-13TeV_Run2.txt --setParameterRanges r=-1,1 --setParameters mask_SR=1,r=0 --cminDefaultMinimizerStrategy 0 -n Hbb_SR3 -m ${i} -v 5 

  fi

done

for i in 700

do

  subrange=SR4

  echo "  "
  echo "  "
  
  eval "mass=${i}"
  

  echo "Computing limits for Ma = $mass GeV"
  
  combineCards.py Run2016_FH=hbb_mbb${mass}_SR3_mssm-13TeV_2016_FH.txt Run2017_FH=hbb_mbb${mass}_SR3_mssm-13TeV_2017_FH.txt Run2018_FH=hbb_mbb${mass}_SR3_mssm-13TeV_2018_FH.txt Run2017_SL=hbb_mbb${mass}_SR3_mssm-13TeV_2017_SL.txt > hbb_mbb${mass}_SR4_mssm-13TeV_Run2.txt

  combine -M AsymptoticLimits -t -1 hbb_mbb${mass}_SR4_mssm-13TeV_Run2.txt --setParameterRanges r=-1,1 --setParameters mask_SR=1,r=0 --cminDefaultMinimizerStrategy 0 -n Hbb_SR4 -m ${i} -v 5 

done


for i in 800 900 1000 

do

  subrange=SR4

  echo "  "
  echo "  "
  
  eval "mass=${i}"
  

  echo "Computing limits for Ma = $mass GeV"
  
  if [ "$i" != 800 -a "$i" != 1000 ]; then
  
    combineCards.py Run2016_FH=hbb_mbb${mass}_SR3_mssm-13TeV_2016_FH.txt Run2017_FH=hbb_mbb${mass}_SR3_mssm-13TeV_2017_FH.txt Run2018_FH=hbb_mbb${mass}_SR3_mssm-13TeV_2018_FH.txt > hbb_mbb${mass}_SR4_mssm-13TeV_Run2.txt

    combine -M AsymptoticLimits -t -1 hbb_mbb${mass}_SR4_mssm-13TeV_Run2.txt --setParameterRanges r=-1,1 --setParameters mask_SR=1,r=0 --cminDefaultMinimizerStrategy 0 -n Hbb_SR4 -m ${i} -v 5 
    
  else
  
  combineCards.py Run2017_FH=hbb_mbb${mass}_SR3_mssm-13TeV_2017_FH.txt Run2018_FH=hbb_mbb${mass}_SR3_mssm-13TeV_2018_FH.txt > hbb_mbb${mass}_SR4_mssm-13TeV_Run2.txt

  #combine -M AsymptoticLimits -t -1 hbb_mbb${mass}_SR4_mssm-13TeV_Run2.txt --setParameterRanges r=-1,1 --setParameters mask_SR=1,r=0 --cminDefaultMinimizerStrategy 0 -n Hbb_SR4 -m ${i} -v 5 

  fi

done


for i in 1000 1200 1400 1600 1800 

do

  subrange=SR5

  echo "  "
  echo "  "
  
  eval "mass=${i}"
  

  echo "Computing limits for Ma = $mass GeV"
  
  combineCards.py Run2017_FH=hbb_mbb${mass}_SR4_mssm-13TeV_2017_FH.txt Run2018_FH=hbb_mbb${mass}_SR4_mssm-13TeV_2018_FH.txt > hbb_mbb${mass}_SR5_mssm-13TeV_Run2.txt
  
  combine -M AsymptoticLimits -t -1 hbb_mbb${mass}_SR5_mssm-13TeV_Run2.txt --setParameterRanges r=-1,1 --setParameters mask_SR=1,r=0 --cminDefaultMinimizerStrategy 0 -n Hbb_SR5 -m ${i} -v 5 

done

