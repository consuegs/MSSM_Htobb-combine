import os
import csv

SR1 = ["SR1"]
SR2 = ["SR2"]
SR3 = ["SR3"]
SR4 = ["SR4"]

masses_SR1 = [300, 350, 400]
masses_SR2 = [400, 450, 500, 600]
masses_SR3 = [600, 700, 800, 900, 1000]
masses_SR4 = [1000, 1200, 1400, 1600, 1800]

template_file = "/nfs/dust/cms/user/leyvaped/Analyses/MSSM/FullRun2/Combine/April_2023/CMSSW_11_3_4/src/Analysis/Combine/Run2018/Inputs_Unc/datacard_Analysis_template_VR.txt"
L1Prefiring_csv_file = "2018_Full_Hadronic_L1Prefiring.csv"
HEMcorrect_csv_file  = "2018_Full_Hadronic_HEMcorrection.csv"
OnlbtagSF_csv_file   = "2018_Full_Hadronic_onlinebtag_syst.csv"
OffbtagSF_csv_file   = "2018_Full_Hadronic_offlinebtag_syst.csv"
JKTESF_csv_file      = "2018_Full_Hadronic_jkte_syst.csv"

with open(HEMcorrect_csv_file, "r") as f:
    reader = csv.DictReader(f)
    HEMcorrect_systematics_dict = {row["Higgs_mass_GeV"]: (row["systematics_symmetric"]) for row in reader}

with open(L1Prefiring_csv_file, "r") as f:
    reader = csv.DictReader(f)
    L1Prefiring_systematics_dict = {row["Higgs_mass_GeV"]: (row["systematics_down"], row["systematics_up"]) for row in reader}

with open(OnlbtagSF_csv_file, "r") as f:
    reader = csv.DictReader(f)
    OnlbtagSF_systematics_dict = {row["Higgs_mass_GeV"]: (row["systematics_down"], row["systematics_up"]) for row in reader}
    
with open(OffbtagSF_csv_file, "r") as f:
    reader = csv.DictReader(f)
    OffbtagSF_systematics_dict = {row["Higgs_mass_GeV"]: (row["systematics_down"], row["systematics_up"]) for row in reader}

with open(JKTESF_csv_file, "r") as f:
    reader = csv.DictReader(f)
    JKTESF_systematics_dict = {row["Higgs_mass_GeV"]: (row["systematics_down"], row["systematics_up"]) for row in reader}



for subrange, masses in [("SR1", masses_SR1), ("SR2", masses_SR2), ("SR3", masses_SR3), ("SR4", masses_SR4)]:
    for mass in masses:
        os.chdir("datacards_vr_bias")
        with open(template_file, "r") as f:
            datacard = f.read().replace("MASS", str(mass)).replace("SUBRANGE", subrange)
            
        if str(mass) in HEMcorrect_systematics_dict:
            HEM_symm = HEMcorrect_systematics_dict[str(mass)]
            datacard = datacard.replace("HEM_SYST", HEM_symm)

        if str(mass) in L1Prefiring_systematics_dict:
            L1_pref_down, L1_pref_up = L1Prefiring_systematics_dict[str(mass)]
            datacard = datacard.replace("PREF_DOWN", L1_pref_down).replace("PREF_UP", L1_pref_up)

        if str(mass) in OnlbtagSF_systematics_dict:
            OnlbtagSF_down, OnlbtagSF_up = OnlbtagSF_systematics_dict[str(mass)]
            datacard = datacard.replace("ONLBTAG_DOWN", OnlbtagSF_down).replace("ONLBTAG_UP", OnlbtagSF_up)

        if str(mass) in OnlbtagSF_systematics_dict:
            OffbtagSF_down, OffbtagSF_up = OffbtagSF_systematics_dict[str(mass)]
            datacard = datacard.replace("OFFBTAG_DOWN", OffbtagSF_down).replace("OFFBTAG_UP", OffbtagSF_up)

        if str(mass) in JKTESF_systematics_dict:
            JKTESF_down, JKTESF_up = JKTESF_systematics_dict[str(mass)]
            datacard = datacard.replace("JKTE_DOWN", JKTESF_down).replace("JKTE_UP", JKTESF_up)

        with open("hbb_mbb{}_{}_mssm-13TeV.txt".format(mass, subrange), "w") as f:
            f.write(datacard)
        print("Datacard {} GeV, {} created".format(mass, subrange))
        os.system("text2workspace.py hbb_mbb{}_{}_mssm-13TeV.txt".format(mass, subrange))
        os.chdir("..")

        