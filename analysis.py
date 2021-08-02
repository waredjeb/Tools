from utils import *
import uproot
import matplotlib.pyplot as plt
import ROOT
import numpy as np
import pandas as pd
import mplhep as hep
plt.style.use([hep.style.ROOT, hep.style.firamath])

fileName = "./data/output_SinglePi.root"
file = uproot.open(fileName)
keys = file.keys()

fileRoot = file[keys[1]]

titleSave = "_"+fileName[22:-5]

pathPlot = './plots/NoSiblings/SinglePions/'
mkdir_p(pathPlot)

sim_rawEnergy = fileRoot['sim_rawEnergy'].array(library='np')
sim_regressedEnergy = fileRoot['sim_regressedEnergy'].array(library='np')
lc_prEnergy = fileRoot['lc_prEnergy'].array(library='np')
lc_allEnergy = fileRoot['lc_allEnergy'].array(library='np')
lc_noPR = fileRoot['lc_noPR'].array(library='np')
lc_noUsed = fileRoot['lc_noUsed'].array(library='np')
EM_raw = fileRoot['EM_raw'].array(library='np')
EM_regressed = fileRoot['EM_regressed'].array(library='np')
HAD_raw = fileRoot['HAD_raw'].array(library='np')
HAD_regressed = fileRoot['HAD_raw'].array(library='np')
CaloParticleEnergy = fileRoot['CaloParticleEnergy'].array(library='np')
lcNoPR_drTrackster = fileRoot['lcNotPR_dRTrackster'].array(library = 'np')
lcNotUsed_dRTrackster = fileRoot['lcNotUsed_dRTrackster'].array(library = 'np')

#Linking
lcNoPR_spectrum = fileRoot['lcNoPR_spectrum'].array(library = 'np')
mask_input = fileRoot['full_mask_input'].array(library = 'np')
mask_output = fileRoot['full_mask_output'].array(library = 'np')
lc_z = fileRoot['lc_z'].array(library = 'np')
lc_E = fileRoot['lc_E'].array(library = 'np')
layers_id = fileRoot['layers_id'].array(library = 'np')
layers_thick = fileRoot['layers_thick'].array(library = 'np')

totalTracksterEnergy_raw = list()
totalTracksterEnergy_regressed = list()

ratioEM_sim = list()
ratioTrkEm_sim = list()
ratioTrk_sim = list()
ratioHad_sim = list()

ratioEM_regressed = list()
ratioTrkEm_regressed = list()
ratioTrk_regressed = list()
ratioHad_regressed = list()

ratioEM_lcPR = list()
ratioTrkEm_lcPR = list()
ratioTrk_lcPR = list()
ratioHad_lcPR = list()

ratioTot_sim_raw = list()
ratioTot_sim_regressed = list()

ratioTot_lcPR = list()
ratioTot_sim_regressed = list()

ratioLCPR_Sim = list()
ratioLCPR_LCNOPR_SIM = list()
ratioLCNOPR_SimRaw = list()

ratioLCNOPR_CP = list()

ratioLinking_sim = list()
ratioLinking_regressed = list()
ratioLinking_lcPR = list()

for i in range(len(sim_rawEnergy)):
    if(sim_regressedEnergy[i][0] > 0 and lc_prEnergy[i][0] > 0):        
        tot_tmp_raw = EM_raw[i][0] +  +  + HAD_raw[i][0]
        tot_tmp_regressed = EM_regressed[i][0] + HAD_regressed[i][0]
        totalTracksterEnergy_raw.append(EM_raw[i][0] +  +  + HAD_raw[i][0])
        totalTracksterEnergy_regressed.append(EM_regressed[i][0] + HAD_regressed[i][0])
        ratioEM_sim.append(EM_raw[i][0] / sim_regressedEnergy[i][0])
        ratioHad_sim.append(HAD_raw[i][0] / sim_regressedEnergy[i][0])
        ratioTot_sim_raw.append(tot_tmp_raw/sim_regressedEnergy[i][0])
        ratioEM_regressed.append(EM_regressed[i][0] / sim_regressedEnergy[i][0])
        ratioHad_regressed.append(HAD_regressed[i][0] / sim_regressedEnergy[i][0])
        ratioTot_sim_regressed.append(tot_tmp_regressed/sim_regressedEnergy[i][0])
        ratioEM_lcPR.append(EM_raw[i][0] / lc_prEnergy[i][0])
        ratioHad_lcPR.append(HAD_raw[i][0] / lc_prEnergy[i][0])
        ratioTot_lcPR.append(tot_tmp_raw / lc_prEnergy[i][0])    
        ratioLCPR_Sim.append(lc_prEnergy[i][0] / sim_regressedEnergy[i][0] )
        ratioLCPR_LCNOPR_SIM.append((lc_prEnergy[i][0] + lc_noPR[i][0]) / sim_regressedEnergy[i][0])
        ratioLCNOPR_SimRaw.append((lc_noPR[i][0]) / sim_rawEnergy[i][0])    
        ratioLCNOPR_CP.append(lc_noPR[i][0] / CaloParticleEnergy[i][0])

zList= list()
for i in layers_id:
    for j in i:
        zList.append(j)
zs = np.unique(zList)[1:]

print(zs)
mask_input_flat = list()
mask_output_flat = list()
lc_id_flat = list()
lc_E_flat = list()
thickness_flat = list()
for ev in range(len(mask_input)):
    for i in range(len(mask_input[ev])):
        mask_input_flat.append(mask_input[ev][i])
        mask_output_flat.append(mask_output[ev][i])
        lc_id_flat.append(layers_id[ev][i])
        lc_E_flat.append(lc_E[ev][i])
        thickness_flat.append(layers_thick[ev][i])

fullLCEnergy = list()
noPRLCEnergy = list()
for z in zs: #all but -999
    fullLCEnergy_z = list()
    noPRLCEnergy_z = list()
    for i in range(len(mask_input_flat)):
        if(lc_id_flat[i] == z): 
            fullLCEnergy_z.append(lc_E_flat[i])#tutti i lcs divisi per z
            if(mask_input_flat[i] == 0):
                noPRLCEnergy_z.append(lc_E_flat[i])
                
    
    fullLCEnergy.append(fullLCEnergy_z)
    noPRLCEnergy.append(noPRLCEnergy_z)

fullLCEnergyHD = list()
noPRLCEnergyHD = list()
for z in zs: #all but -999
    fullLCEnergy_z = list()
    noPRLCEnergy_z = list()
    for i in range(len(mask_input_flat)):
        if(lc_id_flat[i] == z and thickness_flat[i] <= 120): 
            fullLCEnergy_z.append(lc_E_flat[i])#tutti i lcs divisi per z
            if(mask_input_flat[i] == 0):
                noPRLCEnergy_z.append(lc_E_flat[i])
                
    
    fullLCEnergyHD.append(fullLCEnergy_z)
    noPRLCEnergyHD.append(noPRLCEnergy_z)

fullLCEnergyLD = list()
noPRLCEnergyLD = list()
for z in zs: #all but -999
    fullLCEnergy_z = list()
    noPRLCEnergy_z = list()
    for i in range(len(mask_input_flat)):
        if(lc_id_flat[i] == z and thickness_flat[i] > 120): 
            fullLCEnergy_z.append(lc_E_flat[i])#tutti i lcs divisi per z
            if(mask_input_flat[i] == 0):
                noPRLCEnergy_z.append(lc_E_flat[i])
                
    
    fullLCEnergyLD.append(fullLCEnergy_z)
    noPRLCEnergyLD.append(noPRLCEnergy_z)
    
profX = ROOT.TProfile("Profile", "Profile", 50, 0, 50, 0, 1)
for i in range(len(fullLCEnergy)):
    sumZ = 0
    sumZNoPr = 0
    for e in fullLCEnergy[i]:
        sumZ += e
    for e in noPRLCEnergy[i]:
        sumZNoPr += e
    profX.Fill(zs[i], sumZNoPr/ sumZ)
    
profXLD = ROOT.TProfile("ProfileLD", "ProfileLD", 50, 0, 50, 0, 1)
for i in range(len(fullLCEnergyLD)):
    sumZ = 0
    sumZNoPr = 0
    for e in fullLCEnergyLD[i]:
        sumZ += e
    for e in noPRLCEnergyLD[i]:
        sumZNoPr += e
#     print(sumZ)
    if(sumZ == 0):
        profXLD.Fill(zs[i], -99)
    else:
        profXLD.Fill(zs[i], sumZNoPr/ sumZ)
    
profXHD = ROOT.TProfile("ProfileHD", "ProfileHD", 50, 0, 50, 0, 1)
for i in range(len(fullLCEnergyHD)):
    sumZ = 0
    sumZNoPr = 0
    for e in fullLCEnergyHD[i]:
        sumZ += e
    for e in noPRLCEnergyHD[i]:
        sumZNoPr += e
    if(sumZ == 0):
        profXHD.Fill(zs[i], -99)
    else:
        profXHD.Fill(zs[i], sumZNoPr/ sumZ)

c = ROOT.TCanvas("c", "c", 1000,1000)
# profX = histoAll.ProfileX()
profX.SetMarkerColor(4)
profX.SetMarkerSize(1.5)
profX.SetMarkerStyle(20)
profX.SetStats(0000)
profX.SetTitle("LC NHits < 3 / All LCs ")
profX.GetXaxis().SetTitle("HGCal Layer")
profX.GetYaxis().SetTitle("Ratio")
profX.Draw("P")
c.SaveAs(pathPlot + "RatioLCPRAlongLayerID.png", "png")

c2 = ROOT.TCanvas("c2", "c2", 1000,1000)
leg = ROOT.TLegend(0.3,0.7,0.5,0.9)
# profX = histoAll.ProfileX()
profX.GetYaxis().SetRangeUser(0,1)
profX.SetMarkerColor(4)
profX.SetMarkerSize(1.2)
profX.SetMarkerStyle(20)
profX.SetStats(0000)
profX.SetTitle("LC NHits < 3 Energy  / All LCs Energy")
profX.GetXaxis().SetTitle("HGCal Layer")
profX.GetYaxis().SetTitle("Ratio")
profXLD.SetMarkerColor(ROOT.kRed)
profXLD.SetMarkerSize(1.2)
profXLD.SetMarkerStyle(20)
profXLD.SetStats(0000)
profXHD.SetMarkerColor(ROOT.kGreen)
profXHD.SetMarkerSize(1.2)
profXHD.SetMarkerStyle(20)
profXHD.SetStats(0000)
leg.AddEntry(profX, "All")
leg.AddEntry(profXLD, "LD >= 120")
leg.AddEntry(profXHD, "HD < 120")
profX.Draw("P")
profXLD.Draw("SAME P")
profXHD.Draw("SAME P")
leg.Draw("SAME")
c2.SaveAs(pathPlot + "RatioLCPRAlongLayerIDMulti.png", "png")
# c2.Draw()

##############################################################
# LC not in pattern recognition energy spectrumm
##############################################################
lcNoPR_energies = list()
for i in lcNoPR_spectrum:
    for j in i:
        lcNoPR_energies.append(j)

plt.figure(figsize = (13,13))
plt.hist(lcNoPR_energies, bins = 100, range = (0,4))
plt.xlabel('Energy [GeV]')
plt.title("LC Nhits < 3 Total Energy")
plt.savefig(pathPlot + "LCNOPR_energies" + ".png")       

##############################################################
# RECO-TRACKSTER RAW-ENERGY WRT SIMTRACKSTER REGRESSED ENERGY
##############################################################
density = False
bins = 30
plt.figure(figsize = (10,10))
plt.hist(ratioTot_sim_raw, bins = bins, stacked = False,density = density, alpha = 1, range = (0.01, 1),label = "All", histtype='step', lw = 2)
plt.hist(ratioEM_sim, bins = bins, stacked = False,density = density, alpha = 0.4, range = (0.01, 1),label = "EM", histtype='step', lw = 2)
plt.hist(ratioHad_sim, bins = bins, stacked = False,density = density, alpha = 0.4, range = (0.01, 1),label = "HAD", histtype='step', lw = 2)
plt.legend(loc = 'upper left')
plt.xlabel("Ratio")
plt.title("Ratio reco-Trackster SimTrackster Regressed Energy", fontsize = 17)
plt.savefig(pathPlot + "RatioRecoTrackstervsSimTrackster" + ".png")
# plt.show()

##############################################################
# RECO-TRACKSTER RAW-ENERGY WRT SIMTRACKSTER REGRESSED ENERGY
##############################################################

plt.figure(figsize = (10,10))
plt.hist(ratioTot_sim_regressed, bins = bins, stacked = False,density = density, alpha = 1, range = (0.01, 1.5),label = "All", histtype='step', lw = 2)
plt.hist(ratioEM_regressed, bins = bins, stacked = False,density = density, alpha = 0.4, range = (0.01, 1.5),label = "EM", histtype='step', lw = 2)
plt.hist(ratioHad_regressed, bins = bins, stacked = False,density = density, alpha = 0.4, range = (0.01, 1.5),label = "HAD", histtype='step', lw = 2)
plt.legend(loc = 'upper left')
plt.xlabel("Ratio")
plt.title("Ratio reco-Trackster Regressed Energy vs SimTrackster Regressed Energy", fontsize = 17)
plt.savefig(pathPlot + "RatioRecoTracksterRegressedvsSimTrackster" + ".png")
# plt.show()

bins = 50
density = False
plt.figure(figsize = (10,10))
plt.hist(ratioTot_lcPR, bins = bins, stacked = False, density = density, histtype = 'step', range = (0.01,1), lw = 2, label = 'All')
plt.hist(ratioEM_lcPR, bins = bins, stacked = False, density = density, alpha = 0.5, range = (0.01, 1), label = "EM")
plt.hist(ratioHad_lcPR, bins = bins, stacked = False,density = density, alpha = 0.5, range = (0.01, 1),label = "HAD")
plt.legend(loc = 'upper left')
plt.xlabel("Ratio")
plt.ylabel("Entries")
plt.title("Ratio reco-Trackster Raw Energy vs LC NHits > 2 Total Energy", fontsize = 18)
plt.savefig(pathPlot + "RatioRecoTrackstervsLCPR" + ".png")
# plt.show()

##############################################################
# LCs NHits > 2 but not Clusterized
##############################################################

lcNoUsed_flat = []
for i in lc_noUsed:
    lcNoUsed_flat.append(i[0])    
plt.figure(figsize = (10,10))
plt.hist(lcNoUsed_flat, bins = 50, range = (0.01, 200), histtype='step')
plt.title("LCs not clusterized total energy")
plt.xlabel("Energy [GeV]")
plt.ylabel("Entries")
# plt.show()

##############################################################
#Ratio LCs NHits > 2 wrt SimmTrackster
##############################################################
plt.figure(figsize = (10,10))
plt.hist(ratioLCPR_Sim, bins = 30, histtype='step')
plt.title("LC NHits > 2 Energy / SimTrackster Energy")
plt.xlabel("Ratio Energy")
plt.savefig(pathPlot + "RatioLCPRvsSimTrk" + ".png")
# plt.show()