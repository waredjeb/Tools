from utils import *
import uproot
import matplotlib.pyplot as plt
import ROOT
import numpy as np
import pandas as pd
import mplhep as hep

plt.style.use([hep.style.ROOT, hep.style.firamath])

fileName = "./data/output_CloseByPion_NoSiblings_20-25deg.root"
file = uproot.open(fileName)
keys = file.keys()

fileRoot = file[keys[1]]

titleSave = "_" + fileName[22:-5]

pathPlot = "./plots/NoSiblings/CloseByPion_NoSiblings_20-25deg/"
mkdir_p(pathPlot)

sim_rawEnergy = fileRoot["sim_rawEnergy"].array(library="np")
sim_regressedEnergy = fileRoot["sim_regressedEnergy"].array(library="np")
lc_prEnergy = fileRoot["lc_prEnergy"].array(library="np")
lc_allEnergy = fileRoot["lc_allEnergy"].array(library="np")
lc_noPR = fileRoot["lc_noPR"].array(library="np")
lc_noUsed = fileRoot["lc_noUsed"].array(library="np")
EM_raw = fileRoot["EM_raw"].array(library="np")
EM_regressed = fileRoot["EM_regressed"].array(library="np")
HAD_raw = fileRoot["HAD_raw"].array(library="np")
HAD_regressed = fileRoot["HAD_raw"].array(library="np")
CaloParticleEnergy = fileRoot["CaloParticleEnergy"].array(library="np")
lcNoPR_drTrackster = fileRoot["lcNotPR_dRTrackster"].array(library="np")
lcNotUsed_dRTrackster = fileRoot["lcNotUsed_dRTrackster"].array(library="np")

# Linking
lcNoPR_spectrum = fileRoot["lcNoPR_spectrum"].array(library="np")
mask_input = fileRoot["full_mask_input"].array(library="np")
mask_output = fileRoot["full_mask_output"].array(library="np")
lc_z = fileRoot["lc_z"].array(library="np")
lc_E = fileRoot["lc_E"].array(library="np")
layers_id = fileRoot["layers_id"].array(library="np")
layers_thick = fileRoot["layers_thick"].array(library="np")
lc_nhits = fileRoot["LC_nhits"].array(library="np")
n_tracksters = fileRoot["n_tracksters"].array(library="np")
n_trackstersEM = fileRoot["n_trackstersEM"].array(library="np")
n_trackstersHAD = fileRoot["n_trackstersHAD"].array(library="np")

n_lc = fileRoot["n_LCtrackster"].array(library="np")
n_lcEM = fileRoot["n_LCtracksterEM"].array(library="np")
n_lcHAD = fileRoot["n_LCtracksterHAD"].array(library="np")
lc_id_trackster = fileRoot["lc_id_trackster"].array(library="np")
lc_id_tracksterEM = fileRoot["lc_id_tracksterEM"].array(library="np")
lc_id_tracksterHAD = fileRoot["lc_id_tracksterHAD"].array(library="np")


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

N_trackster = list()
N_tracksterEM = list()
N_tracksterHAD = list()

N_lc = list()
N_lcEM = list()
N_lcHAD = list()

for i in range(len(sim_rawEnergy)):
    if sim_regressedEnergy[i][0] > 0 and lc_prEnergy[i][0] > 0:
        tot_tmp_raw = EM_raw[i][0] + HAD_raw[i][0]
        tot_tmp_regressed = EM_regressed[i][0] + HAD_regressed[i][0]
        totalTracksterEnergy_raw.append(EM_raw[i][0] + HAD_raw[i][0])
        totalTracksterEnergy_regressed.append(EM_regressed[i][0] + HAD_regressed[i][0])
        ratioEM_sim.append(EM_raw[i][0] / sim_regressedEnergy[i][0])
        ratioHad_sim.append(HAD_raw[i][0] / sim_regressedEnergy[i][0])
        ratioTot_sim_raw.append(tot_tmp_raw / sim_regressedEnergy[i][0])
        ratioEM_regressed.append(EM_regressed[i][0] / sim_regressedEnergy[i][0])
        ratioHad_regressed.append(HAD_regressed[i][0] / sim_regressedEnergy[i][0])
        ratioTot_sim_regressed.append(tot_tmp_regressed / sim_regressedEnergy[i][0])
        ratioEM_lcPR.append(EM_raw[i][0] / lc_prEnergy[i][0])
        ratioHad_lcPR.append(HAD_raw[i][0] / lc_prEnergy[i][0])
        ratioTot_lcPR.append(tot_tmp_raw / lc_prEnergy[i][0])
        ratioLCPR_Sim.append(lc_prEnergy[i][0] / sim_regressedEnergy[i][0])
        ratioLCPR_LCNOPR_SIM.append(
            (lc_prEnergy[i][0] + lc_noPR[i][0]) / sim_regressedEnergy[i][0]
        )
        ratioLCNOPR_SimRaw.append((lc_noPR[i][0]) / sim_rawEnergy[i][0])
        ratioLCNOPR_CP.append(lc_noPR[i][0] / CaloParticleEnergy[i][0])

zList = list()
for i in layers_id:
    for j in i:
        zList.append(j)
zs = np.unique(zList)[1:]

hgcal_nhits_per_ev = list()
print(len(hgcal_nhits_per_ev), len(zs))

for ev in range(len(lc_nhits)):
    list_per_ev = [0 for i in range(len(zs))]
    for z in range(len(zs)):
        for j in range(len(lc_nhits[ev])):
            if z == layers_id[ev][j]:
                list_per_ev[z] += lc_nhits[ev][j]
    hgcal_nhits_per_ev.append(list_per_ev)


cut_nhits = 10
fullLCEnergy = list()
noPRLCEnergy = list()
for ev in range(len(mask_input)):
    fullLCEnergy_z = [0 for i in range(len(zs))]
    noPRLCEnergy_z = [0 for i in range(len(zs))]
    for z in range(len(zs)):
        for i in range(len(mask_input[ev])):
            if layers_id[ev][i] == zs[z] and hgcal_nhits_per_ev[ev][z] > cut_nhits:
                fullLCEnergy_z[z] += lc_E[ev][i]
                if mask_input[ev][i] == 0:
                    noPRLCEnergy_z[z] += lc_E[ev][i]
    fullLCEnergy.append(fullLCEnergy_z)
    noPRLCEnergy.append(noPRLCEnergy_z)

fullLCEnergyLD = list()
noPRLCEnergyLD = list()
for ev in range(len(mask_input)):
    fullLCEnergy_z = [0 for i in range(len(zs))]
    noPRLCEnergy_z = [0 for i in range(len(zs))]
    for z in range(len(zs)):
        for i in range(len(mask_input[ev])):
            if (
                layers_id[ev][i] == zs[z]
                and hgcal_nhits_per_ev[ev][z] > cut_nhits
                and layers_thick[ev][i] > 120
            ):
                fullLCEnergy_z[z] += lc_E[ev][i]
                if mask_input[ev][i] == 0:
                    noPRLCEnergy_z[z] += lc_E[ev][i]
    fullLCEnergyLD.append(fullLCEnergy_z)
    noPRLCEnergyLD.append(noPRLCEnergy_z)

fullLCEnergyHD = list()
noPRLCEnergyHD = list()
for ev in range(len(mask_input)):
    fullLCEnergy_z = [0 for i in range(len(zs))]
    noPRLCEnergy_z = [0 for i in range(len(zs))]
    for z in range(len(zs)):
        for i in range(len(mask_input[ev])):
            if (
                layers_id[ev][i] == zs[z]
                and hgcal_nhits_per_ev[ev][z] > cut_nhits
                and layers_thick[ev][i] <= 120
            ):
                fullLCEnergy_z[z] += lc_E[ev][i]
                if mask_input[ev][i] == 0:
                    noPRLCEnergy_z[z] += lc_E[ev][i]
    fullLCEnergyHD.append(fullLCEnergy_z)
    noPRLCEnergyHD.append(noPRLCEnergy_z)

histo = ROOT.TH2F("Profile", "Profile", 50, 0, 50, 50, 0, 1)

for ev in range(len(fullLCEnergy)):
    for z in range(len(zs)):
        if fullLCEnergy[ev][z] > 0:
            histo.Fill(zs[z], noPRLCEnergy[ev][z] / fullLCEnergy[ev][z])

histoLD = ROOT.TH2F("ProfileLD", "ProfileLD", 50, 0, 50, 50, 0, 1)

for ev in range(len(fullLCEnergyLD)):
    for z in range(len(zs)):
        if fullLCEnergyLD[ev][z] > 0:
            histoLD.Fill(zs[z], noPRLCEnergyLD[ev][z] / fullLCEnergyLD[ev][z])

histoHD = ROOT.TH2F("ProfileHD", "ProfileHD", 50, 0, 50, 50, 0, 1)

for ev in range(len(fullLCEnergyHD)):
    for z in range(len(zs)):
        if fullLCEnergyHD[ev][z] > 0:
            histoHD.Fill(zs[z], noPRLCEnergyHD[ev][z] / fullLCEnergyHD[ev][z])

for ev in range(len(n_tracksters)):
    for n in n_tracksters[ev]:
        N_trackster.append(n)
    for n in n_trackstersEM[ev]:
        N_tracksterEM.append(n)
    for n in n_trackstersHAD[ev]:
        N_tracksterHAD.append(n)

for ev in range(len(n_lc)):
    for n in n_lc[ev]:
        N_lc.append(n)
    for n in n_lcEM[ev]:
        N_lcEM.append(n)
    for n in n_lcHAD[ev]:
        N_lcHAD.append(n)

trackstersLenght = list()
trackstersSkippedLayers = list()
for ev in lc_id_trackster:
    for t in ev:
        tnp = np.array(t)
        tnp_unique = np.unique(tnp)
        tnpsorted = sorted(tnp_unique)
#         print(tnpsorted)
#         print(tnpsorted[-1] - tnpsorted[0])
        layerLenght = tnpsorted[-1] - tnpsorted[0]
        N_layers_trackster = len(tnpsorted)-1
        skipped_layers = abs(N_layers_trackster - layerLenght)
        if(skipped_layers > 2):
            print(tnpsorted, len(tnpsorted))
        trackstersLenght.append(N_layers_trackster)
        trackstersSkippedLayers.append(skipped_layers)

trackstersLenghtEM = list()
trackstersSkippedLayersEM = list()
for ev in lc_id_tracksterEM:
    for t in ev:
        tnp = np.array(t)
        tnp_unique = np.unique(tnp)
        tnpsorted = sorted(tnp_unique)
#         print(tnpsorted)
#         print(tnpsorted[-1] - tnpsorted[0])
        layerLenght = tnpsorted[-1] - tnpsorted[0]
        N_layers_trackster = len(tnpsorted)-1
        skipped_layers = abs(N_layers_trackster - layerLenght)
        if(skipped_layers > 2):
            print(tnpsorted, len(tnpsorted))
        trackstersLenghtEM.append(N_layers_trackster)
        trackstersSkippedLayersEM.append(skipped_layers)

trackstersLenghtHAD = list()
trackstersSkippedLayersHAD = list()
for ev in lc_id_tracksterHAD:
    for t in ev:
        tnp = np.array(t)
        tnp_unique = np.unique(tnp)
        tnpsorted = sorted(tnp_unique)
#         print(tnpsorted)
#         print(tnpsorted[-1] - tnpsorted[0])
        layerLenght = tnpsorted[-1] - tnpsorted[0]
        N_layers_trackster = len(tnpsorted)-1
        skipped_layers = abs(N_layers_trackster - layerLenght)
        if(skipped_layers > 2):
            print(tnpsorted, len(tnpsorted))
        trackstersLenghtHAD.append(N_layers_trackster)
        trackstersSkippedLayersHAD.append(skipped_layers)

c = ROOT.TCanvas("c", "c", 1000, 1000)
profX = histo.ProfileX()
profX.GetYaxis().SetRangeUser(0, 1)
profX.SetMarkerColor(4)
profX.SetMarkerSize(1.5)
profX.SetMarkerStyle(20)
profX.SetStats(0000)
profX.SetTitle("LC NHits < 3 / All LCs ")
profX.GetXaxis().SetTitle("HGCal Layer")
profX.GetYaxis().SetTitle("Ratio")
profX.Draw("P")
c.Draw()
c.SaveAs(pathPlot + "RatioLCPRAlongLayerIDAll.png", "png")
# c.Show()

c2 = ROOT.TCanvas("c", "c", 1000, 1000)
leg = ROOT.TLegend(0.3, 0.7, 0.5, 0.9)
profXLD = histoLD.ProfileX()
profXHD = histoHD.ProfileX()
profX.GetYaxis().SetRangeUser(0, 1)
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
c2.Draw()
c2.SaveAs(pathPlot + "RatioLCPRAlongLayerIDMulti.png", "png")


##############################################################
# LC not in pattern recognition energy spectrumm
##############################################################
lcNoPR_energies = list()
for i in lcNoPR_spectrum:
    for j in i:
        lcNoPR_energies.append(j)

plt.figure(figsize=(13, 13))
plt.hist(lcNoPR_energies, bins=100, range=(0, 4))
plt.xlabel("Energy [GeV]")
plt.title("LC Nhits < 3 Total Energy")
plt.savefig(pathPlot + "LCNOPR_energies" + ".png")

##############################################################
# RECO-TRACKSTER RAW-ENERGY WRT SIMTRACKSTER REGRESSED ENERGY
##############################################################
density = False
bins = 30
plt.figure(figsize=(10, 10))
plt.hist(
    ratioTot_sim_raw,
    bins=bins,
    stacked=False,
    density=density,
    alpha=1,
    range=(0.01, 1),
    label="All",
    histtype="step",
    lw=2,
)
plt.hist(
    ratioEM_sim,
    bins=bins,
    stacked=False,
    density=density,
    alpha=0.4,
    range=(0.01, 1),
    label="EM",
    histtype="step",
    lw=2,
)
plt.hist(
    ratioHad_sim,
    bins=bins,
    stacked=False,
    density=density,
    alpha=0.4,
    range=(0.01, 1),
    label="HAD",
    histtype="step",
    lw=2,
)
plt.legend(loc="upper left")
plt.xlabel("Ratio")
plt.title("Ratio reco-Trackster SimTrackster Regressed Energy", fontsize=17)
plt.savefig(pathPlot + "RatioRecoTrackstervsSimTrackster" + ".png")
# plt.show()

##############################################################
# RECO-TRACKSTER RAW-ENERGY WRT SIMTRACKSTER REGRESSED ENERGY
##############################################################

plt.figure(figsize=(10, 10))
plt.hist(
    ratioTot_sim_regressed,
    bins=bins,
    stacked=False,
    density=density,
    alpha=1,
    range=(0.01, 1.5),
    label="All",
    histtype="step",
    lw=2,
)
plt.hist(
    ratioEM_regressed,
    bins=bins,
    stacked=False,
    density=density,
    alpha=0.4,
    range=(0.01, 1.5),
    label="EM",
    histtype="step",
    lw=2,
)
plt.hist(
    ratioHad_regressed,
    bins=bins,
    stacked=False,
    density=density,
    alpha=0.4,
    range=(0.01, 1.5),
    label="HAD",
    histtype="step",
    lw=2,
)
plt.legend(loc="upper left")
plt.xlabel("Ratio")
plt.title(
    "Ratio reco-Trackster Regressed Energy vs SimTrackster Regressed Energy",
    fontsize=17,
)
plt.savefig(pathPlot + "RatioRecoTracksterRegressedvsSimTrackster" + ".png")
# plt.show()

bins = 50
density = False
plt.figure(figsize=(10, 10))
plt.hist(
    ratioTot_lcPR,
    bins=bins,
    stacked=False,
    density=density,
    histtype="step",
    range=(0.01, 1),
    lw=2,
    label="All",
)
plt.hist(
    ratioEM_lcPR,
    bins=bins,
    stacked=False,
    density=density,
    alpha=0.5,
    range=(0.01, 1),
    label="EM",
)
plt.hist(
    ratioHad_lcPR,
    bins=bins,
    stacked=False,
    density=density,
    alpha=0.5,
    range=(0.01, 1),
    label="HAD",
)
plt.legend(loc="upper left")
plt.xlabel("Ratio")
plt.ylabel("Entries")
plt.title("Ratio reco-Trackster Raw Energy vs LC NHits > 2 Total Energy", fontsize=18)
plt.savefig(pathPlot + "RatioRecoTrackstervsLCPR" + ".png")
# plt.show()

##############################################################
# LCs NHits > 2 but not Clusterized
##############################################################

lcNoUsed_flat = []
for i in lc_noUsed:
    lcNoUsed_flat.append(i[0])
plt.figure(figsize=(10, 10))
plt.hist(lcNoUsed_flat, bins=50, range=(0.01, 200), histtype="step", lw = 2)
plt.title("LCs not clusterized total energy")
plt.xlabel("Energy [GeV]")
plt.ylabel("Entries")
plt.savefig(pathPlot + "LCNoClusterizedEnergy" + ".png")
# plt.show()

##############################################################
# Ratio LCs NHits > 2 wrt SimmTrackster
##############################################################
plt.figure(figsize=(10, 10))
plt.hist(ratioLCPR_Sim, bins=30, histtype="step", range = (0.01,1), lw = 2)
plt.title("LC NHits > 2 Energy / SimTrackster Energy")
plt.xlabel("Ratio Energy")
plt.savefig(pathPlot + "RatioLCPRvsSimTrk" + ".png")
# plt.show()


plt.figure(figsize=(10, 10))
plt.hist(N_trackster, bins=10, histtype="step", range = (0,10), label = 'MERGE', lw = 2)
plt.hist(N_tracksterHAD, bins = 10, histtype='step', range = (0,10), label = 'HAD', lw = 2)
plt.hist(N_tracksterEM, bins = 10, histtype='step', range = (0,10), label = 'EM', lw = 2)
plt.title(f"Number of Reco-Tracksters: Total {sum(N_trackster)}")
plt.xlabel("Number of Tracksters")
plt.legend(loc = 'upper right')
plt.savefig(pathPlot + "NTracksters" + ".png")


plt.figure(figsize=(10, 10))
plt.hist(N_lc, bins=50, histtype="step", range = (1,max(N_lc)), label = 'MERGE', lw = 2)
plt.hist(N_lcHAD, bins = 50, histtype='step', range = (1,max(N_lc)), label = 'HAD', lw = 2)
plt.hist(N_lcEM, bins = 50, histtype='step', range = (1,max(N_lc)), label = 'EM', lw = 2)
plt.title(f"Number of LC in trackster: Total {sum(N_lc)}")
plt.xlabel("Number of LCs")
plt.legend(loc = 'upper right')
plt.savefig(pathPlot + "NLCs" + ".png")


plt.figure(figsize=(10, 10))
plt.hist(trackstersSkippedLayers, bins=5, histtype="step", range = (0,5), label = 'MERGE', lw = 2)
plt.hist(trackstersSkippedLayersEM, bins=5, histtype="step", range = (0,5), label = 'EM', lw = 2)
plt.hist(trackstersSkippedLayersHAD, bins=5, histtype="step", range = (0,5), label = 'HAD', lw =2)
plt.title("Number of Layers Skipped")
plt.xlabel("Skipped Layers")
plt.legend(loc = 'upper right')
plt.savefig(pathPlot + "SkippedLayers" + ".png")

plt.figure(figsize=(10, 10))
plt.hist(trackstersLenght, bins=50, histtype="step", range = (0,50), label = 'MERGE', lw = 2)
plt.hist(trackstersLenghtEM, bins=50, histtype="step", range = (0,50), label = 'MERGE', lw = 2)
plt.hist(trackstersLenghtHAD, bins=50, histtype="step", range = (0,50), label = 'MERGE', lw = 2)
plt.title("Number of Layers in Trackster")
plt.xlabel("Number of Tracksters Layer")
plt.legend(loc = 'upper right')
plt.savefig(pathPlot + "NLayerTracksters" + ".png")