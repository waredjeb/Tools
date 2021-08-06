from utils import *
import uproot
import matplotlib.pyplot as plt
import ROOT
import numpy as np
import pandas as pd
import mplhep as hep

plt.style.use([hep.style.ROOT, hep.style.firamath])
particle = "Pion"
name = "Siblings_25-30deg"
name2 = "NoSiblings_25-30deg"
fileNames = ["./data/output_CloseBy"+particle+"_"+name+".root","./data/output_CloseBy"+particle+"_"+name2+".root"]
# pathPlot = "./plots/Comparison/CloseBy"+particle+"/Siblings_20-25degVsNoSibling_20-25deg/"
pathPlot = "./plots/Comparison/CloseBy"+particle+"/"+name+"VS"+name2+"/"
mkdir_p(pathPlot)


labels = [name, name2]

profXMulti = list()
lcNoUsed_flatMulti = list()
lcNoPR_energiesMulti = list()
profXHDMulti = list()
profXLDMulti = list()
totalTracksterEnergy_rawMulti = list()
totalTracksterEnergy_regressedMulti = list()

ratioEM_simMulti = list()
ratioHad_simMulti = list()

ratioEM_regressedMulti = list()
ratioTrk_regressedMulti = list()
ratioHad_regressedMulti = list()

ratioEM_lcPRMulti = list()
ratioHad_lcPRMulti = list()

ratioTot_sim_rawMulti = list()
ratioTot_sim_regressedMulti = list()

ratioTot_lcPRMulti = list()
ratioTot_sim_regressedMulti = list()

ratioLCPR_SimMulti = list()
ratioLCPR_LCNOPR_SIMMulti = list()
ratioLCNOPR_SimRawMulti = list()

ratioLCNOPR_CPMulti = list()

ratioLinking_simMulti = list()
ratioLinking_regressedMulti = list()
ratioLinking_lcPRMulti = list()
N_trackstersMulti = list()
N_trackstersMultiEM = list()
N_trackstersMultiHAD = list()

N_lcMulti = list()
N_lcEMMulti = list()
N_lcHADMulti = list()

trackstersLenghtMulti = list()
trackstersSkippedLayersMulti = list()
trackstersLenghtMultiEM = list()
trackstersSkippedLayersMultiEM = list()
trackstersLenghtMultiHAD = list()
trackstersSkippedLayersMultiHAD = list()

histo = ROOT.TH2F("Profile", "Profile", 50, 0, 50, 50, 0, 1)
histo2 = ROOT.TH2F("Profile2", "Profile2", 50, 0, 50, 50, 0, 1)

f_num = 0
for fileName in fileNames:
    file = uproot.open(fileName)
    keys = file.keys()

    fileRoot = file[keys[1]]

    titleSave = "_" + fileName[22:-5]



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
    ratioHad_sim = list()

    ratioEM_regressed = list()
    ratioTrk_regressed = list()
    ratioHad_regressed = list()

    ratioEM_lcPR = list()
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

    N_tracksters = list()
    N_trackstersEM = list()
    N_trackstersHAD = list()

    N_lc = list()
    N_lcEM = list()
    N_lcHAD = list()

    print(f"Total Events {labels[f_num]} {len(sim_rawEnergy)}")
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

    totalTracksterEnergy_rawMulti.append(totalTracksterEnergy_raw)
    totalTracksterEnergy_regressedMulti.append(totalTracksterEnergy_regressed)
    ratioEM_simMulti.append(ratioEM_sim)
    ratioHad_simMulti.append(ratioHad_sim)
    ratioTot_sim_rawMulti.append(ratioTot_sim_raw)
    ratioEM_regressedMulti.append(ratioEM_regressed)
    ratioHad_regressedMulti.append(ratioHad_regressed)
    ratioTot_sim_regressedMulti.append(ratioTot_sim_regressed)
    ratioEM_lcPRMulti.append(ratioEM_lcPR)
    ratioHad_lcPRMulti.append(ratioHad_lcPR)
    ratioTot_lcPRMulti.append(ratioTot_lcPR)
    ratioLCPR_SimMulti.append(ratioLCPR_Sim)
    ratioLCPR_LCNOPR_SIMMulti.append(ratioLCPR_LCNOPR_SIM)
    ratioLCNOPR_SimRawMulti.append(ratioLCNOPR_SimRaw)
    ratioLCNOPR_CPMulti.append(ratioLCNOPR_CP)

    for ev in range(len(n_tracksters)):
        for n in n_tracksters[ev]:
            N_tracksters.append(n)
        for n in n_trackstersHAD[ev]:
            N_trackstersHAD.append(n)
        for n in n_trackstersEM[ev]:
            N_trackstersEM.append(n)

    N_trackstersMulti.append(N_tracksters)
    N_trackstersMultiEM.append(N_trackstersEM)
    N_trackstersMultiHAD.append(N_trackstersHAD)

    for ev in range(len(n_lc)):
        for n in n_lc[ev]:
            N_lc.append(n)
        for n in n_lcEM[ev]:
            N_lcEM.append(n)
        for n in n_lcHAD[ev]:
            N_lcHAD.append(n)

    N_lcMulti.append(N_lc)
    N_lcEMMulti.append(N_lcEM)
    N_lcHADMulti.append(N_lcHAD)


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


    noPrFullEnergyRatioMulti = list()
    for ev in range(len(fullLCEnergy)):
        for z in range(len(zs)):
            if fullLCEnergy[ev][z] > 0:
                if(f_num == 0):
                    histo.Fill(zs[z], noPRLCEnergy[ev][z] / fullLCEnergy[ev][z])
                else:
                    histo2.Fill(zs[z], noPRLCEnergy[ev][z] / fullLCEnergy[ev][z])

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


    lcNoPR_energies = list()
    for i in lcNoPR_spectrum:
        for j in i:
            lcNoPR_energies.append(j)

    lcNoPR_energiesMulti.append(lcNoPR_energies)

    lcNoUsed_flat = []
    for i in lc_noUsed:
        lcNoUsed_flat.append(i[0])

    lcNoUsed_flatMulti.append(lcNoUsed_flat)

    f_num += 1

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

    trackstersLenghtMulti.append(trackstersLenght)
    trackstersSkippedLayersMulti.append(trackstersSkippedLayers)

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

    trackstersLenghtMultiEM.append(trackstersLenghtEM)
    trackstersSkippedLayersMultiEM.append(trackstersSkippedLayersEM)

    

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

    trackstersLenghtMultiHAD.append(trackstersLenghtHAD)
    trackstersSkippedLayersMultiHAD.append(trackstersSkippedLayersHAD)

c = ROOT.TCanvas("c", "c", 1000,1000)
leg = ROOT.TLegend(0.1,0.7,0.3,0.9)

profX1 = histo.ProfileX()
profX2 = histo2.ProfileX()
profX1.GetYaxis().SetRangeUser(0, 1)
profX1.SetMarkerColor(4)
profX1.SetMarkerSize(1.5)
profX1.SetMarkerStyle(20)
profX1.SetStats(0000)
profX2.SetMarkerColor(ROOT.kRed)
profX2.SetMarkerSize(1.5)
profX2.SetMarkerStyle(20)
profX2.SetStats(0000)
leg.AddEntry(profX1, labels[0])
leg.AddEntry(profX2, labels[1])
profX1.SetTitle("LC NHits < 3 / All LCs ")
profX1.GetXaxis().SetTitle("HGCal Layer")
profX1.GetYaxis().SetTitle("Ratio")
profX1.Draw("P")
profX2.Draw("SAME P")
leg.Draw("SAME")
c.SaveAs(pathPlot + "ProfileMulti.png")

##############################################################
# LC not in pattern recognition energy spectrumm
##############################################################
plt.figure(figsize = (10,10))
plt.hist(lcNoPR_energiesMulti, bins = 100, histtype='step', label=labels, lw = 2, range = (0,10))
plt.legend()
plt.xlabel("Energy [GeV]")
plt.ylabel("Entries")
plt.title("LCs Nhits < 3 Energy Spectrum")
plt.savefig(pathPlot + "LCNoPrEnergies.png")

##############################################################
# RECO-TRACKSTER RAW-ENERGY WRT SIMTRACKSTER REGRESSED ENERGY
##############################################################
plt.figure(figsize = (10,10))
plt.hist(ratioTot_sim_rawMulti, bins = 40, histtype='step', label=labels, lw = 2, range = (0.01, 1))
plt.hist(ratioEM_simMulti, bins = 40, histtype='step', label=["EM-" + i for i in labels], lw = 2,range = (0.01, 1))
plt.hist(ratioHad_simMulti, bins = 40, histtype='step', label=["HAD-" + i for i in labels], lw = 2,range = (0.01, 1))
plt.legend(loc = 'upper left')
plt.xlabel("Ratio")
plt.ylabel("Entries")
plt.title("Ratio Reco-Trackster Raw-Energy wrt SimTrackster Regressed Energy", fontsize = 20)
plt.savefig(pathPlot + "RatioRecovsSimRegressed.png")

plt.figure(figsize = (10,10))
plt.hist(ratioTot_sim_rawMulti, bins = 40, histtype='step', label=labels, lw = 2, range = (0.01, 1))
plt.legend(loc = 'upper left')
plt.xlabel("Ratio")
plt.ylabel("Entries")
plt.title("Ratio Reco-Trackster Raw-Energy wrt SimTrackster Regressed Energy", fontsize = 20)
plt.savefig(pathPlot + "RatioRecovsSimRegressedAll.png")

##############################################################
# LC in PR WRT SIMTRACKSTER REGRESSED ENERGY
##############################################################
plt.figure(figsize = (10,10))
plt.hist(ratioTot_lcPRMulti, bins = 40, histtype='step', label=labels, lw = 2, range = (0.01, 1))
plt.hist(ratioEM_lcPRMulti, bins = 40, histtype='step', label=["EM-" + i for i in labels], lw = 2,range = (0.01, 1))
plt.hist(ratioHad_lcPRMulti, bins = 40, histtype='step', label=["HAD-" + i for i in labels], lw = 2,range = (0.01, 1))
plt.legend(loc = 'upper left')
plt.xlabel("Ratio")
plt.ylabel("Entries")
plt.title("Ratio LC Nhits > 2 Total Energy wrt SimTrackster Regressed Energy", fontsize = 20)
plt.savefig(pathPlot + "RatioLCPRvsSimRegressed.png")

plt.figure(figsize = (10,10))
plt.hist(ratioTot_lcPRMulti, bins = 40, histtype='step', label=labels, lw = 2, range = (0.01, 1))
plt.legend(loc = 'upper left')
plt.xlabel("Ratio")
plt.ylabel("Entries")
plt.title("Ratio LC Nhits > 2 Total Energy wrt SimTrackster Regressed Energy", fontsize = 20)
plt.savefig(pathPlot + "RatioLCPRvsSimRegressedAll.png")


##############################################################
# N reco Tracksters
##############################################################

plt.figure(figsize = (10,10))
labels_new = [labels[i] + " [ " + str(sum(N_trackstersMulti[i])) + "]" for i in range(len(labels))]
plt.hist(N_trackstersMulti, bins = 10, histtype='step', label=labels_new, lw = 2, range = (0.0, 10))
plt.legend(loc = 'upper right')
plt.xlabel("Number of Tracksters")
plt.ylabel("Entries")
plt.title("Number of Reco-Tracksters", fontsize = 20)
plt.savefig(pathPlot + "NTracksters.png")


##############################################################
# N reco Tracksters EM
##############################################################

plt.figure(figsize = (10,10))
labels_new = [labels[i] + " [ " + str(sum(N_trackstersMultiEM[i])) + "]" for i in range(len(labels))]
plt.hist(N_trackstersMultiEM, bins = 10, histtype='step', label=labels_new, lw = 2, range = (0.0, 10))
plt.legend(loc = 'upper right')
plt.xlabel("Number of Tracksters")
plt.ylabel("Entries")
plt.title("Number of Reco-Tracksters EM", fontsize = 20)
plt.savefig(pathPlot + "NTrackstersEM.png")

##############################################################
# N reco Tracksters HAD
##############################################################

plt.figure(figsize = (10,10))
labels_new = [labels[i] + " [ " + str(sum(N_trackstersMultiHAD[i])) + "]" for i in range(len(labels))]
plt.hist(N_trackstersMultiHAD, bins = 10, histtype='step', label=labels_new, lw = 2, range = (0.0, 10))
plt.legend(loc = 'upper right')
plt.xlabel("Number of Tracksters")
plt.ylabel("Entries")
plt.title("Number of Reco-Tracksters HAD", fontsize = 20)
plt.savefig(pathPlot + "NTrackstersHAD.png")


##############################################################
# N LC Tracksters
##############################################################

plt.figure(figsize = (10,10))
labels_new = [labels[i] + " [ " + str(sum(N_lcMulti[i])) + "]" for i in range(len(labels))]
plt.hist(N_lcMulti, bins = 50, histtype='step', label=labels_new, lw = 2, range = (0.0, 250))
plt.legend(loc = 'upper right')
plt.xlabel("Number of LC")
plt.ylabel("Entries")
plt.title("Number of LC in Tracksters", fontsize = 20)
plt.savefig(pathPlot + "NLC.png")


##############################################################
# N LC Tracksters EM
##############################################################

plt.figure(figsize = (10,10))
labels_new = [labels[i] + " [ " + str(sum(N_lcEMMulti[i])) + "]" for i in range(len(labels))]
plt.hist(N_lcEMMulti, bins = 50,  histtype='step', label=labels_new, lw = 2, range = (0.0, 250))
plt.legend(loc = 'upper right')
plt.xlabel("Number of LC")
plt.ylabel("Entries")
plt.title("Number of LC in Tracksters EM", fontsize = 20)
plt.savefig(pathPlot + "NLCEM.png")

##############################################################
# N LC Tracksters HAD
##############################################################

plt.figure(figsize = (10,10))
labels_new = [labels[i] + " [ " + str(sum(N_lcHADMulti[i])) + "]" for i in range(len(labels))]
plt.hist(N_lcHADMulti, bins = 50, histtype='step', label=labels_new, lw = 2, range = (0.0, 250))
plt.legend(loc = 'upper right')
plt.xlabel("Number of LC")
plt.ylabel("Entries")
plt.title("Number of LC in Tracksters HAD", fontsize = 20)
plt.savefig(pathPlot + "NLCHad.png")


##############################################################
# Skipped LCs 
##############################################################

plt.figure(figsize = (10,10))
# labels_new = [labels[i] + " [ " + str(sum(N_lcMulti[i])) + "]" for i in range(len(labels))]
plt.hist(trackstersSkippedLayers, bins = 6, histtype='step', label=labels, lw = 2, range = (0.0, 6))
plt.legend(loc = 'upper right')
plt.xlabel("Number of Skipped LCs")
plt.ylabel("Entries")
plt.title("Number Skipped Layers in Tracksters", fontsize = 20)
plt.savefig(pathPlot + "SkippedLCs.png")


##############################################################
# Skipped LCs  EM
##############################################################

plt.figure(figsize = (10,10))
# labels_new = [labels[i] + " [ " + str(sum(N_lcMulti[i])) + "]" for i in range(len(labels))]
plt.hist(trackstersSkippedLayersEM, bins = 6, histtype='step', label=labels, lw = 2, range = (0.0, 6))
plt.legend(loc = 'upper right')
plt.xlabel("Number of Skipped LCs")
plt.ylabel("Entries")
plt.title("Number Skipped Layers in Tracksters EM", fontsize = 20)
plt.savefig(pathPlot + "SkippedLCsEM.png")

##############################################################
# Skipped LCs HAD
##############################################################

plt.figure(figsize = (10,10))
# labels_new = [labels[i] + " [ " + str(sum(N_lcMulti[i])) + "]" for i in range(len(labels))]
plt.hist(trackstersSkippedLayersHAD, bins = 6, histtype='step', label=labels, lw = 2, range = (0.0, 6))
plt.legend(loc = 'upper right')
plt.xlabel("Number of Skipped LCs")
plt.ylabel("Entries")
plt.title("Number Skipped Layers in Tracksters HAD", fontsize = 20)
plt.savefig(pathPlot + "SkippedLCsHAD.png")


##############################################################
# Skipped LCs 
##############################################################

plt.figure(figsize = (10,10))
# labels_new = [labels[i] + " [ " + str(sum(N_lcMulti[i])) + "]" for i in range(len(labels))]
plt.hist(trackstersSkippedLayersMulti, bins = 6, histtype='step', label=labels, lw = 2, range = (0.0, 6))
plt.legend(loc = 'upper right')
plt.xlabel("Number of Skipped LCs")
plt.ylabel("Entries")
plt.title("Number Skipped Layers in Tracksters", fontsize = 20)
plt.savefig(pathPlot + "SkippedLCs.png")


##############################################################
# Skipped LCs  EM
##############################################################

plt.figure(figsize = (10,10))
# labels_new = [labels[i] + " [ " + str(sum(N_lcMulti[i])) + "]" for i in range(len(labels))]
plt.hist(trackstersSkippedLayersMultiEM, bins = 6, histtype='step', label=labels, lw = 2, range = (0.0, 6))
plt.legend(loc = 'upper right')
plt.xlabel("Number of Skipped LCs")
plt.ylabel("Entries")
plt.title("Number Skipped Layers in Tracksters EM", fontsize = 20)
plt.savefig(pathPlot + "SkippedLCsEM.png")

##############################################################
# Skipped LCs HAD
##############################################################

plt.figure(figsize = (10,10))
# labels_new = [labels[i] + " [ " + str(sum(N_lcMulti[i])) + "]" for i in range(len(labels))]
plt.hist(trackstersSkippedLayersMultiHAD, bins = 6, histtype='step', label=labels, lw = 2, range = (0.0, 6))
plt.legend(loc = 'upper right')
plt.xlabel("Number of Skipped LCs")
plt.ylabel("Entries")
plt.title("Number Skipped Layers in Tracksters HAD", fontsize = 20)
plt.savefig(pathPlot + "SkippedLCsHAD.png")

##############################################################
# Skipped LCs 
##############################################################

plt.figure(figsize = (10,10))
# labels_new = [labels[i] + " [ " + str(sum(N_lcMulti[i])) + "]" for i in range(len(labels))]
plt.hist(trackstersLenghtMulti, bins = 50, histtype='step', label=labels, lw = 2, range = (0.0, 50))
plt.legend(loc = 'upper right')
plt.xlabel("Number of LC")
plt.ylabel("Entries")
plt.title("Number of LC in Tracksters", fontsize = 20)
plt.savefig(pathPlot + "TrackstersLength.png")


##############################################################
# Skipped LCs  EM
##############################################################
plt.figure(figsize = (10,10))
# labels_new = [labels[i] + " [ " + str(sum(N_lcMulti[i])) + "]" for i in range(len(labels))]
plt.hist(trackstersLenghtMultiEM, bins = 50, histtype='step', label=labels, lw = 2, range = (0.0, 50))
plt.legend(loc = 'upper right')
plt.xlabel("Number of LC")
plt.ylabel("Entries")
plt.title("Number of LC in Tracksters EM ", fontsize = 20)
plt.savefig(pathPlot + "TrackstersLengthEM.png")

##############################################################
# Skipped LCs HAD
##############################################################
plt.figure(figsize = (10,10))
# labels_new = [labels[i] + " [ " + str(sum(N_lcMulti[i])) + "]" for i in range(len(labels))]
plt.hist(trackstersLenghtMultiHAD, bins = 50, histtype='step', label=labels, lw = 2, range = (0.0, 50))
plt.legend(loc = 'upper right')
plt.xlabel("Number of LC")
plt.ylabel("Entries")
plt.title("Number of LC in Tracksters HAD", fontsize = 20)
plt.savefig(pathPlot + "TrackstersLengthHAD.png")