from builtins import range
import ROOT
import sys
import os 
from DataFormats.FWLite import Events, Handle
import matplotlib.pyplot as plt
import mplhep as hep
from utils import *
import seaborn as sns
import numpy as np


plt.style.use(hep.style.CMS)

# from FWCore.ParameterSet.VarParsing import VarParsing
# options = VarParsing ('python')
# options.parseArguments()

# Events takes either
# - single file name
# - list of file names
# - VarParsing options

# use Varparsing object

#Create new directory
def mkdir_p(mypath):
    '''Function to create a new directory, if it not already exist
        - mypath : directory path
    '''
    from errno import EEXIST
    from os import makedirs,path
    try:
        makedirs(mypath)
    except OSError as exc: 
        if exc.errno == EEXIST and path.isdir(mypath):
            pass
        else: raise

def plotHisto(list, bins, title, xlabel, savename, ylabel = "Entries"):
    plt.figure(figsize = (10,10))
    plt.hist(list, bins = bins )
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig(savename + ".png")

def matchlc(lc1, lc2):
    if(lc1.x() == lc2.x() and lc1.y() == lc2.y() and lc1.z() == lc2.z()):
        return True
    else:
        return False

name_np_missing = sys.argv[1]
# events = Events ("/eos/user/w/wredjeb/HGCAL/TrackstersSeed/Samples/NoOverlapping/singlephoton_closeBy_deltaR_hgcalCenter/step3/step3_singlephoton_Delta1_nopu.root")
# events = Events ("/afs/cern.ch/work/w/wredjeb/public/TICLv4/Seeding/CMSSW_12_3_0_pre1/src/hgcal_private/production/Test/singlepi_closeBy_deltaR_hgcalCenter/step3/step3_singlepi_Delta1_nopu.root")
# events = Events ("/afs/cern.ch/work/w/wredjeb/public/TICLv4/Seeding/CMSSW_12_3_0_pre1/src/hgcal_private/production/Test/new/singlepi_closeBy_deltaR_hgcalCenter/step3/step3_singlepi_Delta1_nopu.root")
# events = Events("/eos/user/w/wredjeb/HGCAL/TrackstersSeed/Samples/singlephoton_closeBy_deltaR_hgcalCenter/step3/step3_singlephoton_Delta10_nopu.root")
# events = Events ("/afs/cern.ch/work/w/wredjeb/public/TICLv4/Seeding/CMSSW_12_3_0_pre1/src/hgcal_private/production/Test/newNoFlatEta/singlepi_closeBy_deltaR_hgcalCenter/step3/step3_singlepi_Delta1_nopu.root")
# events = Events("/afs/cern.ch/work/w/wredjeb/public/TICLv4/Seeding/CMSSW_12_3_0_pre1/src/hgcal_private/production/Test/new2/singlephoton_closeBy_deltaR_hgcalCenter/step3/step3_singlephoton_Delta10_nopu.root")
file_in = "/afs/cern.ch/work/w/wredjeb/public/TICLv4/FineSimTracksters/CMSSW_12_3_0_pre4/src/34693.0_CloseByParticleGun+2026D76+CloseByParticle_Photon_ERZRanges_GenSimHLBeamSpotHGCALCloseBy+DigiTrigger+RecoGlobal+HARVESTGlobal/step3_"+name_np_missing+".root"
print(f"Processing file {file_in}")
events = Events(file_in)

outputdir = './fineSimTracksters_missing_layers/' + name_np_missing + "/"
print(f"OUTPUT DIR {outputdir}")
outputdir3D = outputdir + 'projections/'
mkdir_p(outputdir)
mkdir_p(outputdir3D)
output_np = outputdir + "/output_np/"
mkdir_p(output_np)
caloparticlesH = Handle("std::vector<CaloParticle")
cpslabel = ("mix", "MergedCaloTruth")

gpH = Handle("std::vector<reco::GenParticle")
gpL = ("genParticles")

simTracksterToFineSimtracksterMap_handle = Handle("std::map<unsigned int,std::vector<unsigned int> > ")
simTracksterToFineSimtracksterMap_label = ("ticlFineSimTracksters", "fine")

layerClusters  = Handle("std::vector<reco::CaloCluster>")
layerClustersLabel = ("hgcalLayerClusters")

clue3DLowseedsH = Handle("std::vector<int>")
clue3DHighseedsH = Handle("std::vector<int>")

clue3DLowseedsLabel = ("ticlTrackstersCLUE3DHigh", "tracksterSeeds")
clue3DHighseedsLabel = ("ticlTrackstersCLUE3DHigh", "tracksterSeeds")

fineSimTrackstersSeedL = ("ticlFineSimTracksters", "fine")
fineSimTrackstersSeedDoubletsL = ("ticlFineSimTracksters", "tracksterSeedsDoublets")
fineSimTrackstersSeedH = Handle("std::vector<int>")
# fineSimTrackstersSeedDoubletsH = Handle("std::vector<std::vector<int>>")

clue3DLowseedsDoubletsH = Handle("std::vector<std::vector<int>>")
clue3DHighseedsDoubletsH = Handle("std::vector<std::vector<int>>")
clue3DLowseedsDoubletsLabel = ("ticlTrackstersCLUE3DLow", "tracksterSeedsDoublets")
clue3DHighseedsDoubletsLabel = ("ticlTrackstersCLUE3DHigh", "tracksterSeedsDoublets")

clue3DLowTrackstersH = Handle("std::vector<ticl::Trackster>")
clue3DHighTrackstersH = Handle("std::vector<ticl::Trackster>")
simtrackstersH = Handle("std::vector<ticl::Trackster>")
simtrackstersL = ("ticlSimTracksters")
finesimtrackstersH = Handle("std::vector<ticl::Trackster>")
finesimtrackstersL = ("ticlFineSimTracksters", "fine")
clue3DLowTrackstersLabel = ("ticlTrackstersCLUE3DLow")
clue3DHighTrackstersLabel = ("ticlTrackstersCLUE3DHigh")

CATracksterMergeH = Handle("std::vector<ticl::Trackster>")
CATracksterMergeLabel = ("ticlTrackstersMerge")

number_seedsLow = 0
number_seedsHigh = 0

number_tracksterLow = 0
number_tracksterHigh = 0

number_CAtracksters = 0

seeds_eta = []
seeds_phi = []
seeds_z = []
seeds_energy = []


lcs_eta = []
lcs_phi = []
lcs_z = []
lcs_energy = []

cp_eta = []
cp_phi = []
cp_energy = []

gp_eta = []
gp_phi = []
gp_energy = []

lcs_to_plot = []
lcsE_to_plot = []
simlcs_to_plot = []
finesimlcsE_to_plot = []
seeds_to_plot = []

seeds_merge = []

clue3DTracksterMerge = []

simTrackster_energy = []
fineSimTrackster_energy = []

simTrackster_NLCs = []
finesimTrackster_NLCs = []

n_fineSimTrackster_per_st = []

mergedFineSimTrackster_energy = []
number_of_missing_lcs_per_ev = []
number_of_slimtrackster_per_ev = []
number_of_lcs_per_slimtrackster = []
for i_ev, ev in enumerate(events):
    
    if(i_ev % 1 == 0):
        print(f"-- Processed events {i_ev}")
    ev.getByLabel (gpL, gpH)
    ev.getByLabel (cpslabel, caloparticlesH)
    ev.getByLabel (layerClustersLabel, layerClusters)
    ev.getByLabel (simtrackstersL, simtrackstersH)
    ev.getByLabel (finesimtrackstersL, finesimtrackstersH)
    ev.getByLabel (clue3DLowseedsLabel, clue3DLowseedsH)
    ev.getByLabel (fineSimTrackstersSeedL, fineSimTrackstersSeedH)
    # ev.getByLabel (fineSimTrackstersSeedDoubletsfineSimTrackstersSeedDoubletsLL, fineSimTrackstersSeedDoubletsL)
    ev.getByLabel (clue3DHighseedsLabel, clue3DHighseedsH)
    ev.getByLabel (clue3DLowseedsDoubletsLabel, clue3DLowseedsDoubletsH)
    ev.getByLabel (clue3DHighseedsDoubletsLabel, clue3DHighseedsDoubletsH)
    ev.getByLabel (clue3DLowTrackstersLabel, clue3DLowTrackstersH)
    ev.getByLabel (clue3DHighTrackstersLabel, clue3DHighTrackstersH)
    ev.getByLabel (CATracksterMergeLabel, CATracksterMergeH)
    ev.getByLabel (simTracksterToFineSimtracksterMap_label, simTracksterToFineSimtracksterMap_handle)

    stTofstMap = simTracksterToFineSimtracksterMap_handle.product()
    simTracksters = simtrackstersH.product()
    finesimTracksters = finesimtrackstersH.product()
    lcs = layerClusters.product()
    gps = gpH.product()
    seedsLow = clue3DLowseedsH.product()
    seedsHigh = clue3DHighseedsH.product()
    fineSeeds = fineSimTrackstersSeedH.product()
    s_doubletsLow = clue3DLowseedsDoubletsH.product()
    s_doubletsHigh = clue3DHighseedsDoubletsH.product()
    cps = caloparticlesH.product()
    number_of_slimtrackster_per_ev.append(len(finesimTracksters))
    if(clue3DLowTrackstersH.isValid()):
        tracksterLow = clue3DLowTrackstersH.product()
    else:
        tracksterLow = []

    if(clue3DHighTrackstersH.isValid()):
        tracksterHigh = clue3DHighTrackstersH.product()
    else:
        tracksterHigh = []
    CAtrackster = CATracksterMergeH.product()

    sim_tot_number_of_lcs = 0
    list_lcs_simtracksters = []
    for st in simTracksters:
        simTrackster_energy.append(st.raw_energy())
        simTrackster_NLCs.append(st.vertices().size())
        sim_tot_number_of_lcs += st.vertices().size()
        for i in range(st.vertices().size()):
            list_lcs_simtracksters.append(st.vertices(i))

    fine_tot_number_of_lcs = 0.
    list_lcs_finesimtracksters = []
    lcs_fineSimTracksters = []
    lcsE_fineSimTracksters = []
    lcsETrkIdx_fineSimTracksters = []
    for i_fst, fst in enumerate(finesimTracksters):
        fineSimTrackster_energy.append(fst.raw_energy())
        finesimTrackster_NLCs.append(fst.vertices().size())
        fine_tot_number_of_lcs += fst.vertices().size()
        number_of_lcs_per_slimtrackster.append(fst.vertices().size())
        for i in range(fst.vertices().size()):
            list_lcs_finesimtracksters.append(fst.vertices(i))
            lcs_fineSimTracksters.append(lcs[fst.vertices(i)])
            lcsE_fineSimTracksters.append(lcs[fst.vertices(i)].energy())
            lcsETrkIdx_fineSimTracksters.append(i_fst)

    lcs_missing_to_plot = []
    lcs_missing_to_plotE = []
    fineSeeds_list = []
    fineSeedsE_list = []

    print(f"Difference seeds {len(fineSeeds)} and fineSimTacksters {len(finesimTracksters)} : {len(fineSeeds) - len(finesimTracksters) }")
    for i_s, s in enumerate(fineSeeds):
        fineSeeds_list.append(lcs[s])
        fineSeedsE_list.append(lcs[s].energy())
    diff = sim_tot_number_of_lcs - fine_tot_number_of_lcs
    list_lcs_finesimtracksters_unique = np.unique(np.array(list_lcs_finesimtracksters))
    list_lcs_simtracksters_unique = np.unique(np.array(list_lcs_simtracksters))
    for lc_id in list_lcs_simtracksters_unique:
        if((lc_id not in list_lcs_finesimtracksters_unique)):
            lcs_missing_to_plot.append(lcs[int(lc_id)])    
            lcs_missing_to_plotE.append(lcs[int(lc_id)].energy())    
    print(f"Difference is {len(list_lcs_finesimtracksters_unique) - len(list_lcs_simtracksters_unique)} vs {len(lcs_missing_to_plot)}")
    if(len(lcs_missing_to_plot) >= 10):
        plots3DwithProjectionSeeds3D_2(lcs_fineSimTracksters, lcsE_fineSimTracksters, lcsETrkIdx_fineSimTracksters,lcs_missing_to_plot, lcs_missing_to_plotE, [], fineSeeds_list ,fineSeedsE_list, "SlimSimTracksters and Missing LCs", outputdir3D + 'fineSimtrackster_'+str(i_ev) + '.png', map_ = False)
    number_of_missing_lcs_per_ev.append(len(lcs_missing_to_plot))
    for m in stTofstMap:
        tot_merged_energy = 0.
        n_fineSimTrackster_per_st.append(len(m.second))
        for i_fst in m.second:
            tot_merged_energy += finesimTracksters[i_fst].raw_energy()
        mergedFineSimTrackster_energy.append(tot_merged_energy)
    
    assert len(mergedFineSimTrackster_energy) == len(simTrackster_energy), "Number of mergedFineSimTrackste != Number of simTracksters"

np_number_of_missing_lcs_per_ev = np.array(number_of_missing_lcs_per_ev)
np.save(output_np + "/missing_lcs_per_ev_" + name_np_missing + ".npy",  np_number_of_missing_lcs_per_ev)

np_number_of_slimtrackster_per_ev = np.array(number_of_slimtrackster_per_ev)
np.save(output_np + "/slim_simTrackster_per_ev_" + name_np_missing + ".npy",  np_number_of_slimtrackster_per_ev)

np_number_of_lcs_per_slimtrackster = np.array(number_of_lcs_per_slimtrackster)
np.save(output_np + "/number_of_lcs_per_slimtrackster_" + name_np_missing + ".npy",  np_number_of_lcs_per_slimtrackster)


text = r"2 Pions CloseBy $\Delta$R = 10"
text += "\n"
text += "Energy $\in$ [25,200] GeV"
text += "\n 0PU"

xy_text = (0.7, 0.5)
xy_coord = 'figure fraction'

output_path = "./FineSimTrackstersPlot/"
plt.figure(figsize = (12,10))
plt.hist2d(mergedFineSimTrackster_energy, simTrackster_energy, bins = (30,30))
plt.xlabel("Merged FineSimTrackster Energy")
plt.ylabel("SimTrackster Energy")

plt.savefig(output_path + "FineSimTrackstersMergeVSSimTracksters_Energy.png")

plt.figure(figsize = (12,10))

sns.regplot(simTrackster_energy, mergedFineSimTrackster_energy, x_bins = 15, fit_reg = None)
plt.ylim(0, 210)
plt.ylabel("Merged FineSimTrackster Energy")
plt.xlabel("SimTrackster Energy")
plt.annotate(text, xy= xy_text,  xycoords='figure fraction',
            xytext=xy_text, textcoords='figure fraction', size=10)
plt.savefig(output_path + "FineSimTrackstersMergeVSSimTracksters_Energy_profile.png")

plt.figure(figsize = (12,10))
plt.hist(simTrackster_energy, bins = 80, alpha = 0.5, label = "SimTracksters")
plt.hist(fineSimTrackster_energy, bins = 80, alpha = 0.5,  label  = "FineSimTracksters")
plt.xlabel("Raw Energy")
plt.ylabel("Entries")
plt.title("Energy")
plt.annotate(text, xy= xy_text,  xycoords='figure fraction',
            xytext=xy_text, textcoords='figure fraction', size=10)
plt.legend()
plt.savefig(output_path + "SimTracksterAndFineSimTracksterEnergy.png")

plt.figure(figsize = (12,10))
plt.hist(simTrackster_energy, bins = 80, alpha = 0.5, label = "SimTracksters")
plt.hist(mergedFineSimTrackster_energy, bins = 80, alpha = 0.5,  label  = "MergedFineSimTracksters")
plt.xlabel("Raw Energy")
plt.ylabel("Entries")
plt.title("Energy")
plt.annotate(text, xy= xy_text,  xycoords='figure fraction',
            xytext=xy_text, textcoords='figure fraction', size=10)
plt.legend()
plt.savefig(output_path + "SimTracksterAndMergedFineSimTracksterEnergy.png")

plt.figure(figsize = (12,10))
plt.hist(n_fineSimTrackster_per_st, bins = 80, alpha = 0.5)
plt.xlabel("Number of FineSimTracksters per SimTrackster")
plt.ylabel("Entries")
plt.title("Energy")
plt.annotate(text, xy= xy_text,  xycoords='figure fraction',
            xytext=xy_text, textcoords='figure fraction', size=10)
plt.savefig(output_path + "NumberOfFineSimTracksterPerSimTracksters.png")

plt.figure(figsize = (12,10))
plt.hist(simTrackster_NLCs, bins = 80, alpha = 0.5, label = "SimTracksters")
plt.hist(finesimTrackster_NLCs, bins = 80, alpha = 0.5,  label  = "FineSimTracksters")
plt.xlabel("Number of Layer Clusters in Tracksters")
plt.ylabel("Entries")
plt.title("Number of Layer Clusters in Tracksters")
plt.annotate(text, xy= xy_text,  xycoords='figure fraction',
            xytext=xy_text, textcoords='figure fraction', size=10)
plt.legend()
plt.savefig(output_path + "NumberOfLayerClustersInTracksters.png")



