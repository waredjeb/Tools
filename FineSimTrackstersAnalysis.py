from builtins import range
import ROOT
import sys
import os 
from DataFormats.FWLite import Events, Handle
import matplotlib.pyplot as plt
import mplhep as hep
from utils import *

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

# events = Events ("/eos/user/w/wredjeb/HGCAL/TrackstersSeed/Samples/NoOverlapping/singlephoton_closeBy_deltaR_hgcalCenter/step3/step3_singlephoton_Delta1_nopu.root")
# events = Events ("/afs/cern.ch/work/w/wredjeb/public/TICLv4/Seeding/CMSSW_12_3_0_pre1/src/hgcal_private/production/Test/singlepi_closeBy_deltaR_hgcalCenter/step3/step3_singlepi_Delta1_nopu.root")
# events = Events ("/afs/cern.ch/work/w/wredjeb/public/TICLv4/Seeding/CMSSW_12_3_0_pre1/src/hgcal_private/production/Test/new/singlepi_closeBy_deltaR_hgcalCenter/step3/step3_singlepi_Delta1_nopu.root")
# events = Events("/eos/user/w/wredjeb/HGCAL/TrackstersSeed/Samples/singlephoton_closeBy_deltaR_hgcalCenter/step3/step3_singlephoton_Delta10_nopu.root")
# events = Events ("/afs/cern.ch/work/w/wredjeb/public/TICLv4/Seeding/CMSSW_12_3_0_pre1/src/hgcal_private/production/Test/newNoFlatEta/singlepi_closeBy_deltaR_hgcalCenter/step3/step3_singlepi_Delta1_nopu.root")
# events = Events("/afs/cern.ch/work/w/wredjeb/public/TICLv4/Seeding/CMSSW_12_3_0_pre1/src/hgcal_private/production/Test/new2/singlephoton_closeBy_deltaR_hgcalCenter/step3/step3_singlephoton_Delta10_nopu.root")
events = Events("/afs/cern.ch/work/w/wredjeb/public/TICLv4/FineSimTracksters/CMSSW_12_3_0_pre4/src/34693.0_CloseByParticleGun+2026D76+CloseByParticle_Photon_ERZRanges_GenSimHLBeamSpotHGCALCloseBy+DigiTrigger+RecoGlobal+HARVESTGlobal/step3.root")

outputdir = './fineSimTracksters/'
outputdir3D = outputdir + 'projections/'
mkdir_p(outputdir)
mkdir_p(outputdir3D)

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

clue3DLowseedsLabel = ("ticlTrackstersCLUE3DLow", "tracksterSeeds")
clue3DHighseedsLabel = ("ticlTrackstersCLUE3DHigh", "tracksterSeeds")

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

for i_ev, ev in enumerate(events):
    if(i_ev % 100 == 0):
        print(f"-- Processed events {i_ev}")
    ev.getByLabel (gpL, gpH)
    ev.getByLabel (cpslabel, caloparticlesH)
    ev.getByLabel (layerClustersLabel, layerClusters)
    ev.getByLabel (simtrackstersL, simtrackstersH)
    ev.getByLabel (finesimtrackstersL, finesimtrackstersH)
    ev.getByLabel (clue3DLowseedsLabel, clue3DLowseedsH)
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
    s_doubletsLow = clue3DLowseedsDoubletsH.product()
    s_doubletsHigh = clue3DHighseedsDoubletsH.product()
    cps = caloparticlesH.product()
    if(clue3DLowTrackstersH.isValid()):
        tracksterLow = clue3DLowTrackstersH.product()
    else:
        tracksterLow = []

    if(clue3DHighTrackstersH.isValid()):
        tracksterHigh = clue3DHighTrackstersH.product()
    else:
        tracksterHigh = []
    CAtrackster = CATracksterMergeH.product()

    number_seedsLow += len(seedsLow)
    number_seedsHigh += len(seedsHigh)
    
    number_tracksterLow += len(tracksterLow)  
    number_tracksterHigh += len(tracksterHigh)

    number_CAtracksters += len(CAtrackster)
    
    seeds_merge_ev = []

    lcs_simTracksters = []
    lcs_fineSimTracksters = []
    lcsE_simTracksters = []
    lcsE_fineSimTracksters = []
    lcsTrkIdx_simTracksters = []
    lcsETrkIdx_fineSimTracksters = []
    lcs_CAtrackster = []
    lcsE_CAtrackster = []
    lcsTrkIdx_CAtrackster = []

    for i_t in range(len(tracksterLow)):
        t_l = tracksterLow[i_t]
        # print(len(t_l.vertices()))
        # print(int(t_l.vertices().size()))
        N = int(t_l.vertices().size())
        for k in range(N):
            t_l.vertices(k)
        clue3DTracksterMerge.append(t_l)
    for i_t in range(len(tracksterHigh)):
        t_h = tracksterHigh[i_t]
        # print(int(t_h.vertices().size()))
        N = int(t_h.vertices().size())
        for k in range(N):
            t_h.vertices(k)
        clue3DTracksterMerge.append(t_h)

    for gp in gps:
        gp_eta.append(gp.eta())
        gp_phi.append(gp.phi())
        gp_energy.append(gp.energy())

    for cp in cps:
        cp_eta.append(cp.eta())
        cp_phi.append(cp.phi())
        cp_energy.append(cp.energy())
    
    for lc in lcs:
        lcs_eta.append(lc.eta())
        lcs_phi.append(lc.phi())
        lcs_energy.append(lc.energy())
        lcs_z.append(lc.z())
    final_idx = 0
    print(f"SimTracksters {len(simTracksters)}, FineSimTracksters {len(finesimTracksters)}, CARecoTrackster {len(CAtrackster)}")
    # if(len(simTracksters) != len(finesimTracksters)):
        # print(f"SimTracksters {len(simTracksters)}, FineSimTracksters {len(finesimTracksters)}")
    print(f"SimTracksters {len(simTracksters)}, FineSimTracksters {len(finesimTracksters)}")
    if(len(simTracksters) != len(finesimTracksters)):
        print(f"SimTracksters {len(simTracksters)}, FineSimTracksters {len(finesimTracksters)}")
    for idx, st in enumerate(simTracksters):
        N_lcs = st.vertices().size()
        for i_dx in range(N_lcs):
            sim_lc = lcs[st.vertices(i_dx)]
            lcs_simTracksters.append(sim_lc)
            lcsE_simTracksters.append(sim_lc.energy())
            lcsTrkIdx_simTracksters.append(idx)
        final_idx = idx
    for idx, st in enumerate(finesimTracksters):
        N_lcs = st.vertices().size()
        for i_dx in range(N_lcs):
            sim_lc = lcs[st.vertices(i_dx)]
            lcs_fineSimTracksters.append(sim_lc)
            lcsE_fineSimTracksters.append(sim_lc.energy())
            # print(idx + final_idx)
            lcsETrkIdx_fineSimTracksters.append(idx*2)

    for fs in finesimTracksters:
        for i in range(fs.vertices().size()):
            found = False
            lc_fs = lcs[fs.vertices(i)]
            for ss in simTracksters:
                for j in range(ss.vertices().size()):
                    lc_ss = lcs[ss.vertices(j)]
                    if(matchlc(lc_fs, lc_ss) == True):
                        # print(f"lc_ss {lc_fs.x(), lc_fs.y(), lc_fs.z()} lc_fs {lc_ss.x(), lc_ss.y(), lc_ss.z()}")
                        found = True
            if(found == False):
                print(f"lc_ss {lc_fs.x(), lc_fs.y(), lc_fs.z()}")# lc_fs {lc_ss.x(), lc_ss.y(), lc_ss.z()}")
                print("FOUND NON MATCHING LCS")
        final_idx = idx
    for idx, st in enumerate(CAtrackster):
        N_lcs = st.vertices().size()
        for i_dx in range(N_lcs):
            sim_lc = lcs[st.vertices(i_dx)]
            lcs_CAtrackster.append(sim_lc)
            lcsE_CAtrackster.append(sim_lc.energy())
            lcsTrkIdx_CAtrackster.append(idx*2)
        
    print(len(CAtrackster), lcsTrkIdx_CAtrackster )
    final_sim_lc = lcs_simTracksters
    final_sim_lcE = lcsE_simTracksters
    final_sim_lcIdx = lcsTrkIdx_simTracksters
    # for idx, fs in enumerate(lcs_fineSimTracksters):
    #     final_sim_lc.append(fs)
    #     final_sim_lcE.append(fs.energy())
    #     final_sim_lcIdx.append(lcsETrkIdx_fineSimTracksters[idx])
    # final_sim_lc = [lcs_simTracksters.append(i) for i in lcs_fineSimTracksters]
    # final_sim_lcE = [lcsE_simTracksters.append(i) for i in lcsE_fineSimTracksters]
    # final_sim_lcIdx = [lcsTrkIdx_simTracksters.append(i) for i in lcsETrkIdx_fineSimTracksters]
    # print(final_sim_lc[0])
    # if(i_ev % 1 == 0):
        # plots3DwithProjectionFineTracksters3D(lcs_simTracksters, lcsE_simTracksters, lcsTrkIdx_simTracksters,lcs_fineSimTracksters, lcsE_fineSimTracksters, lcsETrkIdx_fineSimTracksters,[], outputdir3D + 'finetrackster_'+str(i_dx)+"_"+str(i_ev) + '.png')
    for i_sL in seedsLow:
        sL = lcs[i_sL]
        seeds_eta.append(sL.eta())
        seeds_phi.append(sL.phi())
        seeds_energy.append(sL.energy())
        seeds_z.append(sL.z())
        seeds_merge.append(i_sL)
        seeds_merge_ev.append(sL)

    for i_sH in seedsHigh:
        sH = lcs[i_sH]
        seeds_eta.append(sH.eta())
        seeds_phi.append(sH.phi())
        seeds_energy.append(sH.energy())
        seeds_z.append(sH.z())
        seeds_merge.append(i_sH)
        seeds_merge_ev.append(sH)

    lcs_to_plot = []
    lcsE_to_plot = []

    c3dlcs_to_plot = []
    c3dlcsE_to_plot = []
    seeds_to_plot = []    
    index_trackster = []
    index_tracksterc3d = []
    if(i_ev % 1 == 0):
        for i_t,t in enumerate(CAtrackster):
            n = len(t.vertices())
            lcs_trackster = []
            lcsE_trackster = []
            for j in range(n):
                lcs_to_plot.append(lcs[t.vertices(j)])
                lcsE_to_plot.append(lcs[t.vertices(j)].energy())
                index_trackster.append(i_t)

        for i_t in range(len(tracksterLow)):
            t_l = tracksterLow[i_t]
            N = len(t_l.vertices())
            for k in range(N):
                c3dlcs_to_plot.append(lcs[t_l.vertices(k)])
                c3dlcsE_to_plot.append(lcs[t_l.vertices(k)].energy())
                index_tracksterc3d.append(i_t)

        for i_tH in range(len(tracksterHigh)):
            t_h = tracksterHigh[i_tH]
            N = len(t_h.vertices())
            for z in range(N):
                c3dlcs_to_plot.append(lcs[t_h.vertices(z)])
                c3dlcsE_to_plot.append(lcs[t_h.vertices(z)].energy())
                index_tracksterc3d.append(i_tH)

        for s in seeds_merge_ev:
            seeds_to_plot.append(s)
        # for j in range(len(lcs_to_plot)):
        # plots3DwithProjectionSeeds(lcs_to_plot, lcsE_to_plot, index_trackster, seeds_to_plot, outputdir3D + 'trackster_'+str(j)+"_"+str(i_ev) + '.png')
        lc_x = []
        lc_y = []
        lc_z = []
        for l in lcs_simTracksters:
            lc_x.append(l.x())
            lc_y.append(l.y())
            lc_z.append(l.z())
        for l in lcs_fineSimTracksters:
            lc_x.append(l.x())
            lc_y.append(l.y())
            lc_z.append(l.z())
        minx = min(lc_x) -5
        maxx = max(lc_x) + 5
        miny = min(lc_y) - 5
        maxy = max(lc_y) + 5
        minz = min(lc_z) -5 
        maxz = max(lc_z) + 5
        plots3DwithProjectionSeeds3D(lcs_simTracksters, lcsE_simTracksters, lcsTrkIdx_simTracksters,lcs_CAtrackster, lcsE_CAtrackster, lcsTrkIdx_CAtrackster,[], minx,maxx, miny,maxy,minz,maxz, "SimTracksters", outputdir3D + 'Simtrackster_'+str(j)+"_"+str(i_ev) + '.png', map_ = False)
        plots3DwithProjectionSeeds3D(lcs_fineSimTracksters, lcsE_fineSimTracksters, lcsETrkIdx_fineSimTracksters,lcs_CAtrackster, lcsE_CAtrackster, lcsTrkIdx_CAtrackster,[], minx,maxx, miny,maxy,minz,maxz, "SlimSimTracksters", outputdir3D + 'fineSimtrackster_'+str(j)+"_"+str(i_ev) + '.png', map_ = True)


bins = 30
plotHisto(seeds_eta, bins,"seeds eta", "eta", outputdir + "seeds_eta")
plotHisto(seeds_phi, bins, "seeds phi", "phi", outputdir + "seeds_phi")
plotHisto(seeds_energy, 40, "seeds energy", "E", outputdir + "seeds_energy")
plotHisto(seeds_z, 30, "seeds z", "z", outputdir + "seeds_z")

plotHisto(lcs_eta, bins, "lcs eta", "eta", outputdir + "lcs_eta")
plotHisto(lcs_phi, bins, "lcs phi","phi", outputdir + "lcs_phi")
plotHisto(lcs_energy, 40, "lcs energy", "energy", outputdir + "lcs_energy")
plotHisto(lcs_z, 30, "lcs z", "z", outputdir + "lcs_z")

plotHisto(cp_eta, bins, "cp eta", "eta", outputdir + "cp_eta")
plotHisto(cp_phi, bins, "cp phi", "phi", outputdir + "cp_phi")
plotHisto(cp_energy, 30, "cp energy", "energy", outputdir + "cp_energy")


plotHisto(gp_eta, bins, "genPart eta", "eta", outputdir + "gp_eta")
plotHisto(gp_phi, bins, "genPart phi", "phi", outputdir + "gp_phi")
plotHisto(gp_energy, 30, "genPart energy", "energy", outputdir + "gp_energy")




