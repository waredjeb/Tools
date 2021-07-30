path='/data/wredjeb/SingleGammaFlatPt8To150/'
mkdir -p $path
echo 'mkdir -p' $path
for i in {1..10}
do
    subPath=$path"Job"$i"/"
    echo "mkdir -p" $subPath
    mkdir -p $subPath
    outputFileStep1=$subPath'step1_'$i".root"
    echo 'cmsRun SingleGammaFlatPt8To150_cfi_GEN_SIM.py' $outputFileStep1 $outputFileStep1 '-n 8'
    cmsRun SingleGammaFlatPt8To150_cfi_GEN_SIM.py $outputFileStep1 -n 8
    outputFileStep2=$subPath'step2_'$i".root"
    echo 'cmsRun step2_DIGI_L1TrackTrigger_L1_DIGI2RAW_HLT.py' $outputFileStep1 $outputFileStep2 '-n 8'
    cmsRun step2_DIGI_L1TrackTrigger_L1_DIGI2RAW_HLT.py $outputFileStep1 $outputFileStep2 -n 8
    outputFileStep3=$subPath'step3_'$i".root"
    echo 'cmsRun step3_RAW2DIGI_L1Reco_RECO_RECOSIM_PAT_VALIDATION_DQM.py' $outputFileStep2 $outputFileStep3 '-n 8'
    cmsRun step3_RAW2DIGI_L1Reco_RECO_RECOSIM_PAT_VALIDATION_DQM.py $outputFileStep2 $outputFileStep3 -n 8
done

