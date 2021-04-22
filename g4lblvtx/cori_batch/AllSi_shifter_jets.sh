#!/bin/bash

#dir="/global/project/projectdirs/alice/ftorales" #replace per user, often with ~
#dir="~"
dir="/global/homes/r/reynier"
outdir="/project/projectdirs/alice/reynier"
if [[ ! -e $outdir/out_AllSi ]]; then
    mkdir $outdir/out_AllSi
fi

source $dir/Singularity/cvmfs/sphenix.sdcc.bnl.gov/x8664_sl7/opt/sphenix/core/bin/sphenix_setup.sh -n
export MYINSTALL=$dir/Singularity/install
source $dir/Singularity/cvmfs/sphenix.sdcc.bnl.gov/x8664_sl7/opt/sphenix/core/bin/setup_local.sh $MYINSTALL

cd $dir/Singularity/g4lblvtx/macros/
root -b -q "Fun4All_G4_FastMom.C($2, \"$outdir/out_AllSi/out_jets_$1\")"

rm -rf $outdir/out_AllSi/out_jets_$1_G4LBLVtx.root
rm -rf $outdir/out_AllSi/out_jets_$1_g4jet_eval.root
rm -rf $outdir/out_AllSi/out_jets_$1_Beast_FastTrackingEval.root
# rm -rf $outdir/out_AllSi/out_jets_$1_electrons+jets.root
