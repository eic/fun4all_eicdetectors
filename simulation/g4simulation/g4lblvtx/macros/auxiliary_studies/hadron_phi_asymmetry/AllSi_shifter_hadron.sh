#!/bin/bash

dir="/global/homes/r/reynier"
outdir="/project/projectdirs/alice/reynier"
if [[ ! -e $outdir/out_AKiselev ]]; then
    mkdir $outdir/out_AKiselev
fi

source $dir/Singularity/cvmfs/sphenix.sdcc.bnl.gov/x8664_sl7/opt/sphenix/core/bin/sphenix_setup.sh -n
export MYINSTALL=$dir/Singularity/install
source $dir/Singularity/cvmfs/sphenix.sdcc.bnl.gov/x8664_sl7/opt/sphenix/core/bin/setup_local.sh $MYINSTALL

cd $dir/Singularity/g4lblvtx/macros/A_Kiselev_study/Momentum_tutorial/
root -b -q "Fun4All_G4_simple_hadron_GEM.C($2,$3,\"$outdir/out_AKiselev/out_hadron_$1\")"

val2=$((45+$3))

rm -rf $outdir"/out_AKiselev/out_hadron_"$1"_"$3"_"$val2"_G4LBLVtx.root"
