#!/bin/bash

dir="/global/homes/r/reynier"
outdir="/project/projectdirs/alice/reynier"
if [[ ! -e $outdir/out_AllSi ]]; then
    mkdir $outdir/out_AllSi
fi

source $dir/Singularity/cvmfs/sphenix.sdcc.bnl.gov/x8664_sl7/opt/sphenix/core/bin/sphenix_setup.sh -n
export MYINSTALL=$dir/Singularity/install
source $dir/Singularity/cvmfs/sphenix.sdcc.bnl.gov/x8664_sl7/opt/sphenix/core/bin/setup_local.sh $MYINSTALL

cd $dir/Singularity/g4lblvtx/macros/
root -b -q "Fun4All_G4_FastMom_GEM_RICH.C($2, \"$outdir/out_AllSi/out_$3_det$4_10umGEM_RICH_$5x$5_$1\",\"$3\",$4,$5)"

rm -rf $outdir/out_AllSi/out_$3_det$4_10umGEM_RICH_$5x$5_$1_G4LBLVtx.root
