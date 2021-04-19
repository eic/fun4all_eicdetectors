#!/bin/bash

dir="/global/homes/r/reynier"
outdir="/project/projectdirs/alice/reynier"
if [[ ! -e $outdir/out_sPHENIX ]]; then
    mkdir $outdir/out_sPHENIX
fi

source $dir/Singularity/cvmfs/sphenix.sdcc.bnl.gov/x8664_sl7/opt/sphenix/core/bin/sphenix_setup.sh -n
export MYINSTALL=$dir/Singularity/install
source $dir/Singularity/cvmfs/sphenix.sdcc.bnl.gov/x8664_sl7/opt/sphenix/core/bin/setup_local.sh $MYINSTALL

cd $dir/Singularity/macros/macros/g4simulations/
root -b -q "Fun4All_G4_sPHENIX_particle_gen.C($2, \"$outdir/out_sPHENIX/out_$3_$1\",\"$3\")"
