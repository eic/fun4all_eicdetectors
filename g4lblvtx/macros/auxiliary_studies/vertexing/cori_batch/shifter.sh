#!/bin/bash

dir="/global/homes/r/reynier"
outdir="/project/projectdirs/alice/reynier"
if [[ ! -e $outdir/out_AllSi ]]; then
    mkdir $outdir/out_AllSi
fi

source $dir/Singularity/cvmfs/sphenix.sdcc.bnl.gov/x8664_sl7/opt/sphenix/core/bin/sphenix_setup.sh -n
export MYINSTALL=$dir/Singularity/install
source $dir/Singularity/cvmfs/sphenix.sdcc.bnl.gov/x8664_sl7/opt/sphenix/core/bin/setup_local.sh $MYINSTALL

cd $dir/Singularity/g4lblvtx/macros/auxiliary_studies/vertexing/
root -b -q "Fun4All_G4_simple_vertex.C($2,$3,$4,$5,$6,$7,\"$outdir/out_AllSi/out_vtx_$1_$3$4$5_$6XX0_$7um\")"

rm -rf $outdir/out_AllSi/out_vtx_$1_$3$4$5_$6XX0_$7um_G4LBLVtx.root
