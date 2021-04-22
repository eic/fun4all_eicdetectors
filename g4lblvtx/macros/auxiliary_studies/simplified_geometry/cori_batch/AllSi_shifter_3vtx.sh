#!/bin/bash

dir="/global/homes/r/reynier"
outdir="/project/projectdirs/alice/reynier"
if [[ ! -e $outdir/out_AllSi/Simplified_geom/$1 ]]; then
    mkdir $outdir/out_AllSi/Simplified_geom/$1
fi

source $dir/Singularity/cvmfs/sphenix.sdcc.bnl.gov/x8664_sl7/opt/sphenix/core/bin/sphenix_setup.sh -n
export MYINSTALL=$dir/Singularity/install
source $dir/Singularity/cvmfs/sphenix.sdcc.bnl.gov/x8664_sl7/opt/sphenix/core/bin/setup_local.sh $MYINSTALL

cd $dir/Singularity/g4lblvtx/macros/auxiliary_studies/simplified_geometry/
root -b -q "Fun4All_G4_simplified_3vtx.C($3,$4,$5,$6,$7,$8,$9,\"$outdir/out_AllSi/Simplified_geom/$1/out_simp_geom_3vtx_vbd_$4_$5_$6_split_$2\")"

if [ $9 == 1 ]
then
        BLAB="_B_1.5T"
elif [ $9 == 2 ]
then
        BLAB="_B_3.0T"
elif [ $9 == 3 ]
then
        BLAB="_sPHENIX"
else
        BLAB="_Beast"
fi

rm -rf $outdir/out_AllSi/Simplified_geom/$1/out_simp_geom_3vtx_vbd_$4_$5_$6_split_$2$BLAB_G4LBLVtx.root
