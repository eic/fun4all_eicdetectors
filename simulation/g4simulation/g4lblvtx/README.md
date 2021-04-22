This code corresponds to an implementation of the tapered All-Silicon Tracker design described [here](https://indico.bnl.gov/event/7892/contributions/36938/attachments/27856/42740/20200430-EICUG_Tracking_WG_-_eRD16.pdf) into the Fun4All framework.

# Getting Started
1. The code was originally ran inside a Singularity container, which can be found here:
https://github.com/sPHENIX-Collaboration/Singularity
One needs to run ./updatebuild.sh and follow the steps in the README of that repo.

2. The code needs to built and installed with Fun4All, as outlined here: https://wiki.bnl.gov/sPHENIX/index.php/Example_of_using_DST_nodes. References to <sourcedir> for this repository mean `g4lblvtx/src/`, where one should see the autogen.sh file. Create a `build` and `install` directory (I suggest in the same directory that holds this repo) and follow the instructions under the "Building a package" section from the link.

The main steps are:

```
mkdir build
cd build
/global/path/to/g4lblvtx/source/autogen.sh --prefix=$MYINSTALL
make -j 4
make install
source /global/path/to/Singularity/cvmfs/sphenix.sdcc.bnl.gov/x8664_sl7/opt/sphenix/core/bin/sphenix_setup.sh -n
export MYINSTALL=/global/path/to/Singularity/install
source /global/path/to/Singularity/cvmfs/sphenix.sdcc.bnl.gov/x8664_sl7/opt/sphenix/core/bin/setup_local.sh $MYINSTALL
```

# Running a batch job on Cori
From outside the container:
1. `cd cori_batch`.
2. Edit the `run_shared.sh` script (see instructions below).
3. Run `sbatch run_shared.sh`

### Editing the run_shared.sh script

Edit the following lines:
1. `#SBATCH --time=1:00:00` -> Approximate time that running a job will take. A higher value ensures that your job is not stopped because it is taking too long. A lower value ensures that your job starts running sooner.
2. `#SBATCH --array=0-999` -> This means the output will be split into 1000 output root files.
3. `shifter ./AllSi_shifter.sh $SLURM_ARRAY_TASK_ID 100` -> The important parameter to keep in mind here is the number at the end of the line. This is the number of events per root split. 

### Checking the status of your job:

1. On your internet browser, go to: [https://my.nersc.gov/](https://my.nersc.gov/).
2. Login using your Nersc credentials.
3. Go to the tab `Cori Queues`.
