#!/bin/bash
# for dset in 1m0l 1m56 1zoy 2ks1 2wie 3hd7 3rvy;
    
#     do 
#     # scp /home/rsamant2/scratch16-jgray21/rsamant2/ensemble_docking/MPDockTest/Ensemble_docking_fa21/1m0l/docking_fa23.slurm  "$dset"/.; 
#     if [ ! -d ./$dset ]
#     then
#         mkdir -p ./$dset
#     fi

#     scp /home/rsamant2/scr16_jgray21/rsamant2/ensemble_docking/MPDockTest/Ensemble_docking_fa21/"$dset"/complex.pdb ./"$dset"/.; 
#     scp /home/rsamant2/scr16_jgray21/rsamant2/ensemble_docking/MPDockTest/Ensemble_docking_fa21/"$dset"/native.pdb ./"$dset"/.;
#     scp /home/rsamant2/scr16_jgray21/rsamant2/ensemble_docking/MPDockTest/Ensemble_docking_fa21/"$dset"/*.span ./"$dset"/.;
#     scp /home/rsamant2/scr16_jgray21/rsamant2/ensemble_docking/MPDockTest/Ensemble_docking_fa21/"$dset"/partners ./"$dset"/.;
# done

# for dset in 1ehk 1q90 3kcu 1h2s 3chx 4dkl;
    
#     do 
#     # scp /home/rsamant2/scratch16-jgray21/rsamant2/ensemble_docking/MPDockTest/Ensemble_docking_fa21/1m0l/docking_fa23.slurm  "$dset"/.; 
#     if [ ! -d ./$dset ]
#     then
#         mkdir -p ./$dset
#     fi

#     scp /home/rsamant2/scr16_jgray21/rsamant2/ensemble_docking/MPDockTest/Ensemble_docking_fa21/difficult_targets/"$dset"/complex.pdb ./"$dset"/.; 
#     scp /home/rsamant2/scr16_jgray21/rsamant2/ensemble_docking/MPDockTest/Ensemble_docking_fa21/difficult_targets/"$dset"/native.pdb ./"$dset"/.;
#     scp /home/rsamant2/scr16_jgray21/rsamant2/ensemble_docking/MPDockTest/Ensemble_docking_fa21/difficult_targets/"$dset"/*.span ./"$dset"/.;
#     scp /home/rsamant2/scr16_jgray21/rsamant2/ensemble_docking/MPDockTest/Ensemble_docking_fa21/difficult_targets/"$dset"/partners ./"$dset"/.;
# done

# for dset in 1bl8 1e12 2k9j 2qjy 3kly 3oe0;
    
#     do 
#     # scp /home/rsamant2/scratch16-jgray21/rsamant2/ensemble_docking/MPDockTest/Ensemble_docking_fa21/1m0l/docking_fa23.slurm  "$dset"/.; 
#     if [ ! -d ./$dset ]
#     then
#         mkdir -p ./$dset
#     fi

#     scp /home/rsamant2/scr16_jgray21/rsamant2/ensemble_docking/MPDockTest/Ensemble_docking_fa21/medium_targets/"$dset"/complex.pdb ./"$dset"/.; 
#     scp /home/rsamant2/scr16_jgray21/rsamant2/ensemble_docking/MPDockTest/Ensemble_docking_fa21/medium_targets/"$dset"/native.pdb ./"$dset"/.;
#     scp /home/rsamant2/scr16_jgray21/rsamant2/ensemble_docking/MPDockTest/Ensemble_docking_fa21/medium_targets/"$dset"/*.span ./"$dset"/.;
#     scp /home/rsamant2/scr16_jgray21/rsamant2/ensemble_docking/MPDockTest/Ensemble_docking_fa21/medium_targets/"$dset"/partners ./"$dset"/.;
# done

# for dset in 1q90 2nrf 2vt4 2wie ;
    
#     do 
#     # scp /home/rsamant2/scratch16-jgray21/rsamant2/ensemble_docking/MPDockTest/Ensemble_docking_fa21/1m0l/docking_fa23.slurm  "$dset"/.; 
#     if [ ! -d ./$dset ]
#     then
#         mkdir -p ./$dset
#     fi

#     scp /home/rsamant2/data_jgray21/aharmal1/MPDockTest/Benchmark2/Ensemble-fa19/"$dset"/complex.pdb ./"$dset"/.; 
#     scp /home/rsamant2/data_jgray21/aharmal1/MPDockTest/Benchmark2/Ensemble-fa19/"$dset"/native.pdb ./"$dset"/.;
#     scp /home/rsamant2/data_jgray21/aharmal1/MPDockTest/Benchmark2/Ensemble-fa19/"$dset"/*.span ./"$dset"/.;
#     scp /home/rsamant2/data_jgray21/aharmal1/MPDockTest/Benchmark2/Ensemble-fa19/"$dset"/partners ./"$dset"/.;
# done
# for dset in 2r6q 2zxe-ab 2zxe-ag 4hg6 4huq 5a63-ac 5a63-bc 5aww;
    
#     do 
#     # scp /home/rsamant2/scratch16-jgray21/rsamant2/ensemble_docking/MPDockTest/Ensemble_docking_fa21/1m0l/docking_fa23.slurm  "$dset"/.; 
#     if [ ! -d ./$dset ]
#     then
#         mkdir -p ./$dset
#     fi

#     scp /home/rsamant2/data_jgray21/aharmal1/MPDockTest/Benchmark2/Rigid/"$dset"/complex.pdb ./"$dset"/.; 
#     scp /home/rsamant2/data_jgray21/aharmal1/MPDockTest/Benchmark2/Rigid/"$dset"/native.pdb ./"$dset"/.;
#     scp /home/rsamant2/data_jgray21/aharmal1/MPDockTest/Benchmark2/Rigid/"$dset"/*.span ./"$dset"/.;
#     scp /home/rsamant2/data_jgray21/aharmal1/MPDockTest/Benchmark2/Rigid/"$dset"/partners ./"$dset"/.;
# done

for dset in 1bl8 1e12 1ehk 1h2s 1m0l 1m56 1zoy 2k9j 2ks1 2qjy 2wie 3chx 3hd7 3kcu 3kly 3oe0 3rvy 4dkl;
    
    do 
    # scp /home/rsamant2/scratch16-jgray21/rsamant2/ensemble_docking/MPDockTest/Ensemble_docking_fa21/1m0l/docking_fa23.slurm  "$dset"/.; 
    if [ ! -d ./$dset ]
    then
        mkdir -p ./$dset
    fi

    scp /home/rsamant2/data_jgray21/aharmal1/MPDockTest/Benchmark/Structures/"$dset"/bound*.pdb ./"$dset"/.; 
    scp /home/rsamant2/data_jgray21/aharmal1/MPDockTest/Benchmark/Structures/"$dset"/unbound*.pdb ./"$dset"/.;

done

for dset in 1q90 2nrf 2r6q 2vt4 2wie 2zxe-ab 2zxe-ag 4hg6 4huq 5a63-ac 5a63-bc 5aww;
    
    do 
    # scp /home/rsamant2/scratch16-jgray21/rsamant2/ensemble_docking/MPDockTest/Ensemble_docking_fa21/1m0l/docking_fa23.slurm  "$dset"/.; 
    if [ ! -d ./$dset ]
    then
        mkdir -p ./$dset
    fi

    scp /home/rsamant2/data_jgray21/aharmal1/MPDockTest/Benchmark2/Structures/"$dset"/bound*.pdb ./"$dset"/.; 
    scp /home/rsamant2/data_jgray21/aharmal1/MPDockTest/Benchmark2/Structures/"$dset"/unbound*.pdb ./"$dset"/.;
done