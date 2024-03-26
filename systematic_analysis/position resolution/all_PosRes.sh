#!/bin/bash
#loop over sub_volumes
for sub_volume in {20..32}
# for sub_volume in 00
do
    #loop over sub_areas
    for reco in {1..63}
    do
       reco02d=$(printf "%02d" ${reco})
    #    ./calc_dxy2 /mnt/Disk5/F222/zone3/pl${sub_volume}1*/reco${reco02d}*/linked_tracks.root pl${sub_volume}1_reco${reco02d} ${reco02d} --NoAlign&
       if test 0 -eq $((reco%4)) ; then
           wait
       fi
    done
    # for reco in {1..63}
    for reco in 1 10 19 28 37 46 55
    do
        reco02d=$(printf "%02d" ${reco})
        ./fit_deltaXY /mnt/Disk5/F222/zone3/deltaXY/tree_pl${sub_volume}1_reco${reco02d}_BeforeAlign.root 30plates_zone3_division1_xMin11000 ${sub_volume} ${reco} &
        if test 0 -eq $((reco%9)) ; then
            wait
        fi
    done
done
#loop over sub_volumes
for sub_volume in {00..19}
# for sub_volume in 00
do
    # for reco in {1..63}
    for reco in 1 10 19 28 37 46 55
    do
        reco02d=$(printf "%02d" ${reco})
        ./fit_deltaXY /home/administrator/kokui/sub_volume_PosRes/deltaXY/tree_pl${sub_volume}8_reco${reco02d}_BeforeAlign.root 15plates_zone3_division1_xMin11000 ${sub_volume} ${reco} &
        if test 0 -eq $((reco%9)) ; then
            wait
        fi
    done
done
# cat PosRes_zone3_plate_Nodivisions/reco* >>PosRes_zone3_plate_Nodivisions/all.txt