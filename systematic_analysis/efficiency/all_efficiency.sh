#!/bin/bash
#loop over sub_volumes
for sub_volume in {20..32}
# for sub_volume in {00..23}
# for sub_volume in 28
do
    #loop over sub_areas
    # for reco in {1..63}
    for reco in {1..9} 10 18 19 27 28 36 37 45 46 54 {55..63}
    # for reco in 1 10 19 28 37 46 55 2 11 20 29 38 47 56
    # for reco in {3..9} {12..18} {21..27} {30..36} {39..45} {48..54} {57..63}
    do
        reco02d=$(printf "%02d" ${reco})
        ./efficiency_plate_pos /mnt/Disk5/F222/zone3/pl${sub_volume}1*/reco${reco02d}*/linked_tracks.root 30plates_effAreaEachPlate_binWidth1000_noEdgeCut ${reco02d} 1000 &
        # ./efficiency_plate_pos /mnt/Disk4/F222/zone3/pl${sub_volume}8*/reco${reco02d}*/linked_tracks.root 15plates_effAreaEachPlate_binWidth1000 ${reco02d} 1000 &
        if test 0 -eq $((reco%4)) ; then
            wait
        fi
    done
done

#loop over sub_volumes
# for sub_volume in {20..32}
for sub_volume in {00..19}
# for sub_volume in 28
do
    #loop over sub_areas
    # for reco in {1..63}
    for reco in {1..9} 10 18 19 27 28 36 37 45 46 54 {55..63}
    # for reco in 1 10 19 28 37 46 55 2 11 20 29 38 47 56
    # for reco in {3..9} {12..18} {21..27} {30..36} {39..45} {48..54} {57..63}
    do
        reco02d=$(printf "%02d" ${reco})
        # ./efficiency_plate_pos /mnt/Disk5/F222/zone3/pl${sub_volume}1*/reco${reco02d}*/linked_tracks.root 30plates_effAreaEachPlate_binWidth1000 ${reco02d} 1000 &
        ./efficiency_plate_pos /mnt/Disk4/F222/zone3/pl${sub_volume}8*/reco${reco02d}*/linked_tracks.root 15plates_effAreaEachPlate_binWidth1000_noEdgeCut ${reco02d} 1000 &
        if test 0 -eq $((reco%4)) ; then
            wait
        fi
    done
done
echo done