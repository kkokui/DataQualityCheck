#!/bin/bash
divide_align() {
    data="/data/FASER/F222/zone4/temp/TFD/vert32063_pl053_167_new/reco32_065000_050000/v15/linked_tracks.root"
    ./divide_align ${data} binWidth${1}_robustFactor${2} 32 ${1} ${2}
}
export -f divide_align
# parallel -j 5 -u divide_align ::: 5000 2000 1000 500 ::: 1.0 0.{6..9}

measure_momentum() {
    # data="/data/FASER/F222/zone4/temp/TFD/vert32063_pl053_167_new/reco32_065000_050000/v15/linked_tracks.root"
    # ./measure_momentum ${data} before_align_${1} "npl>=100&&Entry$%5=="${1}
    data_dir="/data/Users/kokui/FASERnu/F222/zone4/temp/TFD/vert32063_pl053_167_new/reco32_065000_050000/v15"
    ./measure_momentum ${data_dir}/linked_tracks_after_align_binWidth${1}_robustFactor${2}.root after_align_binWidth${1}_robustFactor${2}_${3} "npl>=100&&Entry$%5=="${3}
}
export -f measure_momentum
# parallel -j 5 -u measure_momentum ::: {0..4}
parallel -j 5 -u measure_momentum ::: 5000 1000 500 ::: 1.0 0.{6..9} ::: {0..4}

calc_pos_res() {
    data_dir="/data/Users/kokui/FASERnu/F222/zone4/temp/TFD/vert32063_pl053_167_new/reco32_065000_050000/v15"
    ./quality_check ${data_dir}/linked_tracks_after_align_binWidth${1}_robustFactor${2}.root after_align_binWidth${1}_robustFactor${2} 65000 50000 ${1}
}
export -f calc_pos_res
parallel -j 5 -u calc_pos_res ::: 2000 5000 1000 500 ::: 1.0 0.{6..9}