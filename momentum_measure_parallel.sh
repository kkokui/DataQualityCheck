#!/bin/bash
data_path="/data/FASER/F222/zone4/temp/TFD/vert32063_pl053_167_new/reco32_065000_050000/v15/linked_tracks.root"
measure_momentum() {
    ./measure_momentum $1 before_align_$2 "npl>=100&&Entry$%5=="$2
}
export -f measure_momentum
seq 0 4 | parallel -j 5 -u measure_momentum ${data_path}