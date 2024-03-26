#!/bin/bash
for i in {1..63}; do
	# for j in {0..19} ; do
	for j in 25 ; do
	# for j in {20..32} ; do
	# sub_v=$(printf "pl%03d_%03d" $(($j*10+8)) $(($j*10+22)))
	# sub_v=$(printf "pl%03d_%03d" $(($j*10+1)) $(($j*10+30)))
	sub_v=$(printf "pl%03d_%03d_ToReprocess" $(($j*10+1)) $(($j*10+30)))
	echo $sub_v
	
	  id=$(printf "%02d" $i)
	  # echo "$id"
	
	  # dir=$(find $shell_path -type d -name "reco$id*")
	  
	#   data_path=$(echo /mnt/Disk4/F222/zone3/$sub_v)
	  data_path=$(echo /mnt/Disk5/F222/zone3/$sub_v)
	  dir=$(find $data_path -type d -name "reco$id*")
	  my_dir=$(echo $dir | cut -d '/' -f7)
	#   root -l -b -q 'getEntriesWithCut.C("'$sub_v'/'$my_dir'/linked_tracks.root")'
	  ./getEntriesWithCut $sub_v/$my_dir/linked_tracks.root
	  # echo $sub_v $my_dir
	done
done
