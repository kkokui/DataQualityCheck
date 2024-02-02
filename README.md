# how to make summary plot
Compile.
```
make quality_check
```
Execute.
```
./quality_check linked_tracks.root title Xcenter Ycenter binWidth
```
Arguments
1. linked_track.root: input file name.
1. title: title of summary plot
1. Xcenter: center value of x in the sub-area
1. Ycenter: center value of y in the sub-area
1. binWidth: division width (unit is $\mathrm{\mu m}$) of divide_align (not used if divide_align is not applied. if so, please put some value like 2000.)
