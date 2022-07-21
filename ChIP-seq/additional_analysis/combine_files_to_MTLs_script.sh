# concatenate the histone mark MIT files from the different patients into one MTL file
cat MIT.Y_H3K9ac.bed MIT.O_H3K9ac.bed MIT.AD_H3K9ac.bed > MTL.H3K9ac.bed
cat MIT.Y_H3K27ac.bed MIT.O_H3K27ac.bed MIT.AD_H3K27ac.bed > MTL.H3K27ac.bed
cat MIT.Y_H3K122ac.bed MIT.O_H3K122ac.bed MIT.AD_H3K122ac.bed > MTL.H3K122ac.bed

# next step would be to sort and merge the MTL files!
# sort:
for i in ./*.bed; do sort -V -k1,1 -k2,2n -i $i -o $i_sorted.bed; done

# merge:
for i in ./*_sorted.bed; do bedtools2/bin/mergeBed -i $i -c 5 -o mean > $i_merged.bed

