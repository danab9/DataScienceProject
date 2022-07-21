# concatenate the different patient files into one file
cat <patient_file1.bed> <patient_file2.bed> ... <patient_fileX.bed> > MIT.Y_H3K9ac.bed
# do this for all 12 cases

# sort all concatenated files (since they were only concatenated) - required for merging
sort -V -k1,1 -k2,2n -i MIT.AD_H3K9ac.bed -o MIT.AD_H3K9ac_new.bed
sort -V -k1,1 -k2,2n -i MIT.Y_H3K9ac.bed -o MIT.Y_H3K9ac_new.bed
sort -V -k1,1 -k2,2n -i MIT.O_H3K9ac.bed -o MIT.O_H3K9ac_new.bed

sort -V -k1,1 -k2,2n -i MIT.O_H3K27ac.bed -o MIT.O_H3K27ac_new.bed
sort -V -k1,1 -k2,2n -i MIT.Y_H3K27ac.bed -o MIT.Y_H3K27ac_new.bed
sort -V -k1,1 -k2,2n -i MIT.AD_H3K27ac.bed -o MIT.AD_H3K27ac_new.bed

sort -V -k1,1 -k2,2n -i MIT.AD_H3K122ac.bed -o MIT.AD_H3K122ac_new.bed
sort -V -k1,1 -k2,2n -i MIT.Y_H3K122ac.bed -o MIT.Y_H3K122ac_new.bed
sort -V -k1,1 -k2,2n -i MIT.O_H3K122ac.bed -o MIT.O_H3K122ac_new.bed

sort -V -k1,1 -k2,2n -i MIT.O_H3K4me1.bed -o MIT.O_H3K4me1_new.bed
sort -V -k1,1 -k2,2n -i MIT.Y_H3K4me1.bed -o MIT.Y_H3K4me1_new.bed
sort -V -k1,1 -k2,2n -i MIT.AD_H3K4me1.bed -o MIT.AD_H3K4me1_new.bed


# merge all overlapping regions from each file, can also be done using bash script
bedtools2/bin/mergeBed -i MIT.AD_H3K9ac_new.bed -c 5 -o mean > MIT.AD_H3K9ac_merged.bed
bedtools2/bin/mergeBed -i MIT.Y_H3K9ac_new.bed -c 5 -o mean > MIT.Y_H3K9ac_merged.bed
bedtools2/bin/mergeBed -i MIT.O_H3K9ac_new.bed -c 5 -o mean > MIT.O_H3K9ac_merged.bed

bedtools2/bin/mergeBed -i MIT.O_H3K27ac_new.bed -c 5 -o mean > MIT.O_H3K27ac_merged.bed
bedtools2/bin/mergeBed -i MIT.Y_H3K27ac_new.bed -c 5 -o mean > MIT.Y_H3K27ac_merged.bed
bedtools2/bin/mergeBed -i MIT.AD_H3K27ac_new.bed -c 5 -o mean > MIT.AD_H3K27ac_merged.bed

bedtools2/bin/mergeBed -i MIT.AD_H3K122ac_new.bed -c 5 -o mean > MIT.AD_H3K122ac_merged.bed
bedtools2/bin/mergeBed -i MIT.Y_H3K122ac_new.bed -c 5 -o mean > MIT.Y_H3K122ac_merged.bed
bedtools2/bin/mergeBed -i MIT.O_H3K122ac_new.bed -c 5 -o mean > MIT.O_H3K122ac_merged.bed

bedtools2/bin/mergeBed -i MIT.O_H3K4me1_new.bed -c 5 -o mean > MIT.O_H3K4me1_merged.bed
bedtools2/bin/mergeBed -i MIT.Y_H3K4me1_new.bed -c 5 -o mean > MIT.Y_H3K4me1_merged.bed
bedtools2/bin/mergeBed -i MIT.AD_H3K4me1_new.bed -c 5 -o mean > MIT.AD_H3K4me1_merged.bed
