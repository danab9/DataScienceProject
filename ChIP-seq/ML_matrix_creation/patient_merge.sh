# merge the different peaks for the individual patients
cat ../../bed_files/Y_H3K9ac/GSM3752815.bed ../../bed_files/Y_H3K27ac/GSM3752845.bed ../../bed_files/Y_H3K122ac/GSM3752875.bed ../../bed_files/Y_H3K4me1/GSM3752902.bed > Y_p1_unsorted.bed


# sort the files of the patients using bedtools
for i in *.bed; do sort -V -k1,1 -k2,2n -i $i -o ../patient_files_sorted/$i; done

# merge the files bedtools (adjust path!)
for i in ./*.bed; do ~/bedtools2/bin/mergeBed -i $i -c 5 -o mean > ../patient_files_merged/$i; done
 

# take the ln of coverage for each file (adjust path!)
for i in ../patient_files_merged/*.bed; do awk -v column=4 ' { $column = log($column); print }' $i > $i; done

# merge all patient files into one matrix -> input for R script filtering the huge matrix
bedtools unionbedg -i ./*.bed > ../healthy_matrix.bed
