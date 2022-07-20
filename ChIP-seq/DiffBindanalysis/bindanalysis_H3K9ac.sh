# filter bed files (required) to use them in unionbedg - get the files into the right format
# REQUIRED PREPARATION: add files of one histone mark (e.g. H3K9ac) of every patient into one folder
# we only used the patients that have data for each histone mark measurement
for i in ./*.bed; do awk '{print $1 "\t" $2 "\t" $3 "\t" $5}' > ../bed_files_filtered/$i; done

# just created new folder such that there is a good separation between the files
cd ../bed_files_filtered

# use the filtered files to create matrices for every individual histone mark (H3K9ac only as example)
bedtools unionbedg -i ./*.bed > ../H3K9ac_mat.txt
