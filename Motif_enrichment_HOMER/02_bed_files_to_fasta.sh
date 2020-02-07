# -------------------------------------------------------------------
# Getting the fasta sequences of a given region from a bed file
# -------------------------------------------------------------------

for file in $(ls *.bed)
do
	base=$(basename ${file} ".bed")     
    bedtools getfasta -fi PCITRI.assembly.v0.fa -bed ${file} > ${base}.fa
done

# Issue with a couple of input beds, one line start < end :S
# all_exons_3plus.bed line 53366 
# unmethylated_exons_3plus.bed line 50031

sed -i '53366d' all_exons_3plus.bed
sed -i '50031d' unmethylated_exons_3plus.bed

bedtools getfasta -fi PCITRI.assembly.v0.fa -bed all_exons_3plus.bed > all_exons_3plus.fa
bedtools getfasta -fi PCITRI.assembly.v0.fa -bed unmethylated_exons_3plus.bed > unmethylated_exons_3plus.fa

