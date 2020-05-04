# Getting the first three exons

grep "exon1;" PCITRI_v0_all_annotation1.txt > exon1.txt
wc -l exon1.txt
# 41192

grep "exon2;" PCITRI_v0_all_annotation1.txt > exon2.txt
wc -l exon2.txt 
# 25122

grep "exon3;" PCITRI_v0_all_annotation1.txt > exon3.txt
wc -l exon3.txt 
# 18039

cat exon1.txt exon2.txt exon3.txt > first_3_exons_annotation.txt