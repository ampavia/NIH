interact
ml bioawk/1.0-GCC-11.2.0
cd /scratch/ac05869/gese_final_yahs/gese_v2.asm/100kb
bioawk -c fastx '{print $name "\t" $seq}' gese_v2.asm.fa > gese_v2.tsv
bioawk -c fastx '{print $name "\t" length($seq) "\t" $seq}' gese_v2.asm.fa > gese_v2.tsv

#Rename contigs in R with ~/NIH/rename_contigs.R

sed 's/\t/\n/g' gese_v2_renamed.tsv > ../final_asm_renamed_contigs/gese_v2.asm.fa
cd ../final_asm_renamed_contigs/
md5sum gese_v2.asm.fa > gese_v2.asm.fa.md5