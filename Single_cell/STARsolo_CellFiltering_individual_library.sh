#interact --mem 20G --cpus-per-task 12 --ntasks 1 --time 02:00:00

ml STAR/2.7.10b-GCC-11.3.0

STAR --runMode soloCellFiltering /scratch/ac05869/nih/kratom/star_results/MIT_AI/MIT_AI_Solo.out/GeneFull/raw/data_for_R \
/scratch/ac05869/nih/kratom/star_results/MIT_AI/MIT_AI_Solo.out/GeneFull/raw/data_for_R/EM_EmptyDrops_Combined/ --soloCellFilter EmptyDrops_CR
