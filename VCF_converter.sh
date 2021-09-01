#!/usr/bin/env bash

#SBATCH --job-name="VCF_converter"
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -o out/VCF_conveter-%N-%j.out
#SBATCH -e err/VCF_converter-%N-%j.err
#SBATCH --mail-user=megahitokiri@hotmail.com
#SBATCH --mail-type=END
#SBATCH --account=pezzolesi-np
#SBATCH --partition=pezzolesi-np

#VCF='/uufs/chpc.utah.edu/common/home/u1123911/EthSEQ/VCF/Low_Filter_retinosis1_RSID.vcf.gz'
#VCF='/uufs/chpc.utah.edu/common/home/u1123911/EthSEQ/VCF/UKS_Final_206_samples.nonRelated.vcf.rsid.gz'
VCF='/uufs/chpc.utah.edu/common/home/u1123911/EthSEQ/VCF/Gencove_imputed_merged.vcf.gz'

tabix -p vcf $VCF

ml bcftools

bcftools view -Oz --max-alleles 2 --exclude-types indels $VCF -o Intermediate.vcf.gz 

tabix -p vcf Intermediate.vcf.gz

bcftools query -H -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t\QUAL\t\FILTER\t\INFO\t\GT[\t%GT]\n' Intermediate.vcf.gz > Intermediate2.vcf

sed '1s/\]/_/g' Intermediate2.vcf | sed '1s/\[/C/g' | sed '1s/\:GT//g' > RP1.vcf

head -n 1 RP1.vcf > RP1.rs.vcf
grep rs RP1.vcf >> RP1.rs.vcf

#mv RP1.rs.vcf Targeted_Seq.rs.vcf

rm Intermediate.vcf.gz
rm Intermediate.vcf.gz.tbi
rm Intermediate2.vcf

ml unload bcftools
