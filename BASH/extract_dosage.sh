#!/bin/bash
#SBATCH --job-name=get_dosages_AL
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=get_dosages.out
#SBATCH --error=get_dosages.err
#SBATCH --mem=20gb

module load plink/1.9-foss-2015b
module load R/3.3.3-foss-2015b
module load RPlus
module load Perl/5.26.2-foss-2015b-bare
module load parallel
module load VCFtools
module load tabix
module load BCFtools/1.9-foss-2018b

inputData='/groups/umcg-griac/tmp01/projects/umcg-cqi/PIAMA/genotype/data/dosage/annot'
outputData='/groups/umcg-griac/tmp01/projects/umcg-alan/forZaid/piama/dosages/LD'
list_dir='/groups/umcg-griac/tmp01/projects/umcg-alan/forZaid/piama/'

cd $inputData
for chrN in $(awk '{print $1}' $list_dir/SNPlistNEW46.txt | sort | uniq) # Only loop over relevant chromosomes
do
echo "Now started on chromosome = " $chrN
vcftools --gzvcf chr$chrN'_dbSNP151.vcf.gz' --positions $list_dir/SNPlistNEW46.txt --out $outputData/chr$chrN'_SNPfilter_asthma' --recode --keep-INFO-all
# If

bgzip $outputData/chr$chrN'_SNPfilter_asthma.recode.vcf'
tabix $outputData/chr$chrN'_SNPfilter_asthma.recode.vcf.gz'
echo $outputData/chr$chrN'_SNPfilter_asthma.recode.vcf.gz' >> $outputData/vcf_merge_list.txt
done

cd $outputData
pwd

# Merge datasets per chromosome which filenames are in vcf_merge_list.txt
bcftools concat --file-list vcf_merge_list.txt -o merged_dosages.vcf.gz

# Recode gzvcf file to standard vcf file
vcftools --gzvcf merged_dosages.vcf.gz --recode --out merged_dosages.vcf

rm -rf merged_dosages.txt || true
bcftools query -f "%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%DS]\n" merged_dosages.vcf.gz --out merged_dosages_scores.txt

grep -E '^#CHROM' merged_dosages.vcf.recode.vcf > headers_temp.txt

sed -i 's/QUAL.*FORMAT/empty/g' headers_temp.txt
sed -i 's/#//g' headers_temp.txt
awk '{$6=""; print $0}' headers_temp.txt > headers.txt

cat headers.txt merged_dosages_scores.txt > merged_dosagesLD.txt
