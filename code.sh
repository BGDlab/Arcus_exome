
########################################################################
# extract the variants with FILTER=PASS and polyphen=probably damaging #
########################################################################
conda activate bcftools

for line in $(<exome_id.txt)
do

export GCS_OAUTH_TOKEN=$(gcloud auth print-access-token)

bcftools query -i'FILTER="PASS" && INFO/CSQ ~ "probably_damaging" && (GT="0/1" || GT="0|1" || GT="1/0" || GT="1|0" || GT="1/1" || GT="1|1" || GT="1/2" || GT="1|2")' -f'%CHROM %POS %ID %REF %ALT %QUAL %FILTER %INFO [%GT]\n' gs://scit1640-clinradomics-55-def/arcus-data/vcf/$line.vcf.gz > /extract_var_v2/$line.txt

done


for line in $(<exome_id.txt)
do

export GCS_OAUTH_TOKEN=$(gcloud auth print-access-token)

bcftools query -i'FILTER="PASS" && INFO/CSQ ~ "probably_damaging" && (GT="0/1" || GT="0|1" || GT="1/0" || GT="1|0" || GT="1/1" || GT="1|1" || GT="1/2" || GT="1|2")' -f'%CHROM %POS %ID %REF %ALT %QUAL %FILTER [%DP %GQ]\n' gs://scit1640-clinradomics-55-def/arcus-data/vcf/$line.vcf.gz > /extract_var_v2_info/$line.txt

done
















