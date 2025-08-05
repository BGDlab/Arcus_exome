Manual for Whole-Exome Sequencing Data Analysis (VCF-Based)

**1. Sample overview.**

The dataset for this project originates from whole-exome sequencing (WES) experiments performed in the Arcus Lab. All sequencing results are provided in VCF (Variant Call Format) files, which capture variant-level data across the exome.

As of May 15, 2025, the dataset includes a total of 4547 WES samples, corresponding to 4333 unique individuals. Due to multiple sequencing runs for some participants, several individuals have more than one sample.

Sample breakdown:
Total WES samples: 4547

Children (primary analysis group): 4182 samples from 3982 unique child participants

Parents: 341 samples

Siblings: 10 samples

(Note: Some children may have multiple samples due to technical replicates or re-sequencing)

In this project, we focus exclusively on the 4182 child samples for downstream genetic analyses.

Important: The list of child sample IDs is provided in the file ABC.txt.

**2. Genetic quality control workflow.**

All exome sequencing samples in Arcus Lab have been pre-annotated using VEP (Variant Effect Predictor). For detailed information on VEP annotations and field definitions, please refer to:
https://useast.ensembl.org/info/docs/tools/vep/index.html

**2.1 Variant filtering criteria** (refer to code.sh)

For downstream analyses, we retain only high-confidence variants and genotypes of interest according to the following criteria:

Filter:
Variants must pass the internal quality filters applied during variant calling:

"FILTER = PASS" in the VCF record.

Genotype (GT) Selection:
We retain genotypes with the following values from the "FORMAT:GT" field:

0/0, 0|0  
0/1, 0|1, 1/0, 1|0  
1/1, 1|1  
1/2, 1|2, 2/1, 2|1

“/” = Unphased genotype (origin of alleles unknown)

“|” = Phased genotype (alleles assigned to parental origin)

0/0: Homozygous for the reference allele

0/1 or 1/0: Heterozygous (REF and ALT)

1/1: Homozygous for the ALT allele

1/2: Heterozygous for two ALT alleles (multiallelic site)

**2.2 Functional prediction filter** (refer to code.sh)

We focus exclusively on variants predicted to be potentially damaging: Include only variants annotated by PolyPhen as:

PolyPhen = probably_damaging
(found in the INFO field of the VCF)

**2.3 Clinical relevance filter** (refer to code.R)

To prioritize clinically significant variants: Include only variants with ClinVar annotations:

CLIN_SIG = likely_pathogenic or pathogenic
(also found in the INFO field of the VCF)

**2.4 Allele frequency filter (rare variants only)** (refer to code.R)

To enrich for rare variants, we exclude common polymorphisms:

Retain variants where:

MAX_AF < 0.01 (i.e., <1% allele frequency)

MAX_AF represents the maximum observed allele frequency across multiple population databases including:

gnomAD

1000 Genomes

ESP

The resulting genetic variants are included in the follow-up genetic analysis.
