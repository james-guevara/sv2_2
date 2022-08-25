# SV2_2
This is a refactored version of SV2. Tested with Python 3.9.7. We recommend creating a separate conda/miniconda (or mamba/micromamba) environment with this version of Python, and then using your created environment's `pip` command to install the requirements, i.e. `pip install --requirements.txt`

Example run command:
```
python run_sv2.py --alignment_file /expanse/projects/sebat1/genomicsdataanalysis/resources/NA12878/illumina_platinum_pedigree/NA12878.alt_bwamem_GRCh38DH.20150706.CEU.illumina_platinum_ped.cram \
                  --reference_fasta /expanse/lustre/projects/ddp195/j3guevar/resources/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa \
                  --snv_vcf_file /home/j3guevar/tests/SV2_nim/src/SV2_nimpkg/refactored_sv2/input_files/NA12878.vcf.gz \
                  --exclude_regions_bed data/excluded_regions_bed_files/hg38_excluded.bed.gz \
                  --sv_bed_file data/test_data/file.nbl.test2.bed \
                  --sex male \
                  --sample_name NA12878
```

As of now, we're targeting only reference genome hg38. We can also make it usable for hg19 and other reference types later.

## To do
- [ ] Test on `HG002` and `NA12878` and some of `REACH` cohort
- [ ] Run from different compute environments/folders
- [ ] Once ROC curves are generated and analyzed, we have to decide on FILTER flag criteria
- [ ] There might be (un)subtle bugs when generating QUAL scores or even genotype likelihoods, need to catch those
- [x] Add SV2 commands to VCF 
- [x] Default output files should include sample name 
- [x] Automatic find the `regions_bed` and `gc_reference_table` (reduce number of command line arguments)
- [x] Make exclude bed file optional? 
- [x] Put requirements.txt file in here 
- [x] Test when run_sv2.py is a symbolic link
