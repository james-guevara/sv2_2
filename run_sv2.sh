#python run_sv2.py --alignment_file /expanse/projects/sebat1/genomicsdataanalysis/resources/NA12878/illumina_platinum_pedigree/NA12878.alt_bwamem_GRCh38DH.20150706.CEU.illumina_platinum_ped.cram \
#                  --reference_fasta /expanse/lustre/projects/ddp195/j3guevar/resources/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa \
#                  --snv_vcf_file /home/j3guevar/tests/SV2_nim/src/SV2_nimpkg/refactored_sv2/input_files/NA12878.vcf.gz \
#                  --exclude_regions_bed data/excluded_regions_bed_files/hg38_excluded.bed.gz \
#                  --sv_bed_file data/test_data/file.nbl.bed \
#                  --sex male \
#                  --sample_name NA12878

#python run_sv2.py --alignment_file /expanse/lustre/projects/ddp195/eiovino/sv_pipeline_nf/REACH000236.chr22.cram \
#                  --reference_fasta /expanse/lustre/projects/ddp195/eiovino/fasta/Homo_sapiens_assembly38.fasta \
#                  --snv_vcf_file /expanse/projects/sebat1/genomicsdataanalysis/REACH_MSSNG/request_202001/REACH000236.vcf.gz \
#                  --exclude_regions_bed data/excluded_regions_bed_files/hg38_excluded.bed.gz \
#                  --sv_bed_file REACH_chr22.bed \
#                  --sex male \
#                  --sample_name REACH000236

#python run_sv2.py --alignment_file /expanse/lustre/projects/ddp195/eiovino/cram_HG002/HG002.100_reads.markdup.recal.cram \
#                  --reference_fasta /expanse/lustre/projects/ddp195/j3guevar/resources/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa \
#                  --snv_vcf_file hg002/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
#                  --exclude_regions_bed data/excluded_regions_bed_files/hg38_excluded.bed.gz \
#                  --sv_bed_file hg002/hg002_results/HG002.merged.sorted.bed \
#                  --sex male \
#                  --sample_name HG002

#python run_sv2.py --alignment_file /expanse/lustre/projects/ddp195/eiovino/cram_HG002/HG002.100_reads.markdup.recal.cram \
#                  --reference_fasta /expanse/lustre/projects/ddp195/j3guevar/resources/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa \
#                  --snv_vcf_file hg002/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
#                  --exclude_regions_bed data/excluded_regions_bed_files/hg38_excluded.bed.gz \
#                  --sv_bed_file hg002/hg002_results/HG002_head.merged.sorted.bed \
#                  --sex male \
#                  --sample_name HG002_head_2 \

#python run_sv2.py --alignment_file /expanse/projects/sebat1/genomicsdataanalysis/resources/NA12878/illumina_platinum_pedigree/NA12878.alt_bwamem_GRCh38DH.20150706.CEU.illumina_platinum_ped.cram \
#                  --reference_fasta /expanse/lustre/projects/ddp195/j3guevar/resources/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa \
#                  --snv_vcf_file /home/j3guevar/tests/SV2_nim/src/SV2_nimpkg/refactored_sv2/input_files/NA12878.vcf.gz \
#                  --exclude_regions_bed data/excluded_regions_bed_files/hg38_excluded.bed.gz \
#                  --sv_bed_file data/test_data/file.nbl.test.bed \
#                  --sex male \
#                  --sample_name NA12878 \

#python run_sv2.py --alignment_file /expanse/projects/sebat1/genomicsdataanalysis/resources/NA12878/illumina_platinum_pedigree/NA12878.alt_bwamem_GRCh38DH.20150706.CEU.illumina_platinum_ped.cram \
#                  --reference_fasta /expanse/lustre/projects/ddp195/j3guevar/resources/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa \
#                  --snv_vcf_file /home/j3guevar/tests/SV2_nim/src/SV2_nimpkg/refactored_sv2/input_files/NA12878.vcf.gz \
#                  --sv_bed_file data/test_data/file.nbl.test.bed \
#                  --sex male \
#                  --sample_name NA12878 \

#python run_sv2.py --alignment_file /expanse/lustre/projects/ddp195/eiovino/cram_HG002/HG002.100_reads.markdup.recal.cram \
#                  --reference_fasta /expanse/lustre/projects/ddp195/j3guevar/resources/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa \
#                  --snv_vcf_file hg002/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
#                  --exclude_regions_bed data/excluded_regions_bed_files/hg38_excluded.bed.gz \
#                  --sv_bed_file hg002/hg002_results/HG002.merged.sorted.bed \
#                  --sex male \
#                  --sample_name HG002

python run_sv2.py --alignment_file /expanse/lustre/projects/ddp195/eiovino/cram_HG002/HG002.100_reads.markdup.recal.cram \
                  --reference_fasta /expanse/lustre/projects/ddp195/j3guevar/resources/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa \
                  --snv_vcf_file hg002/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
                  --sv_bed_file HG002.mod.bed \
                  --sex male \
                  --sample_name HG002 \
                  --exclude_regions_bed data/excluded_regions_bed_files/hg38_excluded.bed.gz
