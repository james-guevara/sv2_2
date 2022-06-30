import argparse
import datetime
import numpy as np
import pandas as pd
from pathlib import Path
from pybedtools import BedTool
import sys
from time import gmtime, strftime

# For step 1
from make_feature_table import make_GC_content_reference_table
from make_feature_table import make_regions_table  
from make_feature_table import make_alignment_preprocessing_table  
from make_feature_table import get_snv_preprocessing_data
from make_feature_table import make_sv_interval_table 
from make_feature_table import make_snv_features_table 
from make_feature_table import make_alignment_features_table
# For step 2
from classify import load_features
from classify import load_features_from_dataframe
from classify import run_highcov_del_gt1kb_classifier 
from classify import run_highcov_del_lt1kb_classifier 
from classify import run_dup_breakpoint_classifier
from classify import run_dup_har_classifier 
from classify import concat_and_sort_pred_dfs
from classify import run_malesexchrom_del_classifier 
from classify import run_malesexchrom_dup_classifier 
from classify import get_genotype 
# For step 3
from make_vcf import make_vcf

parser = argparse.ArgumentParser(description = "SV2 genotyper")
# From make_feature_table.py step
parser.add_argument("--alignment_file", help = "CRAM/BAM file input", required = True)
parser.add_argument("--reference_fasta", help = "Reference fasta file input", required = True)
parser.add_argument("--snv_vcf_file", help = "SNV VCF file input", required = True)
parser.add_argument("--regions_bed", help = "BED file with pre-generated random genomic regions (for estimating coverage per chromosome)")
parser.add_argument("--exclude_regions_bed", help = "BED file with regions to exclude", required = True)
parser.add_argument("--sv_bed_file", help = "SV BED file input", required = True)
parser.add_argument("--preprocessing_table_input", help = "A pre-generated preprocessing table (if SV2 had been run before and you want to skip the preprocessing part of the program")
parser.add_argument("--gc_reference_table", help = "GC content reference table input")
# From classify.py step
parser.add_argument("--sex", help = "Sex of sample: male or female", required = True)
parser.add_argument("--clf_folder", help = "Folder that contains the classifiers, which must be in .pkl format (if not specified, will look for them in the default data folder)")
# From make_vcf.py step
parser.add_argument("--sample_name", help = "Sample name", required = True)
parser.add_argument("--output_folder", help = "Output folder for output files (if not used, then output folder is set to 'sv2_output')")
args = parser.parse_args()

# If regions_bed or gc_reference_table aren't specified, we should look for them using the absolute path of this script (and these files should be found in the data folder, which is inside the same parent folder as this script)
script_folder = Path(__file__).parent.absolute()
regions_bed = str()
if args.regions_bed: regions_bed = args.regions_bed
else: regions_bed = "{}/data/random_regions_hg38.bed".format(script_folder)
gc_reference_table = str()
if args.gc_reference_table: gc_reference_table = args.gc_reference_table
else: gc_reference_table = "{}/data/GC_content_reference.txt".format(script_folder)


# If one of the output args isn't specified, then we should make a folder called "output" to put the output files (if the folder doesn't already exist)
output_folder = Path("sv2_output")
if args.output_folder: output_folder = Path(args.output_folder)
output_folder.mkdir(parents = True, exist_ok = True)

""" make_feature_table.py step """

# Chromosomes to analyze (1,... 22, X, Y)
chroms = list(map(str, np.arange(1, 22 + 1)))
chroms.extend(["X", "Y"])

# SV types to classify
svtypes = ("DEL", "DUP")

# Make the GC content reference table
GC_content_reference_table = make_GC_content_reference_table(gc_reference_table)
# Make the regions table (used to estimate coverage)
regions_table = make_regions_table(regions_bed)

if not args.preprocessing_table_input:
    # Make the CRAM/BAM preprocessing table
    alignment_preprocessing_table = make_alignment_preprocessing_table(args.alignment_file, args.reference_fasta, chroms, regions_table)
    
    # Make the SNV preprocessing table
    snv_preprocessing_table = get_snv_preprocessing_data(args.snv_vcf_file, chroms, regions_table)

    # Merge and output the preprocessing table
    df_alignment_preprocessing_table = pd.DataFrame.from_dict(alignment_preprocessing_table, orient = "index")
    df_snv_preprocessing_table = pd.DataFrame.from_dict(snv_preprocessing_table, orient = "index")
    df_preprocessing_table = df_alignment_preprocessing_table.join(df_snv_preprocessing_table).reset_index(level = 0).rename(columns = {"index": "chrom"})

    # Save preprocessing table to output folder
    current_time = strftime("%Y-%m-%d_%H.%M.%S", gmtime())
    preprocessing_table_filepath =  "{}/{}_sv2_preprocessing_features_{}.tsv".format(output_folder, args.sample_name, current_time)
    df_preprocessing_table.to_csv(preprocessing_table_filepath, sep = "\t", index = False)
else:
    df_preprocessing_table = pd.read_csv(args.preprocessing_table_input, sep = "\t")

# Preprocess the input SV BED file
sv_bed_list = []
with open(args.sv_bed_file, "r") as f:
    for line in f:
        linesplit = line.rstrip().split("\t")
        chrom, start, stop, features = linesplit[0], int(linesplit[1]), int(linesplit[2]), linesplit[3]
        if start > stop: continue
        sv_bed_list.append(line.rstrip())
sv_bed = BedTool(sv_bed_list).filter(lambda x: len(x) > 0).saveas()
exclude_bed = BedTool(args.exclude_regions_bed).merge()
sv_interval_table = make_sv_interval_table(sv_bed, exclude_bed, args.reference_fasta)

# Make SNV features table (for each filtered SV call)
snv_features_table = make_snv_features_table(args.snv_vcf_file, sv_bed, sv_interval_table, svtypes, df_preprocessing_table)

# Make CRAM/BAM features table (for each filtered SV call)
alignment_features_table = make_alignment_features_table(args.alignment_file, args.reference_fasta, sv_bed, df_preprocessing_table, sv_interval_table, svtypes, GC_content_reference_table)

# Merge and output feature table
df_alignment_features_table = pd.DataFrame.from_dict(alignment_features_table, orient = "index")
df_snv_features_table = pd.DataFrame.from_dict(snv_features_table, orient = "index")
df_features_table = df_alignment_features_table.join(df_snv_features_table).reset_index().rename(columns = {"level_0": "chrom", "level_1": "start", "level_2": "end"})

# Save SV features table
current_time = strftime("%Y-%m-%d_%H.%M.%S", gmtime())
features_table_filepath =  "{}/{}_sv2_features_{}.tsv".format(output_folder, args.sample_name, current_time)
df_features_table.to_csv(features_table_filepath, sep = "\t", index = False)

""" classify.py step """

# df, df_male_sex_chromosomes = load_features_from_dataframe(df_features_table, args.sex)
df, df_male_sex_chromosomes = load_features(features_table_filepath, args.sex)

# Set default filepaths for the classifiers
clf_highcov_del_gt1kb_filepath = "data/trained_classifiers/clf_del_gt1kb.pkl"
clf_highcov_del_lt1kb_filepath = "data/trained_classifiers/clf_del_lt1kb.pkl"
clf_dup_breakpoint_filepath = "data/trained_classifiers/clf_dup_breakpoint.pkl"
clf_dup_har_filepath = "data/trained_classifiers/clf_dup_har.pkl"
clf_del_malesexchrom_filepath = "data/trained_classifiers/clf_del_malesexchrom.pkl"
clf_dup_malesexchrom_filepath = "data/trained_classifiers/clf_dup_malesexchrom.pkl"

# If classifier folder is provided as argument, use the ones inside of that folder (must have the same names as default classifiers)
if args.clf_folder:
    clf_highcov_del_gt1kb_filepath = "{}/clf_del_gt1kb.pkl".format(args.clf_folder)
    clf_highcov_del_lt1kb_filepath = "{}/clf_del_lt1kb.pkl".format(args.clf_folder)
    clf_dup_breakpoint_filepath = "{}/clf_dup_breakpoint.pkl".format(args.clf_folder)
    clf_dup_har_filepath = "{}/clf_dup_har.pkl".format(args.clf_folder)
    clf_del_malesexchrom_filepath = "{}/clf_del_malesexchrom.pkl".format(args.clf_folder)
    clf_dup_malesexchrom_filepath = "{}/clf_dup_malesexchrom.pkl".format(args.clf_folder)


# Autosomal SV classifiers
df_highcov_del_gt1kb_preds = run_highcov_del_gt1kb_classifier(df, clf_highcov_del_gt1kb_filepath)
df_highcov_del_lt1kb_preds = run_highcov_del_lt1kb_classifier(df, clf_highcov_del_lt1kb_filepath)
df_dup_breakpoint_preds = run_dup_breakpoint_classifier(df, clf_dup_breakpoint_filepath)
df_dup_har_preds = run_dup_har_classifier(df, clf_dup_har_filepath)

# Concat and sort output diploid dfs
df_preds_concat_sorted = concat_and_sort_pred_dfs([df_highcov_del_gt1kb_preds, df_highcov_del_lt1kb_preds, df_dup_breakpoint_preds, df_dup_har_preds], df)

# Male sex chromosome SV classifiers
df_malesexchrom_del_preds = run_malesexchrom_del_classifier(df_male_sex_chromosomes, clf_del_malesexchrom_filepath) 
df_malesexchrom_dup_preds = run_malesexchrom_dup_classifier(df_male_sex_chromosomes, clf_dup_malesexchrom_filepath) 

# Concat and sort output haploid dfs
df_malesexchrom_preds_concat_sorted = concat_and_sort_pred_dfs([df_malesexchrom_del_preds, df_malesexchrom_dup_preds], df_male_sex_chromosomes)

# Concat the diploid and haploid dfs, and then add quality scores
df_preds_concat_sorted = pd.concat([df_preds_concat_sorted, df_malesexchrom_preds_concat_sorted])
df_preds_concat_sorted["ALT_GENOTYPE_LIKELIHOOD"] = df_preds_concat_sorted["HET_GENOTYPE_LIKELIHOOD"] + df_preds_concat_sorted["HOM_GENOTYPE_LIKELIHOOD"]
df_preds_concat_sorted["REF_QUAL"] = -10.0*np.log10(1.0 - df_preds_concat_sorted["REF_GENOTYPE_LIKELIHOOD"])
df_preds_concat_sorted["ALT_QUAL"] = -10.0*np.log10(df_preds_concat_sorted["REF_GENOTYPE_LIKELIHOOD"])
df_preds_concat_sorted = df_preds_concat_sorted[["chrom", "start", "end", "type", "size", "coverage", "coverage_GCcorrected", "discordant_ratio", "split_ratio", "snv_coverage", "heterozygous_allele_ratio", "snvs", "het_snvs", "ALT_GENOTYPE_LIKELIHOOD", "REF_QUAL", "ALT_QUAL", "HOM_GENOTYPE_LIKELIHOOD", "HET_GENOTYPE_LIKELIHOOD", "REF_GENOTYPE_LIKELIHOOD"]]

# Initialize the GEN column to the missing genotype value ./.
df_preds_concat_sorted["GEN"] = "./."
df_preds_concat_sorted["GEN"] = df_preds_concat_sorted[["REF_GENOTYPE_LIKELIHOOD", "HET_GENOTYPE_LIKELIHOOD", "HOM_GENOTYPE_LIKELIHOOD"]].apply(lambda x: get_genotype(*x), axis = 1)

# Save genotype predictions table
current_time = strftime("%Y-%m-%d_%H.%M.%S", gmtime())
genotype_table_filepath = "{}/{}_genotyping_preds_{}.tsv".format(output_folder, args.sample_name, current_time)
df_preds_concat_sorted.to_csv(genotype_table_filepath, sep = "\t", index = False)

""" make_vcf.py step """
make_vcf(args.sample_name, args.reference_fasta, genotype_table_filepath, output_folder)
