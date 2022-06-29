import pandas as pd
import numpy as np
import joblib
import sys
import argparse

parser = argparse.ArgumentParser(description = "SV2 genotyper")
parser.add_argument("--features_file", help = "Input features .tsv file")
parser.add_argument("--genotype_predictions_output_tsv", help = "Output features .tsv file with genotype predictions")
parser.add_argument("--sex", help = "Sex of sample: male or female")
args = parser.parse_args()

df = pd.read_csv(args.features_file, sep = "\t")
df_male_sex_chromosomes = pd.DataFrame(columns = df.columns)
if args.sex == "male":
    # Filter out the sex chromosome SVs from the main dataframe
    df_male_sex_chromosomes = df[(mask := (df["chrom"].str.contains("X") | df["chrom"].str.contains("Y")))]
    df = df[~mask]

# Autosomal classification

df_highcov_del_gt1kb_values = df[(df["size"] > 1000) & (df["type"] == "DEL")][["chrom", "start", "end", "type", "coverage_GCcorrected", "discordant_ratio", "split_ratio"]].dropna().values
df_highcov_del_lt1kb_values = df[(df["size"] <= 1000) & (df["type"] == "DEL")][["chrom", "start", "end", "type", "coverage_GCcorrected", "discordant_ratio", "split_ratio"]].dropna().values
df_dup_breakpoint_values = df[df["type"] == "DUP"][["chrom", "start", "end", "type", "coverage_GCcorrected", "discordant_ratio", "split_ratio"]].dropna()
df_dup_breakpoint_values = df_dup_breakpoint_values[(df_dup_breakpoint_values["discordant_ratio"] > 0.0) | (df_dup_breakpoint_values["split_ratio"] > 0.0)].values
df_dup_har_values = df[df["type"] == "DUP"][["chrom", "start", "end", "type", "coverage_GCcorrected", "discordant_ratio", "split_ratio", "heterozygous_allele_ratio"]].dropna()
df_dup_har_values = df_dup_har_values[(df_dup_har_values["discordant_ratio"] == 0.0) & (df_dup_har_values["split_ratio"] == 0.0)][["chrom", "start", "end", "type", "coverage_GCcorrected", "heterozygous_allele_ratio"]].values

clf_highcov_del_gt1kb = joblib.load("data/trained_classifiers/clf_del_gt1kb.pkl")
clf_highcov_del_lt1kb = joblib.load("data/trained_classifiers/clf_del_lt1kb.pkl")
clf_dup_breakpoint = joblib.load("data/trained_classifiers//clf_dup_breakpoint.pkl")
clf_dup_har = joblib.load("data/trained_classifiers/clf_dup_har.pkl")

preds_highcov_del_gt1kb = clf_highcov_del_gt1kb.predict_proba(df_highcov_del_gt1kb_values[:, 4:])
preds_highcov_del_lt1kb = clf_highcov_del_lt1kb.predict_proba(df_highcov_del_lt1kb_values[:, 4:])
preds_dup_breakpoint = clf_dup_breakpoint.predict_proba(df_dup_breakpoint_values[:, 4:])
preds_dup_har = clf_dup_har.predict_proba(df_dup_har_values[:, 4:])

df_highcov_del_gt1kb_preds = pd.DataFrame({"chrom": df_highcov_del_gt1kb_values[:, 0], "start": df_highcov_del_gt1kb_values[:, 1], "end": df_highcov_del_gt1kb_values[:, 2], "type": df_highcov_del_gt1kb_values[:, 3], "HOM_GENOTYPE_LIKELIHOOD": preds_highcov_del_gt1kb[:, 0], "HET_GENOTYPE_LIKELIHOOD": preds_highcov_del_gt1kb[:, 1], "REF_GENOTYPE_LIKELIHOOD": preds_highcov_del_gt1kb[:, 2]})
df_highcov_del_lt1kb_preds = pd.DataFrame({"chrom": df_highcov_del_lt1kb_values[:, 0], "start": df_highcov_del_lt1kb_values[:, 1], "end": df_highcov_del_lt1kb_values[:, 2], "type": df_highcov_del_lt1kb_values[:, 3], "HOM_GENOTYPE_LIKELIHOOD": preds_highcov_del_lt1kb[:, 0], "HET_GENOTYPE_LIKELIHOOD": preds_highcov_del_lt1kb[:, 1], "REF_GENOTYPE_LIKELIHOOD": preds_highcov_del_lt1kb[:, 2]})
df_dup_breakpoint_preds = pd.DataFrame({"chrom": df_dup_breakpoint_values[:, 0], "start": df_dup_breakpoint_values[:, 1], "end": df_dup_breakpoint_values[:, 2], "type": df_dup_breakpoint_values[:, 3], "HOM_GENOTYPE_LIKELIHOOD": preds_dup_breakpoint[:, 0], "HET_GENOTYPE_LIKELIHOOD": preds_dup_breakpoint[:, 1], "REF_GENOTYPE_LIKELIHOOD": preds_dup_breakpoint[:, 2]})
df_dup_har_preds = pd.DataFrame({"chrom": df_dup_har_values[:, 0], "start": df_dup_har_values[:, 1], "end": df_dup_har_values[:, 2], "type": df_dup_har_values[:, 3], "HOM_GENOTYPE_LIKELIHOOD": preds_dup_har[:, 0], "HET_GENOTYPE_LIKELIHOOD": preds_dup_har[:, 1], "REF_GENOTYPE_LIKELIHOOD": preds_dup_har[:, 2]})

df_preds_concat_sorted = pd.concat([df_highcov_del_gt1kb_preds, df_highcov_del_lt1kb_preds, df_dup_breakpoint_preds, df_dup_har_preds]).sort_values(["chrom", "start", "end"])
df_preds_concat_sorted = df_preds_concat_sorted.merge(df, on = ["chrom", "start", "end", "type"], how = "left")
df_preds_concat_sorted["ALT_GENOTYPE_LIKELIHOOD"] = df_preds_concat_sorted["HET_GENOTYPE_LIKELIHOOD"] + df_preds_concat_sorted["REF_GENOTYPE_LIKELIHOOD"]
df_preds_concat_sorted["REF_QUAL"] = -10.0*np.log10(1.0 - df_preds_concat_sorted["HOM_GENOTYPE_LIKELIHOOD"])
df_preds_concat_sorted["ALT_QUAL"] = -10.0*np.log10(df_preds_concat_sorted["HOM_GENOTYPE_LIKELIHOOD"])

# Male sex chromosome classifiers
if not df_male_sex_chromosomes.empty:
    df_malesexchrom_del_values = df_male_sex_chromosomes[df_male_sex_chromosomes["type"] == "DEL"][["chrom", "start", "end", "type", "coverage_GCcorrected", "discordant_ratio", "split_ratio"]].dropna().values
    df_malesexchrom_dup_values = df_male_sex_chromosomes[df_male_sex_chromosomes["type"] == "DUP"][["chrom", "start", "end", "type", "coverage_GCcorrected", "discordant_ratio", "split_ratio"]].dropna().values
    
    clf_del_malesexchrom = joblib.load("data/trained_classifiers/clf_del_malesexchrom.pkl")
    clf_dup_malesexchrom = joblib.load("data/trained_classifiers/clf_dup_malesexchrom.pkl")
    
    preds_malesexchrom_del = clf_del_malesexchrom.predict_proba(df_malesexchrom_del_values[:, 4:])
    preds_malesexchrom_dup = clf_dup_malesexchrom.predict_proba(df_malesexchrom_dup_values[:, 4:])
    
    df_malesexchrom_del_preds_values = pd.DataFrame({"chrom": df_malesexchrom_del_values[:, 0], "start": df_malesexchrom_del_values[:, 1], "end": df_malesexchrom_del_values[:, 2], "type": df_malesexchrom_del_values[:, 3], "HOM_GENOTYPE_LIKELIHOOD": preds_malesexchrom_del[:, 0], "REF_GENOTYPE_LIKELIHOOD": preds_malesexchrom_del[:, 1]})
    df_malesexchrom_dup_preds_values = pd.DataFrame({"chrom": df_malesexchrom_dup_values[:, 0], "start": df_malesexchrom_dup_values[:, 1], "end": df_malesexchrom_dup_values[:, 2], "type": df_malesexchrom_dup_values[:, 3], "HOM_GENOTYPE_LIKELIHOOD": preds_malesexchrom_dup[:, 0], "REF_GENOTYPE_LIKELIHOOD": preds_malesexchrom_dup[:, 1]})
    # df_malesexchrom_dup_preds_values["HET_GENOTYPE_LIKELIHOOD"] = float(np.nan) 
    df_malesexchrom_dup_preds_values["HET_GENOTYPE_LIKELIHOOD"] = 0.0
    
    df_malesexchrom_preds_concat_sorted = pd.concat([df_malesexchrom_del_preds_values, df_malesexchrom_dup_preds_values]).sort_values(["chrom", "start", "end"])
    df_malesexchrom_preds_concat_sorted = df_malesexchrom_preds_concat_sorted.merge(df_male_sex_chromosomes, on = ["chrom", "start", "end", "type"], how = "left")
    df_malesexchrom_preds_concat_sorted["ALT_GENOTYPE_LIKELIHOOD"] = df_malesexchrom_preds_concat_sorted["HOM_GENOTYPE_LIKELIHOOD"]
    df_malesexchrom_preds_concat_sorted["REF_QUAL"] = -10.0*np.log10(1.0 - df_malesexchrom_preds_concat_sorted["HOM_GENOTYPE_LIKELIHOOD"])
    df_malesexchrom_preds_concat_sorted["ALT_QUAL"] = -10.0*np.log10(df_malesexchrom_preds_concat_sorted["HOM_GENOTYPE_LIKELIHOOD"])
    
    # Now concat the main dataframe with the male sex chromosome dataframe
    df_preds_concat_sorted = pd.concat([df_preds_concat_sorted, df_malesexchrom_preds_concat_sorted])

df_sort_col = df_preds_concat_sorted[["chrom", "start", "end", "type", "size", "coverage", "coverage_GCcorrected", "discordant_ratio", "split_ratio", "snv_coverage", "heterozygous_allele_ratio", "snvs", "het_snvs", "ALT_GENOTYPE_LIKELIHOOD", "REF_QUAL", "ALT_QUAL", "HOM_GENOTYPE_LIKELIHOOD", "HET_GENOTYPE_LIKELIHOOD", "REF_GENOTYPE_LIKELIHOOD"]]

# Initialize the GEN column to the missing genotype value ./.
df_sort_col["GEN"] = "./."
# Find all rows that fulfill the conditions and set its corresponding genotype
df_sort_col.loc[(df_sort_col['HOM_GENOTYPE_LIKELIHOOD'] > df_sort_col['HET_GENOTYPE_LIKELIHOOD']) &  
       (df_sort_col['HOM_GENOTYPE_LIKELIHOOD'] > df_sort_col['REF_GENOTYPE_LIKELIHOOD']) ,
       'GEN'] = "1/1"
df_sort_col.loc[(df_sort_col['HET_GENOTYPE_LIKELIHOOD'] > df_sort_col['HOM_GENOTYPE_LIKELIHOOD']) &
       (df_sort_col['HET_GENOTYPE_LIKELIHOOD'] > df_sort_col['REF_GENOTYPE_LIKELIHOOD']) ,
       'GEN'] = "0/1"
df_sort_col.loc[(df_sort_col['REF_GENOTYPE_LIKELIHOOD'] > df_sort_col['HOM_GENOTYPE_LIKELIHOOD']) &
       (df_sort_col['REF_GENOTYPE_LIKELIHOOD'] > df_sort_col['HET_GENOTYPE_LIKELIHOOD']) ,
       'GEN'] = "0/0"


if not args.genotype_predictions_output_tsv:
    from time import gmtime, strftime
    current_time = strftime("%Y-%m-%d_%H.%M.%S", gmtime())
    df_sort_col.to_csv("genotyping_preds_{}.tsv".format(current_time), sep = "\t", index = False)
else: df_sort_col.to_csv(args.genotype_predictions_output_tsv, sep = "\t", index = False)

# TODO:
# Reorganize code (such that predictions are done after each df is created first, rather than doing them all at end)
# Some if/else code (for males, empty arrays, etc.)
# Merge dup dataframes?
