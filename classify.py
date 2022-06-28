import pandas as pd
import numpy as np
import joblib
import sys
import argparse

parser = argparse.ArgumentParser(description = "SV2 genotyper")
parser.add_argument("--features_file", help = "Input features file")
parser.add_argument("--sex", help = "Sex of sample: male or female")
args = parser.parse_args()

df = pd.read_csv(args.features_file, sep = "\t")
if sex == "male":
    # Filter out the sex chromosome SVs from the main dataframe
    df = df[~(df["chrom"].str.contains("X") ! df["chrom"].str.contains("Y")]

df_highcov_del_gt1kb_values = df[(df["size"] > 1000) & (df["type"] == "DEL")][["chrom", "start", "end", "type", "coverage_GCcorrected", "discordant_ratio", "split_ratio"]].dropna().values
df_highcov_del_lt1kb_values = df[(df["size"] <= 1000) & (df["type"] == "DEL")][["chrom", "start", "end", "type", "coverage_GCcorrected", "discordant_ratio", "split_ratio"]].dropna().values
df_dup_breakpoint_values = df[df["type"] == "DUP"][["chrom", "start", "end", "type", "coverage_GCcorrected", "discordant_ratio", "split_ratio"]].dropna().values
df_dup_har_values = df[df["type"] == "DUP"][["chrom", "start", "end", "type", "coverage_GCcorrected", "heterozygous_allele_ratio"]].dropna().values

clf_highcov_del_gt1kb = joblib.load("training/clf_del_gt1kb.pkl")
clf_highcov_del_lt1kb = joblib.load("training/clf_del_lt1kb.pkl")
clf_dup_breakpoint = joblib.load("training/clf_dup_breakpoint.pkl")
clf_dup_har = joblib.load("training/clf_dup_har.pkl")

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

from time import gmtime, strftime
current_time = strftime("%Y-%m-%d_%H.%M.%S", gmtime())
df_sort_col.to_csv("genotyping_preds{}.tsv".format(current_time), sep = "\t", index = False)


# TODO:
# Reorganize code (such that predictions are done after each df is created first, rather than doing them all at end)
# Some if/else code (for males, empty arrays, etc.)
# Merge dup dataframes?
