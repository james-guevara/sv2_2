import pandas as pd
import numpy as np
import joblib
import sys
import argparse

# Steps (since I'll be running this from a top-level script, I must use function definitions for these parts)
# def load_features (will create df and df_male_sex_chromosomes) and so forth
# Load dataframe and modify if sex is male 
# Run each classifier in separate function and return as dataframe
# After each classifier is run, concat and sort the output dataframe, and postprocess (e.g. create quality score column and genotype column)
# Perhaps keep male sex chromosome part separate from the rest

def load_features(features_table_filepath, sex):
    df = pd.read_csv(features_table_filepath, sep = "\t")
    df_male_sex_chromosomes = pd.DataFrame(columns = df.columns)
    if sex == "male":
        # Filter out the sex chromosome SVs from the main dataframe
        df_male_sex_chromosomes = df[(mask := (df["chrom"].str.contains("X") | df["chrom"].str.contains("Y")))]
        df = df[~mask]
    return df, df_male_sex_chromosomes

def load_features_from_dataframe(df_original, sex)
    df = df_original 
    df_male_sex_chromosomes = pd.DataFrame(columns = df.columns)
    if sex == "male":
        # Filter out the sex chromosome SVs from the main dataframe
        df_male_sex_chromosomes = df[(mask := (df["chrom"].str.contains("X") | df["chrom"].str.contains("Y")))]
        df = df[~mask]
    return df, df_male_sex_chromosomes

# Run the classifier for the high-coverage deletions greater than 1KB in length
def run_highcov_del_gt1kb_classifier(df, clf_highcov_del_gt1kb_filepath):
    df_highcov_del_gt1kb_values = df[(df["size"] > 1000) & (df["type"] == "DEL")][["chrom", "start", "end", "type", "coverage_GCcorrected", "discordant_ratio", "split_ratio"]].dropna().values
    clf_highcov_del_gt1kb = joblib.load(clf_highcov_del_gt1kb_filepath)
    preds_highcov_del_gt1kb = clf_highcov_del_gt1kb.predict_proba(df_highcov_del_gt1kb_values[:, 4:])
    df_highcov_del_gt1kb_preds = pd.DataFrame({"chrom": df_highcov_del_gt1kb_values[:, 0], "start": df_highcov_del_gt1kb_values[:, 1], "end": df_highcov_del_gt1kb_values[:, 2], "type": df_highcov_del_gt1kb_values[:, 3], "HOM_GENOTYPE_LIKELIHOOD": preds_highcov_del_gt1kb[:, 0], "HET_GENOTYPE_LIKELIHOOD": preds_highcov_del_gt1kb[:, 1], "REF_GENOTYPE_LIKELIHOOD": preds_highcov_del_gt1kb[:, 2]})
    return df_highcov_del_gt1kb_preds

# Likewise for high-coverage deletions less than 1KB in length...
def run_highcov_del_lt1kb_classifier(df, clf_highcov_del_lt1kb_filepath):
    df_highcov_del_lt1kb_values = df[(df["size"] <= 1000) & (df["type"] == "DEL")][["chrom", "start", "end", "type", "coverage_GCcorrected", "discordant_ratio", "split_ratio"]].dropna().values
    clf_highcov_del_lt1kb = joblib.load(clf_highcov_del_lt1kb_filepath)
    preds_highcov_del_lt1kb = clf_highcov_del_lt1kb.predict_proba(df_highcov_del_lt1kb_values[:, 4:])
    df_highcov_del_lt1kb_preds = pd.DataFrame({"chrom": df_highcov_del_lt1kb_values[:, 0], "start": df_highcov_del_lt1kb_values[:, 1], "end": df_highcov_del_lt1kb_values[:, 2], "type": df_highcov_del_lt1kb_values[:, 3], "HOM_GENOTYPE_LIKELIHOOD": preds_highcov_del_lt1kb[:, 0], "HET_GENOTYPE_LIKELIHOOD": preds_highcov_del_lt1kb[:, 1], "REF_GENOTYPE_LIKELIHOOD": preds_highcov_del_lt1kb[:, 2]})
    return df_highcov_del_lt1kb_preds

def run_dup_breakpoint_classifier(df, clf_dup_breakpoint_filepath):
    df_dup_breakpoint_values = df[df["type"] == "DUP"][["chrom", "start", "end", "type", "coverage_GCcorrected", "discordant_ratio", "split_ratio"]].dropna()
    df_dup_breakpoint_values = df_dup_breakpoint_values[(df_dup_breakpoint_values["discordant_ratio"] > 0.0) | (df_dup_breakpoint_values["split_ratio"] > 0.0)].values
    clf_dup_breakpoint = joblib.load(clf_dup_breakpoint_filepath)
    preds_dup_breakpoint = clf_dup_breakpoint.predict_proba(df_dup_breakpoint_values[:, 4:])
    df_dup_breakpoint_preds = pd.DataFrame({"chrom": df_dup_breakpoint_values[:, 0], "start": df_dup_breakpoint_values[:, 1], "end": df_dup_breakpoint_values[:, 2], "type": df_dup_breakpoint_values[:, 3], "HOM_GENOTYPE_LIKELIHOOD": preds_dup_breakpoint[:, 0], "HET_GENOTYPE_LIKELIHOOD": preds_dup_breakpoint[:, 1], "REF_GENOTYPE_LIKELIHOOD": preds_dup_breakpoint[:, 2]})
    return df_dup_breakpoint_preds

def run_dup_har_classifier(df, clf_dup_har_filepath):
    df_dup_har_values = df[df["type"] == "DUP"][["chrom", "start", "end", "type", "coverage_GCcorrected", "discordant_ratio", "split_ratio", "heterozygous_allele_ratio"]].dropna()
    df_dup_har_values = df_dup_har_values[(df_dup_har_values["discordant_ratio"] == 0.0) & (df_dup_har_values["split_ratio"] == 0.0)][["chrom", "start", "end", "type", "coverage_GCcorrected", "heterozygous_allele_ratio"]].values
    clf_dup_har = joblib.load(clf_dup_har_filepath)
    preds_dup_har = clf_dup_har.predict_proba(df_dup_har_values[:, 4:])
    df_dup_har_preds = pd.DataFrame({"chrom": df_dup_har_values[:, 0], "start": df_dup_har_values[:, 1], "end": df_dup_har_values[:, 2], "type": df_dup_har_values[:, 3], "HOM_GENOTYPE_LIKELIHOOD": preds_dup_har[:, 0], "HET_GENOTYPE_LIKELIHOOD": preds_dup_har[:, 1], "REF_GENOTYPE_LIKELIHOOD": preds_dup_har[:, 2]})
    return df_dup_har_preds

def concat_and_sort_pred_dfs(df_list, df):
    df_preds_concat_sorted = pd.concat(df_list).sort_values(["chrom", "start", "end"])
    df_preds_concat_sorted = df_preds_concat_sorted.merge(df, on = ["chrom", "start", "end", "type"], how = "left")
    return df_preds_concat_sorted

def run_malesexchrom_del_classifier(df_male_sex_chromosomes, clf_del_malesexchrom_filepath):
    df_malesexchrom_del_values = df_male_sex_chromosomes[df_male_sex_chromosomes["type"] == "DEL"][["chrom", "start", "end", "type", "coverage_GCcorrected", "discordant_ratio", "split_ratio"]].dropna().values
    clf_del_malesexchrom = joblib.load(clf_del_malesexchrom_filepath)
    preds_malesexchrom_del = clf_del_malesexchrom.predict_proba(df_malesexchrom_del_values[:, 4:]) # Only 2 predicted genotypes (because it's haploid)
    df_malesexchrom_del_preds = pd.DataFrame({"chrom": df_malesexchrom_del_values[:, 0], "start": df_malesexchrom_del_values[:, 1], "end": df_malesexchrom_del_values[:, 2], "type": df_malesexchrom_del_values[:, 3], "HOM_GENOTYPE_LIKELIHOOD": preds_malesexchrom_del[:, 0], "REF_GENOTYPE_LIKELIHOOD": preds_malesexchrom_del[:, 1]})
    df_malesexchrom_del_preds["HET_GENOTYPE_LIKELIHOOD"] = 0.0
    return df_malesexchrom_del_preds

def run_malesexchrom_dup_classifier(df_male_sex_chromosomes, clf_dup_malesexchrom_filepath):
    df_malesexchrom_dup_values = df_male_sex_chromosomes[df_male_sex_chromosomes["type"] == "DUP"][["chrom", "start", "end", "type", "coverage_GCcorrected", "discordant_ratio", "split_ratio"]].dropna().values
    clf_dup_malesexchrom = joblib.load(clf_dup_malesexchrom_filepath)
    preds_malesexchrom_dup = clf_dup_malesexchrom.predict_proba(df_malesexchrom_dup_values[:, 4:])
    df_malesexchrom_dup_preds = pd.DataFrame({"chrom": df_malesexchrom_dup_values[:, 0], "start": df_malesexchrom_dup_values[:, 1], "end": df_malesexchrom_dup_values[:, 2], "type": df_malesexchrom_dup_values[:, 3], "HOM_GENOTYPE_LIKELIHOOD": preds_malesexchrom_dup[:, 0], "REF_GENOTYPE_LIKELIHOOD": preds_malesexchrom_dup[:, 1]})
    df_malesexchrom_dup_preds["HET_GENOTYPE_LIKELIHOOD"] = 0.0
    return df_malesexchrom_dup_preds

def get_genotype(ref_likelihood, het_likelihood, hom_likelihood):
    max_likelihood_index = np.argmax([ref_likelihood, het_likelihood, hom_likelihood])
    if max_likelihood_index == 0: return "0/0"
    elif max_likelihood_index == 1: return "0/1"
    else: return "1/1"

# if not args.genotype_predictions_output_tsv:
#     from time import gmtime, strftime
#     current_time = strftime("%Y-%m-%d_%H.%M.%S", gmtime())
#     df_preds_concat_sorted.to_csv("genotyping_preds_{}.tsv".format(current_time), sep = "\t", index = False)
# else: df_preds_concat_sorted.to_csv(args.genotype_predictions_output_tsv, sep = "\t", index = False)


if __name__ == "__main__":
    # Parse arugments
    parser = argparse.ArgumentParser(description = "SV2 genotyper")
    parser.add_argument("--features_file", help = "Input features .tsv file")
    parser.add_argument("--genotype_predictions_output_tsv", help = "Output features .tsv file with genotype predictions")
    parser.add_argument("--sex", help = "Sex of sample: male or female")
    # Load trained classifiers (.pkl files). In the future, specify the folder where the classifiers are stored (they must have the same names as default classifiers so that the program can load them) 
    parser.add_argument("--clf_folder", help = "Folder that contains the classifiers (classifiers must be in .pkl format)")

    # Set default filepaths for the classifiers
    clf_highcov_del_gt1kb_filepath = "data/trained_classifiers/clf_del_gt1kb.pkl"
    clf_highcov_del_lt1kb_filepath = "data/trained_classifiers/clf_del_lt1kb.pkl"
    clf_dup_breakpoint_filepath = "data/trained_classifiers/clf_dup_breakpoint.pkl"
    clf_dup_har_filepath = "data/trained_classifiers/clf_dup_har.pkl"
    clf_del_malesexchrom_filepath = "data/trained_classifiers/clf_del_malesexchrom.pkl"
    clf_dup_malesexchrom_filepath = "data/trained_classifiers/clf_dup_malesexchrom.pkl"

    args = parser.parse_args()

    # If classifier filepaths are provided as arguments, use the ones provided as arguments
    if args.clf_folder:
        clf_highcov_del_gt1kb_filepath = "{}/clf_del_gt1kb.pkl".format(args.clf_folder)
        clf_highcov_del_lt1kb_filepath = "{}/clf_del_lt1kb.pkl".format(args.clf_folder)
        clf_dup_breakpoint_filepath = "{}/clf_dup_breakpoint.pkl".format(args.clf_folder)
        clf_dup_har_filepath = "{}/clf_dup_har.pkl".format(args.clf_folder)
        clf_del_malesexchrom_filepath = "{}/clf_del_malesexchrom.pkl".format(args.clf_folder)
        clf_dup_malesexchrom_filepath = "{}/clf_dup_malesexchrom.pkl".format(args.clf_folder)

    df, df_male_sex_chromosomes = load_features(args.features_file, args.sex)

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

    if not args.genotype_predictions_output_tsv:
        from time import gmtime, strftime
        current_time = strftime("%Y-%m-%d_%H.%M.%S", gmtime())
        df_preds_concat_sorted.to_csv("genotyping_preds_{}.tsv".format(current_time), sep = "\t", index = False)
    else: df_preds_concat_sorted.to_csv(args.genotype_predictions_output_tsv, sep = "\t", index = False)
