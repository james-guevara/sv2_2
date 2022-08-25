import pandas as pd                                                                                                                                                                                                                                                             
from sklearn.svm import SVC                                                                                                                                                                                                                                                     
import joblib         
import numpy as np
        
# TO DO:
# Feature comparison (between my version of SV2 and Danny's version)
# Create new classifiers and run on my version of SV2

# Compare my features with Danny's.
# Check each copy number value and make sure that the correspondence is the same for mine and Danny's.
# Train classifiers.
# Run on HG002 features.
    
def open_train(filepath):                                                                                                                                                                                                                                                       
    df = pd.read_csv(filepath, sep = "\t")                                                                                                                                                                                                                                      
    return df[df["covr"] < 5.0] # Do not train on outliers                                                                                                                                                                                                                      


# Generate weights for biallelic SVs
# Genotypes must be ordered as [REF, HET, HOM]
def make_biallelic_weights(coverage, genotype, heterozygous_expected_coverage, homozygous_expected_coverage):
    buffer = 0.01
    expected_coverage = 1.0
    if genotype == 1: expected_coverage = heterozygous_expected_coverage
    elif genotype == 0: expected_coverage = homozygous_expected_coverage
    return 1.0/(abs(expected_coverage - coverage) + buffer)

# For male sex chromosomes: Danny creates a weight for each sample as a function of the sample's coverage and genotype
# Genotypes must be ordered as [REF, ALT]
def make_allele_weights(coverage, genotype, alternative_allele_expected_coverage):
    buffer = 0.01
    expected_coverage = 1.0
    # Alternate genotype
    if genotype == 0: expected_coverage = alternative_allele_expected_coverage
    return 1.0/(abs(expected_coverage - coverage) + buffer)
    
def make_snv_weights(coverage, genotype, heterozygous_expected_coverage, homozygous_expected_coverage, heterozygous_allele_ratio, median_heterozygous_ratio_copy_number_2, median_heterozygous_ratio_copy_number_3, median_heterozygous_ratio_copy_number_4):
    buffer = 0.01
    if genotype == 2:
        return 1.0 / np.sqrt(( ((1 - coverage)**2) + ((median_heterozygous_ratio_copy_number_2 - heterozygous_allele_ratio)**2) + buffer ) )
    elif genotype == 1:
        return 1.0 / np.sqrt(( ((heterozygous_expected_coverage - coverage)**2) + ((median_heterozygous_ratio_copy_number_3 - heterozygous_allele_ratio)**2) + buffer ) )
    elif genotype == 0:
        return max(1.0 / np.sqrt(( ((homozygous_expected_coverage - coverage)**2) + ((median_heterozygous_ratio_copy_number_2 - heterozygous_allele_ratio)**2) + buffer ) ),
                   1.0 / np.sqrt(( ((homozygous_expected_coverage - coverage)**2) + ((median_heterozygous_ratio_copy_number_4 - heterozygous_allele_ratio)**2) + buffer ) ) )
                   
MEM = 24000                                                                                                                                                                                                                                                                     

# Deletion > 1000 bp classifier     
df_highcov_del_gt1kb = open_train("data/training_set_files/1kgp_highcov_del_gt1kb.txt")[["covr", "dpe_rat", "sr_rat", "copy_number"]]    
sample_weights = df_highcov_del_gt1kb.apply(lambda sample: make_biallelic_weights(sample["covr"], sample["copy_number"], 0.5, 0.0), axis = 1)
clf_del_gt1kb = SVC(kernel = "rbf" , probability = True, random_state = 42, C = 0.01, gamma = 10, class_weight = "balanced", cache_size = MEM)
clf_del_gt1kb.fit(df_highcov_del_gt1kb.values[:, 0:-1], df_highcov_del_gt1kb.values[:, -1], sample_weight = sample_weights)                                                          
joblib.dump(clf_del_gt1kb, "new_classifiers/clf_del_gt1kb.pkl") 
"""
# Deletion <= 1000 Classifier             
df_highcov_del_lt1kb = open_train("data/training_set_files/1kgp_highcov_del_lt1kb.txt")[["covr", "dpe_rat", "sr_rat", "copy_number"]]
sample_weights = df_highcov_del_lt1kb.apply(lambda sample: make_biallelic_weights(sample["covr"], sample["copy_number"], 0.5, 0.0), axis = 1)
clf_del_lt1kb = SVC(kernel = "rbf" ,probability = True, random_state = 42, C = 1000, gamma = 1, class_weight = "balanced", cache_size = MEM)
clf_del_lt1kb.fit(df_highcov_del_lt1kb.values[:, 0:-1], df_highcov_del_lt1kb.values[:, -1], sample_weight = sample_weights)
joblib.dump(clf_del_lt1kb, "clf_del_lt1kb.pkl")   

# Deletion male sex chromosome classifier
df_del_malesexchrom = open_train("data/training_set_files/1kgp_highcov_del_malesexchrom.txt")[["covr", "dpe_rat", "sr_rat", "copy_number"]]            
sample_weights = df_del_malesexchrom.apply(lambda sample: make_allele_weights(sample["covr"], sample["copy_number"], 0.0), axis = 1)
clf_del_malesexchrom = SVC(kernel = "rbf", probability = True, random_state = 42, C = 1, gamma = 1, class_weight = "balanced", cache_size = MEM)           
clf_del_malesexchrom.fit(df_del_malesexchrom.values[:, 0:-1], df_del_malesexchrom.values[:, -1], sample_weight = sample_weights)
joblib.dump(clf_del_malesexchrom, "clf_del_malesexchrom.pkl")

# Duplication breakpoint classifier
# df_dup_breakpoint["copy_number"].unique() -> array([2, 3, 4])
df_dup_breakpoint = open_train("data/training_set_files/1kgp_lowcov_dup_breakpoint.txt")[["covr", "dpe_rat", "sr_rat", "copy_number"]]
# Danny defines a function called "dup_convert" (for duplications) where he converts "copy_number" 3 to 1, and 4 to 0 (and default is 2)
df_dup_breakpoint["copy_number"] = df_dup_breakpoint["copy_number"].mask(df_dup_breakpoint["copy_number"] == 3, 1)
df_dup_breakpoint["copy_number"] = df_dup_breakpoint["copy_number"].mask(df_dup_breakpoint["copy_number"] == 4, 0)
class_weights = {0: 2, 1: 0.5, 2: 1}    
sample_weights = df_dup_breakpoint.apply(lambda sample: make_biallelic_weights(sample["covr"], sample["copy_number"], 1.5, 2.0), axis = 1)
clf_dup_breakpoint = SVC(kernel = "rbf", probability = True, random_state = 42, C = 1, gamma = 10, class_weight = class_weights, cache_size = MEM)
clf_dup_breakpoint.fit(df_dup_breakpoint.values[:, 0:-1], df_dup_breakpoint.values[:, -1], sample_weight = sample_weights)                                            
joblib.dump(clf_dup_breakpoint, "clf_dup_breakpoint.pkl") 

# Duplication HAR classifier
df_dup_har = open_train("data/training_set_files/1kgp_highcov_dup_har.txt")[["covr", "HET_ratio", "copy_number"]]
df_dup_har["copy_number"] = df_dup_har["copy_number"].mask(df_dup_har["copy_number"] == 3, 1)
df_dup_har["copy_number"] = df_dup_har["copy_number"].mask(df_dup_har["copy_number"] == 4, 0)
median_heterozygous_ratio_copy_number_2 = np.median(df_dup_har[df_dup_har["copy_number"] == 2]["HET_ratio"])
median_heterozygous_ratio_copy_number_3 = np.median(df_dup_har[df_dup_har["copy_number"] == 1]["HET_ratio"])
median_heterozygous_ratio_copy_number_4 = median_heterozygous_ratio_copy_number_2/2.0
sample_weights = df_dup_har.apply(lambda sample: make_snv_weights(sample["covr"], sample["copy_number"], 1.5, 2.0, sample["HET_ratio"], median_heterozygous_ratio_copy_number_2, median_heterozygous_ratio_copy_number_3, median_heterozygous_ratio_copy_number_4), axis = 1)
clf_dup_har = SVC(kernel = "rbf", probability = True, random_state = 42, C = 0.01, gamma = 5000, class_weight = "balanced", cache_size
clf_dup_har.fit(df_dup_har.values[:, 0:-1], df_dup_har.values[:, -1])
joblib.dump(clf_dup_har, "clf_dup_har.pkl")

# Duplication male sex chromosome classifier
# df_dup_malesexchrom["copy_number"].unique() -> array([1, 3])
df_dup_malesexchrom = open_train("data/training_set_files/1kgp_lowcov_dup_malesexchrom.txt")[["covr", "dpe_rat", "sr_rat", "copy_number"]]                                                                                                                                      
# Danny defines a function called "dup_convert_msc" (for duplications on the male sex chromosomes) where he changes "copy_number" to 0 if it's greater than 1
df_dup_malesexchrom["copy_number"] = df_dup_malesexchrom["copy_number"].mask(df_dup_malesexchrom["copy_number"] > 1, 0)
class_weights = {0: 1, 1: 1}    
sample_weights = df_dup_malesexchrom.apply(lambda sample: make_allele_weights(sample["covr"], sample["copy_number"], 2.0), axis = 1)
clf_dup_malesexchrom = SVC(kernel = "rbf", probability = True, random_state = 42, C = 1000, gamma = 0.01, class_weight = class_weights, cache_size = MEM)                                
clf_dup_malesexchrom.fit(df_dup_malesexchrom.values[:, 0:-1], df_dup_malesexchrom.values[:, -1], sample_weight = sample_weights)
joblib.dump(clf_dup_malesexchrom, "clf_dup_malesexchrom.pkl")
"""
