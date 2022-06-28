import pandas as pd
from sklearn.svm import SVC
import joblib

MEM = 12000
# Deletion > 1000 bp classifier
highcov_del_gt1kb  = "1kgp_highcov_del_gt1kb.txt"
df_highcov_del_gt1kb = pd.read_csv(highcov_del_gt1kb, sep = "\t")
df_highcov_del_gt1kb_subset = df_highcov_del_gt1kb[["covr", "dpe_rat", "sr_rat"]]
features_highcov_del_gt1kb = df_highcov_del_gt1kb_subset.values
truth_highcov_del_gt1kb = df_highcov_del_gt1kb[["copy_number"]]

clf_del_gt1kb = SVC(kernel = "rbf" , probability = True, random_state = 42, C = 0.01, gamma = 10, class_weight = "balanced", cache_size = MEM)
clf_del_gt1kb.fit(features_highcov_del_gt1kb, truth_highcov_del_gt1kb)
joblib.dump(clf_del_gt1kb, "clf_del_gt1kb.pkl")

# Deletion <= 1000 Classifier
highcov_del_lt1kb  = "1kgp_highcov_del_lt1kb.txt"
df_highcov_del_lt1kb = pd.read_csv(highcov_del_lt1kb, sep = "\t")
df_highcov_del_lt1kb_subset = df_highcov_del_lt1kb[["covr", "dpe_rat", "sr_rat"]]
features_highcov_del_lt1kb = df_highcov_del_lt1kb_subset.values
truth_highcov_del_lt1kb = df_highcov_del_lt1kb[["copy_number"]]

clf_del_lt1kb = SVC(kernel = "rbf" ,probability = True, random_state = 42, C = 1000, gamma = 1, class_weight = "balanced", cache_size = MEM)
clf_del_lt1kb.fit(features_highcov_del_lt1kb, truth_highcov_del_lt1kb)
joblib.dump(clf_del_lt1kb, "clf_del_lt1kb.pkl")

# Deletion male sex chromosome classifier
del_malesexchrom= "1kgp_highcov_del_malesexchrom.txt"
df_del_malesexchrom = pd.read_csv(del_malesexchrom, sep = "\t")
df_del_malesexchrom_subset = df_del_malesexchrom[["covr", "dpe_rat", "sr_rat"]]
feaures_del_malesexchrom = df_del_malesexchrom_subset.values
truth_del_malesexchrom = df_del_malesexchrom[["copy_number"]]

clf_del_malesexchrom = SVC(kernel = "rbf", probability = True, random_state = 42, C = 1, gamma = 1, class_weight = "balanced", cache_size = MEM)
clf_del_malesexchrom.fit(feaures_del_malesexchrom, truth_del_malesexchrom)
joblib.dump(clf_del_malesexchrom, "clf_del_malesexchrom.pkl")

# Duplication breakpoint classifier
dup_breakpoint = "1kgp_lowcov_dup_breakpoint.txt"
df_dup_breakpoint = pd.read_csv(dup_breakpoint, sep = "\t")
df_dup_breakpoint_subset = df_dup_breakpoint[["covr", "dpe_rat", "sr_rat"]]
features_dup_breakpoint = df_dup_breakpoint_subset.values
truth_dup_breakpoint = df_dup_breakpoint[["copy_number"]]

clf_dup_breakpoint = SVC(kernel = "rbf", probability = True, random_state = 42, C = 1, gamma = 10, class_weight = "balanced", cache_size = MEM)
clf_dup_breakpoint.fit(features_dup_breakpoint,truth_dup_breakpoint)
joblib.dump(clf_dup_breakpoint, "clf_dup_breakpoint.pkl")

# Duplication male sex chromosome classifier
dup_malesexchrom = "1kgp_lowcov_dup_malesexchrom.txt"
df_dup_malesexchrom = pd.read_csv(dup_malesexchrom, sep = "\t")
df_dup_malesexchrom_subset = df_dup_malesexchrom[["covr", "dpe_rat", "sr_rat"]]
features_df_dup_malesexchrom = df_dup_malesexchrom_subset.values
truth_dup_malesexchrom = df_dup_malesexchrom[["copy_number"]]

clf_dup_malesexchrom = SVC(kernel = "rbf", probability = True, random_state = 42, C = 1000, gamma = 0.01, class_weight = "balanced", cache_size = MEM)
clf_dup_malesexchrom.fit(features_df_dup_malesexchrom, truth_dup_malesexchrom)
joblib.dump(clf_dup_malesexchrom, "clf_dup_malesexchrom.pkl")

# Duplication HAR classifier
dup_har = "1kgp_highcov_dup_har.txt"
df_dup_har = pd.read_csv(dup_har, sep = "\t")
df_dup_har_subset = df_dup_har[["covr","HET_ratio"]]
features_dup_har = df_dup_har_subset.values
truth_dup_har = df_dup_har[["copy_number"]]

clf_dup_har = SVC(kernel = "rbf", probability = True, random_state = 42, C = 0.01, gamma = 5000, class_weight = "balanced", cache_size = MEM)
clf_dup_har.fit(features_dup_har, truth_dup_har)
joblib.dump(clf_dup_har, "clf_dup_har.pkl")
