import pandas as pd
from sklearn.svm import SVC
import joblib

def open_train(filepath):
    df = pd.read_csv(filepath, sep = "\t")
    return df[df["covr"] < 5.0] # Do not train on outliers

MEM = 12000

# Deletion > 1000 bp classifier
# df_highcov_del_gt1kb = open_train("data/training_set_files/1kgp_highcov_del_gt1kb.txt")[["covr", "dpe_rat", "sr_rat", "copy_number"]]
# clf_del_gt1kb = SVC(kernel = "rbf" , probability = True, random_state = 42, C = 0.01, gamma = 10, class_weight = "balanced", cache_size = MEM)
# clf_del_gt1kb.fit(df_highcov_del_gt1kb.values[:, 0:-1], df_highcov_del_gt1kb.values[:, -1])
# joblib.dump(clf_del_gt1kb, "clf_del_gt1kb.pkl")

# Deletion <= 1000 Classifier
# df_highcov_del_lt1kb = open_train("data/1kgp_highcov_del_lt1kb.txt")[["covr", "dpe_rat", "sr_rat", "copy_number"]]
# clf_del_lt1kb = SVC(kernel = "rbf" ,probability = True, random_state = 42, C = 1000, gamma = 1, class_weight = "balanced", cache_size = MEM)
# clf_del_lt1kb.fit(df_highcov_del_lt1kb.values[:, 0:-1], df_highcov_del_lt1kb.values[:, -1])
# joblib.dump(clf_del_lt1kb, "clf_del_lt1kb.pkl")

# Deletion male sex chromosome classifier
# df_del_malesexchrom = open_train("data/training_set_files/1kgp_highcov_del_malesexchrom.txt")[["covr", "dpe_rat", "sr_rat", "copy_number"]]
# clf_del_malesexchrom = SVC(kernel = "rbf", probability = True, random_state = 42, C = 1, gamma = 1, class_weight = "balanced", cache_size = MEM)
# clf_del_malesexchrom.fit(df_del_malesexchrom.values[:, 0:-1], df_del_malesexchrom.values[:, -1])
# joblib.dump(clf_del_malesexchrom, "clf_del_malesexchrom.pkl")

# Duplication breakpoint classifier
# df_dup_breakpoint = open_train("data/training_set_files/1kgp_lowcov_dup_breakpoint.txt")[["covr", "dpe_rat", "sr_rat", "copy_number"]] 
# clf_dup_breakpoint = SVC(kernel = "rbf", probability = True, random_state = 42, C = 1, gamma = 10, class_weight = "balanced", cache_size = MEM)
# clf_dup_breakpoint.fit(df_dup_breakpoint.values[:, 0:-1], df_dup_breakpoint.values[:, -1])
# joblib.dump(clf_dup_breakpoint, "clf_dup_breakpoint.pkl")

# Duplication HAR classifier
# df_dup_har = open_train("data/training_set_files/1kgp_highcov_dup_har.txt")[["covr", "HET_ratio", "copy_number"]]
# clf_dup_har = SVC(kernel = "rbf", probability = True, random_state = 42, C = 0.01, gamma = 5000, class_weight = "balanced", cache_size = MEM)
# clf_dup_har.fit(df_dup_har.values[:, 0:-1], df_dup_har.values[:, -1])
# joblib.dump(clf_dup_har, "clf_dup_har.pkl")

# Duplication male sex chromosome classifier
df_dup_malesexchrom = open_train("data/training_set_files/1kgp_lowcov_dup_malesexchrom.txt")[["covr", "dpe_rat", "sr_rat", "copy_number"]] 
clf_dup_malesexchrom = SVC(kernel = "rbf", probability = True, random_state = 42, C = 1000, gamma = 0.01, class_weight = "balanced", cache_size = MEM)
clf_dup_malesexchrom.fit(df_dup_malesexchrom.values[:, 0:-1], df_dup_malesexchrom.values[:, -1])
joblib.dump(clf_dup_malesexchrom, "clf_dup_malesexchrom.pkl")
