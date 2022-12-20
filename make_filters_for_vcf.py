import numpy as np
import pysam
import sys

def gt1kb_del_filter(svlen: int, discordant_ratio: float, split_ratio: float, SQ: float):
    if discordant_ratio > 0 or split_ratio > 0: return SQ >= 8
    else:
        if svlen < 3000: return False
        if svlen < 5000: return SQ >= 20
        return SQ >= 18

def lt1kb_del_filter(svlen: int, discordant_ratio: float, split_ratio: float, SQ: float):
    if discordant_ratio > 0 or split_ratio > 0: return SQ >= 8
    else: return SQ >= 18
     

def del_filter(svlen: int, discordant_ratio: float, split_ratio: float, SQ: float):
    if svlen > 1000: return gt1kb_del_filter(svlen, discordant_ratio, split_ratio, SQ)
    else: return lt1kb_del_filter(svlen, discordant_ratio, split_ratio, SQ)


def dup_filter(svlen: int, discordant_ratio: float, split_ratio: float, SQ: float):
    if discordant_ratio > 0 or split_ratio > 0:
        if svlen <= 1000: return SQ >= 12
        return SQ >= 10
    else: return SQ >= 8

def calculate_filter(svtype: str, svlen: int, discordant_ratio: float, split_ratio: float, SQ: float):
    if svtype == "DEL": return del_filter(svlen, discordant_ratio, split_ratio, SQ)
    else: return dup_filter(svlen, discordant_ratio, split_ratio, SQ)

# DENOVO_FILTER functions
def denovo_dup_filter(svlen: int, discordant_ratio: float, split_ratio: float, SQ: float):
    if discordant_ratio > 0 or split_ratio > 0: return SQ >= 11
    else: return SQ >= 10
    

def denovo_del_filter(svlen: int, discordant_ratio: float, split_ratio: float, SQ: float):
    if discordant_ratio > 0 or split_ratio > 0: return SQ >= 12
    else:
        if svlen <= 1000: return SQ >= 24
        return SQ >= 20

def calculate_denovo_filter(svtype: str, svlen: int, discordant_ratio: float, split_ratio: float, SQ: float):
    if svtype == "DEL": return denovo_del_filter(svlen, discordant_ratio, split_ratio, SQ)
    else: return denovo_dup_filter(svlen, discordant_ratio, split_ratio, SQ)


vcf_iterator = pysam.VariantFile(sys.argv[1])

vcf_out = pysam.VariantFile(sys.argv[2], mode = "w", header = vcf_iterator.header)
vcf_out.header.add_meta(key = "INFO", items = [ ("ID", "ALT_SQ_MEDIAN"), ("Number", 1), ("Type", "Float"), ("Description", "Median of SQs from ALT genotypes") ] )
vcf_out.header.add_meta(key = "INFO", items = [ ("ID", "ALT_STANDARD_FILTER"), ("Number", 1), ("Type", "Float"), ("Description", "Standard filter using median of SQs from ALT genotypes") ] )
vcf_out.header.add_meta(key = "INFO", items = [ ("ID", "ALT_DENOVO_FILTER"), ("Number", 1), ("Type", "Float"), ("Description", "De novo filter using median of SQs from ALT genotypes") ] )

# format_tuple = ("GT", "PE", "SR", "SC", "NS", "HA", "NH", "SQ", "GL")
format_tuple = ("GT", "PE", "SR", "SC", "NS", "HA", "NH", "SQ")
SQs = np.zeros(shape = (len(vcf_iterator.header.samples),))
ALT_SQs = np.nans(shape = (len(vcf_iterator.header.samples),))
for record in vcf_iterator:
    new_record = vcf_out.header.new_record(contig = record.chrom, start = record.start, stop = record.stop, alleles = record.alleles, id = record.id, qual = record.qual, filter = record.filter, info = record.info)
    for i in range(len(record.samples)): new_record.samples[i]["GT"] = record.samples[i]["GT"]
    if record.qual == -1: 
        vcf_out.write(new_record) 
        continue 

    PE = 0.
    SR = 0.
    # Else, SQ should be in the format field (and we use it to get median SQ, which we set as the new QUAL and use for filtering)
    for i in range(len(record.samples)):
        SQs[i] = record.samples[i]["SQ"]
        if record.samples[i]["PE"] not in (".", None): PE += record.samples[i]["PE"]
        if record.samples[i]["SR"] not in (".", None): SR += record.samples[i]["SR"]
        for element in format_tuple: new_record.samples[i][element] = record.samples[i][element]
        if record.samples[i]["GL"] == ".": new_record.samples[i]["GL"] = (0, 0, 0)
        if record.samples[i] == (0, 1):
            ALT_SQs[i] = record.samples[i]["SQ"]

    SQ_median = np.nanmedian(SQs)
    ALT_SQ_median = np.nanmedian(ALT_SQs)
    new_record.qual = SQ_median
    
    if calculate_filter(record.info["SVTYPE"], record.info["SVLEN"], PE, SR, SQ_median): new_record.filter.add("PASS")
    else: new_record.filter.add("FAIL")
    # If variant fails the standard filter, then it also fails DENOVO_FILTER 
    new_record.info["DENOVO_FILTER"] = "FAIL"
    if "FAIL" not in new_record.filter and calculate_denovo_filter(record.info["SVTYPE"], record.info["SVLEN"], PE, SR, SQ_median): new_record.info["DENOVO_FILTER"] = "PASS"


    record.info["ALT_SQ_MEDIAN"] = ALT_SQ_median
    new_record.info["ALT_STANDARD_FILTER"] = "FAIL"
    new_record.info["ALT_DENOVO_FILTER"] = "FAIL"
    if calculate_filter(record.info["SVTYPE"], record.info["SVLEN"], PE, SR, ALT_SQ_median): new_record.info["ALT_STANDARD_FILTER"] = "PASS"
    if new_record.info["ALT_STANDARD_FILTER"] != "FAIL" and calculate_denovo_filter(record.info["SVTYPE"], record.info["SVLEN"], PE, SR, ALT_SQ_median): new_record.info["ALT_DENOVO_FILTER"] = "PASS"

    vcf_out.write(new_record) 
