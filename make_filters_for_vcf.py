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

vcf_iterator = pysam.VariantFile(sys.argv[1])

number_of_samples = len(vcf_iterator.header.samples)
SQs = np.zeros(shape = (number_of_samples,))

vcf_out = pysam.VariantFile(sys.argv[2], mode = "w", header = vcf_iterator.header)

# format_tuple = ("GT", "PE", "SR", "SC", "NS", "HA", "NH", "SQ", "GL")
format_tuple = ("GT", "PE", "SR", "SC", "NS", "HA", "NH", "SQ")
for record in vcf_iterator:
    new_record = vcf_out.header.new_record(contig = record.chrom, start = record.pos, stop = record.stop, alleles = record.alleles, id = record.id, qual = record.qual, filter = record.filter, info = record.info)
    for i in range(len(record.samples)): new_record.samples[0]["GT"] = record.samples[0]["GT"]
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

    SQ_median = np.nanmedian(SQs)

    # If we're dealing with single sample VCF, then this part makes sense. But otherwise, we can't use it...
    # SQ = record.samples[0]["SQ"]
    # PE = record.samples[0]["PE"]
    # SR = record.samples[0]["SR"]
    
    SVLEN = record.info["SVLEN"]
    SVTYPE = record.info["SVTYPE"]

    if calculate_filter(SVTYPE, SVLEN, PE, SR, SQ_median): new_record.filter.add("PASS")
    else: new_record.filter.add("FAIL")

    new_record.qual = SQ_median

    vcf_out.write(new_record) 
