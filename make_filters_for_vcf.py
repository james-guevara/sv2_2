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
vcf_out = pysam.VariantFile(sys.argv[2], mode = "w", header = vcf_iterator.header)

format_tuple = ("GT", "PE", "SR", "SC", "NS", "HA", "NH", "SQ", "GL")
for record in vcf_iterator:
    new_record = vcf_out.header.new_record(contig = record.chrom, start = record.pos, stop = record.stop, alleles = record.alleles, id = record.id, qual = record.qual, filter = record.filter, info = record.info)
    new_record.samples[0]["GT"] = record.samples[0]["GT"]
    if record.qual == -1: 
        vcf_out.write(new_record) 
        continue 
    # Else, SQ should be in the format field
    SQ = record.samples[0]["SQ"]
    PE = record.samples[0]["PE"]
    SR = record.samples[0]["SR"]
    SVLEN = record.info["SVLEN"]
    SVTYPE = record.info["SVTYPE"]

    if calculate_filter(SVTYPE, SVLEN, PE, SR, SQ): new_record.filter.add("PASS")
    else: new_record.filter.add("FAIL")

    for element in format_tuple: new_record.samples[0][element] = record.samples[0][element]
    vcf_out.write(new_record) 
