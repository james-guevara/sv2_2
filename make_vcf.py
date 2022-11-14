# To do:
# Fix positional information (converting from BED/table to VCF...)
# I'm pretty sure the table is in BED format, meaning that we can take the CHROM, START, END, use pybedtools to wrap that as BED formatted data, and then intersect with whatever we want.
#   When we intersect, we can specify new columns/labels and how much each SV overlaps with each type of genomic region.
import argparse
import numpy as np
import pysam
import pysam.bcftools
from time import gmtime, strftime

def gt1kb_del_filter(svlen: int, discordant_ratio: float, split_ratio: float, SQ: float):
    if discordant_ratio > 0 and split_ratio > 0: return SQ >= 8
    else:
        if svlen < 3000: return False
        if svlen < 5000: return SQ >= 20
        return SQ >= 18

def lt1kb_del_filter(svlen: int, discordant_ratio: float, split_ratio: float, SQ: float):
    if discordant_ratio > 0 and split_ratio > 0: return SQ >= 8
    else: return SQ >= 18
     

def del_filter(svlen: int, discordant_ratio: float, split_ratio: float, SQ: float):
    if svlen > 1000: return gt1kb_del_filter(svlen, discordant_ratio, split_ratio)
    else: return lt1kb_del_filter(svlen, discordant_ratio, split_ratio)

def dup_filter(svlen: int, discordant_ratio: float, split_ratio: float, SQ: float):


def calculate_filter(svtype: str, svlen: int, discordant_ratio: float, split_ratio: float, SQ: float):
    if svtype == "DEL": return del_filter(svlen, discordant_ratio, split_ratio)
    else: return dup_filter(svlen, discordant_ratio, split_ratio)


def make_vcf(sample_name, reference_fasta, genotype_predictions_table, output_folder, sv2_command, current_time, sv_bed_list_unfiltered):
    vcf_header = pysam.VariantHeader()
    vcf_header.add_sample(sample_name)
    vcf_header.add_line("##SV2_CMD='{}'".format(sv2_command))
    # FILTER fields
    vcf_header.add_meta("FILTER", items = [ ("ID", "FAIL"), ("Description", "Variant failed standard filters") ] )
    vcf_header.add_meta("FILTER", items = [ ("ID", "PASS"), ("Description", "Variant passed standard filters") ] )
    vcf_header.add_meta("FILTER", items = [ ("ID", "UNGENOTYPED"), ("Description", "Variant was not genotyped (perhaps malformed, or in exclude region)") ] )
    # INFO fields
    vcf_header.add_meta("INFO", items = [ ("ID", "END"), ("Number", 1), ("Type", "Integer"), ("Description", "End position of structural variant") ] )
    vcf_header.add_meta("INFO", items = [ ("ID", "SVTYPE"), ("Number", 1), ("Type", "String"), ("Description", "Type of structural variant") ] )
    vcf_header.add_meta("INFO", items = [ ("ID", "SVLEN"), ("Number", 1), ("Type", "Integer"), ("Description", "Length of structural variant") ] )
    vcf_header.add_meta("INFO", items = [ ("ID", "DENOVO_FILTER"), ("Number", 1), ("Type", "String"), ("Description", "Stringent filter status, recommended for de novo mutation discovery") ] )
    # FORMAT fields
    vcf_header.add_meta("FORMAT", items = [ ("ID", "GT"), ("Number", 1), ("Type", "String"), ("Description", "Genotype") ] )
    vcf_header.add_meta("FORMAT", items = [ ("ID", "PE"), ("Number", 1), ("Type", "Float"), ("Description", "Normalized discordant paired-end count") ] )
    vcf_header.add_meta("FORMAT", items = [ ("ID", "SR"), ("Number", 1), ("Type", "Float"), ("Description", "Normalized split-read count") ] )
    vcf_header.add_meta("FORMAT", items = [ ("ID", "SC"), ("Number", 1), ("Type", "Float"), ("Description", "SNV normalized coverage") ] )
    vcf_header.add_meta("FORMAT", items = [ ("ID", "NS"), ("Number", 1), ("Type", "Integer"), ("Description", "Number of SNVs within locus") ] )
    vcf_header.add_meta("FORMAT", items = [ ("ID", "HA"), ("Number", 1), ("Type", "Float"), ("Description", "Heterozygous allele ratio") ] )
    vcf_header.add_meta("FORMAT", items = [ ("ID", "NH"), ("Number", 1), ("Type", "Integer"), ("Description", "Number of heterozygous SNVs") ] )
    vcf_header.add_meta("FORMAT", items = [ ("ID", "SQ"), ("Number", 1), ("Type", "Float"), ("Description", "Phred-scaled genotype likelihood") ] )
    vcf_header.add_meta("FORMAT", items = [ ("ID", "GL"), ("Number", "G"), ("Type", "Float"), ("Description", "Phred-scaled genotype likelihoods in the order, REF:(0/0), HET:(0/1), HOM:(1/1)") ] )
    # Contigs
    # fasta = pysam.FastaFile(filename = "/expanse/lustre/projects/ddp195/j3guevar/resources/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa")
    fasta = pysam.FastaFile(filename = reference_fasta)
    for contig, contig_length in zip(fasta.references, fasta.lengths):
        vcf_header.contigs.add(contig, length = contig_length)

    vcf_filename = "{}/{}_sv2.vcf".format(output_folder, sample_name)
    vcf = pysam.VariantFile(vcf_filename, "w", header = vcf_header)

    def make_column_indices_table(header_list):
        index_dictionary = {}
        for index, name in enumerate(header_list):
            index_dictionary[name] = index
        return index_dictionary

    # Output all the variants into the VCF, even the ones we don't genotype. This requires the original SV file
    genotype_records = set()
    with open(genotype_predictions_table, "r") as f:
        header_list = f.readline().rstrip().split("\t")
        index_dictionary = make_column_indices_table(header_list)
        for line in f:
            linesplit = line.rstrip().split("\t")
            # SV identifier
            chrom = linesplit[index_dictionary["chrom"]] 
            start = int(linesplit[index_dictionary["start"]])
            end   = int(linesplit[index_dictionary["end"]])
            type_ = linesplit[index_dictionary["type"]] 
            size  = int(linesplit[index_dictionary["size"]])

            # print("{}\t{}\t{}\t{}".format(chrom, start, end, type_))

            genotype_records.add((chrom, start, end, type_))

            # Coverage features
            coverage             = linesplit[index_dictionary["coverage"]] 
            coverage_GCcorrected = linesplit[index_dictionary["coverage_GCcorrected"]] 

            # Read-level features
            discordant_ratio = float(linesplit[index_dictionary["discordant_ratio"]])
            split_ratio      = float(linesplit[index_dictionary["split_ratio"]]) 

            # SNV features
            snv_coverage              = linesplit[index_dictionary["snv_coverage"]]
            if snv_coverage == "": snv_coverage = float("nan")
            elif snv_coverage == "inf": snv_coverage = float("inf")
            else: snv_coverage = float(snv_coverage)
            heterozygous_allele_ratio = linesplit[index_dictionary["heterozygous_allele_ratio"]]
            if heterozygous_allele_ratio == "": heterozygous_allele_ratio = float("nan")
            snvs                      = int(linesplit[index_dictionary["snvs"]])
            het_snvs                  = int(linesplit[index_dictionary["het_snvs"]]) 

            # Genotype likelihoods and quality scores
            REF_GENOTYPE_LIKELIHOOD = float(linesplit[index_dictionary["REF_GENOTYPE_LIKELIHOOD"]])
            HET_GENOTYPE_LIKELIHOOD = float(linesplit[index_dictionary["HET_GENOTYPE_LIKELIHOOD"]]) 
            HOM_GENOTYPE_LIKELIHOOD = float(linesplit[index_dictionary["HOM_GENOTYPE_LIKELIHOOD"]])
            ALT_GENOTYPE_LIKELIHOOD = float(linesplit[index_dictionary["ALT_GENOTYPE_LIKELIHOOD"]])
            REF_QUAL = float(linesplit[index_dictionary["REF_QUAL"]])
            ALT_QUAL = float(linesplit[index_dictionary["ALT_QUAL"]])

            # Classifier name
            Classifier = linesplit[index_dictionary["Classifier"]]

            # Genotype
            GEN = linesplit[index_dictionary["GEN"]]

            # Creating variant record in VCF
            record = vcf.new_record(contig = chrom, start = start, stop = end, alleles = ("N", "<{}>".format(type_)))

            SQ = round(max(REF_QUAL, ALT_QUAL), 2) # SQ and QUAL score

            record.qual = SQ

    
            record.info["SVTYPE"] = type_
            record.info["SVLEN"] = size

            record.samples[sample_name]["GT"] = (int(GEN.split("/")[0]), int(GEN.split("/")[1]))
            record.samples[sample_name]["PE"] = round(discordant_ratio, 2)
            record.samples[sample_name]["SR"] = round(split_ratio, 2)
            record.samples[sample_name]["SC"] = round(snv_coverage, 2)
            record.samples[sample_name]["NS"] = snvs
            record.samples[sample_name]["HA"] = round(float(heterozygous_allele_ratio), 2)
            record.samples[sample_name]["NH"] = het_snvs
            record.samples[sample_name]["SQ"] = SQ 
            record.samples[sample_name]["GL"] = (round(REF_GENOTYPE_LIKELIHOOD, 2), round(HET_GENOTYPE_LIKELIHOOD, 2), round(HOM_GENOTYPE_LIKELIHOOD, 2))


            # Danny gets the median of the alternate allele likelihoods and puts that value in a variable called "non". I'll put it in a variable called "MEDIAN_ALTERNATE_LIKELIHOOD". (He also gets the median of the reference likelihood... why? It should be only 1 value.)
            MEDIAN_ALTERNATE_LIKELIHOOD = np.median([HOM_GENOTYPE_LIKELIHOOD, HET_GENOTYPE_LIKELIHOOD])
            # Danny also adds the split_ratio and discordant_ratio and puts that output in a variable called "breakpoint_feats" (and then later "feats"...). I'll put it in a variable called "breakpoint_features" 
            breakpoint_features = discordant_ratio + split_ratio
            # De novo filter
            DENOVO_FILTER = "FAIL"
            if Classifier in ("clf_highcov_del_gt1kb", "clf_highcov_del_lt1kb"):
                if breakpoint_features > 0 and MEDIAN_ALTERNATE_LIKELIHOOD >= 12 and REF_GENOTYPE_LIKELIHOOD >= 12: DENOVO_FILTER = "PASS"
                elif breakpoint_features == 0 and 1000 <= size < 3000 and MEDIAN_ALTERNATE_LIKELIHOOD >= 24 and REF_GENOTYPE_LIKELIHOOD >= 20: DENOVO_FILTER = "PASS"
                elif breakpoint_features == 0 and 3000 <= size < 5000 and MEDIAN_ALTERNATE_LIKELIHOOD >= 20 and REF_GENOTYPE_LIKELIHOOD >= 20: DENOVO_FILTER = "PASS"
                elif breakpoint_features == 0 and size >= 5000 and MEDIAN_ALTERNATE_LIKELIHOOD >= 20 and REF_GENOTYPE_LIKELIHOOD >= 18: DENOVO_FILTER = "PASS"
            elif Classifier == "clf_dup_breakpoint":
                if breakpoint_features > 0 and MEDIAN_ALTERNATE_LIKELIHOOD >= 11 and REF_GENOTYPE_LIKELIHOOD >= 11: DENOVO_FILTER = "PASS"
                elif breakpoint_features == 0 and 3000 <= size < 100000 and MEDIAN_ALTERNATE_LIKELIHOOD >= 10 and REF_GENOTYPE_LIKELIHOOD >= 15: DENOVO_FILTER = "PASS"
                elif breakpoint_features == 0 and size >= 100000 and MEDIAN_ALTERNATE_LIKELIHOOD >= 10 and REF_GENOTYPE_LIKELIHOOD >= 13: DENOVO_FILTER = "PASS"
            elif Classifier == "clf_dup_har":
                if size < 5000 and REF_GENOTYPE_LIKELIHOOD >= 18 and MEDIAN_ALTERNATE_LIKELIHOOD >= 18: DENOVO_FILTER = "PASS"
                elif 5000 <= size < 100000 and REF_GENOTYPE_LIKELIHOOD >= 15 and MEDIAN_ALTERNATE_LIKELIHOOD >= 10: DENOVO_FILTER = "PASS"
                elif size >= 100000 and REF_GENOTYPE_LIKELIHOOD >= 13 and MEDIAN_ALTERNATE_LIKELIHOOD >= 10: DENOVO_FILTER = "PASS" 
            elif Classifier == "clf_del_malesexchrom":
                if MEDIAN_ALTERNATE_LIKELIHOOD >= 8 and REF_GENOTYPE_LIKELIHOOD >= 8: DENOVO_FILTER = "PASS" 
            elif Classifier == "clf_dup_malesexchrom":
                if size >= 5000 and REF_GENOTYPE_LIKELIHOOD >= 10 and MEDIAN_ALTERNATE_LIKELIHOOD >= 10: DENOVO_FILTER = "PASS" 
            record.info["DENOVO_FILTER"] = DENOVO_FILTER 

            # Standard filter
            FILTER = "PASS"
            if Classifier in ("clf_highcov_del_gt1kb", "clf_highcov_del_lt1kb"):
                if breakpoint_features > 0 and MEDIAN_ALTERNATE_LIKELIHOOD <= 8: FILTER = "FAIL"
                elif breakpoint_features == 0:
                    if "lt1kb" in Classifier:
                        if MEDIAN_ALTERNATE_LIKELIHOOD < 18: FILTER = "FAIL"
                    elif "gt1kb" in Classifier:
                        if size < 3000: FILTER = "FAIL"
                        elif 3000 <= size < 5000 and MEDIAN_ALTERNATE_LIKELIHOOD < 20: FILTER = "FAIL"
                        elif size >= 5000 and MEDIAN_ALTERNATE_LIKELIHOOD < 18: FILTER = "FAIL"
            elif Classifier == "clf_dup_breakpoint":
                if breakpoint_features > 0:
                    if size <  1000 and MEDIAN_ALTERNATE_LIKELIHOOD < 12: FILTER = "FAIL"
                    if size >= 1000 and MEDIAN_ALTERNATE_LIKELIHOOD < 10: FILTER = "FAIL"
                elif breakpoint_features == 0:
                    if size < 3000: FILTER = "FAIL"
                    elif size >= 3000 and MEDIAN_ALTERNATE_LIKELIHOOD < 12: FILTER = "FAIL"
            elif Classifier == "clf_dup_har":
                if MEDIAN_ALTERNATE_LIKELIHOOD < 8: FILTER = "FAIL" 
                elif size < 3000 and MEDIAN_ALTERNATE_LIKELIHOOD < 18: FILTER = "FAIL"
                elif size >= 3000 and MEDIAN_ALTERNATE_LIKELIHOOD < 13: FILTER = "FAIL"
            elif Classifier == "clf_del_malesexchrom":
                if MEDIAN_ALTERNATE_LIKELIHOOD < 8 and breakpoint_features != 0: FILTER = "FAIL"
                elif MEDIAN_ALTERNATE_LIKELIHOOD < 10 and breakpoint_features == 0 and size <= 1000: FILTER = "FAIL"
                elif size < 1000 and breakpoint_features == 0: FILTER = "FAIL" 
            elif Classifier == "clf_dup_malesexchrom":
                if size < 5000: FILTER = "FAIL"
                elif size >= 5000 and MEDIAN_ALTERNATE_LIKELIHOOD < 10: FILTER = "FAIL"
                elif breakpoint_features == 0: FILTER = "FAIL"

            record.filter.add(FILTER)
            vcf.write(record)

    print("Ungenotyped SVs...")
    for sv in sv_bed_list_unfiltered:
        chrom, start, end, type_ = sv.split("\t")[0], int(sv.split("\t")[1]), int(sv.split("\t")[2]), sv.split("\t")[3]
        sv_tuple = (chrom, start, end, type_)
        if sv_tuple not in genotype_records:
            record = vcf.new_record(contig = chrom, start = start, stop = end, alleles = ("N", "<{}>".format(type_)))
            record.filter.add("UNGENOTYPED")
            record.qual = -1
            record.info["SVTYPE"] = type_
            record.info["SVLEN"] = end - start + 1
            record.samples[sample_name]["GT"] = (".", ".")
            vcf.write(record)

    vcf.close()
    pysam.bcftools.sort("-o", "{}.gz".format(vcf_filename), vcf_filename, catch_stdout = False)

if __name__ == "__main__":
    # Argument parsing (when this script is run by itself)
    parser = argparse.ArgumentParser(description = "Make SV2 VCF")
    parser.add_argument("--genotype_predictions_table", help = "Genotype predictions output .tsv file from classify.py step")
    parser.add_argument("--reference_fasta", help = "Reference fasta file input")
    parser.add_argument("--sample_name", help = "Sample name")
    parser.add_argument("--output_folder", help = "Output folder (optional)")
    args = parser.parse_args()

    current_time = strftime("%Y-%m-%d_%H.%M.%S", gmtime())
    make_vcf(args.sample_name, args.reference_fasta, args.genotype_predictions_table, args.output_folder, current_time)
