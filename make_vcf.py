# To do:
# Fix positional information (converting from BED/table to VCF...)
# I'm pretty sure the table is in BED format, meaning that we can take the CHROM, START, END, use pybedtools to wrap that as BED formatted data, and then intersect with whatever we want.
#   When we intersect, we can specify new columns/labels and how much each SV overlaps with each type of genomic region.
# Add sample name. (Should add to table first, and then get that...)
import pysam
import pysam.bcftools
import argparse

def make_vcf(sample_name, reference_fasta, genotype_predictions_table, output_vcf_filepath = None):
    vcf_header = pysam.VariantHeader()
    vcf_header.add_sample(sample_name)
    # FILTER fields
    vcf_header.add_meta("FILTER", items = [ ("ID", "RF"), ("Description", "Variant filed filter due to low RF") ] )
    # INFO fields
    vcf_header.add_meta("INFO", items = [ ("ID", "END"), ("Number", 1), ("Type", "Integer"), ("Description", "End position of structural variant") ] )
    vcf_header.add_meta("INFO", items = [ ("ID", "SVTYPE"), ("Number", 1), ("Type", "String"), ("Description", "Type of structural variant") ] )
    vcf_header.add_meta("INFO", items = [ ("ID", "SVLEN"), ("Number", 1), ("Type", "Integer"), ("Description", "Length of structural variant") ] )
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

    vcf_filename = "{}_sv2.vcf".format(sample_name)
    if output_vcf_filepath: vcf_filename = output_vcf_filepath
    vcf = pysam.VariantFile(vcf_filename, "w", header = vcf_header)

    # with open("genotyping_preds2022-06-03_23.36.29.tsv", "r") as f:
    with open(genotype_predictions_table, "r") as f:
        header = f.readline().rstrip()
        for line in f:
            linesplit = line.rstrip().split("\t")
            chrom, start, end, type_, size, coverage, coverage_GCcorrected, discordant_ratio, split_ratio, snv_coverage, heterozygous_allele_ratio, snvs, het_snvs, ALT_GENOTYPE_LIKELIHOOD, REF_QUAL, ALT_QUAL, HOM_GENOTYPE_LIKELIHOOD, HET_GENOTYPE_LIKELIHOOD, REF_GENOTYPE_LIKELIHOOD, GEN = linesplit[0], linesplit[1], linesplit[2], linesplit[3], linesplit[4], linesplit[5], linesplit[6], linesplit[7], linesplit[8], linesplit[9], linesplit[10], linesplit[11], linesplit[12], linesplit[13], linesplit[14], linesplit[15], linesplit[16], linesplit[17], linesplit[18], linesplit[19]
            if snv_coverage == "": snv_coverage = float("nan") 
            if heterozygous_allele_ratio == "": heterozygous_allele_ratio = float("nan") 

            record = vcf.new_record(contig = chrom, start = int(start), stop = int(end), alleles = ("N", "<{}>".format(type_)))
    
            record.info["SVTYPE"] = type_
            record.info["SVLEN"] = int(size)
    
            record.samples[sample_name]["GT"] = ( GEN.split("/")[0], GEN.split("/")[1] )
            record.samples[sample_name]["PE"] = round(float(discordant_ratio), 1)
            record.samples[sample_name]["SR"] = round(float(split_ratio), 1)
            record.samples[sample_name]["SC"] = round(float(snv_coverage), 1)
            record.samples[sample_name]["NS"] = int(float(snvs))
            record.samples[sample_name]["HA"] = round(float(heterozygous_allele_ratio), 1)
            record.samples[sample_name]["NH"] = int(float(het_snvs))
            record.samples[sample_name]["SQ"] = round(max(float(REF_QUAL), float(ALT_QUAL)), 1)
            record.samples[sample_name]["GL"] = (round(float(REF_GENOTYPE_LIKELIHOOD), 1), round(float(HET_GENOTYPE_LIKELIHOOD), 1), round(float(HOM_GENOTYPE_LIKELIHOOD), 1))
            vcf.write(record)
    
    vcf.close()
    pysam.bcftools.sort("-o", "{}.gz".format(vcf_filename), vcf_filename, catch_stdout = False)

if __name__ == "__main__":
    # Argument parsing (when this script is run by itself)
    parser = argparse.ArgumentParser(description = "Make SV2 VCF")
    parser.add_argument("--genotype_predictions_table", help = "Genotype predictions output .tsv file from classify.py step")
    parser.add_argument("--reference_fasta", help = "Reference fasta file input")
    parser.add_argument("--sample_name", help = "Sample name")
    parser.add_argument("--output_vcf", help = "Output VCF filepath (optional)")
    args = parser.parse_args()

    make_vcf(args.sample_name, args.reference_fasta, args.genotype_predictions_table, args.output_vcf)
