# To do:
# Fix positional information (converting from BED/table to VCF...)
# Get rid of duplicated rows in VCF (using 2 different classifiers for duplications I think? One of them is the HAR classifier...)
    # If split-read and paired-end discordant read features are both 0, then we use the HAR classifier (in classify.py, filter the dataframe based on this...),
# I'm pretty sure the table is in BED format, meaning that we can take the CHROM, START, END, use pybedtools to wrap that as BED formatted data, and then intersect with whatever we want.
#   When we intersect, we can specify new columns/labels and how much each SV overlaps with each type of genomic region.
# Add sample name. (Should add to table first, and then get that...)
# Fix warning message: [W::bcf_update_info] INFO/END=0 is smaller than POS at chr1:1 
##  (I suspect this isn't a big deal and is a minor bug in pysam, which might have already been fixed in a more recent version, but I don't know. EDIT: Doesn't look like it.)
# Get reference allele from FASTA? (And use SVLEN to determine how long the reference allele should be... actually maybe we can just set it to "N"; that's what smoove does.) 

# Fix sex chromosome stuff
    # 1. We've trained a sex chromosome classifier (for DUPs and DELS both I think)
    # 2. In the classify.py step, we have to separate the dataframe for if sex == 1: df[df[chrom.contains("X") | chrom.contains("Y")] (basically, separate out this dataframe if these conditions are met)
    # 3. Then, simply run the sex chromosome classifier on this dataframe
    # 4. Later on, we can filter out the pseudo-autosomal region (PAR) because that is autosomal (so use the regular classifiers for that part; it would depend on reference fasta type)

# For DUP variants
    # 1. If split-read and paired-end features are both 0, then we use the HAR classifier. Else, we use the regular classifier.

import pysam
import pysam.bcftools

vcf_header = pysam.VariantHeader()
vcf_header.add_sample("ahstram")
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
fasta = pysam.FastaFile(filename = "/expanse/lustre/projects/ddp195/j3guevar/resources/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa")
for contig, contig_length in zip(fasta.references, fasta.lengths):
    vcf_header.contigs.add(contig, length = contig_length)

vcf = pysam.VariantFile("example.vcf", "w", header = vcf_header)

with open("genotyping_preds2022-06-03_23.36.29.tsv", "r") as f:
    header = f.readline().rstrip()
    for line in f:
        linesplit = line.rstrip().split("\t")
        chrom, start, end, type_, size, coverage, coverage_GCcorrected, discordant_ratio, split_ratio, snv_coverage, heterozygous_allele_ratio, snvs, het_snvs, ALT_GENOTYPE_LIKELIHOOD, REF_QUAL, ALT_QUAL, HOM_GENOTYPE_LIKELIHOOD, HET_GENOTYPE_LIKELIHOOD, REF_GENOTYPE_LIKELIHOOD, GEN = linesplit[0], linesplit[1], linesplit[2], linesplit[3], linesplit[4], linesplit[5], linesplit[6], linesplit[7], linesplit[8], linesplit[9], linesplit[10], linesplit[11], linesplit[12], linesplit[13], linesplit[14], linesplit[15], linesplit[16], linesplit[17], linesplit[18], linesplit[19]
        if snv_coverage == "": snv_coverage = float("nan") 
        if heterozygous_allele_ratio == "": heterozygous_allele_ratio = float("nan") 
        record = vcf.new_record(contig = chrom, start = int(start), stop = int(end), alleles = ("N", "<{}>".format(type_)))

        record.info["SVTYPE"] = type_
        record.info["SVLEN"] = int(size)

        record.samples["ahstram"]["GT"] = ( GEN.split("/")[0], GEN.split("/")[1] )
        record.samples["ahstram"]["PE"] = round(float(discordant_ratio), 1)
        record.samples["ahstram"]["SR"] = round(float(split_ratio), 1)
        record.samples["ahstram"]["SC"] = round(float(snv_coverage), 1)
        record.samples["ahstram"]["NS"] = int(snvs)
        record.samples["ahstram"]["HA"] = round(float(heterozygous_allele_ratio), 1)
        record.samples["ahstram"]["NH"] = int(het_snvs)
        record.samples["ahstram"]["SQ"] = round(max(float(REF_QUAL), float(ALT_QUAL)), 1)
        record.samples["ahstram"]["GL"] = (round(float(REF_GENOTYPE_LIKELIHOOD), 1), round(float(HET_GENOTYPE_LIKELIHOOD), 1), round(float(HOM_GENOTYPE_LIKELIHOOD), 1))

        vcf.write(record)

vcf.close()

pysam.bcftools.sort("-o", "example.vcf.gz", "example.vcf", catch_stdout = False)
