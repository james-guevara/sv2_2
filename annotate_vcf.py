from pysam import VariantFile
from pybedtools import BedTool

a = BedTool("sv2_output/NA12878_sv2_2022-07-14_23.10.37.vcf.gz")
b = BedTool("data/annotation_files/hg38_cytoband.bed")
# bedtools intersect usage and option summary: https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html
# wao: Write the original A and B entries plus the number of base pairs of overlap between the two features. However, A features w/o overlap are also reported with a NULL B feature and overlap = 0. Restricted by -f and -r.
# f: Minimum overlap required as a fraction of A. Default is 1E-9 (i.e. 1bp).
# F: Minimum overlap required as a fraction of B. Default is 1E-9 (i.e., 1bp).
# wa: Write the original entry in A for each overlap.
# wb: Write the original entry in B for each overlap. Useful for knowing what A overlaps. Restricted by -f and -r.

cytoband_intersection = a.intersect(b, wao = True) # Can output multiple lines for the same variant in the SV2 VCF
# repeatmasker_intersection = a.intersect("data/annotation_files/hg38_repeatmasker.bed.gz", f = 0.8, F = 0.8, wa = True, wb = True)
thousand_genomes_deletions_intersection = a.intersect("data/annotation_files/hg38_1000Genomes_DEL.bed", f = 0.8, F = 0.8, wao = True)
thousand_genomes_duplications_intersection = a.intersect("data/annotation_files/hg38_1000Genomes_DUP.bed", f = 0.8, F = 0.8, wao = True)
# genes_intersection = a.intersect("data/annotation_files/hg38_genes.bed.gz", wa = True, wb = True)

def bedtool_to_dictionary(bed_intersection, annotation_type):
    intersection_table = {}
    for sv in bed_intersection: 
        chrom, pos, svtype, info = sv[0], sv[1], sv[4], sv[7].split(";")
        end, svtype, svlen = info[0].split("=")[1], info[1].split("=")[1], info[2].split("=")[1] # Overwriting the previous svtype but it should be the same thing anyway (and without the "<>" characters)
        key = "{}:{}:{}:{}:{}".format(chrom, pos, end, svtype, svlen)

        if annotation_type == "cytoband": 
            cytoband = sv[-3]
            intersection_table[key] = "{}{}".format(chrom.removeprefix("chr"), cytoband)
        elif annotation_type == "thousand_genomes_deletions":
            if svtype != "DEL": continue
            intersection_table[key] = {"1000G_ID": "NA", "1000G_OVERLAP": 0}
            if sv[-1] == "0": continue
            intersection_table[key]["1000G_ID"] = sv[-2]
            intersection_table[key]["1000G_OVERLAP"] = float(sv[-1])/float(svlen)
        elif annotation_type == "thousand_genomes_duplications":
            if svtype != "DUP": continue
            intersection_table[key] = {"1000G_ID": "NA", "1000G_OVERLAP": 0}
            if sv[-1] == "0": continue
            intersection_table[key]["1000G_ID"] = sv[-2]
            intersection_table[key]["1000G_OVERLAP"] = float(sv[-1])/float(svlen)
        elif annotation_type == "genes":
            intersection_table[key] = sv[-1]

    return intersection_table

cytoband_table = bedtool_to_dictionary(cytoband_intersection, "cytoband")
thousand_genomes_deletions_table = bedtool_to_dictionary(thousand_genomes_deletions_intersection, "thousand_genomes_deletions")
thousand_genomes_duplications_table = bedtool_to_dictionary(thousand_genomes_duplications_intersection, "thousand_genomes_duplications")
genes_table = bedtool_to_dictionary(genes_intersection, "genes")

vcf_in = VariantFile("sv2_output/NA12878_sv2_2022-07-14_23.10.37.vcf.gz")


vcf_out = VariantFile("sv2_output/NA12878_sv2_2022-07-14_23.10.37.annotated.vcf.gz", mode = "w", header = vcf_in.header)
vcf_out.header.add_meta(key = "INFO", items = [ ("ID", "CYTOBAND"), ("Number", 1), ("Type", "String"), ("Description", "Cytoband(s) overlapping the variant") ] )
# Any INFO field that contains "1000G" is reserved (I think) based on looking at this line (line 468 as of 2022-07-15): https://github.com/samtools/htslib/blob/develop/vcf.c
# """
# if ( !strcmp(val,"1000G") ) return 0;
# """
# So I'm renaming them...
vcf_out.header.add_meta(key = "INFO", items = [ ("ID", "THOUSAND_GENOMES_ID_SV2"), ("Number", 1), ("Type", "String"), ("Description", "1000 Genomes Phase 3 integrated SV callset variant identifier") ] )
vcf_out.header.add_meta(key = "INFO", items = [ ("ID", "THOUSAND_GENOMES_OVERLAP_SV2"), ("Number", 1), ("Type", "Float"), ("Description", "Overlap to 1000 Genomes Phase 3 variant, in the range (0,1)") ] )
vcf_out.header.add_meta(key = "INFO", items = [ ("ID", "GENES"), ("Number", 1), ("Type", "String"), ("Description", "Genes overlapping the variant, pipe-separated by transcripts") ] )
# header.add_meta(key = "INFO", items = [ ("ID", "REPEATMASKER"), ("Number", 2), ("Type", "String"), ("Description", "Name and reciprocal overlap of RepeatMasker variant") ] )

format_tuple = ("GT", "PE", "SR", "SC", "NS", "HA", "NH", "SQ", "GL")
for record in vcf_in:
    new_record = vcf_out.header.new_record(contig = record.chrom, start = record.pos, stop = record.stop, alleles = record.alleles, id = record.id, qual = record.qual, filter = record.filter, info = record.info)
    key = "{}:{}:{}:{}:{}".format(record.chrom, record.pos, record.stop, record.info["SVTYPE"], record.info["SVLEN"]) # record.info["END"] doesn't work, so I use this instead which outputs the appropriate info field ("END")

    # Add new annotations
    new_record.info["CYTOBAND"] = cytoband_table[key] 
    if record.info["SVTYPE"] == "DEL":
        new_record.info["THOUSAND_GENOMES_ID_SV2"] = thousand_genomes_deletions_table[key]["1000G_ID"]
        new_record.info["THOUSAND_GENOMES_OVERLAP_SV2"] = thousand_genomes_deletions_table[key]["1000G_OVERLAP"]
    elif record.info["SVTYPE"] == "DUP":
        new_record.info["THOUSAND_GENOMES_ID_SV2"] = thousand_genomes_duplications_table[key]["1000G_ID"]
        new_record.info["THOUSAND_GENOMES_OVERLAP_SV2"] = thousand_genomes_duplications_table[key]["1000G_OVERLAP"]
    # I'm not going to annotate the genes (leave that to another annotation utility)
    # if key not in genes_table: new_record.info["GENES"] = "intergenic"
    # else: new_record.info["GENES"] = genes_table[key]


    # Add format fields for sample (there should only be 1 sample)
    for element in format_tuple: new_record.samples[0][element] = record.samples[0][element]
    vcf_out.write(new_record) 
