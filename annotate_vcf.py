import pysam
from pybedtools import BedTool
import sys

sv_bed_list = []
vcf_iterator = pysam.VariantFile(sys.argv[1])
for record in vcf_iterator:
    chrom, start, stop, svtype = record.chrom, record.start, record.stop, record.info["SVTYPE"]
    sv_bed_list.append("\t".join([chrom, str(start), str(stop), svtype]))
vcf_iterator.close()
sv_bed = BedTool(sv_bed_list)

repeatmasker_bed = BedTool("data/annotation_files/hg38_repeatmasker.bed.gz")
sv2_repeatmasker_intersection = sv_bed.intersect(repeatmasker_bed, f = 0.8, r = True, wa = True, wb = True)
repeatmasker_dict = {}
for interval in sv2_repeatmasker_intersection:
    chrom, start, stop, svtype = interval[0], interval[1], interval[2], interval[3]
    ID = "{}:{}:{}:{}".format(chrom, start, stop, svtype)
    if ID not in repeatmasker_dict: repeatmasker_dict[ID] = {"REPEATMASKER": None, "REPEATMASKER_OVERLAP": 0.}

    name = "{}:{}:{}".format(interval[8], interval[9], interval[7])
    sv_start, sv_end = int(interval[1]), int(interval[2])
    repeatmasker_start, repeatmasker_end = int(interval[5]), int(interval[6])
    overlap = min(sv_end, repeatmasker_end) - max(sv_start, repeatmasker_start)

    reciprocal_overlap = min( float(overlap)/(sv_end - sv_start), float(overlap)/(repeatmasker_end - repeatmasker_start) ) 
    repeatmasker_dict[ID]["REPEATMASKER"] = name
    repeatmasker_dict[ID]["REPEATMASKER_OVERLAP"] = reciprocal_overlap

excluded_regions_bed = BedTool("data/excluded_regions_bed_files/hg38_excluded.bed.gz")
sv2_excluded_regions_intersection = sv_bed.intersect(excluded_regions_bed, wao = True)
excluded_regions_dict = {}
for interval in sv2_excluded_regions_intersection:
    chrom, start, stop, svtype = interval[0], interval[1], interval[2], interval[3]
    ID = "{}:{}:{}:{}".format(chrom, start, stop, svtype)
    if ID not in excluded_regions_dict: excluded_regions_dict[ID] = {"abparts": 0, "centromere": 0, "gap": 0, "segDup": 0, "STR": 0}

    excluded_region_type = interval[-2]
    if excluded_region_type == ".": continue
    sv_basepair_overlap = int(interval[-1])
    excluded_regions_dict[ID][excluded_region_type] += sv_basepair_overlap


vcf_iterator = pysam.VariantFile(sys.argv[1])
vcf_out = pysam.VariantFile(sys.argv[2], mode = "w", header = vcf_iterator.header)

vcf_out.header.add_meta(key = "FILTER", items = [ ("ID", "ABPARTS"), ("Description", "Variant (50% or more) overlaps antibody parts region.") ] )
vcf_out.header.add_meta(key = "FILTER", items = [ ("ID", "GAP"), ("Description", "Variant (50% or more) overlaps gap region.") ] )
vcf_out.header.add_meta(key = "FILTER", items = [ ("ID", "CENTROMERE"), ("Description", "Variant (50% or more) overlaps centromere region.") ] )
vcf_out.header.add_meta(key = "FILTER", items = [ ("ID", "SEGDUP"), ("Description", "Variant (50% or more) overlaps segmental duplication region.") ] )
vcf_out.header.add_meta(key = "FILTER", items = [ ("ID", "STR"), ("Description", "Variant (50% or more) overlaps STR region.") ] )

vcf_out.header.add_meta(key = "INFO", items = [ ("ID", "REPEATMASKER"), ("Number", 2), ("Type", "String"), ("Description", "Name and reciprocal overlap of RepeatMasker variant") ] )
vcf_out.header.add_meta(key = "INFO", items = [ ("ID", "ABPARTS"), ("Number", 1), ("Type", "Float"), ("Description", "Overlap to antibody parts, in the range (0,1)") ] )
vcf_out.header.add_meta(key = "INFO", items = [ ("ID", "CENTROMERE"), ("Number", 1), ("Type", "Float"), ("Description", "Centromere overlap, in the range (0,1)") ] )
vcf_out.header.add_meta(key = "INFO", items = [ ("ID", "GAP"), ("Number", 1), ("Type", "Float"), ("Description", "Overlap to gaps in the reference, in the range (0,1)") ] )
vcf_out.header.add_meta(key = "INFO", items = [ ("ID", "SEGDUP"), ("Number", 1), ("Type", "Float"), ("Description", "Segmental duplication overlap, in the range (0,1)") ] )
vcf_out.header.add_meta(key = "INFO", items = [ ("ID", "STR"), ("Number", 1), ("Type", "Float"), ("Description", "Short tandem repeat overlap, in the range (0,1)") ] )

for record in vcf_iterator:
    new_record = vcf_out.header.new_record(contig = record.chrom, start = record.start, stop = record.stop, alleles = record.alleles, id = record.id, qual = record.qual, filter = record.filter, info = record.info)
    ID = "{}:{}:{}:{}".format(record.chrom, record.start, record.stop, record.info["SVTYPE"])

    if ID in repeatmasker_dict: new_record.info["REPEATMASKER"] = (repeatmasker_dict[ID]["REPEATMASKER"], str(repeatmasker_dict[ID]["REPEATMASKER_OVERLAP"]))

    # Calculate fraction of overlap with each excluded region
    new_record.info["ABPARTS"] = excluded_regions_dict[ID]["abparts"]/float(record.info["SVLEN"])
    new_record.info["CENTROMERE"] = excluded_regions_dict[ID]["centromere"]/float(record.info["SVLEN"])
    new_record.info["GAP"] = excluded_regions_dict[ID]["gap"]/float(record.info["SVLEN"])
    new_record.info["SEGDUP"] = excluded_regions_dict[ID]["segDup"]/float(record.info["SVLEN"])
    new_record.info["STR"] = excluded_regions_dict[ID]["STR"]/float(record.info["SVLEN"])

    # Add filter flags if the variant overlap is greater than 50% (for each excluded region)
    if new_record.info["ABPARTS"] >= 0.5: new_record.filter.add("ABPARTS") 
    if new_record.info["CENTROMERE"] >= 0.5: new_record.filter.add("CENTROMERE") 
    if new_record.info["GAP"] >= 0.5: new_record.filter.add("GAP") 
    if new_record.info["SEGDUP"] >= 0.5: new_record.filter.add("SEGDUP") 
    if new_record.info["STR"] >= 0.5: new_record.filter.add("STR") 

    for i in range(len(record.samples)):
        for format_ in record.samples[i]: new_record.samples[i][format_] = record.samples[i][format_]

    vcf_out.write(new_record)
