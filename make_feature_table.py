import numpy as np
from pybedtools import BedTool
import pysam
import pandas as pd
import sys
import argparse
import datetime

# Load the GC_content_reference.txt into a dictionary:
def make_GC_content_reference_table(GC_content_reference_filepath):
    GC_content_reference_table = {}
    with open(GC_content_reference_filepath, "r") as f:
        for line in f:
            linesplit = line.rstrip().split("\t")
            library_type, GC_content, gc_norm_factor = str(linesplit[0]), int(linesplit[1]), float(linesplit[2])
            GC_content_reference_table[(library_type, GC_content)] = gc_norm_factor
    return GC_content_reference_table

# Load the bed file (used for estimating coverage) into a dictionary:
def make_regions_table(bed_regions_filepath):
    regions_table = {}
    with open(bed_regions_filepath, "r") as f:
        for line in f:
            entry = line.rstrip().split("\t")
            chrom, start, end = str(entry[0]).removeprefix("chr"), int(entry[1]), int(entry[2])
            if chrom not in regions_table: regions_table[chrom] = []
            region = "{}:{}-{}".format(chrom, start, end)
            regions_table[chrom].append(region)
    return regions_table 

# Get the alignment preprocessing data (coverage, median read lengths, median insert sizes, etc.):
def make_alignment_preprocessing_table(alignment_filepath, reference_filepath, chroms, regions_table, threads = 1):
    alignment_preprocessing_table = {}

    genome_coverage_values: list[int] = []
    genome_median_read_lengths: list[int] = [] 
    genome_median_insert_sizes: list[int] = []
    genome_mad_insert_sizes: list[int] = []
    genome_size: int = 0

    # Determine whether the alignment file is a bam file or a cram file:
    mode = "rc"
    if alignment_filepath.endswith(".bam"): mode = "rb"
    alignment_iterator = pysam.AlignmentFile(alignment_filepath, mode = mode, reference_filename = reference_filepath, threads = threads)

    # Later on, I will want to allow the user to specify the index_filename (i.e. the explicit path to the index file)
    # https://readthedocs.org/projects/pysam/downloads/pdf/stable/

    for chrom in alignment_iterator.references: 
        if chrom.removeprefix("chr") not in chroms: continue
        if chrom.removeprefix("chr") not in regions_table: continue

        chrom_size: int = 0
        read_count: int = 0
        read_lengths: list[int] = []
        insert_sizes: list[int] = []

        # The keys to these dictionaries will be "{alignment.query_name}:{alignment.is_read1}"
        query_alignment_lengths: dict[str, int] = {}
        template_lengths: dict[str, int] = {}

        regions = regions_table[chrom.removeprefix("chr")]
        for region in regions:
           start = int(region.split(":")[-1].split("-")[0])
           end   = int(region.split(":")[-1].split("-")[1])
           chrom_size += end - start
           for alignment in alignment_iterator.fetch(contig = chrom, start = start, end = end):
               if (alignment.is_proper_pair == False): continue
               if (alignment.is_qcfail == True): continue
               if (alignment.is_duplicate == True): continue
               if (alignment.mapping_quality < 40): continue

               query_name = alignment.query_name
               is_read1 = alignment.is_read1
               query_name_read1 = "{}-{}".format(query_name, is_read1)
               query_alignment_lengths[query_name_read1] = alignment.query_alignment_length
               template_lengths[query_name_read1] = abs(alignment.template_length)                

        read_count = len(query_alignment_lengths)
        if read_count == 0: continue
        for key in query_alignment_lengths:
            read_lengths.append(query_alignment_lengths[key])
            insert_sizes.append(template_lengths[key])


        normalized_chrom_coverage = np.around( ((read_count/chrom_size) * np.median(read_lengths)), 3)
        median_read_length = np.median(read_lengths)
        median_insert_size = np.median(insert_sizes)
        c = 0.6745
        mad_insert_size = np.around(np.median(np.fabs(insert_sizes - np.median(insert_sizes))/c), 3)

        alignment_preprocessing_table[chrom] = {"normalized_chrom_coverage": normalized_chrom_coverage, "median_read_length": median_read_length, "median_insert_size": median_insert_size, "mad_insert_size": mad_insert_size, "chrom_size": chrom_size}

        genome_coverage_values.append(normalized_chrom_coverage)
        genome_median_read_lengths.append(median_read_length)
        genome_median_insert_sizes.append(median_insert_size)
        genome_mad_insert_sizes.append(mad_insert_size)
        genome_size += chrom_size
    
    
    alignment_preprocessing_table["GENOME"] = {"normalized_chrom_coverage": np.median(genome_coverage_values), "median_read_length": np.median(genome_median_read_lengths), "median_insert_size": np.median(genome_median_insert_sizes), "mad_insert_size": np.median(genome_mad_insert_sizes), "chrom_size": genome_size}
    return alignment_preprocessing_table


# Get the SNV preprocessing data (from the SNV VCF):
def get_snv_preprocessing_data(snv_vcf_filepath, chroms, regions_table, threads = 1):
    snv_preprocessing_table = {} 
    snv_vcf_iterator = pysam.VariantFile(snv_vcf_filepath, mode = "r", threads) 
    for chrom in snv_vcf_iterator.header.contigs:
        if chrom.removeprefix("chr") not in chroms: continue
        if chrom.removeprefix("chr") not in regions_table: continue

        chrom_depths: list[int] = []

        regions = regions_table[chrom.removeprefix("chr")]
        for region in regions:
            start = int(region.split(":")[-1].split("-")[0])
            end   = int(region.split(":")[-1].split("-")[1])
            
            try:
                records = snv_vcf_iterator.fetch(region = "{}:{}-{}".format(chrom, start, end))
            except ValueError:
                print("WARNING: region {}:{}-{} not present in SNV VCF.\n".format(chrom, start, end), file = sys.stderr)
                continue

            for record in records:
                if "AD" not in record.samples[0]: continue
                if "DP" not in record.samples[0]: continue
                if "GT" not in record.samples[0]: continue
                chrom_depths.append(record.samples[0]["DP"])

        snv_preprocessing_table[chrom] = {"median_chrom_depth": 0., "len_chrom_depth": 0}
        if len(chrom_depths) > 0:
            snv_preprocessing_table[chrom]["median_chrom_depth"] = np.nanmedian(chrom_depths) 
            snv_preprocessing_table[chrom]["len_chrom_depth"] = len(chrom_depths) 
    return snv_preprocessing_table


# Make an sv_interval_table that excludes the exclude_bed and contains the GC content per region
def make_sv_interval_table(sv_bed, exclude_bed, reference_fasta):
    def make_sv_interval_table_items(bed, item_name, sv_interval_table):
        for sv in bed:
            ID = sv[3].split(";")[-1].split("#")
            chrom, start, end, svtype = ID[0], int(ID[1]), int(ID[2]), ID[3]
            if (chrom, start, end, svtype) in sv_interval_table:
                sv_interval_table[(chrom, start, end, svtype)][item_name].append((sv[0], int(sv[1]), int(sv[2])))
        return sv_interval_table


    # Add ID to fourth column of BED file 
    sv_list = []
    for index, sv in enumerate(sv_bed): 
        ID = "#".join([sv[0], sv[1], sv[2], sv[3]])
        sv_list.append([sv[0], sv[1], sv[2], sv[3] + ";" + ID])
    sv_list_bed = BedTool(sv_list)
    sv_slop_bed = sv_list_bed.slop(l = 500, r = 500, genome = "hg38")
    sv_slop_df = sv_slop_bed.to_dataframe()
    sv_slop_df["start_plus_1000"] = sv_slop_df["start"] + 1000 + 1
    sv_slop_df["end_minus_1000"] = sv_slop_df["end"] - 1000 
    sv_slop_df["end"] = sv_slop_df["end"] + 1 
    sv_slop_left_flank_df = sv_slop_df[["chrom", "start", "start_plus_1000", "name"]]
    sv_slop_right_flank_df = sv_slop_df[["chrom", "end_minus_1000", "end", "name"]]
    sv_slop_left_flank_bed = BedTool.from_dataframe(sv_slop_left_flank_df)
    sv_slop_right_flank_bed = BedTool.from_dataframe(sv_slop_right_flank_df)

    sv_subtract_exclude_and_nucleotide_content_bed = sv_list_bed.subtract(exclude_bed).nucleotide_content(fi = reference_fasta)
    sv_slop_left_flank_exclude_bed = sv_slop_left_flank_bed.subtract(exclude_bed)
    sv_slop_right_flank_exclude_bed = sv_slop_right_flank_bed.subtract(exclude_bed)

    svtypes = ("DEL", "DUP")
    sv_interval_table = {}
    for index, sv in enumerate(sv_subtract_exclude_and_nucleotide_content_bed):
        ID = sv[3].split(";")[-1].split("#")
        chrom, start, end, svtype = ID[0], int(ID[1]), int(ID[2]), ID[3]
        if svtype not in svtypes: continue
        if (chrom, start, end, svtype) not in sv_interval_table: sv_interval_table[(chrom, start, end, svtype)] = {"SV_SPAN": [], "FLANK_SPAN": [], "WINDOWS": [], "nucleotide_content": []}
        sv_interval_table[(chrom, start, end, svtype)]["SV_SPAN"].append((sv[0], int(sv[1]), int(sv[2])))
        sv_interval_table[(chrom, start, end, svtype)]["nucleotide_content"].append((float(sv[4]), float(sv[5]), int(sv[6]), int(sv[7]), int(sv[8]), int(sv[9]), int(sv[10]), int(sv[11]), int(sv[12])))

    # Get flank spans too
    sv_interval_table = make_sv_interval_table_items(sv_slop_left_flank_exclude_bed,  "FLANK_SPAN", sv_interval_table)
    sv_interval_table = make_sv_interval_table_items(sv_slop_right_flank_exclude_bed, "FLANK_SPAN", sv_interval_table)
    # Get windows too
    sv_interval_table = make_sv_interval_table_items(sv_slop_left_flank_bed,  "WINDOWS", sv_interval_table)
    sv_interval_table = make_sv_interval_table_items(sv_slop_right_flank_bed, "WINDOWS", sv_interval_table)
    return sv_interval_table

def make_snv_features_table(snv_vcf_filepath, sv_bed, sv_interval_table, svtypes, df_preprocessing_table, threads = 1):
    snv_features_table = {}
    snv_vcf_iterator = pysam.VariantFile(snv_vcf_filepath, threads = threads)
    for sv in sv_bed:
        chrom, start, end, svtype = sv[0], int(sv[1]), int(sv[2]), sv[3]
        if (chrom, start, end, svtype) not in sv_interval_table: continue
        if svtype not in svtypes: continue
        svlen = end - start

        sv_span = sv_interval_table[(chrom, start, end, svtype)]["SV_SPAN"]
        flank_span = sv_interval_table[(chrom, start, end, svtype)]["FLANK_SPAN"]
        windows = sv_interval_table[(chrom, start, end, svtype)]["WINDOWS"]
        
        locus_depths: list[int] = []
        locus_HADs: list[int] = [] # heterozygous allele depth
        for locus in sv_span:
            chrom_, start_, end_ = locus[0], locus[1], locus[2]

            try: 
                records = snv_vcf_iterator.fetch(region = "{}:{}-{}".format(chrom, str(start_), str(end_)))
            except ValueError:
                print("WARNING: region {}:{}-{} not present in SNV VCF.\n".format(chrom, str(start_), str(end_)), file = sys.stderr)
                continue

            for record in records:
                if "AD" not in record.samples[0]: continue
                if "DP" not in record.samples[0]: continue
                if "GT" not in record.samples[0]: continue

                locus_depths.append(record.samples[0]["DP"])
                if record.samples[0]["AD"][0] == 0 or record.samples[0]["AD"][1] == 0: continue
                AR = float(record.samples[0]["AD"][0])/float(record.samples[0]["AD"][1])
                if AR > 1: AR = 1/AR
                locus_HADs.append(AR)
       
        snv_features_table[(chrom, start, end)] = {}
        snv_features_table[(chrom, start, end)]["snv_coverage"] = np.nan
        if chrom in df_preprocessing_table["chrom"].values: snv_features_table[(chrom, start, end)]["snv_coverage"] = float(np.nanmedian(locus_depths))/df_preprocessing_table[df_preprocessing_table["chrom"] == chrom]["median_chrom_depth"].values[0]
        snv_features_table[(chrom, start, end)]["heterozygous_allele_ratio"] = np.nanmedian(locus_HADs)
        snv_features_table[(chrom, start, end)]["snvs"] = len(locus_depths)
        snv_features_table[(chrom, start, end)]["het_snvs"] = len(locus_HADs) 

        if np.isnan(snv_features_table[(chrom, start, end)]["snv_coverage"]): snv_features_table[(chrom, start, end)]["snv_coverage"] = "NaN"
        if np.isnan(snv_features_table[(chrom, start, end)]["heterozygous_allele_ratio"]): snv_features_table[(chrom, start, end)]["heterozygous_allele_ratio"] = "NaN"

    return snv_features_table

def make_alignment_features_table(alignment_filepath, reference_filepath, sv_bed, df_preprocessing_table, sv_interval_table, svtypes, GC_content_reference_table, threads = threads):
    def calculate_gc_content_fraction(sv_interval_table_nucleotide_content_list):
        # https://daler.github.io/pybedtools/autodocs/pybedtools.bedtool.BedTool.nucleotide_content.html
        # The elements in each sublist are as follows: 
        # 1) %AT content
        # 2) %GC content
        # 3) Number of As observed
        # 4) Number of Cs observed
        # 5) Number of Gs observed
        # 6) Number of Ts observed
        # 7) Number of Ns observed
        # 8) Number of other bases observed
        # 9) The length of the explored sequence/interval
        # sv_interval_table[(chrom, start, end, svtype)]["nucleotide_content"].append((sv[4], sv[5], sv[6], sv[7], sv[8], sv[9], sv[10], sv[11], sv[12]))
        num_G: int = 0
        num_C: int = 0
        sequence_length: int = 0
        for sublist in sv_interval_table_nucleotide_content_list:
            num_G += sublist[4]
            num_C += sublist[3]
            sequence_length += sublist[8]
        return float(num_G + num_C)/sequence_length 

    alignment_features_table = {}

    ci_insert_size_insert_mad = df_preprocessing_table[df_preprocessing_table["chrom"] == "GENOME"]["median_insert_size"].values[0] + 5.0*df_preprocessing_table[df_preprocessing_table["chrom"] == "GENOME"]["mad_insert_size"].values[0]
    median_read_length = df_preprocessing_table[df_preprocessing_table["chrom"] == "GENOME"]["median_read_length"].values[0]
    c = 0.6745

    def check_read(read):
        if read.mapping_quality >= 40 and read.query_length > median_read_length - 10: return True
        return False

    mode = "rc"
    if alignment_filepath.endswith(".bam"): mode = "rb"
    alignment_iterator = pysam.AlignmentFile(alignment_filepath, mode = mode, reference_filename = reference_filepath, threads = threads)
    for index, sv in enumerate(sv_bed):
        chrom, start, end, svtype = sv[0], int(sv[1]), int(sv[2]), sv[3]
        if (chrom, start, end, svtype) not in sv_interval_table: continue 
        if svtype not in svtypes: continue
        svlen = end - start

        sv_span = sv_interval_table[(chrom, start, end, svtype)]["SV_SPAN"]
        flank_span = sv_interval_table[(chrom, start, end, svtype)]["FLANK_SPAN"]
        windows = sv_interval_table[(chrom, start, end, svtype)]["WINDOWS"]

        # Initialize coverage, which is estimated differently depending on the length of the SV. 
        coverage = np.nan
        coverage_GCcorrected = np.nan
        if chrom not in df_preprocessing_table["chrom"].values: pass
        elif svlen > 1000:
            read_count: int = 0
            alignment_count: int = 0
            basepair_span: int = 0
            for span in sv_span:
                for alignment in alignment_iterator.fetch(contig = span[0], start = span[1], end = span[2]):
                    basepair_span += span[2] - span[1]
                    alignment_count += 1
            coverage = ((float(alignment_count)/svlen)*df_preprocessing_table[df_preprocessing_table["chrom"] == "GENOME"]["median_read_length"].values[0])/df_preprocessing_table[df_preprocessing_table["chrom"] == chrom]["normalized_chrom_coverage"].values[0]
        else: # When SV length is less than 1000, we can use this more cumbersome method to estimate coverage.
            median_depth_of_coverage = np.nan
            positional_depth_of_coverage = {}
            for span in sv_span:
                depth_result = alignment_iterator.count_coverage(contig = span[0], start = span[1] - 1, stop = span[2], read_callback = check_read)
                for i in range(len(depth_result[0])): 
                    positional_depth_of_coverage[span[1] + 1 + i] = depth_result[0][i] + depth_result[1][i] + depth_result[2][i] + depth_result[3][i]
            positional_depth_of_coverage_array = np.array(list(positional_depth_of_coverage.values()))
            median_depth_of_coverage = np.median(positional_depth_of_coverage_array)
            coverage = float(median_depth_of_coverage)/df_preprocessing_table[df_preprocessing_table["chrom"] == chrom]["normalized_chrom_coverage"].values[0]

        # Correcting coverage for GC content...
        gc_content_fraction = calculate_gc_content_fraction(sv_interval_table[(chrom, start, end, svtype)]["nucleotide_content"])
        base = 5
        gc_content = int(base * round((100*gc_content_fraction)/5))
        gc_norm_factor_read_count = 1.0
        gc_norm_factor_depth_of_coverage = 1.0
        if ("PCR-READCOUNT", gc_content) in GC_content_reference_table: gc_norm_factor_read_count = GC_content_reference_table[("PCR-DOC", gc_content)]
        if ("PCR-DOC", gc_content) in GC_content_reference_table: gc_norm_factor_depth_of_coverage = GC_content_reference_table[("PCR-DOC", gc_content)]
        # DEBUGGING
        # Maybe change adjusted coverage value?
        coverage_GCcorrected = gc_norm_factor_read_count * coverage 

        split_read_count: int = 0
        concordant_read_count: int = 0
        discordant_read_count: int = 0

        split_read_ratio: float = 0.0
        discordant_read_ratio: float = 0.0
        for flank in flank_span:
            for alignment in alignment_iterator.fetch(contig = flank[0], start = flank[1], end = flank[2]):
                if alignment.is_qcfail == True: continue
                if alignment.is_duplicate == True: continue
                if alignment.mapping_quality < 10: continue
                if alignment.is_unmapped: continue
                if alignment.is_reverse == alignment.mate_is_reverse: continue
                mate_position = alignment.next_reference_start
                
                if abs(alignment.template_length) >= ci_insert_size_insert_mad:
                    """ insert size is greater than 5 MAD (median absolute deviations)  from median """
                    # Get discordant reads
                    if (windows[0][1] <= alignment.reference_start <= windows[0][2] and windows[1][1] <= mate_position <= windows[1][2]) or (windows[1][1] <= alignment.reference_start <= windows[1][2] and windows[0][1] <= mate_position <= windows[0][2]):
                        discordant_read_count += 1
                    # DEBUGGING
                    # TEST ON original 1000 Genomes file and make sure that split-read ratio values are consistent
                    # Get split reads (but they're all 0 if we do it this way)
                    if alignment.is_supplementary:
                        second_alignment = alignment.get_tag("SA").split(",")
                        if len(second_alignment) != 0:
                            if second_alignment[0] == chrom: # The value c that Danny uses for some reason (and leads to a split_read count of 0...) # c SHOULD BE chromosome instead (the second alignment should be on the same chromosome as the first one)
                                if (windows[0][1] <= alignment.reference_start <= windows[0][2] and windows[1][1] <= int(second_alignment[1]) - 1 <= windows[1][2]) or (windows[1][1] <= alignment.reference_start <= windows[1][2] and windows[0][1] <= int(second_alignment[1]) - 1 <= windows[0][2]):
                                    split_read_count += 1
                
                if (not alignment.is_supplementary) and alignment.is_proper_pair and abs(alignment.template_length) < ci_insert_size_insert_mad:
                    concordant_read_count += 1
                
                if concordant_read_count == 0:
                    discordant_read_ratio = discordant_read_count/1.0
                    split_read_ratio = split_read_count/1.0
                else:
                    discordant_read_ratio = discordant_read_count/float(concordant_read_count)
                    split_read_ratio = split_read_count/float(concordant_read_count)

                # if alignment.is_supplementary: split_read_count += 1
                # if alignment.is_proper_pair: concordant_read_count += 1 
                # else: discordant_read_count += 1


        alignment_features_table[(chrom, start, end)] = {}
        alignment_features_table[(chrom, start, end)]["type"] = svtype
        alignment_features_table[(chrom, start, end)]["size"] = svlen
        alignment_features_table[(chrom, start, end)]["coverage"] = coverage
        alignment_features_table[(chrom, start, end)]["coverage_GCcorrected"] = coverage_GCcorrected
        alignment_features_table[(chrom, start, end)]["discordant_ratio"] = discordant_read_ratio 
        alignment_features_table[(chrom, start, end)]["split_ratio"] = split_read_ratio 

    return alignment_features_table



if __name__ == "__main__":
    # Argument parsing (when this script is run by itself)
    parser = argparse.ArgumentParser(description = "SV2 genotyper")
    parser.add_argument("--alignment_file", help = "CRAM/BAM file input")
    parser.add_argument("--reference_fasta", help = "Reference fasta file input")
    parser.add_argument("--snv_vcf_file", help = "SNV VCF file input")
    parser.add_argument("--regions_bed", help = "BED file with pre-generated random genomic regions (for estimating coverage per chromosome)")
    parser.add_argument("--exclude_regions_bed", help = "BED file with pre-generated random genomic regions (for estimating coverage per chromosome)")
    parser.add_argument("--sv_bed_file", help = "SV BED file input")
    parser.add_argument("--sv_feature_output_tsv", help = "Feature table .tsv output (it will be given the default name if this arugment isn't set)")
    parser.add_argument("--preprocessing_table_input", help = "A pre-generated preprocessing table (if SV2 had been run before and you want to skip the preprocessing part of the program")
    parser.add_argument("--preprocessing_table_output", help = "Preprocessing table .tsv output (it will be given the default name if this argument isn't set)")
    parser.add_argument("--gc_reference_table", help = "GC content reference table input")
    args = parser.parse_args()

    # Chromosomes to analyze (1,... 22, X, Y)
    chroms = list(map(str, np.arange(1, 22 + 1)))
    chroms.extend(["X", "Y"])

    # SV types to classify
    svtypes = ("DEL", "DUP")

    reference_filepath = args.reference_fasta 

    # Make the GC content reference table
    GC_content_reference_table_filepath = args.gc_reference_table
    GC_content_reference_table = make_GC_content_reference_table(GC_content_reference_table_filepath)

    # Make the regions table 
    random_regions_table = args.regions_bed
    regions_table = make_regions_table(random_regions_table)

    alignment_filepath = args.alignment_file
    snv_vcf_filepath = args.snv_vcf_file
    if not args.preprocessing_table_input:
        # Make the CRAM/BAM preprocessing table
        alignment_preprocessing_table = make_alignment_preprocessing_table(alignment_filepath, reference_filepath, chroms, regions_table)
        
        # Make the SNV preprocessing table
        snv_preprocessing_table = get_snv_preprocessing_data(snv_vcf_filepath, chroms, regions_table)
    
        # Merge and output the preprocessing table
        df_alignment_preprocessing_table = pd.DataFrame.from_dict(alignment_preprocessing_table, orient = "index")
        df_snv_preprocessing_table = pd.DataFrame.from_dict(snv_preprocessing_table, orient = "index")
        df_preprocessing_table = df_alignment_preprocessing_table.join(df_snv_preprocessing_table).reset_index(level = 0).rename(columns = {"index": "chrom"})
        df_preprocessing_table.to_csv(args.preprocessing_table_output, sep = "\t", index = False)

    else:
        df_preprocessing_table = pd.read_csv(args.preprocessing_table_input, sep = "\t")

    # Preprocess the input SV BED file
    sv_bed_filepath = args.sv_bed_file
    sv_bed_list = []
    with open(sv_bed_filepath, "r") as f:
        for line in f:
            linesplit = line.rstrip().split("\t")
            chrom, start, stop, features = linesplit[0], int(linesplit[1]), int(linesplit[2]), linesplit[3]
            if start > stop: continue
            sv_bed_list.append(line.rstrip())
    sv_bed = BedTool(sv_bed_list).filter(lambda x: len(x) > 0).saveas()
    exclude_bed = BedTool(args.exclude_regions_bed).merge()
    sv_interval_table = make_sv_interval_table(sv_bed, exclude_bed, args.reference_fasta)

    # Make SNV features table (for each filtered SV call)
    snv_features_table = make_snv_features_table(snv_vcf_filepath, sv_bed, sv_interval_table, svtypes, df_preprocessing_table)
    
    # Make CRAM/BAM features table (for each filtered SV call)
    alignment_features_table = make_alignment_features_table(alignment_filepath, reference_filepath, sv_bed, df_preprocessing_table, sv_interval_table, svtypes, GC_content_reference_table)
    
    # Merge and output feature table
    df_alignment_features_table = pd.DataFrame.from_dict(alignment_features_table, orient = "index")
    df_snv_features_table = pd.DataFrame.from_dict(snv_features_table, orient = "index")
    df_features_table = df_alignment_features_table.join(df_snv_features_table).reset_index().rename(columns = {"level_0": "chrom", "level_1": "start", "level_2": "end"})
    
    if not args.sv_feature_output_tsv: 
        from time import gmtime, strftime
        current_time = strftime("%Y-%m-%d_%H.%M.%S", gmtime())
        df_features_table.to_csv("sv2_features.{}.tsv".format(current_time), sep = "\t", index = False)
    else: df_features_table.to_csv(args.sv_feature_output_tsv, sep = "\t", index = False)
    
    # Sample commands:
    # python sv2_refactor_args.py --reference_fasta /expanse/lustre/projects/ddp195/j3guevar/resources/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa --alignment_file /expanse/projects/sebat1/genomicsdataanalysis/resources/NA12878/illumina_platinum_pedigree/NA12878.alt_bwamem_GRCh38DH.20150706.CEU.illumina_platinum_ped.cram --regions_bed my_random.bed --gc_reference_table GC_content_reference.txt --snv_vcf_file /home/j3guevar/tests/SV2_nim/src/SV2_nimpkg/refactored_sv2/input_files/NA12878.vcf.gz --sv_bed_file /home/j3guevar/tests/sv2_refactored/test_inputs/file.nbl.bed
    # python sv2_refactor_args.py --reference_fasta /expanse/lustre/projects/ddp195/j3guevar/resources/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa --alignment_file /expanse/projects/sebat1/genomicsdataanalysis/resources/NA12878/illumina_platinum_pedigree/NA12878.alt_bwamem_GRCh38DH.20150706.CEU.illumina_platinum_ped.cram --regions_bed my_random.bed --gc_reference_table GC_content_reference.txt --snv_vcf_file /home/j3guevar/tests/SV2_nim/src/SV2_nimpkg/refactored_sv2/input_files/NA12878.vcf.gz --sv_bed_file /home/j3guevar/tests/sv2_refactored/test_inputs/file.nbl.10.bed --preprocessing_table_input preprocessing_table_output_test.txt --exclude_regions_bed /home/j3guevar/tests/SV2_nim/src/SV2_nimpkg/resources/hg38_excluded.bed.gz
