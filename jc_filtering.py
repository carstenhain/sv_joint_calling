import vcf
import numpy as np
import pysam
import os
import mappy as mp

def filter_intersect (original_file_path, filtered_file_path, min_sv_len, filtered_intersect_file_path, bam_file_path):
    vcf_reader = vcf.Reader(open(original_file_path, "r"))
    vcf_writer = vcf.Writer(open(filtered_file_path, "w"), vcf_reader)

    ins_vcf_reader = vcf.Reader(filename=filtered_intersect_file_path)

    samfile = pysam.AlignmentFile(bam_file_path, "rb")

    ins_distance = 150

    svtypes = {}

    for record in vcf_reader:

        if record.INFO["SVTYPE"] == "TRA":
            svlen = 1000000
        else:
            svlen = np.abs(record.INFO["SVLEN"])

        if record.CHROM == "chrM":
            svlen = 1
        if record.INFO["CHR2"] == "chrM":
            svlen = 1

        if svlen >= min_sv_len:

            if record.INFO["SVTYPE"] == "INS":
                record.INFO["SUPP"] = 4
            else:

                # gather SV supporting reads

                reads_at_start = {}
                for read in samfile.fetch(record.CHROM, record.POS - 150, record.POS + 150):
                    if read.has_tag("SA"):
                        if not (read.query_name in reads_at_start):
                            reads_at_start[read.query_name] = read

                sv_reads = []
                for read in samfile.fetch(record.INFO["CHR2"], record.INFO["END"] - 150, record.INFO["END"] + 150):
                    if read.query_name in reads_at_start:
                        sv_reads.append(read)

                if record.INFO["SVTYPE"] in svtypes:
                    svtypes[record.INFO["SVTYPE"]] += 1
                else:
                    svtypes[record.INFO["SVTYPE"]] = 1

                #print(record.ID)

                insertion_in_region = []
                insertion_at_start = False
                insertion_at_end = False

                for unfiltered_variants in ins_vcf_reader.fetch(record.CHROM, record.POS - ins_distance,
                                                                record.POS + ins_distance):
                    if unfiltered_variants.INFO["SVTYPE"] == "INS" and np.abs(unfiltered_variants.INFO["SVLEN"]) > 200:
                        insertion_in_region.append(unfiltered_variants)
                        insertion_at_start = True
                        #print("\t" + unfiltered_variants.ID + " at start")

                for unfiltered_variants in ins_vcf_reader.fetch(record.INFO["CHR2"], record.INFO["END"] - ins_distance,
                                                                record.INFO["END"] + ins_distance):
                    if unfiltered_variants.INFO["SVTYPE"] == "INS" and np.abs(unfiltered_variants.INFO["SVLEN"]) > 200:
                        insertion_in_region.append(unfiltered_variants)
                        insertion_at_end = True
                        #print("\t" + unfiltered_variants.ID + " at end")

                if len(insertion_in_region) == 1:

                    # get alt and ref reference sequences as fasta

                    tmp_bed_fo = open("tmp.bed", "w")

                    start = np.amax([0, insertion_in_region[0].POS - 10000])
                    tmp_bed_fo.write("\t".join(
                        [
                            insertion_in_region[0].CHROM,
                            str(start),
                            str(insertion_in_region[0].POS)
                        ]) + "\n")

                    tmp_bed_fo.write("\t".join(
                        [
                            insertion_in_region[0].CHROM,
                            str(insertion_in_region[0].POS),
                            str(insertion_in_region[0].POS + 10000)
                        ]) + "\n")

                    start = np.amax([0, insertion_in_region[0].POS - 10000])
                    tmp_bed_fo.write("\t".join(
                        [
                            insertion_in_region[0].CHROM,
                            str(start),
                            str(insertion_in_region[0].POS + 10000)
                        ]) + "\n")

                    if not insertion_at_start:
                        start = np.amax([0, record.POS - 10000])
                        tmp_bed_fo.write("\t".join(
                            [
                                record.CHROM,
                                str(start),
                                str(record.POS + 10000)
                            ]) + "\n")
                    if not insertion_at_end:
                        start = np.amax([0, record.INFO["END"] - 10000])
                        tmp_bed_fo.write("\t".join(
                            [
                                record.INFO["CHR2"],
                                str(start),
                                str(record.INFO["END"] + 10000)
                            ]) + "\n")
                    tmp_bed_fo.close()

                    os.system(
                        "bedtools getfasta -bedOut -fi /home/ubuntu/seq/2022_02_WGS_LSK/reference/chm13v2.0.fa -bed tmp.bed > tmp.fa.bed")

                    sequences = []
                    for line in open("tmp.fa.bed"):
                        sequences.append(line.rstrip().split("\t")[-1])

                    fasta_fo = open("tmp.fa", "w")
                    fasta_fo.write(">ALT\n" + sequences[0] + str(insertion_in_region[0].ALT) + sequences[1] + "\n")
                    fasta_fo.write(">REF\n" + sequences[2] + "\n")
                    fasta_fo.write(">BND_WO_INS\n" + sequences[3] + "\n")
                    fasta_fo.close()

                    aligner = mp.Aligner("tmp.fa", preset="map-ont")

                    reads_supporting_alt = 0
                    reads_supporting_ref = 0

                    for read in sv_reads:
                        num_mappings_on_alt = 0
                        num_mappings_on_ref = 0
                        alignments_list = list(aligner.map(read.seq))
                        for aln in alignments_list:
                            if aln.mapq > 30:
                                if aln.ctg == "BND_WO_INS":
                                    num_mappings_on_ref += 1
                                if aln.ctg == "REF":
                                    num_mappings_on_ref += 1
                                if aln.ctg == "ALT":
                                    num_mappings_on_alt += 1

                        if num_mappings_on_ref > num_mappings_on_alt:
                            reads_supporting_ref += 1
                        if num_mappings_on_alt > num_mappings_on_ref:
                            reads_supporting_alt += 1

                    if reads_supporting_alt > reads_supporting_ref:
                        record.FILTER = ["INS"]

                elif len(insertion_in_region) == 0:
                    pass
                else:
                    print("ERROR")
                    print(record.ID)
                    break

                record.INFO["SUPP"] = len(sv_reads)

                if len(sv_reads) < 3:
                    record.FILTER.append("COV")

            vcf_writer.write_record(record)

    print(svtypes)



if __name__ == '__main__':

    samples = [
        ["sample_A", "path/bam_A.bam"],
        ["sample_B", "path/bam_B.bam"],
        ["sample_C", "path/bam_C.bam"],
        ["sample_D", "path/bam_D.bam"],
        ["sample_E", "path/bam_E.bam"]
    ]

    folder = "path/joint_calling/"

    for sample in samples:

        print("Filtering for length for sample: " + sample[0])

        os.system("bedtools sort -header -i " + folder + "intersect/" + sample_name + ".intersect.filtered.vcf > " + folder + "intersect/" + sample_name + ".intersect.filtered.sort.vcf")
        os.system("bgzip -c " + folder + "intersect/" + sample_name + ".intersect.filtered.sort.vcf > " + folder + "intersect/" + sample_name + ".intersect.filtered.vcf.gz")
        os.system("tabix -p vcf " + folder + "intersect/" + sample_name + ".intersect.filtered.vcf.gz")

        filter_intersect(folder + "intersect/" + sample_name + ".final.vcf", folder + "intersect/" + sample_name + ".final.filtered.vcf", 5000, folder + "intersect/" + sample_name + ".intersect.filtered.vcf.gz", sample[1])
