import vcf
import numpy as np

def vcf_subtract (vcf_file_path, vcf_blacklist_file_path, subtracted_vcf_file_path):

    vcf_reader = vcf.Reader(filename=vcf_file_path)
    vcf_writer = vcf.Writer(open(subtracted_vcf_file_path, "w"), vcf_reader)

    n_total = 0
    n_passing = 0

    blacklist_reader = vcf.Reader(filename=vcf_blacklist_file_path)

    for record in vcf_reader:

        if record.CHROM != "chrM":

            n_total += 1

            record_blacklisted = False

            for blacklist_record in blacklist_reader.fetch(record.CHROM, record.POS - 20, record.POS + 20):

                dist_start = 1000000
                dist_end = 1000000

                if record.CHROM == blacklist_record.CHROM:
                    dist_start = np.abs(record.POS - blacklist_record.POS)

                if record.INFO["CHR2"] == blacklist_record.INFO["CHR2"]:
                    dist_end = np.abs(record.INFO["END"] - blacklist_record.INFO["END"])

                if dist_start <= 5 and dist_end <= 5:
                    record_blacklisted = True

            if not record_blacklisted:
                vcf_writer.write_record(record)
                n_passing += 1

    print("Total variants:\t\t" + str(n_total))
    print("Passing variants:\t" + str(n_passing))
