import matplotlib.pyplot as plt
import vcf
import numpy as np
import os
import vcf_subtract
import pysam
import mappy as mp

def filter_sniffles2 (original_sniffles2_file_path, filtered_sniffles2_file_path, dv_cutoff):
    fo = open(filtered_sniffles2_file_path, "w")
    for line in open(original_sniffles2_file_path, "r"):
        if line.startswith("#"):
            fo.write(line)
        else:
            if int(line.split("\t")[-1].split(":")[-1]) >= dv_cutoff:
                fo.write(line)
    fo.close()

def filter_svim (original_svim_file_path, filtered_svim_file_path, dv_cutoff):

    vcf_reader = vcf.Reader(open(original_svim_file_path, "r"))
    vcf_writer = vcf.Writer(open(filtered_svim_file_path, "w"), vcf_reader)

    for record in vcf_reader:
        if record.INFO["SUPPORT"] >= dv_cutoff:
            vcf_writer.write_record(record)

def filter_delly (original_delly_file_path, filtered_delly_file_path, dv_cutoff):

    fo = open(filtered_delly_file_path, "w")
    for line in open(original_delly_file_path, "r"):
        if line.startswith("#"):
            fo.write(line)
        else:

            sr = 0

            for info in line.split("\t")[7].split(";"):
                if info.startswith("SR="):
                    sr = int(info.replace("SR=", ""))

            dv_idx = -1
            rv_idx = -1
            for idx, gt_element in enumerate(line.rstrip().split("\t")[8].split(":")):
                if gt_element == "DV":
                    dv_idx = idx
                if gt_element == "RV":
                    rv_idx = idx

            rv = int(line.rstrip().split("\t")[9].split(":")[rv_idx])

            dv = np.amax([sr, rv])

            A = line.rstrip().split("\t")
            GT = A[-1].split(":")
            GT[dv_idx] = str(dv)

            if dv >= dv_cutoff:

                fo.write("\t".join(A[0:-1]) + "\t" + ":".join(GT) + "\n")

    fo.close()

def preprocess_cuteSV (original_cuteSV_file_path, preprocess_cuteSV_file_path):

    fo = open(preprocess_cuteSV_file_path, "w")
    for line in open(original_cuteSV_file_path, "r"):
        if line.startswith("#"):
            fo.write(line)
        else:
            A = line.rstrip().split("\t")
            if len(A[3]) > 500:
                A[3] = "N"
            fo.write("\t".join(A) + "\n")
            fo.flush()
    fo.close()

def filter_cuteSV (preprocess_cuteSV_file_path, filtered_cuteSV_file_path, dv_cutoff):
    vcf_reader = vcf.Reader(open(preprocess_cuteSV_file_path, "r"))
    vcf_writer = vcf.Writer(open(filtered_cuteSV_file_path, "w"), vcf_reader)

    for record in vcf_reader:
        if record.samples[0]["DV"] >= dv_cutoff:
            vcf_writer.write_record(record)

def filter_intersect (original_intersect_file_path, filtered_intersect_file_path, dv_cutoff):
    vcf_reader = vcf.Reader(open(original_intersect_file_path, "r"))
    vcf_writer = vcf.Writer(open(filtered_intersect_file_path, "w"), vcf_reader)

    for record in vcf_reader:
        max_dv = 0
        for sample in record.samples:
            if sample["DR"][1] > max_dv:
                max_dv = sample["DR"][1]
        if max_dv >= dv_cutoff:
            vcf_writer.write_record(record)

if __name__ == '__main__':

    samples = ["sample_A", sample_B", sample_C", sample_D", sample_E"]

    folder = "path/joint_calling/"

    os.system("mkdir " + folder + "union/")
    os.system("mkdir " + folder + "intersect/")

    survivor_path = "path/./SURVIVOR"

    dgv_sv_vcf_gz_file_path = "reference/DGV_SV_T2Tv2.vcf.gz"

    artifact_fo = open(folder + "artifact.txt", "w")

    for sample_name in samples:

        print("Filtering and merging into union and intersect for sample: " + sample_name)

        filter_sniffles2 (folder + sample_name + ".sniffles2.2.vcf", folder + "union/" + sample_name + ".sniffles2.2.vcf", 2)
        filter_svim (folder + sample_name + ".svim.vcf", folder + "union/" + sample_name + ".svim.2.vcf", 2)
        filter_delly (folder + sample_name + ".delly.vcf", folder + "union/" + sample_name + ".delly.2.vcf", 2)
        preprocess_cuteSV (folder + sample_name + ".cuteSV.2.vcf", folder + sample_name + ".cuteSV.2.processed.vcf")
        filter_cuteSV(folder + sample_name + ".cuteSV.2.processed.vcf", folder + "union/" + sample_name + ".cuteSV.2.filtered.vcf", 2)

        fo = open(folder + "union/" + sample_name + ".callers.txt", "w")
        fo.write(folder + "union/" + sample_name + ".sniffles2.2.vcf\n")
        fo.write(folder + "union/" + sample_name + ".svim.2.vcf\n")
        fo.write(folder + "union/" + sample_name + ".delly.2.vcf\n")
        fo.write(folder + "union/" + sample_name + ".cuteSV.2.filtered.vcf\n")
        fo.close()

        # union
        os.system(survivor_path + " merge " + folder + "union/" + sample_name + ".callers.txt 50 1 0 0 0 30 " + folder + "union/" + sample_name + ".union.vcf")
        artifact_fo.write(folder + "union/" + sample_name + ".union.vcf\n")

        # intersect
        os.system(survivor_path + " merge " + folder + "union/" + sample_name + ".callers.txt 50 2 0 0 0 30 " + folder + "intersect/" + sample_name + ".intersect.vcf")

        # filter for DV
        filter_intersect (folder + "intersect/" + sample_name + ".intersect.vcf", folder + "intersect/" + sample_name + ".intersect.filtered.vcf", 4)

    
    print("Preparing artifacts file for filtering")
    artifact_fo.close()
    os.system(survivor_path + " merge " + folder + "artifact.txt 50 3 0 0 0 30 " + folder + "union/artifacts.3.vcf")
    os.system("bedtools sort -header -i " + folder + "union/artifacts.3.vcf > " + folder + "union/artifacts.3.sort.vcf")
    os.system("bgzip -c " + folder + "/union/artifacts.3.sort.vcf > " + folder + "union/artifacts.3.vcf.gz")
    os.system("tabix -p vcf " + folder + "union/artifacts.3.vcf.gz")
    

    for sample_name in samples:

        print("Subtracting artifacts and DGV_SV from sample: " + sample_name)

        vcf_subtract.vcf_subtract(folder + "intersect/" + sample_name + ".intersect.filtered.vcf", folder + "union/artifacts.3.vcf.gz", "tmp.vcf")
        vcf_subtract.vcf_subtract("tmp.vcf", dgv_sv_vcf_gz_file_path, "tmp2.vcf")
        vcf_subtract.vcf_subtract("tmp2.vcf", folder + "union/19039.T2T.union.sort.vcf.gz" ,folder + "intersect/" + sample_name + ".final.vcf")
