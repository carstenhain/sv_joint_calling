import os

def call_variants (sample_name, bam_file_path, folder, reference):

    os.system("sniffles \
        --input " + bam_file_path + " \
        --vcf " + folder + sample_name + ".sniffles2.2.vcf \
        --reference " + reference + " \
        --non-germline \
        --threads 25 \
        --minsupport 2 \
        --output-rnames")

    os.system("delly \
        lr \
        -g " + reference + " \
        -y ont \
        -o " + folder + sample_name + ".delly.bcf " + bam_file_path)

    os.system("bcftools convert -O v " + folder + sample_name + ".delly.bcf -o " + folder + sample_name + ".delly.vcf")

    os.system("rm -r " + folder + "svim_dir")
    os.system("mkdir " + folder + "svim_dir")

    os.system("svim alignment \
        --min_sv_size 30 \
        --max_sv_size 200000000 \
        --sample " + sample_name + " \
        --read_names \
        " + folder + "svim_dir " + bam_file_path + " " + reference)

    os.system("mv " + folder + "svim_dir/variants.vcf " + folder + sample_name + ".svim.vcf")

    os.system("rm -r " + folder + "cuteSV_dir")
    os.system("mkdir " + folder + "cuteSV_dir")

    os.system("cuteSV \
        --report_readid \
        --max_cluster_bias_INS 100 \
        --diff_ratio_merging_INS 0.3 \
        --max_cluster_bias_DEL 100 \
        --diff_ratio_merging_DEL 0.3 \
        --threads 24 \
        --min_size 30 \
        --min_support 2 \
        --max_size 200000000 " + bam_file_path + " " + reference + " " + folder + sample_name + ".cuteSV.2.vcf \
        " + folder + "cuteSV_dir")

if __name__ == '__main__':

    samples = [
        ["sample_A", "path/bam_A.bam"],
        ["sample_B", "path/bam_B.bam"],
        ["sample_C", "path/bam_C.bam"],
        ["sample_D", "path/bam_D.bam"],
        ["sample_E", "path/bam_E.bam"]
    ]

    for sample in samples:

        call_variants(sample[0], sample[1], "path/joint_calling/", "reference/chm13v2.0.fa")
