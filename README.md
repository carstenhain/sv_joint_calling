# sv_joint_calling
Calling of structural variants (SVs) in multiple samples with multiple callers and subsequent combination and filtering

This is still work in progress and the functionality is currently split into three hard-coded scripts.
In the future this workflow will be combined and sample input via files will be added.

## jc_calling.py

Calls SVs using Sniffles2 (v.2.0.7, Smolka et al., 2022), Delly (v.1.1.6, Rausch et al., 2012) SVIM (v.1.4.2, Heller and Vingrom, 2019) and cuteSV (v.1.0.8, Jiang et al., 2020).

| Tool | Parameters |
|---|---|
| Sniffles2 | --non-germline; --minsupport 2 |
| Delly | no change |
| SVIM | --max_sv_size 200000000 |
| cuteSV | --min_support 2; max_size 200000000 |

## jc_jointing.py

Combines the calls from single callers into a union of all found SVs and an intersection of SVs found in at least 2 callers.
Build a file containing false positive SVs/artifacts and common potential germline SVs.
Removes the SVs from the artifacts file and from the database of genomic variants (MacDonald et al., 2014) from the single sample SV calls. 
Uses SURVIVOR (v.1.0.7, Jeffares et al., 2017)

1. Removes sequence from SVs longer than 500 bases from the cuteSV calls
2. Removes SV calls with less than 2 supporting reads
3. Union of all SVs detected with different callers in a single sample with SURVIVOR
4. Intersection of SVs detected with at least 2 callers in a single sample with SURVIVOR
5. SVs that do not have more than 3 supporting reads with any caller are removed
7. Creation of the artifact and common SV file as intersection of all single samples unions, SVs present in at least 3 samples are added to the artifact file
8. Subtraction of artifacts and common SVs from the single sample intersection calls
    1. Subtraction of all SVs in the artifact file created in step 7
    2. Subtraction of all SVs listed in the 

| Step | SURVIVOR parameters |
|---|---|
| Single sample union | 50 1 0 0 0 30 |
| Single sample intersection | 50 2 0 0 0 30 |
| Artifact intersection | 50 3 0 0 0 30 |

## jc_filtering.py

Removes SVs which results from false mapping of reads adjacent to insertion

1. Detection of close occurance (<150 bp) of a SVs (TRA, DUP, DEL, INV) and a insertion
2. Extractiong of
    1. reference sequences at SV start and end
    2. insertion sequence
    3. SV supporting reads
3. Building of an insertion reference
4. Mapping of the SV supporting reads to a multiple fasta file containing the normal references at SV start and end and potential insertion references
5. Counting of number of reads mapping to the normal reference and to the insertion reference
6. If more reads map to the insertion reference than the normal reference the SV is filtered
7. Filtering of SVs smaller than 5 kb
