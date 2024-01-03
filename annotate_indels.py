#!/usr/bin/env python3
# coding: utf8
# Written By Schaudge King
# Date: 02-05-2023
import pysam
from indelpost import Variant
from indelpost.varaln import VariantAlignment


def main(extended_window_size: int = 50):

    reference = pysam.FastaFile("/yunying/ref/human/b37/b37_Homo_sapiens_assembly19.fasta")
    bam = pysam.AlignmentFile("test/kit_complex_indel.bam")

    #v = Variant("7", 55242467, "AATTAAGAGAAGCAACATCTC", "TCT", reference)
    #valn = VariantAlignment(v, bam, exclude_duplicates=False, downsample_threshold=2000000, base_quality_threshold=5)

    # 1. Count reads:
    #cnt = valn.count_alleles(three_class=True)
    #print(cnt)
    # (4, 0, 4) as (reference, non_reference_non_target, target)

    v2 = Variant("4", 55593611, "TGTTGAGGAGATAAATGGAAACAATTATGTTTACATAGACCCAACACAACTTCCTTATGATC", "CCCTG", reference)
    valn2 = VariantAlignment(v2, bam, window=extended_window_size, exclude_duplicates=False,
                             downsample_threshold=2000000, base_quality_threshold=5)
    cnt2 = valn2.count_alleles(three_class=True)
    print(cnt2)

    # 1.1 forward and reverse
    # fw_rv_cnt = valn.count_alleles(fwrv=True)
    # print(fw_rv_cnt)
    # ((1, 3), (3, 1)) as ((non-target_fw, non-target_rv), (target_fw, target_rv))

    # 1.2 count by fragment
    # f_cnt = valn.count_alleles(by_fragment=True)
    # print(f_cnt)
    # (4, 3)
    # non-supporting fragments = C, D, F, G
    # supporting fragments = A, B, D

    # 2. Fetch reads as a list of pysam’s AlignedSegment:
    # supporting_reads = valn.fetch_reads()
    # non - supporting_reads = valn.fetch_reads(how="non_target")
    # reads_covering_the_locus = valn.fetch_reads(how="covering")

    # for read in supporting_reads:
    #     print(read.mapping_quality)


if __name__ == '__main__':
    from sys import argv
    main(int(argv[1]) if len(argv) > 1 else 50)
