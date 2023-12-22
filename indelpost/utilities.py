import re
import array

import numpy as np
from collections import namedtuple
from operator import mul
from functools import reduce
from pysam.libcbcf import VariantRecord, VariantRecordFilter, VariantFile

from indelpost.local_reference import UnsplicedLocalReference

cigar_ptrn = re.compile(r"[0-9]+[MIDNSHPX=]")


def most_common(lst):
    alst = list(set(lst))
    alst.sort()
    return max(alst, key=lst.count)


def get_gap_ptrn(read):
    return "".join([c for c in read["cigar_list"] if "D" in c or "I" in c])

def get_gap_ptrn2(read):
    ptrn = ""
    pos = read["aln_start"]
    for cigar in read["cigar_list"]:
        event, event_len = cigar[-1], int(cigar[:-1])
        if event in ("M", "X", "="):
            pos += event_len
        elif event in ["I", "D", "N"]:
            ptrn += "{}@{}".format(cigar, pos-1)
            if event == "D":
                pos += event_len
    return ptrn
        
def most_common_gap_pattern(targetpileup):
    ptrns = [get_gap_ptrn(read) for read in targetpileup]
    return most_common(ptrns)
 
def most_common_gap_ptrn(targetpileup):
    ptrns = [get_gap_ptrn2(read) for read in targetpileup]
    return most_common(ptrns)

def to_flat_list(lst_of_lst: list) -> list:
    return [i for lst in lst_of_lst for i in lst]


def to_flat_vcf_records(record: VariantRecord) -> list:
    
    alt = ""
    
    VcfRec = namedtuple(
        "VcfRec", "chrom pos id ref alt qual filter info format samples orig"
    )

    if not record.alts:
        return []
         
    flat_record = [
        VcfRec(
            chrom=record.chrom,
            pos=record.pos,
            id=record.id,
            ref=record.ref,
            alt=alt,
            qual=record.qual,
            filter=record.filter,
            info=record.info,
            format=record.format,
            samples=record.samples,
            orig=record,
        )
        for alt in record.alts
    ]
    
    return flat_record


def to_dict(record: object) -> dict:
    d = {}
    for k, v in record.items():
        if isinstance(v, tuple):
            d[k] = ",".join([str(i) for i in v])
        else:
            d[k] = v

    if d:
        return d


def match_indels(query: object, subject: object, matchby: str, indel_only: bool) -> bool:
    if matchby != "normalization" and indel_only and not query.is_indel:
        return False

    if matchby == "normalization":
        return query == subject

    elif matchby == "locus":
        if query.chrom != subject.chrom:
            return False

        query.normalize(inplace=True)
        subject.normalize(inplace=True)
        
        return query.pos == subject.pos

    elif matchby == "exact":
        return (
            (query.chrom == subject.chrom)
            and (query.pos == subject.pos)
            and (query.ref == subject.ref)
            and (query.alt == subject.alt)
        )
        

def linguistic_complexity(seq: str) -> float:
    n = len(seq)
    if n <= 1:
        return float(n)
    else:
        usage = []
        for i in range(1, n):
            i_mer = [seq[j : j + i] for j in range(n - i + 1)]  
            usage.append(len(set(i_mer)) / min(4 ** i, n - i + 1))
        
        return reduce(mul, usage)


def low_qual_fraction(pileup: list) -> float:
    pileup_vol = 1
    low_qual_vol = 0

    for read in pileup:
        pileup_vol += len(read["read_seq"])
        low_qual_vol += read["low_qual_base_num"]
    
    return low_qual_vol / pileup_vol


def to_minimal_repeat_unit(seq):
    """Find repeat unit in indel sequence
    """
    mid = int(len(seq) / 2)
    min_unit = seq
    
    j = 1
    found = False
    
    while j <= mid and not found:
        tandems = [seq[i : i + j] for i in range(0, len(seq), j)]
        if len(set(tandems)) == 1:
            found = True
            min_unit = list(set(tandems))[0]
        j += 1
    
    return min_unit


def repeat_counter(query_seq, flank_seq):
    """
    """
    qlen, flen = len(query_seq), len(flank_seq)
    count = 0 
    
    if flen < qlen:
        return count
    
    for i in range(0, flen, qlen):
        if flank_seq[i  : i + qlen] == query_seq:
            count += 1
        else:
            break
    
    return count




def get_mapped_subreads(cigarstring: str, aln_start_pos: int, aln_end_pos: int) -> list:
    
    cigar_lst = cigar_ptrn.findall(cigarstring)
    res = []
    
    current_pos = aln_start_pos
    for cigar in cigar_lst:
        event, event_len = cigar[-1], int(cigar[:-1])

        if event in ("M", "X", "="):
            res.append((current_pos, (current_pos + event_len - 1)))
            current_pos += event_len
        elif event in ("I", "S", "H", "P"):
            pass
        else:
            current_pos += event_len 
    
    return res


def get_end_pos(read_start_pos: int, lt_flank: str, cigarstring: str) -> int:
    read_start_pos -= 1

    cigar_lst = cigar_ptrn.findall(cigarstring)
    event_len, i = 0, 0
    flank_len = len(lt_flank)
    
    while flank_len > 0:
        cigar = cigar_lst[i]
        event, event_len = cigar[-1], int(cigar[:-1])
        
        if event == "D" or event == "N":
            read_start_pos += event_len
        elif event == "I":
            flank_len -= event_len
        elif event == "H" or event == "P":
            pass
        else: 
            flank_len -= event_len
            read_start_pos += event_len

        i += 1

    return read_start_pos + flank_len




def split_cigar(cigarstring: str, target_pos: int, start: int) -> tuple:
     
    cigar_lst = cigar_ptrn.findall(cigarstring)
    lt_lst = []
    rt_lst = cigar_lst

    start -= 1
    for cigar in cigar_lst:
        event, event_len = cigar[-1], int(cigar[:-1])
        
        move = 0 if event in ("I", "H", "P") else event_len
        start += move
        rt_lst = rt_lst[1 :]

        if target_pos <= start:
            diff = start - target_pos
            lt_cigar = str(event_len - diff) + event
            lt_lst.append(lt_cigar)
            
            if diff:
                    rt_lst = [str(diff) + event] + rt_lst
            
            return lt_lst, rt_lst
        else:
            lt_lst.append(cigar)


def merge_consecutive_gaps(cigar_lst):

    merged_lst = []
    while cigar_lst:
        c = cigar_lst[0]
        cigar_lst = cigar_lst[1:]

        if "I" in c or "D" in c:
            i = 0
            is_gap = True
            while i < len(cigar_lst) and is_gap:
                tmp = cigar_lst[i]
                is_gap = True if "I" in tmp or "D" in tmp else False
                i += 1

            if i - 1:
                c += "".join(cigar_lst[: i - 1])
                cigar_lst = cigar_lst[i - 1 :]

        merged_lst.append(c)

    return merged_lst


def make_insertion_first(cigarstring):
    
    cigar_lst = cigar_ptrn.findall(cigarstring)    
    
    merged_cigar_lst = merge_consecutive_gaps(cigar_lst)
    new_cigar = []
    for c in merged_cigar_lst:
        if "I" in c and "D" in c:
            c_lst = cigar_ptrn.findall(c)
            if "D" in c_lst[0]:
                swapped = c_lst[::-1]
                new_cigar.append("".join(swapped))
            else:
                new_cigar.append("".join(c_lst))
        else:
            new_cigar.append(c)
    
    return "".join(new_cigar)     


def relative_aln_pos(ref_seq, cigar_lst, aln_start, target_pos, include_clip=False):
    
    current_pos = aln_start - 1
    ref_seq_pos = 0
    for cigar in cigar_lst:
        event, event_len = cigar[-1], int(cigar[:-1])

        event = "M" if include_clip and event == "S" else event

        if event == "M" or event == "D":
            current_pos += event_len
            ref_seq_pos += event_len
        elif event in ("I", "H", "P"):
            pass
        else:
            current_pos += event_len
         
        if current_pos >= target_pos:
            break
    
    ref_seq_pos += (target_pos - current_pos)   
     
    return ref_seq_pos / len(ref_seq)


def split(data: object, cigarstring: str, target_pos: int, string_pos: int, is_for_ref: bool, reverse: bool) -> tuple:
    
    cigar_lst = cigar_ptrn.findall(cigarstring)
    _size = len(cigar_lst)

    
    data_moves = [0] * _size
    genome_moves = [0] * _size
    
    i = 0; j = 0
    
    for cigar in cigar_lst:
        event, event_len = cigar[-1], int(cigar[:-1])
        
        if event == "N":
            d_move = 0
            g_move = event_len
        elif event == "I":
            g_move = 0
            d_move = 0 if is_for_ref else event_len
        elif event == "D":
            g_move = event_len
            d_move = event_len if is_for_ref else 0
        elif event == "H" or event == "P":
            d_move = 0
            g_move = 0    
        else:
            g_move, d_move = event_len, event_len
        
        data_moves[i] = d_move
        genome_moves[i] = g_move
        i += 1

    if reverse:
        string_pos += 1
        data = data[::-1]
        data_moves = data_moves[::-1]
        genome_moves = genome_moves[::-1]
    else:
        string_pos -= 1

    for d_move, g_move in zip(data_moves, genome_moves):
        if reverse:
            if target_pos < string_pos:
                string_pos -= g_move
            else:
                break
        else:
            if string_pos < target_pos:
                string_pos += g_move
            else:
                break
        j += d_move
     
    diff = string_pos - (target_pos + 1)if reverse else target_pos - string_pos
    if reverse:
        lt = data[j + diff :]
        lt = lt[::-1]
        rt = data[: j + diff]
        rt = rt[::-1]
    else:
        lt = data[: j + diff]
        rt = data[j + diff :]
     
    return lt, rt


def get_local_reference(
    target: object,
    pileup: list,
    window: int,
    unspl_loc_ref: UnsplicedLocalReference,
    unspliced: bool = False,
    splice_pattern_only: bool = False,
) -> tuple:

    span = ""

    chrom, pos, reference = target.chrom, target.pos, target.reference
    
    if unspliced:
        splice_patterns = None
    else:
        splice_patterns = [read["splice_pattern"] for read in pileup if read["splice_pattern"] != ("", "")]
   
    ref_len = reference.get_reference_length(chrom)
    
    spl_ptrn = []
    
    if splice_patterns:
        lt_patterns = [ptrn[0] for ptrn in splice_patterns if ptrn[0]]
        if lt_patterns:
            lt_pattern = most_common(lt_patterns)
            lt_spl_pos = []
            for span in lt_pattern.split(":"):
                lt_spl_pos += [int(i) for i in span.split("-")] 
        else:
            lt_spl_pos = []

        rt_patterns = [ptrn[1] for ptrn in splice_patterns if ptrn[1]]
        if rt_patterns:
            rt_pattern = most_common(rt_patterns)
            rt_spl_pos = []
            for span in rt_pattern.split(":"):
                rt_spl_pos += [int(i) for i in span.split("-")]
        else:
            rt_spl_pos = []

        spl_pos = lt_spl_pos + rt_spl_pos
        last_idx = len(spl_pos) - 1
        
        left_len = 0
        first_pass = False
        local_reference = ""
        for i, x in enumerate(spl_pos):
            if i == 0:
                lt_end = max(0, x - window * 2)
                local_reference += reference.fetch(chrom, lt_end, x - 1)
                rt_end = x - 1
                if x + 1 < rt_end:
                    spl_ptrn.append((x + 1, rt_end))
                else:
                    spl_ptrn.append((lt_end, rt_end))
            elif i % 2 == 1 and i != last_idx:
                local_reference += reference.fetch(chrom, x, spl_pos[i+1] - 1)
                rt_end = spl_pos[i+1] - 1
                spl_ptrn.append((x + 1, rt_end))
            elif i % 2 == 0:
                pass
            elif i == last_idx:
                rt_end = min(x + window * 2, ref_len)
                local_reference += reference.fetch(chrom, x, rt_end)
                spl_ptrn.append((x + 1, rt_end))

            if pos <= rt_end and not first_pass:
                left_len = len(local_reference) - (rt_end - pos)
                first_pass = True

    else:
        local_reference = unspl_loc_ref.fetch_ref_seq(pos, window)
        #left_len = unspl_loc_ref.left_len 
        #local_reference = reference.fetch(chrom, max(0, pos - window * 3), min(pos + window * 3, ref_len))
        left_len = pos - max(0, pos - window * 3)

    if splice_pattern_only:
        return tuple(spl_ptrn)

    return local_reference, left_len    
