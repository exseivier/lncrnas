#!/usr/bin/env python

from sys import argv, exit

def load_mef(fin):
    """(STR) -> HASH[HASH]
    Requires the name of the file generated from RNAfold (Vienna) analysis.
    Loads the sequence and structure for every geneID into a hash structure.
    """
    mef = {}
    FIN = open(fin, "r")
    for line in FIN:
        line = line.strip("\n")
        if line[0] == ">":
            head = line
            mef[head] = {}
        elif line[0] in ["A", "T", "U", "G", "C"]:
            mef[head]["sequence"] = line
        elif line[0] in ["|", "(", ")"]:
            mef[head]["structure"] = line
        else:
            mef[head]["unknown"] = "[ERROR!] - Unknown string!"

    FIN.close()

    return mef

def _augc_percent(sequence):
    """(STR) -> FLOAT

    """
    seq_len = len(seq)
    counts = {  "A":0.0,
                "C":0.0,
                "G":0.0,
                "T":0.0}

    for nt in seq:
        if nt == "U":
            nt = "T"
        
        counts[nt] += 1.0

    augc_ratio = (counts["A"] + counts["T"]) / (counts["G"] + counts["C"])
    au_percent = (counts["A"] + counts["T"]) / float(seq_len)
    gc_percent = (counts["G"] + counts["C"]) / float(seq_len)
    return [augc_ratio, au_percent, gc_percent]

def calculate_augc_percent(mef):
    """(HASH) -> HASH

    """
    augc_hash = {}
    for head, items in mef.iteritems():
        # p_augc is an array where augc_ratio, au_percent and gc_percent were alllocated
        p_augc = _augc_percent(items["sequence"])
        # Allocating p_augc to augc_hash table
        augc_hash[head] = p_augc

    return augc_hash

def print_augc(augc_percent):
    """(HASH) -> STDOUT

    """
    for head, percent in augc_percent.iteritems():
        print head\
                + "\t"\
                + "\t".join([str(i) for i in percent])



def main():
    """

    """
    if len(argv) != 2:
        print "[FATAL ERROR!] - Bad number of arguments!"
        exit(1)
    fin =   argv[1]
    mef = load_mef(fin)
    augc_percent = calculate_augc_percent(mef)
    print_augc(augc_percent)


