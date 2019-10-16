import csv
import pprint
import os
import gzip
from collections import Counter
import json

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

def get_sites(r):
    recs = {}
    for row in r:
        sn = row[5]
        # primer = row[7]
        # counts = row[9]
        posid = row[11]

        if sn in recs:
            recs[sn][posid] = {"ids": []}
        else:
            recs[sn] = {}
            recs[sn][posid] = {"ids": []}

    return(recs)

def get_fq_ids(r, recs):
    ignored = []
    for row in r:
        sn = row[7]
        fq_id = row[8]
        # primer = row[9]
        posid = row[10]

        try:
            recs[sn][posid]["ids"].append(fq_id)
        except KeyError:
            ignored.append([sn, posid])

    return(recs)

def scan_fastq(recs, fastq_dir):
    for rep in recs:
        print("working on", rep)

        gz_fn = os.path.join(dirname, fastq_dir, + str(rep) + ".R2.fastq.gz")
        with gzip.open(gz_fn, "rt") as handle:
            fq = SeqIO.to_dict(SeqIO.parse(handle, "fastq"))

            for site in recs[rep]:
                if "U5" in rep:
                    end = 69
                else:
                    end = 62

                c = Counter()
                for i in recs[rep][site]["ids"]:
                    seq = fq[i].seq[0:end]
                    c.update({str(seq): 1})

                recs[rep][site]["seq"] = c.most_common(1)[0]

    with open(counted_json, "w") as output:
        json.dump(recs, output)

    return(recs)

def main():
    dirname = os.path.dirname(__file__)
    fn = os.path.join(dirname, "../data/raw_xofil_df_acute.csv")

    with open(fn) as xofil_csv:
        r = csv.reader(xofil_csv, delimiter=',')
        next(r)

        recs = get_sites(r)

    fn = os.path.join(dirname, "../data/comb_proc_unique_sites_acute.csv")
    with open(fn) as uniq_csv:
        r = csv.reader(uniq_csv, delimiter=',')
        next(r)

        recs = get_fq_ids(r, recs)

    counted_json = "ltr_mapped.json"
    if os.path.exists(counted_json):
        with open(counted_json, "r") as cj:
            recs = json.load(cj)
        
        for rep in recs:
            for site in recs[rep]:
                del recs[rep][site]["ids"]
        
        # pprint.pprint(recs)

    else:
        print("scanning fastq gz files and getting LTR sequences")
        recs = scan_fastq(recs, "../data/fastq_gz/")

    # convert dict to fasta file
    u3_fasta = []
    u5_fasta = []
    for rep in recs:
        for site in recs[rep]:
            rec = "_".join([rep, site])
            seq = Seq(recs[rep][site]['seq'][0], generic_dna)
            seqrec = SeqRecord(seq, id = rec, description = "")

            if "U3" in rep:
                u3_fasta.append(seqrec)
            else:
                u5_fasta.append(seqrec)

    with open("ltr_only_u3.fa", "w") as outfile:
        SeqIO.write(u3_fasta, outfile, "fasta")
    
    with open("ltr_only_u5.fa", "w") as outfile:
        SeqIO.write(u5_fasta, outfile, "fasta")

main()


