import csv
import pprint
import os
import gzip
from collections import Counter

from Bio import SeqIO

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

    for rep in recs:
        print(rep)

        gz_fn = os.path.join(dirname, "../data/fastq_gz/" + str(rep) + ".R2.filt.fastq.gz")
        with gzip.open(gz_fn, "rt") as handle:
            fq = SeqIO.to_dict(SeqIO.parse(handle, "fastq"))

            for site in recs[rep]:
                c = Counter()
                for i in recs[rep][site]["ids"]:
                    seq = fq[i].seq
                    c.update({str(seq): 1})

                pprint.pprint(c)
        break

main()


