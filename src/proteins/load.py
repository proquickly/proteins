import os
import csv
import pickle
import pprint
import sqlite3


DATA = "/Users/andy/data/proteins"
LIMIT = 5000
intact_f = f"{DATA}/intact.txt"
biogrid_f = f"{DATA}/BIOGRID-ALL-4.4.222.tab3.txt"
db = sqlite3.connect(f"{DATA}/proteins.db")


def match():
    intact_dups = []
    intact_processed = {}
    with open(intact_f, "r") as i_f:
        intact = csv.reader(i_f, delimiter=",")
        original_count = 0
        headings = next(intact)
        progress = 0
        for line in intact:
            progress += 1
            if progress > LIMIT:
                break
            line_with_meta = dict(zip(headings[0].split("\t"), line[0].split("\t")))
            original_count += 1
            a = line_with_meta["#ID(s) interactor A"]
            a = a.replace('"', "")
            a = a.replace("uniprotkb:", "")
            a = a.replace("chebi:CHEBI:", "")
            a = a.replace("intact:EBI-", "")
            a = a.replace("ensembl:ENSG00000", "")
            b = line_with_meta["ID(s) interactor B"]
            b = b.replace('"', "")
            b = b.replace("uniprotkb:", "")
            b = b.replace("ensembl:ENSG00000", "")
            b = b.replace("intact:EBI-", "")
            line_with_meta["key"] = (a, b)
            if (a, b) not in intact_dups:
                intact_dups.append((a, b))
                line_with_meta["dup"] = False
            else:
                line_with_meta["dup"] = True
            intact_processed[line_with_meta["key"]] = line_with_meta
            # pprint.pprint(intact_processed[line_with_meta["key"]])
    duplicates = []
    biogrid_matches = {}
    match_number = 0
    progress = 0
    with open(biogrid_f, "r") as b_f:
        biogrid = csv.DictReader(b_f, delimiter="\t", quoting=csv.QUOTE_NONE)
        for line in biogrid:
            progress += 1
            """if progress > LIMIT:
                exit()
                break"""
            sa = line["SWISS-PROT Accessions Interactor A"].replace("-", "")
            ta = line["TREMBL Accessions Interactor A"].replace("-", "")
            sb = line["SWISS-PROT Accessions Interactor B"].replace("-", "")
            tb = line["TREMBL Accessions Interactor B"].replace("-", "")
            # print(sa, sb, ta, tb)
            # print(line["key"][0])
            for a in line["key"][0]:
                for a_i in a:
                    for b in line["key"][1]:
                        for b_i in b:
                            # print(a_i.strip(), b_i.strip())
                            if (a_i, b_i) in intact_processed.keys():
                                print(f"MATCH {match_number}")
                                match_number += 1
                                line["match"] = match_number
                                # intact_processed["key"]["match"] = match_number
                                # print(f"{intact_processed['key']['match']}")
                # print(f"{intact_processed['key']['match']}")
            # continue
            biogrid_matches[line["key"]] = line
            # pprint.pprint(biogrid_matches[line["key"]])

        with open("intact_matches.pkl", "wb") as im:
            pickle.dump(intact_processed, im)

        with open("biogrid_matches.pkl", "wb") as bm:
            pickle.dump(biogrid_matches, bm)
        print("=" * 80)
        count = 0
        for item in intact_processed.items():
            count += 1
            pprint.pprint(item)
            print("-" * 80)
            """if count > 20:
                break"""
        print(f"Intact matches {count}")
        count = 0
        for item in biogrid_matches.items():
            count += 1
            pprint.pprint(item)
            print("-" * 80)
            """if count > 20:
                break"""
        print(f"Biogrid matches {count}")

    """add control total for all catgories"""


def read_files():
    with open("biogrid_matches.pkl", "rb") as bm:
        biogrid_matches = pickle.load(bm)
    with open("intact_matches.pkl", "rb") as im:
        intact_matches = pickle.load(im)
    return intact_matches, biogrid_matches


import pandas as pd

if __name__ == "__main__":
    # intact, biogrid = read_files()
    intact = pd.read_pickle("intact_matches.pkl")
    with open("intact_match.txt", "w") as i:
        count = 0
        for k, v in intact.items():
            i.write(str(k))
            i.write("\n")
            for k2, v2 in v.items():
                i.write(str(k2) + ", " + str(v2))
                i.write("\n")
            count += 1
            if count > 100:
                break
            i.write("\n\n")

    biogrid = pd.read_pickle("biogrid_matches.pkl")
    with open("biogrid_match.txt", "w") as i:
        count = 0
        for k, v in biogrid.items():
            i.write(str(k))
            i.write("\n")
            for k2, v2 in v.items():
                i.write(str(k2) + ", " + str(v2))
                i.write("\n")
            count += 1
            if count > 100:
                break
            i.write("\n\n")
