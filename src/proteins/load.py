import os
import csv
import pickle
import pprint
import sqlite3
import pandas as pd

pd.set_option("display.max_rows", None)
pd.set_option("display.max_columns", None)
pd.set_option("display.width", None)
pd.set_option("display.max_colwidth", None)

DATA = "/Users/andy/data/proteins"
LIMIT = 100
intact_f = f"{DATA}/intact.txt"
biogrid_f = f"{DATA}/BIOGRID-ALL-4.4.222.tab3.txt"
db = sqlite3.connect(f"{DATA}/proteins.db")


def load_raw_data():
    intact_dups = []
    intact_processed = []
    with open(intact_f, "r") as i_f:
        intact = csv.reader(i_f, delimiter=",")
        original_count = 0
        headings = next(intact)
        progress = 0
        for line in intact:
            progress += 1
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
            key = tuple(sorted((a, b)))
            line_with_meta["akey"] = key[0]
            line_with_meta["bkey"] = key[1]
            line_with_meta["type"] = ""
            if key not in intact_dups:
                intact_dups.append(key)
                line_with_meta["dup"] = False
            else:
                line_with_meta["dup"] = True
            line_with_meta["id"] = progress
            line_with_meta["protein_key"] = str(key)
            intact_processed.append(line_with_meta)
            if progress > LIMIT:
                break

            # pprint.pprint(intact_processed[line_with_meta["key"]])
    biogrid_dups = []
    progress = 0
    biogrid_processed = []
    with open(biogrid_f, "r") as b_f:
        biogrid = csv.DictReader(b_f, delimiter="\t", quoting=csv.QUOTE_NONE)
        for line in biogrid:
            line1 = line.copy()
            line2 = line.copy()
            progress += 1
            sa = line1["SWISS-PROT Accessions Interactor A"].replace("-", "")
            sa = sa.strip()
            ta = line1["TREMBL Accessions Interactor A"].replace("-", "")
            ta = ta.strip()
            sb = line1["SWISS-PROT Accessions Interactor B"].replace("-", "")
            sb = sb.strip()
            tb = line1["TREMBL Accessions Interactor B"].replace("-", "")
            sb = sb.strip()
            # print(f"{sa=}, {sb=}, {ta=}, {tb=}")
            if sa + sb != "":
                key = tuple(sorted((sa, sb)))
                # print(key)
                line1["akey"] = key[0]
                line1["bkey"] = key[1]
                line1["type"] = "Swiss-Prot"
                line1["id"] = progress
                key = (key[0], key[1], line1["type"])
                if key not in biogrid_dups:
                    biogrid_dups.append(key)
                    line1["dup"] = False
                else:
                    line1["dup"] = True
                line1["id"] = progress
                line1["protein_key"] = str(key)
                # print(line1, "\n")
                biogrid_processed.append(line1)
            if ta + tb != "":
                key = tuple(sorted((ta, tb)))
                # print(key)
                line2["akey"] = key[0]
                line2["bkey"] = key[1]
                line2["type"] = "Trembl"
                akeys = line2["akey"].split("|")
                bkeys = line2["bkey"].split("|")
                for akey in akeys:
                    for bkey in bkeys:
                        this_line = line2.copy()
                        # key = (key[0], key[1], line2["type"])
                        key = (akey, bkey, this_line["type"])
                        # print(key)
                        if key not in biogrid_dups:
                            biogrid_dups.append(key)
                            this_line["dup"] = False
                        else:
                            this_line["dup"] = True
                        this_line["id"] = progress
                        this_line["protein_key"] = str(key)
                        biogrid_processed.append(this_line)
            # pprint.pprint(biogrid_processed)


            if progress > LIMIT:
                break
        intact = pd.DataFrame(intact_processed)
        biogrid = pd.DataFrame(biogrid_processed)
        print("INTACT")
        print(intact[["akey", "bkey", "dup", "id", "type", "protein_key"]])
        print("BIOGRID")
        print(biogrid[["akey", "bkey", "dup", "id", "type", "protein_key"]])
        print(biogrid.info())
        return intact_processed, biogrid_processed


def match():
    """not needed like this??"""
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
    biogrid_matches[line["id"]] = line


def write_pickles():
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


def display_menu():
    print(
        """
          1. Load raw files
          2. Match protein data
          3. List match data
          4. Save match data
          5. Load match data
          9. Quit
          option = input()
          """
    )
    option = input()
    return option


def main():
    option = display_menu()
    if option == "1":
        intact, biogrid = load_raw_data()
    exit()


if __name__ == "__main__":
    main()
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
