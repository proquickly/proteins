import os
import csv
import pickle
import pprint
import sqlite3
import pandas as pd
import functools
import time

DATA = "/Users/andy/data/proteins"
LIMIT = 500000  # float("inf")  # uncomment to do all

pd.set_option("display.max_rows", None)
pd.set_option("display.max_columns", None)
pd.set_option("display.width", None)
pd.set_option("display.max_colwidth", None)


intact_f = f"{DATA}/intact.txt"
biogrid_f = f"{DATA}/BIOGRID-ALL-4.4.222.tab3.txt"
db = sqlite3.connect(f"{DATA}/proteins.db")


def time_me(func):
    @functools.wraps(func)
    def wrapper_timer(*args, **kwargs):
        start_time = time.perf_counter()
        value = func(*args, **kwargs)
        end_time = time.perf_counter()
        run_time = end_time - start_time
        print(f"Ran {func.__name__!r} in {run_time:.4f} secs")
        return value

    return wrapper_timer


@time_me
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
            line_with_meta["type"] = "intact"
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
                dup_key = (key[0], key[1], line1["type"])
                if dup_key not in biogrid_dups:
                    biogrid_dups.append(dup_key)
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
                        match_key = (akey, bkey)
                        # print(key)
                        if key not in biogrid_dups:
                            biogrid_dups.append(key)
                            this_line["dup"] = False
                        else:
                            this_line["dup"] = True
                        this_line["id"] = progress
                        this_line["protein_key"] = str(match_key)
                        biogrid_processed.append(this_line)

            if progress > LIMIT:
                break

        return intact_processed, biogrid_processed


@time_me
def get_unmatched(intact, biogrid):
    return pd.concat([intact, biogrid]).drop_duplicates(keep=False)


@time_me
def get_matched(intact, biogrid):
    return pd.merge(
        intact, biogrid, how="inner", left_on="protein_key", right_on="protein_key"
    )


def read_files():
    with open(f"{DATA}/biogrid.pkl", "rb") as b:
        biogrid = pickle.load(b)
    with open(f"{DATA}/intact.pkl", "rb") as i:
        intact = pickle.load(i)
    return intact, biogrid


def display_menu():
    print(
        """
          1. Load raw files
          2. Match and list protein data
          9. Quit
          """
    )
    option = input()
    return option


def main():
    while True:
        option = display_menu()
        if option == "9":
            exit(0)
        if option == "1":
            intact, biogrid = load_raw_data()
            intact_dump = pd.DataFrame(intact)
            biogrid_dump = pd.DataFrame(biogrid)
            print("INTACT")
            print(intact_dump[["akey", "bkey", "dup", "id", "type", "protein_key"]])
            print("BIOGRID")
            print(biogrid_dump[["akey", "bkey", "dup", "id", "type", "protein_key"]])
            with open(f"{DATA}/intact.pkl", "wb") as i:
                pickle.dump(intact_dump, i)
            with open(f"{DATA}/biogrid.pkl", "wb") as b:
                pickle.dump(biogrid_dump, b)

        if option == "2":
            intact, biogrid = read_files()
            unmatched = get_unmatched(intact, biogrid)
            matched = get_matched(intact, biogrid)
            with open(f"{DATA}/unmatched.pkl", "wb") as u:
                pickle.dump(unmatched, u)
            with open(f"{DATA}/matched.pkl", "wb") as m:
                pickle.dump(biogrid, m)
            print("MATCHED")
            print(
                matched[
                    [
                        "protein_key",
                        "type_x",
                        "dup_x",
                        "id_x",
                        "type_y",
                        "dup_y",
                        "id_y",
                    ]
                ]
            )
            # print("UNMATCHED")
            # print(unmatched[["protein_key", "type", "dup", "id"]])
            print(f"Number matched is {len(matched.index)}")
            print(f"Number unmatched is {len(unmatched.index)}")
            print(
                f"% matched is {round(len(matched.index) / (len(matched.index) + len(unmatched.index)), 4)}"
            )


if __name__ == "__main__":
    main()
