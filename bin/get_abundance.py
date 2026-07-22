#!/usr/bin/env python

import pandas as pd
from functools import reduce
import requests
import json
import numpy as np
import argparse
import logging

logger = logging.getLogger()

# taxids that were classified but couldn't be resolved to a name from the taxonomy dump -- a sign
# the taxonomy (rankedlineage.dmp) doesn't match the reference databases. Collected and reported loudly.
UNRESOLVED_TAXIDS = set()

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--infile", help="Classification results file", type=str)
    parser.add_argument("--prefix", help="ID of the sample, usually barcode.", type=str)
    parser.add_argument("--outfile", help="Output file name.", type=str, default="rel_abundance")

    return parser.parse_args()

def get_taxname(tax_id, tax_level):
    # Offline-safe. Taxonomy names come from the local dmp (see get_taxname_from_dmp). This
    # previously queried api.unipept.ugent.be, which is fatal on an offline instrument; it now
    # never contacts the internet and just returns a sensible local fallback.
    if str(tax_id) == "nan":
        return "unclassified"
    try:
        return str(int(tax_id))
    except (ValueError, TypeError):
        return "unclassified"


def get_taxname_from_dmp(data, tax_id, tax_level):
    tags = {"S": "species","G": "genus","F": "family", "O": "order"}
    tax_level_tag = tags[tax_level]

    if str(tax_id) == "nan":
        return 'unclassified'
    match = data.loc[data['taxid'] == tax_id]
    if match.empty:
        return 'unclassified'   # taxid not in the classification table (e.g. root/unresolved)
    name = match[tax_level_tag].iloc[0]
    if not isinstance(name, str):
        name = match["name"].iloc[0] if "name" in match.columns else None
    if not isinstance(name, str):
        name = match["sciname"].iloc[0] if "sciname" in match.columns else None
    if isinstance(name, str):
        return name
    # The taxid was classified but has no lineage in the supplied taxonomy: rankedlineage.dmp
    # doesn't match the reference databases. Record it so main() can warn loudly.
    UNRESOLVED_TAXIDS.add(str(tax_id))
    return str(tax_id)


def get_abundance_values(names,paths):
    dfs = []
    for name,path in zip(names,paths):
        data1 = pd.read_csv(path, index_col=False, sep=';').iloc[:,1:]

        total = sum(data1['reads_in_cluster'])
        rel_abundance=[]

        data=choose_classification(data1)

        for index,row in data.iterrows():
            rel_abundance.append(row['reads_in_cluster'] / total * 100)
            
        data['rel_abundance'] = rel_abundance
        dfs.append(pd.DataFrame({'taxid': data['taxid'], 'rel_abundance': rel_abundance, 'reads': data['reads_in_cluster']}))
        data.to_csv("" + name + "_nanoclust_out.txt")

    return dfs, data

def choose_classification(dataframe):
    print(dataframe)
    if len(dataframe.columns)>13:
        chosen_frame=[]
        classification_score={}
        for index, row in dataframe.iterrows():
            print(row['class_level'])
            if row['class_level']=="S":
                chosen_frame.append(row.iloc[:12].tolist())
            else:
                classification_score["kraken2"]=sum(row.notna()[8:12])
                classification_score["blast"]=sum(row.notna()[24:])
                classification_score["seqmatch"]=sum(row.notna()[16:20])
                choice=max(classification_score, key=classification_score.get)
                
                if choice == "kraken2":
                    chosen_frame.append(row.iloc[:12].tolist())
                elif choice == "seqmatch":
                    chosen_frame.append(row.iloc[np.r_[0:4,13:20]].tolist())
                else:
                    chosen_frame.append(row.iloc[np.r_[0:4,20:28]].tolist())
            
        logger.info("Choosing classification")
        logger.debug(chosen_frame)

        chosen_df=pd.DataFrame(chosen_frame, columns=['reads_in_cluster', 'used_for_consensus', 'reads_after_corr', 'draft_id', 'classifier_name', 'taxid', 'stat', 'name', 'species', 'genus', 'family', 'order'])
        logger.debug(len(chosen_df))
        logger.debug(chosen_df)

        return chosen_df
    else:
        return dataframe

def choose_row_classifier(row):
    """Return which classifier won for a single cluster row of the full-mode table.

    Mirrors exactly the per-row decision in choose_classification so the report can be
    filtered to the same winning classifier that produced the reported abundance. The row
    is a positional Series from the id-dropped full-mode dataframe. Column layout:
    0-3 cluster metadata, 4-11 kraken2 block, 12-19 seqmatch block, 20-27 blast block
    (position 6 is the kraken2 class_level; positions 8:12 / 16:20 / 24: are the
    species/genus/family/order columns of each classifier).
    """
    if row.iloc[6] == "S":
        return "kraken2"
    classification_score = {
        "kraken2": sum(row.notna()[8:12]),
        "blast": sum(row.notna()[24:]),
        "seqmatch": sum(row.notna()[16:20]),
    }
    return max(classification_score, key=classification_score.get)


def write_chosen_classifier(infile, prefix):
    """Write a cluster -> winning classifier mapping used to filter the report hit details.

    Only meaningful in 'full' mode, where three classifiers compete per cluster. In the
    single-classifier modes the file is written with a header only (the report then falls
    back to the sole classifier present for each cluster).
    """
    raw = pd.read_csv(infile, index_col=False, sep=';')
    rows = []
    # id column + >13 classifier columns == full mode (matches choose_classification)
    full = raw.shape[1] > 14
    if raw.shape[1] >= 2:
        ids = raw.iloc[:, 0]
        reads_col = pd.to_numeric(raw.iloc[:, 1], errors="coerce").fillna(0)
        total = reads_col.sum() or 1
        data1 = raw.iloc[:, 1:]
        for i, (_, row) in enumerate(data1.iterrows()):
            reads = int(reads_col.iloc[i])
            rows.append({
                "cluster": ids.iloc[i],
                "classifier": choose_row_classifier(row) if full else "",
                "reads": reads,
                "rel_abundance": round(reads / total * 100, 1),
            })
    pd.DataFrame(rows, columns=["cluster", "classifier", "reads", "rel_abundance"]).to_csv(
        prefix + "_chosen_classifier.csv", index=False)


def merge_abundance(dfs, data, tax_level):
    df_final = reduce(lambda left, right: pd.merge(left, right, on='taxid', how='outer').fillna(0), dfs)
    all_tax = []

    for index, row in df_final.iterrows():
        try:
            # Handle collapsing for specific tax IDs at 'S' level
            if tax_level == "S" and row["taxid"] in [1280, 985002, 1654388]:
                all_tax.append("Staphylococcus aureus complex")
            else:
                all_tax.append(get_taxname_from_dmp(data, row["taxid"], tax_level))
        except:
            logger.error("Error getting taxonomic name for tax_id {} in merge_abundance.".format(row["taxid"]))
            if tax_level == "S" and row["taxid"] in [1280, 985002, 1654388]:
                all_tax.append("Staphylococcus aureus complex")
            else:
                all_tax.append("unclassified")   # offline: never contact the taxonomy API

    df_final["taxid"] = all_tax

    # Collapse entries labeled as 'Staphylococcus aureus complex'
    if tax_level == "S":
        df_final.loc[df_final["taxid"] == "Staphylococcus aureus complex", "taxid"] = "Staphylococcus aureus complex"

    # Group by taxid and sum the abundance values
    df_final_grp = df_final.groupby(["taxid"], as_index=False).sum()
    df_final_sorted = df_final_grp.sort_values(by='rel_abundance', ascending=False)

    return df_final_sorted



def get_abundance(names,paths,tax_level, outfile):
    if(not isinstance(paths, list)):
        paths = [paths]
        names = [names]

    dfs, data = get_abundance_values(names,paths)
    df_final_grp = merge_abundance(dfs, data, tax_level)
    df_final_grp.to_csv(outfile + "_"+ names[0] + "_" + tax_level + ".csv", index = False)


def warn_unresolved_taxids(prefix):
    """Loudly flag any classified taxid that the supplied taxonomy couldn't name.

    An empty result is the happy path (nothing written, nothing printed). When taxids are
    unresolved it means rankedlineage.dmp is a different vintage from the reference
    databases, so a real organism can be silently reported as a bare number. We shout about
    it in the log AND drop a published file so it can't be missed.
    """
    if not UNRESOLVED_TAXIDS:
        return
    ids = ", ".join(sorted(UNRESOLVED_TAXIDS, key=lambda x: (len(x), x)))
    bar = "!" * 78
    msg = (
        "\n{bar}\n"
        "!! TAXONOMY MISMATCH: {n} classified taxid(s) are missing from the taxonomy.\n"
        "!! Unresolved taxids: {ids}\n"
        "!! These organisms were identified but could NOT be named, so they appear in the\n"
        "!! results as a bare number instead of a species. This means the taxonomy\n"
        "!! (rankedlineage.dmp) does not match the reference databases. Rebuild the\n"
        "!! taxonomy from the SAME NCBI snapshot as the BLAST/kraken2 databases.\n"
        "{bar}\n"
    ).format(bar=bar, n=len(UNRESOLVED_TAXIDS), ids=ids)
    # stderr so it surfaces in .command.err / the Nextflow log even when stdout is captured
    import sys
    sys.stderr.write(msg)
    print(msg)
    with open(prefix + "_unresolved_taxids.txt", "w") as fh:
        fh.write("Taxids that were classified but are missing from the supplied taxonomy "
                 "(rankedlineage.dmp).\n")
        fh.write("The taxonomy does not match the reference databases; rebuild them as a "
                 "matched set from the same NCBI snapshot.\n\n")
        fh.write("\n".join(sorted(UNRESOLVED_TAXIDS, key=lambda x: (len(x), x))) + "\n")


def main(args):

    write_chosen_classifier(args.infile, args.prefix)

    for level in ["G", "S", "O", "F"]:
        get_abundance(args.prefix, args.infile, level, args.outfile)

    warn_unresolved_taxids(args.prefix)

if __name__=="__main__":
    args = parse_args()

    main(args)