#!/usr/bin/env python

"""Create results report."""

import argparse
from aplanat import report
import pandas as pd
import os
import glob
import re
import logging
import numpy as np
from bokeh.resources import INLINE

logger = logging.getLogger(__name__)

# Number of hits to show per cluster per classifier
TOP_N_HITS = 3

# Marker appended to the classifier's selected (rank 1 / LCA) call
SELECTED_MARK = " *"

# Map the classifier names written by get_abundance.py to the report table labels
CLASSIFIER_LABELS = {
    "blast": "BLAST",
    "seqmatch": "SeqMatch",
    "kraken2": "Kraken2 (LCA)",
}

# Short classifier names used in the hit-details CSV, keyed by the filename token
# (detect_classifier) and by the get_abundance chosen-classifier token
HIT_CSV_CLASSIFIER = {"blastn": "BLAST", "seqmatch": "SeqMatch", "kraken2": "Kraken2"}
CHOSEN_CSV_CLASSIFIER = {"blast": "BLAST", "seqmatch": "SeqMatch", "kraken2": "Kraken2"}


def load_chosen_classifiers(path):
    """Load the cluster -> winning classifier map produced by get_abundance.py.

    Returns {cluster_id (str): report label}. Empty dict if unavailable (e.g. single
    classifier mode), in which case the report shows whatever hits exist per cluster.
    """
    if not path or path in ("none", "unknown") or not os.path.isfile(path):
        return {}
    try:
        df = pd.read_csv(path, dtype=str)
    except Exception:
        logger.exception("Failed to read chosen classifier map %s", path)
        return {}
    mapping = {}
    if "cluster" in df.columns and "classifier" in df.columns:
        for _, r in df.fillna("").iterrows():
            cid = str(r["cluster"]).strip()
            label = CLASSIFIER_LABELS.get(str(r["classifier"]).strip().lower())
            if cid and label:
                mapping[cid] = label
    return mapping


def detect_classifier(filename):
    """Identify the classifier that produced a per-cluster classification CSV."""
    base = os.path.basename(filename)
    if "_blastn_" in base:
        return "blastn"
    if "_seqmatch_" in base:
        return "seqmatch"
    if "_kraken2_" in base:
        return "kraken2"
    return None


def cluster_id_from_filename(filename):
    """Extract the cluster id from a per-cluster classification CSV filename.

    Files are named <prefix>_<cluster>_<tool>_consensus_classification.csv
    """
    m = re.search(
        r"_(\d+)_(?:blastn|seqmatch|kraken2)_consensus_classification",
        os.path.basename(filename),
    )
    return m.group(1) if m else "NA"


def _is_unclassified(name):
    return name is None or str(name).strip() == "" or str(name).strip().lower() == "unclassified"


def blast_records(path, top_n=TOP_N_HITS):
    """Structured top hits from a BLAST per-cluster CSV.

    Columns (no header, ';' separated): sscinames;staxids;evalue;length;pident;bitscore
    Rows are already sorted best-first by the pipeline. Confidence metric is % identity.
    Returns a list of dicts: rank, species, taxid, pct_identity, bitscore, selected.
    """
    df = pd.read_csv(path, sep=";", header=None, dtype=str).fillna("")
    records = []
    for rank, (_, r) in enumerate(df.head(top_n).iterrows(), start=1):
        species = r.get(0, "")
        if _is_unclassified(species):
            records.append({"rank": rank, "species": "unclassified", "taxid": "",
                            "pct_identity": "", "bitscore": "", "selected": False})
            continue
        records.append({"rank": rank, "species": species, "taxid": r.get(1, ""),
                        "pct_identity": r.get(4, ""), "bitscore": r.get(5, ""),
                        "selected": rank == 1})
    return records


def seqmatch_records(path, top_n=TOP_N_HITS):
    """Structured top hits from a SeqMatch per-cluster CSV.

    Columns (no header, ';' separated): sciname;taxid;seqmatch_score
    Rows are already sorted best-first by the pipeline. Confidence metric is the RDP S_ab score.
    Returns a list of dicts: rank, species, taxid, s_ab_score, selected.
    """
    df = pd.read_csv(path, sep=";", header=None, dtype=str).fillna("")
    records = []
    for rank, (_, r) in enumerate(df.head(top_n).iterrows(), start=1):
        species = r.get(0, "")
        if _is_unclassified(species):
            records.append({"rank": rank, "species": "unclassified", "taxid": "",
                            "s_ab_score": "", "selected": False})
            continue
        records.append({"rank": rank, "species": species, "taxid": r.get(1, ""),
                        "s_ab_score": r.get(2, ""), "selected": rank == 1})
    return records


def kraken2_records(path):
    """Structured single LCA call from a kraken2 report.

    kraken2 does lowest-common-ancestor assignment rather than producing a ranked
    candidate list, so there are no runners-up. Report columns (tab separated):
    pct, clade_reads, taxon_reads, rank_code, taxid, name. The directly-assigned node is
    the deepest row with a non-zero self (taxon) read count.
    Returns a single-element list of dicts: rank, species, taxid, lca_reads_pct, selected.
    """
    df = pd.read_csv(path, sep="\t", header=None, dtype=str).fillna("")
    names = ["pct", "clade_reads", "taxon_reads", "rank_code", "taxid", "name"]
    df.columns = names[: df.shape[1]]
    self_reads = pd.to_numeric(df.get("taxon_reads", 0), errors="coerce").fillna(0)
    assigned = df[self_reads > 0]
    if assigned.empty:
        return [{"rank": 1, "species": "unclassified", "taxid": "",
                 "lca_reads_pct": "", "selected": False}]
    row = assigned.iloc[-1]
    name = str(row["name"]).strip()
    return [{"rank": 1, "species": name, "taxid": str(row.get("taxid", "")).strip(),
             "lca_reads_pct": str(row.get("pct", "")).strip(),
             "selected": not _is_unclassified(name)}]


def _mark(record):
    """Species name with the selected-call asterisk appended where applicable."""
    if record["selected"] and not _is_unclassified(record["species"]):
        return record["species"] + SELECTED_MARK
    return record["species"]


def _blank(value):
    return value if str(value).strip() != "" else "-"


def parse_blast_hits(path, top_n=TOP_N_HITS):
    """Display table (Species / % Identity / Bitscore) for a BLAST per-cluster CSV."""
    rows = [{"Species": _mark(r), "% Identity": _blank(r["pct_identity"]),
             "Bitscore": _blank(r["bitscore"])} for r in blast_records(path, top_n)]
    return pd.DataFrame(rows, columns=["Species", "% Identity", "Bitscore"])


def parse_seqmatch_hits(path, top_n=TOP_N_HITS):
    """Display table (Species / S_ab score) for a SeqMatch per-cluster CSV."""
    rows = [{"Species": _mark(r), "S_ab score": _blank(r["s_ab_score"])}
            for r in seqmatch_records(path, top_n)]
    return pd.DataFrame(rows, columns=["Species", "S_ab score"])


def parse_kraken2_hit(path):
    """Display table (Species (LCA) / Reads %) for a kraken2 report."""
    rows = [{"Species (LCA)": _mark(r), "Reads %": _blank(r["lca_reads_pct"])}
            for r in kraken2_records(path)]
    return pd.DataFrame(rows, columns=["Species (LCA)", "Reads %"])


def build_hit_details(hit_details_dir, top_n=TOP_N_HITS):
    """Parse all per-cluster classification CSVs in a directory.

    Returns an ordered dict: {cluster_id: {classifier_label: DataFrame}}.
    """
    pattern = os.path.join(hit_details_dir, "*_consensus_classification.csv")
    files = sorted(glob.glob(pattern))
    clusters = {}
    for f in files:
        classifier = detect_classifier(f)
        if classifier is None:
            continue
        cid = cluster_id_from_filename(f)
        clusters.setdefault(cid, {})
        try:
            if classifier == "blastn":
                clusters[cid]["BLAST"] = parse_blast_hits(f, top_n)
            elif classifier == "seqmatch":
                clusters[cid]["SeqMatch"] = parse_seqmatch_hits(f, top_n)
            elif classifier == "kraken2":
                clusters[cid]["Kraken2 (LCA)"] = parse_kraken2_hit(f)
        except Exception:
            logger.exception("Failed to parse hit details from %s", f)
    return clusters


def add_hit_details_section(reprt, hit_details_dir, chosen_classifier='none', top_n=TOP_N_HITS):
    """Add a per-cluster 'top hits' section to the report if hit details are available.

    When a chosen-classifier map is supplied, each cluster shows only the classifier the
    pipeline selected for it (the same one that produced the reported abundance); otherwise
    every available classifier is shown.
    """
    if not hit_details_dir or hit_details_dir in ("none", "unknown"):
        return
    if not os.path.isdir(hit_details_dir):
        return
    clusters = build_hit_details(hit_details_dir, top_n)
    if not clusters:
        return
    chosen = load_chosen_classifiers(chosen_classifier)

    section = reprt.add_section()
    section.markdown('''
    <br/>
    ### Classifier hit details

    Top {0} hits for the classifier selected for each cluster, with the metric used to rank
    them (% identity for BLAST, S_ab score for SeqMatch). When two candidates are within a
    small margin the call is ambiguous and should be interpreted with care. kraken2 reports a
    single lowest-common-ancestor call rather than a ranked list. Entries marked with an
    asterisk (*) are the classifier's selected call.
    '''.format(top_n))

    # Sort clusters numerically where possible
    def _cluster_key(cid):
        return (0, int(cid)) if str(cid).isdigit() else (1, str(cid))

    for cid in sorted(clusters, key=_cluster_key):
        tables = clusters[cid]
        winner = chosen.get(str(cid))
        if winner is not None and winner in tables:
            tables = {winner: tables[winner]}
        section.markdown("**Cluster {0}**".format(cid))
        for classifier_label, table in tables.items():
            section.markdown("_{0}_".format(classifier_label))
            section.table(table, classes="larger-first-column")

def read_patient_info(file, barcode):
    """
    Read patient info file and return relevant row for the barcode
    
    Args:
        file (str): path to the patient info file
        barcode (str): barcode identifier for the patient sample
        
    Returns:
        list: list of pandas Series with patient info
    """
    #funtion parsing CSV patient file and looking up info for relevant barcode
    info=pd.read_excel(file, usecols=range(0,11))
    #return row with a barcode as a Series
    if barcode=="discontinued":
        relevant_row=[]
        relevant_rows=info.loc[info['Status'] == barcode]
        for index,row in relevant_rows.iterrows():
            single_row=row.transpose()
            single_df=single_row.reset_index()
            single_df.columns=['Metadata', 'Sample Information']
            relevant_row.append(single_df)
    else:
        relevant_rows=info.loc[info['Barcode'] == barcode]
        #move row names into a column
        relevant_row=relevant_rows.transpose()
        relevant_row.index.name = 'Metadata'
        relevant_row.reset_index(inplace=True)
        #rename columns
        relevant_row.columns=['Metadata', 'Sample Information']
    return relevant_row

def read_abundance_results(file):
    """
    Read abundance results csv file and return top 3 results

    Args:
        file (path): path to the abundance results file    

    Returns:
        pandas.DataFrame: top 3 abundance results
    """
    # The abundance results file is a CSV file
    abundance_results = pd.read_csv(file)
    
    # Rename the columns
    abundance_results.columns = ['Detected Species', 'Relative Abundance (%)', 'Number of Reads']
    
    # Round the abundance results
    abundance_results['Relative Abundance (%)'] = abundance_results['Relative Abundance (%)'].apply(np.around)
    
    # Return the top 3 abundance results
    abundance_results_t3 = abundance_results.head(n=3)
    
    return abundance_results_t3

def process_controls(ctrl):
    """
    Process controls files

    Args:
        control (path): Path to the positive or negative control file

    Returns:
        dataframe: Pandas dataframe with the control results
    """
    ctrl = ctrl[1:-1]
    if ctrl == '' or ctrl == 'None':
        control = None
    else:
        control = read_abundance_results(ctrl)
            

    return control

def write_hit_details_csv(hit_details_dir, outpath, chosen_classifier='none', top_n=TOP_N_HITS):
    """Write a consolidated hit-details CSV covering every cluster and classifier.

    Unlike the report (which shows only the winning classifier per cluster), this is the full
    record: top {top_n} hits for BLAST and SeqMatch and the single kraken2 LCA call, one row per
    hit. Returns the path written, or None if there were no per-cluster files.
    """
    if not hit_details_dir or hit_details_dir in ("none", "unknown"):
        return None
    if not os.path.isdir(hit_details_dir):
        return None
    files = sorted(glob.glob(os.path.join(hit_details_dir, "*_consensus_classification.csv")))
    if not files:
        return None

    # cluster -> classifier the pipeline selected (short name), for the cluster_winner column
    reverse_label = {label: token for token, label in CLASSIFIER_LABELS.items()}
    winners = {}
    for cid, label in load_chosen_classifiers(chosen_classifier).items():
        winners[cid] = CHOSEN_CSV_CLASSIFIER.get(reverse_label.get(label, ""), "")

    parsers = {"blastn": blast_records, "seqmatch": seqmatch_records, "kraken2": kraken2_records}
    rows = []
    for f in files:
        classifier = detect_classifier(f)
        if classifier is None:
            continue
        cid = cluster_id_from_filename(f)
        try:
            records = parsers[classifier](f, top_n) if classifier != "kraken2" else parsers[classifier](f)
        except Exception:
            logger.exception("Failed to parse hit details from %s", f)
            continue
        for r in records:
            rows.append({
                "cluster": cid,
                "classifier": HIT_CSV_CLASSIFIER[classifier],
                "cluster_winner": winners.get(cid, ""),
                "rank": r["rank"],
                "species": r["species"],
                "taxid": r.get("taxid", ""),
                "pct_identity": r.get("pct_identity", ""),
                "bitscore": r.get("bitscore", ""),
                "s_ab_score": r.get("s_ab_score", ""),
                "lca_reads_pct": r.get("lca_reads_pct", ""),
                "top_hit": r["selected"],
            })

    cols = ["cluster", "classifier", "cluster_winner", "rank", "species", "taxid",
            "pct_identity", "bitscore", "s_ab_score", "lca_reads_pct", "top_hit"]
    out = pd.DataFrame(rows, columns=cols)
    # Sort by numeric cluster, then classifier, then rank
    out["_c"] = pd.to_numeric(out["cluster"], errors="coerce")
    out = out.sort_values(["_c", "classifier", "rank"], kind="stable").drop(columns="_c")
    out.to_csv(outpath, index=False)
    return outpath


# ----------------------------------------------------------------------------
# Per-cluster QC traffic lights
# ----------------------------------------------------------------------------
# Tunable thresholds (see the report design notes)
QC_BLAST_SCORE_GREEN = 99.0     # >= green, [amber, green) amber, < amber red
QC_BLAST_SCORE_AMBER = 98.0
QC_SEQMATCH_SCORE_GREEN = 0.95
QC_SEQMATCH_SCORE_AMBER = 0.90
QC_BLAST_CLOSE_AMBER = 1.0      # top1 - top2 (% identity): <= amber, <= red red
QC_BLAST_CLOSE_RED = 0.3
QC_SEQMATCH_CLOSE_AMBER = 0.02  # top1 - top2 (S_ab)
QC_SEQMATCH_CLOSE_RED = 0.005
POS_CONTROL_SPECIES = "Marinobacter nauticus"   # spiked into every sample
SHOW_TOP_N_CLUSTERS = 3   # matches the previous report's top-3 abundance display; None = show all

_WIN_TOKEN = {"BLAST": "blast", "SeqMatch": "seqmatch", "Kraken2 (LCA)": "kraken2"}
_LEVEL_RANK = {"green": 0, "amber": 1, "red": 2}


def _qc_num(x):
    try:
        return float(str(x).strip())
    except (TypeError, ValueError):
        return None


def _species_key(name):
    """Normalise a species name to 'genus species' for comparison."""
    if _is_unclassified(name):
        return ""
    toks = re.sub(r"[^A-Za-z0-9 ]", " ", str(name)).lower().split()
    return " ".join(toks[:2])


def _worst(*levels):
    return max(levels, key=lambda lvl: _LEVEL_RANK[lvl])


def load_cluster_info(path):
    """Load cluster -> {classifier label, reads, rel_abundance} from chosen_classifier.csv."""
    info = {}
    if not path or path in ("none", "unknown") or not os.path.isfile(path):
        return info
    try:
        df = pd.read_csv(path, dtype=str).fillna("")
    except Exception:
        logger.exception("Failed to read chosen classifier map %s", path)
        return info
    for _, r in df.iterrows():
        cid = str(r.get("cluster", "")).strip()
        if not cid:
            continue
        info[cid] = {
            "classifier": CLASSIFIER_LABELS.get(str(r.get("classifier", "")).strip().lower(), ""),
            "reads": r.get("reads", ""),
            "rel_abundance": r.get("rel_abundance", ""),
        }
    return info


def collect_cluster_records(hit_details_dir, top_n=TOP_N_HITS):
    """{cluster_id: {'blast': [recs], 'seqmatch': [recs], 'kraken2': [recs]}}."""
    out = {}
    if not hit_details_dir or hit_details_dir in ("none", "unknown") or not os.path.isdir(hit_details_dir):
        return out
    parsers = {"blastn": ("blast", blast_records),
               "seqmatch": ("seqmatch", seqmatch_records),
               "kraken2": ("kraken2", kraken2_records)}
    for f in sorted(glob.glob(os.path.join(hit_details_dir, "*_consensus_classification.csv"))):
        classifier = detect_classifier(f)
        if classifier is None:
            continue
        token, parser = parsers[classifier]
        try:
            recs = parser(f) if classifier == "kraken2" else parser(f, top_n)
        except Exception:
            logger.exception("QC: failed to parse %s", f)
            continue
        out.setdefault(cluster_id_from_filename(f), {})[token] = recs
    return out


def _rank1(recs):
    return recs[0] if recs else None


def _close_light(recs, metric, amber, red):
    """Near-tie only counts when the runner-up is a DIFFERENT species.

    Strain-level ties (e.g. several 'Klebsiella pneumoniae, NBRC ...' rows) collapse to the
    same species key, so they are not flagged. We compare the top hit against the first
    runner-up that is a genuinely different species.
    """
    if not recs:
        return "green"
    top_key = _species_key(recs[0].get("species"))
    top_val = _qc_num(recs[0].get(metric))
    if not top_key or top_val is None:
        return "green"
    for r in recs[1:]:
        k = _species_key(r.get("species"))
        if k and k != top_key:
            v = _qc_num(r.get(metric))
            if v is None:
                return "green"
            gap = top_val - v
            if gap <= red:
                return "red"
            if gap <= amber:
                return "amber"
            return "green"
    return "green"   # all kept hits are the same species


def assess_cluster(recs, winner_token, neg_keys):
    """Return the four per-cluster QC light levels + the called species."""
    blast, seq, krak = recs.get("blast", []), recs.get("seqmatch", []), recs.get("kraken2", [])

    # 1. sequencer agreement: compare each classifier's top species
    keys = [_species_key(r["species"]) for r in (_rank1(blast), _rank1(seq), _rank1(krak)) if r]
    keys = [k for k in keys if k]
    if len(keys) <= 1 or len(set(keys)) == 1:
        agreement = "green"
    elif len(set(keys)) == len(keys):
        agreement = "red"           # no two agree
    else:
        agreement = "amber"         # a dissenter

    # 2. close top hits: worst near-tie across the ranked classifiers
    close = _worst(
        _close_light(blast, "pct_identity", QC_BLAST_CLOSE_AMBER, QC_BLAST_CLOSE_RED),
        _close_light(seq, "s_ab_score", QC_SEQMATCH_CLOSE_AMBER, QC_SEQMATCH_CLOSE_RED),
    )

    # 3. low absolute score: prefer BLAST % identity, else SeqMatch S_ab
    score = "amber"
    b1 = _rank1(blast)
    if b1 and not _is_unclassified(b1["species"]) and _qc_num(b1["pct_identity"]) is not None:
        p = _qc_num(b1["pct_identity"])
        score = "green" if p >= QC_BLAST_SCORE_GREEN else "amber" if p >= QC_BLAST_SCORE_AMBER else "red"
    else:
        s1 = _rank1(seq)
        if s1 and not _is_unclassified(s1["species"]) and _qc_num(s1["s_ab_score"]) is not None:
            v = _qc_num(s1["s_ab_score"])
            score = "green" if v >= QC_SEQMATCH_SCORE_GREEN else "amber" if v >= QC_SEQMATCH_SCORE_AMBER else "red"

    # the "call" = winning classifier's top species (fallback to any available)
    call_recs = recs.get(winner_token) if winner_token else None
    if not call_recs:
        for t in ("blast", "seqmatch", "kraken2"):
            if recs.get(t):
                call_recs = recs[t]
                break
    call = _rank1(call_recs)
    call_species = call["species"] if call else "unclassified"

    # 4. matches negative control (simple: called species seen in neg control)
    ck = _species_key(call_species)
    negctrl = "red" if ck and ck in neg_keys else "green"

    return {"agreement": agreement, "close": close, "score": score,
            "negctrl": negctrl, "call_species": call_species}


def marinobacter_present(clusters, name=POS_CONTROL_SPECIES):
    """True if the positive-control spike appears in any cluster / classifier."""
    key = _species_key(name)
    for recs in clusters.values():
        for token in ("blast", "seqmatch", "kraken2"):
            for r in recs.get(token, []):
                if _species_key(r["species"]) == key:
                    return True
    return False


QC_CSS = """
<style>
.qc{--g:#16a34a;--a:#d97706;--r:#dc2626}
.qc .lamp{display:inline-block;width:12px;height:12px;border-radius:50%;vertical-align:middle}
.qc .g .lamp{background:var(--g)} .qc .a .lamp{background:var(--a)} .qc .r .lamp{background:var(--r)}
.qc a.lamp-link{text-decoration:none} .qc a.lamp-link:hover .lamp{transform:scale(1.3)}
.qc .qc-legend{font-size:12px;color:#5c6773;margin:0 0 12px}
.qc .qc-legend i{width:9px;height:9px;border-radius:50%;display:inline-block;vertical-align:middle;margin:0 4px 0 14px}
.qc table{border-collapse:collapse;margin:0 0 1rem;font-size:14px}
.qc>table{width:100%}
.qc th,.qc td{padding:.5rem .75rem;border:1px solid #dee2e6;text-align:left;vertical-align:middle}
.qc thead th{background-color:var(--brand-primary,#0084a9);color:#fff;font-weight:700;border-color:var(--brand-primary,#0084a9)}
.qc thead th a{color:#fff;text-decoration:underline}
.qc th.c,.qc td.c{text-align:center}
.qc-explain h4{margin:14px 0 3px;font-size:14px} .qc-explain p{margin:0 0 4px}
.qc-explain{font-size:13.5px}
.qc .sci{font-style:italic}
.qc .qc-detail{display:none}
.qc .qc-detail:target{display:table-row}
.qc .qc-detail>td{background:#f6f8fa}
.qc .qc-detail table{margin:8px 0 2px}
.qc .reason{margin:2px 0 8px} .qc .reason .lamp{margin-right:8px}
</style>
"""


def _lamp(level, href=None):
    cls = {"green": "g", "amber": "a", "red": "r"}[level]
    if href and level != "green":
        return '<a class="lamp-link {0}" href="#{1}"><span class="lamp"></span></a>'.format(cls, href)
    return '<span class="{0}"><span class="lamp"></span></span>'.format(cls)


def _esc(s):
    return (str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))


def _cluster_sort_key(cid):
    return (0, int(cid)) if str(cid).isdigit() else (1, str(cid))


def _cluster_reads(cid, cluster_info):
    try:
        return float(cluster_info.get(str(cid), {}).get("reads") or 0)
    except (TypeError, ValueError):
        return 0.0


def build_qc_html(clusters, cluster_info, neg_species, top_clusters=SHOW_TOP_N_CLUSTERS):
    """Build the QC traffic-light section HTML from parsed cluster records."""
    neg_keys = {_species_key(s) for s in (neg_species or []) if _species_key(s)}

    # show the most abundant clusters first, limited to the previous top-N display
    ordered = sorted(clusters, key=lambda c: -_cluster_reads(c, cluster_info))
    shown = ordered[:top_clusters] if top_clusters else ordered

    body = []
    for cid in shown:
        recs = clusters[cid]
        info = cluster_info.get(str(cid), {})
        winner_token = _WIN_TOKEN.get(info.get("classifier", ""))
        a = assess_cluster(recs, winner_token, neg_keys)
        lights = {k: a[k] for k in ("agreement", "close", "score", "negctrl")}
        flagged = any(v != "green" for v in lights.values())
        anchor = "qc-{0}".format(cid)
        reads, pct = info.get("reads", ""), info.get("rel_abundance", "")
        reads_txt = ("{0} ({1}%)".format(reads, pct)
                     if str(reads).strip() not in ("", "nan") else "-")

        def cell(level):
            return '<td class="c">{0}</td>'.format(_lamp(level, anchor if flagged else None))

        body.append(
            '<tr><td>{cid}</td><td><span class="sci">{sp}</span></td>'
            '<td class="c">{reads}</td>{ag}{cl}{sc}{ng}</tr>'.format(
                cid=_esc(cid), sp=_esc(a["call_species"]), reads=reads_txt,
                ag=cell(lights["agreement"]), cl=cell(lights["close"]),
                sc=cell(lights["score"]), ng=cell(lights["negctrl"])))

        if flagged:
            body.append(_build_detail(cid, anchor, recs, a, lights))   # dropdown right below its row

    if not body:
        return ""

    legend = ('<p class="qc-legend">Per-cluster checks &mdash; '
              '<i style="background:#16a34a"></i>confident '
              '<i style="background:#d97706"></i>interpret with care '
              '<i style="background:#dc2626"></i>unreliable / QC concern. '
              'Click an amber/red light for detail; click a column heading for how it is scored.</p>')

    table = (
        '<table>'
        '<thead><tr><th>Cluster</th><th>Call</th><th class="c">Reads (%)</th>'
        '<th class="c"><a href="#qcx-agreement">Agreement</a></th>'
        '<th class="c"><a href="#qcx-close">Close hits</a></th>'
        '<th class="c"><a href="#qcx-score">Abs. score</a></th>'
        '<th class="c"><a href="#qcx-negctrl">Neg. control</a></th></tr></thead>'
        '<tbody>{0}</tbody></table>'.format("".join(body)))

    return QC_CSS + '<div class="qc">' + legend + table + '</div>'


def _map3(level):
    return {"green": "g", "amber": "a", "red": "r"}[level]


def _build_detail(cid, anchor, recs, a, lights):
    reasons = []
    if lights["agreement"] != "green":
        calls = []
        for token, label in (("blast", "BLAST"), ("kraken2", "kraken2"), ("seqmatch", "SeqMatch")):
            r = _rank1(recs.get(token, []))
            if r and not _is_unclassified(r["species"]):
                calls.append("{0}: <i>{1}</i>".format(label, _esc(r["species"])))
        reasons.append((lights["agreement"],
                        "<b>Classifiers disagree.</b> " + "; ".join(calls)))
    if lights["score"] != "green":
        b1 = _rank1(recs.get("blast", []))
        detail = ""
        if b1 and _qc_num(b1["pct_identity"]) is not None:
            detail = " Best BLAST hit {0}% identity.".format(_esc(b1["pct_identity"]))
        reasons.append((lights["score"], "<b>Low absolute score.</b>" + detail))
    if lights["close"] != "green":
        reasons.append((lights["close"],
                        "<b>Close top hits.</b> Runner-up is within the near-tie margin."))
    if lights["negctrl"] != "green":
        reasons.append((lights["negctrl"],
                        "<b>Also in the negative control.</b> <i>{0}</i> appears in the neg-control "
                        "results &mdash; possible contamination.".format(_esc(a["call_species"]))))

    reason_html = "".join(
        '<p class="reason {0}"><span class="lamp"></span>{1}</p>'.format(_map3(lvl), txt)
        for lvl, txt in reasons)

    # a small per-classifier breakdown table (unclassified hits are omitted; kraken2 has no
    # comparable score so it shows "-")
    mini_rows = []
    for token, label in (("blast", "BLAST"), ("seqmatch", "SeqMatch"), ("kraken2", "kraken2")):
        real = [r for r in recs.get(token, []) if not _is_unclassified(r["species"])]
        for i, r in enumerate(real):
            if str(r.get("pct_identity", "")).strip():
                score = _esc(r["pct_identity"]) + "%"
            elif str(r.get("s_ab_score", "")).strip():
                score = _esc(r["s_ab_score"])
            else:
                score = "-"
            mini_rows.append(
                '<tr><td>{lab}</td><td><span class="sci">{sp}</span></td>'
                '<td class="c">{sc}</td></tr>'.format(
                    lab=label if i == 0 else "", sp=_esc(r["species"]), sc=score))
    mini = ('<table><thead><tr><th>Classifier</th><th>Hit</th>'
            '<th class="c">Score</th></tr></thead><tbody>{0}</tbody></table>'.format("".join(mini_rows))
            if mini_rows else "")

    return ('<tr class="qc-detail" id="{a}"><td colspan="7">'
            '<b>Cluster {c}</b>{reasons}{mini}</td></tr>'.format(
                a=anchor, c=_esc(cid), reasons=reason_html, mini=mini))


def build_qc_explanations():
    """Reference block explaining each QC check; the table column headings link here."""
    return (
        '<div class="qc-explain">'
        '<h4 id="qcx-agreement">Agreement</h4>'
        '<p>Do the three classifiers name the same species for this cluster? '
        '<b>Green</b> = all agree; <b>amber</b> = one dissents; <b>red</b> = all differ. '
        'Compared at species level, so different strains of one species (for example several '
        '<i>Klebsiella pneumoniae</i> reference strains) count as agreement.</p>'
        '<h4 id="qcx-close">Close hits</h4>'
        '<p>Is the top hit clearly ahead of the next <i>different</i> species? '
        '<b>Amber/red</b> when a different species sits within the near-tie margin, so the call '
        'could plausibly be either. Margins: BLAST within about 1% identity (red within 0.3%); '
        'SeqMatch within about 0.02 S_ab (red within 0.005).</p>'
        '<h4 id="qcx-score">Absolute score</h4>'
        '<p>How strong is the best hit on its own? BLAST % identity: '
        '<b>green</b> at least 99%, <b>amber</b> 98&ndash;99%, <b>red</b> below 98%. '
        'SeqMatch S_ab: <b>green</b> at least 0.95, <b>amber</b> 0.90&ndash;0.95, <b>red</b> below 0.90.</p>'
        '<h4 id="qcx-negctrl">Negative control</h4>'
        '<p>Is the called species also detected in the negative control for this run? '
        '<b>Red</b> if it is &mdash; likely reagent or environmental contamination rather than a true finding.</p>'
        '<h4>How each classifier scores</h4>'
        '<p><b>BLAST</b> &mdash; % identity of the best alignment to the 16S database '
        '(higher is better; a confident species match is usually at least 99%). '
        '<b>SeqMatch (RDP)</b> &mdash; S_ab similarity score from 0 to 1 (higher is better). '
        '<b>kraken2</b> &mdash; a lowest-common-ancestor assignment built from shared k-mers; '
        'it returns a single taxon with no numeric score, so it often stops at genus and shows '
        '&ldquo;-&rdquo; in the score column.</p>'
        '</div>')


def internal_control_html(hit_details_dir, name=POS_CONTROL_SPECIES):
    """Internal spike-control status line for the Run QC section."""
    clusters = collect_cluster_records(hit_details_dir)
    if not clusters:
        return ""
    present = marinobacter_present(clusters, name)
    color = "#16a34a" if present else "#dc2626"
    dot = ('<span style="display:inline-block;width:12px;height:12px;border-radius:50%;'
           'background:{0};vertical-align:middle;margin-right:8px"></span>'.format(color))
    if present:
        msg = "Internal spike control <i>{0}</i> detected.".format(_esc(name))
    else:
        msg = ('Internal spike control <i>{0}</i> <font color="red">NOT detected &mdash; the run '
               'may have failed and results should be treated with caution</font>.'.format(_esc(name)))
    return '<br/>\n<b>INTERNAL CONTROL</b><br/>\n{0}{1}'.format(dot, msg)


def add_qc_section(reprt, hit_details_dir, chosen_classifier="none", neg_species=None, top_n=TOP_N_HITS):
    """Add the per-cluster QC traffic-light section to the report."""
    clusters = collect_cluster_records(hit_details_dir, top_n)
    if not clusters:
        return
    cluster_info = load_cluster_info(chosen_classifier)
    html = build_qc_html(clusters, cluster_info, neg_species)
    if not html:
        return
    section = reprt.add_section()
    section.markdown("<br/>\n### Sequence identification & QC\n")
    section.markdown(html)
    section.markdown(build_qc_explanations())


def parse_args():
    """Run the entry point."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--infile", default='unknown',
        help='Table file with classification and abundance results')
    parser.add_argument(
        "--output", default='unknown',
        help="Report output file name")
    parser.add_argument(
        "--barcode", default='discontinued',
        help="barcode identifier for the patient sample")
    parser.add_argument(
        "--info", required=True,
        help="Experiment information file mapping patient ID with a barcode")
    parser.add_argument(
        "--revision", default='unknown',
        help="git branch/tag of the executed workflow")
    parser.add_argument(
        "--commit", default='unknown',
        help="git commit of the executed workflow")
    parser.add_argument(
        "--demux", default='unknown',
        help='demultiplexing method and software version')
    parser.add_argument(
        "--clustering_size", default='unknown',
        help="Amount of reads used for UMAP HDBSCAN clustering")
    parser.add_argument(
        "--positive", default='unknown',
        help="File path for positive control results")
    parser.add_argument(
        "--negative", default='unknown',
        help="File path for negative control results")
    parser.add_argument(
        "--reads_count", default='0',
        help="Reads count after quality control")
    parser.add_argument(
        "--kit", default='unknown',
        help="Kit used for barcoding and demultiplexing")
    parser.add_argument(
        "--report_template",
        help="path to report template")
    parser.add_argument(
        "--logo",
        help="custom logo")
    parser.add_argument(
        "--run_id", default='unknown',
        help="Run ID of the sequencing run")
    parser.add_argument(
        "--seq_start", default='unknown',
        help="Start time of the sequencing run")
    parser.add_argument(
        "--hit_details", default='none',
        help="Directory containing per-cluster *_consensus_classification.csv files "
             "used to render the top hits per cluster")
    parser.add_argument(
        "--chosen_classifier", default='none',
        help="CSV mapping cluster -> winning classifier (from get_abundance.py). When "
             "given, each cluster shows only the classifier the pipeline selected for it")

    args = parser.parse_args()

    return(args)

def main(args):
    
    if args.barcode=="discontinued":

        metadata_table_list=read_patient_info(args.info, args.barcode)
                
        for patient in metadata_table_list:
            # Restructure the metadata table
            restructured=[]
            for index, row in patient.iloc[[0,2,1,4,6,7,10]].iterrows():
                restructured.append(": ".join([str(row['Metadata']), str(row['Sample Information'])]))
            restructured.insert(4, " ".join(["Sequencing start:", args.seq_start]))
            rest_df=pd.DataFrame(list(zip(restructured[:4],restructured[4:])), columns=['Sample Information', 'Time Stamps'])

            # Generate the report
            title="Patient " + patient.iloc[0,1] + " Report"
            reprt = report.UoSReport(
                title=title, report_template=args.report_template, about=False, style='UoS', logo=args.logo)

            section=reprt.add_section()
            section.markdown('''
            ### Sample Information
            ''')

            section.table(rest_df)

            section=reprt.add_section()

            assay_type=patient.iloc[4]['Sample Information']
            if assay_type == '16S':
                assay_info = 'Bacterial 16s'
            else:
                assay_info = 'Fungal ITS2'

            section.markdown('''
            <br/>
            ### Results

            <font color="red">**{0} rRNA NOT detected**</font>
            '''.format(assay_info))

            reprt.write("patient_report_" + str(patient.iloc[1,1]) + ".html")
    else:        
        metadata_table=read_patient_info(args.info, args.barcode)
        restructured=[]
        for index, row in metadata_table.iloc[[0,2,1,4,6,7,10]].iterrows():
            restructured.append(": ".join([str(row['Metadata']), str(row['Sample Information'])]))
        restructured.insert(4, " ".join(["Sequencing start:", args.seq_start]))
        rest_df=pd.DataFrame(list(zip(restructured[:4],restructured[4:])), columns=['Sample Information', 'Time Stamps'])

        # Create the title for the report
        title="Patient " + metadata_table.iloc[0,1] + " Report"

        if args.infile == "input.1":
            results_table = None
        else:
            results_table=read_abundance_results(args.infile)

        positive=process_controls(args.positive)
        negative=process_controls(args.negative)

        reprt = report.UoSReport(
            title=title, workflow="NanoCLUST", report_template=args.report_template,
            revision=args.revision, commit=args.commit, style='UoS', logo=args.logo)

        section=reprt.add_section()
        section.markdown('''
        ### Sample Information
        ''')

        section.table(rest_df)

        section=reprt.add_section()

        section.markdown('''
        <br/>
        ### Results

        Total reads in this sample: {0}
        '''.format(args.reads_count))

        assay_type=metadata_table.loc[metadata_table['Metadata'] == 'Assay', 'Sample Information'].iloc[0]

        if assay_type == '16S':
            database_info="16s bacterial sequencing results were compared against 16S & 18S database, build 18 Jan 2022."
            infection_type='bacterial'
            assay_info='Bacterial 16s'
        else:
            database_info="ITS2 fungal sequencing results were compared against ITS2 database, build 15 Mar 2022."
            infection_type='fungal'
            assay_info='Fungal ITS2'

        # The per-cluster QC table below is now the results table, so we no longer show the
        # separate top-3 abundance table (they duplicated each other). Only the "not detected"
        # message remains for the empty case.
        if results_table is None:
            section.markdown('''
            <font color="red">**{0} rRNA NOT detected**</font>
            '''.format(assay_info))

        # Per-cluster QC traffic lights (agreement, close hits, absolute score, neg-control
        # match) plus the Marinobacter nauticus positive-control spike check
        neg_species = list(negative['Detected Species']) if negative is not None else []
        add_qc_section(reprt, args.hit_details, args.chosen_classifier, neg_species)

        # Full machine-readable record: top hits for every classifier, all clusters
        write_hit_details_csv(
            args.hit_details,
            "hit_details_" + str(metadata_table.iloc[1, 1]) + ".csv",
            args.chosen_classifier)

        section=reprt.add_section()
        section.markdown('''
        <br/>
        ### Run QC


        **NEGATIVE CONTROL**
        ''')
        comment=0
        if negative is not None:
            comment+=10
            section.markdown('''
            Total reads in negative control: {0} 
            '''.format(negative['Number of Reads'].sum()))

            section.table(negative, classes='larger-first-column')
        else:
            section.markdown('''
            No species detected in negative control.
            ''')

        section.markdown('''
        <br/>
        **POSITIVE CONTROL**
        ''')

        if positive is not None:
            comment+=1
            section.markdown('''
            Total reads in positive control: {0}
            '''.format(positive['Number of Reads'].sum()))
            
            section.table(positive, classes='larger-first-column')
        else:
            comment+=2
            section.markdown('''
            No species detected in positive control.
            ''')

        # Internal spike control (Marinobacter nauticus) added to every sample
        internal = internal_control_html(args.hit_details)
        if internal:
            section.markdown(internal)

        if comment == 1:
            section.markdown('''
            comment: <font color="green">QC for this sample was **successful**</font>
            ''')
        elif comment == 11:
            section.markdown('''
            comment: <font color="orange">QC for this sample shows the presence of {0} reads in the negative control. Check if there is overlap with any detected pathogen in the sample that may indicate contamination.</font>
            '''.format(infection_type))
        else:
            section.markdown('''
            comment: <font color="red">QC for this sample has **failed** and results cannot be validated</font>
            ''')

        section=reprt.add_section()
        run_id=args.run_id
        barcoding_kit=args.kit
        print(barcoding_kit)
        demux_method=args.demux
        species_database=metadata_table['Sample Information'].iloc[4]
        clustering_size=args.clustering_size

        run_params=pd.DataFrame(list(zip(['Run ID: '+str(run_id), 'Barcoding kit: '+str(barcoding_kit), 'Demultiplex method: '+str(demux_method)], ['Species database: '+str(species_database), 'Clustering size: '+str(clustering_size), 'Sample barcode: '+str(metadata_table.iloc[1,1])])), columns=['GridIon properties', 'NanoCLUST properties'])
        #run_params.reset_index(drop=True, inplace=True)

        section.markdown('''
        ### Run parameters
        ''')
        section.table(run_params)
        section.markdown('''
        Sample was sequenced on a ONT GridION Mk1. 
        Sequencing data was processed and analysed using a custom nanoclust pipeline.
        {6}

        '''.format(run_id, barcoding_kit, demux_method, species_database, clustering_size, metadata_table.iloc[1,1], database_info))

        #write report
        reprt.write(args.output + "_" + str(metadata_table.iloc[1,1]) + ".html")


if __name__ == "__main__":
    args = parse_args()

    main(args)