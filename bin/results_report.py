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
    if not recs or len(recs) < 2:
        return "green"
    a, b = _qc_num(recs[0].get(metric)), _qc_num(recs[1].get(metric))
    if a is None or b is None:
        return "green"
    gap = a - b
    if gap <= red:
        return "red"
    if gap <= amber:
        return "amber"
    return "green"


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
.qc{--g:#16a34a;--a:#d97706;--r:#dc2626;font-size:14px}
.qc .lamp{display:inline-block;width:13px;height:13px;border-radius:50%;vertical-align:middle}
.qc .g .lamp{background:var(--g)} .qc .a .lamp{background:var(--a)} .qc .r .lamp{background:var(--r)}
.qc a.lamp-link{text-decoration:none} .qc a.lamp-link:hover .lamp{transform:scale(1.3)}
.qc .banner{display:flex;flex-wrap:wrap;gap:12px;margin:6px 0}
.qc .card{flex:1 1 220px;border:1px solid #d7dde3;border-radius:10px;padding:11px 14px;display:flex;gap:12px;align-items:center}
.qc .card .lamp{width:15px;height:15px;flex:none}
.qc .card b{font-size:13px} .qc .card small{display:block;color:#5c6773;font-size:12px}
.qc .legend{font-size:12px;color:#5c6773;margin:6px 0 10px}
.qc .legend i{width:9px;height:9px;border-radius:50%;display:inline-block;vertical-align:middle;margin:0 4px 0 12px}
.qc table{border-collapse:collapse;width:100%;font-size:13.5px}
.qc thead th{text-align:left;font-size:11px;letter-spacing:.04em;text-transform:uppercase;color:#5c6773;padding:8px 10px;border-bottom:1px solid #d7dde3}
.qc thead th.c{text-align:center}
.qc tbody td{padding:9px 10px;border-bottom:1px solid #e6eaee}
.qc tbody td.c{text-align:center}
.qc .sci{font-style:italic}
.qc .cid{font-weight:700;color:#0f766e}
.qc .reads{color:#5c6773;font-variant-numeric:tabular-nums}
.qc .qc-detail{display:none;background:#f6f8fa}
.qc .qc-detail:target{display:table-row}
.qc .qc-detail td{padding:12px 16px}
.qc .reason{margin:0 0 8px} .qc .reason .lamp{margin-right:8px}
.qc .mini{border-collapse:collapse;margin-top:6px;font-size:13px}
.qc .mini th,.qc .mini td{border-bottom:1px solid #e6eaee;padding:5px 10px;text-align:left}
.qc .mini td.n{text-align:right;font-variant-numeric:tabular-nums}
.qc .pill{font-size:11px;font-weight:600;color:#16a34a;margin-left:6px}
</style>
"""


def _lamp(level, href=None):
    cls = {"green": "g", "amber": "a", "red": "r"}[level]
    dot = '<span class="{0}"><span class="lamp"></span></span>'.format(cls)
    if href and level != "green":
        return '<a class="lamp-link {0}" href="#{1}"><span class="lamp"></span></a>'.format(cls, href)
    return dot


def _esc(s):
    return (str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))


def _cluster_sort_key(cid):
    return (0, int(cid)) if str(cid).isdigit() else (1, str(cid))


def build_qc_html(clusters, cluster_info, neg_species):
    """Build the QC traffic-light section HTML from parsed cluster records."""
    neg_keys = {_species_key(s) for s in (neg_species or []) if _species_key(s)}

    rows_html, detail_html = [], []
    flagged_total = 0
    for cid in sorted(clusters, key=_cluster_sort_key):
        recs = clusters[cid]
        info = cluster_info.get(str(cid), {})
        winner_label = info.get("classifier", "")
        winner_token = _WIN_TOKEN.get(winner_label)
        a = assess_cluster(recs, winner_token, neg_keys)
        lights = {k: a[k] for k in ("agreement", "close", "score", "negctrl")}
        flagged = any(v != "green" for v in lights.values())
        anchor = "qc-{0}".format(cid)
        reads = info.get("rel_abundance", "")
        reads_txt = "{0}%".format(reads) if str(reads).strip() not in ("", "nan") else "-"

        def cell(level):
            return '<td class="c">{0}</td>'.format(_lamp(level, anchor if flagged else None))

        rows_html.append(
            '<tr><td><span class="cid">{cid}</span></td>'
            '<td><span class="sci">{sp}</span></td>'
            '<td class="c reads">{reads}</td>{ag}{cl}{sc}{ng}</tr>'.format(
                cid=_esc(cid), sp=_esc(a["call_species"]), reads=reads_txt,
                ag=cell(lights["agreement"]), cl=cell(lights["close"]),
                sc=cell(lights["score"]), ng=cell(lights["negctrl"])))

        if flagged:
            flagged_total += 1
            detail_html.append(_build_detail(cid, anchor, recs, a, lights))

    if not rows_html:
        return ""

    n = len(rows_html)
    overall = "green" if flagged_total == 0 else "amber" if flagged_total < n else "red"
    overall_txt = ("all clusters confident" if flagged_total == 0
                   else "{0} of {1} cluster(s) flagged".format(flagged_total, n))
    pos = marinobacter_present(clusters)
    pos_level = "green" if pos else "red"
    pos_txt = ("<i>{0}</i> present".format(_esc(POS_CONTROL_SPECIES)) if pos
               else "<i>{0}</i> NOT detected".format(_esc(POS_CONTROL_SPECIES)))

    banner = (
        '<div class="banner">'
        '<div class="card {ol}"><span class="lamp"></span><div><b>Overall</b>'
        '<small>{ot}</small></div></div>'
        '<div class="card {pl}"><span class="lamp"></span><div><b>Positive-control spike</b>'
        '<small>{pt}</small></div></div></div>'.format(ol=_map3(overall), ot=overall_txt,
                                                       pl=_map3(pos_level), pt=pos_txt))

    legend = ('<p class="legend">Per-cluster checks &mdash; '
              '<i style="background:#16a34a"></i>confident '
              '<i style="background:#d97706"></i>interpret with care '
              '<i style="background:#dc2626"></i>unreliable / QC concern. '
              'Click an amber/red light for detail.</p>')

    table = (
        '<table><thead><tr><th>Cluster</th><th>Call</th><th class="c">Reads</th>'
        '<th class="c">Agreement</th><th class="c">Close hits</th>'
        '<th class="c">Abs. score</th><th class="c">Neg. control</th></tr></thead>'
        '<tbody>{0}</tbody></table>'.format("".join(rows_html))
        + "".join(detail_html))

    return QC_CSS + '<div class="qc">' + banner + legend + table + '</div>'


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

    # a small per-classifier breakdown table
    mini_rows = []
    for token, label in (("blast", "BLAST"), ("seqmatch", "SeqMatch"), ("kraken2", "kraken2")):
        for i, r in enumerate(recs.get(token, [])):
            metric = (r.get("pct_identity") or r.get("s_ab_score") or r.get("lca_reads_pct") or "")
            unit = ("%" if r.get("pct_identity") else "")
            mini_rows.append(
                '<tr><td>{lab}</td><td><span class="sci">{sp}</span></td>'
                '<td class="n">{m}{u}</td></tr>'.format(
                    lab=label if i == 0 else "", sp=_esc(r["species"]),
                    m=_esc(metric) if str(metric).strip() else "-", u=unit))
    mini = ('<table class="mini"><thead><tr><th>Classifier</th><th>Hit</th>'
            '<th class="n">Score</th></tr></thead><tbody>{0}</tbody></table>'.format("".join(mini_rows)))

    return ('<tr class="qc-detail" id="{a}"><td colspan="7">'
            '<b>Cluster {c}</b>{reasons}{mini}</td></tr>'.format(
                a=anchor, c=_esc(cid), reasons=reason_html, mini=mini))


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

        if results_table is not None:
            section.table(results_table, classes=['highlighted', 'larger-first-column'])
        else:
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