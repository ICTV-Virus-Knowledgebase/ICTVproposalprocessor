#!/usr/bin/env python3

"""Native Python proposal validator converted from merge_proposal_zips.R.

The first version preserves the R command-line and output contracts while
intentionally omitting reference and proposal caches. Correctness and explicit
pandas mutation are preferred over premature optimization.
"""

from __future__ import annotations

import argparse
import csv
import re
import subprocess
import sys
import tempfile
import warnings
import zipfile
from pathlib import Path
from typing import Any, Iterable, Iterator, Mapping, Sequence
from xml.etree import ElementTree

from openpyxl import load_workbook

try:
    import pandas as pd
except ImportError as exc:
    raise SystemExit(
        "ERROR: pandas is required. Create the environment with "
        "`conda env create -f environment.python.yml`."
    ) from exc


# ---------------------------------------------------------------------------
# R-compatible process state
# ---------------------------------------------------------------------------
#
# These module-level objects intentionally retain the R script names. Proposal
# processing mutates them in stages: reference load, proposal QC, change
# ordering/application, and output generation.

actionOrder = 0
cvList: dict[str, Any] = {}
dbCvList: dict[str, Any] = {}
dbCvMapList: dict[str, Any] = {}
scAbbrevNameMap: dict[str, str] = {}
oldMSLs = pd.DataFrame()
curMSL = pd.DataFrame()
newMSL = pd.DataFrame()
vmrDf = pd.DataFrame()
vmrAccessions = pd.DataFrame()
abolishedAccessions = pd.DataFrame(columns=["accession", "species", "virus_name"])
proposalsDf = pd.DataFrame()
changeList: dict[str, pd.DataFrame] = {}
allChangeDf = pd.DataFrame()

ERROR_COLUMNS = [
    "subcommittee", "code", "docx", "xlsx", "row", "change", "rank",
    "taxon", "level", "error", "message", "notes", "scAbbrev",
    "validator_version", "order",
]
SUMMARY_METADATA_COLUMNS = ["docpath", "path", "file", "basename"]
allErrorDf = pd.DataFrame(columns=ERROR_COLUMNS)


# ---------------------------------------------------------------------------
# Command-line compatibility
# ---------------------------------------------------------------------------
#
# Params preserves both attribute access and dotted R option names. parse_args
# keeps the established flags so existing shell and container entrypoints can
# invoke the Python implementation with minimal changes.


class Params(dict):
    """Dictionary with attribute access and support for exact dotted R keys."""

    def __getattr__(self, name: str) -> Any:
        try:
            return self[name]
        except KeyError as exc:
            raise AttributeError(name) from exc

    def __setattr__(self, name: str, value: Any) -> None:
        self[name] = value


params = Params()


def parse_args(argv: Sequence[str] | None = None) -> Params:
    parser = argparse.ArgumentParser(
        description="Load and QC ICTV taxonomic proposals", allow_abbrev=False
    )
    parser.add_argument("-v", "--verbose", action="store_true", default=False,
                        help="Print extra output")
    parser.add_argument("-t", "--tmi", action="store_true", default=False,
                        help="Print lots of extra output (Too Much Information)")
    parser.add_argument("-q", "--quiet", action="store_false", dest="verbose",
                        help="Print no output")
    parser.add_argument("-d", "--debug", action="store_true", default=False,
                        dest="debug_on_error")
    parser.add_argument("--noInfo", action="store_false", default=True,
                        dest="show.xlsx.code_miss")
    # Accepted for bash compatibility; caching is intentionally inactive.
    parser.add_argument("-c", "--useCache", action="store_true", default=False,
                        dest="use_cache")
    parser.add_argument("-u", "--updateCache", action="store_true", default=False,
                        dest="update_cache")
    parser.add_argument("-C", "--loadProposalCache", action="store_true",
                        default=False, dest="load_proposal_cache")
    parser.add_argument("-U", "--saveProposalCache", action="store_true",
                        default=False, dest="save_proposal_cache")
    parser.add_argument("--msl", action="store_true", default=False,
                        dest="export_msl")
    parser.add_argument("-D", "--outputDeltas", action="store_true",
                        default=False, dest="output_change_report")
    parser.add_argument("-i", "--proposalsDir", default="proposalsTest",
                        dest="proposals_dir")
    parser.add_argument("-o", "--outDir", default="results", dest="out_dir")
    parser.add_argument("-r", "--refDir", default="current_msl/msl41v1",
                        dest="ref_dir")
    parser.add_argument("--mslTsv", default="msl.tsv", dest="msl_tsv")
    parser.add_argument("--docxSummary", default="QC.docx_summary.tsv",
                        dest="proposals_meta")
    parser.add_argument("--mslLoadSql", default="msl_load.sql",
                        dest="sql_load_filename")
    parser.add_argument("--sqlInsertBatch", default="200",
                        dest="sql_insert_batch_size")
    parser.add_argument("--taxnodeIdDelta", default="100000",
                        dest="taxnode_delta")
    parser.add_argument("--qcPrettySummary", default="QC.pretty_summary.all.xlsx",
                        dest="qc_summary_fname")
    parser.add_argument("--qcTsvSummary", default="QC.summary.tsv",
                        dest="qc_summary_tsv_fname")
    parser.add_argument("--qcTsvRegression", default="QC.regression.new.tsv",
                        dest="qc_regression_tsv_fname")
    parser.add_argument("--dbSourceHost", default="taxonomy_host_source.utf8.txt",
                        dest="db_host_source_fname")
    parser.add_argument("--dbRanks", default="taxonomy_level.utf8.txt",
                        dest="db_rank_fname")
    parser.add_argument("--dbMolecules", default="taxonomy_molecule.utf8.txt",
                        dest="db_molecule_fname")
    parser.add_argument("--dbTaxa", default="taxonomy_node_export.utf8.txt",
                        dest="db_taxonomy_node_fname")
    parser.add_argument("--cvTemplate", default="TP_Template.xlsx",
                        dest="template_xlsx_fname")
    parser.add_argument("--cvTemplateSheet", default="Menu Items (Do not change)",
                        dest="template_xlsx_sheet")
    parser.add_argument("--vmr", default="species_isolates.utf8.txt",
                        dest="vmr_fname")
    parser.add_argument("--abolishedAccessions", default="abolished_accessions.tsv",
                        dest="abolished_accessions_fname")
    parser.add_argument("--templateURL",
                        default="https://ictv.global/taxonomy/templates",
                        dest="template_url")
    parser.add_argument("--cacheFile", default=".RData", dest="cache_fname")
    parser.add_argument(
        "--infileSupplPat",
        default=r"(^|[_. -])(suppl|viridic|table|matrix|sup)([_. -]|$)",
        dest="infile_suppl_pat",
    )
    parser.add_argument("--mode", default="validate", dest="processing_mode")
    parser.add_argument("--version_file", default="version_git.txt",
                        dest="version_file")
    parser.add_argument("--version", dest="version")
    parser.add_argument("--newMslName", default="YYYY", dest="msl_name")
    parser.add_argument(
        "--newMslNotes",
        default=("Provisional EC ##, Online meeting, July YYYY; "
                 "Email ratification March YYYY (MSL ###)"),
        dest="msl_notes",
    )
    result = Params(vars(parser.parse_args(argv)))
    result.debug = result.debug_on_error
    if result.tmi:
        result.verbose = True
    return result


# ---------------------------------------------------------------------------
# R and pandas compatibility helpers
# ---------------------------------------------------------------------------
#
# These small adapters centralize NA semantics, R-like vector recycling, frame
# appends, and POSIX-regex translation used throughout the conversion.


def is_na(value: Any) -> bool:
    if value is None:
        return True
    if isinstance(value, (list, tuple, pd.Series, pd.Index, pd.DataFrame)):
        return False
    try:
        result = pd.isna(value)
        return bool(result) if not hasattr(result, "__len__") else False
    except (TypeError, ValueError):
        return False


def text(value: Any, na: str = "NA") -> str:
    if is_na(value):
        return na
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    return str(value)


def first(value: Any, default: Any = pd.NA) -> Any:
    if isinstance(value, pd.DataFrame):
        return value.iloc[0, 0] if not value.empty else default
    if isinstance(value, (pd.Series, pd.Index, list, tuple)):
        return value[0] if len(value) else default
    return value


def values(value: Any) -> list[Any]:
    if isinstance(value, pd.DataFrame):
        return value.to_numpy().ravel().tolist()
    if isinstance(value, (pd.Series, pd.Index)):
        return value.tolist()
    if isinstance(value, (list, tuple)):
        return list(value)
    return [value]


def broadcast(*items: Any) -> Iterator[tuple[Any, ...]]:
    seqs = [values(item) for item in items]
    size = max((len(seq) for seq in seqs), default=0)
    for i in range(size):
        yield tuple(seq[i % len(seq)] if seq else pd.NA for seq in seqs)


def append_frames(left: pd.DataFrame, right: pd.DataFrame) -> pd.DataFrame:
    if left.empty:
        return right.copy()
    if right.empty:
        return left.copy()
    return pd.concat([left, right], ignore_index=True, sort=False)


def normalize_posix_regex(pattern: str) -> str:
    for source, target in {
        "[:alnum:]": "A-Za-z0-9", "[:alpha:]": "A-Za-z",
        "[:print:]": r"\x20-\x7E",
    }.items():
        pattern = pattern.replace(source, target)
    return pattern


# ---------------------------------------------------------------------------
# Input and output format adapters
# ---------------------------------------------------------------------------
#
# Database exports require a dedicated parser so unquoted \N/NULL, quoted empty
# strings, embedded control characters, and MySQL backslash escapes remain
# distinct. Excel readers preserve the fixed A:AO proposal-template contract.


def _parse_mysql_record(line: str) -> list[Any]:
    """Parse one MySQL tab export record, including MySQL backslash escapes."""
    escape_map = {
        "0": "\0", "b": "\b", "n": "\n", "r": "\r", "t": "\t",
        "Z": "\x1a", "\\": "\\", '"': '"', "'": "'",
    }
    fields: list[Any] = []
    value: list[str] = []
    quoted = False
    was_quoted = False
    index = 0
    line = line.rstrip("\r\n")
    while index < len(line):
        char = line[index]
        if char == '"' and not value and not quoted:
            quoted = True
            was_quoted = True
        elif char == '"' and quoted:
            # MySQL uses backslash escaping, but accept CSV-style doubled quotes too.
            if index + 1 < len(line) and line[index + 1] == '"':
                value.append('"')
                index += 1
            else:
                quoted = False
        elif char == "\\" and index + 1 < len(line):
            following = line[index + 1]
            # Preserve the unquoted NULL sentinel for field-final handling.
            if not was_quoted and not value and following == "N" and (
                    index + 2 == len(line) or line[index + 2] == "\t"):
                value.extend(["\\", "N"])
            else:
                value.append(escape_map.get(following, following))
            index += 1
        elif char == "\t" and not quoted:
            raw = "".join(value)
            fields.append(pd.NA if (not was_quoted and raw in {"\\N", "NULL"}) else raw)
            value = []
            was_quoted = False
        else:
            value.append(char)
        index += 1
    raw = "".join(value)
    fields.append(pd.NA if (not was_quoted and raw in {"\\N", "NULL"}) else raw)
    return fields


def read_mysql_tsv(path: str | Path) -> pd.DataFrame:
    """Read old and new MySQL TSV formats without conflating empty and NULL."""
    source = Path(path)
    with source.open("r", encoding="utf-8", newline="") as stream:
        header_line = stream.readline()
        if not header_line:
            return pd.DataFrame()
        header = [text(value, "") for value in _parse_mysql_record(header_line)]
        rows: list[list[Any]] = []
        for line_number, line in enumerate(stream, 2):
            parsed = _parse_mysql_record(line)
            if len(parsed) != len(header):
                raise RuntimeError(
                    f"ERROR: {source}:{line_number} has {len(parsed)} fields; expected {len(header)}"
                )
            rows.append(parsed)
    return pd.DataFrame.from_records(rows, columns=header)

def read_excel_raw(path: str | Path, sheet: str) -> pd.DataFrame:
    source = Path(path)
    engine = "xlrd" if source.suffix.lower() == ".xls" else "openpyxl"
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            frame = pd.read_excel(path, sheet_name=sheet, header=None,
                                  engine=engine, dtype=object)
    except ImportError as exc:
        raise RuntimeError("Legacy .xls input requires the xlrd package") from exc
    # readxl range=cell_cols("A:AO") retains 41 columns and overrides skip.
    frame = frame.reindex(columns=range(41))
    frame = frame.replace({
        "Please select": pd.NA, "[Please select]": pd.NA,
        "[Please\u00a0select]": pd.NA,
    })
    # The current template derives E1 from the workbook filename with a
    # formula. Some files contain a stale/error cached result, which pandas
    # sees instead of the formula Excel displays. Retain the formula marker so
    # qc_proposal can trust the already-validated filename-derived code.
    frame.attrs["formula_cells"] = {}
    if source.suffix.lower() != ".xls":
        workbook = load_workbook(source, read_only=True, data_only=False)
        try:
            formula = workbook[sheet]["E1"].value
            if isinstance(formula, str) and formula.startswith("="):
                frame.attrs["formula_cells"]["E1"] = formula
        finally:
            workbook.close()
    return frame


def write_tsv(frame: pd.DataFrame, path: str | Path, na_rep: str = "NA") -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False, encoding="utf-8",
                 na_rep=na_rep, lineterminator="\n")


def write_xlsx_sheets(sheets: Mapping[str, pd.DataFrame], path: str | Path) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        for name, frame in sheets.items():
            frame.to_excel(writer, sheet_name=str(name)[:31], index=False)


# ---------------------------------------------------------------------------
# Validator versioning and QC report collection
# ---------------------------------------------------------------------------
#
# All validation paths append to allErrorDf through log_error. The summary
# writer then creates the consolidated all-results, issues-only, TSV, and
# regression artifacts using the same public column contract as the R tools.


def load_version() -> None:
    if params.version:
        print(f"VERSION: {params.version}")
        return
    version_path = Path(params.version_file)
    if not version_path.exists():
        print(f"BUILD_VERSION_GIT: {version_path} does not exists. Rebuilding with ./version_git.sh")
        script = Path(__file__).resolve().with_name("version_git.sh")
        if script.exists():
            subprocess.run([str(script)], check=False)
    if not version_path.exists():
        raise RuntimeError(f"BUILD_VERSION_GIT: {version_path} still does not exist, after running ./version_git.sh")
    params.version = version_path.read_text(encoding="utf-8").split()[0]
    print(f"VERSION: {params.version}")


def _proposal_value(code: Any, column: str, default: Any = pd.NA) -> Any:
    if proposalsDf.empty or is_na(code) or code not in proposalsDf.index:
        return default
    return first(proposalsDf.loc[code, column], default) if column in proposalsDf else default


def log_error(code: Any, linenum: Any, action: Any, rank: Any, taxon: Any,
              levelStr: Any, errorCode: Any, errorStr: Any, notes: Any = "",
              actionOrder: Any = "") -> None:
    global allErrorDf
    rows: list[dict[str, Any]] = []
    for args in broadcast(code, linenum, action, rank, taxon, levelStr,
                          errorCode, errorStr, notes, actionOrder):
        (one_code, one_line, one_action, one_rank, one_taxon, one_level,
         one_error, one_message, one_notes, one_order) = args
        if is_na(one_order) or one_order == "":
            one_order = globals()["actionOrder"]
        try:
            row_value = int(float(str(one_line).split(":")[-1]))
        except (TypeError, ValueError):
            row_value = pd.NA
        docx = _proposal_value(one_code, "docx", pd.NA)
        xlsx = _proposal_value(one_code, "xlsx", pd.NA)
        rows.append({
            "subcommittee": _proposal_value(one_code, "subcommittee", pd.NA),
            "code": one_code, "docx": "MISSING" if is_na(docx) else docx,
            "xlsx": "MISSING" if is_na(xlsx) else xlsx, "row": row_value,
            "change": one_action, "rank": one_rank, "taxon": one_taxon,
            "level": one_level, "error": one_error, "message": one_message,
            "notes": one_notes,
            "scAbbrev": _proposal_value(one_code, "scAbbrev", pd.NA),
            "validator_version": params.get("version", pd.NA), "order": one_order,
        })
    if rows:
        allErrorDf = append_frames(allErrorDf, pd.DataFrame(rows))


def log_change_error(curChangeDf: pd.Series, levelStr: str, errorCode: str,
                     errorStr: str, notes: str = "") -> None:
    log_error(curChangeDf.get(".code"), curChangeDf.get(".linenum"),
              curChangeDf.get("change"), curChangeDf.get("rank"),
              curChangeDf.get(".changeTaxon"), levelStr, errorCode, errorStr,
              notes, curChangeDf.get(".order", ""))


def _sort_errors(frame: pd.DataFrame) -> pd.DataFrame:
    result = frame.copy()
    if result.empty:
        return result
    result["_row_n"] = pd.to_numeric(result["row"], errors="coerce")
    result["_order_n"] = pd.to_numeric(result["order"], errors="coerce")
    result = result.sort_values(["code", "_row_n", "_order_n"],
                                kind="stable", na_position="last")
    return result.drop(columns=["_row_n", "_order_n"])


def _issue_action_item(row: pd.Series) -> str:
    """Return the same action text assigned by qc_extract_issues.R."""
    error = text(row.get("error"), "").upper()
    actions = {
        "DOCX_MISSING": "Add the missing proposal .docx file (matching the proposal code) and rerun.",
        "DOCX_BAD_EXT": "Rename/convert the document to .docx and rerun.",
        "DOCX_UNKNONW_TEMPLATE": "Use the current ICTV .docx template and rerun.",
        "DOCX_TITLE_MISSING": "Add the proposal title to the .docx (required section) and rerun.",
        "DOCX_AUTHORS_MISSING": "Add authors to the .docx (required section) and rerun.",
        "DOCX_CORR_AUTHOR_MISSING": "Add a corresponding author to the .docx (required section) and rerun.",
        "DOCX_ABSTRACT_MISSING": "Add an abstract to the .docx (required section) and rerun.",
        "XLSX.MISSING_COLUMN": "Add the missing required column(s) using the current template and rerun.",
        "XLSX.INVALID_TERM": "Replace the invalid term with an allowed controlled-vocabulary value (template dropdowns) and rerun.",
        "XLSX.NO_SRC_DEST": "Fill in the 'current taxonomy' and/or 'proposed taxonomy' columns for this row and rerun.",
        "XLSX.CODE_BAD": "Fix the proposal code mismatch (filename vs sheet header) and rerun.",
        "XLSX.CODE_MISS": "Fix the proposal code (filename/header rows) and rerun.",
        "XLSX.EMPTY": "Add at least one change row (or remove the empty proposal file) and rerun.",
        "XLSX_EXTRA_SHEETS": "Remove extra worksheets or keep the template sheet as the first sheet; rerun.",
        "XLSX_MULTI_PROPOSAL": "Split into one proposal per .xlsx (one code per file) and rerun.",
        "XLSX_NOT_PROPOSAL": "Remove/rename the non-proposal .xlsx so it is not picked up by the scanner; rerun.",
        "RANK.UNK": "Fix the rank value to a supported rank (matches the template/reference) and rerun.",
        "ACTION.UNK": "Fix the 'Change' value to a supported action verb (matches the template) and rerun.",
        "MISSING_RANK_SUFFIX": "Fix the taxon name suffix to match the rank (per ICTV naming conventions) and rerun.",
        "SUFFIX_RANK_MISMATCH": "Fix the taxon name suffix to match the rank (per ICTV naming conventions) and rerun.",
        "TAXON_NAME_CASE": "Fix capitalization for the taxon name (non-species taxa capitalized) and rerun.",
        "SRC.LINEAGE_WRONG": "Update the 'current taxonomy' lineage to match the current MSL for this taxon and rerun.",
        "CREATE.PARENT_NO_EXIST": "Ensure the proposed parent taxon exists (or is created earlier in this proposal set) and rerun.",
        "MOVE.PARENT_NO_EXIST": "Ensure the destination parent taxon exists (or is created earlier) and rerun.",
    }
    if error in actions:
        return actions[error]
    if error.startswith("XLSX."):
        return "Fix the proposal .xlsx per the message/notes, then rerun."
    if error.startswith("DOCX_"):
        return "Fix the proposal .docx per the message/notes, then rerun."
    return "Review the message/notes; update the proposal row(s) and rerun the validator."


def write_error_summary(errorDf: pd.DataFrame, final: bool = False) -> None:
    """Write R-compatible all, issues-only, TSV, and regression summaries.

    ``final`` remains in the signature for compatibility. Per-subcommittee
    workbooks are intentionally no longer generated.
    """
    del final
    frame = errorDf.copy()
    for column in ERROR_COLUMNS + SUMMARY_METADATA_COLUMNS:
        if column not in frame:
            frame[column] = pd.NA
    frame["validator_version"] = params.get("version", pd.NA)
    known_columns = ERROR_COLUMNS + SUMMARY_METADATA_COLUMNS
    extra_columns = [column for column in frame.columns if column not in known_columns]
    output_columns = known_columns + extra_columns
    frame = frame.reindex(columns=output_columns)
    sorted_errors = _sort_errors(frame)

    pretty_rows: list[dict[str, Any]] = []
    previous = None
    for _, row in sorted_errors.iterrows():
        code = row.get("code")
        if not is_na(code) and code != previous:
            previous = code
            pretty_rows.append({"subcommittee": row.get("subcommittee")})
            pretty_rows.append({
                "subcommittee": row.get("subcommittee"), "code": code,
                "xlsx": row.get("docx") if is_na(row.get("xlsx")) else row.get("xlsx"),
            })
        pretty_rows.append(row.to_dict())
    pretty = pd.DataFrame(pretty_rows).reindex(columns=output_columns)
    public_columns = [
        column for column in output_columns
        if column not in {"subcommittee", "code", "docx"}
    ]

    xlsx_path = Path(params.out_dir) / params.qc_summary_fname
    write_xlsx_sheets({"QC_Report": pretty.loc[:, public_columns]}, xlsx_path)
    if params.verbose or params.tmi:
        print(f"Wrote: {xlsx_path} with 1 sheets")
        print(f"     sheet 'QC_Report {len(pretty)} rows)")

    issue_mask = sorted_errors["level"].astype("string").isin(["ERROR", "WARNING"])
    issues = sorted_errors.loc[issue_mask].copy()
    issues["action_item"] = issues.apply(_issue_action_item, axis=1)
    issue_columns = list(public_columns)
    insert_at = issue_columns.index("notes") + 1 if "notes" in issue_columns else len(issue_columns)
    issue_columns.insert(insert_at, "action_item")
    issues_path = Path(params.out_dir) / "QC.pretty_summary.issues.xlsx"
    write_xlsx_sheets({"QC_Issues": issues.loc[:, issue_columns]}, issues_path)
    if params.verbose or params.tmi:
        print(f"Wrote: {issues_path} with {len(issues)} issue rows")

    summary_path = Path(params.out_dir) / params.qc_summary_tsv_fname
    write_tsv(sorted_errors.loc[:, output_columns], summary_path)
    if params.verbose or params.tmi:
        print(f"Wrote: {summary_path} ({len(frame)} rows)")
    else:
        print(f"Wrote: OUTDIR/{params.qc_summary_tsv_fname} ({len(frame)} rows)")

    regression_columns = [column for column in output_columns if "version" not in column]
    regression_path = Path(params.out_dir) / params.qc_regression_tsv_fname
    write_tsv(sorted_errors.loc[:, regression_columns], regression_path)
    if params.verbose or params.tmi:
        print(f"Wrote: {regression_path} ({len(frame)} rows)")


# ---------------------------------------------------------------------------
# Proposal template schema and validation rules
# ---------------------------------------------------------------------------
#
# The column maps support the historical and current workbook layouts. The
# controlled vocabularies and validation tables below are data-driven mirrors
# of the R configuration rather than action-specific conditionals.

xlsx_change_ranks = [
    "Realm", "Subrealm", "Kingdom", "Subkingdom", "Phylum", "Subphylum",
    "Class", "Subclass", "Order", "Suborder", "Family", "Subfamily",
    "Genus", "Subgenus", "Species",
]
xlsx_change_src_colnames = [f"src{name}" for name in xlsx_change_ranks]
xlsx_change_dest_colnames = [name.lower() for name in xlsx_change_ranks]
xlsx_change_other = [
    "exemplarAccession", "exemplarName", "Abbrev", "exemplarIsolate",
    "genomeCoverage", "molecule", "hostSource", "change", "rank", "comments",
]
xlsx_change_colnames = (
    xlsx_change_src_colnames + xlsx_change_dest_colnames + xlsx_change_other
)

XLSX_V1_CHANGE_COLS = list(range(0, 15)) + list(range(16, 31)) + list(range(31, 41))
XLSX_V2_CHANGE_COLS = list(range(1, 16)) + list(range(16, 31)) + list(range(31, 41))
XLSX_2023_CHANGE_COLS = list(range(0, 15)) + list(range(15, 30)) + list(range(30, 40))

XLSX_2023_ROW3 = [
    "CURRENT TAXONOMY", *([pd.NA] * 14), "PROPOSED TAXONOMY",
    *([pd.NA] * 14), "DESCRIPTIVES", *([pd.NA] * 6), "ACTION", pd.NA,
    "COMMENTS",
]
XLSX_2023_ROW4 = xlsx_change_ranks + xlsx_change_ranks + [
    "Exemplar GenBank Accession Number", "Exemplar virus name",
    "Virus name abbreviation", "Exemplar isolate designation",
    "Genome coverage", "Genome composition", "Host/Source", "Change",
    "Proposed Rank", "Comments",
]

VALUE_VALIDATION: list[dict[str, Any]] = []
for col in [c for c in xlsx_change_src_colnames + xlsx_change_dest_colnames
            if "pecies" not in c.lower()]:
    VALUE_VALIDATION.append({
        "col": col, "regex": r"([^[:alnum:]-]+)", "type": "replace",
        "class": "ERROR", "code": "XLSX.NON_ALPHA-NUMERIC",
        "warn": "non-(AlphaNumeric or hyphen) characters removed", "replace": "",
    })
for col in ["srcSpecies", "species"]:
    VALUE_VALIDATION.append({
        "col": col, "regex": r"([^[:alnum:] -]+)", "type": "replace",
        "class": "ERROR", "code": "XLSX.NON_ALPHA-NUMERIC-SPACE",
        "warn": "non-(AlphaNumeric or hyphen) characters removed", "replace": "",
    })
VALUE_VALIDATION.extend([
    {"col": "species", "regex": r"^([[:alnum:]-]+ [[:alnum:]-]+)$",
     "type": "required", "class": "ERROR", "code": "XLSX.SPECIES_BAD_NAME",
     "warn": "Species name must be 'genus[space]species' binomial naming", "replace": ""},
    {"col": "change", "regex": r"([^[:alnum:] ;-]+)", "type": "replace",
     "class": "INFO", "code": "XLSX.CHANGE_CV_REMOVE_NON_ALPHA-NUMERIC-SPACE-SEMI",
     "warn": "Change col:non-(AlphaNumeric,hyphen,space,dash) characters removed", "replace": " "},
    {"col": "Abbrev", "regex": r"([^[:alnum:];,_/. -]+)", "type": "replace",
     "class": "ERROR", "code": "XLSX.ABBREV.REMOVE_ILLEGAL_CHARS",
     "warn": "Should be semicolon-separated list", "replace": ""},
    {"col": "exemplarAccession", "regex": r"( +and +)", "type": "replace",
     "class": "WARNING", "code": "XLSX.ACCESSION.REPLACE_AND_WITH_SEMI",
     "warn": "Should be semicolon-separated list, with optional colon-separated name prefixes", "replace": ";"},
    {"col": "exemplarAccession", "regex": r"([A-Z0-9.]+) *\(([^\)]+) segment\)",
     "type": "replace", "class": "ERROR", "code": "XLSX.ACCESSION.REPLACE_PAREN_WITH_PREFIX",
     "warn": "Should be semicolon-separated list, with optional colon-separated name prefixes", "replace": r"\2:\1"},
    {"col": "exemplarAccession", "regex": r"(,)", "type": "replace",
     "class": "ERROR", "code": "XLSX.ACCESSION.CONVERT_COMMA_TO_SEMICOLON",
     "warn": "commas: should be a semicolon-separated list, with optional colon-separated name prefixes: convert commas to semicolons", "replace": ";"},
    {"col": "exemplarAccession", "regex": r"([^[:alnum:]:;_/. -]+)",
     "type": "replace", "class": "ERROR", "code": "XLSX.ACCESSION.REMOVED_ILLEGAL_CHARS",
     "warn": "unexpected characters in accession list: should be a semi-colon separated list, with optional colon-separated name prefixes", "replace": ""},
])
for col in ["exemplarIsolate", "genomeCoverage", "molecule", "hostSource"]:
    VALUE_VALIDATION.append({
        "col": col, "regex": r"([^[:print:]]+)", "type": "replace",
        "class": "ERROR", "code": "XLSX.NON_PRINTABLE_REMOVED",
        "warn": "unprintable/non-ASCII removed", "replace": " ",
    })


ACTION_MAP = {
    "create": "new", "create new": "new", "rename": "rename",
    "abolish": "abolish", "move": "move", "move; rename": "move",
    "promote": "promote", "demote": "demote", "merge": "merge",
    "split": "split",
}
XLSX2DB_MAP = {
    "exemplarAccession": "genbank_accession_csv",
    "exemplarName": "exemplar_name", "Abbrev": "abbrev_csv",
    "exemplarIsolate": "isolate_csv", "genomeCoverage": "genome_coverage",
    "molecule": "molecule_id", "hostSource": "host_source",
    "rank": "level_id", "comments": "notes",
}


# ---------------------------------------------------------------------------
# Taxonomy lineage helpers
# ---------------------------------------------------------------------------
#
# Proposal rows represent a lineage as one column per rank. These helpers
# derive the active name, rank, parent, and semicolon-delimited lineage used by
# both QC and the working MSL.


def taxon_get_name(realmSpecies: Sequence[Any] | pd.Series) -> Any:
    items = list(realmSpecies)[:15]
    present = [value for value in items if not is_na(value) and str(value) != ""]
    return present[-1] if present else pd.NA


def taxon_get_parent_name(realmSpecies: Sequence[Any] | pd.Series) -> Any:
    items = list(realmSpecies)[:15]
    present = [value for value in items if not is_na(value) and str(value) != ""]
    if len(present) > 1:
        return present[-2]
    if len(present) == 1:
        return params.msl_name
    return pd.NA


def taxon_get_lineage(realmSpecies: Sequence[Any] | pd.Series,
                       masks: Sequence[str] = ()) -> Any:
    excluded = {str(item).lower() for item in masks}
    kept = [value for rank, value in zip(xlsx_change_ranks, list(realmSpecies)[:15])
            if rank.lower() not in excluded and not is_na(value) and str(value) != ""]
    return ";".join(str(value) for value in kept) if kept else pd.NA


def taxon_get_parent_lineage(realmSpecies: Sequence[Any] | pd.Series) -> Any:
    rank = taxon_get_rank(realmSpecies)
    return taxon_get_lineage(realmSpecies, [rank]) if not is_na(rank) else pd.NA


def taxon_get_rank(realmSpecies: Sequence[Any] | pd.Series) -> Any:
    items = list(realmSpecies)[:15]
    present = [i for i, value in enumerate(items)
               if not is_na(value) and str(value) != ""]
    return xlsx_change_ranks[present[-1]].lower() if present else pd.NA


def diff_strings(s1: Any, s2: Any) -> str:
    left, right = text(s1), text(s2)
    prefix = 0
    while prefix < min(len(left), len(right)) and left[prefix] == right[prefix]:
        prefix += 1
    suffix = 0
    while (suffix < len(left) - prefix and suffix < len(right) - prefix
           and left[-suffix - 1] == right[-suffix - 1]):
        suffix += 1
    if prefix + suffix < 4:
        return f"[{left}//{right}]"
    left_mid = left[prefix:len(left) - suffix if suffix else None]
    right_mid = right[prefix:len(right) - suffix if suffix else None]
    return f"{left[:prefix]}[{left_mid}//{right_mid}]{left[len(left)-suffix:] if suffix else ''}"


def diff_lineages(lin1: Any, lin2: Any) -> str:
    left = text(lin1).split(";")
    right = text(lin2).split(";")
    size = max(len(left), len(right))
    left.extend([""] * (size - len(left)))
    right.extend([""] * (size - len(right)))
    parts = [a if a == b else diff_strings(a, b) for a, b in zip(left, right)]
    return ("[identical]" if left == right else "") + ";".join(parts)


# ---------------------------------------------------------------------------
# GenBank accession parsing and reference lookup
# ---------------------------------------------------------------------------
#
# species_isolates is checked for every listed owner, including rows labelled
# abolished. abolished_accessions.tsv supplements it with historical removals
# that predate tracking in species_isolates. Source taxnode exclusions permit
# an action to retain its own accession while still exposing other owners.


def normalize_accession_value(accession: Any) -> Any:
    if is_na(accession):
        return pd.NA
    value = str(accession).strip()
    if not value:
        return pd.NA
    value = value.rsplit(":", 1)[-1].strip()
    return re.sub(r" \([0-9]+\.[0-9]+\)$", "", value).strip() or pd.NA


def _expand_accession_range(token: str) -> list[str]:
    """Expand AB000001-AB000003 (or AB000001-000003) with a safe limit."""
    value = normalize_accession_value(token)
    if is_na(value):
        return []
    match = re.fullmatch(r"([A-Za-z_]+)(\d+)(?:\.\d+)?\s*[-–—]\s*([A-Za-z_]*)(\d+)(?:\.\d+)?", str(value))
    if not match:
        return [str(value)]
    prefix1, start_s, prefix2, end_s = match.groups()
    prefix2 = prefix2 or prefix1
    start, end = int(start_s), int(end_s)
    if prefix1 != prefix2 or end < start or end - start > 10000:
        return [str(value)]
    width = max(len(start_s), len(end_s))
    return [f"{prefix1}{number:0{width}d}" for number in range(start, end + 1)]


def split_accessions(accession: Any) -> list[str]:
    if is_na(accession):
        return []
    result: list[str] = []
    for token in re.split(r"\s*;\s*", str(accession)):
        if token:
            result.extend(_expand_accession_range(token))
    return list(dict.fromkeys(result))


def parseGenbankAccessionsV1(frame: pd.DataFrame) -> pd.DataFrame:
    return parseGenbankAccessions(
        frame, "genbank_accessions",
        ["isolate_id", "taxnode_id", "species_name", "isolate_type", "isolate_names"],
    )["validAccessions"]


def parseGenbankAccessions(frame: pd.DataFrame, accessionCol: str,
                           extraCols: Sequence[str]) -> dict[str, pd.DataFrame]:
    if accessionCol not in frame:
        raise RuntimeError("The specified accession column does not exist in the data frame.")
    missing = [column for column in extraCols if column not in frame]
    if missing:
        raise RuntimeError("The following extra columns are missing in the data frame: " + ", ".join(missing))
    valid: list[dict[str, Any]] = []
    invalid: list[dict[str, Any]] = []
    for _, source in frame.iterrows():
        raw = source.get(accessionCol)
        if is_na(raw) or str(raw).strip() == "":
            continue
        index = 0
        for item in str(raw).split(";"):
            index += 1
            item = item.strip()
            if not re.fullmatch(r"[a-zA-Z0-9._-]*:?[^:()]+(?: \([0-9]+\.[0-9]+\))?", item):
                row = {column: source.get(column) for column in extraCols}
                row.update({accessionCol: item, "accession_index": index})
                invalid.append(row)
                continue
            segment = item.split(":", 1)[0].strip() if ":" in item else pd.NA
            for accession in _expand_accession_range(item):
                row = {column: source.get(column) for column in extraCols}
                range_match = re.search(r"\(([0-9]+\.[0-9]+)\)$", item)
                row.update({"segment_name": segment, "accession": accession,
                            "accession_range": range_match.group(1) if range_match else pd.NA,
                            "accession_index": index})
                valid.append(row)
    return {"validAccessions": pd.DataFrame(valid),
            "invalidAccessions": pd.DataFrame(invalid)}


def find_vmr_accession_matches(accession: Any,
                               exclude_taxnode_ids: Sequence[Any] = ()) -> list[str]:
    accessions = split_accessions(accession)
    if vmrAccessions.empty or not accessions:
        return []
    matches = vmrAccessions["accession"].isin(accessions)
    excluded = {text(value) for value in exclude_taxnode_ids if not is_na(value)}
    if excluded and "taxnode_id" in vmrAccessions:
        matches &= ~vmrAccessions["taxnode_id"].map(text).isin(excluded)
    return list(dict.fromkeys(vmrAccessions.loc[matches, "species_name"].dropna().astype(str)))


def load_abolished_accessions() -> None:
    global abolishedAccessions
    candidates = [Path(params.abolished_accessions_fname),
                  Path(params.ref_dir) / params.abolished_accessions_fname,
                  Path(__file__).resolve().parent / params.abolished_accessions_fname]
    filename = next((candidate for candidate in candidates if candidate.exists()), None)
    if filename is None:
        if params.verbose:
            print(f"# SKIP ABOLISHED ACCESSIONS LOAD: {candidates[0]}")
        return
    frame = read_mysql_tsv(filename)
    if "accession" not in frame:
        raise RuntimeError(f"ERROR: {filename} is missing required column: accession")
    for column in ["species", "virus_name"]:
        if column not in frame:
            frame[column] = pd.NA
    frame["accession"] = frame["accession"].map(normalize_accession_value)
    abolishedAccessions = frame.dropna(subset=["accession"])[["accession", "species", "virus_name"]].drop_duplicates()
    if params.verbose:
        print(f"AbolishedAccessions:  {abolishedAccessions.shape} from {filename}")


def find_abolished_accession_matches(accession: Any) -> pd.DataFrame:
    if abolishedAccessions.empty:
        return abolishedAccessions.iloc[0:0].copy()
    return abolishedAccessions.loc[
        abolishedAccessions["accession"].isin(split_accessions(accession))
    ].drop_duplicates()


def format_abolished_accession_matches(matches: pd.DataFrame) -> str:
    labels: list[str] = []
    for _, row in matches.iterrows():
        label = row.get("species") if not is_na(row.get("species")) and row.get("species") != "" else row.get("virus_name")
        if not is_na(label) and label != "" and str(label) not in labels:
            labels.append(str(label))
    return "; ".join(labels)


def is_pending_accession(accession: Any) -> bool:
    return not is_na(accession) and bool(re.search(r"\bpending\b", str(accession), re.I))


# ---------------------------------------------------------------------------
# Reference MSL construction and controlled-vocabulary loading
# ---------------------------------------------------------------------------
#
# load_reference builds curMSL/newMSL from database exports, loads rank,
# molecule, host/source vocabularies, and expands species isolate accessions to
# one normalized accession per row.


def create_new_msl(current: pd.DataFrame, prev_msl: int, dest_msl: int,
                   taxnode_delta: int) -> pd.DataFrame:
    result = current.loc[pd.to_numeric(current["msl_release_num"], errors="coerce") == prev_msl].copy(deep=True)
    result[".prev_taxnode_id"] = result["taxnode_id"]
    result["prev_proposals"] = pd.NA
    for column in ["exemplar_name", "genome_coverage", "host_source"]:
        if column not in result:
            result[column] = pd.NA
    for column in ["taxnode_id", "parent_id", "tree_id"]:
        result[column] = pd.to_numeric(result[column], errors="coerce").astype("Int64") + int(taxnode_delta)
    result["msl_release_num"] = dest_msl
    for column in [c for c in result if re.match(r"^(in|out)_", c)]:
        result[column] = pd.NA
    for column in [".emptyReported", ".otherLineage", ".otherLineageProposal",
                   ".otherLineageAction"]:
        result[column] = pd.NA
    root = result["rank"].astype(str) == "tree"
    result.loc[root, "name"] = params.msl_name
    result.loc[root, "notes"] = params.msl_notes
    return result.reset_index(drop=True)


def _clean_cv_name(value: Any) -> str:
    return re.sub(r"[^A-Za-z]+", " ", text(value, "")).strip()


def load_reference() -> None:
    global oldMSLs, curMSL, newMSL, dbCvList, dbCvMapList, cvList
    global scAbbrevNameMap, vmrDf, vmrAccessions
    taxonomy_filename = Path(params.ref_dir) / params.db_taxonomy_node_fname
    taxonomy = read_mysql_tsv(taxonomy_filename)
    numeric_columns = ["taxnode_id", "parent_id", "tree_id", "msl_release_num",
                       "level_id", "ictv_id", "molecule_id"]
    for column in numeric_columns:
        if column in taxonomy:
            taxonomy[column] = pd.to_numeric(taxonomy[column], errors="coerce").astype("Int64")
    # data.table fread(key="taxnode_id") sorts the R reference at load time;
    # preserving that order keeps msl.tsv useful for byte-oriented regressions.
    taxonomy = taxonomy.sort_values("taxnode_id", kind="stable").reset_index(drop=True)
    print(f"Previous taxa: {taxonomy.shape[0]} {taxonomy.shape[1]} from {taxonomy_filename}")
    if "host_source" not in taxonomy:
        print("WARNING: no host_source column in taxonomy_node dump!!! (Adding)")
        taxonomy["host_source"] = pd.NA
    last_msl = int(taxonomy["msl_release_num"].max())
    oldMSLs = taxonomy.loc[taxonomy["msl_release_num"] < last_msl].copy()
    curMSL = taxonomy.loc[taxonomy["msl_release_num"] == last_msl].copy().reset_index(drop=True)
    curMSL["out_updated"] = False
    newMSL = create_new_msl(curMSL, last_msl, last_msl + 1, int(params.taxnode_delta))

    rank_filename = Path(params.ref_dir) / params.db_rank_fname
    molecule_filename = Path(params.ref_dir) / params.db_molecule_fname
    host_filename = Path(params.ref_dir) / params.db_host_source_fname
    rank_cv = read_mysql_tsv(rank_filename)
    molecule_cv = read_mysql_tsv(molecule_filename)
    host_cv = read_mysql_tsv(host_filename)
    for frame in [rank_cv, molecule_cv]:
        if "id" in frame:
            frame["id"] = pd.to_numeric(frame["id"], errors="coerce").astype("Int64")
    dbCvList = {"rank": rank_cv, "molecule": molecule_cv,
                "hostSource": host_cv.get("host_source", pd.Series(dtype=object))}
    dbCvMapList = {
        "rank": dict(zip(rank_cv["name"].astype(str), rank_cv["id"])),
        "molecule": dict(zip(molecule_cv["abbrev"].astype(str), molecule_cv["id"])),
        "hostSource": {value: value for value in host_cv["host_source"].dropna().astype(str)},
    }
    if params.verbose:
        print(f"RankCV: {rank_cv.shape} from {rank_filename}")
        print(f"MoleculeCV: {molecule_cv.shape} from {molecule_filename}")
        print(f"HostSourceCV: {host_cv.shape} from {host_filename}")

    template_filename = Path(params.ref_dir) / params.template_xlsx_fname
    template = read_excel_raw(template_filename, params.template_xlsx_sheet)
    if params.verbose:
        print(f"ProposalTemplate[{params.template_xlsx_sheet}]: {template.shape} from {template_filename}")
    template.iloc[0, 5] = "Subcommittee Abbrev"
    template.iloc[0, 6] = "Subcommittee Name"
    name_map = {
        "Genome coverage": "genomeCoverage", "Genome composition": "molecule",
        "Host Source": "hostSource", "Change": "change",
        "Proposed Rank": "rank", "Subcommittee Abbrev": "scAbbrev",
        "Subcommittee Name": "scName",
    }
    loaded: dict[str, list[Any]] = {}
    for column in template.columns:
        cv_name = _clean_cv_name(template.iloc[0, column])
        mapped = name_map.get(cv_name)
        if not mapped:
            continue
        vals = []
        for value in template.iloc[1:, column].tolist():
            if is_na(value):
                continue
            cleaned = re.sub(r"[^\w\s\-+;()/.,]", " ", str(value).replace("\u00a0", " ")).strip()
            if cleaned:
                vals.append(cleaned)
        loaded[mapped] = list(dict.fromkeys(vals))
    loaded[".change2action"] = ACTION_MAP.copy()
    for required in ["change", "rank", "genomeCoverage", "molecule", "hostSource"]:
        loaded.setdefault(required, [])
    loaded["change"] = [v for v in loaded["change"] if re.sub(r"\W", "", v).lower() != "pleaseselect"]
    loaded["rank"] = [v for v in loaded["rank"] if re.sub(r"\W", "", v).lower() != "pleaseselect"]
    sc_abbrev = loaded.pop("scAbbrev", [])
    sc_name = loaded.pop("scName", [])
    scAbbrevNameMap = dict(zip(sc_abbrev, sc_name))
    cvList = loaded
    for value in cvList.get("molecule", []):
        compact = value.replace(" ", "")
        if value not in dbCvMapList["molecule"] and compact in dbCvMapList["molecule"]:
            dbCvMapList["molecule"][value] = dbCvMapList["molecule"][compact]
    dbCvMapList["molecule"]["multiple"] = pd.NA
    cvList["hostSource"] = list(dict.fromkeys(
        cvList.get("hostSource", []) + host_cv["host_source"].dropna().astype(str).tolist()
    ))

    vmr_filename = Path(params.ref_dir) / params.vmr_fname
    if vmr_filename.exists():
        vmrDf = read_mysql_tsv(vmr_filename)
        parsed = parseGenbankAccessions(
            vmrDf, "genbank_accessions",
            ["isolate_id", "taxnode_id", "species_name", "isolate_type", "isolate_names"],
        )
        vmrAccessions = parsed["validAccessions"]
        invalid = parsed["invalidAccessions"]
        if not invalid.empty:
            for _, bad in invalid.iterrows():
                print(f"ERROR: VMR accession badly formated: species:{bad.get('species_name')} iso_id:{bad.get('isolate_id')} value: {bad.get('genbank_accessions')}")
            raise RuntimeError(f"ERROR: {vmr_filename} contains {len(invalid)} badly formatted accessions")
    else:
        print(f"# SKIP VMR LOAD: {vmr_filename}")


# ---------------------------------------------------------------------------
# Proposal discovery and input-file classification
# ---------------------------------------------------------------------------
#
# The scanner recursively classifies proposal, document, supplemental, and
# unknown files; associates them by proposal code; and repairs spaces in names
# when the directory is writable.


def _append_scan_error(row: Mapping[str, Any]) -> None:
    global allErrorDf
    allErrorDf = append_frames(allErrorDf, pd.DataFrame([dict(row)]))


def _proposal_code(filename: str) -> Any:
    match = re.match(r"^([0-9]+\.[0-9]+[A-Z]X*)\.", Path(filename).stem)
    return match.group(1) if match else pd.NA


def _repair_spaced_filenames(root: Path) -> Path:
    """Rename files containing spaces when collision-free; otherwise warn."""
    files = [root] if root.is_file() else sorted(p for p in root.rglob("*") if p.is_file())
    for source in files:
        if " " not in source.name:
            continue
        target = source.with_name(re.sub(r" +", "_", source.name))
        code = _proposal_code(source.name)
        if not target.exists():
            try:
                source.rename(target)
            except OSError as exc:
                _append_scan_error({
                    "code": code, "level": "WARNING", "error": "FILENAME_SPACES",
                    "message": "filename contains a space: please replace with _ or -",
                    "notes": f"Could not rename {source.name}: {exc}",
                })
                continue
            _append_scan_error({
                "code": code, "level": "WARNING", "error": "FILENAME_SPACES_RENAMED",
                "message": "filename contains a space: replaced with _",
                "notes": f"{source.name} -> {target.name}",
                "docx": target.name if target.suffix.lower() in {".doc", ".docx"} else pd.NA,
                "xlsx": target.name if target.suffix.lower() in {".xls", ".xlsx"} else pd.NA,
            })
            if params.verbose:
                print(f"WARNING: renamed file containing spaces: {source} -> {target}")
            if source == root:
                root = target
        else:
            _append_scan_error({
                "code": code, "level": "WARNING", "error": "FILENAME_SPACES",
                "message": "filename contains a space: please replace with _ or -",
                "notes": f"Could not rename {source.name}: {target.name} already exists",
            })
    return root


def _input_file_record(path: Path) -> dict[str, Any]:
    stem = re.sub(r"\.(doc|xls)x?$", "", path.name, flags=re.I)
    code = _proposal_code(path.name)
    abbreviation = pd.NA
    if not is_na(code):
        match = re.match(r"^[0-9]{4}\.[0-9]{3}([A-Z])X*$", str(code))
        abbreviation = match.group(1) if match else pd.NA
    return {"docpath": str(path), "path": str(path.parent), "file": path.name,
            "basename": stem, "code": code, "scAbbrev": abbreviation}


def scan_for_proposals() -> pd.DataFrame:
    """Discover proposal files and return one metadata row per proposal code."""
    global proposalsDf, allErrorDf
    root = Path(params.proposals_dir)
    if params.processing_mode == "final":
        filename_regex = re.compile(r"^[0-9]{4}\.[0-9]{3}[A-Z]X*\.[^ ]*", re.I)
        filename_message = "final:####[A-Z].###[A-Z][X].____"
    elif params.processing_mode == "draft":
        filename_regex = re.compile(r"^[0-9]{4}\.[0-9]{3}[A-Z]X*\.[A-Za-z]+\.v[0-9]+\.[^ ]*", re.I)
        filename_message = "draft:####[A-Z].###[A-Z][X].[A-Z]+*.v#.____"
    elif params.processing_mode == "validate":
        filename_regex = re.compile(r"^.*", re.I)
        filename_message = "validate:*.____"
    else:
        raise RuntimeError(f"ERROR: --mode='{params.processing_mode}' is not a valid option: validate, draft, or final")
    if params.verbose:
        print(f"(PROCESSING) MODE      : {params.processing_mode}")
    if not root.exists():
        _append_scan_error({"level": "ERROR", "error": "PROPOSAL_DIR_NO_EXIST",
                            "message": "path not found", "notes": f"Input folder='{root}/'"})
        write_error_summary(allErrorDf)
        raise RuntimeError(f"# ERROR: PROPOSAL_DIR_NO_EXIST:{root}")

    root = _repair_spaced_filenames(root)
    candidates = [root] if root.is_file() else sorted(
        p for p in root.rglob("*")
        if p.is_file() and re.search(r"\.(doc|xls)x?$", p.name, re.I)
        and not p.name.startswith(("~", "."))
    )
    records = [_input_file_record(path) for path in candidates]
    inputs = pd.DataFrame(records)
    if inputs.empty:
        _append_scan_error({"level": "ERROR", "error": "NO_INPUT_FILES",
                            "message": "found no .xls[x] or .doc[x] files",
                            "notes": f"Input folder='{root}/'"})
        write_error_summary(allErrorDf)
        raise RuntimeError(f"# ERROR: NO_INPUT_FILES found in {root}")
    if params.tmi:
        print(f"# xls|doc(x) files found: N={len(inputs)}")
        for item in inputs["docpath"]:
            print(f"\t{item}")

    def filter_inputs(mask: pd.Series, level: str, error: str, message: str) -> None:
        nonlocal inputs
        global allErrorDf
        if mask.any():
            errors = inputs.loc[mask].copy()
            errors["level"], errors["error"], errors["message"] = level, error, message
            allErrorDf = append_frames(allErrorDf, errors)
            write_error_summary(allErrorDf)
            inputs = inputs.loc[~mask].copy()

    if params.processing_mode in {"draft", "final"}:
        mask = inputs["file"].str.contains(r"\.Ud\.", case=False, regex=True)
        filter_inputs(mask, "INFO", "IGNORE_FNAME_.Ud.", "All files with '.Ud.' in filename are ignored")
    mask = inputs["file"].map(
        lambda value: bool(re.search(params.infile_suppl_pat, str(value), re.I))
    )
    filter_inputs(mask, "INFO", "IGNORE_FNAME_SUPPL",
                  f"All files with '{params.infile_suppl_pat}' in filename are ignored")
    mask = inputs["file"].str.contains("appendix", case=False, regex=False)
    filter_inputs(mask, "INFO", "IGNORE_FNAME_APPENDIXL",
                  "All files with 'appendix' in filename are ignored")

    docxs = inputs.loc[inputs["file"].str.contains(r"\.docx?$", case=False, regex=True)].copy()
    docxs = docxs.rename(columns={"docpath": "docxpath", "file": "docx"})
    xlsxs = inputs.loc[inputs["file"].str.contains(r"\.xlsx?$", case=False, regex=True)].copy()
    xlsxs = xlsxs.rename(columns={"docpath": "xlsxpath", "file": "xlsx"})

    for frame, file_col, error_code in [
        (docxs, "docx", "DUPCODE.DOCX"), (xlsxs, "xlsx", "XLSX_DUPCODE")
    ]:
        duplicate = frame["code"].notna() & frame["code"].duplicated(keep=False)
        if duplicate.any():
            errors = frame.loc[duplicate, ["scAbbrev", "code", file_col]].copy()
            errors["level"] = "ERROR"
            errors["error"] = error_code
            errors["message"] = "duplicate proposal ID"
            allErrorDf = append_frames(allErrorDf, errors)
            write_error_summary(allErrorDf)
            if frame is xlsxs:
                raise RuntimeError("ERROR: can not proceed with duplicated proposal IDs: " +
                                   ",".join(frame.loc[duplicate, "code"].astype(str).unique()))

    for frame, file_col, kind in [(docxs, "docx", "DOCX"), (xlsxs, "xlsx", "XLSX")]:
        if params.processing_mode in {"draft", "final"}:
            good = frame[file_col].map(lambda value: bool(filename_regex.match(str(value))))
            if (~good).any():
                errors = frame.loc[~good, ["scAbbrev", "code", file_col]].copy()
                errors["level"] = "ERROR"
                errors["error"] = f"{kind}.BAD_FILENAME_FORMAT"
                errors["message"] = f"Should be '{filename_message}.{kind.lower()}'"
                errors["notes"] = "####[A-Z]=year/study_section, ###[A-Z]=index/type, [A-Z]+=status, v#=version"
                allErrorDf = append_frames(allErrorDf, errors)

    if xlsxs.empty:
        _append_scan_error({"level": "ERROR", "error": "NO_INPUT_FILES",
                            "message": "found no .xls[x] files that passed filename QC",
                            "notes": f"Input folder='{root}/'"})
        write_error_summary(allErrorDf)
        raise RuntimeError(f"# ERROR: NO_INPUT_FILES found in {root}")

    # Give files without a proposal code a stable independent key.
    for frame in [docxs, xlsxs]:
        for idx in frame.index[frame["code"].isna()]:
            frame.at[idx, "code"] = str(idx + 1)
    merge_keys = ["code", "path", "scAbbrev"]
    proposals = pd.merge(xlsxs, docxs, on=merge_keys, how="outer",
                         suffixes=(".xlsx", ".docx"), sort=True)
    proposals["basename"] = proposals.get("basename.xlsx").combine_first(proposals.get("basename.docx"))
    proposals = proposals.drop(columns=[c for c in ["basename.xlsx", "basename.docx"] if c in proposals])
    proposals["cleanbase"] = proposals["basename"].map(
        lambda value: re.sub(r"^([0-9]+\.[0-9]+[A-Z]X*)(\.[A-Z]+)(\.v[0-9]+)*(\.fix)*(\..*)$", r"\1\5", str(value))
    )
    missing = proposals["xlsx"].isna()
    if missing.any():
        errors = proposals.loc[missing, ["code", "docx"]].copy()
        errors["xlsx"], errors["row"], errors["level"] = pd.NA, pd.NA, "ERROR"
        errors["error"], errors["message"] = "XLSX.MISSING", "DOCX has no matching XLSX"
        errors["notes"] = "Suggestions: contact corresponding author"
        allErrorDf = append_frames(allErrorDf, errors)
    proposals["scAbbrev"] = proposals["code"].map(
        lambda code: (re.match(r"^[0-9]{4}\.[0-9]{3}([A-Z])X*", str(code)).group(1)
                      if re.match(r"^[0-9]{4}\.[0-9]{3}([A-Z])X*", str(code)) else pd.NA)
    )
    bad = ~proposals["scAbbrev"].isin(scAbbrevNameMap)
    if bad.any() and params.processing_mode in {"draft", "final"}:
        errors = proposals.loc[bad, ["scAbbrev", "code", "xlsx", "docx"]].copy()
        errors["level"], errors["error"] = "ERROR", "CODE_BAD_SC_ABBREV"
        errors["message"] = "Last letter of CODE not a valid ICTV Subcommittee letter"
        errors["notes"] = errors["scAbbrev"].map(
            lambda value: f"'{text(value)}' not in [{','.join(scAbbrevNameMap)}]"
        )
        allErrorDf = append_frames(allErrorDf, errors)
    proposals["subcommittee"] = proposals.apply(
        lambda row: (scAbbrevNameMap.get(row["scAbbrev"])
                     if row["scAbbrev"] in scAbbrevNameMap
                     else ("unspecified" if is_na(row["scAbbrev"])
                           else f"unknown-{row['scAbbrev']}")), axis=1
    )
    proposalsDf = proposals.set_index("code", drop=False)
    return proposalsDf


# ---------------------------------------------------------------------------
# XLS/XLSX proposal loading
# ---------------------------------------------------------------------------
#
# Workbook loading selects the proposal sheet, normalizes legacy/current
# layouts into the shared schema, and retains original Excel row numbers for
# diagnostics.


def load_proposal(code: str, xlsxpath_override: Any = pd.NA) -> dict[str, Any]:
    """Load one proposal workbook into the normalized change-row schema."""
    print(f"# LOAD_PROPOSAL({code},{text(xlsxpath_override)})")
    xlsxpath = _proposal_value(code, "xlsxpath") if is_na(xlsxpath_override) else xlsxpath_override
    try:
        engine = "xlrd" if Path(str(xlsxpath)).suffix.lower() == ".xls" else "openpyxl"
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            excel = pd.ExcelFile(xlsxpath, engine=engine)
        sheet_names = excel.sheet_names
    except Exception as exc:
        log_error(code, 0, "OPEN_XLSX", "", "", "ERROR", "XLSX.OPEN_FAILED",
                  "XLS file could not be opened", str(exc), actionOrder)
        return {"proposalDf": None}
    expected_2022 = {"Proposals Template", "Menu Items (Do not change)"}
    expected_2023 = {"Instructions", "Proposal Template", "Menu Items (Do not change)"}
    expected = expected_2023 if expected_2023.issubset(sheet_names) else expected_2022
    extra = [name for name in sheet_names if name not in expected]
    if extra and is_na(xlsxpath_override):
        log_error(code, 0, "OPEN_XLSX", "", "", "WARNING", "XLSX_EXTRA_SHEETS",
                  "XLS file has additional sheets not present in template",
                  "Worksheet(s) named '" + "','".join(extra) + "' were added", actionOrder)
    proposal_sheets = [name for name in sheet_names if re.search(r"proposal*", name, re.I)]
    if not proposal_sheets:
        if is_na(xlsxpath_override):
            log_error(code, 0, "OPEN_XLSX", "", "", "ERROR", "XLSX_NOT_PROPOSAL",
                      "XLS no acceptably named 'proposal' sheet",
                      "A worksheet name matching 'proposal*' was not found, only: '" + "','".join(sheet_names) + "'", actionOrder)
        return {"proposalDf": None}
    if len(proposal_sheets) > 1:
        if is_na(xlsxpath_override):
            log_error(code, 0, "OPEN_XLSX", "", "", "ERROR", "XLSX_MULTI_PROPOSAL",
                      "XLS file has more than one sheet named 'proposal'",
                      "More than 1 worksheet named 'proposal*' were found: '" + "','".join(proposal_sheets) + "'", actionOrder)
        return {"proposalDf": None}
    return {"proposalDf": read_excel_raw(xlsxpath, proposal_sheets[0])}


# ---------------------------------------------------------------------------
# DOCX proposal metadata loading
# ---------------------------------------------------------------------------
#
# DOCX files are read directly from their XML payload. Required title, author,
# corresponding-author, and abstract sections are reported through the common
# QC logger.


def _docx_text(path: str | Path) -> list[str]:
    namespace = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
    with zipfile.ZipFile(path) as archive:
        root = ElementTree.fromstring(archive.read("word/document.xml"))
    lines: list[str] = []
    for paragraph in root.iter(namespace + "p"):
        pieces: list[str] = []
        for node in paragraph.iter():
            if node.tag == namespace + "t" and node.text:
                pieces.append(node.text)
            elif node.tag == namespace + "tab":
                pieces.append("\t")
            elif node.tag == namespace + "br":
                pieces.append("\n")
        value = "".join(pieces)
        if value:
            lines.append(value)
    return lines


def load_proposal_docx(code: str) -> pd.DataFrame:
    print(f"# LOAD_PROPOSAL_DOCX({code})")
    meta = {"code": code}
    path = _proposal_value(code, "docxpath")
    if is_na(path):
        if params.processing_mode in {"draft", "final"}:
            log_error(code, 0, "OPEN_DOCX", "", "", "ERROR", "DOCX_MISSING",
                      ".docx file not available", "", actionOrder)
        return pd.DataFrame([meta])
    extension = Path(str(path)).suffix.lower().lstrip(".")
    if extension != "docx":
        log_error(code, 0, "OPEN_DOCX", "", "", "WARNING", "DOCX_BAD_EXT",
                  f"File type .{extension} not (yet) supported. Only .docx",
                  f"DOCX_FILENAME={_proposal_value(code, 'docx')}", actionOrder)
        return pd.DataFrame([meta])
    try:
        lines = _docx_text(path)
    except Exception as exc:
        log_error(code, 0, "OPEN_DOCX", "", "", "ERROR", "DOCX_OPEN_FAILED",
                  ".docx file could not be read", str(exc), actionOrder)
        return pd.DataFrame([meta])
    if not any(re.search("Part 1: TITLE, AUTHORS, APPROVALS, etc", line, re.I) for line in lines):
        log_error(code, 0, "OPEN_DOCX", "", "", "INFO", "DOCX_UNKNOWN_TEMPLATE",
                  "Unrecognized version of DOCX proposal template. Not (yet) supported.",
                  f"DOCX_FILENAME={_proposal_value(code, 'docx')}", actionOrder)
        return pd.DataFrame([meta])

    def after(pattern: str, key: str, missing_code: str, missing_message: str,
              same_line: bool = False) -> None:
        indices = [i for i, line in enumerate(lines) if re.search(pattern, line, re.I)]
        if indices:
            index = indices[0]
            if same_line:
                value = re.split(pattern, lines[index], maxsplit=1, flags=re.I)[-1]
            else:
                value = lines[index + 1] if index + 1 < len(lines) else ""
            meta[key] = re.sub(r"\t+", " ", re.sub(r"\n+", "; ", value)).strip()
        else:
            log_error(code, 0, "OPEN_DOCX", "", "", "WARNING", missing_code,
                      missing_message, f"Line containing expression '{pattern}' not found", actionOrder)
    after(r"^Short title\s*:\s*", "title", "DOCX_TITLE_MISSING", "Title not found", True)
    after(r"^Author.* and email address.*", "authorsEmails", "DOCX_AUTHORS_MISSING", "Authors list not found")
    after(r"^Corresponding author", "correspondingAuthor", "DOCX_CORR_AUTHOR_MISSING", "Corresponding Author not found")
    after(r"^Abstract", "abstract", "DOCX_ABSTRACT_MISSING", "Abstract not found")
    return pd.DataFrame([meta])



# ---------------------------------------------------------------------------
# Proposal row normalization and QC
# ---------------------------------------------------------------------------
#
# qc_proposal validates workbook identity, controlled vocabulary values,
# taxonomy fields, naming conventions, descriptives, and action prerequisites.
# It also derives the dot-prefixed fields consumed by change application.


def _normalized_header(value: Any) -> str:
    return re.sub(r"[^A-Za-z0-9]", "", text(value, "")).lower()


def _excel_column(number: int) -> str:
    result = ""
    value = number + 1
    while value:
        value, remainder = divmod(value - 1, 26)
        result = chr(65 + remainder) + result
    return result


def _log_change_rows(code: str, rows: pd.DataFrame, level: str, error: str,
                     message: str, note_builder: Any = "") -> None:
    if rows.empty:
        return
    notes = (rows.apply(note_builder, axis=1).tolist()
             if callable(note_builder) else note_builder)
    log_error(code, rows[".linenum"], rows.get("change", pd.NA),
              rows.get("rank", pd.NA), rows.get(".changeTaxon", pd.NA),
              level, error, message, notes, actionOrder)


def _derive_change_taxonomy(change: pd.DataFrame) -> None:
    """Refresh taxonomy helper columns after proposal-cell normalization."""
    change[".srcTaxon"] = change[xlsx_change_src_colnames].apply(taxon_get_name, axis=1)
    change[".srcRank"] = change[xlsx_change_src_colnames].apply(taxon_get_rank, axis=1)
    change[".srcLineage"] = change[xlsx_change_src_colnames].apply(taxon_get_lineage, axis=1)
    change[".srcParentName"] = change[xlsx_change_src_colnames].apply(taxon_get_parent_name, axis=1)
    change[".destTaxon"] = change[xlsx_change_dest_colnames].apply(taxon_get_name, axis=1)
    change[".destRank"] = change[xlsx_change_dest_colnames].apply(taxon_get_rank, axis=1)
    change[".destLineage"] = change[xlsx_change_dest_colnames].apply(taxon_get_lineage, axis=1)
    change[".destParentName"] = change[xlsx_change_dest_colnames].apply(taxon_get_parent_name, axis=1)
    change[".destParentLineage"] = change[xlsx_change_dest_colnames].apply(taxon_get_parent_lineage, axis=1)
    change[".changeTaxon"] = change[".srcTaxon"].combine_first(change[".destTaxon"])


def qc_proposal(code: str, proposalDf: pd.DataFrame) -> dict[str, Any]:
    """Validate and normalize every change row from one proposal workbook."""
    global proposalsDf
    template_version = "error"
    version_cell = proposalDf.iloc[1, 0] if len(proposalDf) > 1 else pd.NA
    if not is_na(version_cell) and re.match(r"version 202[34]\.", str(version_cell), re.I):
        full_version = str(version_cell)[8:].strip()
        if full_version.startswith("2024.") and full_version != "2024.1":
            proposalsDf.loc[code, "templateVersion"] = "error"
            log_error(code, 2, "OPEN_XLSX", "", "", "ERROR", "XLSX.TEMPLATE_UNK",
                      "XLSX template version", f"Cell A2='{version_cell}'", actionOrder)
            return {}
        template_version = "2023."
        # Compare only expected populated headings; empty merged cells are
        # represented slightly differently by Excel engines.
        for row_number, expected in [(3, XLSX_2023_ROW3), (4, XLSX_2023_ROW4)]:
            observed = proposalDf.iloc[row_number - 1, :len(expected)].tolist()
            mismatches = [i for i, (actual, wanted) in enumerate(zip(observed, expected))
                          if (not is_na(wanted) and _normalized_header(actual) != _normalized_header(wanted))]
            if mismatches:
                details = "; ".join(
                    f"column {_excel_column(i)}='{text(observed[i])}' instead of '{text(expected[i])}'"
                    for i in mismatches
                )
                log_error(code, 2, "OPEN_XLSX", "", "", "ERROR",
                          f"XLSX.TEMPLATE_LINE{row_number}", "XLSX template version",
                          f"row{row_number} mismatches template at: {details}", actionOrder)
                proposalsDf.loc[code, "templateVersion"] = "error"
                return {}
    else:
        # Locate the old template header rather than depending on pandas/readxl
        # differences in trailing blank cells.
        header_row = None
        for idx in range(min(8, len(proposalDf) - 1)):
            values_norm = [_normalized_header(v) for v in proposalDf.iloc[idx].tolist()]
            if "currenttaxonomy" in values_norm and "proposedtaxonomy" in values_norm:
                header_row = idx + 1
                break
        if header_row is None:
            log_error(code, 3, "OPEN_XLSX", "", "", "ERROR", "XLSX.TEMPLATE_UNK",
                      "XLSX template version", "ROW2/ROW3 is unrecognized", actionOrder)
            return {}
        second = proposalDf.iloc[header_row].tolist()
        template_version = "v2" if _normalized_header(second[0]).startswith("code") else "v1"
    proposalsDf.loc[code, "templateVersion"] = template_version
    if template_version == "v1":
        log_error(code, 3, "OPEN_XLSX", "", "", "INFO", "XLSX.OLD_TEMPLATE_V1",
                  "XLSX template version",
                  f"You are using version {template_version}. Please get the latest version from {params.template_url}",
                  actionOrder)

    code_value, code_cell, code_row = "missing", "undefined", pd.NA
    if template_version == "v1":
        code_value, code_cell, code_row = proposalDf.iloc[0, 0], "A1", 1
    elif template_version == "v2":
        code_value, code_cell, code_row = proposalDf.iloc[2, 0], "A3", 3
    elif template_version == "2023.":
        code_value, code_cell, code_row = proposalDf.iloc[0, 4], "E1", 1
        if "E1" in proposalDf.attrs.get("formula_cells", {}):
            # E1 is generated from this already-validated filename. Ignore a
            # stale cached formula value such as #VALUE! or 2025.000B.
            code_value = code
    if text(code_value) != code:
        if text(code_value).startswith("Code"):
            if params["show.xlsx.code_miss"]:
                log_error(code, code_row, "OPEN_XLSX", "", "", "INFO", "XLSX.CODE_MISS",
                          "XLSX code missing",
                          f"XLSX cell {code_cell} is '{text(code_value)}'; replace with the actual code: '{code}'",
                          actionOrder)
        else:
            log_error(code, code_row, "OPEN_XLSX", "", "", "ERROR", "XLSX.CODE_BAD",
                      "XLSX code wrong",
                      f"XLSX cell {code_cell} does not match proposal code from filename: '{text(code_value)}' should be '{code}' ",
                      actionOrder)
    if params.verbose:
        print(f"      {code} XLSX template {template_version}")

    if template_version == "v1":
        first_data_row, selected = 4, XLSX_V1_CHANGE_COLS
    elif template_version == "v2":
        first_data_row, selected = 4, XLSX_V2_CHANGE_COLS
    else:
        first_data_row, selected = 5, XLSX_2023_CHANGE_COLS
    if len(proposalDf) < first_data_row:
        log_error(code, first_data_row, "OPEN_XLSX", "", "", "ERROR", "XLSX.EMPTY",
                  "XLSX no change rows found", "", actionOrder)
        return {}
    change = proposalDf.iloc[first_data_row - 1:, selected].copy()
    change.columns = xlsx_change_colnames
    change[".noErrors"] = True
    change[".errors"] = pd.NA
    change[".code"] = code
    change[".linenum"] = change.index + 1
    has_data = change[xlsx_change_src_colnames + xlsx_change_dest_colnames].notna().any(axis=1)
    change = change.loc[has_data].copy()
    if change.empty:
        log_error(code, first_data_row, "OPEN_XLSX", "", "", "ERROR", "XLSX.EMPTY",
                  "XLSX no change rows found", "", actionOrder)
        return {}
    proposalsDf.loc[code, "nChanges"] = len(change)

    _derive_change_taxonomy(change)

    data_columns = [column for column in change if not column.startswith(".")]
    for column in data_columns:
        series = change[column]
        cleaned = series.map(lambda value: value if is_na(value) else str(value))
        replacements = [
            (r"0xCA|\u00a0", " "), (r"[\u2012\u2013\u2014\u2015]+", "-"),
            (r"[\u201C\u201D\u201E\u201F]+", '"'),
            (r"[\u2018\u2019\u201A\u201B]+", "'"), (r"[\r\n]+", ";"),
            (r"^[ \t]+", ""), (r"[ \t]+$", ""), (r" {2,}", " "),
        ]
        for pattern, replacement in replacements:
            cleaned = cleaned.map(
                lambda value: value if is_na(value) else re.sub(pattern, replacement, str(value))
            )
        if column != "comments":
            quote_mask = cleaned.map(lambda value: False if is_na(value) else bool(re.search(r'"+', str(value))))
            if quote_mask.any():
                rows = change.loc[quote_mask]
                _log_change_rows(code, rows, "INFO", "XLSX.QUOTES_REMOVED",
                                 "XLSX has quote",
                                 lambda row: f"{column}:{text(row[column])} (replacing with '')")
                cleaned = cleaned.map(lambda value: value if is_na(value) else re.sub(r'"+', "", str(value)))
        change[column] = cleaned

    for validation in VALUE_VALIDATION:
        column = validation["col"]
        pattern = normalize_posix_regex(validation["regex"])
        if column not in change:
            continue
        if validation["type"] == "replace":
            unicode_printable = "[:print:]" in validation["regex"]
            matches = change[column].map(
                lambda value: False if is_na(value) else (
                    any(not char.isprintable() for char in str(value))
                    if unicode_printable else bool(re.search(pattern, str(value)))
                )
            )
            if matches.any():
                rows = change.loc[matches]
                _log_change_rows(code, rows, validation["class"], validation["code"],
                                 f"XLSX has {validation['warn']}",
                                 lambda row: f"{column} : {text(row[column])}")
                if unicode_printable:
                    replacement = validation["replace"]
                    change[column] = change[column].map(
                        lambda value: value if is_na(value) else "".join(
                            char if char.isprintable() else replacement
                            for char in str(value)
                        )
                    )
                else:
                    change[column] = change[column].map(
                        lambda value: value if is_na(value) else re.sub(
                            pattern, validation["replace"], str(value)
                        )
                    )
        else:
            failed = change[column].map(
                lambda value: False if is_na(value) else re.fullmatch(pattern, str(value)) is None
            )
            if failed.any():
                _log_change_rows(code, change.loc[failed], validation["class"],
                                 validation["code"], validation["warn"],
                                 lambda row: f"{column}:{text(row[column])}")

    # Convert any remaining please-select variants to missing.
    for column in data_columns:
        change[column] = change[column].map(
            lambda value: pd.NA if (not is_na(value) and re.search(r"please.*select", str(value), re.I)) else value
        )

    # Controlled vocabulary corrections and errors.
    for cv_name, allowed in cvList.items():
        if cv_name.startswith(".") or cv_name not in change:
            continue
        allowed_values = [str(value) for value in allowed if not is_na(value)]
        exact = set(allowed_values)
        close_map: dict[str, list[str]] = {}
        for value in allowed_values:
            close_map.setdefault(re.sub(r"[^A-Za-z;+-]", "", value).lower(), []).append(value)
        bad_indices: list[Any] = []
        for idx, value in change[cv_name].items():
            if is_na(value):
                # NA is allowed for descriptives; required checks are below.
                continue
            if str(value) in exact:
                continue
            key = re.sub(r"[^A-Za-z;+-]", "", str(value)).lower()
            candidates = close_map.get(key, [])
            if len(candidates) == 1:
                if params.tmi:
                    _log_change_rows(code, change.loc[[idx]], "TMI", "XLSX.TYPO_TERM",
                                     f"fixed term with typo (space,caps) in column {cv_name}",
                                     f"Term '{value}' replaced with '{candidates[0]}'.")
                change.at[idx, cv_name] = candidates[0]
            else:
                bad_indices.append(idx)
        if bad_indices:
            bad_rows = change.loc[bad_indices]
            change.loc[bad_indices, ".errors"] = [
                f"CV '{cv_name}' incorrect value [{text(v)}]. Valid terms: [{','.join(allowed_values)}]"
                for v in bad_rows[cv_name]
            ]
            change.loc[bad_indices, ".noErrors"] = False
            _log_change_rows(code, bad_rows, "ERROR", "XLSX.INVALID_TERM",
                             f"XLSX incorrect term in column {cv_name}",
                             lambda row: f"XLSX incorrect value [{text(row[cv_name])}]. Valid terms: [{','.join(allowed_values)}]")

    # Taxon helper values were initially needed for cleaning diagnostics.
    # Recompute them from the cleaned cells so invisible surrounding spaces do
    # not cause false source, parent, binomial, suffix, or case errors.
    _derive_change_taxonomy(change)

    change[".action"] = change["change"].map(
        lambda value: ACTION_MAP.get(str(value).lower()) if not is_na(value) else pd.NA
    )
    bad_action = change[".action"].isna()
    if bad_action.any():
        change.loc[bad_action, ".noErrors"] = False
        _log_change_rows(code, change.loc[bad_action], "ERROR", "ACTION.UNK",
                         "XLSX incorrect term in column Change",
                         f"Valid terms: [{','.join(ACTION_MAP)}]")

    # Genome Coverage is required only for an explicitly new species. Existing
    # taxa moved, renamed, merged, split, promoted, or demoted retain it.
    species_rows = change[".destRank"].eq("species") & change[".action"].eq("new")
    missing_coverage = species_rows & (
        change["genomeCoverage"].isna() | change["genomeCoverage"].astype("string").str.strip().eq("")
    )
    if missing_coverage.any():
        change.loc[missing_coverage, ".noErrors"] = False
        _log_change_rows(code, change.loc[missing_coverage], "ERROR",
                         "XLSX.GENOME_COVERAGE_MISSING",
                         "Proposed species must have a Genome Coverage value",
                         "Genome Coverage is empty")

    # Rank suffix and capitalization QC.
    rank_cv = dbCvList.get("rank", pd.DataFrame())
    suffix_columns = [column for column in rank_cv if "suffix" in column]
    all_suffixes: list[tuple[str, str]] = []
    for _, rank_row in rank_cv.iterrows():
        for column in suffix_columns:
            suffix = rank_row.get(column)
            if not is_na(suffix) and str(suffix):
                all_suffixes.append((str(rank_row.get("name")), str(suffix)))
    for idx, row in change.iterrows():
        row_action = row.get(".action")
        if (not is_na(row_action) and row_action == "abolish") or is_na(row.get(".destTaxon")):
            continue
        dest_name, dest_rank = str(row[".destTaxon"]), str(row[".destRank"])
        matching_rank = rank_cv.loc[rank_cv["name"].astype(str) == dest_rank]
        valid_suffixes = [str(matching_rank.iloc[0][column]) for column in suffix_columns
                          if not matching_rank.empty and not is_na(matching_rank.iloc[0][column])]
        if dest_rank != "species" and valid_suffixes and not any(dest_name.lower().endswith(s.lower()) for s in valid_suffixes):
            _log_change_rows(code, change.loc[[idx]], "ERROR", "MISSING_RANK_SUFFIX",
                             "Taxon name does not end with a valid suffix for it's rank",
                             f"New name: {dest_name}; Valid suffixes for rank [{dest_rank}]: {','.join(valid_suffixes)}")
        if dest_rank != "species" and not valid_suffixes:
            wrong = next(((rank, suffix) for rank, suffix in all_suffixes
                          if dest_name.lower().endswith(suffix.lower())), None)
            if wrong:
                _log_change_rows(code, change.loc[[idx]], "ERROR", "SUFFIX_RANK_MISMATCH",
                                 "Taxon name does not end with a valid suffix for it's rank",
                                 f"New name: '{dest_name}' has the suffix of a '{wrong[0]}' but is a '{dest_rank}'")
        if dest_rank != "species" and not re.fullmatch(r"[A-Z][a-z0-9-]+", dest_name):
            _log_change_rows(code, change.loc[[idx]], "ERROR", "TAXON_NAME_CASE",
                             "Non-species taxa must be Capitalized.",
                             f"{dest_name}: Non-species taxa must start with a capital, followed by one or more lower-case letters, numbers or hyphens.")
    return {"changeDf": change}


# ---------------------------------------------------------------------------
# Proposal-set loading and QC orchestration
# ---------------------------------------------------------------------------
#
# Documents and workbooks are validated proposal by proposal. Results remain
# separated in changeList until cross-proposal dependency ordering is possible.


def load_and_qc_proposals(proposals: pd.DataFrame,
                          existing: dict[str, pd.DataFrame]) -> dict[str, pd.DataFrame]:
    global proposalsDf
    result = dict(existing)
    codes = list(proposals.index)
    for position, code in enumerate(codes, 1):
        if code in result:
            continue
        if is_na(_proposal_value(code, "xlsx")):
            print(f"# SKIPPED: {code} (no .xlsx)")
            continue
        before = len(allErrorDf)
        meta = load_proposal_docx(code)
        for column in meta:
            proposalsDf.loc[code, column] = first(meta[column])
        print(f"# LOADED: {code} DOCX with {len(allErrorDf) - before} errors/warnings")
        before = len(allErrorDf)
        raw = load_proposal(code).get("proposalDf")
        print(f"# LOADED: {code} XLS with {len(allErrorDf) - before} errors/warnings")
        if raw is None:
            print(f"SKIP: {code}: proposal could not be loaded")
            continue
        print(f"# QC start: {code}: proposal loaded ({position} out of {len(codes)})")
        checked = qc_proposal(code, raw)
        print(f"# QCed:      {code} with {len(allErrorDf) - before} errors/warnings")
        if "changeDf" in checked:
            result[code] = checked["changeDf"]
        write_error_summary(allErrorDf)
    print("changeList: " + " ".join(result))
    write_error_summary(allErrorDf)
    return result


# ---------------------------------------------------------------------------
# Cross-proposal change dependency ordering
# ---------------------------------------------------------------------------
#
# A stable topological order lets proposals refer to taxa created or changed by
# other rows in the same run. Descendants move before branch removal, while
# parent creation precedes actions that target the new parent.

proposedAccessions: dict[str, list[dict[str, Any]]] = {}


def merge_and_order_changes(changes: Mapping[str, pd.DataFrame]) -> pd.DataFrame:
    """Merge proposal changes and return a stable dependency-aware order."""
    """Combine changes and order both R-style ranks and proposal dependencies."""
    frames = [frame.copy() for frame in changes.values() if frame is not None and not frame.empty]
    if not frames:
        return pd.DataFrame()
    merged = pd.concat(frames, ignore_index=True, sort=False)
    merged[".codeLine"] = merged[".code"].map(text) + ":" + merged[".linenum"].map(text)

    # Mirror the R tentative ordering: destructive actions bottom-up, then
    # constructive actions top-down, with creates before moves at each rank.
    action_group = {
        "abolish": -1, "merge": -1,
        "new": 1, "promote": 1, "demote": 1, "split": 1,
        "rename": 1, "move": 1,
    }
    sub_action = {
        "abolish": 0, "merge": 0,
        "new": 1, "promote": 1, "demote": 1, "split": 1,
        "rename": 2, "move": 2,
    }
    rank_order = {rank.lower(): index for index, rank in enumerate(xlsx_change_ranks)}
    effective_rank = merged[".destRank"].combine_first(merged[".srcRank"])
    merged[".actionGroup"] = merged[".action"].map(action_group).fillna(99)
    raw_rank = effective_rank.astype("string").str.lower().map(rank_order).fillna(99)
    merged[".rankPriority"] = raw_rank * merged[".actionGroup"].map(
        lambda value: -1 if value == -1 else 1
    )
    merged[".subActionPriority"] = merged[".action"].map(sub_action).fillna(99)
    merged = merged.sort_values(
        [".actionGroup", ".rankPriority", ".subActionPriority", ".code", ".linenum"],
        kind="stable", na_position="last",
    ).reset_index(drop=True)

    # Build a dependency graph on top of that stable tentative order.
    count = len(merged)
    outgoing: dict[int, set[int]] = {i: set() for i in range(count)}
    incoming = [0] * count

    def add_edge(before: int, after: int) -> None:
        if before == after or after in outgoing[before]:
            return
        outgoing[before].add(after)
        incoming[after] += 1

    destinations: dict[str, list[int]] = {}
    source_users: dict[str, list[int]] = {}
    parent_users: dict[str, list[int]] = {}
    for index, row in merged.iterrows():
        destination = row.get(".destTaxon")
        source = row.get(".srcTaxon")
        parent = row.get(".destParentName")
        if not is_na(destination):
            destinations.setdefault(str(destination), []).append(index)
        if not is_na(source):
            source_users.setdefault(str(source), []).append(index)
        if not is_na(parent):
            parent_users.setdefault(str(parent), []).append(index)

    # A taxon created/renamed in the proposal must exist before it is used as
    # either a source or a destination parent.
    for name, creators in destinations.items():
        for before in creators:
            for after in source_users.get(name, []) + parent_users.get(name, []):
                add_edge(before, after)

    def lineage_parts(value: Any) -> list[str]:
        if is_na(value):
            return []
        return [part for part in str(value).split(";") if part]

    # Before removing a taxon, first apply proposal rows that move or remove
    # its descendants. This lets promotion/reorganization proposals empty the
    # old branch before their final Abolish/Merge rows.
    removal_indices = merged.index[merged[".action"].isin(["abolish", "merge"])]
    for removal in removal_indices:
        removed_name = merged.at[removal, ".srcTaxon"]
        if is_na(removed_name):
            continue
        removed_name = str(removed_name)
        for before, row in merged.iterrows():
            if before == removal:
                continue
            source_parts = lineage_parts(row.get(".srcLineage"))
            if removed_name not in source_parts[:-1]:
                continue
            destination_parts = lineage_parts(row.get(".destLineage"))
            if removed_name not in destination_parts:
                add_edge(before, removal)

    available = [i for i, degree in enumerate(incoming) if degree == 0]
    ordered: list[int] = []
    while available:
        node = available.pop(0)
        ordered.append(node)
        for target in sorted(outgoing[node]):
            incoming[target] -= 1
            if incoming[target] == 0:
                available.append(target)
                available.sort()
    if len(ordered) != count:
        cyclic = [i for i in range(count) if i not in ordered]
        for index in cyclic:
            row = merged.iloc[index]
            log_change_error(row, "ERROR", "CHAIN.HIDDEN",
                             "Circular or hidden proposal dependency",
                             f"Could not order {row.get('.srcTaxon')} -> {row.get('.destTaxon')}")
            ordered.append(index)
    merged[".originalOrder"] = range(1, count + 1)
    merged = merged.iloc[ordered].reset_index(drop=True)
    merged[".chainOrder"] = range(1, count + 1)
    merged[".changeOrder"] = merged[".chainOrder"]
    return merged.drop(columns=[
        ".actionGroup", ".rankPriority", ".subActionPriority",
    ])


# ---------------------------------------------------------------------------
# Working-taxonomy lookup and lineage maintenance
# ---------------------------------------------------------------------------
#
# These helpers resolve unique taxa and recursively propagate lineage changes
# through descendants after moves, renames, promotions, or demotions.


def _matching_taxa(frame: pd.DataFrame, name: Any) -> pd.Index:
    if is_na(name) or "name" not in frame:
        return pd.Index([])
    return frame.index[frame["name"].astype("string").eq(str(name)).fillna(False)]


def _resolve_single(frame: pd.DataFrame, name: Any, change: pd.Series,
                    role: str) -> Any:
    matches = _matching_taxa(frame, name)
    if len(matches) == 1:
        return matches[0]
    if len(matches) == 0:
        log_change_error(change, "ERROR", f"{role}.NOT_FOUND",
                         f"{role.title()} taxon not found in current working MSL",
                         text(name))
    else:
        log_change_error(change, "ERROR", f"{role}.MULTIPLE",
                         f"{role.title()} taxon is not unique in current working MSL",
                         f"{text(name)} matched {len(matches)} rows")
    return None


def _lineage_from_parent(parent_index: Any, name: str) -> str:
    if parent_index is None:
        return name
    if text(newMSL.at[parent_index, "rank"], "") == "tree":
        return name
    parent_lineage = newMSL.at[parent_index, "lineage"]
    if is_na(parent_lineage) or not str(parent_lineage):
        return name
    return f"{parent_lineage};{name}"


def update_lineage(parent_id: Any, parent_lineage: Any,
                   parent_otherLineage: Any = pd.NA,
                   parent_otherLineageProposal: Any = pd.NA,
                   parent_otherLineageAction: Any = pd.NA) -> int:
    """Recursively update descendant lineages and R-style provenance."""
    parent_rows = newMSL.index[newMSL["taxnode_id"].eq(parent_id)]
    if len(parent_rows):
        parent_row = parent_rows[0]
        if ".otherLineage" in newMSL:
            newMSL.at[parent_row, ".otherLineage"] = parent_otherLineage
        for column, value in [
            (".otherLineageProposal", parent_otherLineageProposal),
            (".otherLineageAction", parent_otherLineageAction),
        ]:
            if column not in newMSL:
                continue
            previous = newMSL.at[parent_row, column]
            newMSL.at[parent_row, column] = (
                value if is_na(previous) else f"{previous};{value}"
            )

    count = 0
    children = newMSL.index[newMSL["parent_id"].eq(parent_id)]
    for child in children:
        child_lineage = f"{parent_lineage};{newMSL.at[child, 'name']}"
        newMSL.at[child, "lineage"] = child_lineage
        if ".otherLineage" in newMSL and not is_na(parent_otherLineage):
            newMSL.at[child, ".otherLineage"] = (
                f"{parent_otherLineage};{newMSL.at[child, 'name']}"
            )
        child_other_lineage = (
            newMSL.at[child, ".otherLineage"]
            if ".otherLineage" in newMSL else pd.NA
        )
        count += 1 + update_lineage(
            newMSL.at[child, "taxnode_id"],
            child_lineage,
            child_other_lineage,
            parent_otherLineageProposal,
            parent_otherLineageAction,
        )
    return count


# ---------------------------------------------------------------------------
# Cache compatibility stubs
# ---------------------------------------------------------------------------
#
# Cache flags remain accepted for shell compatibility, but the first Python
# implementation deliberately performs uncached loads. These boundaries keep a
# later cache implementation isolated from validation logic.


def load_reference_cache() -> bool:
    """Compatibility stub: reference caching is deferred in version one."""
    return False


def save_reference_cache() -> bool:
    """Compatibility stub: reference caching is deferred in version one."""
    return False


# ---------------------------------------------------------------------------
# Taxonomy mutation metadata helpers
# ---------------------------------------------------------------------------
#
# Mutation helpers update descriptives and the in_*/out_* audit columns while
# preserving links to the previous MSL taxnode and proposal source.


def _refresh_descendant_lineages(parent_taxnode_id: Any) -> None:
    children = newMSL.index[newMSL["parent_id"].eq(parent_taxnode_id)]
    for child in children:
        parent_matches = newMSL.index[newMSL["taxnode_id"].eq(newMSL.at[child, "parent_id"])]
        parent = parent_matches[0] if len(parent_matches) else None
        newMSL.at[child, "lineage"] = _lineage_from_parent(parent, str(newMSL.at[child, "name"]))
        _refresh_descendant_lineages(newMSL.at[child, "taxnode_id"])


def _set_descriptives(index: Any, change: pd.Series) -> None:
    species_only = {
        "exemplarAccession", "exemplarName", "Abbrev", "exemplarIsolate",
        "genomeCoverage",
    }
    for source, destination in XLSX2DB_MAP.items():
        if destination not in newMSL or source not in change:
            continue
        # The R implementation only persists exemplar/isolate fields on
        # species. Templates may repeat the accession on a preceding genus row.
        if source in species_only and change.get(".destRank") != "species":
            continue
        value = change.get(source)
        if source in {"molecule", "hostSource", "rank"} and not is_na(value):
            value = dbCvMapList.get(source, {}).get(str(value), value)
        if not is_na(value) and str(value) != "":
            newMSL.at[index, destination] = value


def _proposal_zip(change: pd.Series) -> Any:
    basename = _proposal_value(change.get(".code"), "basename", pd.NA)
    return pd.NA if is_na(basename) else f"{basename}.zip"


def _set_in_fields(index: Any, change: pd.Series, action: str) -> None:
    newMSL.at[index, "in_change"] = action
    newMSL.at[index, "in_target"] = change.get(".destLineage")
    newMSL.at[index, "in_filename"] = _proposal_zip(change)
    newMSL.at[index, "in_notes"] = f"xlsx_row={text(change.get('.linenum'), '')}"
    newMSL.at[index, "prev_proposals"] = change.get(".code")


def _append_in_fields(index: Any, change: pd.Series) -> None:
    filename = _proposal_zip(change)
    old_filename = newMSL.at[index, "in_filename"]
    if not is_na(filename):
        newMSL.at[index, "in_filename"] = (
            filename if is_na(old_filename) or not str(old_filename)
            else f"{old_filename};{filename}"
        )
    old_notes = newMSL.at[index, "in_notes"]
    note = f"linenum={text(change.get('.linenum'), '')}"
    newMSL.at[index, "in_notes"] = (
        note if is_na(old_notes) or not str(old_notes) else f"{old_notes};{note}"
    )
    newMSL.at[index, "in_target"] = change.get(".destLineage")


def _set_out_fields(index: Any, change: pd.Series, action: str,
                    target: Any) -> None:
    filename = _proposal_zip(change)
    previous = newMSL.at[index, ".prev_taxnode_id"] if ".prev_taxnode_id" in newMSL else pd.NA
    matches = curMSL.index[curMSL["taxnode_id"].eq(previous)].tolist() if not is_na(previous) else []
    for match in matches:
        old_filename = curMSL.at[match, "out_filename"]
        curMSL.at[match, "out_change"] = action
        curMSL.at[match, "out_target"] = target
        curMSL.at[match, "out_filename"] = (
            filename if is_na(old_filename) or not str(old_filename)
            else f"{old_filename};{filename}"
        )
        curMSL.at[match, "out_notes"] = f"linenum={text(change.get('.linenum'), '')}"
        if "out_updated" in curMSL:
            curMSL.at[match, "out_updated"] = True

# ---------------------------------------------------------------------------
# Per-change validation and proposed-accession tracking
# ---------------------------------------------------------------------------
#
# Accession checks combine species_isolates, historical abolished accessions,
# and accessions already accepted earlier in this proposal run. Results are
# aggregated into one R-style diagnostic per change row.


def _check_accessions(change: pd.Series,
                      exclude_taxnode_ids: Sequence[Any] = ()) -> bool:
    """Validate every accession and report one R-style result per change row."""
    raw = change.get("exemplarAccession")
    accessions = split_accessions(raw)
    if not accessions:
        return True
    okay = True
    if is_pending_accession(raw):
        log_change_error(change, "WARNING", "ACCESSION.PENDING",
                         "Exemplar accession is pending", text(raw))
        return True
    abolished = find_abolished_accession_matches(raw)
    if not abolished.empty:
        okay = False
        log_change_error(
            change, "ERROR", "ACCESSION.ABOLISHED",
            "Exemplar accession is abolished",
            f"{';'.join(accessions)}: {format_abolished_accession_matches(abolished)}",
        )

    duplicate_owners: list[str] = []
    for accession in accessions:
        duplicate_owners.extend(
            find_vmr_accession_matches(accession, exclude_taxnode_ids)
        )
        prior = proposedAccessions.get(accession, [])
        current_taxon = text(change.get(".destTaxon"), "")
        conflicting = [item for item in prior if item["taxon"] != current_taxon]
        duplicate_owners.extend(
            f"{item['taxon']} ({item['codeLine']})" for item in conflicting
        )
    if duplicate_owners:
        okay = False
        owners = list(dict.fromkeys(duplicate_owners))
        action = text(change.get(".action"), "").lower()
        error_code = (
            "CREATE.DUP_ACC" if action in {"new", "split"}
            else "MOVE.DUP_ACC"
        )
        log_change_error(
            change, "ERROR", error_code,
            f"Change={action.upper()}, a species with this accession number already exists",
            f"accession={text(raw)}, existingSpecies={'; '.join(owners)}",
        )
    return okay


def _register_accessions(change: pd.Series) -> None:
    if (change.get(".destRank") != "species"
            or is_pending_accession(change.get("exemplarAccession"))):
        return
    for accession in split_accessions(change.get("exemplarAccession")):
        proposedAccessions.setdefault(accession, []).append({
            "taxon": text(change.get(".destTaxon"), ""),
            "codeLine": text(change.get(".codeLine"), ""),
        })


def _validate_destination(change: pd.Series, allow_index: Any = None) -> bool:
    destination = change.get(".destTaxon")
    matches = _matching_taxa(newMSL, destination)
    if allow_index is not None:
        matches = matches[matches != allow_index]
    if len(matches):
        log_change_error(change, "ERROR", "DEST.EXISTS",
                         "Proposed destination taxon already exists",
                         f"{text(destination)} matched {len(matches)} row(s)")
        return False
    return True


def _validate_species_fields(change: pd.Series, require_coverage: bool = True) -> bool:
    if change.get(".destRank") != "species":
        return True
    okay = True
    name = text(change.get(".destTaxon"), "")
    genus = text(change.get("genus"), "")
    if not re.fullmatch(r"\S+ \S+", name):
        log_change_error(change, "ERROR", "CREATE.SPECIES_BINOMIAL",
                         "Species name must use binomial naming", name)
        okay = False
    elif genus and name.split()[0] != genus:
        log_change_error(change, "ERROR", "CREATE.SPECIES_GENUS",
                         "Species first word must match proposed genus",
                         f"species={name}; genus={genus}")
        okay = False
    if require_coverage and (
            is_na(change.get("genomeCoverage"))
            or not text(change.get("genomeCoverage"), "").strip()):
        # Normally logged during XLSX QC; retained here for direct callers.
        log_change_error(change, "ERROR", "XLSX.GENOME_COVERAGE_MISSING",
                         "Proposed species must have a Genome Coverage value", "Genome Coverage is empty")
        okay = False
    return okay


def _lineage_metadata(index: Any, column: str) -> Any:
    """Read optional R-style lineage provenance from the working MSL."""
    return newMSL.at[index, column] if column in newMSL else pd.NA


def _original_parent_lineage(parent_index: Any) -> Any:
    """Find a working parent's pre-proposal lineage in the reference MSL."""
    previous_id = _lineage_metadata(parent_index, ".prev_taxnode_id")
    if is_na(previous_id) or "taxnode_id" not in curMSL or "lineage" not in curMSL:
        return pd.NA
    matches = curMSL.index[curMSL["taxnode_id"].eq(previous_id)]
    return curMSL.at[matches[0], "lineage"] if len(matches) == 1 else pd.NA


def _log_create_parent_lineage(change: pd.Series, parent_index: Any) -> None:
    """Mirror the R create/split diagnostic for an unexpected parent lineage."""
    expected = change.get(".destParentLineage")
    observed = newMSL.at[parent_index, "lineage"]
    if is_na(expected) or is_na(observed) or str(expected) == str(observed):
        return

    other_lineage = _lineage_metadata(parent_index, ".otherLineage")
    original_lineage = _original_parent_lineage(parent_index)
    previous_proposals = _lineage_metadata(parent_index, "prev_proposals")
    parent_created_during_run = (
        is_na(_lineage_metadata(parent_index, ".prev_taxnode_id"))
        and not is_na(previous_proposals)
    )
    # R records .otherLineage while applying earlier proposals. Until Python's
    # complete recursive provenance model is ported, the reference MSL gives
    # the same classification: if it matched the proposal before this run, the
    # mismatch was introduced by an earlier applied proposal and is INFO.
    parent_changed_during_run = (
        not is_na(original_lineage) and str(original_lineage) == str(expected)
    )
    if (not is_na(other_lineage)
            or parent_changed_during_run
            or parent_created_during_run):
        lineage_proposal = _lineage_metadata(parent_index, ".otherLineageProposal")
        if is_na(lineage_proposal):
            lineage_proposal = previous_proposals
        lineage_action = _lineage_metadata(parent_index, ".otherLineageAction")
        if is_na(lineage_action):
            lineage_action = "earlier proposal change"
        log_change_error(
            change, "INFO", "CREATE.PARENT_RENAMED",
            f"Change={text(change.get('.action'), '').upper()}, proposed parent taxon "
            "exists, but not with expected name/lineage, using observed lineage",
            f"otherProposal(s)={text(lineage_proposal)} did a "
            f"{text(lineage_action)}, "
            f"PROPOSED//OBSERVED={diff_lineages(expected, observed)}",
        )
        return

    log_change_error(
        change, "WARNING", "CREATE.PARENT_LINEAGE",
        f"Change={text(change.get('.action'), '').upper()}, proposed parent taxon "
        "exists, but not with expected lineage, using observed lineage",
        f"PROPOSED//OBSERVED={diff_lineages(expected, observed)}, "
        f"otherProposals={text(previous_proposals)}",
    )


def _log_move_parent_lineage(change: pd.Series, parent_index: Any) -> None:
    """Mirror the R move-family diagnostic for an unexpected parent lineage."""
    expected = change.get(".destParentLineage")
    observed = newMSL.at[parent_index, "lineage"]
    if is_na(expected) or is_na(observed) or str(expected) == str(observed):
        return

    other_lineage = _lineage_metadata(parent_index, ".otherLineage")
    if not is_na(other_lineage):
        log_change_error(
            change, "INFO", "MOVE.PROPOSED_PARENT_LINEAGE_CHANGED",
            "PROPOSED parent lineage already modified",
            f"proposal(s)="
            f"{text(_lineage_metadata(parent_index, '.otherLineageProposal'))} did a "
            f"{text(_lineage_metadata(parent_index, '.otherLineageAction'))}; "
            f"OBSERVED//PROPOSED={diff_lineages(observed, expected)}",
        )
        return

    parent_name = change.get(".destParentName")
    previous = _matching_taxa(curMSL, parent_name)
    previous_lineage = (
        curMSL.at[previous[0], "lineage"]
        if len(previous) == 1 and "lineage" in curMSL else pd.NA
    )
    if not is_na(previous_lineage) and str(previous_lineage) != str(expected):
        log_change_error(
            change, "WARNING", "MOVE.PARENT_LINEAGE",
            "PROPOSED parent taxon exists, but not with expected lineage",
            f",PROPOSED//CUR={diff_lineages(expected, observed)}, "
            f"PROPOSED={text(expected)}, CUR={text(observed)}, "
            f"otherProposals={text(_lineage_metadata(parent_index, 'prev_proposals'))}",
        )
    else:
        log_change_error(
            change, "INFO", "MOVE.PARENT_LINEAGE_UPDATED2",
            "PROPOSED parent taxon exists, but with updated lineage",
            f",PROPOSED//CUR={diff_lineages(expected, observed)}, "
            f"PROPOSED={text(expected)}, CUR={text(observed)}",
        )


def _merge_into_notes(source_index: Any, destination_index: Any) -> str:
    """Build the R MERGE_INTO.OK diagnostic before the source row is removed."""
    source_rank = text(newMSL.at[source_index, "rank"], "")
    source_name = text(newMSL.at[source_index, "name"], "")
    source_lineage = newMSL.at[source_index, "lineage"]
    destination_rank = text(newMSL.at[destination_index, "rank"], "")
    destination_name = text(newMSL.at[destination_index, "name"], "")
    destination_lineage = newMSL.at[destination_index, "lineage"]
    return (
        f"MERGE {source_rank} named '{source_name}' into {destination_rank} named "
        f"'{destination_name}' CUR//PROPOSED="
        f"{diff_lineages(source_lineage, destination_lineage)}"
    )


# ---------------------------------------------------------------------------
# Atomic taxonomy action implementations
# ---------------------------------------------------------------------------
#
# Create, remove, move, and rename operations mutate newMSL only after their
# action-specific prerequisites pass. Promotions, demotions, splits, and some
# merges reuse these primitives through the dispatcher below.


def _create_taxon(change: pd.Series, action: str,
                  exclude_taxnode_ids: Sequence[Any] = ()) -> bool:
    global newMSL
    destination = change.get(".destTaxon")
    if is_na(destination) or is_na(change.get(".destRank")):
        log_change_error(change, "ERROR", "CREATE.MISSING_DEST",
                         "New taxon requires a proposed name and rank", "")
        return False
    if (not _validate_destination(change)
            or not _validate_species_fields(change, require_coverage=action == "new")):
        return False
    if not _check_accessions(change, exclude_taxnode_ids):
        return False
    parent_name = change.get(".destParentName")
    parent_index = _resolve_single(newMSL, parent_name, change, "PARENT")
    if parent_index is None:
        return False
    _log_create_parent_lineage(change, parent_index)
    parent_taxnode_id = newMSL.at[parent_index, "taxnode_id"]
    row = {column: pd.NA for column in newMSL.columns}
    row["taxnode_id"] = int(pd.to_numeric(newMSL["taxnode_id"], errors="coerce").max()) + 1
    row["parent_id"] = parent_taxnode_id
    row["tree_id"] = newMSL.at[parent_index, "tree_id"]
    row["msl_release_num"] = newMSL.at[parent_index, "msl_release_num"]
    row["level_id"] = dbCvMapList.get("rank", {}).get(str(change.get(".destRank")), pd.NA)
    row["rank"] = change.get(".destRank")
    row["name"] = destination
    row["lineage"] = _lineage_from_parent(parent_index, str(destination))
    row["ictv_id"] = row["taxnode_id"]
    row["is_hidden"] = 0
    row[".prev_taxnode_id"] = pd.NA
    parent_other_lineage = _lineage_metadata(parent_index, ".otherLineage")
    if not is_na(parent_other_lineage):
        row[".otherLineage"] = f"{parent_other_lineage};{destination}"
        row[".otherLineageProposal"] = _lineage_metadata(
            parent_index, ".otherLineageProposal"
        )
        row[".otherLineageAction"] = _lineage_metadata(
            parent_index, ".otherLineageAction"
        )
    index = len(newMSL)
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore", message="The behavior of DataFrame concatenation with empty or all-NA entries",
            category=FutureWarning,
        )
        newMSL.loc[index] = pd.Series(row)
    _set_descriptives(index, change)
    _set_in_fields(index, change, "new" if action == "new" else action)
    _register_accessions(change)
    return True


def _remove_taxon(index: Any, change: pd.Series, action: str,
                  target: Any) -> bool:
    global newMSL
    taxnode_id = newMSL.at[index, "taxnode_id"]
    children = newMSL.index[newMSL["parent_id"].eq(taxnode_id)]
    if len(children):
        log_change_error(change, "ERROR", "ACTION.HAS_CHILDREN",
                         f"Cannot {action} a taxon that still has children",
                         ", ".join(newMSL.loc[children, "name"].astype(str)))
        return False
    _set_out_fields(index, change, action, target)
    newMSL = newMSL.drop(index=index).reset_index(drop=True)
    return True


def _move_or_rename_taxon(index: Any, change: pd.Series, action: str) -> bool:
    destination = change.get(".destTaxon")
    destination_rank = change.get(".destRank")
    if is_na(destination):
        destination = newMSL.at[index, "name"]
    if is_na(destination_rank):
        destination_rank = newMSL.at[index, "rank"]
    if str(destination) != str(newMSL.at[index, "name"]):
        if not _validate_destination(change, allow_index=index):
            return False
    parent_name = change.get(".destParentName")
    parent_index = _resolve_single(newMSL, parent_name, change, "PARENT")
    if parent_index is None or parent_index == index:
        if parent_index == index:
            log_change_error(change, "ERROR", "PARENT.SELF",
                             "Taxon cannot be its own parent", text(destination))
        return False
    exclude = [newMSL.at[index, ".prev_taxnode_id"]] if ".prev_taxnode_id" in newMSL else []
    if not _check_accessions(change, exclude):
        return False
    if action in {"move", "promote", "demote", "split"}:
        _log_move_parent_lineage(change, parent_index)
    previous_name = newMSL.at[index, "name"]
    previous_lineage = newMSL.at[index, "lineage"]
    out_target = destination if action == "rename" else change.get(".destLineage")
    _set_out_fields(index, change, action, out_target)
    newMSL.at[index, "name"] = destination
    newMSL.at[index, "rank"] = destination_rank
    newMSL.at[index, "level_id"] = dbCvMapList.get("rank", {}).get(str(destination_rank), newMSL.at[index, "level_id"])
    newMSL.at[index, "parent_id"] = newMSL.at[parent_index, "taxnode_id"]
    newMSL.at[index, "lineage"] = _lineage_from_parent(parent_index, str(destination))
    if action != "rename":
        _set_descriptives(index, change)
    previous_taxnode_id = newMSL.at[index, ".prev_taxnode_id"] if ".prev_taxnode_id" in newMSL else pd.NA
    if is_na(previous_taxnode_id):
        _append_in_fields(index, change)

    lineage_proposal = (
        f"{text(_proposal_zip(change))}:{text(change.get('.linenum'), '')}"
    )
    if action == "rename":
        alternate_lineage = previous_lineage
        lineage_action = (
            f"{action} {text(change.get('.srcRank'))} {text(previous_name)} "
            f"to {text(destination)}"
        )
    else:
        alternate_lineage = (
            f"{text(_lineage_metadata(parent_index, '.otherLineage'))};"
            f"{text(destination)}"
        )
        lineage_action = (
            f"{action} {text(change.get('.srcRank'))} "
            f"{diff_lineages(change.get('.srcLineage'), newMSL.at[index, 'lineage'])}"
        )
    update_lineage(
        newMSL.at[index, "taxnode_id"],
        newMSL.at[index, "lineage"],
        alternate_lineage,
        lineage_proposal,
        lineage_action,
    )
    _register_accessions(change)
    if params.tmi:
        print(f"  {action}: {previous_name} -> {destination}")
    return True


# ---------------------------------------------------------------------------
# Action result logging and dispatch
# ---------------------------------------------------------------------------
#
# apply_change maps each normalized action to an atomic operation and records a
# success only after mutation. apply_changes preserves the chosen dependency
# order and records whether each row was applied.


def _log_apply_success(change: pd.Series, action: str,
                       error_code: str | None = None,
                       notes_override: str | None = None) -> None:
    success_codes = {
        "new": "CREATE.OK", "split": "SPLIT.OK",
        "rename": "RENAME.OK", "abolish": "ABOLISH.OK",
        "move": "MOVE.OK", "promote": "PROMOTE.OK",
        "demote": "DEMOTE.OK", "merge": "MERGE_INTO.OK",
    }
    source = text(change.get(".srcTaxon"), "")
    destination = text(change.get(".destTaxon"), "")
    rank = text(change.get(".destRank"), text(change.get(".srcRank"), ""))
    if action == "new":
        notes = f"Create {rank} of '{text(change.get('.destLineage'), '')}'"
    elif action == "rename":
        notes = (f"RENAME {rank} from '{source}' to '{destination}' in "
                 f"'{text(change.get('.destParentLineage'), '')}'")
    elif action == "abolish":
        notes = f"ABOLISH {text(change.get('.srcRank'), '')} named {source}"
    else:
        notes = f"{action.upper()} {source} to {destination}"
    log_change_error(
        change, "SUCCESS",
        error_code or success_codes.get(action, f"{action.upper()}.OK"),
        f"Change={action.upper()}, applied successfully",
        notes_override if notes_override is not None else notes,
    )


def apply_change(change: pd.Series) -> bool:
    """Validate dispatch prerequisites and apply one ordered taxonomy change."""
    global actionOrder
    actionOrder = int(change.get(".changeOrder", actionOrder + 1))
    if not bool(change.get(".noErrors", True)):
        log_change_error(change, "ERROR", "ACTION.SKIPPED_QC",
                         "Change not applied because its XLSX row has errors",
                         text(change.get(".errors"), ""))
        return False
    action = change.get(".action")
    if is_na(action):
        return False
    action = str(action)
    source = change.get(".srcTaxon")
    destination = change.get(".destTaxon")
    success_error_code = None
    success_notes = None

    if action == "new":
        success = _create_taxon(change, action)
    elif action == "split" and (is_na(source) or str(source) != str(destination)):
        source_index = _resolve_single(newMSL, source, change, "SOURCE")
        if source_index is None:
            return False
        exclude = [newMSL.at[source_index, ".prev_taxnode_id"]] if ".prev_taxnode_id" in newMSL else []
        success = _create_taxon(change, action, exclude)
        if success:
            _set_out_fields(source_index, change, "split", change.get(".destLineage"))
    else:
        source_index = _resolve_single(newMSL, source, change, "SOURCE")
        if source_index is None:
            return False
        if action == "abolish":
            success = _remove_taxon(source_index, change, action, destination)
        elif action == "merge" and not is_na(destination):
            destination_matches = _matching_taxa(newMSL, destination)
            destination_matches = destination_matches[destination_matches != source_index]
            if len(destination_matches) == 1:
                success_error_code = "MERGE_INTO.OK"
                success_notes = _merge_into_notes(
                    source_index, destination_matches[0]
                )
                destination_taxnode_id = newMSL.at[destination_matches[0], "taxnode_id"]
                success = _remove_taxon(source_index, change, action, destination)
                if success:
                    destination_index = newMSL.index[newMSL["taxnode_id"].eq(destination_taxnode_id)][0]
                    previous_taxnode_id = newMSL.at[destination_index, ".prev_taxnode_id"]
                    if is_na(previous_taxnode_id):
                        _append_in_fields(destination_index, change)
                    _set_descriptives(destination_index, change)
                    _register_accessions(change)
            else:
                success = _move_or_rename_taxon(source_index, change, action)
        elif action in {"rename", "move", "promote", "demote", "split", "merge"}:
            success = _move_or_rename_taxon(source_index, change, action)
        else:
            log_change_error(change, "ERROR", "ACTION.UNIMPLEMENTED",
                             "Proposal action is not implemented", action)
            success = False
    if success:
        _log_apply_success(
            change, action, success_error_code, success_notes
        )
    if success and params.verbose:
        print(f"APPLIED {change.get('.codeLine')}: {action} {text(source)} -> {text(destination)}")
    return success


def apply_changes(changes: pd.DataFrame) -> pd.DataFrame:
    applied: list[bool] = []
    for _, change in changes.iterrows():
        applied.append(apply_change(change))
    result = changes.copy()
    result[".applied"] = applied
    return result


# ---------------------------------------------------------------------------
# Output schemas and file generation
# ---------------------------------------------------------------------------
#
# Output functions preserve the R filenames and database column order. MSL/SQL
# artifacts are optional under --msl; consolidated QC outputs are always
# written.

TSV_COLUMNS = [
    "msl_release_num", "level_id", "name", "ictv_id", "molecule_id",
    "abbrev_csv", "genbank_accession_csv", "genbank_refseq_accession_csv",
    "refseq_accession_csv", "isolate_csv", "notes", "in_change",
    "in_target", "in_filename", "in_notes", "out_change", "out_target",
    "out_filename", "out_notes", "lineage", "host_source", "exemplar_name",
    "genome_coverage", "notes",
]
SQL_COLUMNS = [
    "taxnode_id", "parent_id", "tree_id", "msl_release_num", "level_id",
    "name", "ictv_id", "molecule_id", "abbrev_csv", "genbank_accession_csv",
    "genbank_refseq_accession_csv", "refseq_accession_csv", "isolate_csv",
    "notes", "is_hidden", "in_change", "in_target", "in_filename",
    "in_notes", "out_change", "out_target", "out_filename", "out_notes",
    "host_source", "exemplar_name", "genome_coverage",
]


def _ensure_columns(frame: pd.DataFrame, columns: Sequence[str]) -> pd.DataFrame:
    result = frame.copy()
    for column in dict.fromkeys(columns):
        if column not in result:
            result[column] = pd.NA
    return result


def _short_proposal_filename(value: Any) -> Any:
    if is_na(value):
        return value
    return re.sub(r"^([0-9]+\.[0-9A-Z]+)\..*", r"\1...", str(value))


# ---------------------------------------------------------------------------
# MSL TSV and proposal-metadata exports
# ---------------------------------------------------------------------------


def export_msl_tsv() -> None:
    frame = _ensure_columns(newMSL, TSV_COLUMNS)
    selected = frame.loc[:, TSV_COLUMNS].copy()
    selected["in_filename"] = selected["in_filename"].map(_short_proposal_filename)
    selected["out_filename"] = selected["out_filename"].map(_short_proposal_filename)
    filename = Path(params.out_dir) / params.msl_tsv
    print(f"Writing {filename}")
    selected.to_csv(filename, sep="\t", index=False, header=True, na_rep="NULL",
                    quoting=csv.QUOTE_NONE, escapechar="\\", lineterminator="\n")
    print(f"WROTE   {filename} ({len(selected)} rows)")


def export_docx_summary() -> None:
    columns = ["subcommittee", "code", "docx", "xlsx", "title",
               "authorsEmails", "correspondingAuthor", "abstract"]
    frame = proposalsDf.reset_index(drop=True).copy()
    frame = _ensure_columns(frame, columns).loc[:, columns]
    filename = Path(params.out_dir) / params.proposals_meta
    print(f"Writing {filename}")
    frame.to_csv(filename, sep="\t", index=False, na_rep="",
                 quoting=csv.QUOTE_NONE, escapechar="\\", lineterminator="\n")
    print(f"WROTE   {filename} ({len(frame)})")


# ---------------------------------------------------------------------------
# MariaDB load-script generation
# ---------------------------------------------------------------------------
#
# SQL output batches taxonomy inserts, updates previous-MSL out_* fields,
# rebuilds indexes/deltas, refreshes VMR, and only then runs database QC.


def sql_value_nullable(value: Any) -> str:
    if is_na(value):
        return "NULL"
    return "'" + str(value).replace("'", "''") + "'"


def _sql_index_rebuild(tree_id: int) -> str:
    ranks = [rank.lower() for rank in xlsx_change_ranks]
    identifiers = ",\n    ".join(f"@{rank}_id" for rank in ranks)
    descendant_counts = ",\n    ".join(f"@{rank}_desc_ct" for rank in ranks)
    variables = ", ".join(f"@{rank}_id := NULL" for rank in ranks)
    count_variables = ", ".join(f"@{rank}_desc_ct := 0" for rank in ranks)
    return f"""
-- set for taxonomy_node_compute_indexes
SET SESSION max_sp_recursion_depth = 255;
SET @right_idx := NULL;
SET {variables}, @inher_molecule_id := NULL, @lineage := NULL;
SET {count_variables};

-- rebuilt indexes
CALL `taxonomy_node_compute_indexes`(
    {tree_id}, 1, @right_idx, 1,
    {identifiers},
    {descendant_counts},
    @inher_molecule_id, @lineage
);
"""


def export_sql() -> None:
    """Write the ordered MariaDB load and post-load maintenance script."""
    frame = _ensure_columns(newMSL, SQL_COLUMNS + ["lineage"])
    sql_msl_num = int(pd.to_numeric(frame["msl_release_num"], errors="coerce").dropna().iloc[0])
    tree_id = int(pd.to_numeric(frame["tree_id"], errors="coerce").dropna().max())
    filename = Path(params.out_dir) / params.sql_load_filename
    print(f"Writing {filename}")
    lines = [
        "-- begin transaction\n-- rollback transaction\n",
        "SET SESSION sql_mode = CONCAT(@@SESSION.sql_mode, ',NO_BACKSLASH_ESCAPES');\n\n",
        f"insert into `taxonomy_toc` (`tree_id`,`msl_release_num`,`comments`) values ({tree_id},{sql_msl_num},NULL);\n",
    ]
    ordered = frame.assign(_level=pd.to_numeric(frame["level_id"], errors="coerce")).sort_values("_level", kind="stable")
    batch_size = max(1, int(params.sql_insert_batch_size))
    for offset in range(0, len(ordered), batch_size):
        batch = ordered.iloc[offset:offset + batch_size]
        lines.append("insert into `taxonomy_node` (`" + "`,`".join(SQL_COLUMNS) + "`)\n values \n")
        row_lines = []
        for _, row in batch.iterrows():
            values_sql = ",".join(sql_value_nullable(row.get(column)) for column in SQL_COLUMNS)
            row_lines.append(f"({values_sql}) -- lineage={text(row.get('lineage'), '')}")
        lines.append(",\n".join(row_lines) + ";\n")

    current_columns = ["out_change", "out_target", "out_filename", "out_notes"]
    current = _ensure_columns(curMSL, current_columns)
    current = current.loc[current["out_change"].notna()]
    lines.append("\n-- updates for previous MSL\n")
    for _, row in current.iterrows():
        assignments = ",".join(
            f"`{column}`={sql_value_nullable(row.get(column))}" for column in current_columns
        )
        lines.append(f"update `taxonomy_node` set {assignments} where `taxnode_id`={int(row['taxnode_id'])};\n")
    lines.append(_sql_index_rebuild(tree_id))
    lines.append("""
-- build deltas and merge/split relationships
CALL `rebuild_delta_nodes`(NULL);
CALL `rebuild_node_merge_split`();

-- Refresh VMR from the newly loaded MSL before QC. If a deployment defines
-- parameters for this procedure, adjust this single call to that signature.
CALL `VMR_update_from_MSL`();

-- QC SPs (must remain after VMR_update_from_MSL)
CALL `QC_run_modules`(NULL);
""")
    filename.write_text("".join(lines), encoding="utf-8")
    print(f"WROTE   {filename}")


# ---------------------------------------------------------------------------
# Optional change report and final output orchestration
# ---------------------------------------------------------------------------


def export_change_report(changes: pd.DataFrame) -> None:
    if not params.output_change_report:
        return
    filename = Path(params.out_dir) / "changes.tsv"
    public_columns = [column for column in changes if not column.startswith(".void")]
    write_tsv(changes.loc[:, public_columns], filename)
    if params.verbose:
        print(f"WROTE   {filename} ({len(changes)} rows)")


def export_outputs(changes: pd.DataFrame) -> None:
    # Remove per-subcommittee workbooks produced by older Python runs. The
    # consolidated all/issue reports are the supported workbook outputs.
    for stale in Path(params.out_dir).glob("QC.pretty_summary.?_*.xlsx"):
        stale.unlink(missing_ok=True)
    if params.export_msl:
        export_msl_tsv()
        export_docx_summary()
        export_sql()
    export_change_report(changes)
    write_error_summary(allErrorDf, final=True)


# ---------------------------------------------------------------------------
# Program entry point and end-to-end lifecycle
# ---------------------------------------------------------------------------
#
# Lifecycle: parse options -> load references -> discover/QC proposals -> order
# and apply changes -> export MSL/SQL/QC artifacts.


def main(argv: Sequence[str] | None = None) -> int:
    """Run one complete proposal validation and MSL-generation workflow."""
    global params, changeList, allChangeDf, allErrorDf, proposedAccessions
    # 1. Initialize invocation-specific options and mutable result state.
    params = parse_args(argv)
    Path(params.out_dir).mkdir(parents=True, exist_ok=True)
    allErrorDf = pd.DataFrame(columns=ERROR_COLUMNS)
    proposedAccessions = {}
    if params.use_cache or params.update_cache or params.load_proposal_cache or params.save_proposal_cache:
        print("WARNING: cache flags are accepted for compatibility but caching is not implemented in this first Python version")
    # 2. Load validator/reference data before reading proposal content.
    load_version()
    load_reference()
    load_abolished_accessions()

    # 3. Discover, load, and independently QC each proposal pair.
    proposals = scan_for_proposals()
    changeList = load_and_qc_proposals(proposals, {})
    # 4. Resolve cross-proposal dependencies, then mutate the working MSL.
    allChangeDf = merge_and_order_changes(changeList)
    if allChangeDf.empty:
        log_error(pd.NA, 0, "", "", "", "ERROR", "NO_CHANGES",
                  "No valid proposal changes were loaded", "")
    else:
        allChangeDf = apply_changes(allChangeDf)
    # 5. Write requested database artifacts and the always-on QC reports.
    export_outputs(allChangeDf)
    print("# COMPLETED.")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        raise SystemExit(130)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        if params.get("debug_on_error"):
            raise
        raise SystemExit(1)
