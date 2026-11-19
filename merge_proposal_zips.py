#!/usr/bin/env python3

"""Native Python proposal validator converted from merge_proposal_zips.R.

The R script is the reference implementation. This file follows its stages,
function names, error codes, messages and notes so that QC output can be
compared with R row by row:

    load_reference -> scan_for_proposals -> load_and_qc_proposals
    -> merge/sort/chain the changes -> apply_changes -> exports

Reference and proposal caches are intentionally not implemented.

Behavior that exists only here, or that deliberately does not copy an R
defect, is listed in PYTHON_ONLY below. Each item is on by default and can be
switched off for equivalence testing with the environment variable
MERGE_PROPOSAL_ZIPS_R_COMPAT (see r_compat()).
"""

from __future__ import annotations

import argparse
import datetime
import os
import re
import subprocess
import sys
import unicodedata
import warnings
import zipfile
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence
from xml.etree import ElementTree

from openpyxl import Workbook

try:
    import pandas as pd
except ImportError as exc:
    raise SystemExit(
        "ERROR: pandas is required. Create the environment with "
        "`conda env create -f environment.python.yml`."
    ) from exc


# ---------------------------------------------------------------------------
# Behavior that differs from merge_proposal_zips.R on purpose
# ---------------------------------------------------------------------------
#
# Every entry is a decision for the maintainer, not a porting gap. With the
# default (all on) Python keeps the behavior it had before the 2026-11 sync.
# For equivalence testing set, for example:
#
#     MERGE_PROPOSAL_ZIPS_R_COMPAT=all
#     MERGE_PROPOSAL_ZIPS_R_COMPAT=in_change_text,move_lineage_text
#
# and the named items behave exactly as the R script does.

PYTHON_ONLY = {
    # Supplemental, nonblocking QC: R skips E1 for modern templates and
    # accession validation for rename actions. Neither check alters taxonomy.
    "modern_header_code": True,
    "rename_duplicate_accession": True,
    # msl.tsv / the SQL lineage comment get the lineage rebuilt from parent_id.
    # R leaves the lineage text of a moved/promoted/demoted taxon and of its
    # descendants stale (open defect), and keeps the tree name in the lineage
    # of a taxon renamed directly under the root.
    "move_lineage_text": True,
    # A "Move; rename" onto a name another taxon already has is rejected with
    # DEST.EXISTS. R skips every name check for a move (open defect).
    "move_rename_name_check": True,
    # Input files whose names contain spaces are renamed in place. R only logs
    # DOCX_FILENAME_SPACES / XLSX_FILENAME_SPACES.
    "rename_spaced_filenames": True,
    # msl.tsv in_change holds the word ("new", "split"). R writes the factor
    # level number instead (1, 2, ...), an R defect in the TSV export only.
    "in_change_text": True,
    # A "please select" cell becomes NA by itself. In R one such cell (in a
    # spelling read_excel does not already drop) blanks its whole column.
    "please_select_per_cell": True,
    # A space that is a rich-text run of its own is kept ("caudata (long-tailed").
    # readxl drops such runs, so R reads "caudata(long-tailed".
    "rich_text_spaces": True,
}


def r_compat() -> None:
    """Switch PYTHON_ONLY items off as MERGE_PROPOSAL_ZIPS_R_COMPAT asks."""
    requested = os.environ.get("MERGE_PROPOSAL_ZIPS_R_COMPAT", "").strip()
    if not requested:
        return
    names = [name.strip() for name in requested.split(",") if name.strip()]
    if "all" in names:
        names = list(PYTHON_ONLY)
    unknown = [name for name in names if name not in PYTHON_ONLY]
    if unknown:
        raise RuntimeError(
            "MERGE_PROPOSAL_ZIPS_R_COMPAT: unknown item(s) " + ",".join(unknown)
            + "; known: all," + ",".join(PYTHON_ONLY))
    for name in names:
        PYTHON_ONLY[name] = False
    print("R_COMPAT: Python-only behavior switched off: " + ",".join(names))


# ---------------------------------------------------------------------------
# R-compatible process state
# ---------------------------------------------------------------------------
#
# These module-level objects intentionally retain the R script names.

actionOrder = 0
cvList: dict[str, Any] = {}
dbCvList: dict[str, Any] = {}
dbCvMapList: dict[str, Any] = {}
scAbbrevNameMap: dict[Any, Any] = {}
rankNames: list[str] = []
proposalsDf: dict[str, dict[str, Any]] = {}
proposalsDfColumns: list[str] = []
changeList: dict[str, list[dict[str, Any]]] = {}
allChangeDf: list[dict[str, Any]] = []
vmrAccessions: dict[str, list[tuple[int, Any, Any]]] = {}
abolishedAccessions: list[dict[str, Any]] = []
inChangeLevels: list[str] = []

ERROR_COLUMNS = [
    "subcommittee", "code", "docx", "xlsx", "row", "change", "rank",
    "taxon", "level", "error", "message", "notes", "scAbbrev",
    "validator_version", "order",
]


# ---------------------------------------------------------------------------
# Command-line compatibility
# ---------------------------------------------------------------------------


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
        # R f84c1c5: "suppl" became "suppl[a-z0-9]*"
        default=r"(^|[_. -])(suppl[a-z0-9]*|viridic|table|matrix|sup)([_. -]|$)",
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
# R value semantics
# ---------------------------------------------------------------------------
#
# Missing values are None everywhere in this file. rstr()/paste0() format a
# value the way R's paste0() does, so notes read "NA" where R's do.


def is_na(value: Any) -> bool:
    if value is None:
        return True
    if isinstance(value, float):
        return value != value
    if isinstance(value, (str, int, bool, list, tuple, dict)):
        return False
    try:
        return bool(pd.isna(value))
    except (TypeError, ValueError):
        return False


def rstr(value: Any) -> str:
    """as.character() for one value."""
    if is_na(value):
        return "NA"
    if isinstance(value, bool):
        return "TRUE" if value else "FALSE"
    if isinstance(value, float):
        if value.is_integer() and abs(value) < 1e15:
            return str(int(value))
        return format(value, ".15g")
    return str(value)


def paste0(*parts: Any) -> str:
    return "".join(rstr(part) for part in parts)


def text(value: Any, na: str = "NA") -> str:
    return na if is_na(value) else rstr(value)


def r_strsplit(value: str, sep: str = ";") -> list[str]:
    """strsplit(value, sep)[[1]]: one trailing empty piece is dropped."""
    parts = value.split(sep)
    if parts and parts[-1] == "":
        parts.pop()
    return parts


_ICU_PUNCTUATION = "_-,;:!?.'\"()[]{}@*/\\&#%`^+<=>|~$"


def _collation_is_c() -> bool:
    """R compares strings bytewise only in the C/POSIX locale."""
    for name in ("LC_ALL", "LC_COLLATE", "LANG"):
        value = os.environ.get(name)
        if value:
            return value in ("C", "POSIX")
    return True


def r_collate_key(value: Any) -> tuple:
    """Sort key for the order R's sort()/order()/list.files() give strings.

    Outside the C locale R collates with ICU: spaces, then punctuation, then
    digits, then letters without regard to case or accents; ties go to fewer
    accents, then lower case. NA sorts last.
    """
    if is_na(value):
        return (1,)
    string = rstr(value)
    if _collation_is_c():
        return (0, string)
    primary: list[tuple[int, Any]] = []
    accents: list[int] = []
    cases: list[int] = []
    for char in unicodedata.normalize("NFD", string):
        if unicodedata.combining(char):
            if accents:
                accents[-1] += 1
            continue
        if char.isspace():
            primary.append((0, 0))
        elif char in _ICU_PUNCTUATION:
            primary.append((1, _ICU_PUNCTUATION.index(char)))
        elif char.isdigit():
            primary.append((3, char))
        elif char.isalpha():
            primary.append((4, char.lower()))
        else:
            primary.append((2, char))
        accents.append(0)
        cases.append(1 if char.isupper() else 0)
    return (0, primary, accents, cases, string)


def normalize_posix_regex(pattern: str) -> str:
    """Translate the POSIX classes R patterns use to Python's re syntax."""
    for source, target in {
        "[:alnum:]": "A-Za-z0-9", "[:alpha:]": "A-Za-z", "[:digit:]": "0-9",
        "[:space:]": r"\s", "[:upper:]": "A-Z", "[:lower:]": "a-z",
        "[:print:]": r"\x20-\x7E", "[:punct:]": r"!-/:-@\[-`{-~",
    }.items():
        pattern = pattern.replace(source, target)
    return pattern


# ---------------------------------------------------------------------------
# Tables: rows are dicts, missing columns read as NA
# ---------------------------------------------------------------------------


class Table:
    """A small ordered table (list of dict rows) standing in for a data.frame."""

    def __init__(self, columns: Iterable[str], rows: list[dict[str, Any]]):
        self.columns = list(columns)
        self.rows = rows

    def add_column(self, name: str) -> None:
        if name not in self.columns:
            self.columns.append(name)


class NewMsl(Table):
    """The working MSL, with the lookups apply_changes makes kept as indexes.

    Rows keep the order R's newMSL has: the copied MSL first, created taxa
    appended. Every index returns rows in that order. Names, parents and
    membership must only change through set_name/set_parent/append/remove.
    """

    def __init__(self, columns: Iterable[str], rows: list[dict[str, Any]]):
        super().__init__(columns, rows)
        self._next_pos = 0
        self.by_name: dict[Any, list[dict[str, Any]]] = {}
        self.by_lower: dict[Any, list[dict[str, Any]]] = {}
        self.by_id: dict[Any, list[dict[str, Any]]] = {}
        self.by_parent: dict[Any, list[dict[str, Any]]] = {}
        self.by_prev: dict[Any, list[dict[str, Any]]] = {}
        for row in rows:
            self._index(row)

    @staticmethod
    def _lower(name: Any) -> Any:
        return None if is_na(name) else str(name).lower()

    @staticmethod
    def _insert(index: dict, key: Any, row: dict[str, Any]) -> None:
        if is_na(key):
            return
        bucket = index.setdefault(key, [])
        bucket.append(row)
        if len(bucket) > 1 and bucket[-2]["_pos"] > row["_pos"]:
            bucket.sort(key=lambda item: item["_pos"])

    @staticmethod
    def _drop(index: dict, key: Any, row: dict[str, Any]) -> None:
        bucket = index.get(key)
        if bucket is None:
            return
        for i, item in enumerate(bucket):
            if item is row:
                del bucket[i]
                break
        if not bucket:
            del index[key]

    def _index(self, row: dict[str, Any]) -> None:
        row["_pos"] = self._next_pos
        self._next_pos += 1
        self._insert(self.by_name, row.get("name"), row)
        self._insert(self.by_lower, self._lower(row.get("name")), row)
        self._insert(self.by_id, row.get("taxnode_id"), row)
        self._insert(self.by_parent, row.get("parent_id"), row)
        self._insert(self.by_prev, row.get(".prev_taxnode_id"), row)

    def append(self, row: dict[str, Any]) -> None:
        self.rows.append(row)
        self._index(row)

    def remove(self, row: dict[str, Any]) -> None:
        self.rows = [item for item in self.rows if item is not row]
        self._drop(self.by_name, row.get("name"), row)
        self._drop(self.by_lower, self._lower(row.get("name")), row)
        self._drop(self.by_id, row.get("taxnode_id"), row)
        self._drop(self.by_parent, row.get("parent_id"), row)
        self._drop(self.by_prev, row.get(".prev_taxnode_id"), row)

    def set_name(self, row: dict[str, Any], name: Any) -> None:
        self._drop(self.by_name, row.get("name"), row)
        self._drop(self.by_lower, self._lower(row.get("name")), row)
        row["name"] = name
        self._insert(self.by_name, name, row)
        self._insert(self.by_lower, self._lower(name), row)

    def set_parent(self, row: dict[str, Any], parent_id: Any) -> None:
        self._drop(self.by_parent, row.get("parent_id"), row)
        row["parent_id"] = parent_id
        self._insert(self.by_parent, parent_id, row)

    def named(self, name: Any) -> list[dict[str, Any]]:
        return [] if is_na(name) else list(self.by_name.get(str(name), []))

    def with_id(self, taxnode_id: Any) -> list[dict[str, Any]]:
        return [] if is_na(taxnode_id) else list(self.by_id.get(taxnode_id, []))

    def kids(self, taxnode_id: Any) -> list[dict[str, Any]]:
        return [] if is_na(taxnode_id) else list(self.by_parent.get(taxnode_id, []))

    def from_prev(self, prev_id: Any) -> list[dict[str, Any]]:
        return [] if is_na(prev_id) else list(self.by_prev.get(prev_id, []))


oldMSLs = Table([], [])
curMSL = Table([], [])
newMSL = NewMsl([], [])


# ---------------------------------------------------------------------------
# Input format adapters
# ---------------------------------------------------------------------------
#
# Database exports require a dedicated parser so unquoted \N/NULL, quoted empty
# strings and MySQL backslash escapes remain distinct. Excel cells are read the
# way readxl::read_excel(col_names=FALSE, trim_ws=TRUE) returns them.


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
            fields.append(None if (not was_quoted and raw in {"\\N", "NULL"}) else raw)
            value = []
            was_quoted = False
        else:
            value.append(char)
        index += 1
    raw = "".join(value)
    fields.append(None if (not was_quoted and raw in {"\\N", "NULL"}) else raw)
    return fields


def read_mysql_tsv(path: str | Path) -> Table:
    """Read old and new MySQL TSV formats without conflating empty and NULL."""
    source = Path(path)
    with source.open("r", encoding="utf-8", newline="") as stream:
        header_line = stream.readline()
        if not header_line:
            return Table([], [])
        header = [text(value, "") for value in _parse_mysql_record(header_line)]
        rows: list[dict[str, Any]] = []
        for line_number, line in enumerate(stream, 2):
            if line.strip("\r\n") == "":
                continue
            parsed = _parse_mysql_record(line)
            if len(parsed) != len(header):
                raise RuntimeError(
                    f"ERROR: {source}:{line_number} has {len(parsed)} fields; expected {len(header)}"
                )
            rows.append(dict(zip(header, parsed)))
    return Table(header, rows)


def _cell_text(value: Any, na_strings: Sequence[str]) -> Any:
    """One worksheet cell as readxl returns it in a text column."""
    if is_na(value):
        return None
    if isinstance(value, bool):
        return "TRUE" if value else "FALSE"
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        return rstr(value)
    if isinstance(value, (datetime.datetime, datetime.date, datetime.time)):
        return str(value)
    result = str(value)
    # readxl: a cell of only whitespace is blank; trim_ws trims only spaces
    # and tabs, so a line break at either end is kept
    if result.strip(" \t\r\n\f\v") == "":
        return None
    result = result.strip(" \t")
    if result in na_strings:
        return None
    return result


# the text of a <t> element in sharedStrings.xml
_SHARED_STRING_TEXT = re.compile(rb"(<(?:\w+:)?t(?:\s[^>]*)?>)(.*?)(</(?:\w+:)?t>)", re.S)


def _readxl_shared_strings(path: Path, drop_space_runs: bool) -> dict[str, str]:
    """Shared strings that readxl returns differently from openpyxl.

    Keys are the strings as openpyxl returns them, values as readxl does.
    openpyxl's XML parser turns a stored CRLF (or CR) into LF, as XML
    requires; readxl keeps the CR. The CRs are recovered by parsing a copy in
    which each CR inside a <t> element is the character reference &#13;.

    With drop_space_runs, a rich-text run that holds only whitespace is also
    dropped, even one marked xml:space="preserve", as readxl does. Used only
    when PYTHON_ONLY["rich_text_spaces"] is off.
    """
    namespace = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
    try:
        with zipfile.ZipFile(path) as archive:
            raw = archive.read("xl/sharedStrings.xml")
    except (KeyError, OSError, zipfile.BadZipFile):
        return {}
    if b"\r" not in raw and not drop_space_runs:
        return {}
    stored = _SHARED_STRING_TEXT.sub(
        lambda match: match.group(1) + match.group(2).replace(b"\r", b"&#13;") + match.group(3), raw)
    try:
        openpyxl_root = ElementTree.fromstring(raw)
        readxl_root = ElementTree.fromstring(stored)
    except ElementTree.ParseError:
        return {}

    def runs(item: ElementTree.Element) -> list[str]:
        texts = [node.text or "" for node in item.findall(namespace + "t")]
        for run in item.findall(namespace + "r"):
            node = run.find(namespace + "t")
            if node is not None:
                texts.append(node.text or "")
        return texts

    # every readxl value seen for each openpyxl string; a string that two
    # different stored strings share is left as openpyxl returns it
    values: dict[str, set[str]] = {}
    for seen, item in zip(openpyxl_root.iter(namespace + "si"), readxl_root.iter(namespace + "si")):
        texts = runs(item)
        if drop_space_runs:
            texts = [text for text in texts if text.strip(" \t\r\n") != ""]
        values.setdefault("".join(runs(seen)), set()).add("".join(texts))
    return {key: next(iter(found)) for key, found in values.items()
            if len(found) == 1 and key not in found}


def excel_sheets(path: str | Path) -> list[str]:
    engine = "xlrd" if Path(str(path)).suffix.lower() == ".xls" else "openpyxl"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with pd.ExcelFile(path, engine=engine) as excel:
            return list(excel.sheet_names)


def read_excel_raw(path: str | Path, sheet: str, na_strings: Sequence[str] = (),
                   columns: int | None = None) -> list[list[Any]]:
    """Rows of a worksheet as text cells; row i is Excel row i+1."""
    source = Path(path)
    engine = "xlrd" if source.suffix.lower() == ".xls" else "openpyxl"
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            # cells such as "NA" and "N/A" are text, not missing values
            frame = pd.read_excel(path, sheet_name=sheet, header=None,
                                  engine=engine, dtype=object,
                                  keep_default_na=False, na_filter=False)
    except ImportError as exc:
        raise RuntimeError("Legacy .xls input requires the xlrd package") from exc
    width = columns if columns is not None else frame.shape[1]
    as_readxl = ({} if engine != "openpyxl"
                 else _readxl_shared_strings(source, drop_space_runs=not PYTHON_ONLY["rich_text_spaces"]))
    rows: list[list[Any]] = []
    for values in frame.itertuples(index=False, name=None):
        cells = [_cell_text(as_readxl.get(value, value) if isinstance(value, str) else value,
                            na_strings) for value in values[:width]]
        cells.extend([None] * (width - len(cells)))
        rows.append(cells)
    # readxl stops at the last row that has data
    while rows and all(cell is None for cell in rows[-1]):
        rows.pop()
    return rows


def _tsv_field(value: Any) -> str:
    """One field as readr::write_delim(delim='\\t') writes it."""
    if is_na(value):
        return "NA"
    field = rstr(value)
    if any(char in field for char in '\t"\n\r'):
        return '"' + field.replace('"', '""') + '"'
    return field


def write_tsv(columns: Sequence[str], rows: Iterable[Mapping[str, Any]],
              path: str | Path) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as out:
        out.write("\t".join(_tsv_field(column) for column in columns) + "\n")
        for row in rows:
            out.write("\t".join(_tsv_field(row.get(column)) for column in columns) + "\n")


def write_xlsx_sheet(sheet_name: str, columns: Sequence[str],
                     rows: Iterable[Mapping[str, Any]], path: str | Path) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    workbook = Workbook(write_only=True)
    sheet = workbook.create_sheet(title=str(sheet_name)[:31])
    sheet.append(list(columns))
    illegal = re.compile(r"[\x00-\x08\x0b\x0c\x0e-\x1f]")
    for row in rows:
        cells = []
        for column in columns:
            value = row.get(column)
            if is_na(value):
                cells.append(None)
            elif isinstance(value, str):
                cells.append(illegal.sub("", value))
            else:
                cells.append(value)
        sheet.append(cells)
    workbook.save(path)


# ---------------------------------------------------------------------------
# Validator versioning and QC report collection
# ---------------------------------------------------------------------------


class ErrorLog:
    """allErrorDf: rows plus the column order rbindlist(fill=TRUE) builds."""

    def __init__(self) -> None:
        self.columns = list(ERROR_COLUMNS)
        self.rows: list[dict[str, Any]] = []
        # factor levels of the level column, in R's order of appearance
        self.levels = ["ERROR", "WARNING", "INFO"]

    def append(self, row: Mapping[str, Any]) -> None:
        for column in row:
            if column not in self.columns:
                self.columns.append(column)
        level = row.get("level")
        if not is_na(level) and level not in self.levels:
            self.levels.append(level)
        self.rows.append(dict(row))

    def __len__(self) -> int:
        return len(self.rows)


allErrorDf = ErrorLog()
_last_summary_state: tuple[int, int] | None = None


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


def _proposal_value(code: Any, column: str) -> Any:
    if is_na(code):
        return None
    return proposalsDf.get(str(code), {}).get(column)


def _as_numeric(value: Any) -> Any:
    """as.numeric() of a line number: NA when it is not a number."""
    if is_na(value):
        return None
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return int(value) if float(value).is_integer() else value
    try:
        number = float(str(value).strip())
    except ValueError:
        return None
    return int(number) if number.is_integer() else number


def _vector(value: Any) -> list[Any]:
    return list(value) if isinstance(value, (list, tuple)) else [value]


def log_error(code: Any, linenum: Any, action: Any, rank: Any, taxon: Any,
              levelStr: Any, errorCode: Any, errorStr: Any, notes: Any = "",
              actionOrder: Any = "") -> None:
    """Append one row per element of the longest argument (R recycling)."""
    vectors = [_vector(item) for item in
               (code, linenum, action, rank, taxon, levelStr, errorCode, errorStr, notes)]
    if any(len(vector) == 0 for vector in vectors):
        return
    size = max(len(vector) for vector in vectors)
    order = actionOrder
    if is_na(order) or order == "":
        order = globals()["actionOrder"]
    for i in range(size):
        (one_code, one_line, one_action, one_rank, one_taxon, one_level,
         one_error, one_message, one_notes) = (vector[i % len(vector)] for vector in vectors)
        docx = _proposal_value(one_code, "docx")
        xlsx = _proposal_value(one_code, "xlsx")
        allErrorDf.append({
            "subcommittee": _proposal_value(one_code, "subcommittee"),
            "code": one_code,
            "row": _as_numeric(one_line),
            "order": order,
            "change": one_action, "rank": one_rank, "taxon": one_taxon,
            "docx": "MISSING" if is_na(docx) else docx,
            "xlsx": "MISSING" if is_na(xlsx) else xlsx,
            "level": one_level, "error": one_error, "message": one_message,
            "notes": one_notes,
        })


def log_change_error(curChangeDf: Mapping[str, Any], levelStr: str, errorCode: str,
                     errorStr: Any, notes: Any = "") -> None:
    log_error(curChangeDf.get(".code"), curChangeDf.get(".linenum"),
              curChangeDf.get("change"), curChangeDf.get("rank"),
              curChangeDf.get(".changeTaxon"), levelStr, errorCode, errorStr,
              notes, curChangeDf.get(".order"))


def _sorted_errors() -> list[dict[str, Any]]:
    """order(code, as.integer(row), as.integer(order)), NA last, stable."""
    def key(row: dict[str, Any]) -> tuple:
        code, line, order = row.get("code"), _as_numeric(row.get("row")), _as_numeric(row.get("order"))
        return (r_collate_key(code),
                is_na(line), 0 if is_na(line) else int(line),
                is_na(order), 0 if is_na(order) else int(order))
    return sorted(allErrorDf.rows, key=key)


def _issue_action_item(row: Mapping[str, Any]) -> str:
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


def write_error_summary(errorDf: ErrorLog | None = None, final: bool = False) -> None:
    """Write the pretty workbook, QC.summary.tsv and the regression TSV.

    Called often, as in R, so an artifact exists if a later step fails. When
    no row or column was added since the last call the files are already
    current and are not rewritten; the "Wrote:" line is still printed.
    Per-subcommittee workbooks (R's final mode) are intentionally not written.
    The issues-only workbook exists only in Python and is written with the
    final summary.
    """
    global _last_summary_state
    del errorDf
    columns = list(allErrorDf.columns)
    state = (len(allErrorDf.rows), len(columns))
    summary_path = Path(params.out_dir) / params.qc_summary_tsv_fname
    xlsx_path = Path(params.out_dir) / params.qc_summary_fname
    regression_path = Path(params.out_dir) / params.qc_regression_tsv_fname
    public_columns = [column for column in columns
                      if not re.search(r"(code|docx|subcommittee)", column)]
    if state != _last_summary_state or final:
        version = params.get("version")
        ordered = _sorted_errors()
        for row in ordered:
            row["validator_version"] = version
        pretty_rows: list[Mapping[str, Any]] = []
        previous = ""
        for row in ordered:
            code = row.get("code")
            if not is_na(code) and code != previous:
                previous = code
                pretty_rows.append({"subcommittee": row.get("subcommittee")})
                pretty_rows.append({
                    "subcommittee": row.get("subcommittee"), "code": code,
                    "xlsx": row.get("docx") if is_na(row.get("xlsx")) else row.get("xlsx"),
                })
            pretty_rows.append(row)
        write_xlsx_sheet("QC_Report", public_columns, pretty_rows, xlsx_path)
        write_tsv(columns, ordered, summary_path)
        write_tsv([column for column in columns if "version" not in column],
                  ordered, regression_path)
        if final:
            issues = []
            for row in ordered:
                if row.get("level") in ("ERROR", "WARNING"):
                    issue = dict(row)
                    issue["action_item"] = _issue_action_item(row)
                    issues.append(issue)
            issue_columns = list(public_columns)
            position = issue_columns.index("notes") + 1 if "notes" in issue_columns else len(issue_columns)
            issue_columns.insert(position, "action_item")
            write_xlsx_sheet("QC_Issues", issue_columns, issues,
                             Path(params.out_dir) / "QC.pretty_summary.issues.xlsx")
        _last_summary_state = state
    if params.verbose or params.tmi:
        print(f"Wrote: {xlsx_path} with 1 sheets")
        print(f"Wrote: {summary_path} ({state[0]} rows)")
        print(f"Wrote: {regression_path} ({state[0]} rows)")
    else:
        print(f"Wrote: OUTDIR/{params.qc_summary_tsv_fname} ({state[0]} rows)")


# ---------------------------------------------------------------------------
# Proposal template schema and validation rules
# ---------------------------------------------------------------------------

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
xlsx_change_srcDest_colnames = xlsx_change_src_colnames + xlsx_change_dest_colnames

# 0-based worksheet columns for each template layout
XLSX_V1_CHANGE_COLS = list(range(0, 15)) + list(range(16, 31)) + list(range(31, 41))
XLSX_V2_CHANGE_COLS = list(range(1, 16)) + list(range(16, 31)) + list(range(31, 41))
XLSX_2023_CHANGE_COLS = list(range(0, 15)) + list(range(15, 30)) + list(range(30, 40))

# The 2023.x, 2024.1, and 2026.1 templates share this worksheet layout.
# Later revisions must be added explicitly after their headers are verified.
XLSX_2023_STYLE_EXACT_VERSIONS = ("2024.1", "2026.1")

_NA = None
XLSX_V1_ROW2 = (["CURRENT TAXONOMY"] + [_NA] * 15 + ["PROPOSED TAXONOMY"] + [_NA] * 21
                + ["SPECIFY PROPOSED CHANGE", _NA, "COMMENTS"] + [_NA] * 8)
XLSX_V2_ROW2 = (["CURRENT TAXONOMY"] + [_NA] * 14 + ["PROPOSED TAXONOMY"] + [_NA] * 21
                + ["SPECIFY PROPOSED CHANGE", _NA, "COMMENTS"])
XLSX_V1_ROW3 = xlsx_change_ranks + ["Exemplar GenBank Accession Number"] + xlsx_change_ranks + [
    "Exemplar GenBank Accession Number", "Exemplar \r\nvirus name", "Virus name abbrevn",
    "Exemplar\r\nisolate designation", "Genome coverage", "Genome composition",
    "Host/Source", "Change", "Rank",
]
XLSX_V2_ROW3 = [_NA] + xlsx_change_ranks + xlsx_change_ranks + [
    "Exemplar GenBank Accession Number", "Exemplar \r\nvirus name", "Virus name abbreviation",
    "Exemplar\r\nisolate designation", "Genome coverage", "Genome composition",
    "Host/Source", "Change", "Rank", "Comments",
]
XLSX_2023_ROW3 = (["CURRENT TAXONOMY"] + [_NA] * 14 + ["PROPOSED TAXONOMY"] + [_NA] * 14
                  + ["DESCRIPTIVES"] + [_NA] * 6 + ["ACTION", _NA, "COMMENTS"])
XLSX_2023_ROW4 = xlsx_change_ranks + xlsx_change_ranks + [
    "Exemplar GenBank Accession Number", "Exemplar virus name",
    "Virus name abbreviation", "Exemplar isolate designation",
    "Genome coverage", "Genome composition", "Host/Source", "Change",
    "Proposed Rank", "Comments",
]
XLSX_COLUMN_LETTERS = [chr(65 + i) for i in range(26)] + ["A" + chr(65 + i) for i in range(26)]

# R's [:alnum:] is Unicode-aware in a UTF-8 locale: "not alphanumeric" is
# [\W_] here, and the characters each rule also allows are excluded from it.
VALUE_VALIDATION: list[dict[str, Any]] = []
for _col in [c for c in xlsx_change_srcDest_colnames if "pecies" not in c.lower()]:
    VALUE_VALIDATION.append({
        "col": _col, "regex": r"((?:[^\w\-]|_)+)", "type": "replace",
        "class": "ERROR", "code": "XLSX.NON_ALPHA-NUMERIC",
        "warn": "non-(AlphaNumeric or hyphen) characters removed", "replace": "",
    })
for _col in ["srcSpecies", "species"]:
    VALUE_VALIDATION.append({
        "col": _col, "regex": r"((?:[^\w \-]|_)+)", "type": "replace",
        "class": "ERROR", "code": "XLSX.NON_ALPHA-NUMERIC-SPACE",
        "warn": "non-(AlphaNumeric or hyphen) characters removed", "replace": "",
    })
VALUE_VALIDATION.extend([
    {"col": "species", "regex": r"^((?:[^\W_]|-)+ (?:[^\W_]|-)+)$",
     "type": "required", "class": "ERROR", "code": "XLSX.SPECIES_BAD_NAME",
     "warn": "Species name must be 'genus[space]species' binomial naming", "replace": ""},
    {"col": "change", "regex": r"((?:[^\w ;\-]|_)+)", "type": "replace",
     "class": "INFO", "code": "XLSX.CHANGE_CV_REMOVE_NON_ALPHA-NUMERIC-SPACE-SEMI",
     "warn": "Change col:non-(AlphaNumeric,hyphen,space,dash) characters removed", "replace": " "},
    {"col": "Abbrev", "regex": r"([^\w;,/. \-]+)", "type": "replace",
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
    {"col": "exemplarAccession", "regex": r"([^\w:;/. \-]+)",
     "type": "replace", "class": "ERROR", "code": "XLSX.ACCESSION.REMOVED_ILLEGAL_CHARS",
     "warn": "unexpected characters in accession list: should be a semi-colon separated list, with optional colon-separated name prefixes", "replace": ""},
    # a period ending an accession, at the end of the cell or before a ";":
    # "MG029287." -> "MG029287" (2026.010M row 16). R's regex engine has no
    # lookahead, so the ";" is matched and put back with \3, as in R.
    {"col": "exemplarAccession", "regex": r"(([.]+) *(;|\Z))",
     "type": "replace", "class": "WARNING", "code": "XLSX.ACCESSION.TRAILING_PERIOD_REMOVED",
     "warn": "accession ends with a period: period removed", "replace": r"\3"},
    # each ;-separated item is one [name:]accession[(range)]: a range
    # ("A10004-A10006") or a list separated with spaces or dashes is an error.
    # \Z, not $: R's $ doesn't match before a trailing line break.
    {"col": "exemplarAccession",
     "regex": r"^ *([^:;]+: *)?[A-Za-z0-9_./]+( ?[(][0-9]+[.][0-9]+[)])? *"
              r"(; *([^:;]+: *)?[A-Za-z0-9_./]+( ?[(][0-9]+[.][0-9]+[)])? *)*;? *\Z",
     "type": "required", "class": "ERROR", "code": "XLSX.ACCESSION.NOT_A_LIST",
     "warn": "accession must be a semicolon-separated list of single accessions, with optional colon-separated name prefixes: no ranges (A10004-A10006) or lists separated by spaces or dashes", "replace": ""},
])
for _col in ["exemplarIsolate", "genomeCoverage", "molecule", "hostSource"]:
    VALUE_VALIDATION.append({
        "col": _col, "regex": "([^\\x20-\\x7E\\u00A1-\\uFFFF\\U00010000-\\U0010FFFF]+)",
        "type": "replace", "class": "ERROR", "code": "XLSX.NON_PRINTABLE_REMOVED",
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
# A proposal row gives a lineage as one cell per rank. As in R, a cell counts
# as filled when it is not NA; an empty string left by cleaning still counts.


def _filled(realmSpecies: Sequence[Any]) -> list[int]:
    return [i for i, value in enumerate(list(realmSpecies)[:15]) if not is_na(value)]


def taxon_get_name(realmSpecies: Sequence[Any]) -> Any:
    present = _filled(realmSpecies)
    return realmSpecies[present[-1]] if present else None


def taxon_get_parent_name(realmSpecies: Sequence[Any]) -> Any:
    present = _filled(realmSpecies)
    if len(present) > 1:
        return realmSpecies[present[-2]]
    if len(present) == 1:
        # only one name means the parent is the tree root
        return params.msl_name
    return None


def taxon_get_lineage(realmSpecies: Sequence[Any], masks: Sequence[Any] = ()) -> Any:
    excluded = {str(item).lower() for item in masks if not is_na(item)}
    kept = [realmSpecies[i] for i in _filled(realmSpecies)
            if xlsx_change_ranks[i].lower() not in excluded]
    return ";".join(rstr(value) for value in kept) if kept else None


def taxon_get_rank(realmSpecies: Sequence[Any]) -> Any:
    present = _filled(realmSpecies)
    if not present:
        return None
    # dbCvList[["rank"]]$name[max(which(nonNAs))+1]: the rank table starts with "tree"
    position = present[-1] + 1
    return rankNames[position] if position < len(rankNames) else None


def taxon_get_parent_lineage(realmSpecies: Sequence[Any]) -> Any:
    if not _filled(realmSpecies):
        return None
    return taxon_get_lineage(realmSpecies, [taxon_get_rank(realmSpecies)])


def diff_strings(s1: Any, s2: Any) -> str:
    """Mark the differing middle of two strings, exactly as R does.

    The common prefix and suffix are measured independently, so they may
    overlap (R: "abab" vs "ab" gives "ab[//]ab").
    """
    left, right = rstr(s1), rstr(s2)
    longest = max(len(left), len(right))
    lmatch = 0
    while lmatch < longest and left[:lmatch + 1] == right[:lmatch + 1]:
        lmatch += 1
    rmatch = 0
    while (rmatch < longest
           and left[max(0, len(left) - rmatch - 1):] == right[max(0, len(right) - rmatch - 1):]):
        rmatch += 1
    if lmatch + rmatch < 4:
        return f"[{left}//{right}]"
    left_mid = left[lmatch:max(lmatch, len(left) - rmatch)]
    right_mid = right[lmatch:max(lmatch, len(right) - rmatch)]
    suffix = left[max(0, len(left) - rmatch):] if rmatch else ""
    return f"{left[:lmatch]}[{left_mid}//{right_mid}]{suffix}"


def diff_lineages(lin1: Any, lin2: Any) -> str:
    """Rank-by-rank diff of two lineages, with R's handling of NA lineages."""
    left = [None] if is_na(lin1) else r_strsplit(rstr(lin1))
    right = [None] if is_na(lin2) else r_strsplit(rstr(lin2))
    size = max(len(left), len(right))
    left = left + [""] * (size - len(left))
    right = right + [""] * (size - len(right))
    any_na = False
    all_same = True
    parts: list[str] = []
    for a, b in zip(left, right):
        if a is None or b is None:
            any_na = True
            parts.append("NA")
        elif a == b:
            parts.append(a)
        else:
            all_same = False
            parts.append(diff_strings(a, b))
    prefix = "NA" if any_na else ("[identical]" if all_same else "")
    return prefix + ";".join(parts)


# ---------------------------------------------------------------------------
# GenBank accession parsing and reference lookup
# ---------------------------------------------------------------------------
#
# Accessions are compared one at a time and without version suffixes (R
# fbafc8e): "DNA-A: MN013784.1; DNA-B: MN013785" has the keys MN013784 and
# MN013785. Tokens without a digit ("pending", "NA") are not accessions.


def normalize_accession_value(accession: Any) -> Any:
    if is_na(accession):
        return None
    value = rstr(accession).strip()
    if value == "":
        return None
    match = re.search(r"[^:]+$", value)
    if not match:
        return None
    # a range may follow with or without a space: "AH009795(1.3949)" (R 6ac8228)
    return re.sub(r" ?[(][0-9]+[.][0-9]+[)]$", "", match.group(0), count=1).strip()


def accession_key(accession: str) -> str:
    """Version-independent accession: MN013784.1 and MN013784. match MN013784."""
    return re.sub(r"[.][0-9]*$", "", accession, count=1)


def accession_keys(accessionList: Any) -> list[str]:
    if is_na(accessionList):
        return []
    keys: list[str] = []
    for token in re.split(r"[;,]", rstr(accessionList)):
        value = normalize_accession_value(token)
        if is_na(value):
            continue
        key = accession_key(value)
        if re.search(r"[0-9]", key) and key not in keys:
            keys.append(key)
    return keys


def accession_list_overlaps(accessionLists: Sequence[Any], accession: Any) -> list[bool]:
    """For each list: TRUE if it equals accession or shares an accession key."""
    if is_na(accession):
        return [False] * len(accessionLists)
    target = rstr(accession)
    keys = accession_keys(accession)
    key_set = set(keys)
    result = []
    for value in accessionLists:
        if is_na(value):
            result.append(False)
        elif value == target:
            result.append(True)
        elif keys and any(key in value for key in keys):
            # substring pre-filter, so only a few candidate lists get parsed
            result.append(any(key in key_set for key in accession_keys(value)))
        else:
            result.append(False)
    return result


def parseGenbankAccessions(table: Table, accessionCol: str,
                           extraCols: Sequence[str]) -> dict[str, list[dict[str, Any]]]:
    """One row per accession of a ;-separated list with optional "name:" prefixes."""
    if accessionCol not in table.columns:
        raise RuntimeError("The specified accession column does not exist in the data frame.")
    missing = [column for column in extraCols if column not in table.columns]
    if missing:
        raise RuntimeError("The following extra columns are missing in the data frame: " + ", ".join(missing))
    valid: list[dict[str, Any]] = []
    invalid: list[dict[str, Any]] = []
    pattern = re.compile(r"[a-zA-Z0-9._-]*:?[^:()]+(?: ?\([0-9]+\.[0-9]+\))?")
    for source in table.rows:
        raw = source.get(accessionCol)
        if is_na(raw) or str(raw) == "":
            continue
        for index, item in enumerate(r_strsplit(str(raw)), 1):
            item = item.strip()
            row = {column: source.get(column) for column in extraCols}
            if not pattern.fullmatch(item):
                row.update({accessionCol: item, "accession_index": index})
                invalid.append(row)
                continue
            full = re.search(r"[^:]+$", item).group(0)
            range_match = re.search(r"(?<=\()[0-9]+\.[0-9]+(?=\))", full)
            segment = re.match(r"^[^:]+(?=:)", item)
            accession = re.sub(r" ?\([0-9]+\.[0-9]+\)", "", full, count=1).strip()
            row.update({"segment_name": segment.group(0).strip() if segment else None,
                        "accession": accession,
                        "accession_range": range_match.group(0) if range_match else None,
                        "accession_index": index})
            valid.append(row)
    return {"validAccessions": valid, "invalidAccessions": invalid}


def _taxnode_int(value: Any) -> Any:
    if is_na(value):
        return None
    try:
        return int(float(str(value)))
    except ValueError:
        return value


def _vmr_matches(accession: Any) -> list[tuple[int, Any, Any]]:
    matches: list[tuple[int, Any, Any]] = []
    for key in accession_keys(accession):
        matches.extend(vmrAccessions.get(key, []))
    return sorted(set(matches), key=lambda item: item[0])


def find_vmr_accession_matches(accession: Any,
                               exclude_taxnode_ids: Sequence[Any] = ()) -> list[str]:
    excluded = {_taxnode_int(value) for value in exclude_taxnode_ids if not is_na(value)}
    names: list[str] = []
    for _, taxnode_id, species in _vmr_matches(accession):
        if taxnode_id in excluded or is_na(species):
            continue
        if rstr(species) not in names:
            names.append(rstr(species))
    return names


def vmr_taxon_has_accession(accession: Any, taxnode_ids: Sequence[Any]) -> bool:
    wanted = {_taxnode_int(value) for value in taxnode_ids if not is_na(value)}
    if not wanted:
        return False
    return any(taxnode_id in wanted for _, taxnode_id, _ in _vmr_matches(accession))


def load_abolished_accessions() -> None:
    global abolishedAccessions
    abolishedAccessions = []
    filename = params.abolished_accessions_fname
    if not os.path.exists(filename) and not filename.startswith("/"):
        candidate = os.path.join(params.ref_dir, filename)
        if os.path.exists(candidate):
            filename = candidate
    if not os.path.exists(filename) and not filename.startswith("/"):
        candidate = str(Path(__file__).resolve().parent / filename)
        if os.path.exists(candidate):
            filename = candidate
    if not os.path.exists(filename):
        if params.verbose:
            print(f"# SKIP ABOLISHED ACCESSIONS LOAD: {filename}")
        return
    table = read_mysql_tsv(filename)
    if "accession" not in table.columns:
        raise RuntimeError(f"ERROR: {filename} is missing required column: accession")
    seen = set()
    for source in table.rows:
        row = {column: (None if is_na(source.get(column)) or source.get(column) in ("", "NA", "NULL")
                        else source.get(column))
               for column in ("accession", "species", "virus_name")}
        # one row per accession key, so lists and versioned values match like VMR accessions
        for key in accession_keys(row["accession"]):
            entry = (rstr(row["accession"]).strip(), row["species"], row["virus_name"], key)
            if entry in seen:
                continue
            seen.add(entry)
            abolishedAccessions.append({"accession": entry[0], "species": entry[1],
                                        "virus_name": entry[2], "accession_key": key})
    if params.verbose:
        print(f"AbolishedAccessions:  {len(abolishedAccessions)} 4  from  {filename} ")


def find_abolished_accession_matches(accession: Any) -> list[dict[str, Any]]:
    keys = set(accession_keys(accession))
    if not keys:
        return []
    result: list[dict[str, Any]] = []
    seen = set()
    for row in abolishedAccessions:
        if row["accession_key"] in keys:
            entry = (row["accession"], row["species"], row["virus_name"])
            if entry not in seen:
                seen.add(entry)
                result.append(row)
    return result


def format_abolished_accession_matches(matches: Sequence[Mapping[str, Any]]) -> str:
    labels: list[str] = []
    for row in matches:
        species = row.get("species")
        label = species if (not is_na(species) and species != "") else row.get("virus_name")
        if not is_na(label) and label != "" and rstr(label) not in labels:
            labels.append(rstr(label))
    return "; ".join(labels)


def is_pending_accession(accession: Any) -> bool:
    return not is_na(accession) and bool(re.search(r"\bpending\b", rstr(accession), re.I))


# ---------------------------------------------------------------------------
# Reference MSL construction and controlled-vocabulary loading
# ---------------------------------------------------------------------------

TAXONOMY_INTEGER_COLUMNS = ["taxnode_id", "parent_id", "tree_id", "msl_release_num",
                            "level_id", "ictv_id", "molecule_id"]


def create_new_msl(current: Table, prev_msl: int, dest_msl: int,
                   taxnode_delta: int) -> NewMsl:
    columns = list(current.columns)
    rows = [dict(row) for row in current.rows if row.get("msl_release_num") == prev_msl]
    for column in [".prev_taxnode_id", "prev_proposals", "exemplar_name", "genome_coverage",
                   ".emptyReported", ".otherLineage", ".otherLineageProposal",
                   ".otherLineageAction"]:
        if column not in columns:
            columns.append(column)
    clear = [column for column in columns if re.match(r"^(in|out)_", column)]
    for row in rows:
        row[".prev_taxnode_id"] = row.get("taxnode_id")
        row["prev_proposals"] = None
        for column in ["taxnode_id", "parent_id", "tree_id"]:
            if not is_na(row.get(column)):
                row[column] = row[column] + int(taxnode_delta)
        row["msl_release_num"] = dest_msl
        for column in clear:
            row[column] = None
        for column in [".emptyReported", ".otherLineage", ".otherLineageProposal",
                       ".otherLineageAction"]:
            row[column] = None
        if row.get("rank") == "tree":
            row["name"] = params.msl_name
            row["notes"] = params.msl_notes
    return NewMsl(columns, rows)


def load_reference() -> None:
    global oldMSLs, curMSL, newMSL, dbCvList, dbCvMapList, cvList
    global scAbbrevNameMap, vmrAccessions, rankNames, inChangeLevels

    #### load taxonomy_node ####
    taxonomy_filename = os.path.join(params.ref_dir, params.db_taxonomy_node_fname)
    taxonomy = read_mysql_tsv(taxonomy_filename)
    for row in taxonomy.rows:
        for column, value in row.items():
            # fread: unquoted fields are stripped, and "", NULL, \N are NA
            if isinstance(value, str):
                value = value.strip(" ")
                row[column] = None if value in ("", "NULL", "\\N") else value
        for column in TAXONOMY_INTEGER_COLUMNS:
            if not is_na(row.get(column)):
                row[column] = int(float(row[column]))
    # fread(key="taxnode_id") sorts the reference at load time
    taxonomy.rows.sort(key=lambda row: row.get("taxnode_id") if not is_na(row.get("taxnode_id")) else -1)
    print(f"Previous taxa: {len(taxonomy.rows)} {len(taxonomy.columns)}  from  {taxonomy_filename} ")
    if "host_source" not in taxonomy.columns:
        print("WARNING: no host_source column in taxonomy_node dump!!! (Adding)")
        taxonomy.add_column("host_source")
    # in_change is a factor in R; its sorted levels decide the numbers R writes to msl.tsv
    inChangeLevels = sorted({row["in_change"] for row in taxonomy.rows if not is_na(row.get("in_change"))})

    last_msl = max(row["msl_release_num"] for row in taxonomy.rows)
    oldMSLs = Table(taxonomy.columns, [row for row in taxonomy.rows if row["msl_release_num"] < last_msl])
    curMSL = Table(taxonomy.columns, [row for row in taxonomy.rows if row["msl_release_num"] == last_msl])
    curMSL.add_column("out_updated")
    for row in curMSL.rows:
        row["out_updated"] = False
    newMSL = create_new_msl(curMSL, last_msl, last_msl + 1, int(params.taxnode_delta))

    #### load CVs ####
    rank_filename = os.path.join(params.ref_dir, params.db_rank_fname)
    rank_cv = read_mysql_tsv(rank_filename)
    for row in rank_cv.rows:
        row["id"] = _taxnode_int(row.get("id"))
    molecule_filename = os.path.join(params.ref_dir, params.db_molecule_fname)
    molecule_cv = read_mysql_tsv(molecule_filename)
    for row in molecule_cv.rows:
        row["id"] = _taxnode_int(row.get("id"))
    rankNames = [rstr(row.get("name")) for row in rank_cv.rows]
    dbCvList = {"rank": rank_cv, "molecule": molecule_cv}
    dbCvMapList = {
        "rank": {rstr(row.get("name")): row.get("id") for row in rank_cv.rows},
        # R builds aliases for terms that differ only in spaces, and "multiple",
        # in a local copy that is discarded: only the database abbreviations map.
        "molecule": {rstr(row.get("abbrev")): row.get("id") for row in molecule_cv.rows},
    }
    if params.verbose:
        print(f"RankCV:  {len(rank_cv.rows)} {len(rank_cv.columns)}  from  {rank_filename} ")
        print(f"MoleculeCV:  {len(molecule_cv.rows)} {len(molecule_cv.columns)}  from  {molecule_filename} ")

    ##### CVs from Proposals Template #####
    template_filename = os.path.join(params.ref_dir, params.template_xlsx_fname)
    template = read_excel_raw(template_filename, params.template_xlsx_sheet)
    width = max((len(row) for row in template), default=0)
    if params.verbose:
        print(f"ProposalTemplate[ {params.template_xlsx_sheet} ]:  {len(template)} {width}  from  {template_filename} ")
    if template:
        # patch up study section, which is a map, not just a CV list
        if width > 5:
            template[0][5] = "Subcommittee Abbrev"
        if width > 6:
            template[0][6] = "Subcommittee Name"
    name_map = {
        "Genome coverage": "genomeCoverage", "Genome composition": "molecule",
        "Host Source": "hostSource", "Change": "change",
        "Proposed Rank": "rank", "Subcommittee Abbrev": "scAbbrev",
        "Subcommittee Name": "scName",
    }
    loaded: dict[str, list[Any]] = {}
    for column in range(width):
        header = template[0][column]
        cv_name = "NA" if is_na(header) else re.sub(r"[^A-Za-z]+", " ", header)
        mapped = name_map.get(cv_name)
        if not mapped:
            continue
        values: list[Any] = []
        for row in template[1:]:
            value = row[column]
            if is_na(value):
                continue
            value = value.replace(" ", " ")
            # gsub("[^[:alnum:][:punct:]]+", " ", value)
            values.append(re.sub(r"[^\w!-/:-@\[-`{-~]+", " ", value))
        # R appends NA: an empty cell is a valid term for every CV...
        loaded[mapped] = values + [None]
    for required in ["genomeCoverage", "molecule", "hostSource", "change", "rank"]:
        loaded.setdefault(required, [None])
    # ...except "change" and "rank", where NA and "please select" are removed
    for cv in ["change", "rank"]:
        loaded[cv] = [value for value in loaded[cv]
                      if not is_na(value) and re.sub(r"[^A-Za-z0-9]", "", value).lower() != "pleaseselect"]
    sc_abbrev = loaded.pop("scAbbrev", [None])
    sc_name = loaded.pop("scName", [None])
    scAbbrevNameMap = dict(zip(sc_abbrev, sc_name))

    host_filename = os.path.join(params.ref_dir, params.db_host_source_fname)
    host_cv = read_mysql_tsv(host_filename)
    host_terms = [row.get("host_source") for row in host_cv.rows]
    dbCvList["hostSource"] = host_terms
    dbCvMapList["hostSource"] = {value: value for value in host_terms if not is_na(value)}
    if params.verbose:
        print(f"HostSourceCV:  {len(host_cv.rows)} {len(host_cv.columns)}  from  {host_filename} ")
    for value in host_terms:
        if value not in loaded["hostSource"]:
            loaded["hostSource"].append(value)
    cvList = {name: loaded[name] for name in
              ["genomeCoverage", "molecule", "hostSource", "change", "rank"]}
    cvList[".change2action"] = ACTION_MAP.copy()

    #####  VMR (species_isolates) for Genbank Accessions  #####
    vmrAccessions = {}
    vmr_filename = os.path.join(params.ref_dir, params.vmr_fname)
    if os.path.exists(vmr_filename):
        vmr = read_mysql_tsv(vmr_filename)
        parsed = parseGenbankAccessions(
            vmr, "genbank_accessions",
            ["isolate_id", "taxnode_id", "species_name", "isolate_type", "isolate_names"],
        )
        invalid = parsed["invalidAccessions"]
        if invalid:
            print("".join(
                f"ERROR: VMR accession badly formated: species:{rstr(bad.get('species_name'))} "
                f"iso_id:{rstr(bad.get('isolate_id'))} value: {rstr(bad.get('genbank_accessions'))}"
                for bad in invalid))
            raise RuntimeError(f"ERROR: {vmr_filename} contains {len(invalid)} badly formatted accessions")
        for position, row in enumerate(parsed["validAccessions"]):
            vmrAccessions.setdefault(accession_key(row["accession"]), []).append(
                (position, _taxnode_int(row.get("taxnode_id")), row.get("species_name")))
    else:
        print(f'[1] "# SKIP VMR LOAD: {vmr_filename}"')


def load_reference_cache() -> bool:
    """Compatibility stub: reference caching is deferred in version one."""
    return False


def save_reference_cache() -> bool:
    """Compatibility stub: reference caching is deferred in version one."""
    return False


# ---------------------------------------------------------------------------
# Proposal discovery and input-file classification
# ---------------------------------------------------------------------------


class RStop(RuntimeError):
    """A condition where the R script stops (quit(status=1) or stop())."""


def _scan_error(row: Mapping[str, Any]) -> None:
    allErrorDf.append(dict(row))


def _repair_spaced_filenames(root: Path) -> None:
    """PYTHON_ONLY rename_spaced_filenames: rename when collision-free, else warn."""
    files = [root] if root.is_file() else sorted(p for p in root.rglob("*") if p.is_file())
    for source in files:
        if " " not in source.name:
            continue
        target = source.with_name(re.sub(r" +", "_", source.name))
        match = re.match(r"^([0-9]+\.[0-9]+[A-Z]X*)\.", source.stem)
        code = match.group(1) if match else None
        if target.exists():
            _scan_error({
                "code": code, "level": "WARNING", "error": "FILENAME_SPACES",
                "message": "filename contains a space: please replace with _ or -",
                "notes": f"Could not rename {source.name}: {target.name} already exists",
            })
            continue
        try:
            source.rename(target)
        except OSError as exc:
            _scan_error({
                "code": code, "level": "WARNING", "error": "FILENAME_SPACES",
                "message": "filename contains a space: please replace with _ or -",
                "notes": f"Could not rename {source.name}: {exc}",
            })
            continue
        _scan_error({
            "code": code, "level": "WARNING", "error": "FILENAME_SPACES_RENAMED",
            "message": "filename contains a space: replaced with _",
            "notes": f"{source.name} -> {target.name}",
            "docx": target.name if target.suffix.lower() in {".doc", ".docx"} else None,
            "xlsx": target.name if target.suffix.lower() in {".xls", ".xlsx"} else None,
        })
        if params.verbose:
            print(f"WARNING: renamed file containing spaces: {source} -> {target}")


def _list_input_files(root: str) -> list[str]:
    """list.files(root, "^[^~.].*\\.(doc|xls)x*$", recursive=TRUE, full.names=TRUE)."""
    found: list[str] = []
    for directory, subdirs, files in os.walk(root):
        subdirs[:] = [name for name in subdirs if not name.startswith(".")]
        relative = os.path.relpath(directory, root)
        for name in files:
            if re.match(r"^[^~.].*\.(doc|xls)x*$", name):
                found.append(name if relative == "." else relative.replace(os.sep, "/") + "/" + name)
    return [root + "/" + name for name in sorted(found, key=r_collate_key)]


def _select(row: Mapping[str, Any], columns: Sequence[str]) -> dict[str, Any]:
    return {column: row.get(column) for column in columns}


def scan_for_proposals() -> dict[str, dict[str, Any]]:
    """Discover proposal files and return one metadata row per proposal code."""
    global proposalsDf, proposalsDfColumns
    proposals_dir = params.proposals_dir
    mode = params.processing_mode
    ##### filename regex #####
    if mode == "final":
        filename_regex = r"^[0-9][0-9][0-9][0-9]\.[0-9][0-9][0-9][A-Z]X*\.[^ ]*"
        filename_message = "final:####[A-Z].###[A-Z][X].____"
    elif mode == "draft":
        filename_regex = r"^[0-9][0-9][0-9][0-9]\.[0-9][0-9][0-9][A-Z]X*\.[A-Za-z]+\.v[0-9]+\.[^ ]*"
        filename_message = "draft:####[A-Z].###[A-Z][X].[A-Z]+*.v#.____"
    elif mode == "validate":
        filename_regex = r"^.*"
        filename_message = "validate:*.____"
    else:
        print(f"ERROR: --mode='{mode}' is not a valid option: validate, draft, or final")
        raise RStop("")
    if params.verbose:
        print(f"(PROCESSING) MODE      : {mode} ")

    ##### scan proposal_dir #####
    if os.path.isdir(proposals_dir):
        if PYTHON_ONLY["rename_spaced_filenames"]:
            _repair_spaced_filenames(Path(proposals_dir))
        paths = _list_input_files(proposals_dir)
    elif os.path.isfile(proposals_dir):
        paths = [proposals_dir]
    else:
        _scan_error({"notes": f"Input folder='{proposals_dir}/'", "level": "ERROR",
                     "error": "PROPOSAL_DIR_NO_EXIST", "message": "path not found"})
        write_error_summary()
        print(f"# ERROR:  PROPOSAL_DIR_NO_EXIST : {proposals_dir} ")
        raise RStop("")
    if params.tmi:
        print(f"# xls|doc(x) files found: N= {len(paths)} ")
        for path in paths:
            print(f"\t {path} ")
    if not paths:
        _scan_error({"notes": f"Input folder='{proposals_dir}/'", "level": "ERROR",
                     "error": "NO_INPUT_FILES", "message": "found no .xls[x] or .doc[x] files"})
        write_error_summary()
        print(f"# ERROR: NO_INPUT_FILES found in  {proposals_dir} ")
        raise RStop("")

    ##### parse file paths/names #####
    inputs: list[dict[str, Any]] = []
    for number, docpath in enumerate(paths, 1):
        name = os.path.basename(docpath)
        basename = re.sub(r"(.*).(doc|xls)x*$", r"\1", name)
        record = {"_rowname": str(number), "docpath": docpath,
                  "path": os.path.dirname(docpath) or ".", "file": name,
                  "basename": basename, "code": None, "scAbbrev": None}
        # code may not be in early versions of filenames
        if re.match(r"^[0-9][0-9][0-9][0-9]\.[0-9][0-9][0-9][A-Z]X*\.", basename):
            record["code"] = re.sub(r"^([0-9]+\.[0-9]+[A-Z]X*).*", r"\1", basename, count=1)
            record["scAbbrev"] = re.sub(r"^[0-9][0-9][0-9][0-9]\.[0-9][0-9][0-9]([A-Z])X*$", r"\1",
                                        record["code"], count=1)
        inputs.append(record)
    input_columns = ["docpath", "path", "file", "basename", "code", "scAbbrev"]

    ##### filter filenames #####
    def filter_inputs(matcher: Any, error: str, message: str) -> None:
        nonlocal inputs
        filtered = [record for record in inputs if matcher(record["file"])]
        if filtered:
            for record in filtered:
                row = _select(record, input_columns)
                row.update({"level": "INFO", "error": error, "message": message})
                _scan_error(row)
            write_error_summary()
        inputs = [record for record in inputs if not matcher(record["file"])]

    if mode in {"draft", "final"}:
        filter_inputs(lambda name: bool(re.search(r".*\.Ud\..*", name, re.I)),
                      "IGNORE_FNAME_.Ud.", "All files with '.Ud.' in filename are ignored")
    suppl = re.compile(normalize_posix_regex(params.infile_suppl_pat), re.I)
    filter_inputs(lambda name: bool(suppl.search(name)), "IGNORE_FNAME_SUPPL",
                  f"All files with '{params.infile_suppl_pat}' in filename are ignored")
    filter_inputs(lambda name: "appendix" in name.lower(), "IGNORE_FNAME_APPENDIXL",
                  "All files with 'appendix' in filename are ignored")

    ##### SECTION scan DOCX #####
    docxs = [dict(record) for record in inputs if re.search(r"\.docx*$", record["file"], re.I)]
    for record in docxs:
        record["docxpath"], record["docx"] = record.pop("docpath"), record.pop("file")
    seen: set[Any] = set()
    dup_codes: set[Any] = set()
    for record in docxs:
        if record["code"] in seen:
            dup_codes.add(record["code"])
        seen.add(record["code"])
    if dup_codes:
        dup_rows = [dict(_select(record, ["scAbbrev", "code", "docx"]), level="ERROR",
                         error="DUPCODE.DOCX", message="duplicate proposal ID")
                    for record in docxs if record["code"] in dup_codes]
        for row in dup_rows:
            _scan_error(row)
        write_xlsx_sheet("Sheet1", ["scAbbrev", "code", "docx", "level", "error", "message"],
                         dup_rows, Path(params.out_dir) / "QC01.docx_duplicate_ids.xlsx")
        write_error_summary()
    # use codes for row names, but fall back to ordinals, if no code in filename
    for record in docxs:
        if is_na(record["code"]):
            record["code"] = record["_rowname"]
    docx_codes = [record["code"] for record in docxs]
    if len(set(docx_codes)) != len(docx_codes):
        duplicated = sorted({code for code in docx_codes if docx_codes.count(code) > 1})
        # R: rownames(docxs) = ... fails
        raise RStop("duplicate 'row.names' are not allowed (non-unique values when setting "
                    "'row.names': " + ", ".join(f"'{code}'" for code in duplicated) + ")")
    for record in docxs:
        if " " in record["docx"]:
            row = _select(record, ["scAbbrev", "code", "docx"])
            row.update({"level": "WARNING", "error": "DOCX_FILENAME_SPACES",
                        "message": "filename contains a space: please replace with _ or -",
                        "notes": record["docx"].replace(" ", "[ ]")})
            _scan_error(row)
    if mode in {"draft", "final"}:
        for record in docxs:
            if not re.search(filename_regex + ".docx$", record["docx"], re.I):
                row = _select(record, ["scAbbrev", "code", "docx"])
                row.update({"level": "ERROR", "error": "DOCX.BAD_FILENAME_FORMAT",
                            "message": f"Should be '{filename_message}.docx'",
                            "notes": "####[A-Z]=year/study_section, ###[A-Z]=index/type, [A-Z]+=status, v#=version"})
                _scan_error(row)

    ##### SECTION: scan XLSX #####
    xlsxs = [dict(record) for record in inputs if re.search(r"\.xlsx*$", record["file"], re.I)]
    for record in xlsxs:
        record["xlsxpath"], record["xlsx"] = record.pop("docpath"), record.pop("file")
    seen = set()
    dup_list: list[Any] = []
    for record in xlsxs:
        if not is_na(record["code"]):
            if record["code"] in seen:
                dup_list.append(record["code"])
            seen.add(record["code"])
    if dup_list:
        all_dups = {id(record) for record in xlsxs if record["code"] in dup_list}
        error_rows = [dict(_select(record, ["scAbbrev", "code", "xlsx"]), level="ERROR",
                           error="XLSX_DUPCODE") for record in xlsxs if id(record) in all_dups]
        dups_ok = True
        filter_rows: set[int] = set()
        for code in dup_list:
            for record in [item for item in xlsxs if item["code"] == code]:
                load_check = load_proposal(code, record["xlsxpath"])
                if load_check.get("proposalDf") is None:
                    # report the non-proposal non-suppl-named xlsx file
                    row = _select(record, ["xlsxpath", "path", "xlsx", "basename", "code", "scAbbrev"])
                    row.update({"docpath": record["xlsxpath"], "file": record["xlsx"],
                                "level": "ERROR", "error": "NON_PROPOSAL_NOT_SUPPL",
                                "message": "Note a proposal XLSX: All supplemental XLSX files should have '_Suppl.' in the name!"})
                    _scan_error(row)
                    write_error_summary()
                    all_dups.discard(id(record))
                    filter_rows.add(id(record))
            # check if there was only one parsable xlsx
            if len([item for item in xlsxs if item["code"] == code and id(item) in all_dups]) > 1:
                dups_ok = False
            for row in error_rows:
                if row["code"] == code:
                    row["message"] = "duplicate proposal ID: " + ",".join(
                        item["xlsx"] for item in xlsxs if item["code"] == code)
        for row in error_rows:
            _scan_error(row)
        write_error_summary()
        if not dups_ok:
            print("ERROR: can not proceed with duplicated proposal IDs: "
                  + ",".join(sorted(set(dup_list), key=r_collate_key)) + " ")
            print(" ".join(f"\t {item['xlsxpath']} \n" for item in xlsxs
                           if id(item) not in filter_rows), end="")
            raise RStop("1")
        xlsxs = [record for record in xlsxs if id(record) not in filter_rows]
    for record in xlsxs:
        if is_na(record["code"]):
            record["code"] = record["_rowname"]
    if mode in {"draft", "final"}:
        bad_format = [record for record in xlsxs
                      if not re.search(filename_regex + ".xlsx*$", record["xlsx"], re.I)]
        for record in bad_format:
            row = _select(record, ["scAbbrev", "code", "xlsx"])
            row.update({"level": "ERROR", "error": "XLSX.BAD_FILENAME_FORMAT",
                        "message": f"Should be '{filename_message}.xlsx'",
                        "notes": "####[A-Z]=year/study_section, ###[A-Z]=index/type, [A-Z]+=status, v#=version"})
            _scan_error(row)
        if bad_format:
            write_error_summary()
    spaced = [record for record in xlsxs if " " in record["xlsx"]]
    if spaced:
        for record in spaced:
            row = _select(record, ["scAbbrev", "code", "xlsx"])
            row.update({"level": "ERROR", "error": "XLSX_FILENAME_SPACES",
                        "message": "filename contains a space: please replace with _ or -",
                        "notes": record["xlsx"].replace(" ", "[ ]")})
            _scan_error(row)
        write_error_summary()
    if not xlsxs:
        _scan_error({"notes": f"Input folder='{proposals_dir}/'", "level": "ERROR",
                     "error": "NO_INPUT_FILES",
                     "message": "found no .xls[x] files that passed filename QC"})
        write_error_summary()
        print(f"# ERROR: NO_INPUT_FILES found in  {proposals_dir} ")
        raise RStop("")

    ##### merge XLSX list into DOCX list #####
    def merge_key(record: Mapping[str, Any]) -> tuple:
        return (record["code"], record["path"], record.get("scAbbrev"))
    merged: dict[tuple, dict[str, Any]] = {}
    for record in xlsxs:
        merged.setdefault(merge_key(record), {})["xlsx"] = record
    for record in docxs:
        merged.setdefault(merge_key(record), {})["docx"] = record
    proposals: dict[str, dict[str, Any]] = {}
    ordered_keys = sorted(merged, key=lambda key: tuple(r_collate_key(part) for part in key))
    for key in ordered_keys:
        pair = merged[key]
        xlsx, docx = pair.get("xlsx"), pair.get("docx")
        code = key[0]
        if code in proposals:
            raise RStop("duplicate 'row.names' are not allowed (non-unique value when setting "
                        f"'row.names': '{code}')")
        # prioritize the xlsx basename for zip file naming and reporting
        basename = xlsx["basename"] if xlsx else docx["basename"]
        proposals[code] = {
            "code": code, "path": key[1], "scAbbrev": key[2],
            "xlsxpath": xlsx["xlsxpath"] if xlsx else None,
            "xlsx": xlsx["xlsx"] if xlsx else None,
            "docxpath": docx["docxpath"] if docx else None,
            "docx": docx["docx"] if docx else None,
            "basename": basename,
            # strip off version, workflow status and .fix, to get final, production filename
            "cleanbase": re.sub(r"^([0-9]+\.[0-9]+[A-Z]X*)(\.[A-Z]+)(\.v[0-9]+)*(\.fix)*(\..*)$",
                                r"\1\5", basename),
        }
    # QC - missing xlsx file
    missing = [row for row in proposals.values() if is_na(row["xlsx"])]
    if missing:
        for proposal in missing:
            _scan_error({"code": proposal["code"], "docx": proposal["docx"], "xlsx": None,
                         "row": None, "level": "ERROR", "error": "XLSX.MISSING",
                         "message": "DOCX has no matching XLSX",
                         "notes": "Suggestions: contact corresponding author"})
        write_error_summary()

    # get SC names from last letter of code
    bad_abbrevs = []
    for proposal in proposals.values():
        proposal["scAbbrev"] = None
        if re.search(r"[0-9][0-9][0-9][0-9]\.[0-9][0-9][0-9]([A-Z])X*", proposal["code"]):
            proposal["scAbbrev"] = re.sub(r"[0-9][0-9][0-9][0-9]\.[0-9][0-9][0-9]([A-Z])X*", r"\1",
                                          proposal["code"])
        # the template's abbreviation list ends with NA, so a code without a
        # letter counts as known and gets the subcommittee NA
        if proposal["scAbbrev"] not in scAbbrevNameMap:
            bad_abbrevs.append(proposal)
            proposal["subcommittee"] = ("unspecified" if is_na(proposal["scAbbrev"])
                                        else f"unknown-{proposal['scAbbrev']}")
        else:
            proposal["subcommittee"] = scAbbrevNameMap[proposal["scAbbrev"]]
    if bad_abbrevs:
        if mode in {"draft", "final"}:
            known = ",".join(rstr(name) for name in scAbbrevNameMap)
            for proposal in bad_abbrevs:
                row = _select(proposal, ["scAbbrev", "code", "xlsx", "docx"])
                row.update({"level": "ERROR", "error": "CODE_BAD_SC_ABBREV",
                            "message": "Last letter of CODE not a valid ICTV Subcommittee letter",
                            "notes": f"'{rstr(proposal['scAbbrev'])}' not in [{known}]"})
                _scan_error(row)
        write_error_summary()
    proposalsDf = proposals
    proposalsDfColumns = ["code", "path", "scAbbrev", "xlsxpath", "xlsx", "docxpath", "docx",
                          "basename", "cleanbase", "subcommittee"]
    return proposalsDf


# ---------------------------------------------------------------------------
# XLS/XLSX proposal loading
# ---------------------------------------------------------------------------

PROPOSAL_NA_STRINGS = ("Please select", "[Please select]", "[Please select]")


def load_proposal(code: str, xlsxpath_override: Any = None) -> dict[str, Any]:
    """Load one proposal workbook: 41 text columns, row i is Excel row i+1."""
    print(f"# LOAD_PROPOSAL( {code} , {rstr(xlsxpath_override)} )")
    xlsxpath = _proposal_value(code, "xlsxpath") if is_na(xlsxpath_override) else xlsxpath_override
    try:
        sheet_names = excel_sheets(xlsxpath)
    except Exception as exc:  # R stops here; Python reports the file and goes on
        if is_na(xlsxpath_override):
            log_error(code, 0, "OPEN_XLSX", "", "", "ERROR", "XLSX.OPEN_FAILED",
                      "XLS file could not be opened", str(exc), actionOrder)
        return {"proposalDf": None}
    sheets2022 = ["Proposals Template", "Menu Items (Do not change)"]
    sheets2023 = ["Instructions", "Proposal Template", "Menu Items (Do not change)"]
    extra: list[str] = []
    if all(name in sheet_names for name in sheets2022):
        extra = [name for name in sheet_names if name not in sheets2022]
    elif all(name in sheet_names for name in sheets2023):
        extra = [name for name in sheet_names if name not in sheets2023]
    if extra:
        print(f"code {code}: Extra worksheets found in xlsx: '" + "','".join(extra) + "'")
        if is_na(xlsxpath_override):
            log_error(code, 0, "OPEN_XLSX", "", "", "WARNING", "XLSX_EXTRA_SHEETS",
                      "XLS file has additional sheets not present in template",
                      "Worksheet(s) named '" + "','".join(extra) + "' were added ", actionOrder)
    regex = "proposal*"
    proposal_sheets = [name for name in sheet_names if re.search(regex, name, re.I)]
    if not proposal_sheets:
        print(f"code {code}: A worksheet name matching '{regex}' was not found: logging error")
        if is_na(xlsxpath_override):
            log_error(code, 0, "OPEN_XLSX", "", "", "ERROR", "XLSX_NOT_PROPOSAL",
                      "XLS no acceptably named 'proposal' sheet",
                      f"A worksheet name matching '{regex}' was not found, only: '"
                      + "','".join(sheet_names) + "'", actionOrder)
        return {"proposalDf": None}
    if len(proposal_sheets) > 1:
        print(f"code {code}: ERROR: More than 1 worksheet named like '{regex}' found: '"
              + "','".join(proposal_sheets) + "'")
        if is_na(xlsxpath_override):
            log_error(code, 0, "OPEN_XLSX", "", "", "ERROR", "XLSX_MULTI_PROPOSAL",
                      "XLS file has more than one sheet named 'proposal'",
                      f"More than 1 worksheet named '{regex}' were found: '"
                      + "','".join(proposal_sheets) + "'", actionOrder)
        return {"proposalDf": None}
    try:
        rows = read_excel_raw(xlsxpath, proposal_sheets[0], PROPOSAL_NA_STRINGS, columns=41)
    except Exception as exc:
        if is_na(xlsxpath_override):
            log_error(code, 0, "OPEN_XLSX", "", "", "ERROR", "XLSX.OPEN_FAILED",
                      "XLS file could not be opened", str(exc), actionOrder)
        return {"proposalDf": None}
    return {"proposalDf": rows}


# ---------------------------------------------------------------------------
# DOCX proposal metadata loading
# ---------------------------------------------------------------------------

_W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"


def _docx_root(path: str | Path) -> ElementTree.Element:
    with zipfile.ZipFile(path) as archive:
        return ElementTree.fromstring(archive.read("word/document.xml"))


def read_docx(path: str | Path) -> list[str]:
    """qdapTools::read_docx: the text of every non-empty paragraph."""
    lines = ["".join(paragraph.itertext()) for paragraph in _docx_root(path).iter(_W + "p")]
    return [line for line in lines if line != ""]


def read_proposal_docx_tables(docxpath: str | Path) -> dict[str, str]:
    """Metadata of the table-based 2024-2026 forms (R ff66cad).

    Cells are read rather than flattened paragraphs, so blank cells and
    multi-paragraph values keep their place.
    """
    root = _docx_root(docxpath)

    def cell_text(cell: ElementTree.Element) -> str:
        paragraphs = ["".join(p.itertext()) for p in cell.findall(_W + "p")]
        paragraphs = [re.sub(r"[\t\r\n]+", " ", p).strip(" \t\r\n") for p in paragraphs]
        return "; ".join(p for p in paragraphs if p)

    meta = {"title": "", "authorsEmails": "", "correspondingAuthor": "", "abstract": ""}
    authors: list[str] = []
    corresponding: list[str] = []
    taxonomy_abstracts: list[str] = []
    general_abstracts: list[str] = []
    for table in root.iter(_W + "tbl"):
        rows = [[cell_text(cell) for cell in row.findall(_W + "tc")]
                for row in table.findall(_W + "tr")]
        for i, cells in enumerate(rows):
            if not cells:
                continue
            last = i == len(rows) - 1
            if re.search(r"^Title *:$", cells[0], re.I) and len(cells) > 1:
                meta["title"] = "; ".join(cell for cell in cells[1:] if cell)

            # 2024 uses Name; 2025/2026 use Given name and Surname.
            name_cols = [j for j, cell in enumerate(cells)
                         if re.search(r"^(Name$|Given name\b|Surname$)", cell, re.I)]
            email_cols = [j for j, cell in enumerate(cells) if re.search(r"^Email address", cell, re.I)]
            corr_cols = [j for j, cell in enumerate(cells)
                         if re.search(r"^(Corresponding author|Corr\. author)", cell, re.I)]
            if name_cols and len(email_cols) == 1 and len(corr_cols) == 1 and not last:
                for author_cells in rows[i + 1:]:
                    if len(author_cells) <= max(name_cols + email_cols + corr_cols):
                        continue
                    name = " ".join(author_cells[j] for j in name_cols if author_cells[j])
                    email = author_cells[email_cols[0]]
                    if not name and not email:
                        continue
                    authors.append(" ".join(part for part in (name, email) if part))
                    if author_cells[corr_cols[0]]:
                        corresponding.append(name if name else email)

            # Some forms retain an empty General Proposal abstract before the
            # completed Taxonomy Proposal abstract. Prefer the taxonomy section.
            is_taxonomy = bool(re.search(r"^Abstract of Taxonomy Proposal", cells[0], re.I))
            is_general = bool(re.search(r"^Abstract for General Proposal", cells[0], re.I))
            if (is_taxonomy or is_general) and not last:
                content = [cell for row in rows[i + 1:] for cell in row]
                abstract = "; ".join(cell for cell in content if cell)
                if abstract:
                    (taxonomy_abstracts if is_taxonomy else general_abstracts).append(abstract)
    meta["authorsEmails"] = "; ".join(authors)
    meta["correspondingAuthor"] = "; ".join(corresponding)
    meta["abstract"] = "; ".join(taxonomy_abstracts if taxonomy_abstracts else general_abstracts)
    return meta


def load_proposal_docx(code: str) -> dict[str, Any]:
    print(f"# LOAD_PROPOSAL_DOCX( {code} )")
    meta: dict[str, Any] = {"code": code}
    path = _proposal_value(code, "docxpath")
    if is_na(path):
        if params.processing_mode in {"draft", "final"}:
            log_error(code, 0, "OPEN_DOCX", "", "", "ERROR", "DOCX_MISSING",
                      ".docx file not available", "", actionOrder)
        return meta
    extension_match = re.search(r"\.([A-Za-z0-9]+)$", str(path))
    extension = extension_match.group(1) if extension_match else ""
    if extension != "docx":
        log_error(code, 0, "OPEN_DOCX", "", "", "WARNING", "DOCX_BAD_EXT",
                  f"File type .{extension} not (yet) supported. Only .docx",
                  f"DOCX_FILENAME={rstr(_proposal_value(code, 'docx'))}", actionOrder)
        return meta
    try:
        lines = read_docx(path)
    except Exception as exc:  # R stops here; Python reports the file and goes on
        log_error(code, 0, "OPEN_DOCX", "", "", "ERROR", "DOCX_OPEN_FAILED",
                  ".docx file could not be read", str(exc), actionOrder)
        return meta

    # Table-based forms replaced the older Part 1 / Short title layout.
    if any(re.search(r"^ *Taxonomy Proposal Form, (2024|2025|2026) *$", line, re.I) for line in lines):
        fields = read_proposal_docx_tables(path)
        missing_codes = {"title": "DOCX_TITLE_MISSING", "authorsEmails": "DOCX_AUTHORS_MISSING",
                         "correspondingAuthor": "DOCX_CORR_AUTHOR_MISSING",
                         "abstract": "DOCX_ABSTRACT_MISSING"}
        missing_messages = {"title": "Title not found", "authorsEmails": "Authors list not found",
                            "correspondingAuthor": "Corresponding Author not found",
                            "abstract": "Abstract not found"}
        for field, value in fields.items():
            if value:
                meta[field] = value
            else:
                log_error(code, 0, "OPEN_DOCX", "", "", "WARNING", missing_codes[field],
                          missing_messages[field],
                          f"Field not found in DOCX metadata tables: {field}", actionOrder)
        return meta

    # check original template version
    if not any(re.search("Part 1: TITLE, AUTHORS, APPROVALS, etc", line, re.I) for line in lines):
        log_error(code, 0, "OPEN_DOCX", "", "", "INFO", "DOCX_UNKNOWN_TEMPLATE",
                  "Unrecognized version of DOCX proposal template. Not (yet) supported.",
                  f"DOCX_FILENAME={rstr(_proposal_value(code, 'docx'))}", actionOrder)
        return meta

    def clean(value: str) -> str:
        return re.sub(r"\t+", " ", re.sub(r"\n+", "; ", value))

    def following(pattern: str, key: str, missing_code: str, missing_message: str) -> None:
        indices = [i for i, line in enumerate(lines) if re.search(pattern, line, re.I)]
        if not indices:
            log_error(code, 0, "OPEN_DOCX", "", "", "WARNING", missing_code, missing_message,
                      f"Line containing expression '{pattern}' not found", actionOrder)
            return
        if key == "title":
            parts = re.split(r"[sS]hort title *: *", lines[indices[0]])
            meta[key] = clean(parts[1]) if len(parts) > 1 else None
        else:
            meta[key] = clean(lines[indices[0] + 1]) if indices[0] + 1 < len(lines) else None

    following(r"^Short title:", "title", "DOCX_TITLE_MISSING", "Title not found")
    following(r"^Author.* and email address.*", "authorsEmails", "DOCX_AUTHORS_MISSING",
              "Authors list not found")
    following(r"^Corresponding author", "correspondingAuthor", "DOCX_CORR_AUTHOR_MISSING",
              "Corresponding Author not found")
    following(r"^Abstract", "abstract", "DOCX_ABSTRACT_MISSING", "Abstract not found")
    return meta


# ---------------------------------------------------------------------------
# Proposal row normalization and QC
# ---------------------------------------------------------------------------


def _row_cell(rows: list[list[Any]], row_number: int, column: int) -> Any:
    """proposalDf[row_number, column+1]; NA outside the sheet."""
    if 1 <= row_number <= len(rows) and 0 <= column < len(rows[row_number - 1]):
        return rows[row_number - 1][column]
    return None


def _old_template_row(rows: list[list[Any]], row_number: int, code: str, first: int,
                      expected: Sequence[Any]) -> tuple[int, list[int], list[Any]]:
    """Compare a sheet row with a v1/v2 header row the way R's == does.

    proposalDf has 42 columns in R (A:AO plus "code"); positions past the
    sheet compare as NA. Returns the match count and the mismatch positions.

    Header cells such as "Exemplar \\r\\nvirus name" are stored with a CRLF.
    readxl returns it as stored, but openpyxl's XML parser turns it into LF,
    so line breaks are compared as LF on both sides.
    """
    observed = [(_row_cell(rows, row_number, column) if column < 41 else code)
                for column in range(first, min(len(expected) + first, 42))]
    matches = 0
    mismatches: list[int] = []
    for position, (actual, wanted) in enumerate(zip(observed, expected)):
        if is_na(actual) or is_na(wanted):
            continue
        if actual.replace("\r\n", "\n") == wanted.replace("\r\n", "\n"):
            matches += 1
        else:
            mismatches.append(position)
    return matches, mismatches, observed


def _template_mismatch_text(prefix: str, mismatches: Sequence[int], shown: Sequence[Any],
                            expected: Sequence[Any]) -> str:
    if mismatches:
        items = [f"column {XLSX_COLUMN_LETTERS[i]}='{rstr(shown[i] if i < len(shown) else None)}' "
                 f"instead of '{rstr(expected[i] if i < len(expected) else None)}'"
                 for i in mismatches]
    else:
        items = ["column ='' instead of ''"]
    return f"{prefix}  " + "; ".join(items)


def qc_modern_header_code(code: str, value: Any) -> None:
    """Check E1 identity and explicit version; descriptions/status are optional.

    Modern templates often calculate just the proposal ID. Do not require a
    version in those headers, or treat a descriptive suffix as part of the ID.
    This is advisory: header metadata must not suppress valid change rows.
    """
    if is_na(value):
        return
    header = rstr(value).strip()
    if header.startswith("Code"):
        if params["show.xlsx.code_miss"]:
            log_error(code, 1, "OPEN_XLSX", "", "", "INFO", "XLSX.CODE_MISS",
                      "XLSX code missing", f"XLSX cell E1 is '{header}'; expected '{code}'",
                      actionOrder)
        return
    identity = re.match(r"^(\d{4}\.\d{3}[A-Z]X*)(?=$|[._\s])", header)
    filename = rstr(_proposal_value(code, "xlsx"))
    version_pattern = r"(?:^|[._])v(\d+)(?=$|[._\s])"
    header_version = re.search(version_pattern, header)
    file_version = re.search(version_pattern, filename)
    wrong_id = identity is None or identity.group(1) != code
    wrong_version = (header_version is not None and file_version is not None
                     and int(header_version.group(1)) != int(file_version.group(1)))
    if wrong_id or wrong_version:
        expected = code + (f" (version v{file_version.group(1)})" if file_version else "")
        log_error(code, 1, "OPEN_XLSX", "", "", "WARNING", "XLSX.CODE_BAD",
                  "XLSX header code/version does not match filename",
                  f"XLSX cell E1 is '{header}'; expected {expected}; filename='{filename}'",
                  actionOrder)


def qc_proposal(code: str, proposalDf: list[list[Any]]) -> dict[str, Any]:
    """Validate and normalize every change row from one proposal workbook."""
    templateVersion = "error"
    nrow = len(proposalDf)
    version_cell = _row_cell(proposalDf, 2, 0)

    #### guess template version ####
    if not is_na(version_cell) and version_cell[:13] in (
            "version 2023.", "version 2024.", "version 2026."):
        templateVersionX = version_cell[8:]
        templateVersion = templateVersionX[:5]
        if templateVersion in ("2024.", "2026."):
            # These verified releases use the same layout as the 2023.x series (R 214291b).
            if templateVersionX in XLSX_2023_STYLE_EXACT_VERSIONS:
                templateVersion = "2023."
            else:
                log_error(code, 2, "OPEN_XLSX", "", "", "ERROR", "XLSX.TEMPLATE_UNK",
                          "XLSX template version", f"Cell A2='{version_cell}'", actionOrder)
                return {}

        def header_mismatches(row_number: int, expected: Sequence[Any]) -> list[int]:
            result = []
            for i, wanted in enumerate(expected):
                if is_na(wanted):
                    continue
                actual = re.sub(r"[\W_]", "", rstr(_row_cell(proposalDf, row_number, i))).lower()
                if actual != re.sub(r"[\W_]", "", wanted).lower():
                    result.append(i)
            return result

        def header_error(row_number: int, mismatches: Sequence[int], expected: Sequence[Any]) -> str:
            # R shows the cell of row 2, not of the mismatching row
            return f"row{row_number} mismatches template at:  " + "; ".join(
                f"column {XLSX_COLUMN_LETTERS[i]}='{rstr(_row_cell(proposalDf, 2, i))}' "
                f"instead of '{rstr(expected[i])}'" for i in mismatches)

        row3 = header_mismatches(3, XLSX_2023_ROW3)
        row4 = header_mismatches(4, XLSX_2023_ROW4)
        if row3 or row4:
            if row3:
                log_error(code, 2, "OPEN_XLSX", "", "", "ERROR", "XLSX.TEMPLATE_LINE3",
                          "XLSX template version", header_error(3, row3, XLSX_2023_ROW3), actionOrder)
            if row4:
                log_error(code, 2, "OPEN_XLSX", "", "", "ERROR", "XLSX.TEMPLATE_LINE4",
                          "XLSX template version", header_error(4, row4, XLSX_2023_ROW4), actionOrder)
            return {}
        if PYTHON_ONLY["modern_header_code"]:
            qc_modern_header_code(code, _row_cell(proposalDf, 1, 4))
    else:
        # figure out v1 or v2 templates from contents of rows 2 & 3
        v1_count, v1_bad, v1_seen = _old_template_row(proposalDf, 2, code, 0, XLSX_V1_ROW2)
        v2_count, v2_bad, _ = _old_template_row(proposalDf, 2, code, 1, XLSX_V2_ROW2)
        line2_error = ""
        if v1_count == sum(not is_na(value) for value in XLSX_V1_ROW2):
            line2_version = "v1"
        elif v2_count == sum(not is_na(value) for value in XLSX_V2_ROW2):
            line2_version = "v2"
        else:
            line2_version = "unrecognized"
            if v1_count > v2_count:
                line2_error = _template_mismatch_text("similar to v1, but with", v1_bad,
                                                      v1_seen, XLSX_V1_ROW2)
            if v1_count < v2_count:
                # R indexes row 2 and the v1 header by the v2 positions
                line2_error = _template_mismatch_text("similar to v2, but with", v2_bad,
                                                      v1_seen, XLSX_V1_ROW2)
        v1_count, v1_bad, v1_seen = _old_template_row(proposalDf, 3, code, 0, XLSX_V1_ROW3)
        v2_count, v2_bad, v2_seen = _old_template_row(proposalDf, 3, code, 0, XLSX_V2_ROW3)
        line3_error = ""
        if v1_count == sum(not is_na(value) for value in XLSX_V1_ROW3):
            line3_version = "v1"
        elif v2_count == sum(not is_na(value) for value in XLSX_V2_ROW3):
            line3_version = "v2"
        else:
            line3_version = "unrecognized"
            if v1_count > v2_count:
                line3_error = _template_mismatch_text("similar to v1, but with", v1_bad,
                                                      v1_seen, XLSX_V1_ROW3)
            if v1_count < v2_count:
                line3_error = _template_mismatch_text("similar to v2, but with", v2_bad,
                                                      v2_seen, XLSX_V2_ROW3)
        if "unrecognized" in (line2_version, line3_version):
            log_error(code, 3, "OPEN_XLSX", "", "", "ERROR", "XLSX.TEMPLATE_UNK",
                      "XLSX template version",
                      f"ROW2 is {line2_version} ({line2_error}), ROW3 is {line3_version} ({line3_error})",
                      actionOrder)
            return {}
        templateVersion = line2_version
        if templateVersion == "v1":
            # R reads params$templateURL, which does not exist, so no URL is shown
            log_error(code, 3, "OPEN_XLSX", "", "", "INFO", "XLSX.OLD_TEMPLATE_V1",
                      "XLSX template version",
                      f"You are using version {templateVersion}. Please get the latest version from ",
                      actionOrder)
        # QC Proposal ID
        if templateVersion == "v1":
            code_value, code_cell, code_row = _row_cell(proposalDf, 1, 0), "A1", 1
        else:
            code_value, code_cell, code_row = _row_cell(proposalDf, 3, 0), "A3", 3
        if not is_na(code_value) and code_value != code:
            if code_value.startswith("Code"):
                if params["show.xlsx.code_miss"]:
                    log_error(code, code_row, "OPEN_XLSX", "", "", "INFO", "XLSX.CODE_MISS",
                              "XLSX code missing",
                              f"XLSX cell {code_cell} is '{code_value}'; replace with the actual code: '{code}'",
                              actionOrder)
            else:
                log_error(code, code_row, "OPEN_XLSX", "", "", "ERROR", "XLSX.CODE_BAD",
                          "XLSX code wrong",
                          f"XLSX cell {code_cell} does not match proposal code from filename: "
                          f"'{code_value}' should be '{code}' ", actionOrder)
    if params.verbose:
        print(f"      {code}  XLSX template  {templateVersion} ")

    # extract & standardize cols & rows with change data
    if templateVersion == "v1":
        firstDataRow, selected = 4, XLSX_V1_CHANGE_COLS
    elif templateVersion == "v2":
        firstDataRow, selected = 4, XLSX_V2_CHANGE_COLS
    else:
        firstDataRow, selected = 5, XLSX_2023_CHANGE_COLS
    # R's firstDataRow:nrow(proposalDf) counts down when the sheet is shorter
    row_numbers = (range(firstDataRow, nrow + 1) if nrow >= firstDataRow
                   else range(firstDataRow, nrow - 1, -1))
    changeDf: list[dict[str, Any]] = []
    for row_number in row_numbers:
        row = {name: _row_cell(proposalDf, row_number, column)
               for name, column in zip(xlsx_change_colnames, selected)}
        # find non-empty rows
        if all(is_na(row[name]) for name in xlsx_change_srcDest_colnames):
            continue
        row[".noErrors"] = True
        row[".errors"] = None
        row[".code"] = code
        row[".linenum"] = str(row_number)
        changeDf.append(row)
    if not changeDf:
        log_error(code, firstDataRow, "OPEN_XLSX", "", "", "ERROR", "XLSX.EMPTY",
                  "XLSX no change rows found", "", actionOrder)
        return {}

    #### extract src/dest taxon names for error reporting ####
    # As in R this is done before the cells are cleaned and not repeated, so
    # these helper columns hold the text as read from the workbook.
    for row in changeDf:
        src = [row[name] for name in xlsx_change_src_colnames]
        dest = [row[name] for name in xlsx_change_dest_colnames]
        row[".srcTaxon"] = taxon_get_name(src)
        row[".srcRank"] = taxon_get_rank(src)
        row[".srcLineage"] = taxon_get_lineage(src)
        row[".destTaxon"] = taxon_get_name(dest)
        row[".destRank"] = taxon_get_rank(dest)
        row[".destLineage"] = taxon_get_lineage(dest)
        row[".destParentName"] = taxon_get_parent_name(dest)
        row[".destParentLineage"] = taxon_get_parent_lineage(dest)
        row[".changeTaxon"] = row[".srcTaxon"] if not is_na(row[".srcTaxon"]) else row[".destTaxon"]

    def log_rows(rows: list[dict[str, Any]], level: str, error: str, message: str,
                 notes: list[str]) -> None:
        if rows:
            log_error(code, [row[".linenum"] for row in rows], [row["change"] for row in rows],
                      [row["rank"] for row in rows], [row[".changeTaxon"] for row in rows],
                      level, error, message, notes, actionOrder)

    #### QC spaces, quotes, etc ####
    silent_fixes = [
        (r"0xCA", " "), (" ", " "), ("([‒–—―]+)", "-"),
        ("([“”„‟]+)", '"'), ("([‘’‚‛]+)", "'"),
        (r"([\r\n]+)", ";"), (r"^([ \t]+)", ""), (r"([ \t]+)$", ""), (r"( [ ]+)", " "),
    ]
    quotes = re.compile(r'(["]+)')
    for column in xlsx_change_colnames:
        for pattern, replacement in silent_fixes:
            for row in changeDf:
                if not is_na(row[column]):
                    row[column] = re.sub(pattern, replacement, row[column])
        if column != "comments":
            matches = [row for row in changeDf if not is_na(row[column]) and '"' in row[column]]
            if matches:
                if params.verbose:
                    print(f"INFO: {code} has {len(matches)} cells with quote in column {column} ")
                log_rows(matches, "INFO", "XLSX.QUOTES_REMOVED", "XLSX has quote",
                         [column + ":" + quotes.sub(r"[\1]", row[column]) + " (replacing with '')"
                          for row in matches])
                for row in changeDf:
                    if not is_na(row[column]):
                        row[column] = quotes.sub("", row[column])

    # check regex's for each column
    for validation in VALUE_VALIDATION:
        column = validation["col"]
        pattern = re.compile(validation["regex"])
        if validation["type"] == "replace":
            matches = [row for row in changeDf if not is_na(row[column]) and pattern.search(row[column])]
            if matches:
                if params.verbose:
                    print(f"{validation['class']}:{validation['code']} has {len(matches)} cells with "
                          f"{validation['warn']} in column '{column}'")
                log_rows(matches, validation["class"], validation["code"],
                         "XLSX has " + validation["warn"],
                         [column + " : " + pattern.sub(r"[\1]", row[column]) + " " for row in matches])
                for row in changeDf:
                    if not is_na(row[column]):
                        row[column] = pattern.sub(validation["replace"], row[column])
        else:
            matches = [row for row in changeDf if not is_na(row[column]) and not pattern.search(row[column])]
            if matches:
                log_rows(matches, validation["class"], validation["code"], validation["warn"],
                         [column + ":" + row[column] for row in matches])

    #### QC controlled vocabularies ####
    # convert "please select" to NA
    if PYTHON_ONLY["please_select_per_cell"]:
        for row in changeDf:
            for column in xlsx_change_colnames:
                if not is_na(row[column]) and re.search(r"please.*select", row[column], re.I):
                    row[column] = None
    else:
        # R tests each column as one deparsed string, so a match blanks the column
        for column in list(changeDf[0].keys()):
            joined = ", ".join(rstr(row.get(column)) for row in changeDf)
            if re.search(r"please.*select", joined, re.I | re.S):
                for row in changeDf:
                    row[column] = None

    def cv_key(value: Any) -> Any:
        return None if is_na(value) else re.sub(r"[^A-Za-z;+-]", "", value).lower()

    for cv, terms in cvList.items():
        if cv.startswith("."):
            continue
        term_keys = [cv_key(term) for term in terms]
        perfect = [row[cv] in terms for row in changeDf]
        close = [cv_key(row[cv]) in term_keys for row in changeDf]
        corrections: dict[int, Any] = {}
        for i, row in enumerate(changeDf):
            if perfect[i] or not close[i]:
                continue
            found = [j for j, key in enumerate(term_keys)
                     if key is not None and key == cv_key(row[cv])]
            if len(found) > 1:
                close[i] = False  # matching multiple things solves nothing
            elif found:
                corrections[i] = terms[found[0]]
        bad_terms = [row for i, row in enumerate(changeDf) if not close[i]]
        if bad_terms:
            valid = ",".join(rstr(term) for term in terms)
            if params.verbose:
                print(f"ERROR: {code} illegal term in CV column ' {cv} ': [ "
                      + ",".join(f"{row['.linenum']}:{rstr(row[cv])}" for row in bad_terms)
                      + " ] on rows ( " + ",".join(row[".linenum"] for row in bad_terms) + " )")
            notes = [f"XLSX incorrect value [ {rstr(row[cv])} ]. Valid terms: [ {valid} ]"
                     for row in bad_terms]
            for row in bad_terms:
                row[".errors"] = f"CV ' {cv} ' incorrect value [ {rstr(row[cv])} ]. Valid terms: [ {valid} ]"
            log_rows(bad_terms, "ERROR", "XLSX.INVALID_TERM", f"XLSX incorrect term in column {cv}", notes)
            for row in bad_terms:
                if cv in ("genomeCoverage", "molecule", "hostSource"):
                    row[cv] = None  # optional CV: drop the value
                else:
                    row[".noErrors"] = False  # obligatory CV: skip the line
            # R: next -- terms with only spacing/caps problems are left
            # uncorrected in this column when any term in it is unknown
            continue
        for i, row in enumerate(changeDf):
            if close[i] and not perfect[i]:
                if params.tmi:
                    log_rows([row], "TMI", "XLSX.TYPO_TERM",
                             f"fixed term with typo (space,caps) in column  {cv}",
                             [f"Term '{rstr(row[cv])}' replaced with '{rstr(corrections.get(i))}'. "
                              f"Valid terms: [{','.join(rstr(term) for term in terms)}]"])
                row[cv] = corrections.get(i)

    # clean change -> action
    for row in changeDf:
        row[".action"] = None if is_na(row["change"]) else ACTION_MAP.get(row["change"].lower())
    bad_rows = [row for row in changeDf if is_na(row[".action"])]
    if bad_rows:
        valid = ",".join(ACTION_MAP)
        for row in bad_rows:
            row[".noErrors"] = False
        log_rows(bad_rows, "ERROR", "ACTION.UNK", "XLSX incorrect term in column Change",
                 [f"XLSX incorrect value [ {rstr(row['change'])} ]. Valid terms: [ {valid} ]"
                  for row in bad_rows])
        # remove bad rows from further processing
        changeDf = [row for row in changeDf if not is_na(row[".action"])]

    # A species without a genome coverage is not dropped here: apply_changes
    # logs the soft error CREATE.SPECIES_NO_GENOME_COVERAGE and still creates
    # it, as R does (882245b).

    #### QC suffixes ####
    rank_rows = dbCvList["rank"].rows
    suffix_columns = [column for column in dbCvList["rank"].columns if "suffix" in column]

    def ends_with(name: str, suffix: Any) -> bool:
        return re.search(rstr(suffix) + "$", name, re.I) is not None

    missing_suffix: list[dict[str, Any]] = []
    wrong_suffix: list[tuple[dict[str, Any], str]] = []
    for row in changeDf:
        if row["change"] == "Abolish" or is_na(row[".destRank"]) or is_na(row[".destTaxon"]):
            continue
        cv_row = next((item for item in rank_rows if item.get("name") == row[".destRank"]), None)
        if cv_row is None:
            continue
        if not is_na(cv_row.get("suffix")):
            # this rank has suffixes defined: one of them must be used
            if not any(not is_na(cv_row.get(column)) and ends_with(row[".destTaxon"], cv_row[column])
                       for column in suffix_columns):
                missing_suffix.append(row)
        else:
            # this rank has no suffixes (species): no other rank's suffix may be used
            suffix_rank = None
            for item in rank_rows:
                suffixes = [item.get(column) for column in suffix_columns]
                if any(not is_na(suffix) for suffix in suffixes):
                    if any(ends_with(row[".destTaxon"], suffix) for suffix in suffixes):
                        suffix_rank = rstr(item.get("name"))
            if suffix_rank is not None:
                wrong_suffix.append((row, suffix_rank))
    if missing_suffix:
        notes = []
        for row in missing_suffix:
            cv_row = next(item for item in rank_rows if item.get("name") == row[".destRank"])
            valid = ",".join(rstr(cv_row[column]) for column in suffix_columns if not is_na(cv_row.get(column)))
            notes.append(f"New name:  {row['.destTaxon']} ; Valid suffixes for rank [ {row['.destRank']} ]:  {valid}")
        log_rows(missing_suffix, "ERROR", "MISSING_RANK_SUFFIX",
                 "Taxon name does not end with a valid suffix for it's rank", notes)
    if wrong_suffix:
        notes = []
        for row, suffix_rank in wrong_suffix:
            # every suffix of the other ranks, column by column
            invalid = ",".join(rstr(item[column]) for column in suffix_columns for item in rank_rows
                               if item.get("name") != row[".destRank"] and not is_na(item.get(column)))
            notes.append(f"New name: '{row['.destTaxon']}' has the suffix of a '{suffix_rank}' "
                         f"but is a '{row['.destRank']}'; *IN*valid suffixes for rank "
                         f"[{row['.destRank']}]: {invalid}")
        log_rows([row for row, _ in wrong_suffix], "ERROR", "SUFFIX_RANK_MISMATCH",
                 "Taxon name does not end with a valid suffix for it's rank", notes)

    #### QC capitalization ####
    bad_case = [row for row in changeDf
                if row["change"] != "Abolish" and not is_na(row[".destRank"])
                and not is_na(row[".destTaxon"]) and row[".destRank"] != "species"
                and not re.search(r"^[A-Z][a-z0-9-]+$", row[".destTaxon"])]
    log_rows(bad_case, "ERROR", "TAXON_NAME_CASE", "Non-species taxa must be Capitalized.",
             [f"{row['.destTaxon']}: Non-species taxa must start with a capital, followed by one "
              "or more lower-case letters, numbers or hyphens." for row in bad_case])
    return {"changeDf": changeDf}


# ---------------------------------------------------------------------------
# Proposal-set loading and QC orchestration
# ---------------------------------------------------------------------------


def load_and_qc_proposals(proposals: dict[str, dict[str, Any]],
                          existing: dict[str, list[dict[str, Any]]]) -> dict[str, list[dict[str, Any]]]:
    result = dict(existing)
    codes = list(proposals)
    for position, code in enumerate(codes, 1):
        if code in result:
            continue
        if is_na(proposals[code].get("xlsx")):
            print(f"# SKIPPED:  {code}  (no .xlsx)")
            continue
        # try to load proposal metadata from DOCX
        meta = load_proposal_docx(code)
        # R counts every row logged so far, not only this proposal's
        print(f"# LOADED:  {code}  DOCX with  {len(allErrorDf)}  errors/warnings")
        for column, value in meta.items():
            if column not in proposalsDfColumns:
                proposalsDfColumns.append(column)
            proposalsDf[code][column] = value
        # try to load proposal from XLSX
        raw = load_proposal(code).get("proposalDf")
        print(f"# LOADED:  {code}  XLS with  {len(allErrorDf)}  errors/warnings")
        if len(allErrorDf) > 0:
            write_error_summary()
        print(f"# LOADED:  {code} ")
        if raw is None:
            print(f"SKIP:  {code} : proposal could not be loaded")
            continue
        print(f"# QC start: {code}: proposal loaded ({position} out of {len(codes)})")
        checked = qc_proposal(code, raw)
        print(f"# QCed:      {code}  with  {len(allErrorDf)}  errors/warnings")
        if "changeDf" in checked:
            result[code] = checked["changeDf"]
    print("changeList:  " + " ".join(result) + " ")
    write_error_summary()
    return result


# ---------------------------------------------------------------------------
# Merge the proposals' changes, sort them, and order chained changes
# ---------------------------------------------------------------------------
#
# Abolish and merge run first, from species up; then create/split, then
# rename/move, from the root down. Changes linked through a taxon (A->B, then
# B->C) are then re-ordered along their chain.

ACTION_ORDER = {
    # destructive: run first, species->root by rank
    "abolish": (-1, 0), "merge": (-1, 0),
    # constructive: root -> species by rank; first create (w/in a rank)
    "new": (1, 1), "promote": (1, 1), "demote": (1, 1), "split": (1, 1),
    # then move (w/in a rank)
    "rename": (1, 2), "move": (1, 2),
}


def merge_and_order_changes(changes: Mapping[str, list[dict[str, Any]]]) -> list[dict[str, Any]]:
    """Return every change without QC errors, in the order they are applied."""
    global allChangeDf
    skipped = merged = 0
    codes = list(proposalsDf)
    rows: list[dict[str, Any]] = []
    for position, code in enumerate(codes, 1):
        changeDf = changes.get(code)
        if changeDf is None:
            skipped += 1
            continue
        clean = [dict(row) for row in changeDf if row.get(".noErrors") is True]
        print(f"# MERGE PROC: {code} with {len(clean)} changes  ({position} out of {len(codes)})")
        rows.extend(clean)
        merged += 1
    for row in rows:
        row[".codeLine"] = f"{row['.code']}:{row['.linenum']}"
    print(f"# DONE. MERGED {merged}  proposal(s), skipped {skipped} ; total changes: {len(rows)} ")
    allChangeDf = rows
    if not rows:
        return []

    ##### sort merged actions #####
    def sort_key(item: tuple[int, dict[str, Any]]) -> tuple:
        index, row = item
        action_order, sub_rank = ACTION_ORDER.get(row.get(".action"), (None, None))
        rank_id = None if is_na(row.get("rank")) else dbCvMapList["rank"].get(row["rank"])
        rank_order = (None if (is_na(rank_id) or action_order is None)
                      else rank_id * (1 if action_order > 0 else -1))
        return (action_order is None, action_order or 0, rank_order is None, rank_order or 0,
                sub_rank is None, sub_rank or 0, index)
    order = [index for index, _ in sorted(enumerate(rows), key=sort_key)]
    for position, index in enumerate(order, 1):
        rows[index][".changeOrder"] = position

    ##### detect chaining #####
    by_line = {row[".codeLine"]: row for row in rows}
    for row in rows:
        row[".chainOrder"] = 0
        row[".chainStart"] = ""
        row[".originalOrder"] = row[".changeOrder"]
        row[".srcRankTaxon"] = paste0(row.get(".srcRank"), ":", row.get(".srcTaxon"))
        row[".destRankTaxon"] = paste0(row.get(".destRank"), ":", row.get(".destTaxon"))
    # links: prev.destTaxon = next.srcTaxon (and not self)
    sources: dict[str, list[dict[str, Any]]] = {}
    for row in rows:
        if row[".srcRankTaxon"] != "NA:NA":
            sources.setdefault(row[".srcRankTaxon"], []).append(row)
    links: list[dict[str, Any]] = []
    for prev in rows:
        key = prev[".destRankTaxon"]
        if key == "NA:NA":
            continue
        for nxt in sources.get(key, []):
            if prev[".codeLine"] == nxt[".codeLine"]:
                continue
            links.append({
                "key": key, "prev": prev[".codeLine"], "src": prev[".srcRankTaxon"],
                "changePrev": prev["change"], "orderPrev": prev[".changeOrder"],
                "chainOrderPrev": 0, "chainStartPrev": "",
                "next": nxt[".codeLine"], "destNext": nxt[".destRankTaxon"],
                "changeNext": nxt["change"], "orderNext": nxt[".changeOrder"],
                "chainOrderNext": 0, "chainStartNext": "",
            })
    links.sort(key=lambda link: r_collate_key(link["key"]))  # merge() sorts by the key; stable
    if not links:
        if params.verbose:
            print("# CHAIN: NO chained changes found")
        return [rows[index] for index in order]
    if params.verbose:
        print(f"# CHAIN: Found {len(links)} links!")

    def r_number(value: float) -> float:
        """A number after R stored it as text and read it back."""
        return float(format(value, ".15g"))

    chain_depth = 1
    # mark first link in each chain
    keys = {link["key"] for link in links}
    keys_of_real_changes = {link["key"] for link in links if link["src"] != link["key"]}
    for link in links:
        # the link's source is not the DEST of another link
        simple = link["src"] not in keys
        # links [A,A] (src=dest), where A is not the DEST of any other link
        no_op = link["src"] == link["key"] and link["src"] not in keys_of_real_changes
        if simple or no_op:
            link["start"] = True
            link["chainOrderPrev"] = float(f"{link['orderPrev']}.{chain_depth}")
            link["chainStartPrev"] = link["prev"]
    active = [dict(link, _ref=link) for link in links if link.get("start")]
    if params.verbose:
        print(f"# CHAIN: level  {chain_depth}  has  {len({a['prev'] for a in active})}  "
              f"startLinks in  {len(active)}  chains ")
    for snapshot in active:
        by_line[snapshot["prev"]][".chainOrder"] = snapshot["chainOrderPrev"]
        by_line[snapshot["prev"]][".chainStart"] = snapshot["chainStartPrev"]

    # get and mark next link in each chain
    while active and chain_depth < 5000:
        chain_depth += 1
        # the greater of the chain head's order and the next link's own order
        for snapshot in active:
            link = snapshot["_ref"]
            chain_order_next = int(link["chainOrderPrev"]) + chain_depth / 10
            order_next = int(link["orderNext"]) + chain_depth / 10
            link["chainOrderNext"] = chain_order_next if chain_order_next > order_next else order_next
            link["chainStartNext"] = link["chainStartPrev"]
        for snapshot in active:
            link = snapshot["_ref"]
            by_line[link["next"]][".chainOrder"] = r_number(link["chainOrderNext"])
            by_line[link["next"]][".chainStart"] = link["chainStartNext"]
        # copy outgoing (Next) .chainOrder to incoming (Prev) side of any following links
        for snapshot in active:
            link = snapshot["_ref"]
            for onward in links:
                if onward["prev"] == link["next"]:
                    onward["chainOrderPrev"] = r_number(link["chainOrderNext"])
                    onward["chainStartPrev"] = link["chainStartNext"]
        # find next links in this chain: they connect to the Next of the active
        # links, and their chain was not on the active links already
        next_lines = {snapshot["next"] for snapshot in active}
        seen_starts = {snapshot["chainStartNext"] for snapshot in active}
        active = [dict(link, _ref=link) for link in links
                  if link["prev"] in next_lines and link["chainStartPrev"] not in seen_starts]
        if params.verbose:
            print(f"# CHAIN: level  {chain_depth}  has  {len({a['prev'] for a in active})}  "
                  f"startLinks in  {len(active)}  chains ")

    # record final ordering
    def chained_order(item: tuple[int, dict[str, Any]]) -> tuple:
        index, row = item
        value = row[".changeOrder"] if row[".chainOrder"] == 0 else row[".chainOrder"]
        return (float(value), index)
    order = [index for index, _ in sorted(enumerate(rows), key=chained_order)]
    for position, index in enumerate(order, 1):
        rows[index][".changeOrder"] = position

    # errors for any links that didn't get traversed
    bad_links = [link for link in links if link["chainOrderPrev"] == 0 or link["chainOrderNext"] == 0]
    for link in bad_links:
        prev = by_line[link["prev"]]
        # R passes prevTaxon$actionOrder, a column that does not exist: order is NA
        allErrorDf.append({
            "subcommittee": _proposal_value(prev[".code"], "subcommittee"),
            "code": prev[".code"], "row": _as_numeric(prev[".linenum"]), "order": None,
            "change": prev["change"], "rank": prev["rank"], "taxon": prev[".changeTaxon"],
            "docx": text(_proposal_value(prev[".code"], "docx"), "MISSING"),
            "xlsx": text(_proposal_value(prev[".code"], "xlsx"), "MISSING"),
            "level": "ERROR", "error": "CHAIN.HIDDEN",
            "message": "Chained change order not resolved",
            "notes": (f"Change links with  {link['next']} :  {link['key']}  "
                      f"[ {rstr(link['changeNext'])} ]  {link['destNext']}"),
        })
    if bad_links:
        write_error_summary()
    return [rows[index] for index in order]


# ---------------------------------------------------------------------------
# Applying the changes to the working MSL
# ---------------------------------------------------------------------------


def update_lineage(parent_id: Any, parent_lineage: Any, parent_otherLineage: Any,
                   parent_otherLineageProposal: Any, parent_otherLineageAction: Any) -> int:
    """Set the lineage of every descendant of parent_id, and how it changed."""
    count = 0
    stack = [(parent_id, parent_lineage, parent_otherLineage)]
    while stack:
        node_id, lineage, other = stack.pop()
        for row in newMSL.with_id(node_id):
            row[".otherLineage"] = other
            for column, value in ((".otherLineageProposal", parent_otherLineageProposal),
                                  (".otherLineageAction", parent_otherLineageAction)):
                row[column] = value if is_na(row.get(column)) else paste0(row[column], ";", value)
        kids = newMSL.kids(node_id)
        for kid in kids:
            kid["lineage"] = paste0(lineage, ";", kid.get("name"))
            if not is_na(other):
                kid[".otherLineage"] = paste0(other, ";", kid.get("name"))
        count += len(kids)
        # depth-first, in row order, like R's recursion
        for kid in reversed(kids):
            stack.append((kid.get("taxnode_id"), kid.get("lineage"), kid.get(".otherLineage")))
    return count


def newmsl_lineage_from_parents(taxnode_id: Any) -> str:
    """Lineage of a newMSL taxon, built from parent_id up to the tree root.

    The lineage column is not always current: MOVE (and PROMOTE/DEMOTE) update
    a taxon's parent_id but not its own lineage.
    """
    names: list[str] = []
    current = taxnode_id
    for _ in range(50):
        rows = newMSL.with_id(current)
        # stop at the tree root (level 100), or if the parent is missing
        if len(rows) != 1 or is_na(rows[0].get("level_id")) or rows[0]["level_id"] <= 100:
            break
        names.insert(0, rstr(rows[0].get("name")))
        parent = rows[0].get("parent_id")
        if is_na(parent) or parent == current:
            break
        current = parent
    return ";".join(names)


def newmsl_resolve_merged_id(taxnode_id: Any) -> Any:
    """taxnode_id if still in newMSL, else the taxon its merges lead to, else NA."""
    seen: list[Any] = []
    current = taxnode_id
    while not is_na(current) and current not in seen:
        if newMSL.with_id(current):
            return current
        seen.append(current)
        merged = [row for row in curMSL.rows
                  if row.get("out_change") == "merge" and row.get(".merge_src_taxnode_id") == current]
        if len(merged) != 1:
            break
        current = merged[0].get(".out_taxnode_id")
    return None


def _append_filename(old: Any, proposalZip: str) -> str:
    return proposalZip if is_na(old) else paste0(old, ";", proposalZip)


_cur_by_id: dict[Any, list[dict[str, Any]]] = {}


def _cur_rows(taxnode_id: Any) -> list[dict[str, Any]]:
    return [] if is_na(taxnode_id) else _cur_by_id.get(taxnode_id, [])


def abolish_taxon(curChangeDf: dict[str, Any], srcNewTarget: dict[str, Any]) -> None:
    """Remove a taxon from newMSL, and record the abolish on curMSL."""
    proposalZip = paste0(_proposal_value(curChangeDf[".code"], "basename"), ".zip")
    src = [curChangeDf.get(name) for name in xlsx_change_src_colnames]
    dest = [curChangeDf.get(name) for name in xlsx_change_dest_colnames]
    previous = _cur_rows(srcNewTarget.get(".prev_taxnode_id"))
    newMSL.remove(srcNewTarget)
    for row in previous:
        row["out_updated"] = True
        row["out_change"] = "abolish"
        row["out_filename"] = _append_filename(row.get("out_filename"), proposalZip)
        row["out_target"] = taxon_get_name(dest)
        row["out_notes"] = paste0("linenum=", curChangeDf[".linenum"])
        row[".actionOrder"] = curChangeDf.get(".order")
    log_change_error(curChangeDf, "SUCCESS", "ABOLISH.OK", "Change=ABOLISH, applied successfully",
                     paste0("ABOLISH ", taxon_get_rank(src), " named ", taxon_get_name(src)))


def _find_genus(parent: Mapping[str, Any]) -> tuple[str, Mapping[str, Any]]:
    """The genus a species under `parent` belongs to: ("ok"|"subgenus"|"none", taxon)."""
    if parent.get("rank") == "genus":
        return "ok", parent
    if parent.get("rank") == "subgenus":
        grand = newMSL.with_id(parent.get("parent_id"))
        if grand and grand[0].get("rank") == "genus":
            return "ok", grand[0]
        return "subgenus", grand[0] if grand else {}
    return "none", parent


def apply_changes(changesDf: list[dict[str, Any]]) -> None:
    """Apply every change, in order, to newMSL and record it on curMSL."""
    global actionOrder, _cur_by_id
    # Successful rename/move outcomes keyed by proposal change row. A later
    # rename can replace the name again, so final names do not prove whether
    # an earlier operation succeeded on the merge's intended target.
    appliedRenameMoves: dict[str, dict[str, Any]] = {}
    # genera to scan for binomial issues when we're done
    renamedGenera: list[dict[str, Any]] = []
    changes_by_row = {row[".codeLine"]: row for row in changesDf}

    #### admin columns ####
    for table in (curMSL, newMSL):
        for row in table.rows:
            row[".split"] = False
            row[".split_kept"] = False
            for column in (".split_code", ".split_linenum", ".split_actionOrder", ".split_acc_used"):
                row[column] = None
    # ABOLISH of a taxon with kids is deferred until all changes are applied
    for row in newMSL.rows:
        row[".abolish"] = False
        row[".abolish_row"] = None
        row[".abolish_actionOrder"] = None

    # the current and older MSLs do not change: index their names once
    cur_by_name: dict[Any, list[dict[str, Any]]] = {}
    cur_by_lower: dict[Any, list[dict[str, Any]]] = {}
    _cur_by_id = {}
    for row in curMSL.rows:
        if not is_na(row.get("name")):
            cur_by_name.setdefault(row["name"], []).append(row)
            cur_by_lower.setdefault(row["name"].lower(), []).append(row)
        _cur_by_id.setdefault(row.get("taxnode_id"), []).append(row)
    old_by_name: dict[Any, list[dict[str, Any]]] = {}
    old_by_lower: dict[Any, list[dict[str, Any]]] = {}
    for row in oldMSLs.rows:
        if not is_na(row.get("name")):
            old_by_name.setdefault(row["name"], []).append(row)
            old_by_lower.setdefault(row["name"].lower(), []).append(row)

    def newest(rows: list[dict[str, Any]]) -> dict[str, Any]:
        """matchDf[order(msl_release_num, decreasing=T)[1],]"""
        best = rows[0]
        for row in rows[1:]:
            if row.get("msl_release_num") > best.get("msl_release_num"):
                best = row
        return best

    def case_only(index: Mapping[Any, list[dict[str, Any]]], destName: Any, excluded: Any) -> list[dict[str, Any]]:
        # rows whose name equals destName apart from capitalization; only renames
        # may exclude their source, because they preserve the taxon's identity
        return [row for row in index.get(str(destName).lower(), [])
                if is_na(excluded) or row.get("name") != excluded]

    def rank_id(rank: Any) -> Any:
        return dbCvMapList["rank"].get(rank) if not is_na(rank) else None

    def molecule_id(molecule: Any) -> Any:
        return None if is_na(molecule) else dbCvMapList["molecule"].get(molecule)

    def binomial_ok(species: Any, genus: Mapping[str, Any]) -> bool:
        return paste0(genus.get("name"), " ") in rstr(species)

    def prior_change_notes(rows: list[dict[str, Any]]) -> list[str]:
        return [paste0("proposal(s)=", row.get("out_filename"), " also did a '", row.get("out_change"),
                       "' to '", row.get("out_target"), "'") for row in rows]

    def other_lineage_notes(rows: list[dict[str, Any]]) -> list[str]:
        notes = [paste0("proposal(s)=", row.get(".otherLineageProposal"), " did a '",
                        row.get(".otherLineageAction"), "' ") for row in rows]
        return notes or ["proposal(s)= did a '' "]

    def kid_names(kids: list[dict[str, Any]]) -> list[str]:
        return [paste0(kid.get("rank"), ":", kid.get("name")) for kid in kids]

    def root_rows() -> list[dict[str, Any]]:
        return [row for row in newMSL.rows if row.get("taxnode_id") == row.get("tree_id")]

    def duplicate_accession_rows(accession: Any, skip: Any = None) -> list[dict[str, Any]]:
        if is_na(accession):
            return []
        flags = accession_list_overlaps([row.get("genbank_accession_csv") for row in newMSL.rows],
                                        accession)
        return [row for row, flag in zip(newMSL.rows, flags) if flag and row is not skip]

    prevDestParent: list[dict[str, Any]] | None = None  # as in R, kept between changes
    code = None
    order = 0
    #### FOR(row) ####
    for change in changesDf:
        if change.get(".noErrors") is not True:
            continue
        ##### change definition #####
        order += 1
        row_key = change[".codeLine"]
        cur = dict(change)
        action = cur.get("change")
        code = cur[".code"]
        linenum = cur[".linenum"]
        proposalZip = paste0(_proposal_value(code, "basename"), ".zip")
        cur[".order"] = order
        act = cur.get(".action")
        ACT = rstr(act).upper()
        CHANGE = rstr(action).upper()

        rankClean = None
        if is_na(cur.get("rank")) or cur["rank"].lower() not in rankNames:
            # (R stops here: the call that logs RANK.UNK has a wrong argument name)
            log_change_error(cur, "ERROR", "RANK.UNK", "Rank name unknown",
                             paste0(", rank=", cur.get("rank"), "; knownRanks=", ",".join(rankNames[1:])))
        else:
            rankClean = cur["rank"].lower()

        srcRealmSpecies = [cur.get(name) for name in xlsx_change_src_colnames]
        destRealmSpecies = [cur.get(name) for name in xlsx_change_dest_colnames]
        srcTaxonName = taxon_get_name(srcRealmSpecies)
        srcTaxonRank = taxon_get_rank(srcRealmSpecies)
        srcLineage = taxon_get_lineage(srcRealmSpecies)
        destTaxonName = taxon_get_name(destRealmSpecies)
        destTaxonRank = taxon_get_rank(destRealmSpecies)
        destParentName = taxon_get_parent_name(destRealmSpecies)
        destLineage = taxon_get_lineage(destRealmSpecies)
        destParentLineage = taxon_get_lineage(destRealmSpecies, [destTaxonRank])

        # check where this name may or may not exist already
        srcNewMatches = newMSL.named(srcTaxonName)

        # check that at least one of dest/src taxa were given
        if is_na(srcTaxonName) and is_na(destTaxonName):
            log_change_error(cur, "ERROR", "XLSX.NO_SRC_DEST",
                             f"Change={CHANGE}, but neither 'current taxonomy', nor 'proposed taxonomy', were specified")
            continue
        # check if src was already abolished
        # (ABOLISH of a taxon with kids is deferred, so it's still in newMSL)
        srcAbolished = [row for row in srcNewMatches if row.get(".abolish")]
        if srcAbolished:
            abolishChange = changes_by_row[srcAbolished[0][".abolish_row"]]
            log_change_error(cur, "ERROR", f"{ACT}.SRC_ALREADY_CHANGED",
                             f"Change={CHANGE}, but taxon already abolished",
                             paste0("proposal(s)=", _proposal_value(abolishChange[".code"], "basename"),
                                    ".zip line ", abolishChange[".linenum"],
                                    " also did a 'abolish' to '", srcTaxonName, "'"))
            continue
        # check if lineage not correct in "current"
        if len(srcNewMatches) == 1:
            srcNewTaxon = srcNewMatches[0]
            if cur.get(".srcLineage") != srcNewTaxon.get("lineage"):
                if not is_na(srcNewTaxon.get(".otherLineage")):
                    # this is a known change, done by another proposal
                    log_change_error(cur, "INFO", "SRC.LINEAGE_CHANGED",
                                     "CURRENT lineage changed by another proposal",
                                     paste0("proposal(s)=", srcNewTaxon.get(".otherLineageProposal"),
                                            " did a ", srcNewTaxon.get(".otherLineageAction"),
                                            " PROPOSAL//NEW_MSL=",
                                            diff_lineages(cur.get(".srcLineage"), srcNewTaxon.get("lineage"))))
                else:
                    # unexplained lineage mismatch
                    log_change_error(cur, "WARNING", "SRC.LINEAGE_WRONG",
                                     "'current taxonomy' in proposal doesn't match current MSL. Typo?",
                                     paste0("MSL//PROPOSED=",
                                            diff_lineages(cur.get(".srcLineage"), cur.get(".destLineage")),
                                            ", PROPOSAL=", cur.get(".srcLineage"),
                                            "; MSL=", srcNewTaxon.get("lineage")))

        # make sure new taxon does not exist already
        #   move/merge are not checked; split can be to the same name, if src=dest
        if (not is_na(destTaxonName) and act not in ("move", "merge")
                and not (act == "split" and not is_na(srcTaxonName) and destTaxonName == srcTaxonName)):
            def existing_notes(match: Mapping[str, Any]) -> str:
                return paste0("proposed ", destTaxonRank, "=", destTaxonName, "; existing ",
                              match.get("rank"), "=", "MSL", match.get("msl_release_num"), ":",
                              match.get("lineage"))
            destCurMatches = cur_by_name.get(destTaxonName, [])
            destOldMatches = old_by_name.get(destTaxonName, [])
            destNewMatches = newMSL.named(destTaxonName)
            # check if exists in the current MSL
            if destCurMatches:
                log_change_error(cur, "ERROR", "DEST.IN_CUR",
                                 f"Change={CHANGE}, but taxon name already exists",
                                 existing_notes(newest(destCurMatches)))
                continue
            # check if exists in a historical MSL (not the most recent)
            if destOldMatches:
                log_change_error(cur, "ERROR", "DEST.IN_OLD",
                                 f"Change={CHANGE},  but taxon name existed historically",
                                 existing_notes(newest(destOldMatches)))
                continue
            # check if already created in newMSL
            if destNewMatches:
                match = newest(destNewMatches)
                log_change_error(cur, "ERROR", "DEST.IN_NEW",
                                 f"Change={CHANGE}, but taxon name already created in new MSL",
                                 paste0(existing_notes(match), "; otherProposal=", match.get("in_filename")))
                continue
            # Reject capitalization collisions before applying the change (R 8824d74).
            # A split can introduce a new taxon, so its source must also be checked.
            excluded = srcTaxonName if act == "rename" else None
            destCurCaseMatches = case_only(cur_by_lower, destTaxonName, excluded)
            destOldCaseMatches = case_only(old_by_lower, destTaxonName, excluded)
            destNewCaseMatches = case_only(newMSL.by_lower, destTaxonName, excluded)
            if destCurCaseMatches:
                log_change_error(cur, "ERROR", "DEST.IN_CUR_CASE",
                                 f"Change={CHANGE}, but taxon name already exists with different capitalization",
                                 existing_notes(newest(destCurCaseMatches)))
                continue
            if destOldCaseMatches:
                log_change_error(cur, "ERROR", "DEST.IN_OLD_CASE",
                                 f"Change={CHANGE}, but taxon name existed historically with different capitalization",
                                 existing_notes(newest(destOldCaseMatches)))
                continue
            if destNewCaseMatches:
                match = newest(destNewCaseMatches)
                log_change_error(cur, "ERROR", "DEST.IN_NEW_CASE",
                                 f"Change={CHANGE}, but taxon name already created in new MSL with different capitalization",
                                 paste0(existing_notes(match), "; otherProposal=", match.get("in_filename")))
                continue
            # check if taxon rank matches change rank
            if rstr(destTaxonRank).lower() != rstr(cur.get("rank")).lower():
                log_change_error(cur, "ERROR", "DEST.RANK_MISMATCH",
                                 f"Change={CHANGE} but proposed taxon rank does not match [rank] column",
                                 paste0("rankColumn=", cur.get("rank"), ", proposedTaxonRank=", destTaxonRank,
                                        ", proposedTaxonomy=", destLineage))

        if act == "new" or (act == "split" and (is_na(srcTaxonName) or is_na(destTaxonName)
                                                or srcTaxonName != destTaxonName)):
            #### CREATE/SPLIT ####
            # SPLIT (srcName != destName); split (same name) is handled by "move" code.
            # If no split line keeps the original name by the end, the original is
            # removed from newMSL:
            #   .split_kept = T means the MOVE code saw a srcName = destName split
            #   .split = T means this code saw a srcName != destName split
            srcCurTaxa: list[dict[str, Any]] = []
            if act == "new":
                # should have NO source taxon
                if not is_na(srcTaxonName):
                    # report on the proposed taxon: this may really be a split
                    cur[".changeTaxon"] = cur.get(".destTaxon")
                    log_change_error(cur, "ERROR", "CREATE.W_SRC",
                                     "Change=CREATE, but 'current taxonomy' is not empty; perhaps you meant SPLIT",
                                     paste0("CURRENT//PROPOSED=",
                                            diff_lineages(cur.get(".srcLineage"), cur.get(".destLineage"))))
            else:
                if is_na(srcTaxonName):
                    log_change_error(cur, "ERROR", "SPLIT.NO_SRC",
                                     "Change=SPLIT, but 'current taxonomy' columns are empty",
                                     paste0("proposedTaxonomy=", cur.get(".destLineage")))
                    continue
                if is_na(destTaxonName):
                    log_change_error(cur, "ERROR", "SPLIT.NO_DEST",
                                     "Change=SPLIT, but 'proposed taxonomy' columns are empty",
                                     paste0("currentTaxonomy=", cur.get(".srcLineage")))
                    continue
                # has a source taxon, let's find it
                srcCurTaxa = cur_by_name.get(srcTaxonName, [])
                if not srcCurTaxa:
                    log_change_error(cur, "ERROR", "SPLIT.SRC_NO_EXIST",
                                     "Change=SPLIT, but 'current taxonomy' doesn't exist",
                                     paste0("currentTaxonomy=", cur.get(".srcLineage")))
                    continue
                # remember this taxon has been split - may need to obsolete it,
                # later, if no split directive keeps the old name
                for row in srcCurTaxa:
                    row[".split"] = True
                    row[".split_code"] = code
                    row[".split_linenum"] = linenum
                    row[".split_actionOrder"] = order
                    row[".actionOrder"] = order
                # mark the newMSL entry we may delete, if not already marked
                if not srcCurTaxa[0].get(".split_kept"):
                    srcNewTaxa = [row for src in srcCurTaxa for row in newMSL.from_prev(src.get("taxnode_id"))]
                    if srcNewTaxa:
                        for row in srcNewTaxa:
                            row[".split"] = True
                            row[".split_code"] = code
                            row[".split_linenum"] = linenum
                            row[".split_actionOrder"] = order
                    else:
                        log_change_error(cur, "ERROR", "SPLIT.SRC_NOT_IN_NEW",
                                         "Change=SPLIT, but 'current taxonomy' can not be found in newMSL",
                                         paste0("currentTaxonomy=", cur.get(".srcLineage")))
                        continue

            # check DEST taxon specified
            if is_na(destTaxonName):
                log_change_error(cur, "ERROR", "CREATE.NO_DEST",
                                 f"Change={ACT}, but 'proposed taxonomy' columns are empty",
                                 paste0("CURRENT//PROPOSED=",
                                        diff_lineages(cur.get(".srcLineage"), cur.get(".destLineage")),
                                        ",CURRENT=", cur.get(".srcLineage"),
                                        ", PROPOSED=", cur.get(".destLineage")))
                continue

            accession = cur.get("exemplarAccession")
            # check if accession number is pending
            if is_pending_accession(accession):
                log_change_error(cur, "ERROR", "CREATE.PENDING_ACC",
                                 f"Change={ACT}, accession number is pending",
                                 paste0("accession=", accession))
                continue
            # check if accession number was used in a previous VMR but is no longer current
            abolishedMatches = find_abolished_accession_matches(accession)
            if abolishedMatches:
                log_change_error(cur, "ERROR", "CREATE.ABOLISHED_ACC",
                                 f"Change={ACT}, accession number was used previously and is no longer a current VMR accession",
                                 paste0("accession=", accession, ", previousSpecies=",
                                        format_abolished_accession_matches(abolishedMatches)))
                continue

            # check if same accession number already exists
            splitReusesSourceAccession = False
            if act == "split" and srcCurTaxa and any(row.get(".split") for row in srcCurTaxa):
                splitReusesSourceAccession = (
                    any(accession_list_overlaps([row.get("genbank_accession_csv") for row in srcCurTaxa],
                                                accession))
                    or vmr_taxon_has_accession(accession, [row.get("taxnode_id") for row in srcCurTaxa]))
            dupRows = duplicate_accession_rows(accession)
            vmrDupAccessionTaxa = [] if dupRows else find_vmr_accession_matches(accession)
            if dupRows or vmrDupAccessionTaxa:
                if splitReusesSourceAccession:
                    # a split may re-use an accession under a new or same name, but just once
                    if not is_na(srcCurTaxa[0].get(".split_acc_used")):
                        log_change_error(cur, "ERROR", "SPLIT.REUSE_ACC",
                                         f"Change={ACT}, accession re-used more than once in a split",
                                         paste0("accession=", accession, "; already reused on ",
                                                srcCurTaxa[0].get(".split_acc_used")))
                        continue
                    for row in srcCurTaxa:
                        row[".split_acc_used"] = paste0(code, ":", linenum)
                    if params.tmi:
                        log_change_error(cur, "INFO", "SPLIT.REUSE_ACC",
                                         f"Change={ACT}, accession re-used in a split",
                                         paste0("accession=", accession, "; currentName=", srcTaxonName,
                                                "; proposedName=", destTaxonName))
                else:
                    existing = ("; ".join(rstr(row.get("lineage")) for row in dupRows) if dupRows
                                else "; ".join(vmrDupAccessionTaxa))
                    log_change_error(cur, "ERROR", "CREATE.DUP_ACC",
                                     f"Change={ACT}, a species with this accession number already exists",
                                     paste0("accession=", accession, ", existingSpecies=", existing))
                    # hard error in FINAL mode
                    if params.processing_mode == "final":
                        continue

            # verify that PARENT taxon exists already in newMSL
            if is_na(destParentName):
                parents = root_rows()
            else:
                parents = newMSL.named(destParentName)
                if not parents:
                    # look for the name in the current MSL, and assume it was renamed
                    prevDestParent = cur_by_name.get(destParentName, [])
                    if len(prevDestParent) == 1:
                        parents = newMSL.with_id(prevDestParent[0].get(".out_taxnode_id"))
                        log_change_error(cur, "INFO", "CREATE.PARENT_ALREADY_CHANGED",
                                         "PROPOSED parent taxon already modified",
                                         other_lineage_notes(parents))
                    else:
                        # look for another change that already changed the name in newMSL
                        changedParent = [row for row in allChangeDf if row.get(".srcTaxon") == destParentName]
                        if len(changedParent) == 1:
                            log_change_error(cur, "INFO", "CREATE.PARENT_ALREADY_CHANGED_2",
                                             "PROPOSED parent taxon already modified",
                                             other_lineage_notes(newMSL.named(changedParent[0].get(".destTaxon"))))
                        else:
                            # just missing
                            log_change_error(cur, "ERROR", "CREATE.PARENT_NO_EXIST",
                                             f"Change={ACT}, but parent rank taxon does not exist",
                                             paste0("parentTaxon=", destParentName,
                                                    ", proposedTaxonomy=", destLineage))
            if params.verbose:
                label = "CREATE(SPLIT):     " if act == "split" else "CREATE:            "
                print(f"{label}{rstr(cur.get('rank')).upper()} code:{code} line:{linenum} "
                      f"'{rstr(destTaxonName)}' findParent({rstr(destParentName)})={len(parents)} ")
            # no parent found: skip record, otherwise, get parent
            if len(parents) != 1:
                continue
            destParentTaxon = parents[0]

            # check if taxon rank matches change rank
            if rstr(destTaxonRank).lower() != rstr(cur.get("rank")).lower():
                log_change_error(cur, "ERROR", "CREATE.RANK_MISMATCH",
                                 f"Change={ACT}, but proposed taxon rank does not match [rank] column",
                                 paste0("ankColumn=", cur.get("rank"), ", proposedTaxonRank=", destTaxonRank,
                                        ", proposedTaxonomy=", destLineage))

            ###### binomial check ######
            if destTaxonRank == "species":
                # check hostSource
                if is_na(cur.get("hostSource")) or re.search(r"please.*select", cur["hostSource"], re.I):
                    log_change_error(cur, "ERROR", "CREATE.SPECIES_NO_HOST_SOURCE",
                                     f"Change={ACT}, but proposed species must have a host/source value")
                # check exemplar virus name, genome coverage and genome composition (R 882245b)
                # Soft errors: the species is still created. These cells fill the new
                # species' VMR isolate. A cell with no letters or digits counts as empty.
                for column, label, error in (
                        ("exemplarName", "Exemplar virus name", "CREATE.SPECIES_NO_EXEMPLAR_NAME"),
                        ("genomeCoverage", "Genome coverage", "CREATE.SPECIES_NO_GENOME_COVERAGE"),
                        ("molecule", "Genome composition", "CREATE.SPECIES_NO_GENOME_COMPOSITION")):
                    value = cur.get(column)
                    if (is_na(value) or not re.search(r"[^\W_]", value)
                            or re.search(r"please.*select", value, re.I)):
                        log_change_error(cur, "ERROR", error,
                                         f"Change={ACT}, but proposed species is missing '{label}'")
                # check if parent is a genus
                status, genus = _find_genus(destParentTaxon)
                if status == "subgenus":
                    log_change_error(cur, "ERROR", "CREATE.SPECIES_SUBGENUS_NO_GENUS",
                                     f"Change={ACT}, can not find parent genus",
                                     paste0("subgenus '", destParentTaxon.get("name"), "' is not in a genus;",
                                            " it's parent '", genus.get("name"), "' is a ", genus.get("rank")))
                    continue
                if status == "none":
                    log_change_error(cur, "ERROR", "CREATE.SPECIES_NO_GENUS",
                                     f"Change={ACT}, can not find parent genus",
                                     paste0("parent '", destParentTaxon.get("name"), "' is a ",
                                            destParentTaxon.get("rank")))
                    continue
                # check binomial naming
                if not binomial_ok(destTaxonName, genus):
                    log_change_error(cur, "ERROR", "CREATE.SPECIES_BINOMIAL_MISMATCH",
                                     f"Change={ACT}, but proposed species names does not start with 'genus[space]' per binomial naming convention",
                                     paste0("parent genus name =", genus.get("name")))
                    continue

            #####  create new taxon #####
            # A created taxon starts as a copy of its parent's row.
            newTaxon = {key: value for key, value in destParentTaxon.items() if key != "_pos"}
            # WARN if PARENT_LINEAGE is not expected, AND USE PARENT LINEAGE
            if not is_na(destParentLineage) and destParentTaxon.get("lineage") != destParentLineage:
                if not is_na(destParentTaxon.get(".otherLineage")):
                    # (the parent itself being renamed is reported above as
                    # CREATE.PARENT_ALREADY_CHANGED; this is its ancestors changing)
                    if cur.get(".destParentName") == destParentTaxon.get("name"):
                        log_change_error(cur, "INFO", "CREATE.PARENT_RENAMED",
                                         f"Change={ACT}, proposed parent taxon exists, but not with expected name/lineage, using observed lineage",
                                         paste0("otherProposal(s)=", destParentTaxon.get(".otherLineageProposal"),
                                                " did a ", destParentTaxon.get(".otherLineageAction"),
                                                ", PROPOSED//OBSERVED=",
                                                diff_lineages(destParentLineage, destParentTaxon.get("lineage"))))
                else:
                    log_change_error(cur, "WARNING", "CREATE.PARENT_LINEAGE",
                                     f"Change={ACT}, proposed parent taxon exists, but not with expected lineage, using observed lineage",
                                     paste0("PROPOSED//OBSERVED=",
                                            diff_lineages(destParentLineage, destParentTaxon.get("lineage")),
                                            ", otherProposals=", destParentTaxon.get("prev_proposals")))
                destLineage = paste0(destParentTaxon.get("lineage"), ";", destTaxonName)

            newTaxon["taxnode_id"] = max(newMSL.by_id) + 1
            newTaxon["parent_id"] = destParentTaxon.get("taxnode_id")
            newTaxon["ictv_id"] = newTaxon["taxnode_id"]
            newTaxon["name"] = destTaxonName
            newTaxon["cleaned_name"] = destTaxonName
            newTaxon["level_id"] = rank_id(destTaxonRank)
            newTaxon["rank"] = destTaxonRank
            newTaxon["in_change"] = act
            newTaxon["in_filename"] = proposalZip
            newTaxon["in_notes"] = paste0("xlsx_row=", linenum)
            newTaxon["in_target"] = destLineage
            # clear any fields that should NOT be "inherited"
            newTaxon["out_change"] = None
            newTaxon["out_target"] = ""
            newTaxon["out_filename"] = ""
            newTaxon["out_notes"] = ""
            newTaxon["out_updated"] = False
            for column in ("is_ref", "is_offical", "is_hidden", "is_deleted", "is_typo", "is_obsolete"):
                newTaxon[column] = "0"
            newTaxon[".abolish"] = False
            # a new taxon is never the original of a split, even if its parent is (R 992902e)
            newTaxon[".split"] = False
            newTaxon[".split_kept"] = False
            for column in (".split_code", ".split_linenum", ".split_actionOrder", ".split_acc_used"):
                newTaxon[column] = None
            # propagate alternate names/lineages
            if not is_na(destParentTaxon.get(".otherLineage")):
                newTaxon[".otherLineage"] = paste0(destParentTaxon[".otherLineage"], ";", destTaxonName)
                newTaxon[".otherLineageProposal"] = destParentTaxon.get(".otherLineageProposal")
                newTaxon[".otherLineageAction"] = destParentTaxon.get(".otherLineageAction")
            newTaxon["molecule_id"] = molecule_id(cur.get("molecule"))
            newTaxon["host_source"] = cur.get("hostSource")
            newTaxon["notes"] = cur.get("comments")
            newTaxon[".actionOrder"] = order
            if act == "split":
                # a split line that is not the original name: keep it, and record the split
                newTaxon[".split_kept"] = True
                newTaxon[".split_code"] = code
                newTaxon[".split_actionOrder"] = order
                newTaxon[".split_linenum"] = linenum
                newTaxon["in_change"] = act
                newTaxon["in_filename"] = _append_filename(newTaxon["in_filename"], proposalZip)
                newTaxon["in_target"] = srcLineage
                newTaxon["in_notes"] = paste0("linenum=", linenum)
            if destTaxonRank == "species":
                for column in ("exemplarAccession", "exemplarName", "Abbrev", "exemplarIsolate",
                               "genomeCoverage"):
                    newTaxon[XLSX2DB_MAP[column]] = cur.get(column)
            newTaxon["lineage"] = destLineage
            newTaxon[".prev_taxnode_id"] = None
            newTaxon["prev_proposals"] = paste0(_proposal_value(code, "basename"), ".zip:", linenum)
            for column in newTaxon:
                newMSL.add_column(column)
            if act not in inChangeLevels:
                inChangeLevels.append(act)
            newMSL.append(newTaxon)
            log_change_error(cur, "SUCCESS", "CREATE.OK" if act == "new" else "SPLIT.OK",
                             f"Change={ACT}, applied successfully",
                             paste0("Create ", newTaxon["rank"], " of '", newTaxon["lineage"], "'"))

        elif act == "rename":
            #### RENAME ####
            if is_na(srcTaxonName):
                log_change_error(cur, "ERROR", "RENAME.WO_SRC",
                                 "Change=RENAME, but 'current taxonomy' columns are empty",
                                 paste0("currentTaxonomy=", srcLineage, ", destTaxonomy=", destLineage))
                continue
            # check if renamed taxon has the same name
            if not is_na(destTaxonName) and srcTaxonName == destTaxonName:
                log_change_error(cur, "ERROR", "RENAME.SAME_NAME",
                                 "Change=RENAME, but new name is the same; might this be a move?",
                                 paste0("currentVsProposed=", diff_lineages(srcLineage, destLineage),
                                        ", currentTaxonomy=", srcLineage, ", destTaxonomy=", destLineage))
            # find the target to rename in NEW
            targets = newMSL.named(srcTaxonName)
            if params.verbose:
                print(f"RENAME:            {rstr(cur.get('rank')).upper()} code:{code} line:{linenum} "
                      f"'{rstr(destTaxonName)}' findTarget({srcTaxonName})={len(targets)} ")
            # In the next three errors R passes the texts in the other order:
            # the message column holds the detail and notes holds the summary.
            if not targets:
                prevSrcTaxon = cur_by_name.get(srcTaxonName, [])
                if not prevSrcTaxon:
                    log_change_error(cur, "ERROR", "RENAME.NO_EXIST",
                                     paste0("taxon=", srcTaxonName, ", lineage=", srcLineage,
                                            ", proposedTaxon=", destTaxonName),
                                     "Change=RENAME, but taxon does not exist")
                else:
                    log_change_error(cur, "ERROR", "RENAME.SRC_ALREADY_CHANGED",
                                     prior_change_notes(prevSrcTaxon),
                                     "Change=RENAME, but taxon already modified")
                continue
            if len(targets) > 1:
                log_change_error(cur, "ERROR", "RENAME.MANY",
                                 paste0("taxon=", srcTaxonName, ", lineage=", srcLineage,
                                        ", proposedTaxon=", destTaxonName, " matches ", len(targets)),
                                 "Change=RENAME, multiple taxa exist with parent name")
                continue
            target = targets[0]
            srcPrevTargets = _cur_rows(target.get(".prev_taxnode_id"))
            srcNewParents = newMSL.with_id(target.get("parent_id"))
            parent: Mapping[str, Any] = srcNewParents[0] if srcNewParents else {}
            ###### binomial check ######
            if destTaxonRank == "species":
                status, genus = _find_genus(parent)
                if status == "subgenus":
                    log_change_error(cur, "ERROR", "RENAME.SPECIES_SUBGENUS_NO_GENUS",
                                     "Change=RENAME, can not find parent genus",
                                     paste0("subgenus '", parent.get("name"), "' is not in a genus;",
                                            " it's parent '", genus.get("name"), "' is a ", genus.get("rank")))
                    continue
                if status == "none":
                    log_change_error(cur, "ERROR", "CREATE.SPECIES_NO_GENUS",
                                     "Change=CREATE, can not find parent genus",
                                     paste0("parent '", parent.get("name"), "' is a ", parent.get("rank")))
                    continue
                if not binomial_ok(destTaxonName, genus):
                    log_change_error(cur, "ERROR", "RENAME.SPECIES_BINOMIAL_MISMATCH",
                                     "Change=RENAME, but proposed species names does not start with 'genus[space]' per binomial naming convention",
                                     paste0("parent genus name =", genus.get("name")))
                    continue
            if PYTHON_ONLY["rename_duplicate_accession"] and srcTaxonRank == "species":
                accession = cur.get("exemplarAccession")
                dupRows = duplicate_accession_rows(accession, skip=target)
                vmrDuplicates = find_vmr_accession_matches(
                    accession, [target.get(".prev_taxnode_id")])
                existing = list(dict.fromkeys([rstr(row.get("name")) for row in dupRows] + vmrDuplicates))
                if existing:
                    # Rename does not change exemplar metadata. Flag the sheet
                    # for correction, but preserve R's taxonomy/dependency behavior.
                    log_change_error(cur, "WARNING", "RENAME.DUP_ACC",
                                     "Change=RENAME, a species with this accession number already exists",
                                     paste0("accession=", accession, ", existingSpecies=",
                                            "; ".join(existing)))
            ##### rename apply changes  #####
            oldLineage = target.get("lineage")
            newMSL.set_name(target, destTaxonName)
            target["lineage"] = paste0(parent.get("lineage"), ";", destTaxonName)
            target[".actionOrder"] = order
            # set the lineage of the renamed taxon's kids
            update_lineage(target.get("taxnode_id"), target["lineage"], oldLineage,
                           paste0(proposalZip, ":", linenum),
                           paste0(act, " ", cur.get(".srcRank"), " ", cur.get(".srcTaxon"),
                                  " to ", cur.get(".destTaxon")))
            # add to list of genera to check for binomial problems when all changes are applied
            if destTaxonRank == "genus":
                renamedGenera.append({"name": destTaxonName, "code": code, "linenum": linenum,
                                      "actionOrder": order})
            target["prev_proposals"] = paste0(
                "" if is_na(target.get("prev_proposals")) else paste0(target["prev_proposals"], ","),
                proposalZip, ":", linenum)
            # put out_* changes on curMSL
            for row in srcPrevTargets:
                row["out_updated"] = True
                row["out_change"] = "rename"
                row["out_filename"] = _append_filename(row.get("out_filename"), proposalZip)
                row["out_target"] = destTaxonName
                row["out_notes"] = paste0("linenum=", linenum)
                row[".out_taxnode_id"] = target.get("taxnode_id")
                row[".actionOrder"] = order
            appliedRenameMoves[row_key] = {"taxnode_id": target.get("taxnode_id"),
                                           "name": rstr(target.get("name")),
                                           "rank": rstr(target.get("rank"))}
            log_change_error(cur, "SUCCESS", "RENAME.OK", "Change=RENAME, applied successfully",
                             paste0("RENAME ", srcTaxonRank, " from '", srcTaxonName, "' to '",
                                    target.get("name"), "' in '", parent.get("lineage"), "'"))

        elif act == "abolish":
            ##### ABOLISH #####
            if is_na(srcTaxonName):
                log_change_error(cur, "ERROR", "ABOLISH.WO_SRC",
                                 "Change=ABOLISH, but 'current taxonomy' columns are empty",
                                 paste0("destTaxonomy=", destLineage))
                continue
            targets = newMSL.named(srcTaxonName)
            if not targets:
                prevSrcTaxon = cur_by_name.get(srcTaxonName, [])
                if not prevSrcTaxon:
                    # not modified, just missing
                    log_change_error(cur, "ERROR", "ABOLISH.NO_EXIST",
                                     "Change=ABOLISH, but taxon does not exist",
                                     paste0("taxon=", srcTaxonName, ", lineage=", srcLineage,
                                            ", proposedTaxon=", destTaxonName))
                else:
                    # previous record - must have been modified already this MSL
                    log_change_error(cur, "ERROR", "ABOLISH.SRC_ALREADY_CHANGED",
                                     "Change=ABOLISH, but taxon already modified",
                                     prior_change_notes(prevSrcTaxon))
                continue
            if len(targets) > 1:
                log_change_error(cur, "ERROR", "ABOLISH.MANY",
                                 "Change=ABOLISH, multiple taxa exist with name",
                                 paste0("taxon=", srcTaxonName, ", lineage=", srcLineage,
                                        ", proposedTaxon=", destTaxonName, " matches ", len(targets)))
                continue
            target = targets[0]
            # check for kids - can't abolish before kids. Abolish is sorted before
            # move, so the kids may still be moved out by a later change: mark it,
            # and re-check once all changes are applied (R f0e638e)
            if newMSL.kids(target.get("taxnode_id")):
                target[".abolish"] = True
                target[".abolish_row"] = row_key
                target[".abolish_actionOrder"] = order
                continue
            abolish_taxon(cur, target)

        elif (act in ("move", "promote", "demote")
              or (act in ("split", "merge") and not is_na(srcTaxonName) and not is_na(destTaxonName)
                  and srcTaxonName == destTaxonName)):
            ##### MOVE, PROMOTE/DEMOTE, SPLIT=, MERGE= #####
            # these all require a SRC and DEST, and change structure
            if is_na(srcTaxonName) and act in ("promote", "demote"):
                # (R d6484d0: this used to stop the run)
                log_change_error(cur, "ERROR", "CHANGE_RANK.NO_SRC",
                                 f"Change={ACT}, but 'current taxonomy' columns are empty",
                                 paste0("destTaxonomy=", destLineage))
                continue
            if is_na(srcTaxonName) and act == "move":
                # they didn't specify the source, but it's a MOVE: try to find it
                guess = newMSL.named(destTaxonName)
                if len(guess) == 1:
                    log_change_error(cur, "ERROR", "MOVE.NO_SRC_CAN_GUESS",
                                     "Change=MOVE, but 'current taxonomy' columns are empty",
                                     paste0("destTaxonomy=", destLineage,
                                            "; we guess you meant=", guess[0].get("lineage")))
                    srcTaxonName = guess[0].get("name")
                    srcTaxonRank = guess[0].get("rank")
                    srcLineage = guess[0].get("lineage")
                else:
                    log_change_error(cur, "ERROR", "MOVE.NO_SRC",
                                     "Change=MOVE, but 'current taxonomy' columns are empty",
                                     paste0("destTaxonomy=", destLineage,
                                            "; no existing taxon '", destTaxonName, "' found"))
                    continue

            ##### validate rank #####
            if is_na(destTaxonRank) or destTaxonRank != rankClean:
                # 1. consistent: dest rank doesn't match action rank
                log_change_error(cur, "ERROR", "MOVE.DEST_WRONG_RANK",
                                 f"Change={ACT}, but 'proposed taxonomy' ranks doesn't match 'proposed rank' column",
                                 paste0("rank=", rankClean, ", srcTaxonomy=[", srcTaxonRank, "]", srcLineage))
                continue
            if act == "move" and srcTaxonRank != destTaxonRank:
                # 2. MOVE: src/dest ranks NOT the same
                log_change_error(cur, "ERROR", "MOVE.DIFF_RANK",
                                 "Change=MOVE, but 'current taxonomy' and 'proposed taxonomy' ranks differ; use promote or demote, not move, to change rank",
                                 paste0("CUR//PROPOSED=[", srcTaxonRank, "//", destTaxonRank, "] ",
                                        diff_lineages(srcLineage, destLineage)))
                continue
            if act in ("promote", "demote") and srcTaxonRank == destTaxonRank:
                # 3. PROMOTE/DEMOTE src/dest ranks ARE the same
                log_change_error(cur, "ERROR", "CHANGE_RANK.SAME_RANK",
                                 f"Change={ACT}, but 'current taxonomy' and 'proposed taxonomy' ranks are the same; use MOVE",
                                 paste0("rank=", rankClean, ", CUR//PROPOSED=",
                                        diff_lineages(srcLineage, destLineage)))
                continue

            # find the SRC target to MOVE in NEW and CUR
            targets = newMSL.named(srcTaxonName)
            if not targets:
                prevSrcTaxon = cur_by_name.get(srcTaxonName, [])
                if not prevSrcTaxon:
                    log_change_error(cur, "ERROR", "MOVE.NO_EXIST",
                                     f"Change={ACT}, but taxon does not exist",
                                     paste0("taxon=", srcTaxonName, ",CUR//PROPOSED=",
                                            diff_lineages(srcLineage, destLineage),
                                            ", CUR=", srcLineage, ", PROPOSED=", destLineage))
                else:
                    log_change_error(cur, "ERROR", "MOVE.SRC_ALREADY_CHANGED",
                                     f"Change={ACT}, but taxon already modified",
                                     prior_change_notes(prevSrcTaxon))
                continue
            if len(targets) > 1:
                log_change_error(cur, "ERROR", "MOVE.MANY",
                                 f"Change={ACT}, multiple taxa exist with name",
                                 paste0("taxon=", srcTaxonName, ", lineage=", srcLineage,
                                        ", proposedTaxon=", destTaxonName, " matches ", len(targets)))
                continue
            target = targets[0]
            srcPrevTaxnodeID = target.get(".prev_taxnode_id")
            srcPrevTargets = _cur_rows(srcPrevTaxnodeID)
            if params.verbose:
                print(f"MOVE:              {rstr(cur.get('rank')).upper()} code:{code} line:{linenum} "
                      f"'{srcTaxonName}' findTarget({srcTaxonName}): new={len(targets)}/cur={len(srcPrevTargets)} ")

            # taxon not linked to curMSL: non-fatal ERROR
            if is_na(srcPrevTaxnodeID):
                log_change_error(cur, "ERROR", "ALTER_NEW_TAXON",
                                 f"Change={ACT}, but taxon was created in this MSL",
                                 paste0(srcTaxonRank, ":", srcTaxonName, " was ",
                                        rstr(target.get("in_change")).upper(), " by ",
                                        target.get("prev_proposals")))

            if (PYTHON_ONLY["move_rename_name_check"] and act == "move"
                    and destTaxonName != srcTaxonName):
                others = [row for row in newMSL.named(destTaxonName) if row is not target]
                if others:
                    log_change_error(cur, "ERROR", "DEST.EXISTS",
                                     "Proposed destination taxon already exists",
                                     f"{destTaxonName} matched {len(others)} row(s)")
                    continue

            accession = cur.get("exemplarAccession")
            # check if accession number was used in a previous VMR but is no longer current
            abolishedMatches = find_abolished_accession_matches(accession)
            if abolishedMatches:
                log_change_error(cur, "ERROR", "MOVE.ABOLISHED_ACC",
                                 f"Change={ACT}, accession number was used previously and is no longer a current VMR accession",
                                 paste0("accession=", accession, ", previousSpecies=",
                                        format_abolished_accession_matches(abolishedMatches)))
                continue
            # check if same accession number already exists (not checking against ourself)
            dupRows = duplicate_accession_rows(accession, skip=target)
            vmrDupAccessionTaxa = [] if dupRows else find_vmr_accession_matches(accession, [srcPrevTaxnodeID])
            if dupRows or vmrDupAccessionTaxa:
                existing = ("; ".join(rstr(row.get("lineage")) for row in dupRows) if dupRows
                            else "; ".join(vmrDupAccessionTaxa))
                log_change_error(cur, "ERROR", "MOVE.DUP_ACC",
                                 f"Change={ACT}, a species with this accession number already exists",
                                 paste0("accession=", accession, ", existingSpecies=", existing))
                continue

            # verify that PARENT taxon exists already in newMSL
            parents = root_rows() if is_na(destParentName) else newMSL.named(destParentName)
            if params.verbose:
                print(f"MOVE:              {rstr(cur.get('rank')).upper()} code:{code} line {linenum} "
                      f"'{rstr(destTaxonName)}' findParent({rstr(destParentName)})={len(parents)} ")
            if not parents:
                # check if it used to exist
                prevDestParent = cur_by_name.get(destParentName, [])
                if not prevDestParent:
                    log_change_error(cur, "ERROR", "MOVE.PARENT_NO_EXIST",
                                     f"Change={ACT}, but parent rank taxon does not exist",
                                     paste0("parentTaxon=", destParentName, ", proposedTaxonomy=", destLineage))
                else:
                    log_change_error(cur, "ERROR", "MOVE.PARENT_ALREADY_CHANGED",
                                     f"Change={ACT}, but parent taxon already modified",
                                     prior_change_notes(prevDestParent))
                continue
            if len(parents) > 1:
                log_change_error(cur, "ERROR", "MOVE.PARENT_MANY",
                                 f"Change={ACT}, multiple taxa exist with parent name",
                                 paste0("parentTaxon=", destParentName, ", proposedTaxonomy=", destLineage))
                continue
            destParentTaxon = parents[0]

            # remember to check later if a genus name change broke binomial
            if destTaxonRank == "genus":
                renamedGenera.append({"name": destTaxonName, "code": code, "linenum": linenum,
                                      "actionOrder": order})
            # check that promote/demote change name
            if act in ("promote", "demote") and srcTaxonName == destTaxonName:
                # (R stops here: the notes use a variable that does not exist)
                log_change_error(cur, "ERROR", "CHANGE_RANK.SAME_NAME",
                                 f"Change={ACT}, but current and proposed names are the same.",
                                 paste0("current name=[", srcTaxonRank, "]", srcTaxonName,
                                        " and proposed name [", destTaxonRank, "]", destTaxonName))

            # WARN if PARENT_LINEAGE is not expected
            if not is_na(destParentLineage) and destParentTaxon.get("lineage") != destParentLineage:
                if not is_na(destParentTaxon.get(".otherLineage")):
                    # a parent of the destParent was renamed, moved, etc.
                    log_change_error(cur, "INFO", "MOVE.PROPOSED_PARENT_LINEAGE_CHANGED",
                                     "PROPOSED parent lineage already modified",
                                     paste0("proposal(s)=", destParentTaxon.get(".otherLineageProposal"),
                                            " did a ", destParentTaxon.get(".otherLineageAction"),
                                            "; OBSERVED//PROPOSED=",
                                            diff_lineages(destParentTaxon.get("lineage"), destParentLineage)))
                elif prevDestParent and prevDestParent[0].get("lineage") != destParentLineage:
                    # prevDestParent is whatever an earlier change last looked up,
                    # as in R: it is not reset for each change
                    log_change_error(cur, "WARNING", "MOVE.PARENT_LINEAGE",
                                     "PROPOSED parent taxon exists, but not with expected lineage",
                                     paste0(",PROPOSED//CUR=",
                                            diff_lineages(destParentLineage, destParentTaxon.get("lineage")),
                                            ", PROPOSED=", destParentLineage,
                                            ", CUR=", destParentTaxon.get("lineage"),
                                            ", otherProposals=", destParentTaxon.get("prev_proposals")))
                else:
                    log_change_error(cur, "INFO", "MOVE.PARENT_LINEAGE_UPDATED2",
                                     "PROPOSED parent taxon exists, but with updated lineage",
                                     paste0(",PROPOSED//CUR=",
                                            diff_lineages(destParentLineage, destParentTaxon.get("lineage")),
                                            ", PROPOSED=", destParentLineage,
                                            ", CUR=", destParentTaxon.get("lineage")))

            ##### binomial check (species) #####
            if destTaxonRank == "species":
                status, genus = _find_genus(destParentTaxon)
                if status == "subgenus":
                    log_change_error(cur, "ERROR", "MOVE.SPECIES_SUBGENUS_NO_GENUS",
                                     f"Change={ACT}, can not find parent genus",
                                     paste0("subgenus '", destParentTaxon.get("name"), "' is not in a genus;",
                                            " it's parent '", genus.get("name"), "' is a ", genus.get("rank")))
                    continue
                if status == "none":
                    log_change_error(cur, "ERROR", "MOVE.SPECIES_NO_GENUS",
                                     f"Change={ACT}, can not find parent genus",
                                     paste0("parent '", destParentTaxon.get("name"), "' is a ",
                                            destParentTaxon.get("rank")))
                    continue
                if not binomial_ok(destTaxonName, genus):
                    log_change_error(cur, "ERROR", "MOVE.SPECIES_BINOMIAL_MISMATCH",
                                     f"Change={ACT}, but proposed species names does not start with 'genus[space]' per binomial naming convention",
                                     paste0("parent genus name =", genus.get("name")))
                    continue

            # MOVE the taxon. Its own lineage text is not updated (R open defect);
            # see PYTHON_ONLY["move_lineage_text"] for the export.
            newMSL.set_name(target, destTaxonName)
            target["level_id"] = rank_id(destTaxonRank)
            target["rank"] = destTaxonRank
            newMSL.set_parent(target, destParentTaxon.get("taxnode_id"))
            target[".actionOrder"] = order
            if not is_na(cur.get("molecule")):
                target["molecule_id"] = molecule_id(cur["molecule"])
            if not is_na(cur.get("hostSource")):
                target["host_source"] = cur["hostSource"]
            if not is_na(cur.get("comments")):
                target["notes"] = cur["comments"]
            if destTaxonRank == "species":
                for column in ("exemplarAccession", "exemplarName", "Abbrev", "exemplarIsolate",
                               "genomeCoverage"):
                    if not is_na(cur.get(column)):
                        target[XLSX2DB_MAP[column]] = cur[column]
            # set the lineage of the moved taxon's kids
            moveDiff = diff_lineages(cur.get(".srcLineage"), target.get("lineage"))
            update_lineage(target.get("taxnode_id"), target.get("lineage"),
                           paste0(destParentTaxon.get(".otherLineage"), ";", target.get("name")),
                           paste0(proposalZip, ":", linenum),
                           paste0(act, " ", cur.get(".srcRank"), " ", moveDiff))

            created = is_na(srcPrevTaxnodeID)
            if created:
                # the taxon was CREATED in newMSL: just extend its in_*
                target["in_filename"] = paste0(target.get("in_filename"), ";", proposalZip)
                target["in_notes"] = paste0(target.get("in_notes"), ";linenum=", linenum)
                target["in_target"] = destLineage
            elif act == "split":
                # SPLIT= : a split line kept the original name, so the original is not abolished
                for row in srcPrevTargets:
                    row[".split_kept"] = True
                target[".split_kept"] = True
                target[".split_code"] = code
                target[".split_actionOrder"] = order
                target[".split_linenum"] = linenum
                target["in_change"] = act
                target["in_filename"] = _append_filename(target.get("in_filename"), proposalZip)
                target["in_target"] = srcLineage
                target["in_notes"] = paste0("linenum=", linenum)
                if act not in inChangeLevels:
                    inChangeLevels.append(act)
            elif act == "merge":
                # MERGE= : another taxon is merged into this one, which keeps its name
                target[".merge_kept"] = True
                target[".merge_code"] = code
                target[".merge_actionOrder"] = order
                target[".merge_linenum"] = linenum
                for row in srcPrevTargets:
                    row[".merge_kept"] = True
                    row["out_change"] = act
                    row["out_filename"] = _append_filename(row.get("out_filename"), proposalZip)
                    row["out_target"] = destLineage
                    row["out_notes"] = paste0("linenum=", linenum)
                    row[".actionOrder"] = order
            else:
                # MOVE, PROMOTE, DEMOTE
                for row in srcPrevTargets:
                    row["out_updated"] = True
                    row["out_change"] = act
                    row["out_filename"] = _append_filename(row.get("out_filename"), proposalZip)
                    row["out_target"] = destLineage
                    row["out_notes"] = paste0("linenum=", linenum)
                    row[".out_taxnode_id"] = target.get("taxnode_id")
                    row[".actionOrder"] = order
            success_code = {"move": "MOVE.OK", "promote": "PROMOTE.OK", "demote": "DEMOTE.OK",
                            "split": "SPLIT=.OK", "merge": "MERGE=.OK"}[act]
            if act == "move":
                appliedRenameMoves[row_key] = {"taxnode_id": target.get("taxnode_id"),
                                               "name": rstr(target.get("name")),
                                               "rank": rstr(target.get("rank"))}
            if not created:
                previous: Mapping[str, Any] = srcPrevTargets[0] if srcPrevTargets else {}
                log_change_error(cur, "SUCCESS", success_code, f"Change={ACT}, applied successfully",
                                 paste0(ACT, " ", previous.get("rank"), " named '", previous.get("name"),
                                        "' to ", target.get("rank"), " named '", target.get("name"),
                                        "' PROPOSED//CUR=",
                                        diff_lineages(target.get("lineage"), previous.get("lineage"))))
            else:
                log_change_error(cur, "SUCCESS", success_code, f"Change={ACT}, applied successfully",
                                 paste0(ACT, " ", cur.get(".srcRank"), ":'", cur.get(".srcTaxon"),
                                        "' to ", target.get("rank"), ":'", target.get("name"),
                                        "' PROPOSED//CUR=",
                                        diff_lineages(cur.get(".destLineage"), cur.get(".srcLineage"))))

        elif act == "merge":
            ###### MERGE != ######
            # curMSL    newMSL  case
            #   A   ->   B      (B exists; just set out_*)
            #   A   ->   C      (C does not exist; like rename: set out_*)
            if is_na(srcTaxonName):
                log_change_error(cur, "ERROR", "MERGE.WO_SRC",
                                 "Change=MERGE, but 'current taxonomy' columns are empty",
                                 paste0("currentTaxonomy=", srcLineage, ", destTaxonomy=", destLineage))
                continue

            ##### find SRC #####
            # first look in curMSL
            srcCurTaxa = cur_by_name.get(srcTaxonName, [])
            if not srcCurTaxa:
                # not found in curMSL - see if it's been renamed, look in newMSL
                srcNewTaxa = newMSL.named(srcTaxonName)
                if not srcNewTaxa:
                    # (R passes these two texts in the other order)
                    log_change_error(cur, "ERROR", "MERGE.SRC_NOT_FOUND",
                                     paste0("Looked for=", srcTaxonName, "; merge(CURRENT//PROPOSED)=",
                                            diff_lineages(srcLineage, destLineage)),
                                     "Change=MERGE, but taxon not found")
                    continue
                # found in newMSL, but not in curMSL; go find original curMSL name
                srcCurTaxa = _cur_rows(srcNewTaxa[0].get(".prev_taxnode_id"))
                if not srcCurTaxa:
                    # created in this MSL, so there is no current taxon to merge
                    log_change_error(cur, "ERROR", "MERGE.SRC_NOT_IN_CUR",
                                     "Change=MERGE, but taxon does not exist in the current MSL",
                                     paste0("Looked for=", srcTaxonName,
                                            "; found only in new MSL, from proposal(s)=",
                                            srcNewTaxa[0].get("in_filename")))
                    continue
                log_change_error(cur, "WARNING", "MERGE.SRC_ALREADY_CHANGED",
                                 "Change=MERGE, but current taxon already changed",
                                 paste0("Looked for=", srcTaxonName, "; ORIGNAL//FOUND=",
                                        diff_lineages(srcCurTaxa[0].get("lineage"), srcNewTaxa[0].get("lineage")),
                                        "; proposal(s)=", srcCurTaxa[0].get("out_filename"),
                                        " already did a '", srcCurTaxa[0].get("out_change"),
                                        "' to '", srcCurTaxa[0].get("out_target"), "'"))
            else:
                # found in curMSL: use its copy in newMSL, by taxnode_id
                srcNewTaxa = newMSL.from_prev(srcCurTaxa[0].get("taxnode_id"))
                if len(srcNewTaxa) != 1:
                    log_change_error(cur, "ERROR", "MERGE.SRC_NOT_IN_NEW",
                                     "Change=MERGE, but taxon is no longer in the new MSL",
                                     paste0("Looked for=", srcTaxonName, "; matches in new MSL=", len(srcNewTaxa),
                                            "; proposal(s)=", srcCurTaxa[0].get("out_filename"),
                                            " already did a '", srcCurTaxa[0].get("out_change"),
                                            "' to '", srcCurTaxa[0].get("out_target"), "'"))
                    continue
            srcCurTaxon = srcCurTaxa[0]
            srcNewTaxon = srcNewTaxa[0]

            ##### tag SRC in newMSL #####
            srcCurTaxon[".merge"] = True
            srcCurTaxon[".merge_code"] = code
            srcCurTaxon[".merge_linenum"] = linenum
            srcCurTaxon[".merge_actionOrder"] = order
            srcCurTaxon[".actionOrder"] = order
            srcNewTaxon[".merge"] = True
            srcNewTaxon[".merge_code"] = code
            srcNewTaxon[".merge_linenum"] = linenum
            srcNewTaxon[".merge_actionOrder"] = order

            ##### look for DEST #####
            foundMergeTarget = False
            pendingDestChange = None  # set if a later rename/move in this proposal gives the target its proposed name
            destNewTaxa = newMSL.named(destTaxonName)
            if len(destNewTaxa) == 1:
                foundMergeTarget = True
            elif not destNewTaxa:
                # check if someone else already modified it
                destCurTaxa = cur_by_name.get(destTaxonName, []) if not is_na(destTaxonName) else []
                if not destCurTaxa:
                    # Not in the current or new MSL. Merges run before renames/moves, so
                    # the proposed name may come from a rename/move later in this proposal
                    # (eg 2026.009F). If so, merge into the taxon that change renames/moves:
                    # it keeps its taxnode_id, and out_target is set to its final name
                    # after all changes are applied (R 971f6dc).
                    pending = [row for row in allChangeDf
                               if row.get(".code") == code and row.get(".action") in ("rename", "move")
                               and not is_na(destTaxonRank) and row.get(".destRank") == destTaxonRank
                               and not is_na(destTaxonName) and row.get(".destTaxon") == destTaxonName
                               and not is_na(row.get(".srcTaxon"))
                               and row.get(".codeLine") != cur.get(".codeLine")]
                    if len(pending) == 1:
                        pendingDestChange = pending[0]
                        destNewTaxa = newMSL.named(pendingDestChange.get(".srcTaxon"))
                    if len(pending) == 1 and len(destNewTaxa) == 1:
                        # found it, under its current name
                        pendingDestNotes = paste0(
                            "line ", pendingDestChange.get(".linenum"), " does a '",
                            pendingDestChange.get("change"), "' of ", pendingDestChange.get(".srcRank"),
                            " '", pendingDestChange.get(".srcTaxon"), "' to '",
                            pendingDestChange.get(".destTaxon"), "'")
                        if destNewTaxa[0].get("taxnode_id") == srcNewTaxon.get("taxnode_id"):
                            # MERGE A->C and RENAME A->C: the merge source itself is renamed
                            log_change_error(cur, "ERROR", "MERGE.DEST_IS_SRC",
                                             "Change=MERGE, but the proposed taxon is the current taxon, renamed/moved later in this proposal; a taxon can not be merged into itself",
                                             pendingDestNotes)
                            continue
                        foundMergeTarget = True
                        log_change_error(cur, "INFO", "MERGE.DEST_PENDING_RENAME",
                                         "Change=MERGE, proposed taxon is renamed/moved to that name later in this proposal; merged into that taxon",
                                         paste0(pendingDestNotes, "; merge into ", destNewTaxa[0].get("lineage")))
                    elif pending:
                        # A change in this proposal gives a taxon the proposed name, but that
                        # taxon can't be identified, so that change will fail too.
                        lines = "; ".join(
                            paste0("line ", row.get(".linenum"), " does a '", row.get("change"),
                                   "' of ", row.get(".srcRank"), " '", row.get(".srcTaxon"), "'")
                            for row in pending)
                        log_change_error(cur, "ERROR", "MERGE.DEST_PENDING_RENAME_FAILED",
                                         "Change=MERGE, proposed taxon is renamed/moved to that name later in this proposal, but the taxon that change renames/moves was not found",
                                         paste0(lines, "; matches in new MSL for the renamed/moved taxon: ",
                                                len(destNewTaxa) if len(pending) == 1
                                                else "not checked, several changes"))
                        continue
                    else:
                        # just not there - we'll be creating it!
                        pendingDestChange = None
                        destNewTaxa = []
                elif len(destCurTaxa) == 1:
                    # previous record - must have been modified already this MSL
                    destCurTaxon = destCurTaxa[0]
                    destNewTaxa = newMSL.from_prev(destCurTaxon.get("taxnode_id"))
                    if len(destNewTaxa) == 1:
                        foundMergeTarget = True
                        # (R passes these two texts in the other order)
                        log_change_error(cur, "WARNING", "MERGE.DEST_ALREADY_CHANGED",
                                         paste0("proposal(s)=", destCurTaxon.get("out_filename"),
                                                " also did a '", destCurTaxon.get("out_change"),
                                                "' to '", destCurTaxon.get("out_target"), "'",
                                                "; PROPOSED//OBSERVED=",
                                                diff_lineages(destLineage, destCurTaxon.get("lineage"))),
                                         "Change=MERGE, but proposed taxon already modified")
                    elif not destNewTaxa and not is_na(destCurTaxon.get(".merge_src_taxnode_id")):
                        # the proposed taxon was merged into another taxon already (eg MERGE
                        # B->C ran before MERGE A->B): merge into the taxon it was merged into.
                        mergedDestId = newmsl_resolve_merged_id(destCurTaxon.get(".out_taxnode_id"))
                        mergedDestNotes = paste0(
                            "proposal(s)=", destCurTaxon.get("out_filename"), " line ",
                            destCurTaxon.get(".merge_linenum"), " merged '", destCurTaxon.get("name"),
                            "' into '", destCurTaxon.get("out_target"), "'")
                        if is_na(mergedDestId):
                            log_change_error(cur, "ERROR", "MERGE.DEST_NOT_IN_NEW",
                                             "Change=MERGE, but the proposed taxon was merged into a taxon that is no longer in the new MSL",
                                             mergedDestNotes)
                            continue
                        if mergedDestId == srcNewTaxon.get("taxnode_id"):
                            # MERGE A->B, then MERGE B->A: A is already part of B
                            log_change_error(cur, "ERROR", "MERGE.DEST_IS_SRC",
                                             "Change=MERGE, but the proposed taxon was already merged into the current taxon; a taxon can not be merged into itself",
                                             mergedDestNotes)
                            continue
                        destNewTaxa = newMSL.with_id(mergedDestId)
                        foundMergeTarget = True
                        log_change_error(cur, "WARNING", "MERGE.DEST_MERGED",
                                         "Change=MERGE, but the proposed taxon is itself merged into another taxon; merged into that taxon",
                                         paste0(mergedDestNotes, "; merge into ", destNewTaxa[0].get("lineage")))
                    elif len(destNewTaxa) > 1:
                        log_change_error(cur, "ERROR", "MERGE.DEST_MANY",
                                         paste0("proposal(s)=", destNewTaxa[0].get("in_filename"),
                                                " already did a SPLIT? "),
                                         "Change=MERGE, but proposed taxon already split")
                        continue
                else:
                    log_change_error(cur, "ERROR", "MERGE.DEST_MANY_2",
                                     paste0("proposal(s)=", destCurTaxa[0].get("out_filename"),
                                            " already did a MERGED? "),
                                     "Change=MERGE, but multiple taxa")
                    continue

            # a taxon can not be merged into itself: merging would remove it from the new MSL
            if foundMergeTarget and any(row.get("taxnode_id") == srcNewTaxon.get("taxnode_id")
                                        for row in destNewTaxa):
                log_change_error(cur, "ERROR", "MERGE.DEST_IS_SRC",
                                 "Change=MERGE, but the proposed taxon is the current taxon; a taxon can not be merged into itself",
                                 paste0("current taxon '", srcNewTaxon.get("lineage"),
                                        "', proposed taxon found as '", destNewTaxa[0].get("lineage"), "'"))
                continue

            ##### rename/move or merge #####
            if foundMergeTarget:
                ###### . merge into existing ######
                destNewTaxon = destNewTaxa[0]
                srcCurTaxon["out_updated"] = True
                srcCurTaxon["out_change"] = "merge"
                srcCurTaxon["out_filename"] = proposalZip
                srcCurTaxon["out_target"] = destNewTaxon.get("lineage")
                srcCurTaxon["out_notes"] = paste0("linenum=", linenum)
                srcCurTaxon[".out_taxnode_id"] = destNewTaxon.get("taxnode_id")
                srcCurTaxon[".actionOrder"] = order
                # this change, to report on it after all changes
                srcCurTaxon[".merge_row"] = row_key
                srcCurTaxon[".merge_order"] = order
                destNewTaxon[".actionOrder"] = order
                # remove the now abandoned record for src in newMSL: check for kids
                kids = newMSL.kids(srcNewTaxon.get("taxnode_id"))
                if kids:
                    log_change_error(cur, "ERROR", "MERGE_INTO.KIDS.UNIMP",
                                     "Change=MERGE, merging child taxa not implemented (yet) ",
                                     [paste0("taxon='", srcTaxonName, "', lineage='", srcLineage,
                                             "', kids: N=", len(kids), ", NAMES=[", name, "]")
                                      for name in kid_names(kids)])
                    continue
                newMSL.remove(srcNewTaxon)
                # record which newMSL taxon was removed, so a merge into it can be followed
                srcCurTaxon[".merge_src_taxnode_id"] = srcNewTaxon.get("taxnode_id")
                # Remember the exact dependent rename/move and intended identity, so its
                # successful outcome can be checked even after further renames or merges.
                if pendingDestChange is not None:
                    srcCurTaxon[".merge_pending_name"] = rstr(destTaxonName)
                    srcCurTaxon[".merge_pending_rank"] = rstr(destTaxonRank)
                    srcCurTaxon[".merge_pending_linenum"] = rstr(pendingDestChange.get(".linenum"))
                    srcCurTaxon[".merge_pending_row"] = pendingDestChange.get(".codeLine")
                    srcCurTaxon[".merge_pending_taxnode_id"] = destNewTaxon.get("taxnode_id")
                logNotes = paste0(ACT, " ", srcCurTaxon.get("rank"), " named '", srcCurTaxon.get("name"),
                                  "' into ", destNewTaxon.get("rank"), " named '", destNewTaxon.get("name"),
                                  "' CUR//PROPOSED=",
                                  diff_lineages(srcCurTaxon.get("lineage"), destNewTaxon.get("lineage")))
                # merging child nodes is not implemented, so only SPECIES
                if srcCurTaxon.get("rank") != "species":
                    log_change_error(cur, "ERROR", "MERGE_INTO.SRC_NOT_SPECIES.UNIMP",
                                     "Change=MERGE only species merging has been implemented", logNotes)
                log_change_error(cur, "SUCCESS", "MERGE_INTO.OK",
                                 f"Change={ACT}, applied successfully", logNotes)
            else:
                ###### . move/rename taxon ######
                parents = root_rows() if is_na(destParentName) else newMSL.named(destParentName)
                if not parents:
                    # look for the name in curMSL, and assume it got renamed already in newMSL
                    parentDestCurTaxa = cur_by_name.get(destParentName, []) if not is_na(destParentName) else []
                    if len(parentDestCurTaxa) == 1:
                        parents = newMSL.from_prev(parentDestCurTaxa[0].get("taxnode_id"))
                        log_change_error(cur, "INFO", "MERGE_MOVE_RENAME.PARENT_ALREADY_CHANGED",
                                         "PROPOSED parent taxon already modified", other_lineage_notes(parents))
                    elif not parentDestCurTaxa:
                        changedParent = [row for row in allChangeDf if row.get(".srcTaxon") == destParentName]
                        if len(changedParent) == 1:
                            parents = newMSL.named(changedParent[0].get(".destTaxon"))
                            log_change_error(cur, "INFO", "MERGE_MOVE_RENAME.PARENT_ALREADY_CHANGED_2",
                                             "PROPOSED parent taxon already modified",
                                             other_lineage_notes(parents))
                        else:
                            log_change_error(cur, "ERROR", "MERGE_MOVE_RENAME.PARENT_NOT_FOUND_2",
                                             "PROPOSED parent taxon not found",
                                             paste0("'", destParentName, "' : ", destParentLineage))
                            continue
                    else:
                        log_change_error(cur, "ERROR", "MERGE_MOVE_RENAME.PARENT_MULTIPLE_MATCH",
                                         "Probably parent was split by another proposal",
                                         paste0("parent '", destParentName, "' matches ", len(parentDestCurTaxa),
                                                " taxa in the current MSL: ",
                                                ", ".join(rstr(row.get("lineage")) for row in parentDestCurTaxa)))
                        continue
                if not parents:
                    # (R stops here: it reads the rank of a parent it did not find)
                    log_change_error(cur, "ERROR", "MERGE_MOVE_RENAME.PARENT_NOT_FOUND_2",
                                     "PROPOSED parent taxon not found",
                                     paste0("'", destParentName, "' : ", destParentLineage))
                    continue
                parent = parents[0]
                ###### . create binomial check ######
                if destTaxonRank == "species":
                    status, genus = _find_genus(parent)
                    if status == "subgenus":
                        log_change_error(cur, "ERROR", "MERGE_CREATE.SPECIES_SUBGENUS_NO_GENUS",
                                         "Change=MERGE_CREATE, can not find parent genus",
                                         paste0("subgenus '", parent.get("name"), "' is not in a genus;",
                                                " it's parent '", genus.get("name"), "' is a ", genus.get("rank")))
                        continue
                    if status == "none":
                        log_change_error(cur, "ERROR", "MERGE_CREATE.SPECIES_NO_GENUS",
                                         "Change=MERGE_CREATE, can not find parent genus",
                                         paste0("parent '", parent.get("name"), "' is a ", parent.get("rank")))
                        continue
                    if not binomial_ok(destTaxonName, genus):
                        log_change_error(cur, "ERROR", "CREATE_RENAME.SPECIES_BINOMIAL_MISMATCH",
                                         "Change=RENAME, but proposed species names does not start with 'genus[space]' per binomial naming convention",
                                         paste0("GENUS//SPECIES =", diff_strings(genus.get("name"), destTaxonName)))
                        continue
                # add to list of genera to check for binomial problems when all changes are applied
                if destTaxonRank == "genus":
                    renamedGenera.append({"name": destTaxonName, "code": code, "linenum": linenum,
                                          "actionOrder": order})
                srcNewTaxon["prev_proposals"] = paste0(
                    "" if is_na(srcNewTaxon.get("prev_proposals")) else paste0(srcNewTaxon["prev_proposals"], ","),
                    proposalZip, ":", linenum)
                # update (rename/move) dest in newMSL
                newMSL.set_name(srcNewTaxon, destTaxonName)
                srcNewTaxon["lineage"] = paste0(parent.get("lineage"), ";", destTaxonName)
                srcNewTaxon[".actionOrder"] = order
                # put out_* changes on curMSL
                srcCurTaxon["out_updated"] = True
                srcCurTaxon["out_change"] = "merge"
                srcCurTaxon["out_filename"] = proposalZip
                srcCurTaxon["out_target"] = srcNewTaxon["lineage"]
                srcCurTaxon["out_notes"] = paste0("linenum=", linenum)
                srcCurTaxon[".out_taxnode_id"] = srcNewTaxon.get("taxnode_id")
                srcCurTaxon[".actionOrder"] = order
                srcCurTaxon[".merge_row"] = row_key
                srcCurTaxon[".merge_order"] = order
                log_change_error(cur, "SUCCESS", "MERGE_MOVE_RENAME.OK", "Change=MERGE, applied successfully",
                                 paste0("MERGE ", srcTaxonRank, " from '", srcTaxonName, "' to '",
                                        srcNewTaxon.get("name"), "' in '", parent.get("lineage"),
                                        "', CURRENT//PROPOSED=",
                                        diff_lineages(srcCurTaxon.get("lineage"), srcNewTaxon["lineage"])))
        else:
            ##### UNKNOWN ACTION not implemented. #####
            log_change_error(cur, "ERROR", "CHANGE.UNIMP",
                             f"Change={CHANGE} is NOT (yet) implemented",
                             paste0("lineageChange=", diff_lineages(srcLineage, destLineage)))
            continue

    #### POST-QC ####
    actionOrder += 1

    ##### abolish taxa - deferred / split taxa - delete #####
    # ABOLISH of a taxon that still had kids was deferred, as the kids may have
    # been moved out by a later change. If a taxon is split, and no split keeps
    # its original name, the original name is removed from the new MSL. Either
    # may be the last kid of the other, so remove both together, lowest rank
    # first, so kids are removed before their parents are checked.
    #
    # A split source that another change in the same proposal renamed, moved,
    # promoted or demoted still exists under its new name/rank/parent, so that
    # change keeps it (R 971f6dc).
    for row in newMSL.rows:
        if row.get(".split") and not row.get(".split_kept"):
            previous_rows = _cur_rows(row.get(".prev_taxnode_id"))
            if len(previous_rows) != 1:
                continue
            splitZip = paste0(_proposal_value(row.get(".split_code"), "basename"), ".zip")
            out_files = ([] if is_na(previous_rows[0].get("out_filename"))
                         else str(previous_rows[0]["out_filename"]).split(";"))
            if (previous_rows[0].get("out_change") in ("rename", "move", "promote", "demote")
                    and splitZip in out_files):
                row[".split_kept"] = True
    splitDelete = [row for row in newMSL.rows if row.get(".split") and not row.get(".split_kept")]
    if splitDelete:
        actionOrder += 1
        log_error([row.get(".split_code") for row in splitDelete],
                  [row.get(".split_linenum") for row in splitDelete],
                  "split_abolish", [rstr(row.get("rank")) for row in splitDelete],
                  [row.get("name") for row in splitDelete], "INFO", "SPLIT.IMPLICIT_ABOLISH",
                  "Change=SPLIT, but no split line kept original name, so remove original name from new MSL",
                  [paste0("removed ", row.get("rank"), " ", row.get("name"), " from ", row.get("lineage"))
                   for row in splitDelete], order)
    # check for kids - by taxnode_id for stability across row deletions
    split_ids = {id(row) for row in splitDelete}
    removals = [row for row in newMSL.rows if row.get(".abolish") or id(row) in split_ids]
    removals.sort(key=lambda row: -(row.get("level_id") or 0))
    for row in removals:
        kids = newMSL.kids(row.get("taxnode_id"))
        if row.get(".abolish") is True:
            # deferred ABOLISH: restore the original ABOLISH change
            cur = dict(changes_by_row[row[".abolish_row"]])
            cur[".order"] = row.get(".abolish_actionOrder")
            if kids:
                src = [cur.get(name) for name in xlsx_change_src_colnames]
                log_change_error(cur, "ERROR", "ABOLISH.KIDS",
                                 "Change=ABOLISH, taxon still has un-abolished children",
                                 [paste0("taxon=", taxon_get_name(src), ", lineage=", taxon_get_lineage(src),
                                         ", kids: N=", len(kids), ", NAMES=[", name, "]")
                                  for name in kid_names(kids)])
            else:
                abolish_taxon(cur, row)
        else:
            # SPLIT implicit abolish
            actionOrder += 1
            if kids:
                log_error(row.get(".split_code"), row.get(".split_linenum"), "split_abolish",
                          rstr(row.get("rank")), row.get("name"), "ERROR",
                          "SPLIT.IMPLICIT_ABOLISH_WITH_KIDS",
                          "Change=ABOLISH, taxon still has un-abolished/moved children",
                          [paste0("taxon=", row.get("name"), ", lineage=", row.get("lineage"),
                                  ", kids: N=", len(kids), ", NAMES=[", name, "]")
                           for name in kid_names(kids)], order)
            else:
                newMSL.remove(row)

    # remaining split taxa that kept their name
    splitKeep = [row for row in newMSL.rows if row.get(".split") and row.get(".split_kept")]
    if splitKeep and params.verbose:
        log_error([row.get(".split_code") for row in splitKeep],
                  [row.get(".split_linenum") for row in splitKeep], "split_keep",
                  [row.get("rank") for row in splitKeep], [row.get("name") for row in splitKeep],
                  "INFO", "SPLIT.KEEP", "Change=SPLIT, but one split directive kept the original name",
                  [paste0("split and keep ", row.get("rank"), " '", row.get("name"), "' from ",
                          row.get("lineage")) for row in splitKeep], order)

    ##### empty taxa #####
    # did proposal create empty (non-species) taxa? (global scan)
    # It is reported under the last change's proposal, as in R; the notes say
    # which proposal emptied the taxon.
    emptyTaxa = [row for row in newMSL.rows
                 if not newMSL.kids(row.get("taxnode_id")) and row.get("level_id") != 600
                 and is_na(row.get(".emptyReported"))]
    for emptyTaxon in emptyTaxa:
        errLineNum, errChange = "", ""
        errText = "check if lower rank taxa were moved elsewhere?"
        if not is_na(emptyTaxon.get(".split_linenum")):
            errLineNum = emptyTaxon[".split_linenum"]
            errChange = "Split"
            errText = paste0(_proposal_value(emptyTaxon.get(".split_code"), "basename"), ":",
                             emptyTaxon[".split_linenum"], " did a split")
        elif not is_na(emptyTaxon.get("prev_proposals")):
            errChange = "previous"
            errText = emptyTaxon["prev_proposals"]
        elif not is_na(emptyTaxon.get(".otherLineageProposal")):
            errChange = "previous"
            errText = paste0(emptyTaxon[".otherLineageProposal"], " ", emptyTaxon.get(".otherLineageAction"))
        log_error(code, errLineNum, errChange, emptyTaxon.get("rank"), emptyTaxon.get("name"),
                  "ERROR", "PROPOSAL.EMPTY_TAXA", "Proposal created empty (non-species) taxa",
                  paste0("the ", emptyTaxon.get("rank"), " '", emptyTaxon.get("name"),
                         "' is empty - it does not contain any lower rank taxons; ", errText), order)
        # mark empty taxa so we don't re-report them
        for row in newMSL.named(emptyTaxon.get("name")):
            row[".emptyReported"] = _proposal_value(code, "xlsx")

    ##### renamed genera - binomial #####
    for entry in renamedGenera:
        actionOrder += 1
        genera = newMSL.named(entry["name"])
        if not genera:
            continue
        genus_ids = [row.get("taxnode_id") for row in genera]
        subgenera = [kid for genus_id in genus_ids for kid in newMSL.kids(genus_id)
                     if kid.get("level_id") == 550]
        parent_ids = set(genus_ids) | {row.get("taxnode_id") for row in subgenera}
        species = [row for row in newMSL.rows
                   if row.get("parent_id") in parent_ids and row.get("level_id") == 600]
        prefix = rstr(genera[0].get("name")) + " "
        violations = sorted((rstr(row.get("name")) for row in species
                             if not rstr(row.get("name")).startswith(prefix)), key=r_collate_key)
        if violations:
            log_error(entry["code"], entry["linenum"], "rename_genus", "species", violations, "ERROR",
                      "RENAME_GENUS.SPECIES_BINOMIAL_MISMATCH",
                      "Change=RENAME_GENUS, but species name does not start with 'genus[space]' per binomial naming convention",
                      paste0("new genus name: ", genera[0].get("name")), entry["actionOrder"])

    ##### merged taxa - final merge target #####
    # A merge records its target's lineage in out_target when it is applied.
    # Merges run before renames/moves, so a later change (eg renaming the target
    # or its genus) can leave out_target naming a taxon that is not in the new
    # MSL; the database load (rebuild_delta_nodes, which matches out_target to
    # the new MSL's lineage or name) then can not find the merge target.
    # If the target's lineage changed after the merge, record its final name.
    for merged in curMSL.rows:
        if merged.get("out_change") != "merge" or is_na(merged.get(".out_taxnode_id")):
            continue
        # the merge's change line, to report on (None if not recorded)
        mergeChange = None
        if not is_na(merged.get(".merge_row")):
            mergeChange = dict(changes_by_row[merged[".merge_row"]])
            mergeChange[".order"] = merged.get(".merge_order")
        mergeNotes = paste0("merged ", merged.get("rank"), " '", merged.get("name"),
                            "' into '", merged.get("out_target"), "'")
        # follow merges of the target
        recordedTargetId = merged[".out_taxnode_id"]
        mergeTargetId = newmsl_resolve_merged_id(recordedTargetId)
        if is_na(mergeTargetId):
            # target no longer in newMSL: leave out_target as is
            if mergeChange is not None:
                log_change_error(mergeChange, "ERROR", "MERGE.DEST_NOT_IN_NEW",
                                 "Change=MERGE, but the taxon it was merged into is no longer in the new MSL",
                                 mergeNotes)
            continue
        mergeTarget = newMSL.with_id(mergeTargetId)[0]
        if mergeTargetId != recordedTargetId:
            if mergeChange is not None:
                log_change_error(mergeChange, "WARNING", "MERGE.DEST_MERGED",
                                 "Change=MERGE, but the proposed taxon is itself merged into another taxon; merged into that taxon",
                                 paste0(mergeNotes, ", which was then merged into '",
                                        newmsl_lineage_from_parents(mergeTargetId), "'"))
            merged[".out_taxnode_id"] = mergeTargetId
        # check a merge into a taxon renamed/moved later in its proposal: by the
        # recorded outcome of that change, not by the final name, which later
        # changes may replace
        if not is_na(merged.get(".merge_pending_name")):
            pendingName = merged[".merge_pending_name"]
            pendingRank = merged.get(".merge_pending_rank")
            appliedChange = appliedRenameMoves.get(merged.get(".merge_pending_row"))
            pendingSucceeded = (appliedChange is not None
                                and appliedChange["taxnode_id"] == merged.get(".merge_pending_taxnode_id")
                                and appliedChange["name"] == pendingName
                                and appliedChange["rank"] == pendingRank)
            if not pendingSucceeded and mergeChange is not None:
                log_change_error(mergeChange, "ERROR", "MERGE.DEST_PENDING_RENAME_FAILED",
                                 "Change=MERGE, merged into a taxon that a later rename/move in this proposal was to give the proposed name, but that change did not",
                                 paste0("line ", merged.get(".merge_pending_linenum"),
                                        " was to rename/move it to ", pendingRank, " '", pendingName, "'; ",
                                        mergeNotes, ", which is now ", mergeTarget.get("rank"), " '",
                                        newmsl_lineage_from_parents(mergeTargetId), "'"))
        if rstr(merged.get("out_target")) != newmsl_lineage_from_parents(mergeTarget.get("taxnode_id")):
            merged["out_target"] = rstr(mergeTarget.get("name"))


def report_applied(changes: list[dict[str, Any]]) -> None:
    """The per-proposal and total "# APPLIED" lines R prints after applying."""
    def level_summary(rows: Iterable[Mapping[str, Any]]) -> str:
        counts: dict[str, int] = {}
        for row in rows:
            if not is_na(row.get("level")):
                counts[row["level"]] = counts.get(row["level"], 0) + 1
        return " ".join(f"{level}:{counts[level]}" for level in allErrorDf.levels if counts.get(level))

    applied = [row for row in allErrorDf.rows
               if not is_na(_as_numeric(row.get("row"))) and _as_numeric(row.get("row")) > 0]
    codes = sorted({row[".code"] for row in changes}, key=r_collate_key)
    for code in codes:
        count = sum(1 for row in changes if row[".code"] == code)
        summary = level_summary(row for row in applied if row.get("code") == code)
        print(f"# APPLIED: {count:3d} changes from {rstr(_proposal_value(code, 'basename'))}  with {summary} ")
    print(f"# TOTAL  : {len(changes):3d} changes  from {len(codes)} proposals  "
          f"with {level_summary(allErrorDf.rows)} ")


# ---------------------------------------------------------------------------
# Output generation
# ---------------------------------------------------------------------------
#
# The files are written byte for byte as R writes them, including the spaces
# cat() puts between its arguments.

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


def _short_proposal_filename(value: Any) -> Any:
    if is_na(value):
        return value
    return re.sub(r"^([0-9]+\.[0-9A-Z]+)\..*", r"\1...", rstr(value), count=1, flags=re.S)


def _export_lineage(row: Mapping[str, Any]) -> Any:
    """The lineage text to export for a newMSL row."""
    if not PYTHON_ONLY["move_lineage_text"]:
        return row.get("lineage")
    if is_na(row.get("level_id")) or row["level_id"] <= 100:
        return row.get("lineage")
    return newmsl_lineage_from_parents(row.get("taxnode_id"))


def export_msl_tsv() -> None:
    filename = os.path.join(params.out_dir, "msl.tsv")
    print(f"Writing  {filename} ")
    with open(filename, "w", encoding="utf-8", newline="") as out:
        out.write("\t".join(TSV_COLUMNS) + " \n")
        for row in newMSL.rows:
            values = []
            for column in TSV_COLUMNS:
                value = row.get(column)
                if column in ("in_filename", "out_filename"):
                    value = _short_proposal_filename(value)
                elif column == "lineage":
                    value = _export_lineage(row)
                elif column == "in_change" and not is_na(value) and not PYTHON_ONLY["in_change_text"]:
                    # R writes the factor level number, not the word
                    value = inChangeLevels.index(value) + 1 if value in inChangeLevels else value
                # map NA's to NULL's
                values.append("NULL" if is_na(value) else rstr(value))
            out.write("\t".join(values) + " \n")
    print(f"WROTE   {filename} ({len(newMSL.rows)} rows)")


def export_docx_summary() -> None:
    columns = ["subcommittee", "code", "docx", "xlsx", "title",
               "authorsEmails", "correspondingAuthor", "abstract"]
    filename = os.path.join(params.out_dir, params.proposals_meta)
    print(f"Writing  {filename} ")
    with open(filename, "w", encoding="utf-8", newline="") as out:
        out.write("\t".join(columns) + " \n")
        for proposal in proposalsDf.values():
            # a column no DOCX filled is ""; a column this DOCX did not fill is NA
            out.write("\t".join("" if column not in proposalsDfColumns else rstr(proposal.get(column))
                                for column in columns) + " \n")
    print(f"WROTE   {filename} ({len(proposalsDf)})")


def sql_value_nullable(value: Any) -> str:
    if is_na(value):
        return "NULL"
    return "'" + rstr(value).replace("'", "''") + "'"


def export_sql() -> None:
    """Write the MariaDB load and post-load maintenance script."""
    sql_msl_num = newMSL.rows[0].get("msl_release_num")
    msl = rstr(sql_msl_num)
    previous_msl = rstr(sql_msl_num - 1)
    sql_tree_id = sorted({rstr(row.get("tree_id")) for row in newMSL.rows
                          if not is_na(row.get("tree_id"))}, reverse=True)[0]
    filename = os.path.join(params.out_dir, params.sql_load_filename)
    print(f"Writing  {filename} ")
    out: list[str] = []
    out.append("-- begin transaction\n")
    out.append("-- rollback transaction\n")
    # MariaDB/MySQL treat backslash as an escape character by default, so a value
    # ending with "\" can accidentally escape the closing quote and break parsing.
    # Force standard single-quote escaping (''), and make backslashes literal.
    out.append("SET SESSION sql_mode = CONCAT(@@SESSION.sql_mode, ',NO_BACKSLASH_ESCAPES');\n\n")

    ##### SQL to add new MSL to `taxonomy_toc` #####
    out.append("insert into `taxonomy_toc` (`tree_id`,`msl_release_num`,`comments`)  values ( "
               + f"{sql_tree_id},{msl},NULL" + " );\n")

    ##### SQL to insert newMSL into `taxonomy_node` #####
    batch_size = int(params.sql_insert_batch_size)
    ordered = sorted(newMSL.rows, key=lambda row: (is_na(row.get("level_id")), row.get("level_id") or 0))
    for row_count, row in enumerate(ordered, 1):
        # insert several rows per batch insert statement: fewer trigger calls
        if row_count % batch_size == 1:
            if row_count > 1:
                out.append(";\n")
            out.append("insert into `taxonomy_node` (`" + "`,`".join(SQL_COLUMNS) + "`)\n values \n")
        else:
            out.append(",")
        values = ",".join(sql_value_nullable(row.get(column)) for column in SQL_COLUMNS)
        out.append("(" + values + ")  -- lineage= " + rstr(_export_lineage(row)) + " \n")
    out.append(";\n")

    # QC queries (the odd spacing is what R's cat()/paste() write)
    out.append("-- QC queries \n")
    out.append("\n".join([
        "",
        "  select level_id,",
        "      count(case when in_change like 'new' then 1 end) as new_ct, ",
        "    \tcount(case when in_change<>'new' then 1 end) as other_ct,",
        "    \tcount(case when in_change is null then 1 end) as null_ct,",
        "    \tcount(*) as total_ct",
        "  from `taxonomy_node`",
        "  where msl_release_num =  " + msl + " ",
        "  group by level_id",
        "  order by level_id;",
        "  ",
    ]))
    out.append("\n")
    out.append("\n".join([
        "-- QC Query",
        "  select ",
        "       level_id,in_change, count(*) as ct",
        "  from `taxonomy_node`",
        "  where msl_release_num =  " + msl + " ",
        "  group by level_id,in_change",
        "  order by level_id,in_change;",
        "  ",
    ]))
    out.append("\n")

    ##### SQL to update prevMSL.out in `taxonomy_node` #####
    print("Writing out_* updates for prevMSL")
    current_columns = ["out_change", "out_target", "out_filename", "out_notes"]
    out.append("-- begin transaction\n")
    out.append("-- rollback transaction\n")
    for row in curMSL.rows:
        if is_na(row.get("out_change")):
            continue
        assignments = ",".join(f"`{column}`={sql_value_nullable(row.get(column))}"
                               for column in current_columns)
        out.append(f"update `taxonomy_node` set {assignments} where `taxnode_id`="
                   + rstr(row.get("taxnode_id")) + " ;\n")

    ranks = [rank.lower() for rank in xlsx_change_ranks]
    out.append(
        "\n-- set for taxonomy_node_compute_indexes\n"
        "SET SESSION max_sp_recursion_depth = 255;\n\n"
        "SET @right_idx := NULL;\n\n"
        "SET @realm_id := NULL, @subrealm_id := NULL, @kingdom_id := NULL, @subkingdom_id := NULL,\n"
        "    @phylum_id := NULL, @subphylum_id := NULL, @class_id := NULL, @subclass_id := NULL,\n"
        "    @order_id := NULL, @suborder_id := NULL, @family_id := NULL, @subfamily_id := NULL,\n"
        "    @genus_id := NULL, @subgenus_id := NULL, @species_id := NULL,\n"
        "    @inher_molecule_id := NULL, @lineage := NULL;\n\n"
        "SET @realm_desc_ct := 0, @subrealm_desc_ct := 0, @kingdom_desc_ct := 0, @subkingdom_desc_ct := 0,\n"
        "    @phylum_desc_ct := 0, @subphylum_desc_ct := 0, @class_desc_ct := 0, @subclass_desc_ct := 0,\n"
        "    @order_desc_ct := 0, @suborder_desc_ct := 0, @family_desc_ct := 0, @subfamily_desc_ct := 0,\n"
        "    @genus_desc_ct := 0, @subgenus_desc_ct := 0, @species_desc_ct := 0;\n\n"
        "-- rebuilt indexes\n"
        "CALL `taxonomy_node_compute_indexes`(\n"
        "    " + sql_tree_id + ",\n"
        "    1,\n"
        "    @right_idx,\n"
        "    1,\n"
        + "".join(f"    @{rank}_id,\n" for rank in ranks)
        + "".join(f"    @{rank}_desc_ct,\n" for rank in ranks)
        + "    @inher_molecule_id,\n"
        "    @lineage\n"
        ");\n"
    )
    out.append("\n".join([
        "",
        "-- build deltas (~7min) ",
        "",
        "  CALL `rebuild_delta_nodes`(NULL);",
        "  ",
        "  -- rebuild merge/split (seconds) ",
        "",
        "  CALL `rebuild_node_merge_split`(); ",
        "  ",
    ]))
    out.append("\n".join([
        "",
        "  -- NOW check if all newMSL have delta in, and prevMSL have delta out",
        "  SELECT 'prevMSL w/o delta to new', COUNT(*) FROM `taxonomy_node`",
        "  WHERE msl_release_num =  " + previous_msl + " ",
        "  AND taxnode_id NOT IN (SELECT prev_taxid FROM taxonomy_node_delta)",
        "  UNION ALL",
        "  SELECT  'newMSL w/o delta to prev',  COUNT(*) FROM `taxonomy_node`",
        "  WHERE msl_release_num =  " + msl + " ",
        "  AND taxnode_id NOT IN (SELECT new_taxid FROM taxonomy_node_delta);",
        "  ",
    ]))
    out.append("\n")
    # R 984fb27: refresh the VMR from the new MSL, then its sort columns, before QC
    out.append("\n".join([
        "-- QC Query",
        "  SELECT ",
        "       level_id,out_change, COUNT(*) as ct",
        "  FROM `taxonomy_node`",
        "  WHERE msl_release_num =  " + previous_msl + " ",
        "  GROUP BY level_id,out_change",
        "  ORDER BY level_id,out_change;",
        "  ",
        "  CALL `VMR_update_from_new_MSL`();",
        "  CALL `species_isolates_update_sorts`();",
        "  ",
        "  -- QC SPs",
        "  CALL `QC_run_modules`(NULL);",
        "  ",
        "  ",
    ]))
    with open(filename, "w", encoding="utf-8", newline="") as stream:
        stream.write("".join(out))
    print(f"WROTE    {filename} ")


def export_change_report(changes: list[dict[str, Any]]) -> None:
    if not params.output_change_report or not changes:
        return
    filename = Path(params.out_dir) / "changes.tsv"
    columns: list[str] = []
    for row in changes:
        for column in row:
            if column not in columns:
                columns.append(column)
    write_tsv(columns, changes, filename)
    if params.verbose:
        print(f"WROTE   {filename} ({len(changes)} rows)")


# ---------------------------------------------------------------------------
# Program entry point and end-to-end lifecycle
# ---------------------------------------------------------------------------


def main(argv: Sequence[str] | None = None) -> int:
    """Run one complete proposal validation and MSL-generation workflow."""
    global params, changeList, allErrorDf, actionOrder, _last_summary_state
    params = parse_args(argv)
    r_compat()
    if params.verbose:
        print(f"PROPOSALS_DIR:  {params.proposals_dir}/")
        print(f"OUT_DIR:        {params.out_dir}/")
    Path(params.out_dir).mkdir(parents=True, exist_ok=True)
    if params.verbose:
        print(f"REF_DIR:        {params.ref_dir}/")
    allErrorDf = ErrorLog()
    actionOrder = 0
    _last_summary_state = None
    if params.use_cache or params.update_cache or params.load_proposal_cache or params.save_proposal_cache:
        print("WARNING: cache flags are accepted for compatibility but caching is not implemented in this first Python version")
    load_version()
    load_reference()
    load_abolished_accessions()

    proposals = scan_for_proposals()
    changeList = load_and_qc_proposals(proposals, {})
    changes = merge_and_order_changes(changeList)
    if changes:
        apply_changes(changes)
        report_applied(changes)
    write_error_summary(final=True)
    for stale in Path(params.out_dir).glob("QC.pretty_summary.?_*.xlsx"):
        stale.unlink(missing_ok=True)
    if params.export_msl:
        export_msl_tsv()
        export_docx_summary()
        export_sql()
    export_change_report(changes)
    print("# COMPLETED.")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        raise SystemExit(130)
    except RStop as exc:
        if str(exc):
            print(f"Error: {exc}", file=sys.stderr)
        raise SystemExit(1)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        if params.get("debug_on_error"):
            raise
        raise SystemExit(1)
