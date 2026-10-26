"""lncRNA identifiers -> HGNC symbols, for the PubMed search.

The network stores lncRNAs as Ensembl transcript ids (ENST...). PubMed abstracts
name lncRNAs by symbol (MAGI2-AS3, PVT1), so each transcript is mapped
    transcript -> gene      via app_v1/lncrna_reference.db  (species, transcript_id, gene)
    gene       -> symbol    via app_v1/lncrna_symbols.tsv   (species, gene, symbol)
The TSV is built by reference_mapping/build_lncrna_symbols.py from the Ensembl
ncRNA dumps' `gene_symbol:` field. An lncRNA with no row in either table keeps
its Ensembl id and is reported as unsearchable on its two arms.
"""
import csv
import os
import sqlite3

from app_v1.lncrna_reference import normalize_lncrna_id

_HERE = os.path.dirname(os.path.abspath(__file__))
_APP = os.path.dirname(_HERE)


def default_symbols_path():
    return os.environ.get("ISOTAR_LNCRNA_SYMBOLS_TSV") or os.path.join(_APP, "lncrna_symbols.tsv")


def default_reference_db():
    env = os.environ.get("ISOTAR_LNCRNA_REFERENCE_DB")
    if env:
        return env
    local = os.path.join(_APP, "lncrna_reference.db")
    return local if os.path.exists(local) else "/app_v1/lncrna_reference.db"


def species_for_genome(genome):
    """Species key used by both tables: hg19/hg38 -> hsa_hg19/hsa_hg38, else the code."""
    if genome in ("hg19", "hg38"):
        return "hsa_" + genome
    return genome


class LncrnaSymbols(object):
    def __init__(self, species, db_path=None, tsv_path=None):
        self.species = species
        self.db_path = db_path or default_reference_db()
        self.tsv_path = tsv_path or default_symbols_path()
        self._gene_symbol = None  # lazily loaded {gene: symbol} for this species
        self._tx_gene = {}

    def _load_symbols(self):
        """The symbol table is keyed by organism ("hsa"), not assembly: ENSG ids
        are stable across hg19/hg38, and the GRCh37 dump carries no symbols."""
        self._gene_symbol = {}
        if not os.path.exists(self.tsv_path):
            return
        organism = self.species.split("_")[0]
        with open(self.tsv_path, newline="", encoding="utf-8") as f:
            for row in csv.DictReader(f, delimiter="\t"):
                if row.get("species") == organism and row.get("symbol"):
                    self._gene_symbol[normalize_lncrna_id(row["gene"])] = row["symbol"].strip()

    def gene_of(self, transcript_id):
        tx = normalize_lncrna_id(transcript_id)
        if tx in self._tx_gene:
            return self._tx_gene[tx]
        gene = None
        if os.path.exists(self.db_path):
            conn = sqlite3.connect(self.db_path)
            try:
                row = conn.execute(
                    "SELECT gene FROM lncrna_reference WHERE species = ? AND transcript_id = ? LIMIT 1",
                    (self.species, tx)).fetchone()
            finally:
                conn.close()
            gene = row[0] if row and row[0] else None
        self._tx_gene[tx] = gene
        return gene

    def symbol_for(self, lncrna_id):
        """HGNC symbol for a transcript or gene id, or None when unknown."""
        if self._gene_symbol is None:
            self._load_symbols()
        norm = normalize_lncrna_id(lncrna_id)
        if norm in self._gene_symbol:
            return self._gene_symbol[norm]
        gene = self.gene_of(norm)
        if gene and gene in self._gene_symbol:
            return self._gene_symbol[gene]
        return None
