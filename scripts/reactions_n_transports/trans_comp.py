import re
import sqlite3
import argparse
import shutil
import logging
import sys
import requests

# Database for compund names for the transport reactions
CREATE_TABLE = """
CREATE TABLE IF NOT EXISTS transport_compound_map (
    go_id            TEXT NOT NULL,
    go_name          TEXT NOT NULL,
    substrate_guess  TEXT,
    compound_id      TEXT,
    compound_name    TEXT,
    match_tier       TEXT,
    ambiguous        INTEGER NOT NULL DEFAULT 0,
    ambiguous_reason TEXT,
    chebi_id         TEXT,
    PRIMARY KEY (go_id),
    FOREIGN KEY (compound_id) REFERENCES compounds(compound_id)
)
"""

CREATE_CHEBI_CACHE = """
CREATE TABLE IF NOT EXISTS chebi_synonym_cache (
    normalized_name  TEXT PRIMARY KEY,
    compound_id      TEXT NOT NULL
)
"""

# Suffix and prefix list

_TRANSPORT_SUFFIXES = [
    " transmembrane transporter activity",
    " transporter activity",
    " transmembrane transport",
    " transport",
    " uptake",
    " efflux",
    " import",
    " export",
    " secretion",
    " absorption",
    " import across plasma membrane",
    " export across plasma membrane",
    " import into cell",
    " export from cell",
    " import across the plasma membrane",
    " export across the plasma membrane",
    " floppase activity",
    " transfer activity",
    " translocation",
    " releasing activity",
    " oxidase activity",
    " synthase activity, rotational mechanism",
    " atpase activity, rotational mechanism",
    " phosphotransferase activity",
    " phosphotransferase system transporter activity",
    " phosphotransferase system",
    " oxidoreductase activity",
    " dehydrogenase activity",
    " transhydrogenase activity",
    " transfer",
    "-transporter activity",
    " by the type iii secretion system",
    " dehydrogenase (ubiquinone) activity",
    " involved in transformation",
]

_TRANSPORT_PREFIXES = [
    "ABC-type ",
    "P-type ",
    "ATPase-coupled ",
    "proton-coupled ",
    "sodium-coupled ",
    "energy-coupled ",
    "ATP-dependent ",
    "proton-translocating ",
    "oxidoreduction-driven active ",
    "protein-n(pi)-phosphohistidine-",
    "protein-phosphocysteine-",
]

_BLOCKLIST_SUBSTRATES = {
    "dna import into cell",
    "protein secretion",
    "intermembrane phospholipid",
    "lipopolysaccharide",
    "phospholipid",

    "transmembrane",
    "atpase-coupled",
    "abc-type",
    "p-type ion",
    "oxidoreduction-driven active",
    "monoatomic cation",
    "transition metal ion",
    "metal ion",
    "carbohydrate",
    "lipid",
    "xenobiotic",
    "polyamine",
    "peptide",
    "solute",
    "ion",
    "sulfur compound",
    "quaternary ammonium group",
    "electron transport coupled proton",
    "fructose phosphotransferase system",
    "mannose phosphotransferase system",
    "sorbose phosphotransferase system",
    "d-fructose-phosphotransferase system",
    "galactitol-phosphotransferase system",
    "l-ascorbate-phosphotransferase system",
    "n,n'-diacetylchitobiose phosphotransferase system",
    "sugar phosphotransferase system",
    "cytochrome-c",
    "cytochrome bo3 ubiquinol",
    "nadh",
    "dna",
    "proton-transporting atp",
    "proton-transporting atpase",
    "proton-transporting",
    "lipoprotein",
    "protein secretion by the type iii secretion system",
    "polar amino acid",
    "amino acid",
    "l-alpha-amino acid",
    "hexose",
    "monosaccharide",
    "alkanesulfonate transporter",
    "taurine transporter",
    "siderophore uptake",
    "siderophore",
    "siderophore-iron",
    "xenobiotic detoxification by transmembrane",
    "phosphoenolpyruvate-dependent sugar",
    "amino-acid betaine",
    "autoinducer ai-2",
    "nad(p)+",
    "ferric hydroxamate",
    "ferric-hydroxamate",
    "lipid-linked peptidoglycan",
    "monovalent copper",
    "organic phosphonate",
    "polymyxin",
}

# Manuel liste
_MANUAL_SUBSTRATE_MAP = {
    # ions not in DB by common name
    "zinc":             "C00038",   # Zinc cation
    "zinc ion":         "C00038",
    "magnesium":        "C00305",   # Magnesium cation
    "magnesium ion":    "C00305",
    "nickel":           "C00291",   # Nickel
    "nickel cation":    "C00291",
    "cadmium":          "C00076",   # Cadmium
    "cadmium ion":      "C00076",
    "copper ion":       "C00070",   # Copper
    "potassium":        "C00238",   # Potassium
    "potassium ion":    "C00238",
    "lead ion":         None,       # Not in standard KEGG metabolic models
    # sugars
    "maltose":          "C00208",
    "methylgalactoside":"C01353",
    "n-acetylglucosamine": "C00140",
    # other
    "phosphonate":      "C06701",
}
_BLOCKLISTED = object()

# Extractor
def extract_substrate(go_name: str) -> str | None | object:
    name = go_name.lower().strip()

    matched = False
    for suffix in sorted(_TRANSPORT_SUFFIXES, key=len, reverse=True):
        if name.endswith(suffix.lower()):
            name = name[: -len(suffix)].strip()
            matched = True
            break

    if not matched:
        return None

    for prefix in sorted(_TRANSPORT_PREFIXES, key=len, reverse=True):
        if name.startswith(prefix.lower()):
            name = name[len(prefix):].strip()
            break

    if not name:
        return None

    if name in _BLOCKLIST_SUBSTRATES:
        log.debug(f"  BLOCKLIST  {go_name!r} → {name!r}")
        return _BLOCKLISTED

    return name

# Normalizer for the compound name
def normalize_name(name: str) -> str:
    name = re.sub(r"<[^>]+>", "", name)   # strip HTML tags from ChEBI
    name = name.lower().strip()
    name = re.sub(r"^(alpha|beta|gamma|delta)-", "", name)
    name = re.sub(r"^[dlsn]-(?=[a-z])", "", name)
    name = re.sub(r"\s+(ion|cation|anion)$", "", name)
    name = re.sub(r"\s*\(\d+[\+\-]\)$", "", name)
    return name.strip()

# --- ChEBI functions ---
class ChebiClient:
    MAX_CONSECUTIVE_ERRORS = 3
    BASE_URL = "https://www.ebi.ac.uk/chebi/backend/api/public"

    def __init__(self):
        self.consecutive_errors = 0
        self.network_disabled = False
        self.errors = []

    def _record_error(self, label: str, reason: str):
        self.consecutive_errors += 1

        self.errors.append({
            "accession": label,
            "reason": reason,
        })

        if self.consecutive_errors >= self.MAX_CONSECUTIVE_ERRORS:
            self.network_disabled = True

    def _success(self):
        self.consecutive_errors = 0

    def get_compound_detail(self, chebi_accession: str) -> dict | None:
        if self.network_disabled:
            return None

        try:
            r = requests.get(
                f"{self.BASE_URL}/compound/{chebi_accession}/",
                timeout=5,
            )
            r.raise_for_status()
            self._success()
            return r.json()

        except requests.exceptions.Timeout:
            self._record_error(chebi_accession, "timeout")

        except requests.exceptions.HTTPError as e:
            status = e.response.status_code if e.response else "unknown"
            self._record_error(
                chebi_accession,
                f"HTTP {status}",
            )

        except requests.exceptions.ConnectionError:
            self._record_error(
                chebi_accession,
                "connection error",
            )

        except Exception as e:
            self._record_error(
                chebi_accession,
                type(e).__name__,
            )

        return None

def build_chebi_kegg_map(cur: sqlite3.Cursor,) -> tuple[dict[str, str], dict]:
    """
    Fetch ChEBI compound details and build a normalized-name → KEGG map.

    Returns:
        (
            {normalized_name: kegg_compound_id},
            {
                "found": number of successful API lookups,
                "failed": number of failed API lookups,
                "errors": list of network errors,
                "network_disabled": whether enrichment was stopped,
            }
        )
    """

    cur.execute("""
        SELECT cx.compound_id, cx.external_id
        FROM compound_xrefs cx
        WHERE cx.namespace = 'chebi'
    """)

    chebi_rows = cur.fetchall()
    client = ChebiClient()
    enriched = {}
    found = 0
    failed = 0

    for compound_id, chebi_accession in chebi_rows:

        if client.network_disabled:
            break

        detail = client.get_compound_detail(chebi_accession)

        if not detail:
            failed += 1
            continue

        found += 1

        # Collect all name strings from all name types
        name_strings = []

        for name_type, name_list in detail.get("names", {}).items():
            for entry in name_list:
                raw = entry.get("ascii_name") or entry.get("name", "")
                clean = re.sub(r"<[^>]+>", "", raw).strip()

                if clean:
                    name_strings.append(clean)

        # Also add top-level name and ascii_name
        for field in ("name", "ascii_name"):
            val = detail.get(field, "")
            if val:
                name_strings.append(val)

        # Get KEGG C-number from MANUAL_X_REF.
        # Fall back to our own compound_id.
        kegg_id = compound_id

        for xref in detail.get(
            "database_accessions", {}
        ).get("MANUAL_X_REF", []):

            if (
                xref.get("prefix") == "kegg.compound"
                and xref.get("accession_number", "").startswith("C")
            ):
                kegg_id = xref["accession_number"]
                break

        for name in name_strings:
            norm = normalize_name(name)

            if norm and len(norm) > 2 and norm not in enriched:
                enriched[norm] = kegg_id

    chebi_stats = {
        "found": found,
        "failed": failed,
        "errors": client.errors,
        "network_disabled": client.network_disabled,
        "mappings": len(enriched),
        "total": len(chebi_rows),
    }

    return enriched, chebi_stats

def load_or_build_chebi_cache(cur: sqlite3.Cursor, con: sqlite3.Connection) -> tuple[dict[str, str], dict]:
    cur.execute("SELECT COUNT(*) FROM chebi_synonym_cache")
    count = cur.fetchone()[0]

    if count > 0:
        cur.execute("""
            SELECT normalized_name, compound_id
            FROM chebi_synonym_cache
        """)

        cached = dict(cur.fetchall())

        return cached, {
            "total": 0,
            "found": 0,
            "failed": 0,
            "errors": [],
            "network_disabled": False,
            "mappings": len(cached),
            "from_cache": True,
        }

    enriched, chebi_stats = build_chebi_kegg_map(cur)

    cur.executemany("""
        INSERT OR IGNORE INTO chebi_synonym_cache
            (normalized_name, compound_id)
        VALUES (?, ?)
    """, enriched.items())

    con.commit()

    return enriched, chebi_stats

def build_synonym_map(cur: sqlite3.Cursor) -> tuple[dict[str, str], dict[str, str]]:
    cur.execute("SELECT compound_id, name FROM compounds WHERE name IS NOT NULL")
    rows = cur.fetchall()

    synonym_map:    dict[str, str] = {}
    normalized_map: dict[str, str] = {}

    for compound_id, name in rows:
        synonym_map[name.lower()] = compound_id
        norm = normalize_name(name)
        if norm and norm not in normalized_map:
            normalized_map[norm] = compound_id

    cur.execute("SELECT compound_id, external_id FROM compound_xrefs WHERE namespace = 'bigg'")
    for compound_id, ext_id in cur.fetchall():
        readable = ext_id.lower().replace("__", " ").replace("_", " ").strip()
        if readable and len(readable) > 2:
            synonym_map.setdefault(readable, compound_id)
            norm = normalize_name(readable)
            if norm:
                normalized_map.setdefault(norm, compound_id)
        synonym_map.setdefault(ext_id.lower(), compound_id)

    cur.execute("""
        SELECT compound_id, external_id FROM compound_xrefs
        WHERE namespace IN ('hmdb', 'mnx')
    """)
    for compound_id, external_id in cur.fetchall():
        cleaned = external_id.lower().replace("_", " ").replace("-", " ").strip()
        if len(cleaned) > 4:
            synonym_map.setdefault(cleaned, compound_id)
            norm = normalize_name(cleaned)
            if norm:
                normalized_map.setdefault(norm, compound_id)

    log.info(f"  Local synonym map: {len(synonym_map)} entries, {len(normalized_map)} normalized.")
    return synonym_map, normalized_map

def resolve_local(
    substrate: str,
    synonym_map: dict[str, str],
    normalized_map: dict[str, str],
) -> tuple[str, str] | None:
    # tier 1 — exact
    if substrate in synonym_map:
        return synonym_map[substrate], "exact"

    # tier 2 — normalized
    norm = normalize_name(substrate)
    if norm in normalized_map:
        return normalized_map[norm], "normalized"

    # tier 3 — whole-word regex on normalized map
    pattern = re.compile(
        r"(^|(?<=\s))" + re.escape(norm) +
        r"(?!\s*\d)(?!-phosphate)(?!-sulfate)(?!-diphosphate)(?![\w-])",
        re.IGNORECASE,
    )
    matches = {
        cid for norm_name, cid in normalized_map.items()
        if pattern.search(norm_name)
    }
    if len(matches) == 1:
        return next(iter(matches)), "regex"
    elif len(matches) > 1:
        log.debug(f"  Regex ambiguous for {substrate!r}: {len(matches)} candidates")

    return None

def build_transport_compound_map(db_path: str, reset: bool = False, reset_chebi: bool = False) -> None:
    log.info("=== Building transport_compound_map ===")

    con = sqlite3.connect(db_path)
    cur = con.cursor()

    if reset:
        log.warning("  --reset flag set: dropping transport_compound_map.")
        cur.execute("DROP TABLE IF EXISTS transport_compound_map")
        con.commit()

    if reset_chebi:
        log.warning("  --reset-chebi flag set: dropping chebi_synonym_cache.")
        cur.execute("DROP TABLE IF EXISTS chebi_synonym_cache")
        con.commit()

    cur.executescript(CREATE_TABLE)
    cur.executescript(CREATE_CHEBI_CACHE)
    con.commit()

    log.info("  Table transport_compound_map created/verified.")

    cur.execute("""
        SELECT DISTINCT go_id, go_name
        FROM transport_annotations
        WHERE go_aspect IN ('biological_process', 'molecular_function')
    """)
    go_terms = cur.fetchall()
    log.info(f"  Found {len(go_terms)} distinct GO terms to map.")

    # build local maps
    log.info("Building local synonym map...")
    synonym_map, normalized_map = build_synonym_map(cur)

    # enrich with ChEBI — one API call per compound, not per GO term
    log.info("Enriching via ChEBI compound details (one call per DB compound)...")
    chebi_enriched = load_or_build_chebi_cache(cur, con)
    added = 0
    for norm_name, kegg_id in chebi_enriched.items():
        if norm_name not in normalized_map:
            normalized_map[norm_name] = kegg_id
            added += 1
    log.info(f"  Merged {added} new entries from ChEBI into normalized map.")

    stats = {
        "exact": 0, "normalized": 0, "regex": 0,
        "chebi": 0, "manual":0, "ambiguous": 0, "no_suffix": 0,
    }

    for go_id, go_name in go_terms:
        cur.execute("SELECT 1 FROM transport_compound_map WHERE go_id = ?", (go_id,))
        if cur.fetchone():
            log.debug(f"  SKIP (already in table) {go_id}")
            continue

        substrate = extract_substrate(go_name)

        if substrate is None:
            reason = "GO name did not match any transport suffix"
            log.debug(f"  NO SUFFIX  {go_id}  {go_name!r}")
        elif substrate is _BLOCKLISTED:
            reason = "GO term describes mechanism, not substrate (blocklisted)"
            substrate = None
            log.debug(f"  BLOCKLISTED  {go_id}  {go_name!r}")
        else:
            reason = None

        if reason:
            cur.execute("""
                INSERT OR REPLACE INTO transport_compound_map
                    (go_id, go_name, substrate_guess, ambiguous, ambiguous_reason)
                VALUES (?, ?, NULL, 1, ?)
            """, (go_id, go_name, reason))
            stats["no_suffix"] += 1
            continue

        result = resolve_local(substrate, synonym_map, normalized_map)

        # Check explicit manual overrides first before testing general local rules
        if substrate in _MANUAL_SUBSTRATE_MAP:
            log.debug(f"  MANUAL HIT  {substrate!r}")
            kegg_id = _MANUAL_SUBSTRATE_MAP[substrate]
            if kegg_id:
                result = (kegg_id, "manual")
            else:
                cur.execute("""
                    INSERT OR REPLACE INTO transport_compound_map
                        (go_id, go_name, substrate_guess, ambiguous, ambiguous_reason)
                    VALUES (?, ?, ?, 1, 'known ambiguous compound class — no single KEGG mapping')
                """, (go_id, go_name, substrate))
                stats["ambiguous"] += 1
                log.info(f"  KNOWN-AMBIGUOUS  {go_id}  {substrate!r}")
                continue
        else:
            result = resolve_local(substrate, synonym_map, normalized_map)

        if result:
            compound_id, tier = result
            # if this came from ChEBI enrichment, label it as chebi tier
            norm = normalize_name(substrate)
            if tier == "normalized" and norm in chebi_enriched:
                tier = "chebi"
            stats[tier] += 1

            cur.execute("SELECT name FROM compounds WHERE compound_id = ?", (compound_id,))
            row = cur.fetchone()
            compound_name = row[0] if row else None

            cur.execute("""
                INSERT OR REPLACE INTO transport_compound_map
                    (go_id, go_name, substrate_guess, compound_id, compound_name,
                     match_tier, ambiguous)
                VALUES (?, ?, ?, ?, ?, ?, 0)
            """, (go_id, go_name, substrate, compound_id, compound_name, tier))
            log.info(f"  {tier.upper():10s}  {go_id}  {substrate!r} → {compound_id} ({compound_name})")
            continue

        # all tiers failed
        cur.execute("""
            INSERT OR REPLACE INTO transport_compound_map
                (go_id, go_name, substrate_guess, ambiguous, ambiguous_reason)
            VALUES (?, ?, ?, 1, 'no match after local + ChEBI lookup')
        """, (go_id, go_name, substrate))
        log.info(f"  AMBIGUOUS   {go_id}  {substrate!r}")
        stats["ambiguous"] += 1

    con.commit()
    con.close()

def write_report(
    report_path: str,
    stats: dict,
    chebi_stats: dict,
):
    with open(report_path, "w", encoding="utf-8") as f:

        f.write("=== Transport Compound Mapping Report ===\n\n")

        f.write("=== Mapping Results ===\n")
        f.write(f"  Exact      : {stats['exact']}\n")
        f.write(f"  Normalized : {stats['normalized']}\n")
        f.write(f"  Regex      : {stats['regex']}\n")
        f.write(f"  ChEBI      : {stats['chebi']}\n")
        f.write(f"  Manual     : {stats['manual']}\n")
        f.write(f"  No suffix  : {stats['no_suffix']}\n")
        f.write(f"  Ambiguous  : {stats['ambiguous']}\n")

        f.write("\n=== ChEBI Enrichment ===\n")

        if chebi_stats.get("from_cache"):
            f.write("  Source     : cache\n")
        else:
            f.write("  Source     : ChEBI API\n")
            f.write(f"  Total      : {chebi_stats['total']}\n")
            f.write(f"  Found      : {chebi_stats['found']}\n")
            f.write(f"  Failed     : {chebi_stats['failed']}\n")

        f.write(f"  Mappings   : {chebi_stats['mappings']}\n")

        if chebi_stats["network_disabled"]:
            f.write("  Network    : DISABLED after consecutive errors\n")
        else:
            f.write("  Network    : OK\n")

        if chebi_stats["errors"]:
            f.write("\n  ChEBI errors:\n")

            for error in chebi_stats["errors"]:
                f.write(
                    f"    {error['accession']}: "
                    f"{error['reason']}\n"
                )

        f.write("\n=== Review ===\n")
        f.write("Review ambiguous entries with:\n")
        f.write(
            "SELECT go_id, go_name, substrate_guess "
            "FROM transport_compound_map "
            "WHERE ambiguous = 1;\n"
        )



if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Map GO transport terms to KEGG compound IDs and write to transport_compound_map table."
    )
    parser.add_argument("--input_db", required=True)
    parser.add_argument("--output_db", required=True)
    args = parser.parse_args()

    build_transport_compound_map(args.input_db, args_output_db)