#!/usr/bin/env python3
"""Fetch the open data BioPathAgent is built on.

Sources, all free and redistributable for research:
  Rhea   (EMBL-EBI/SIB)  curated reactions + reaction SMILES + EC + xrefs
  KEGG   REST            pathway <-> reaction membership (whole-database LINK queries)
  RetroRules            SMARTS reaction templates at several radii

Deliberately NOT fetched: MetaCyc (needs a licence) and KEGG per-entry records
(the REST service is not a bulk endpoint). PLAN.md 6 explains the substitution:
KEGG supplies pathway membership only, and the chemistry comes from Rhea.

Usage
-----
    python preprocessing/fetch_data.py --what rhea kegg      # small, ~50 MB
    python preprocessing/fetch_data.py --what retrorules     # larger, asks first
    python preprocessing/fetch_data.py --what all --yes
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

import requests

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from biopathagent import config  # noqa: E402

# A generic client name. Never put the user's identity in an outbound header.
HEADERS = {"User-Agent": "BioPathAgent/0.1 (research pipeline; python-requests)"}
CHUNK = 1 << 20


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(CHUNK), b""):
            h.update(block)
    return h.hexdigest()


def download(url: str, dest: Path, *, timeout: int = 120, retries: int = 3) -> dict:
    """Stream a URL to dest. Returns a provenance record.

    Skips the download when dest already exists, so the script is restartable and
    cheap to re-run; delete the file to force a refresh.
    """
    if dest.exists() and dest.stat().st_size > 0:
        return {
            "url": url,
            "path": str(dest),
            "bytes": dest.stat().st_size,
            "sha256": _sha256(dest),
            "status": "cached",
        }

    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    last_err: Exception | None = None

    for attempt in range(1, retries + 1):
        try:
            with requests.get(url, headers=HEADERS, stream=True, timeout=timeout) as r:
                r.raise_for_status()
                with tmp.open("wb") as fh:
                    for chunk in r.iter_content(CHUNK):
                        if chunk:
                            fh.write(chunk)
            tmp.replace(dest)
            return {
                "url": url,
                "path": str(dest),
                "bytes": dest.stat().st_size,
                "sha256": _sha256(dest),
                "status": "downloaded",
                "fetched_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            }
        except Exception as exc:  # noqa: BLE001 - we retry on anything transient
            last_err = exc
            tmp.unlink(missing_ok=True)
            if attempt < retries:
                time.sleep(2 ** attempt)

    return {"url": url, "path": str(dest), "status": "FAILED", "error": repr(last_err)}


def fetch_rhea(outdir: Path) -> list[dict]:
    records = []
    for name in config.RHEA_FILES:
        rec = download(f"{config.RHEA_TSV_BASE}/{name}", outdir / "rhea" / name)
        rec["role"] = config.RHEA_FILES[name]
        print(f"  rhea/{name:<34} {rec['status']:<11} {rec.get('bytes', 0):>12,} B")
        records.append(rec)
    return records


def fetch_kegg(outdir: Path) -> list[dict]:
    records = []
    for name, url in config.KEGG_LINKS.items():
        rec = download(url, outdir / "kegg" / name)
        print(f"  kegg/{name:<34} {rec['status']:<11} {rec.get('bytes', 0):>12,} B")
        records.append(rec)
        time.sleep(1)  # courtesy to a free academic service
    return records


def fetch_retrorules(outdir: Path) -> list[dict]:
    """RetroRules has moved hosts; try the known locations in order."""
    for url in config.RETRORULES_CANDIDATES:
        suffix = ".tar.gz" if url.endswith(".tar.gz") else ".tar.gz"
        dest = outdir / "retrorules" / ("retrorules_rr02_rp2_hs" + suffix)
        rec = download(url, dest, timeout=600)
        if rec["status"] != "FAILED":
            print(f"  retrorules from {url}\n    {rec['status']} {rec.get('bytes', 0):,} B")
            return [rec]
        print(f"  retrorules: {url} failed ({rec.get('error', '')[:80]})")
    return [{"url": "retrorules", "status": "FAILED", "error": "all candidate URLs failed"}]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--what", nargs="+", default=["rhea", "kegg"],
                    choices=["rhea", "kegg", "retrorules", "all"])
    ap.add_argument("--outdir", type=Path, default=None,
                    help="default: config.DATA_RAW")
    ap.add_argument("--yes", action="store_true",
                    help="do not prompt before the larger RetroRules download")
    args = ap.parse_args()

    outdir = args.outdir or config.DATA_RAW
    config.ensure_dirs()
    outdir.mkdir(parents=True, exist_ok=True)

    what = set(args.what)
    if "all" in what:
        what = {"rhea", "kegg", "retrorules"}

    print(f"destination: {outdir}\n")
    provenance: list[dict] = []

    if "rhea" in what:
        print("Rhea (EMBL-EBI/SIB) ...")
        provenance += fetch_rhea(outdir)
    if "kegg" in what:
        print("KEGG REST (whole-database LINK queries only) ...")
        provenance += fetch_kegg(outdir)
    if "retrorules" in what:
        if not args.yes:
            print("\nRetroRules is a few hundred MB. Re-run with --yes to fetch it.")
        else:
            print("RetroRules ...")
            provenance += fetch_retrorules(outdir)

    manifest = outdir / "provenance.json"
    old = json.loads(manifest.read_text()) if manifest.exists() else []
    by_url = {r["url"]: r for r in old}
    by_url.update({r["url"]: r for r in provenance})
    manifest.write_text(json.dumps(sorted(by_url.values(), key=lambda r: r["url"]), indent=2))

    failed = [r for r in provenance if r.get("status") == "FAILED"]
    print(f"\nmanifest: {manifest}")
    if failed:
        print(f"FAILED: {len(failed)}")
        for r in failed:
            print(f"  {r['url']}  {r.get('error', '')[:120]}")
        return 1
    print("all requested sources present")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
