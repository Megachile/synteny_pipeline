#!/usr/bin/env python3
"""
Phase 6d: Validate Exonerate-extracted proteins to flag likely false positives.

Adds post-Exonerate QC using intrinsic low-complexity heuristics and optional
domain HMM checks (hmmscan against Pfam or a custom HMM library).

Outputs a TSV of per-sequence QC metrics and optional FASTAs of flagged sequences.

Usage examples:
  python 06d_validate_exonerate_outputs.py \
    --fasta outputs/<FAMILY>/06_extracted_sequences/all_extracted_genes.faa \
    --output-dir outputs/<FAMILY>/06_extracted_sequences \
    --pfam-db /path/to/Pfam-A.hmm \
    --min-domain-coverage 0.35 --max-evalue 1e-5

If hmmscan is not available or --pfam-db is omitted, only intrinsic complexity
filters are applied.
"""

from pathlib import Path
import argparse
import shutil
import subprocess
import tempfile
import math
import csv
from typing import Dict, List, Tuple

from Bio import SeqIO


def shannon_entropy(seq: str) -> float:
    """Compute global Shannon entropy (bits) for amino-acid composition."""
    if not seq:
        return 0.0
    counts: Dict[str, int] = {}
    length = 0
    for aa in seq:
        if aa == '*' or aa == '-':
            # Ignore stop and gap in entropy
            continue
        counts[aa] = counts.get(aa, 0) + 1
        length += 1
    if length == 0:
        return 0.0
    ent = 0.0
    for c in counts.values():
        p = c / length
        ent -= p * math.log2(p)
    return ent


def windowed_entropy_fraction(seq: str, w: int = 30, thresh: float = 2.2) -> float:
    """Fraction of sliding windows with entropy below threshold (lower = better)."""
    if len(seq) < w:
        return 1.0 if shannon_entropy(seq) < thresh else 0.0
    low = 0
    total = 0
    for i in range(0, len(seq) - w + 1):
        window = seq[i : i + w]
        total += 1
        if shannon_entropy(window) < thresh:
            low += 1
    return low / total if total else 0.0


def longest_homopolymer(seq: str) -> int:
    """Return the longest run of identical residues (simple low-complexity proxy)."""
    if not seq:
        return 0
    best = 1
    cur = 1
    last = seq[0]
    for a in seq[1:]:
        if a == last:
            cur += 1
            best = max(best, cur)
        else:
            cur = 1
            last = a
    return best


def run_hmmscan(hmmscan_bin: str, hmm_db: Path, fasta: Path, domtblout: Path, cpus: int = 2) -> None:
    cmd = [
        hmmscan_bin,
        "--noali",
        "--cut_ga",
        "--domtblout",
        str(domtblout),
        str(hmm_db),
        str(fasta),
    ]
    env = dict(**dict(), **{"HMMER_NCPU": str(cpus)})
    subprocess.run(cmd, check=True, text=True, env=env)


def run_hmmbuild(hmmbuild_bin: str, aln_fasta: Path, out_hmm: Path) -> None:
    cmd = [hmmbuild_bin, str(out_hmm), str(aln_fasta)]
    subprocess.run(cmd, check=True, text=True)


def run_hmmsearch(hmmsearch_bin: str, profile_hmm: Path, fasta: Path, domtblout: Path, cpus: int = 2) -> None:
    cmd = [
        hmmsearch_bin,
        "--noali",
        "--domtblout",
        str(domtblout),
        str(profile_hmm),
        str(fasta),
    ]
    env = dict(**dict(), **{"HMMER_NCPU": str(cpus)})
    subprocess.run(cmd, check=True, text=True, env=env)


def parse_domtblout_hmmscan(domtblout: Path) -> Dict[str, List[Dict[str, str]]]:
    """Parse HMMER domtblout and return hits per query id.

    Returns a mapping: query_id -> list of hits with keys including:
      hmm, evalue, qlen, ali_from, ali_to
    """
    hits: Dict[str, List[Dict[str, str]]] = {}
    with open(domtblout) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.strip().split()
            # HMMER domtblout columns (hmmscan):
            # target name, target acc, tlen, query name, query acc, qlen, E-value, score, bias,
            # #, of, c-Evalue, i-Evalue, score, bias,
            # hmm from, hmm to, ali from, ali to, env from, env to, acc, description...
            try:
                hmm_name = parts[0]
                q_name = parts[3]
                q_len = int(parts[5])
                i_eval = float(parts[12])
                ali_from = int(parts[17])
                ali_to = int(parts[18])
            except Exception:
                continue
            rec = {
                "hmm": hmm_name,
                "qname": q_name,
                "qlen": q_len,
                "i_evalue": i_eval,
                "ali_from": ali_from,
                "ali_to": ali_to,
            }
            hits.setdefault(q_name, []).append(rec)
    return hits


def parse_domtblout_hmmsearch(domtblout: Path) -> Dict[str, List[Dict[str, str]]]:
    """Parse HMMER domtblout from hmmsearch; key by target (sequence) id.

    Returns mapping: target_id -> list of hits with keys including:
      hmm, tname, tlen, i_evalue, ali_from, ali_to
    """
    hits: Dict[str, List[Dict[str, str]]] = {}
    with open(domtblout) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.strip().split()
            try:
                t_name = parts[0]
                t_len = int(parts[2])
                hmm_name = parts[3]
                i_eval = float(parts[12])
                ali_from = int(parts[17])
                ali_to = int(parts[18])
            except Exception:
                continue
            rec = {
                "hmm": hmm_name,
                "tname": t_name,
                "tlen": t_len,
                "i_evalue": i_eval,
                "ali_from": ali_from,
                "ali_to": ali_to,
            }
            hits.setdefault(t_name, []).append(rec)
    return hits


def compute_best_domain_metrics(hits: List[Dict[str, str]]) -> Tuple[float, float]:
    """Return (best_evalue, best_coverage) across domains for a sequence."""
    if not hits:
        return (1.0, 0.0)
    best_e = min(h["i_evalue"] for h in hits)
    best_cov = 0.0
    for h in hits:
        qlen = max(1, int(h["qlen"]))
        cov = (int(h["ali_to"]) - int(h["ali_from"]) + 1) / qlen
        if cov > best_cov:
            best_cov = cov
    return (best_e, best_cov)


def validate_sequences(
    fasta: Path,
    output_dir: Path,
    hmm_db: Path = None,
    hmmscan_bin: str = "hmmscan",
    profile_from_fasta: Path = None,
    mafft_bin: str = "mafft",
    hmmbuild_bin: str = "hmmbuild",
    hmmsearch_bin: str = "hmmsearch",
    min_len: int = 30,
    max_stop: int = 0,
    entropy_thresh: float = 2.4,
    low_entropy_window: int = 30,
    low_entropy_frac_thresh: float = 0.5,
    homopolymer_thresh: int = 12,
    max_evalue: float = 1e-5,
    min_domain_cov: float = 0.35,
    cpus: int = 2,
    profile_max_evalue: float = 1e-5,
    profile_min_target_cov: float = 0.35,
) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    qc_tsv = output_dir / "exonerate_validation.tsv"

    # Write sequences to a temp for hmmscan batch
    tmp_fasta = None
    domtblout = None
    hits_by_query: Dict[str, List[Dict[str, str]]] = {}

    do_hmm = False
    if hmm_db is not None:
        if shutil.which(hmmscan_bin):
            do_hmm = True
        else:
            print(f"WARNING: hmmscan not found in PATH; skipping HMM domain checks")

    profile_hits_by_target: Dict[str, List[Dict[str, str]]] = {}
    do_profile = False
    if profile_from_fasta is not None:
        if not profile_from_fasta.exists():
            print(f"WARNING: profile_from_fasta not found: {profile_from_fasta}; skipping profile HMM")
        elif not shutil.which(hmmbuild_bin) or not shutil.which(hmmsearch_bin):
            print("WARNING: hmmbuild/hmmsearch not found in PATH; skipping profile HMM")
        else:
            do_profile = True

    if do_hmm:
        tmp_fasta = Path(tempfile.mkstemp(prefix="hmmscan_q_")[1])
        with open(tmp_fasta, "w") as out_fa:
            for rec in SeqIO.parse(str(fasta), "fasta"):
                SeqIO.write(rec, out_fa, "fasta")
        domtblout = Path(tempfile.mkstemp(prefix="hmmscan_domtbl_")[1])
        run_hmmscan(hmmscan_bin, hmm_db, tmp_fasta, domtblout, cpus=cpus)
        hits_by_query = parse_domtblout_hmmscan(domtblout)

    if do_profile:
        # Build alignment if 2+ sequences and mafft available; else use raw as alignment
        seqs = list(SeqIO.parse(str(profile_from_fasta), "fasta"))
        if len(seqs) == 0:
            print(f"WARNING: profile_from_fasta contained 0 sequences; skipping profile HMM")
            do_profile = False
        else:
            tmp_aln = Path(tempfile.mkstemp(prefix="profile_aln_")[1])
            tmp_src = Path(tempfile.mkstemp(prefix="profile_src_")[1])
            SeqIO.write(seqs, tmp_src, "fasta")
            if len(seqs) >= 2 and shutil.which(mafft_bin):
                # mafft --auto src > aln
                with open(tmp_aln, "w") as out_aln:
                    subprocess.run([mafft_bin, "--auto", str(tmp_src)], check=True, text=True, stdout=out_aln)
            else:
                # Single-seq alignment fallback
                shutil.copyfile(tmp_src, tmp_aln)

            profile_hmm = Path(tempfile.mkstemp(prefix="profile_")[1])
            run_hmmbuild(hmmbuild_bin, tmp_aln, profile_hmm)

            domtbl_profile = Path(tempfile.mkstemp(prefix="hmmsearch_domtbl_")[1])
            run_hmmsearch(hmmsearch_bin, profile_hmm, fasta, domtbl_profile, cpus=cpus)
            profile_hits_by_target = parse_domtblout_hmmsearch(domtbl_profile)

    rows = []
    failed_records = []
    lcr_records = []

    for rec in SeqIO.parse(str(fasta), "fasta"):
        seq = str(rec.seq)
        length = len(seq)
        stop_codons = seq.count("*")
        ent = shannon_entropy(seq)
        low_frac = windowed_entropy_fraction(seq, w=low_entropy_window, thresh=entropy_thresh)
        hom = longest_homopolymer(seq)

        hmm_best_eval = None
        hmm_best_cov = None
        if do_hmm:
            hits = hits_by_query.get(rec.id, [])
            be, bc = compute_best_domain_metrics(hits)
            hmm_best_eval, hmm_best_cov = be, bc

        profile_best_eval = None
        profile_best_cov = None
        if do_profile:
            phits = profile_hits_by_target.get(rec.id, [])
            if phits:
                # Coverage on target length
                best_e = min(h["i_evalue"] for h in phits)
                best_cov = 0.0
                for h in phits:
                    tlen = max(1, int(h["tlen"]))
                    cov = (int(h["ali_to"]) - int(h["ali_from"]) + 1) / tlen
                    if cov > best_cov:
                        best_cov = cov
                profile_best_eval, profile_best_cov = best_e, best_cov

        low_complexity = (ent < entropy_thresh) or (low_frac > low_entropy_frac_thresh) or (hom >= homopolymer_thresh)

        domain_ok = True
        if do_hmm:
            domain_ok = (hmm_best_eval is not None) and (hmm_best_eval <= max_evalue) and (hmm_best_cov is not None) and (hmm_best_cov >= min_domain_cov)

        profile_ok = True
        if do_profile:
            profile_ok = (profile_best_eval is not None) and (profile_best_eval <= profile_max_evalue) and (profile_best_cov is not None) and (profile_best_cov >= profile_min_target_cov)

        passes = (
            (length >= min_len)
            and (stop_codons <= max_stop)
            and (not low_complexity)
        )

        rows.append({
            "id": rec.id,
            "length": length,
            "stop_codons": stop_codons,
            "entropy": round(ent, 3),
            "low_entropy_frac": round(low_frac, 3),
            "longest_homopolymer": hom,
            "low_complexity": low_complexity,
            "domain_ok": domain_ok if do_hmm else "NA",
            "pfam_best_evalue": ("{:.2e}".format(hmm_best_eval) if hmm_best_eval is not None else ""),
            "pfam_best_cov": (round(hmm_best_cov, 3) if hmm_best_cov is not None else ""),
            "profile_best_evalue": ("{:.2e}".format(profile_best_eval) if profile_best_eval is not None else ""),
            "profile_best_cov": (round(profile_best_cov, 3) if profile_best_cov is not None else ""),
            "profile_ok": profile_ok if do_profile else "NA",
            "status": "pass" if passes else "flag",
        })

        if not passes:
            failed_records.append(rec)
        elif low_complexity:
            lcr_records.append(rec)

    # Save TSV
    with open(qc_tsv, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "id",
                "length",
                "stop_codons",
                "entropy",
                "low_entropy_frac",
                "longest_homopolymer",
                "pfam_best_evalue",
                "pfam_best_cov",
                "profile_best_evalue",
                "profile_best_cov",
                "low_complexity",
                "domain_ok",
                "profile_ok",
                "status",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)

    # Save flagged sets for quick inspection
    if failed_records:
        SeqIO.write(failed_records, output_dir / "exonerate_validation_failed.faa", "fasta")
    if lcr_records:
        SeqIO.write(lcr_records, output_dir / "exonerate_low_complexity.faa", "fasta")

    # Cleanup temps
    if domtblout and domtblout.exists():
        domtblout.unlink(missing_ok=True)  # type: ignore[arg-type]
    if tmp_fasta and Path(tmp_fasta).exists():
        Path(tmp_fasta).unlink(missing_ok=True)

    return qc_tsv


def parse_args():
    p = argparse.ArgumentParser(description="Validate Exonerate outputs with low-complexity and optional HMM domain checks")
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--fasta", type=Path, help="Protein FASTA to validate (e.g., all_extracted_genes.faa)")
    g.add_argument("--extracted-dir", type=Path, help="06_extracted_sequences directory; uses all_extracted_genes.faa inside")

    p.add_argument("--output-dir", type=Path, help="Directory to write validation outputs", required=False)
    p.add_argument("--pfam-db", type=Path, help="Pfam-A.hmm (or custom HMM library) for hmmscan", required=False)
    p.add_argument("--hmmscan", type=str, default="hmmscan", help="hmmscan binary (default: hmmscan in PATH)")
    p.add_argument("--threads", type=int, default=2, help="CPUs for hmmscan")

    p.add_argument("--profile-from", type=Path, help="FASTA of family input/query proteins to build an HMM profile (hmmbuild+hmmsearch)", required=False)
    p.add_argument("--mafft", type=str, default="mafft", help="mafft binary (default: mafft in PATH)")
    p.add_argument("--hmmbuild", type=str, default="hmmbuild", help="hmmbuild binary")
    p.add_argument("--hmmsearch", type=str, default="hmmsearch", help="hmmsearch binary")

    p.add_argument("--min-length", type=int, default=30, help="Minimum protein length")
    p.add_argument("--max-stop", type=int, default=0, help="Maximum internal stop codons allowed")
    p.add_argument("--entropy-thresh", type=float, default=2.4, help="Global/window entropy threshold (bits)")
    p.add_argument("--low-entropy-window", type=int, default=30, help="Window size for entropy fraction")
    p.add_argument("--low-entropy-frac", type=float, default=0.5, help="Max fraction of low-entropy windows allowed")
    p.add_argument("--homopolymer-thresh", type=int, default=12, help="Max allowed homopolymer length")
    p.add_argument("--max-evalue", type=float, default=1e-5, help="Domain i-Evalue cutoff for hmmscan")
    p.add_argument("--min-domain-coverage", type=float, default=0.35, help="Minimum domain coverage on query (0-1)")

    return p.parse_args()


def main():
    args = parse_args()

    if args.extracted_dir:
        fasta = args.extracted_dir / "all_extracted_genes.faa"
        if not fasta.exists():
            raise SystemExit(f"Cannot find FASTA: {fasta}")
        output_dir = args.output_dir or args.extracted_dir
    else:
        fasta = args.fasta
        if not fasta.exists():
            raise SystemExit(f"Cannot find FASTA: {fasta}")
        output_dir = args.output_dir or fasta.parent

    hmm_db = args.pfam_db if args.pfam_db else None

    print("=" * 80)
    print("PHASE 6d: VALIDATE EXONERATE OUTPUTS")
    print("=" * 80)
    print(f"FASTA: {fasta}")
    print(f"Output: {output_dir}")
    if hmm_db:
        print(f"HMM DB: {hmm_db}")
        print(f"hmmscan: {args.hmmscan}")
    else:
        print("HMM DB: [none] (domain checks disabled)")
    if args.profile_from:
        print(f"Profile from: {args.profile_from}")
        print(f"mafft: {args.mafft} | hmmbuild: {args.hmmbuild} | hmmsearch: {args.hmmsearch}")

    qc_path = validate_sequences(
        fasta=fasta,
        output_dir=output_dir,
        hmm_db=hmm_db,
        hmmscan_bin=args.hmmscan,
        profile_from_fasta=args.profile_from,
        mafft_bin=args.mafft,
        hmmbuild_bin=args.hmmbuild,
        hmmsearch_bin=args.hmmsearch,
        min_len=args.min_length,
        max_stop=args.max_stop,
        entropy_thresh=args.entropy_thresh,
        low_entropy_window=args.low_entropy_window,
        low_entropy_frac_thresh=args.low_entropy_frac,
        homopolymer_thresh=args.homopolymer_thresh,
        max_evalue=args.max_evalue,
        min_domain_cov=args.min_domain_coverage,
        cpus=args.threads,
        profile_max_evalue=args.max_evalue,
        profile_min_target_cov=args.min_domain_coverage,
    )

    print(f"\nValidation written: {qc_path}")
    print("Done.")


if __name__ == "__main__":
    main()
