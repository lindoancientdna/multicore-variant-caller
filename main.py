#!/usr/bin/env python3
import argparse
import gzip
import os
import random
from datetime import date
from concurrent.futures import ProcessPoolExecutor, as_completed

import pysam
from tqdm import tqdm


# -------- BED loading --------
def load_bed_sites(bed_file):
    """
    Load BED sites (plain, gzipped, or Tabix-indexed .bed.gz).
    Expected format (1 bp intervals):
      chrom  start  end  rsid:REF/ALT  [label...]
    Returns: {chrom: {pos: {"rsid": str, "ref": "A", "alt": "G"}}}
    """
    bed_sites = {}

    if os.path.exists(bed_file + ".tbi"):
        # Tabix-indexed: iterate contigs efficiently
        tbx = pysam.TabixFile(bed_file)
        for chrom in tbx.contigs:
            for line in tbx.fetch(chrom):
                parts = line.strip().split()
                if len(parts) < 4:
                    continue
                c, start, end, rsinfo = parts[:4]
                pos = int(start) + 1  # BED 0-based -> VCF 1-based
                rsid, alleles = rsinfo.split(":") if ":" in rsinfo else (rsinfo, "N/N")
                ref, alt = alleles.split("/") if "/" in alleles else ("N", "N")
                bed_sites.setdefault(c, {})[pos] = {"rsid": rsid, "ref": ref, "alt": alt}
        tbx.close()
    else:
        opener = gzip.open if bed_file.endswith(".gz") else open
        with opener(bed_file, "rt") as f:
            for line in f:
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.strip().split()
                if len(parts) < 4:
                    continue
                c, start, end, rsinfo = parts[:4]
                pos = int(start) + 1
                rsid, alleles = rsinfo.split(":") if ":" in rsinfo else (rsinfo, "N/N")
                ref, alt = alleles.split("/") if "/" in alleles else ("N", "N")
                bed_sites.setdefault(c, {})[pos] = {"rsid": rsid, "ref": ref, "alt": alt}

    return bed_sites


# -------- Chrom-name harmonization --------
def harmonize_chrom_names(bed_sites, bam_file, verbose=False):
    """
    Ensure chromosome naming (with or without 'chr' prefix) matches between BED and BAM.
    Returns a harmonized copy of bed_sites (keys renamed if needed).
    """
    bam = pysam.AlignmentFile(bam_file, "rb")
    bam_contigs = set(bam.references)
    bam.close()

    bam_has_chr = any(c.startswith("chr") for c in bam_contigs)
    if verbose:
        print(f"[INFO] Detected BAM naming: {'chr-prefixed' if bam_has_chr else 'numeric-only'}")

    harmonized = {}
    for chrom in bed_sites:
        new_chrom = chrom
        if bam_has_chr and not chrom.startswith("chr"):
            new_chrom = "chr" + chrom
        elif not bam_has_chr and chrom.startswith("chr"):
            new_chrom = chrom.replace("chr", "", 1)
        harmonized.setdefault(new_chrom, {}).update(bed_sites[chrom])

    if verbose and harmonized.keys() != bed_sites.keys():
        print("[INFO] BED chromosome names harmonized to match BAM.")
    return harmonized


# -------- Pileup & calling --------
def call_pseudohaploid_per_sample(bam_file, chrom, pos, min_baseq=20, min_mapq=30):
    """
    Return a single pseudo-haploid basecall (A/C/G/T/N) for one BAM at one position.
    Randomly draws from reads passing mapping/base quality filters.
    """
    basecall = "N"
    bam = pysam.AlignmentFile(bam_file, "rb")
    try:
        for pileup in bam.pileup(
            chrom,
            pos - 1,
            pos,
            truncate=True,
            min_base_quality=min_baseq,
            min_mapping_quality=min_mapq,
        ):
            if pileup.pos + 1 != pos:
                continue
            bases = []
            for pr in pileup.pileups:
                if pr.is_del or pr.is_refskip:
                    continue
                qpos = pr.query_position
                if qpos is None:
                    continue
                aln = pr.alignment
                if aln.mapping_quality >= min_mapq:
                    # guard against missing query_qualities
                    if aln.query_qualities is not None and qpos < len(aln.query_qualities):
                        if aln.query_qualities[qpos] >= min_baseq:
                            base = aln.query_sequence[qpos]
                            bases.append(base.upper())
            if bases:
                basecall = random.choice(bases)
            break
    except ValueError:
        # Chromosome not found in BAM
        pass
    finally:
        bam.close()
    return basecall


# -------- Per-chromosome worker --------
def process_chromosome(chrom, sites, bams, min_baseq, min_mapq, verbose=False):
    """
    Process all sites for a single chromosome; return list of VCF data lines (without header).
    Each line is tab-separated: CHROM POS ID REF ALT . PASS . GT <sample GTs...>
    """
    lines = []
    iterator = sorted(sites.items())
    pbar = tqdm(iterator, desc=f"Chr {chrom}", total=len(sites), dynamic_ncols=True, disable=not verbose)
    for pos, meta in pbar:
        ref, alt, rsid = meta["ref"], meta["alt"], meta["rsid"]
        gts = []
        n_called = 0
        for bam_file in bams:
            base = call_pseudohaploid_per_sample(bam_file, chrom, pos, min_baseq, min_mapq)
            if base == ref:
                gts.append("0")
                n_called += 1
            elif base == alt:
                gts.append("1")
                n_called += 1
            else:
                gts.append(".")
        # Skip site entirely if no sample had a valid call
        if n_called == 0:
            continue
        lines.append(f"{chrom}\t{pos}\t{rsid}\t{ref}\t{alt}\t.\tPASS\t.\tGT\t" + "\t".join(gts))
    return lines


# -------- Sorting helpers --------
# Standard human order: 1..22, X (23), Y (24), MT (25)
def _chrom_rank(chrom_str):
    """
    Map chromosome name (with or without 'chr') to a sortable rank and cleaned name for ordering.
    """
    c = chrom_str
    if c.startswith("chr"):
        c = c[3:]
    c_upper = c.upper()
    # Normalize common mitochondrial aliases to MT for ordering (output stays original)
    if c_upper in {"MT", "M", "CHRMT"}:
        return 25, c
    if c_upper == "X":
        return 23, c
    if c_upper == "Y":
        return 24, c
    # numeric?
    try:
        n = int(c)
        if 1 <= n <= 22:
            return n, c
    except ValueError:
        pass
    # Everything else goes after MT in lexicographic bucket
    return 10_000, c

def vcf_sort_key(line):
    chrom, pos = line.split("\t", 2)[:2]
    rank, _ = _chrom_rank(chrom)
    return (rank, int(pos))


# -------- Main --------
def main():
    ap = argparse.ArgumentParser(
        description="Multi-sample pseudo-haploid caller restricted to exact BED sites (parallel + sorted VCF)."
    )
    ap.add_argument("--bams", nargs="+", required=True, help="Input BAMs (sorted, indexed). One per sample.")
    ap.add_argument("--samples", nargs="+", required=True, help="Sample names (must match number of BAMs).")
    ap.add_argument("--bed", required=True, help="BED file (plain, .gz, or tabix-indexed .bed.gz+.tbi).")
    ap.add_argument("--out", required=True, help="Output VCF file (sorted by 1..22,X,Y,MT).")
    ap.add_argument("--threads", type=int, default=4, help="# chromosomes to process in parallel.")
    ap.add_argument("--min-baseq", type=int, default=20, help="Minimum base (Phred) quality.")
    ap.add_argument("--min-mapq", type=int, default=30, help="Minimum mapping quality.")
    ap.add_argument("--verbose", action="store_true", help="Print progress and naming info.")
    args = ap.parse_args()

    if len(args.bams) != len(args.samples):
        raise ValueError("Number of BAMs must equal number of sample names.")

    # Load and harmonize BED sites
    bed_sites = load_bed_sites(args.bed)
    bed_sites = harmonize_chrom_names(bed_sites, args.bams[0], verbose=args.verbose)
    total_sites = sum(len(bed_sites[c]) for c in bed_sites)
    if args.verbose:
        print(f"[INFO] Loaded {total_sites:,} target sites across {len(bed_sites)} chromosomes.")

    # Parallel per chromosome
    results = []
    with ProcessPoolExecutor(max_workers=max(1, args.threads)) as ex:
        futures = {
            ex.submit(process_chromosome, chrom, bed_sites[chrom], args.bams, args.min_baseq, args.min_mapq, args.verbose): chrom
            for chrom in bed_sites
        }
        for f in as_completed(futures):
            chrom = futures[f]
            try:
                results.extend(f.result())
            except Exception as e:
                print(f"[ERROR] Failed processing {chrom}: {e}")

    # Sort all lines by chrom/pos in standard human order (1..22,X,Y,MT)
    results.sort(key=vcf_sort_key)

    # Write VCF (sorted)
    with open(args.out, "w") as out:
        out.write("##fileformat=VCFv4.1\n")
        out.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        out.write("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Number of high-quality bases\">\n")

        #Add standard contigs
        for chrom in list(map(str, range(1, 23))) + ["X", "Y"]:
            out.write(f"##contig=<ID={chrom}>\n")

        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(args.samples) + "\n")
        out.flush()
        for line in results:
            out.write(line + "\n")
            out.flush()

    if args.verbose:
        print(f"[INFO] Completed. Wrote {len(results):,} variant sites to {args.out} (sorted).")
        print("[INFO] Optional: compress & index for downstream tools:")
        print(f"       bgzip {args.out} && bcftools index {args.out}.gz")


if __name__ == "__main__":
    main()