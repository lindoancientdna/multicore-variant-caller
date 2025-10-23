#!/usr/bin/env python3
import pysam
import random
import argparse
import gzip
import os
from datetime import date
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed

def load_bed_sites(bed_file):
    """
    Load BED sites (plain, gzipped, or Tabix-indexed .bed.gz).
    Expected format:
      chrom  start  end  rsid:REF/ALT  label
    Returns dict: {chrom: {pos: {"rsid": str, "ref": "G", "alt": "A"}}}
    """
    bed_sites = {}

    if os.path.exists(bed_file + ".tbi"):
        tbx = pysam.TabixFile(bed_file)
        for chrom in tbx.contigs:
            for line in tbx.fetch(chrom):
                parts = line.strip().split()
                if len(parts) < 4:
                    continue
                chrom, start, end, rsinfo = parts[:4]
                start = int(start)
                pos = start + 1
                rsid, alleles = rsinfo.split(":") if ":" in rsinfo else (rsinfo, "N/N")
                ref, alt = alleles.split("/") if "/" in alleles else ("N", "N")
                bed_sites.setdefault(chrom, {})[pos] = {"rsid": rsid, "ref": ref, "alt": alt}
        tbx.close()
    else:
        opener = gzip.open if bed_file.endswith(".gz") else open
        with opener(bed_file, 'rt') as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.strip().split()
                if len(parts) < 4:
                    continue
                chrom, start, end, rsinfo = parts[:4]
                start = int(start)
                pos = start + 1
                rsid, alleles = rsinfo.split(":") if ":" in rsinfo else (rsinfo, "N/N")
                ref, alt = alleles.split("/") if "/" in alleles else ("N", "N")
                bed_sites.setdefault(chrom, {})[pos] = {"rsid": rsid, "ref": ref, "alt": alt}

    return bed_sites


def harmonize_chrom_names(bed_sites, bam_file, verbose=False):
    """Ensure chromosome naming (with or without 'chr' prefix) matches between BED and BAM."""
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
        harmonized[new_chrom] = bed_sites[chrom]

    if verbose and harmonized != bed_sites:
        print("[INFO] BED chromosome names harmonized to match BAM.")
    return harmonized


def call_pseudohaploid_per_sample(bam_file, chrom, pos, min_baseq=20, min_mapq=30):
    """Return a single pseudo-haploid basecall (A/C/G/T/N) for a given BAM and position."""
    bam = pysam.AlignmentFile(bam_file, "rb")
    basecall = "N"

    try:
        for pileup in bam.pileup(chrom, pos - 1, pos,
                                 truncate=True,
                                 min_base_quality=min_baseq,
                                 min_mapping_quality=min_mapq):
            if pileup.pos + 1 != pos:
                continue
            bases = []
            for read in pileup.pileups:
                if read.is_del or read.is_refskip:
                    continue
                base = read.alignment.query_sequence[read.query_position]
                if (read.alignment.mapping_quality >= min_mapq and
                    read.alignment.query_qualities[read.query_position] >= min_baseq):
                    bases.append(base.upper())
            if bases:
                basecall = random.choice(bases)
            break
    except ValueError:
        pass  # chromosome not found in BAM

    bam.close()
    return basecall


def process_chromosome(chrom, sites, bams, samples, min_baseq, min_mapq, verbose=False):
    """Process all sites for a single chromosome."""
    output_lines = []
    for pos, meta in tqdm(sorted(sites.items()),
                          desc=f"Chr {chrom}",
                          total=len(sites),
                          disable=not verbose,
                          dynamic_ncols=True):
        ref, alt = meta["ref"], meta["alt"]
        rsid = meta["rsid"]

        gts = []
        n_called = 0
        for bam_file in bams:
            base = call_pseudohaploid_per_sample(bam_file, chrom, pos, min_baseq, min_mapq)
            if base == ref:
                gt = "0"
                n_called += 1
            elif base == alt:
                gt = "1"
                n_called += 1
            else:
                gt = "."
            gts.append(gt)

        if n_called > 0:
            line = f"{chrom}\t{pos}\t{rsid}\t{ref}\t{alt}\t.\tPASS\t.\tGT\t" + "\t".join(gts)
            output_lines.append(line)
    return output_lines


def main():
    parser = argparse.ArgumentParser(
        description="Multi-sample pseudo-haploid caller (parallelized per chromosome)."
    )
    parser.add_argument("--bams", nargs="+", required=True,
                        help="Input BAM files (sorted, indexed). One per sample.")
    parser.add_argument("--samples", nargs="+", required=True,
                        help="Sample names (must match number of BAMs).")
    parser.add_argument("--bed", required=True,
                        help="BED file (plain, gz, or tabix-indexed .bed.gz)")
    parser.add_argument("--out", required=True,
                        help="Output VCF file")
    parser.add_argument("--cores", type=int, default=4,
                        help="Number of chromosomes to process in parallel.")
    parser.add_argument("--min-baseq", type=int, default=20,
                        help="Minimum base quality")
    parser.add_argument("--min-mapq", type=int, default=30,
                        help="Minimum mapping quality")
    parser.add_argument("--verbose", action="store_true",
                        help="Print progress and chromosome naming info.")
    args = parser.parse_args()

    if len(args.bams) != len(args.samples):
        raise ValueError("Number of BAM files must match number of sample names.")

    bed_sites = load_bed_sites(args.bed)
    bed_sites = harmonize_chrom_names(bed_sites, args.bams[0], verbose=args.verbose)
    total_sites = sum(len(bed_sites[c]) for c in bed_sites)
    if args.verbose:
        print(f"[INFO] Loaded {total_sites:,} target sites across {len(bed_sites)} chromosomes.")

    # --- Parallel execution ---
    results = []
    with ProcessPoolExecutor(max_workers=args.cores) as executor:
        future_to_chrom = {
            executor.submit(process_chromosome, chrom, bed_sites[chrom],
                            args.bams, args.samples, args.min_baseq, args.min_mapq, args.verbose): chrom
            for chrom in bed_sites
        }

        for future in as_completed(future_to_chrom):
            chrom = future_to_chrom[future]
            try:
                chrom_results = future.result()
                results.extend(chrom_results)
            except Exception as e:
                print(f"[ERROR] Failed processing {chrom}: {e}")

    # --- Write final VCF ---
    with open(args.out, "w") as out:
        out.write("##fileformat=VCFv4.2\n")
        out.write(f"##fileDate={date.today()}\n")
        out.write("##source=bam-caller_bedfiltered_multisample_vcf_parallel.py\n")
        out.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(args.samples) + "\n")
        for line in sorted(results):
            out.write(line + "\n")
            out.flush()

    if args.verbose:
        print(f"[INFO] Completed successfully. Wrote {len(results):,} variant sites to {args.out}.")


if __name__ == "__main__":
    main()