# Multicore Variant Caller  
**Parallel pseudo-haploid variant caller for ancient DNA**

Created and maintained by the Lindo Ancient DNA Lab at Emory University

## Overview  

`multicore-variant-caller` is a **high-performance, multicore pseudo-haploid SNP caller** for ancient DNA data.  
It efficiently processes large numbers of BAM files against dense SNP panels (e.g. 1240k) using **per-chromosome parallelization**.  
Each process samples one allele per SNP per individual, generating pseudo-haploid VCFs ideal for low-coverage aDNA analyses.

### Features
- Multi-core processing via `ProcessPoolExecutor`
- Accepts multiple BAM files and sample names  
- Reads BED or Tabix-indexed SNP panels  
- Automatically harmonizes chromosome naming 
- Pseudo-haploid base calling with adjustable quality filters  
- Outputs clean multi-sample VCF  
- Scales cores for efficient genome-wide processing  

---

## Dependencies

This tool requires **Python 3.8+** and a few external libraries.

### Required packages

| Library | Purpose | Install via pip |
|----------|----------|----------------|
| `pysam` | Reading and parsing BAM/CRAM files; pileup iteration | `pip install pysam` |
| `tqdm` | Visual progress bars for chromosome/site processing | `pip install tqdm` |

### Built-in Python modules (no installation needed)

These are included in the standard Python library:
- `argparse`
- `gzip`
- `os`
- `random`
- `datetime`
- `concurrent.futures`

### Quick installation

Install all required dependencies in one command:

```bash
pip install pysam tqdm
```

Clone and install dependencies:

```bash
git clone https://github.com/lindoancientdna/multicore-variant-caller.git
```

Example bash script:

    python main.py \
    --bams Sample1 Sample2 \
    --samples sample1.bam sample2.bam \
    --bed 1240k.gz \
    --out output.vcf \
    --cores 48 \
    --verbose

