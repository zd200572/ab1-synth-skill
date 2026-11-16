# ab1-synth-skill

Turn a reference sequence (ASV / OTU / consensus FASTA) plus its supporting
reads into a **synthetic Sanger `.ab1` chromatogram** that opens in Chromas,
FinchTV, SnapGene and Biopython.

FASTQ files contain no fluorescence channels, so the trace is reconstructed:
read variation supplies minor peaks, Phred scores drive peak heights, and peak
spacing/width/baseline are modelled to look like a real electropherogram.
**The result is simulated data, not instrument output.**

## Requirements

- Python 3.9+
- Biopython (`pip install biopython`)

## Usage

```bash
python scripts/make_ab1.py \
  --fasta asv_rep.fasta \
  --id ASV1 \
  --fastq trim.L9.fq.gz \
  --out ASV1.ab1
```

| Option | Meaning |
| --- | --- |
| `--fasta` | reference FASTA |
| `--id` | sequence id to use (default: first record) |
| `--fastq` | reads for peak variation (FASTQ/FASTA, `.gz` supported); optional |
| `--out` | output `.ab1` |
| `--sample` | `SMPL1` sample name |
| `--seed` | random seed for reproducible traces |
| `--chunk-size` | split the reference into reads of N bp, one `.ab1` per chunk |
| `--overlap` | overlap between consecutive chunks when splitting |
| `--both-strands` | write `<out>_F.ab1` and `<out>_R.ab1` (reverse complement) |
| `--open` | open the result in Chromas on Windows |

### Split a long sequence into several reads

```bash
python scripts/make_ab1.py --fasta ref.fasta --fastq reads.fq.gz \
  --chunk-size 800 --out sample.ab1
# -> sample_part1.ab1 (800 bp), sample_part2.ab1 (800 bp, overlaps the tail)
```

Chunks are consecutive and non-overlapping; `--overlap 50` adds an overlap.
Every part is exactly `--chunk-size` long: a short tail is covered by shifting
the last window back, so it overlaps the previous read. Each part is treated
as an independent read: quality is tapered at both ends and the random seed
differs per part.

### Short sequences: forward + reverse

When the reference is shorter than `--chunk-size` (or with `--both-strands`),
two full-length reads are written instead: `sample_F.ab1` and `sample_R.ab1`
(reverse complement). Read evidence is mirrored, so the reverse trace carries
the same variation and quality profile.

```bash
python scripts/make_ab1.py --fasta asv_rep.fasta --id ASV1 \
  --fastq trim.L9.fq.gz --chunk-size 800 --out ASV1.ab1
# -> ASV1_F.ab1, ASV1_R.ab1
```


Output is validated by re-reading it with Biopython: base calls must match the
reference and all eight trace channels must have equal length.

## How it works

1. Unique 17-mers from the reference vote for a read offset; reads must cover at
   least half the reference and stay below a 4% mismatch rate.
2. Per-position base counts and quality sums are accumulated.
3. Peaks are drawn with base-dependent spacing and height, asymmetric Gaussian
   shape, spectral pull-up, baseline drift, noise and a few dye bumps.
4. A valid ABIF container is written.

## ABIF gotchas

- Header is **34 bytes** (it includes a 4-byte handle). A 30-byte header makes
  Chromas report *"not a recognised file type or the data is corrupted"*.
- Payload starts at 128; every data offset and the directory offset must be
  4-byte aligned.
- Directory entries are sorted by tag name, then number.
- `FWO_1 = GATC` defines the channel order of `DATA1..4`.

