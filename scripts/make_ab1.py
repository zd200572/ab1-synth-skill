#!/usr/bin/env python3
"""Build a synthetic Sanger .ab1 chromatogram from a reference sequence and reads.

The output is a valid ABIF container (34-byte header, sorted directory, 4-byte
aligned payloads) that opens in Chromas, SnapGene, FinchTV and Biopython.

Example:
    python make_ab1.py --fasta asv_rep.fasta --fastq trim.fq.gz \
        --id ASV1 --out ASV1.ab1
"""

from __future__ import annotations

import argparse
import gzip
import math
import os
import random
import struct
import sys
from collections import Counter, defaultdict

try:
    from Bio import SeqIO
except ImportError:  # pragma: no cover
    sys.exit("Biopython is required: pip install biopython")

BASES = "GATC"  # ABIF filter-wheel order used by FWO_1
PULLUP = {"A": "G", "G": "A", "C": "T", "T": "C"}  # neighbouring dye channels
COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}
BASE_SPACING = 12.0
SPACING_BIAS = {"A": 0.0, "C": -1.1, "G": 1.3, "T": 0.5}
HEIGHT_BIAS = {"A": 1.02, "C": 0.82, "G": 1.18, "T": 0.94}


def revcomp(seq: str) -> str:
    return "".join(COMPLEMENT[b] for b in reversed(seq))


def open_text(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def load_target(fasta: str, seq_id: str | None) -> tuple[str, str]:
    for record in SeqIO.parse(fasta, "fasta"):
        if seq_id is None or record.id == seq_id:
            return record.id, str(record.seq).upper()
    if seq_id is None:
        raise SystemExit(f"No sequences found in {fasta}")
    raise SystemExit(f"Sequence id {seq_id!r} not found in {fasta}")


def map_reads(fastq: str, target: str, k: int, max_mismatch_rate: float):
    """Return per-position base counts, quality sums, and the number of reads used."""
    if len(target) < k:
        k = max(9, len(target) // 2)
    seeds: dict[str, int] = {}
    for i in range(len(target) - k + 1):
        seeds.setdefault(target[i : i + k], i)
    unique = {s: p for s, p in seeds.items() if target.count(s) == 1}

    counts = [Counter() for _ in target]
    qsum = [0] * len(target)
    used = 0

    fmt = "fastq" if fastq.lower().endswith((".fq", ".fastq", ".fq.gz", ".fastq.gz")) else "fasta"
    with open_text(fastq) as handle:
        for record in SeqIO.parse(handle, fmt):
            seq = str(record.seq).upper()
            if len(seq) < k:
                continue
            qual = record.letter_annotations.get("phred_quality")
            votes: Counter = Counter()
            for j in range(len(seq) - k + 1):
                p = unique.get(seq[j : j + k])
                if p is not None:
                    votes[p - j] += 1
            if not votes:
                continue
            offset = votes.most_common(1)[0][0]
            lo, hi = max(0, offset), min(len(target), offset + len(seq))
            overlap = hi - lo
            if overlap < max(50, len(target) // 2):
                continue
            mismatches = sum(seq[p - offset] != target[p] for p in range(lo, hi))
            if mismatches > max(4, int(overlap * max_mismatch_rate)):
                continue
            used += 1
            for p in range(lo, hi):
                rp = p - offset
                base = seq[rp]
                if base in "ACGT":
                    counts[p][base] += 1
                    qsum[p] += max(0, (qual[rp] if qual else 40))
    return counts, qsum, used


def gaussian(channel, center, height, sigma_left, sigma_right, length):
    lo = max(0, int(center - 4 * sigma_left))
    hi = min(length, int(center + 4 * sigma_right) + 1)
    for x in range(lo, hi):
        sigma = sigma_left if x < center else sigma_right
        channel[x] += height * math.exp(-0.5 * ((x - center) / sigma) ** 2)


def build_trace(target, counts, qsum, seed=20261109, taper=False):
    rng = random.Random(seed)
    centers, heights, sigmas, quals = [], [], [], []
    pos = 26.0
    for i, called in enumerate(target):
        depth = sum(counts[i].values())
        q = round(qsum[i] / depth) if depth >= 5 else 18
        q = max(12, min(55, q))
        if taper:
            # a real read loses quality near both ends
            edge = min(i, len(target) - 1 - i, 60)
            q = min(q, int(11 + edge * 0.62))
        pos += BASE_SPACING + SPACING_BIAS[called] + rng.gauss(0, 1.0)
        height = 620 * HEIGHT_BIAS[called] * (0.62 + q / 55 * 0.55)
        height *= math.exp(rng.gauss(0, 0.16))
        centers.append(pos)
        heights.append(height)
        sigmas.append(rng.uniform(2.9, 3.6))
        quals.append(q)

    length = int(centers[-1] + 60)
    channels = {b: [0.0] * length for b in BASES}

    for i, called in enumerate(target):
        depth = sum(counts[i].values())
        c, h, s = centers[i], heights[i], sigmas[i]
        if depth >= 5:
            for b in BASES:
                fraction = counts[i][b] / depth
                if fraction > 0:
                    gaussian(channels[b], c + rng.gauss(0, 0.25),
                             max(8.0, h * fraction), s * 0.92, s * 1.3, length)
        else:
            gaussian(channels[called], c, h, s * 0.92, s * 1.3, length)
        # spectral pull-up between neighbouring dye channels
        gaussian(channels[PULLUP[called]], c + 0.4,
                 h * rng.uniform(0.02, 0.05), s, s, length)

    for b in BASES:
        phase = rng.uniform(0, 2 * math.pi)
        for x in range(length):
            channels[b][x] += 8 + 3 * math.sin(x / 90 + phase) + rng.gauss(0, 2.2)
    for _ in range(6):
        b = rng.choice(BASES)
        c = rng.uniform(40, max(41, length - 40))
        gaussian(channels[b], c, rng.uniform(15, 45), rng.uniform(6, 12), rng.uniform(6, 12), length)

    traces = {b: [max(0, min(32000, int(round(v)))) for v in channels[b]] for b in BASES}
    ploc = [int(round(c)) for c in centers]
    return traces, ploc, quals, length


def write_abif(path, target, traces, ploc, quals, sample_name):
    entries = []

    def add(name, number, elem_type, elem_size, data):
        entries.append([name.encode("ascii"), number, elem_type, elem_size, data])

    for number, base in zip((1, 2, 3, 4), BASES):
        add("DATA", number, 4, 2, struct.pack(">" + "h" * len(traces[base]), *traces[base]))
    for number, base in zip((9, 10, 11, 12), BASES):
        add("DATA", number, 4, 2, struct.pack(">" + "h" * len(traces[base]), *traces[base]))
    add("FWO_", 1, 2, 1, BASES.encode())
    add("MODL", 1, 2, 1, b"3730")
    add("PBAS", 1, 2, 1, target.encode("ascii"))
    add("PBAS", 2, 2, 1, target.encode("ascii"))
    add("PCON", 1, 2, 1, bytes(quals))
    add("PCON", 2, 2, 1, bytes(quals))
    add("PLOC", 1, 4, 2, struct.pack(">" + "h" * len(ploc), *ploc))
    add("PLOC", 2, 4, 2, struct.pack(">" + "h" * len(ploc), *ploc))
    add("SCAN", 1, 5, 4, struct.pack(">i", len(traces["G"])))
    encoded = sample_name.encode("ascii")[:250]
    add("SMPL", 1, 18, 1, bytes([len(encoded)]) + encoded)
    entries.sort(key=lambda e: (e[0], e[1]))

    header_size = 34  # 4-byte ABIF magic + 30-byte header (includes handle)
    payload = bytearray(b"\0" * (128 - header_size))
    directory = []
    for name, number, elem_type, elem_size, data in entries:
        if len(data) <= 4:
            offset = int.from_bytes(data.ljust(4, b"\0"), "big")
        else:
            payload.extend(b"\0" * ((-(header_size + len(payload))) % 4))
            offset = header_size + len(payload)
            payload.extend(data)
        directory.append((name, number, elem_type, elem_size,
                          len(data) // elem_size, len(data), offset, 0))
    payload.extend(b"\0" * ((-(header_size + len(payload))) % 4))
    dir_offset = header_size + len(payload)
    dir_bytes = b"".join(struct.pack(">4sI2H4I", *d) for d in directory)
    root = struct.pack(">H4sI2H3I4x", 101, b"tdir", 1, 1023, 28,
                       len(directory), len(dir_bytes), dir_offset)
    with open(path, "wb") as handle:
        handle.write(b"ABIF")
        handle.write(root)
        handle.write(payload)
        handle.write(dir_bytes)


def validate(path, target):
    record = SeqIO.read(path, "abi")
    if str(record.seq) != target:
        raise SystemExit("Validation failed: base calls differ from the reference")
    raw = record.annotations["abif_raw"]
    length = len(raw["DATA1"])
    for key in ("DATA2", "DATA3", "DATA4", "DATA9", "DATA10", "DATA11", "DATA12"):
        if len(raw[key]) != length:
            raise SystemExit(f"Validation failed: {key} length mismatch")
    return record


def open_in_chromas(path):
    if not sys.platform.startswith("win") and not os.path.exists("/mnt/c/Windows/System32/WindowsPowerShell/v1.0/powershell.exe"):
        print("Chromas launch skipped (Windows not detected)")
        return
    win_path = path.replace("/mnt/", "").replace("/", "\\")
    if win_path[1:2] != ":":
        drive = win_path[0].upper()
        win_path = f"{drive}:{win_path[1:]}"
    exe = r"C:\Program Files\Chromas\Chromas.exe"
    os.system(f'powershell.exe -NoProfile -Command "Start-Process -FilePath \'{exe}\' -ArgumentList \'{win_path}\'"')


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--fasta", required=True, help="reference FASTA (ASV/consensus)")
    parser.add_argument("--fastq", help="reads (FASTQ/FASTA, .gz ok) used for peak variation")
    parser.add_argument("--id", dest="seq_id", help="sequence id in the FASTA (default: first record)")
    parser.add_argument("--out", required=True, help="output .ab1 path")
    parser.add_argument("--sample", help="SMPL1 sample name (default: sequence id)")
    parser.add_argument("--seed", type=int, default=20261109, help="random seed")
    parser.add_argument("--kmer", type=int, default=17, help="seed k-mer length for read mapping")
    parser.add_argument("--max-mismatch-rate", type=float, default=0.04,
                        help="maximum mismatch fraction for a read to be used")
    parser.add_argument("--open", action="store_true", help="open the result in Chromas on Windows")
    parser.add_argument("--chunk-size", type=int, default=0, metavar="N",
                        help="split the reference into reads of N bp (e.g. 800); one .ab1 per chunk")
    parser.add_argument("--overlap", type=int, default=0,
                        help="overlap between consecutive chunks when splitting")
    parser.add_argument("--both-strands", action="store_true",
                        help="write <out>_F.ab1 and <out>_R.ab1 (reverse complement)")
    args = parser.parse_args(argv)

    seq_id, target = load_target(args.fasta, args.seq_id)
    if not target:
        raise SystemExit("Reference sequence is empty")

    if args.fastq:
        counts, qsum, used = map_reads(args.fastq, target, args.kmer, args.max_mismatch_rate)
    else:
        counts, qsum, used = [Counter() for _ in target], [0] * len(target), 0
    print(f"Reference: {seq_id} ({len(target)} bp); reads used: {used}")

    if args.chunk_size > 0:
        step = args.chunk_size - args.overlap
        if step <= 0:
            raise SystemExit("--overlap must be smaller than --chunk-size")
        spans = []
        start = 0
        while start < len(target):
            end = min(start + args.chunk_size, len(target))
            if end == len(target) and end - start < args.chunk_size:
                # keep every read at full length: shift the last window back
                start = max(0, end - args.chunk_size)
            spans.append((start, end))
            if end == len(target):
                break
            start += step
    else:
        spans = [(0, len(target))]

    base, ext = os.path.splitext(args.out)
    both_strands = args.both_strands or (args.chunk_size > 0 and len(target) < args.chunk_size)
    if both_strands and len(spans) > 1:
        raise SystemExit("--both-strands cannot be combined with multi-chunk output")

    jobs = []
    if both_strands:
        if len(target) < args.chunk_size:
            print(f"Reference ({len(target)} bp) is shorter than --chunk-size; "
                  "writing forward and reverse-complement reads instead")
        rev_counts = [
            Counter({COMPLEMENT[b]: n for b, n in counts[len(target) - 1 - i].items()})
            for i in range(len(target))
        ]
        rev_qsum = list(reversed(qsum))
        jobs.append((target, counts, qsum, f"{base}_F{ext or '.ab1'}",
                     f"{args.sample or seq_id}_F", False))
        jobs.append((revcomp(target), rev_counts, rev_qsum, f"{base}_R{ext or '.ab1'}",
                     f"{args.sample or seq_id}_R", False))
    elif len(spans) == 1:
        jobs.append((target, counts, qsum, args.out, args.sample or seq_id, False))
    else:
        for index, (start, end) in enumerate(spans):
            jobs.append((target[start:end], counts[start:end], qsum[start:end],
                         f"{base}_part{index + 1}{ext or '.ab1'}",
                         f"{args.sample or seq_id}_{index + 1}", True))

    for index, (seq, seq_counts, seq_qsum, out_path, sample, taper) in enumerate(jobs):
        traces, ploc, quals, length = build_trace(
            seq, seq_counts, seq_qsum, seed=args.seed + index, taper=taper)
        write_abif(out_path, seq, traces, ploc, quals, sample)
        validate(out_path, seq)
        print(f"Wrote {out_path} ({os.path.getsize(out_path)} bytes, "
              f"{len(seq)} bp, {length} trace samples)")
    print("ABIF validation passed")
    if args.open:
        for job in jobs:
            open_in_chromas(job[3])


if __name__ == "__main__":
    main()
