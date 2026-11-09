#!/usr/bin/env python3
"""
dupcaller_realign.py - local re-assembly + realignment pre-step for DupCaller.

Sits between "align (BWA-MEM)" and "DupCaller call".  It looks for regions where
reads disagree with the reference (mismatch clusters, indels, soft-clips), builds
a small de Bruijn graph there (the same idea as GATK HaplotypeCaller / Mutect2
local re-assembly), extracts the well-supported haplotypes, and re-aligns reads
that fit one of those haplotypes much better than they fit their current
alignment.  Everything else is left byte-for-byte alone.

Design rules that matter for a duplex-somatic caller:
  * Support is counted in MOLECULES (UMI tag, or read-name if no tag), not reads,
    so PCR copies of one molecule cannot create a haplotype by themselves.
  * Haplotypes are discovered jointly across all BAMs given (tumour + matched
    normal), so both samples are realigned against the same candidate set.
  * A read is only rewritten if the new alignment beats the old one by
    --min-gain score points. Reads that already fit well are never touched.
  * Output uses the same reference/coordinates; changed reads get tags
    OA (original alignment) and ZR:i:1 so you can audit every change.

Usage:
  python dupcaller_realign.py --bam tumour.bam normal.bam --ref ref.fa \
      --out-dir realigned/ --umi-tag RX --report regions.tsv

Requires: python>=3.9, numpy, pysam.
"""
import argparse
import heapq
import os
import sys
from collections import defaultdict

import numpy as np
import pysam

# ---------------------------------------------------------------- scoring ----
# BWA-MEM-like defaults: match +1, mismatch -4, gap of length n costs 6 + n, clip 5
MATCH, MISMATCH, GAP_OPEN, GAP_EXT, CLIP_PEN = 1, -4, 6, 1, 5
NEG = -10 ** 6


def log(*a):
    print(*a, file=sys.stderr, flush=True)


# ------------------------------------------------------ affine-gap alignment --
def affine_align(q, t, local=True):
    """
    Affine-gap alignment of query q against target t (numpy, row-vectorised).
      local=True : Smith-Waterman, free clipping on both sequences (used read->target)
      local=False: global (Needleman-Wunsch)               (used haplotype->window)
    Returns dict(score, ops, q_start, q_end, t_start, t_end) where ops is a list of
    (op, length), op in 'M','I','D' (I = extra base in q, D = extra base in t),
    covering only the aligned part (clips are implied by q_start/q_end).
    Gaps are pushed as far LEFT as possible (diagonal preferred in traceback).
    """
    qa = np.frombuffer(q.encode(), dtype=np.uint8)
    ta = np.frombuffer(t.encode(), dtype=np.uint8)
    n, m = len(qa), len(ta)
    N = ord("N")
    H = np.zeros((n + 1, m + 1), dtype=np.int32)
    E = np.full((n + 1, m + 1), NEG, dtype=np.int32)  # consumes target (deletion)
    F = np.full((n + 1, m + 1), NEG, dtype=np.int32)  # consumes query  (insertion)
    jj = np.arange(m + 1, dtype=np.int32)
    if not local:
        H[0, 1:] = -(GAP_OPEN + GAP_EXT * jj[1:])
        E[0, 1:] = H[0, 1:]
        H[1:, 0] = -(GAP_OPEN + GAP_EXT * np.arange(1, n + 1))
        F[1:, 0] = H[1:, 0]
    tN = ta == N
    for i in range(1, n + 1):
        sub = np.where(ta == qa[i - 1], MATCH, MISMATCH).astype(np.int32)
        if qa[i - 1] == N:
            sub[:] = 0
        else:
            sub[tN] = 0
        F[i, 1:] = np.maximum(F[i - 1, 1:] - GAP_EXT, H[i - 1, 1:] - GAP_OPEN - GAP_EXT)
        pre = np.empty(m + 1, dtype=np.int32)
        pre[0] = H[i, 0]
        pre[1:] = np.maximum(H[i - 1, :-1] + sub, F[i, 1:])
        if local:
            pre = np.maximum(pre, 0)
        # E[j] = max_{j'<j} pre[j'] - open - ext*(j-j')   (vectorised with cummax)
        g = pre + GAP_EXT * jj
        cm = np.maximum.accumulate(g)
        E[i, 1:] = cm[:-1] - GAP_OPEN - GAP_EXT * jj[1:]
        H[i] = np.maximum(pre, E[i])
        if not local:
            H[i, 0] = -(GAP_OPEN + GAP_EXT * i)

    if local:
        last = int(np.argmax(H[n]))
        s_full = int(H[n, last])
        flat = int(np.argmax(H))
        bi, bj = divmod(flat, m + 1)
        s_max = int(H[bi, bj])
        if s_full > 0 and s_full >= s_max - CLIP_PEN:
            i, j = n, last
        else:
            i, j = bi, bj
        if int(H[i, j]) <= 0:
            return None
    else:
        i, j = n, m
    q_end, t_end = i, j
    score = int(H[i, j])

    ops = []
    state = "H"
    while True:
        if state == "H":
            if local and H[i, j] == 0:
                break
            if i == 0 and j == 0:
                break
            if i > 0 and j > 0:
                a, b = qa[i - 1], ta[j - 1]
                s = 0 if (a == N or b == N) else (MATCH if a == b else MISMATCH)
                if H[i, j] == H[i - 1, j - 1] + s:
                    ops.append("M"); i -= 1; j -= 1
                    continue
            if j > 0 and H[i, j] == E[i, j]:
                state = "E"
            elif i > 0 and H[i, j] == F[i, j]:
                state = "F"
            elif j == 0:
                state = "F"
            else:
                state = "E"
        elif state == "E":
            ops.append("D")
            if j >= 2 and E[i, j] == E[i, j - 1] - GAP_EXT:
                pass
            else:
                state = "H"
            j -= 1
        else:
            ops.append("I")
            if i >= 2 and F[i, j] == F[i - 1, j] - GAP_EXT:
                pass
            else:
                state = "H"
            i -= 1
    ops.reverse()
    rle = []
    for o in ops:
        if rle and rle[-1][0] == o:
            rle[-1][1] += 1
        else:
            rle.append([o, 1])
    return dict(score=score, ops=[tuple(x) for x in rle], q_start=i, q_end=q_end,
                t_start=j, t_end=t_end)


# ----------------------------------------------- alignment bookkeeping helpers --
def cols_from_read(read, ref_off):
    """pysam read -> per-column list [(op,q,r)], r relative to window start (ref_off).
    Returns None for CIGARs we do not handle (N, H, P)."""
    cols, q, r = [], 0, read.reference_start - ref_off
    for op, ln in read.cigartuples:
        if op in (0, 7, 8):
            for _ in range(ln):
                cols.append(("M", q, r)); q += 1; r += 1
        elif op == 1:
            for _ in range(ln):
                cols.append(("I", q, None)); q += 1
        elif op == 2:
            for _ in range(ln):
                cols.append(("D", None, r)); r += 1
        elif op == 4:
            q += ln
        else:
            return None
    return cols


def score_cols(cols, seq, win, clip_l, clip_r):
    s, prev = 0, None
    for op, q, r in cols:
        if op == "M":
            a, b = seq[q], win[r]
            s += 0 if (a == "N" or b == "N") else (MATCH if a == b else MISMATCH)
        else:
            s -= GAP_EXT + (GAP_OPEN if prev != op else 0)
        prev = op
    return s, s - CLIP_PEN * ((clip_l > 0) + (clip_r > 0))


def compose_to_ref(q_start, t_start, ops, h2r):
    """Read->target alignment ops  +  target->window map  =>  read->window columns."""
    emis, qi, tj = [], q_start, t_start
    for op, ln in ops:
        for _ in range(ln):
            if op == "M":
                r = tj if h2r is None else h2r[tj]
                emis.append((qi, r if r >= 0 else None)); qi += 1; tj += 1
            elif op == "I":
                emis.append((qi, None)); qi += 1
            else:
                tj += 1
    cols, prev = [], None
    for qi, r in emis:
        if r is None:
            cols.append(["I", qi, None])
        else:
            if prev is not None and r > prev + 1:
                for rr in range(prev + 1, r):
                    cols.append(["D", None, rr])
            cols.append(["M", qi, r]); prev = r
    while cols and cols[0][0] != "M":
        cols.pop(0)
    while cols and cols[-1][0] != "M":
        cols.pop()
    return cols


def left_align(cols, seq, win):
    """Shift every indel to its leftmost equivalent position (GATK/BWA convention)."""
    s = 1
    while s < len(cols):
        if cols[s][0] in "ID":
            kind = cols[s][0]
            e = s
            while e < len(cols) and cols[e][0] == kind:
                e += 1
            while s >= 2 and cols[s - 1][0] == "M":
                m = cols[s - 1]
                qp, rp = m[1], m[2]
                if seq[qp] != win[rp]:
                    break
                if kind == "D":
                    if win[cols[s][2] - 1] != win[cols[e - 1][2]]:
                        break
                    blk = [["D", None, c[2] - 1] for c in cols[s:e]]
                    newm = ["M", qp, cols[e - 1][2]]
                else:
                    if seq[qp] != seq[cols[e - 1][1]]:
                        break
                    blk = [["I", c[1] - 1, None] for c in cols[s:e]]
                    newm = ["M", cols[e - 1][1], rp]
                cols[s - 1:e] = blk + [newm]
                s -= 1; e -= 1
            s = e
        else:
            s += 1
    return cols


def finish_alignment(cols, seq, win, win_start):
    """cols -> (ref_start_abs, cigartuples, NM, MD, AS).  Adds soft clips."""
    q_first, q_last = cols[0][1], cols[-1][1]
    n = len(seq)
    ops = []
    if q_first:
        ops.append([4, q_first])
    code = {"M": 0, "I": 1, "D": 2}
    for op, _, _ in cols:
        c = code[op]
        if ops and ops[-1][0] == c:
            ops[-1][1] += 1
        else:
            ops.append([c, 1])
    if n - 1 - q_last:
        ops.append([4, n - 1 - q_last])
    nm, md, run, i = 0, [], 0, 0
    while i < len(cols):
        op, q, r = cols[i]
        if op == "M":
            if seq[q] == win[r] or seq[q] == "N":
                run += 1
            else:
                md.append(str(run)); md.append(win[r]); run = 0; nm += 1
        elif op == "I":
            nm += 1
        else:
            md.append(str(run)); run = 0
            dele = ""
            while i < len(cols) and cols[i][0] == "D":
                dele += win[cols[i][2]]; nm += 1; i += 1
            md.append("^" + dele)
            continue
        i += 1
    md.append(str(run))
    sc, _ = score_cols(cols, seq, win, 0, 0)
    return win_start + cols[0][2], [tuple(x) for x in ops], nm, "".join(md), sc


# ------------------------------------------------------- active region scan --
def read_ok(r, min_mapq):
    return not (r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_qcfail
                or r.mapping_quality < min_mapq or r.cigartuples is None
                or r.query_sequence is None)


def find_active_regions(bams, fasta, contig, c_start, c_end, a):
    """Positions where >= min_alt_reads reads (and >= min_alt_frac of coverage) show a
    mismatch / indel / soft-clip, merged into cores, then padded into windows."""
    flagged = []
    for cs in range(c_start, c_end, a.chunk):
        ce = min(cs + a.chunk, c_end)
        ev = np.zeros(ce - cs, dtype=np.int32)
        cov = np.zeros(ce - cs + 1, dtype=np.int32)
        ref = fasta.fetch(contig, max(0, cs - 2000), ce + 2000).upper()
        ref_off = max(0, cs - 2000)
        for bam in bams:
            for r in bam.fetch(contig, cs, ce):
                if not read_ok(r, a.min_mapq):
                    continue
                rs, re_ = r.reference_start, r.reference_end
                cov[max(rs, cs) - cs] += 1
                cov[min(re_, ce) - cs] -= 1
                seq, qual = r.query_sequence, r.query_qualities
                pos, q, hits = rs, 0, set()
                for op, ln in r.cigartuples:
                    if op in (0, 7, 8):
                        for k in range(ln):
                            rb = ref[pos + k - ref_off]
                            if seq[q + k] != rb and rb != "N" and seq[q + k] != "N" \
                                    and qual[q + k] >= a.min_baseq:
                                hits.add(pos + k)
                        pos += ln; q += ln
                    elif op == 1:
                        hits.add(pos); q += ln
                    elif op == 2:
                        hits.update(range(pos, pos + ln)); pos += ln
                    elif op == 4:
                        if ln >= a.min_clip:
                            hits.add(pos)
                        q += ln
                    elif op == 3:
                        pos += ln
                for h in hits:
                    if cs <= h < ce:
                        ev[h - cs] += 1
        cov = np.cumsum(cov)[:-1]
        act = (ev >= a.min_alt_reads) & (ev >= a.min_alt_frac * np.maximum(cov, 1))
        flagged.extend((np.nonzero(act)[0] + cs).tolist())
    if not flagged:
        return []
    cores, s, p = [], flagged[0], flagged[0]
    for x in flagged[1:]:
        if x - p > a.merge_dist:
            cores.append((s, p + 1)); s = x
        p = x
    cores.append((s, p + 1))
    out = []
    for cs, ce in cores:
        if ce - cs > a.max_region:
            log(f"  skip core {contig}:{cs}-{ce}: longer than --max-region (likely SV/repeat)")
            continue
        ws, we = max(c_start, cs - a.pad), min(c_end, ce + a.pad)
        if out and ws <= out[-1][2]:       # overlapping windows -> merge
            out[-1] = (contig, out[-1][1], we, out[-1][3], ce)
        else:
            out.append((contig, ws, we, cs, ce))
    return out


# ------------------------------------------------------------- re-assembly ---
def prune_graph(support, ref_edges, min_pruning):
    """GATK-style chain pruning: drop a linear chain only if EVERY edge is weakly supported."""
    outdeg, indeg = defaultdict(int), defaultdict(int)
    for (a_, b_) in support:
        outdeg[a_] += 1; indeg[b_] += 1
    parent = {e: e for e in support}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    by_src, by_dst = defaultdict(list), defaultdict(list)
    for e in support:
        by_src[e[0]].append(e); by_dst[e[1]].append(e)
    for node in by_dst:
        if indeg[node] == 1 and outdeg[node] == 1:
            parent[find(by_dst[node][0])] = find(by_src[node][0])
    comp = defaultdict(list)
    for e in support:
        comp[find(e)].append(e)
    keep = set()
    for edges in comp.values():
        if any(e in ref_edges for e in edges) or \
                max(len(support[e]) for e in edges) >= min_pruning:
            keep.update(edges)
    return keep


def best_paths(weights, source, sink, max_haps, max_expansions=60000):
    out = defaultdict(list)
    rev = defaultdict(list)
    for (a_, b_), w in weights.items():
        out[a_].append((b_, w)); rev[b_].append(a_)
    reach, stack = {sink}, [sink]
    while stack:
        for p in rev[stack.pop()]:
            if p not in reach:
                reach.add(p); stack.append(p)
    if source not in reach:
        return []
    tot = {a_: sum(w for _, w in lst) for a_, lst in out.items()}
    heap, cnt, res, exp = [(0.0, 0, (source,))], 1, [], 0
    while heap and len(res) < max_haps and exp < max_expansions:
        cost, _, path = heapq.heappop(heap)
        node = path[-1]
        if node == sink:
            res.append(path); continue
        exp += 1
        seen = set(path)
        for dst, w in out.get(node, ()):
            if dst in seen or dst not in reach:
                continue
            heapq.heappush(heap, (cost - np.log(w / tot[node]), cnt, path + (dst,))); cnt += 1
    return res


def assemble(window, rinfo, a):
    """Return (haplotype list excluding the reference window, list of k actually used)."""
    ks = sorted(a.kmers)

    def ref_unique(k):
        km = [window[i:i + k] for i in range(len(window) - k + 1)]
        return len(km) > 1 and len(set(km)) == len(km)
    good = [k for k in ks if ref_unique(k)]
    kmax, tries = ks[-1], 0
    while not good and tries < 6:                        # same fallback as GATK
        kmax += 10; tries += 1
        if kmax < len(window) - 5 and ref_unique(kmax):
            good = [kmax]
    if not good:
        return [], []
    haps = []
    for k in good:
        ref_k = [window[i:i + k] for i in range(len(window) - k + 1)]
        ref_edges = set(zip(ref_k, ref_k[1:]))
        support = defaultdict(set)
        for seq, qual, mol in rinfo:
            prev = None
            for i in range(len(seq) - k + 1):
                km = seq[i:i + k]
                if "N" in km or qual[i:i + k].min() < a.min_baseq:
                    prev = None; continue
                if prev is not None:
                    support[(prev, km)].add(mol)
                prev = km
        for e in ref_edges:
            support.setdefault(e, set())
        keep = prune_graph(support, ref_edges, a.min_pruning)
        weights = {e: len(support[e]) + 0.01 for e in keep}
        for p in best_paths(weights, ref_k[0], ref_k[-1], a.max_haps):
            h = p[0] + "".join(n[-1] for n in p[1:])
            if h != window and h not in haps:
                haps.append(h)
    return haps, good


# ------------------------------------------------------------ region worker --
def molecule_id(r, umi_tag):
    """Key that identifies the original DNA molecule. ADAPT THIS to how your pipeline
    stores barcodes (DupCaller's trimming step may put them in the read name)."""
    if umi_tag and r.has_tag(umi_tag):
        return (r.get_tag(umi_tag), min(r.reference_start,
                r.next_reference_start if r.is_paired and r.next_reference_id == r.reference_id
                else r.reference_start) // 50)
    return r.query_name


def process_region(bams, fasta, region, a, changes):
    contig, ws, we, cs, ce = region
    window = fasta.fetch(contig, ws, we).upper()
    reads = []
    for bi, bam in enumerate(bams):
        for r in bam.fetch(contig, ws, we):
            if read_ok(r, a.min_mapq):
                reads.append((bi, r))
    rinfo = [(r.query_sequence, np.asarray(r.query_qualities, dtype=np.int16),
              molecule_id(r, a.umi_tag)) for _, r in reads]
    haps, ks = assemble(window, rinfo, a)
    rec = dict(contig=contig, start=ws, end=we, core_start=cs, core_end=ce,
               n_reads=len(reads), k=",".join(map(str, ks)), n_haps=len(haps), n_realigned=0)
    # haplotype -> window maps (global alignment, like GATK's SW step 6)
    hmaps = []
    for h in haps:
        al = affine_align(h, window, local=False)
        h2r, qi, tj = [], 0, 0
        for op, ln in al["ops"]:
            for _ in range(ln):
                if op == "M":
                    h2r.append(tj); tj += 1
                elif op == "I":
                    h2r.append(-1)
                else:
                    tj += 1
        hmaps.append(h2r)

    for bi, r in reads:
        key = (bi, r.query_name, r.flag, r.reference_id, r.reference_start, r.cigarstring)
        if key in changes:
            continue
        cig = r.cigartuples
        lc = cig[0][1] if cig[0][0] == 4 else 0
        rc = cig[-1][1] if cig[-1][0] == 4 else 0
        if r.reference_end + rc <= cs or r.reference_start - lc >= ce:
            continue                                       # not near the active core
        if r.reference_start < ws + 5 or r.reference_end > we - 5:
            continue
        cols0 = cols_from_read(r, ws)
        if cols0 is None:
            continue
        seq = r.query_sequence
        _, orig = score_cols(cols0, seq, window, lc, rc)
        if len(seq) - orig < a.min_gain:                   # already fits well
            continue
        best_ref = affine_align(seq, window, local=True)
        s_ref = -10 ** 9
        if best_ref:
            s_ref = best_ref["score"] - CLIP_PEN * ((best_ref["q_start"] > 0) + (best_ref["q_end"] < len(seq)))
        best_hap, s_hap, hap_i = None, -10 ** 9, -1
        for hi, h in enumerate(haps):
            al = affine_align(seq, h, local=True)
            if not al or al["t_start"] == 0 or al["t_end"] == len(h):
                continue
            s = al["score"] - CLIP_PEN * ((al["q_start"] > 0) + (al["q_end"] < len(seq)))
            if s > s_hap:
                best_hap, s_hap, hap_i = al, s, hi
        if best_hap is not None and s_hap - max(orig, s_ref) >= a.min_gain:
            cols = compose_to_ref(best_hap["q_start"], best_hap["t_start"], best_hap["ops"], hmaps[hap_i])
        elif best_ref and s_ref - orig >= a.min_gain:
            cols = compose_to_ref(best_ref["q_start"], best_ref["t_start"], best_ref["ops"], None)
        else:
            continue
        if len(cols) < 20:
            continue
        cols = left_align(cols, seq, window)
        new_start, cigt, nm, md, sc = finish_alignment(cols, seq, window, ws)
        if abs(new_start - r.reference_start) > a.max_shift:
            continue
        if new_start == r.reference_start and cigt == tuple(map(tuple, cig)):
            continue
        changes[key] = (new_start, cigt, nm, md, sc)
        rec["n_realigned"] += 1
    return rec


# ------------------------------------------------------------------- driver --
def parse_region(s, lengths):
    if ":" not in s:
        return s, 0, lengths[s]
    c, rng = s.split(":")
    lo, hi = rng.replace(",", "").split("-")
    return c, int(lo) - 1, min(int(hi), lengths[c])


def write_output(path_in, path_out, bi, changes, fix_mates, tmpdir):
    with pysam.AlignmentFile(path_in, "rb") as fin:
        tmp = os.path.join(tmpdir, os.path.basename(path_out) + ".tmp.bam")
        with pysam.AlignmentFile(tmp, "wb", template=fin) as fout:
            for r in fin.fetch(until_eof=True):
                key = (bi, r.query_name, r.flag, r.reference_id, r.reference_start, r.cigarstring)
                ch = changes.get(key)
                if ch:
                    new_start, cigt, nm, md, sc = ch
                    oa = f"{fin.get_reference_name(r.reference_id)},{r.reference_start + 1}," \
                         f"{'-' if r.is_reverse else '+'},{r.cigarstring},{r.mapping_quality},"
                    oa += f"{r.get_tag('NM') if r.has_tag('NM') else ''};"
                    r.reference_start = new_start
                    r.cigartuples = list(cigt)
                    r.set_tag("NM", nm, "i"); r.set_tag("MD", md, "Z"); r.set_tag("AS", sc, "i")
                    r.set_tag("OA", oa, "Z"); r.set_tag("ZR", 1, "i")
                fout.write(r)
    if fix_mates:
        ns, fx = tmp + ".nsort.bam", tmp + ".fix.bam"
        pysam.sort("-n", "-o", ns, tmp)
        pysam.fixmate("-O", "bam", ns, fx)
        pysam.sort("-o", path_out, fx)
        for f in (ns, fx):
            os.remove(f)
    else:
        pysam.sort("-o", path_out, tmp)
    os.remove(tmp)
    pysam.index(path_out)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bam", nargs="+", required=True, help="coordinate-sorted, indexed BAMs (tumour + normal)")
    p.add_argument("--ref", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--region", help="chr or chr:start-end (1-based); default whole genome")
    p.add_argument("--umi-tag", default=None, help="BAM tag holding the molecule barcode (e.g. RX)")
    p.add_argument("--report", default=None, help="TSV summary of regions processed")
    p.add_argument("--kmers", type=int, nargs="+", default=[15, 25, 35])
    p.add_argument("--min-pruning", type=int, default=2, help="min supporting MOLECULES per non-reference chain")
    p.add_argument("--max-haps", type=int, default=16, help="per k")
    p.add_argument("--min-gain", type=int, default=8, help="score gain needed to rewrite a read")
    p.add_argument("--min-alt-reads", type=int, default=3)
    p.add_argument("--min-alt-frac", type=float, default=0.05)
    p.add_argument("--min-clip", type=int, default=5)
    p.add_argument("--min-mapq", type=int, default=10)
    p.add_argument("--min-baseq", type=int, default=15)
    p.add_argument("--merge-dist", type=int, default=50)
    p.add_argument("--pad", type=int, default=200, help="flank added around active core (>= read length)")
    p.add_argument("--max-region", type=int, default=1500)
    p.add_argument("--max-shift", type=int, default=100)
    p.add_argument("--chunk", type=int, default=1_000_000)
    p.add_argument("--no-fix-mates", action="store_true", help="skip samtools fixmate after realignment")
    a = p.parse_args()

    os.makedirs(a.out_dir, exist_ok=True)
    fasta = pysam.FastaFile(a.ref)
    bams = [pysam.AlignmentFile(b, "rb") for b in a.bam]
    lengths = dict(zip(bams[0].references, bams[0].lengths))
    contigs = [parse_region(a.region, lengths)] if a.region else [(c, 0, l) for c, l in lengths.items()]

    changes, report = {}, []
    for contig, s, e in contigs:
        regions = find_active_regions(bams, fasta, contig, s, e, a)
        log(f"{contig}: {len(regions)} active region(s)")
        for reg in regions:
            rec = process_region(bams, fasta, reg, a, changes)
            report.append(rec)
            log(f"  {contig}:{rec['start']}-{rec['end']}  reads={rec['n_reads']} k={rec['k']} "
                f"haplotypes={rec['n_haps']} realigned={rec['n_realigned']}")
    for b in bams:
        b.close()
    for bi, path in enumerate(a.bam):
        out = os.path.join(a.out_dir, os.path.basename(path).replace(".bam", "") + ".realigned.bam")
        write_output(path, out, bi, changes, not a.no_fix_mates, a.out_dir)
        log("wrote", out)
    if a.report:
        cols = ["contig", "start", "end", "core_start", "core_end", "n_reads", "k", "n_haps", "n_realigned"]
        with open(a.report, "w") as fh:
            fh.write("\t".join(cols) + "\n")
            for rec in report:
                fh.write("\t".join(str(rec[c]) for c in cols) + "\n")
        log("report:", a.report)


if __name__ == "__main__":
    main()