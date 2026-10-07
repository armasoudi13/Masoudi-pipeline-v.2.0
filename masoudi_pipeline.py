#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
masoudi-pipeline : Oxford Nanopore (ONT) reference-based pipeline

INPUT   one or many FASTQ(.gz) files (or barcode folders) + a custom reference FASTA
STEP 1  FastQC -> Trim Galore -> FastQC (validation) -> MultiQC
STEP 2  minimap2 | BWA-MEM2 -> SAM -> filter/sort -> BAM + BAI -> samtools flagstat
STEP 3  samtools markdup (duplicate removal) -> LoFreq + bcftools(mpileup) SNP calling
STEP 4  bcftools consensus (low-depth positions masked with N)
OUTPUT  VCF, consensus FASTA, QC summary (TSV + MultiQC), QC / depth / VCF plots
"""
import argparse
import glob
import gzip
import os
import re
import shlex
import shutil
import subprocess
import sys
import time
import traceback
from array import array
from collections import OrderedDict

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

__version__ = "2.0.0"

# --------------------------------------------------------------------------- #
#  Pac-Man header
# --------------------------------------------------------------------------- #
PACMAN = r"""
            .-~~~~~~~-.
         .'   _        '.
        /    (o)    _.-~'
       |          .-'              .    .    .    .    .    .    .
       |          '-.
        \           '-._
         '.             .'
           '-.._____..-'
"""

_L = {
    "M": [" __  __ ", "|  \\/  |", "| |\\/| |", "| |  | |", "|_|  |_|"],
    "A": ["    _    ", "   / \\   ", "  / _ \\  ", " / ___ \\ ", "/_/   \\_\\"],
    "S": [" ____  ", "/ ___| ", "\\___ \\ ", " ___) |", "|____/ "],
    "O": ["  ___  ", " / _ \\ ", "| | | |", "| |_| |", " \\___/ "],
    "U": [" _   _ ", "| | | |", "| | | |", "| |_| |", " \\___/ "],
    "D": [" ____  ", "|  _ \\ ", "| | | |", "| |_| |", "|____/ "],
    "I": [" ___ ", "|_ _|", " | | ", " | | ", "|___|"],
}


_NOCOLOR = bool(os.environ.get("NO_COLOR")) or os.environ.get("TERM") == "dumb"
_ANSI = dict(BLUE="\033[94m", LBLUE="\033[96m", DBLUE="\033[34m", GREEN="\033[92m", YELLOW="\033[93m",
             RED="\033[91m", MAGENTA="\033[95m", WHITE="\033[97m", CYAN="\033[96m", BOLD="\033[1m", RESET="\033[0m")
_ASCII_MAP = {"═": "=", "║": "|", "╔": "+", "╗": "+", "╚": "+", "╝": "+", "█": "#", "─": "-", "│": "|",
              "┌": "+", "┐": "+", "└": "+", "┘": "+", "├": "+", "┤": "+", "↳": ">", "▸": ">", "●": "*",
              "◉": "o", "·": ".", "→": ">"}


def banner():
    c = {k: ("" if _NOCOLOR else v) for k, v in _ANSI.items()}
    BLUE, LBLUE, DBLUE, GREEN, YELLOW = c["BLUE"], c["LBLUE"], c["DBLUE"], c["GREEN"], c["YELLOW"]
    RED, MAGENTA, WHITE, CYAN, BOLD, RESET = c["RED"], c["MAGENTA"], c["WHITE"], c["CYAN"], c["BOLD"], c["RESET"]
    out = [""]

    # ---- Pac-Man logo ----
    for row in ("        ████████                          ",
                "    ████████████████                      ",
                "  ████████████████████                    ",
                " ████████████████████                     ",
                "██████████████████                        ",
                "████████████████    %s◉%s                   " % (WHITE, YELLOW),
                "████████                              ",
                "████████████████    ·  ·  ·  ·  ·  ·      ",
                "██████████████████                        ",
                " ████████████████████                     ",
                "  ████████████████████                    ",
                "    ████████████████                      ",
                "        ████████                          "):
        out.append("                %s%s%s" % (YELLOW, row, RESET))
    out.append("              %s%s         MASOUDI PIPELINE  🦠  Oxford Nanopore%s" % (BOLD, WHITE, RESET))
    out.append("")

    # ---- Welcome box ----
    def box(segs, width=51):
        vis = sum(len(t) + sum(1 for ch in t if ord(ch) > 0xFFFF) for t, _ in segs)  # emoji = 2 columns
        pad = max(0, width - vis)
        left = pad // 2 if segs[0][1] == "center" else 0
        body = "".join("%s%s%s" % ("" if col == "center" else col, t, RESET) for t, col in segs)
        return "%s%s     ║%s%s%s%s%s%s║%s" % (BOLD, BLUE, RESET, " " * left, body, " " * (pad - left), BOLD, BLUE, RESET)

    out.append("%s%s     ╔═══════════════════════════════════════════════════╗%s" % (BOLD, BLUE, RESET))
    out.append(box([("  ", ""), ("🧬 ", GREEN), ("Welcome to Masoudi Pipeline", LBLUE), (" v.2.0", YELLOW),
                    (" :)", YELLOW), ("  ", ""), ("🦠", MAGENTA)]))
    out.append(box([("Virology Department", "center"), (" ", YELLOW), ("●", RED), (" ", ""), ("SPH-TUMS", GREEN)]))
    out.append(box([("QC > Trim > Map > Dedup > Call > Consensus", "center")]))
    out.append("%s%s     ╚═══════════════════════════════════════════════════╝%s" % (BOLD, BLUE, RESET))
    out.append("")

    # ---- Pipeline flow ----
    out.append("  %s%s%s" % (DBLUE, "═" * 64, RESET))
    out.append("  %s%s        🔬  ONT SEQUENCING ANALYSIS PIPELINE FLOW  🔬%s" % (BOLD, LBLUE, RESET))
    out.append("  %s%s%s" % (DBLUE, "═" * 64, RESET))
    out.append("")
    steps = [
        ("STEP 1", "📋", "QC & Trimming",           "FastQC > Trim Galore > FastQC > MultiQC", MAGENTA),
        ("STEP 2", "🗺 ", "Read Mapping",            "minimap2 | BWA-MEM2 > BAM > flagstat",    LBLUE),
        ("STEP 3", "🔍", "Variant Calling",          "markdup > LoFreq + bcftools mpileup",     YELLOW),
        ("STEP 4", "🧬", "Consensus Generation",     "bcftools consensus + low-depth N mask",   GREEN),
        ("PLOTS",  "📈", "Plots & Summary",          "QC · Depth · Coverage · Variant PNGs",    BLUE),
    ]
    for i, (num, icon, name, desc, col) in enumerate(steps):
        out.append("        %s%s%s" % (DBLUE, ("┌" if i == 0 else "├") + "─" * 54 + ("┐" if i == 0 else "┤"), RESET))
        out.append("        %s│%s  %s%s%s%s  %s  %s%s%s%s  %s%s│%s" % (DBLUE, RESET, col, BOLD, num.ljust(6), RESET, icon, BOLD, WHITE,
                                                                       name.ljust(22), RESET, " " * 16, DBLUE, RESET))
        out.append("        %s│%s          %s↳  %s%s%s│%s" % (DBLUE, RESET, LBLUE, desc.ljust(41), RESET, DBLUE, RESET))
    out.append("        %s└%s┘%s" % (DBLUE, "─" * 54, RESET))
    out.append("")

    # ---- Output summary ----
    out.append("  %s%s%s" % (DBLUE, "─" * 64, RESET))
    out.append("  %s%s  📁  OUTPUT FILES%s" % (BOLD, YELLOW, RESET))
    out.append("  %s%s%s" % (DBLUE, "─" * 64, RESET))
    for fname, fdesc, col in [
        ("07_variants/",  "LoFreq + bcftools filtered VCFs",               YELLOW),
        ("08_consensus/", "Consensus FASTA (low-depth positions = N)",     GREEN),
        ("11_summary/",   "QC summary TSV + run overview",                 MAGENTA),
        ("10_plots/",     "QC, depth, coverage & variant plots (PNG)",     GREEN),
    ]:
        out.append("    %s▸%s %s%s%s%s  %s%s%s" % (col, RESET, BOLD, WHITE, fname.ljust(15), RESET, LBLUE, fdesc, RESET))
    out.append("  %s%s%s" % (DBLUE, "─" * 64, RESET))
    out.append("")

    txt = "\n".join(out) + "\n"
    try:  # consoles that cannot print box characters / emoji (e.g. cp1252): fall back to plain ASCII
        txt.encode(getattr(sys.stdout, "encoding", None) or "utf-8")
    except (UnicodeEncodeError, LookupError):
        txt = "".join(_ASCII_MAP.get(ch, ch) for ch in txt).encode("ascii", "ignore").decode("ascii")
    return txt


# --------------------------------------------------------------------------- #
#  Globals / helpers
# --------------------------------------------------------------------------- #
A = None          # parsed args
DIRS = {}         # output directories
LOGDIR = ""
CMD_FILE = ""
MAIN_LOG = ""
FQ_EXT = (".fastq", ".fq", ".fastq.gz", ".fq.gz")
q = shlex.quote


def info(msg):
    line = "[%s] %s" % (time.strftime("%H:%M:%S"), msg)
    print(("\033[1m\033[96m%s\033[0m" % line) if (msg.startswith("=== STEP") and not _NOCOLOR) else line, flush=True)
    if MAIN_LOG:
        with open(MAIN_LOG, "a") as fh:
            fh.write(line + "\n")


def die(msg):
    print("[ERROR] %s" % msg, file=sys.stderr)
    sys.exit(1)


def run(cmd, logname):
    """Run a shell command (bash, pipefail); stdout/stderr go to logs/<logname>.log."""
    with open(CMD_FILE, "a") as fh:
        fh.write(cmd + "\n")
    lp = os.path.join(LOGDIR, logname + ".log")
    with open(lp, "a") as lf:
        lf.write("$ %s\n" % cmd)
        lf.flush()
        p = subprocess.run(["bash", "-o", "pipefail", "-c", cmd], stdout=lf, stderr=subprocess.STDOUT)
    if p.returncode != 0:
        with open(lp) as lf:
            tail = "".join(lf.readlines()[-15:])
        raise RuntimeError("Command failed (exit %d): %s\n--- log tail (%s) ---\n%s" % (p.returncode, cmd, lp, tail))


def need(path):
    """True if `path` must be (re)generated."""
    if A.resume and os.path.exists(path) and os.path.getsize(path) > 0:
        return False
    return True


def sanitize(name):
    return re.sub(r"[^A-Za-z0-9._-]", "_", name)


def strip_fq_ext(fname):
    for e in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if fname.endswith(e):
            return fname[: -len(e)]
    return fname


TOOL_VERSION = OrderedDict([
    ("fastqc", "fastqc --version"),
    ("trim_galore", "trim_galore --version | grep -i version"),
    ("cutadapt", "cutadapt --version"),
    ("multiqc", "multiqc --version"),
    ("minimap2", "minimap2 --version"),
    ("bwa-mem2", "bwa-mem2 version"),
    ("samtools", "samtools --version | head -1"),
    ("bcftools", "bcftools --version | head -1"),
    ("bgzip", "bgzip --version | head -1"),
    ("lofreq", "lofreq version | head -2 | tail -1"),
])


def check_tools():
    needed = ["fastqc", "trim_galore", "cutadapt", "multiqc", "samtools", "bcftools", "bgzip"]
    needed.append("minimap2" if A.mapper == "minimap2" else "bwa-mem2")
    if "lofreq" in A.callers:
        needed.append("lofreq")
    missing = [t for t in needed if shutil.which(t) is None]
    if missing:
        die("Missing tools in PATH: %s\n        Activate the environment:  conda activate masoudi-pipeline" % ", ".join(missing))
    info("Tool versions:")
    for t in needed:
        try:
            out = subprocess.run(TOOL_VERSION[t], shell=True, capture_output=True, text=True).stdout.strip().splitlines()
            info("   %-12s %s" % (t, out[0].strip() if out else "ok"))
        except Exception:
            info("   %-12s ok" % t)


# --------------------------------------------------------------------------- #
#  Input handling
# --------------------------------------------------------------------------- #
def collect_inputs(paths, merge, merged_name):
    expanded = []
    for p in paths:
        hits = glob.glob(p)
        expanded.extend(sorted(hits) if hits else [p])
    groups = OrderedDict()
    for p in expanded:
        if os.path.isdir(p):
            files = sorted(os.path.join(p, f) for f in os.listdir(p) if f.endswith(FQ_EXT))
            if not files:
                die("No FASTQ files found in directory: %s" % p)
            name = sanitize(os.path.basename(os.path.normpath(p)))
            if name in groups:
                die("Duplicate sample name '%s'. Rename inputs or use --merge-fastq." % name)
            groups[name] = files
        elif os.path.isfile(p):
            if not p.endswith(FQ_EXT):
                die("Not a FASTQ file (.fastq/.fq/.gz): %s" % p)
            name = sanitize(strip_fq_ext(os.path.basename(p)))
            if name in groups and not merge:
                die("Duplicate sample name '%s' (%s). Rename inputs or use --merge-fastq." % (name, p))
            groups.setdefault(name, []).append(p)
        else:
            die("Input not found: %s" % p)
    if merge:
        allf = [f for fl in groups.values() for f in fl]
        return OrderedDict([(sanitize(merged_name), allf)])
    return groups


def prepare_fastq(sample, files):
    out = os.path.join(DIRS["input"], sample + ".fastq.gz")
    if not need(out):
        return out
    if os.path.lexists(out):
        os.remove(out)
    if len(files) == 1 and files[0].endswith(".gz"):
        try:
            os.symlink(os.path.abspath(files[0]), out)
        except OSError:
            shutil.copyfile(files[0], out)
        return out
    with open(out, "wb") as o:
        for f in files:
            if f.endswith(".gz"):
                with open(f, "rb") as fin:
                    shutil.copyfileobj(fin, o)      # concatenated gzip members are valid
            else:
                with open(f, "rb") as fin, gzip.GzipFile(fileobj=o, mode="wb") as g:
                    shutil.copyfileobj(fin, g)
    return out


# --------------------------------------------------------------------------- #
#  Parsers / statistics
# --------------------------------------------------------------------------- #
_ERR = np.ones(256)
_ERR[33:] = 10 ** (-(np.arange(33, 256) - 33) / 10.0)


def n50(lengths):
    if len(lengths) == 0:
        return 0
    s = np.sort(lengths)[::-1]
    c = np.cumsum(s)
    return int(s[np.searchsorted(c, c[-1] / 2.0)])


def fastq_stats(path):
    """Read length + mean read quality (error-probability based) for every read."""
    opener = gzip.open if path.endswith(".gz") else open
    lens, quals = array("I"), array("f")
    with opener(path, "rb") as fh:
        while True:
            h = fh.readline()
            if not h:
                break
            seq = fh.readline().rstrip(b"\r\n")
            fh.readline()
            qual = fh.readline().rstrip(b"\r\n")
            lens.append(len(seq))
            if not qual:
                quals.append(0.0)
                continue
            e = _ERR[np.frombuffer(qual, dtype=np.uint8)].mean()
            quals.append(float(-10.0 * np.log10(max(e, 1e-10))))
    L = np.frombuffer(lens, dtype=np.uint32).astype(np.int64) if len(lens) else np.zeros(0, dtype=np.int64)
    Q = np.frombuffer(quals, dtype=np.float32).astype(float) if len(quals) else np.zeros(0)
    return {
        "reads": int(len(L)), "bases": int(L.sum()),
        "mean_len": float(L.mean()) if len(L) else 0.0,
        "median_len": float(np.median(L)) if len(L) else 0.0,
        "n50": n50(L), "mean_q": float(Q.mean()) if len(Q) else 0.0,
        "lens": L, "quals": Q,
    }


def parse_flagstat(path):
    t = open(path).read()

    def g(pat):
        m = re.search(pat, t, re.M)
        return int(m.group(1)) if m else 0

    return {
        "total": g(r"^(\d+) \+ \d+ in total"),
        "primary": g(r"^(\d+) \+ \d+ primary$"),
        "mapped": g(r"^(\d+) \+ \d+ mapped \("),
        "primary_mapped": g(r"^(\d+) \+ \d+ primary mapped \("),
        "dups": g(r"^(\d+) \+ \d+ primary duplicates$") or g(r"^(\d+) \+ \d+ duplicates$"),
    }


def read_fai(path):
    out = OrderedDict()
    with open(path) as fh:
        for line in fh:
            f = line.split("\t")
            out[f[0]] = int(f[1])
    return out


def load_depth(path, fai):
    arrs = OrderedDict((c, np.zeros(l, dtype=np.int32)) for c, l in fai.items())
    with open(path) as fh:
        for line in fh:
            c, p, d = line.rstrip("\n").split("\t")
            if c in arrs:
                arrs[c][int(p) - 1] = int(d)
    return arrs


def read_vcf(path):
    recs = []
    if not os.path.exists(path):
        return recs
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            chrom, pos, ref, alt, qual, infos = f[0], int(f[1]), f[3], f[4].split(",")[0], f[5], f[7]
            info_d = {}
            for kv in infos.split(";"):
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    info_d[k] = v
                else:
                    info_d[kv] = True
            dp = float(info_d.get("DP", 0) or 0)
            af = None
            if "AF" in info_d:
                try:
                    af = float(str(info_d["AF"]).split(",")[0])
                except ValueError:
                    af = None
            if af is None and "DP4" in info_d:
                d4 = [float(x) for x in str(info_d["DP4"]).split(",")]
                tot = sum(d4)
                af = (d4[2] + d4[3]) / tot if tot > 0 else 0.0
            if af is None:
                af = 1.0
            if len(ref) == len(alt):
                vt = "SNV"
            elif len(ref) < len(alt):
                vt = "INS"
            else:
                vt = "DEL"
            try:
                qv = float(qual)
            except ValueError:
                qv = 0.0
            recs.append({"chrom": chrom, "pos": pos, "ref": ref, "alt": alt, "qual": qv, "dp": dp, "af": af, "type": vt})
    return recs


def read_fasta(path):
    seqs, name = OrderedDict(), None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name is not None:
                seqs[name].append(line)
    return OrderedDict((k, "".join(v)) for k, v in seqs.items())


# --------------------------------------------------------------------------- #
#  Plots
# --------------------------------------------------------------------------- #
C_RAW, C_TRIM = "#8d99ae", "#e63946"
CALLER_COL = {"lofreq": "#1d3557", "bcftools": "#e9a100"}


def plot_qc(sample, st_raw, st_trim, fate, out):
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
    mx = max(st_raw["lens"].max() if len(st_raw["lens"]) else 1, 1)
    bins = np.logspace(1, np.log10(mx + 10), 60)
    ax[0].hist(st_raw["lens"], bins=bins, color=C_RAW, alpha=0.8, label="raw (N50 %d)" % st_raw["n50"])
    if len(st_trim["lens"]):
        ax[0].hist(st_trim["lens"], bins=bins, color=C_TRIM, alpha=0.6, label="trimmed (N50 %d)" % st_trim["n50"])
    ax[0].set_xscale("log")
    ax[0].set(xlabel="Read length (bp)", ylabel="Reads", title="Read length distribution")
    ax[0].legend()
    qb = np.linspace(0, max(30, np.ceil(st_raw["quals"].max()) if len(st_raw["quals"]) else 30), 60)
    ax[1].hist(st_raw["quals"], bins=qb, color=C_RAW, alpha=0.8, label="raw (mean Q%.1f)" % st_raw["mean_q"])
    if len(st_trim["quals"]):
        ax[1].hist(st_trim["quals"], bins=qb, color=C_TRIM, alpha=0.6, label="trimmed (mean Q%.1f)" % st_trim["mean_q"])
    ax[1].set(xlabel="Mean read quality (Phred)", ylabel="Reads", title="Read quality distribution")
    ax[1].legend()
    labels = list(fate.keys())
    vals = list(fate.values())
    bars = ax[2].bar(labels, vals, color=["#8d99ae", "#e63946", "#457b9d", "#1d3557"][: len(labels)])
    for b, v in zip(bars, vals):
        ax[2].text(b.get_x() + b.get_width() / 2, v, "{:,}".format(int(v)), ha="center", va="bottom", fontsize=8)
    ax[2].set(ylabel="Reads", title="Read fate")
    ax[2].tick_params(axis="x", rotation=20)
    fig.suptitle("QC summary - %s" % sample, fontweight="bold")
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


def moving_avg(a, w):
    """Simple centred moving average (edge-padded)."""
    if w <= 1 or len(a) == 0:
        return a.astype(float)
    pad = w // 2
    ap = np.pad(a.astype(float), (pad, pad), mode="edge")
    return np.convolve(ap, np.ones(w) / w, mode="valid")[: len(a)]


def plot_depth(sample, arrs, mind, out):
    """Depth plot per sample in the VIRIONT style: one coverage line per contig, grey facet strip
    (sample name / contig name) on top, light-grey panel with white grid. Results only."""
    from matplotlib.ticker import AutoMinorLocator, FuncFormatter, MaxNLocator
    n = max(1, len(arrs))
    W, strip_h, panel_h, gap = 10.0, 0.45, 1.0, 0.12   # inches
    left, right, bottom, top = 0.75, 0.15, 0.62, 0.05
    H = top + n * (strip_h + panel_h) + (n - 1) * gap + bottom
    fig = plt.figure(figsize=(W, H))
    for i, (c, a) in enumerate(arrs.items()):
        y_top = H - top - i * (strip_h + panel_h + gap)
        ax_s = fig.add_axes([left / W, (y_top - strip_h) / H, (W - left - right) / W, strip_h / H])
        ax = fig.add_axes([left / W, (y_top - strip_h - panel_h) / H, (W - left - right) / W, panel_h / H])
        # facet strip: sample name + contig name
        ax_s.set_facecolor("#D9D9D9")
        ax_s.set_xticks([]); ax_s.set_yticks([])
        for sp in ax_s.spines.values():
            sp.set_visible(False)
        ax_s.text(0.5, 0.74, sample, ha="center", va="center", fontsize=8, color="#1A1A1A", transform=ax_s.transAxes)
        ax_s.text(0.5, 0.26, c, ha="center", va="center", fontsize=8, color="#1A1A1A", transform=ax_s.transAxes)
        # coverage line (1-based genomic position); bin only for very long contigs
        L = len(a)
        if L > 100000:
            win = int(np.ceil(L / 100000.0))
            m = (L // win) * win
            y = a[:m].reshape(-1, win).mean(axis=1)
            x = np.arange(y.size) * win + win / 2.0
        else:
            y = a.astype(float)
            x = np.arange(1, L + 1)
        ax.set_facecolor("#EBEBEB")
        ax.plot(x, y, color="#F1B0AB", lw=1.0, solid_joinstyle="round")
        ax.set_xlim(-0.05 * L + 1, 1.05 * L)
        ymax = float(y.max()) if y.size else 1.0
        ax.set_ylim(-0.05 * max(ymax, 1.0), 1.05 * max(ymax, 1.0))
        ax.yaxis.set_major_locator(MaxNLocator(nbins=4, integer=True, steps=[1, 2, 4, 5, 10]))
        ax.xaxis.set_major_locator(MaxNLocator(nbins=6, integer=True, steps=[1, 2, 5, 10]))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda v, p: "%d" % v))
        ax.xaxis.set_major_formatter(FuncFormatter(lambda v, p: "%d" % v))
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        ax.grid(which="major", color="white", lw=1.0)
        ax.grid(which="minor", color="white", lw=0.5)
        ax.set_axisbelow(True)
        for sp in ax.spines.values():
            sp.set_visible(False)
        ax.tick_params(which="major", length=2.5, width=0.6, color="#333333", labelsize=8, labelcolor="#4D4D4D")
        ax.tick_params(which="minor", length=0)
        ax.set_ylabel("Coverage", fontsize=11)
        if i < n - 1:
            ax.set_xticklabels([])
        else:
            ax.set_xlabel("Genomic position", fontsize=11)
    fig.savefig(out, dpi=200)
    plt.close(fig)


def plot_vcf(sample, calls, fai, out):
    offs, o = {}, 0
    for c, l in fai.items():
        offs[c] = o
        o += l
    total = max(o, 1)
    fig, ax = plt.subplots(2, 2, figsize=(14, 9))
    mk = {"SNV": "o", "INS": "^", "DEL": "v"}
    any_var = any(len(v) for v in calls.values())
    for caller, recs in calls.items():
        for vt, m in mk.items():
            sub = [r for r in recs if r["type"] == vt]
            if sub:
                ax[0, 0].scatter([offs.get(r["chrom"], 0) + r["pos"] for r in sub], [r["af"] for r in sub], s=22, marker=m,
                                 color=CALLER_COL.get(caller, "k"), alpha=0.75, label="%s %s" % (caller, vt))
    ax[0, 0].set(xlim=(0, total), ylim=(0, 1.05), xlabel="Genome position (bp)", ylabel="Allele frequency", title="Variants along the genome")
    if any_var:
        ax[0, 0].legend(fontsize=7, ncol=2)
    else:
        ax[0, 0].text(0.5, 0.5, "No variants passed filters", ha="center", transform=ax[0, 0].transAxes)
    bins = np.linspace(0, 1, 21)
    for caller, recs in calls.items():
        if recs:
            ax[0, 1].hist([r["af"] for r in recs], bins=bins, histtype="step", lw=2, color=CALLER_COL.get(caller, "k"), label=caller)
    ax[0, 1].set(xlabel="Allele frequency", ylabel="Variants", title="Allele-frequency spectrum")
    if any_var:
        ax[0, 1].legend()
    types = ["SNV", "INS", "DEL"]
    w = 0.8 / max(len(calls), 1)
    for i, (caller, recs) in enumerate(calls.items()):
        cnt = [sum(1 for r in recs if r["type"] == t) for t in types]
        ax[1, 0].bar(np.arange(3) + i * w, cnt, w, color=CALLER_COL.get(caller, "k"), label=caller)
    ax[1, 0].set_xticks(np.arange(3) + w * (len(calls) - 1) / 2.0)
    ax[1, 0].set_xticklabels(types)
    ax[1, 0].set(ylabel="Variants", title="Variant types")
    if any_var:
        ax[1, 0].legend()
    subs = [a + ">" + b for a in "ACGT" for b in "ACGT" if a != b]
    for i, (caller, recs) in enumerate(calls.items()):
        cnt = [sum(1 for r in recs if r["type"] == "SNV" and r["ref"].upper() + ">" + r["alt"].upper() == s) for s in subs]
        ax[1, 1].bar(np.arange(len(subs)) + i * w, cnt, w, color=CALLER_COL.get(caller, "k"), label=caller)
    ax[1, 1].set_xticks(np.arange(len(subs)) + w * (len(calls) - 1) / 2.0)
    ax[1, 1].set_xticklabels(subs, rotation=60, fontsize=8)
    ax[1, 1].set(ylabel="SNVs", title="Substitution spectrum")
    fig.suptitle("Variants (filtered) - %s" % sample, fontweight="bold")
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


def plot_overview(rows, mind, out):
    if not rows:
        return
    names = [r["sample"] for r in rows]
    x = np.arange(len(names))
    fig, ax = plt.subplots(2, 3, figsize=(16, 8))
    ax = ax.ravel()
    w = 0.4
    ax[0].bar(x - w / 2, [r["reads_raw"] for r in rows], w, color=C_RAW, label="raw")
    ax[0].bar(x + w / 2, [r["reads_trimmed"] for r in rows], w, color=C_TRIM, label="trimmed")
    ax[0].set(title="Reads", ylabel="Reads")
    ax[0].legend()
    ax[1].bar(x, [r["pct_aligned"] for r in rows], color="#457b9d")
    ax[1].set(title="% trimmed reads aligned", ylim=(0, 100))
    ax[2].bar(x, [r["mean_depth"] for r in rows], color="#1d3557")
    ax[2].set(title="Mean depth (x)")
    ax[3].bar(x, [r["breadth_min_depth"] for r in rows], color="#2a9d8f")
    ax[3].set(title="%% genome >= %dx" % mind, ylim=(0, 100))
    ax[4].bar(x - w / 2, [r["lofreq_variants"] for r in rows], w, color=CALLER_COL["lofreq"], label="lofreq")
    ax[4].bar(x + w / 2, [r["bcftools_variants"] for r in rows], w, color=CALLER_COL["bcftools"], label="bcftools")
    ax[4].set(title="Filtered variants")
    ax[4].legend()
    ax[5].bar(x, [r["consensus_N_pct"] for r in rows], color="#6c757d")
    ax[5].set(title="% N in consensus")
    for a in ax:
        a.set_xticks(x)
        a.set_xticklabels(names, rotation=45, ha="right", fontsize=8)
    fig.suptitle("masoudi-pipeline - run overview", fontweight="bold")
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


# --------------------------------------------------------------------------- #
#  Extra plots ported from ont_analysis_pipeline.py: coverage (bar / scatter /
#  heatmap / depth-of-coverage / per-reference) and variant plots.
#  numpy + matplotlib only (no pandas / seaborn / pysam needed).
#  Files are written to 10_plots/ as <sample>.<name>.png (variant plots: <sample>.<caller>.<name>.png)
# --------------------------------------------------------------------------- #
CATEGORY_COLORS = {"Human": "#4C72B0", "Virus": "#C44E52", "Other": "#8C8C8C"}


def classify_reference(ref_id, description=""):
    """Heuristic Human / Virus / Other classification from the FASTA header (offline keyword match)."""
    text = ("%s %s" % (ref_id, description)).lower()
    human_patterns = ["homo sapiens", "human", "grch", "hg19", "hg38", "chromosome"]
    chrom_regex = re.compile(r"^(chr)?([1-9]|1[0-9]|2[0-2]|x|y|mt|m)$", re.IGNORECASE)
    if any(p in text for p in human_patterns) or chrom_regex.match(ref_id.strip()):
        return "Human"
    virus_patterns = [
        "virus", "viral", "viridae", "phage", "sars", "cov-2", "coronavirus",
        "influenza", "hiv", "hepatitis", "hbv", "hcv", "papillomavirus", "hpv",
        "norovirus", "rotavirus", "adenovirus", "rsv", "monkeypox", "mpox",
        "dengue", "zika", "ebola", "rabies", "poliovirus", "measles", "mumps",
        "rubella", "enterovirus", "rhinovirus", "parvovirus", "herpesvirus",
        "cytomegalovirus", "segment",
    ]
    if any(p in text for p in virus_patterns):
        return "Virus"
    return "Other"


def build_ref_info(ref, fai):
    """id -> {description, category, length} from the reference FASTA headers."""
    ri = OrderedDict()
    with open(ref) as fh:
        for line in fh:
            if line.startswith(">"):
                hdr = line[1:].strip()
                rid = hdr.split()[0] if hdr else ""
                ri[rid] = {"description": hdr, "category": classify_reference(rid, hdr), "length": fai.get(rid, 0)}
    return ri


def _safe_ref(rn):
    return rn.replace("|", "_").replace("/", "_").replace(" ", "_")[:50]


def read_coverage_tsv(path):
    """(rname list, meandepth list) from `samtools coverage` output."""
    names, means = [], []
    if not os.path.exists(path):
        return names, means
    with open(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        if "#rname" not in head or "meandepth" not in head:
            return names, means
        i_n, i_m = head.index("#rname"), head.index("meandepth")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) > max(i_n, i_m):
                names.append(f[i_n])
                means.append(float(f[i_m]))
    return names, means


def plot_coverage_annotated(pre, covtsv, arrs, ref_info, mind):
    """Bar plot + scatter plot of mean depth per reference (coloured Human/Virus/Other),
    then the genome-wide heatmap and the depth-of-coverage plots."""
    info("[%s] plotting coverage (genome-type annotated)" % os.path.basename(pre).rstrip("."))
    rnames, meandepth = read_coverage_tsv(covtsv)
    if not rnames:
        info("No coverage data to plot")
        return
    cats = [ref_info.get(r, {}).get("category", "Other") for r in rnames]
    bar_colors = [CATEGORY_COLORS[c] for c in cats]

    plt.figure(figsize=(max(12, len(rnames) * 0.5), 6))
    plt.bar(range(len(rnames)), meandepth, color=bar_colors, edgecolor="black", alpha=0.85)
    plt.xlabel("Reference Sequence")
    plt.ylabel("Mean Depth")
    plt.title("Coverage Across Reference Sequences (annotated by genome type)")
    plt.xticks(range(len(rnames)), ["%s [%s]" % (r, c) for r, c in zip(rnames, cats)], rotation=45, ha="right")
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in CATEGORY_COLORS.values()]
    plt.legend(handles, CATEGORY_COLORS.keys(), title="Genome Type")
    plt.tight_layout()
    plt.savefig(pre + "coverage_bar_plot.png", dpi=300)
    plt.close()

    if len(rnames) > 1:
        plt.figure(figsize=(10, 6))
        plt.scatter(range(len(rnames)), meandepth, s=100, alpha=0.7, edgecolors="black", c=bar_colors)
        plt.xlabel("Reference Sequence Index")
        plt.ylabel("Mean Depth")
        plt.title("Mean Depth Distribution")
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(pre + "coverage_scatter.png", dpi=300)
        plt.close()

    plot_coverage_heatmap(pre, arrs, ref_info)
    plot_depth_of_coverage(pre, arrs, ref_info, mind)


def plot_coverage_heatmap(pre, arrs, ref_info, n_bins=100):
    """Heatmap of mean depth binned along each reference; rows grouped / colour-labelled by genome type."""
    if not arrs:
        info("No depth data available for coverage heatmap")
        return
    contigs = sorted(arrs.keys(), key=lambda c: (ref_info.get(c, {}).get("category", "Other"), c))
    matrix, row_labels = [], []
    for contig in contigs:
        a = arrs[contig]
        length = ref_info.get(contig, {}).get("length", len(a)) or len(a)
        pos = np.arange(1, len(a) + 1)
        bins = np.linspace(0, max(length, 1), n_bins + 1)
        bin_idx = np.clip(np.digitize(pos, bins) - 1, 0, n_bins - 1)
        depth_sum = np.zeros(n_bins)
        depth_count = np.zeros(n_bins)
        np.add.at(depth_sum, bin_idx, a)
        np.add.at(depth_count, bin_idx, 1)
        mean_depth = np.divide(depth_sum, depth_count, out=np.zeros_like(depth_sum), where=depth_count > 0)
        matrix.append(mean_depth)
        row_labels.append("%s [%s]" % (contig, ref_info.get(contig, {}).get("category", "Other")))
    log_matrix = np.log1p(np.array(matrix))

    fig, ax = plt.subplots(figsize=(14, max(4, 0.4 * len(contigs) + 2)))
    im = ax.imshow(log_matrix, cmap="viridis", aspect="auto", interpolation="nearest")
    fig.colorbar(im, ax=ax).set_label("log(1 + Mean Depth)")
    ax.set_yticks(range(len(contigs)))
    ax.set_yticklabels(row_labels, rotation=0)
    ax.set_xticks([])
    ax.set_xlabel("Genome Position (binned, 5' -> 3')")
    ax.set_ylabel("Reference Sequence")
    ax.set_title("Coverage Heatmap Across Reference Genomes (annotated by type)")
    for tick_label, contig in zip(ax.get_yticklabels(), contigs):
        tick_label.set_color(CATEGORY_COLORS.get(ref_info.get(contig, {}).get("category", "Other"), "black"))
    fig.tight_layout()
    fig.savefig(pre + "coverage_heatmap.png", dpi=300)
    plt.close(fig)


def plot_depth_of_coverage(pre, arrs, ref_info, mind, max_individual_plots=25):
    """Per-position depth for each reference (one PNG each) + a combined overview grid."""
    if not arrs:
        info("No depth data available for depth-of-coverage plots")
        return
    contigs = sorted(arrs.keys())
    if len(contigs) > max_individual_plots:
        info("%d reference sequences; plotting depth-of-coverage for the first %d" % (len(contigs), max_individual_plots))
        contigs = contigs[:max_individual_plots]

    for contig in contigs:
        a = arrs[contig]
        pos = np.arange(1, len(a) + 1)
        category = ref_info.get(contig, {}).get("category", "Other")
        color = CATEGORY_COLORS.get(category, "#8C8C8C")
        plt.figure(figsize=(12, 5))
        plt.fill_between(pos, a, color=color, alpha=0.4)
        plt.plot(pos, a, color=color, linewidth=0.8)
        plt.axhline(mind, color="red", linestyle="--", linewidth=1, label="Min depth threshold (%dx)" % mind)
        plt.axhline(a.mean(), color="black", linestyle=":", linewidth=1, label="Mean depth (%.1fx)" % a.mean())
        plt.xlabel("Genome Position (bp)")
        plt.ylabel("Read Depth")
        plt.title("Depth of Coverage: %s [%s]" % (contig, category))
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(pre + "depth_coverage_%s.png" % re.sub(r"[^A-Za-z0-9_.-]+", "_", contig), dpi=300)
        plt.close()

    n = len(contigs)
    ncols = 2 if n > 1 else 1
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(14, 3.2 * nrows), squeeze=False)
    for i, contig in enumerate(contigs):
        ax = axes[i // ncols][i % ncols]
        a = arrs[contig]
        pos = np.arange(1, len(a) + 1)
        category = ref_info.get(contig, {}).get("category", "Other")
        color = CATEGORY_COLORS.get(category, "#8C8C8C")
        ax.fill_between(pos, a, color=color, alpha=0.4)
        ax.plot(pos, a, color=color, linewidth=0.6)
        ax.axhline(mind, color="red", linestyle="--", linewidth=0.8)
        ax.set_title("%s [%s]" % (contig, category), fontsize=9)
        ax.set_xlabel("Position (bp)", fontsize=8)
        ax.set_ylabel("Depth", fontsize=8)
        ax.tick_params(labelsize=7)
    for j in range(n, nrows * ncols):
        axes[j // ncols][j % ncols].axis("off")
    fig.suptitle("Depth-of-Coverage Overview (all reference sequences)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(pre + "depth_coverage_overview.png", dpi=300)
    plt.close(fig)


def plot_coverage_per_reference(pre, arrs, mind):
    """Per-reference coverage line plot, binned heatmap and depth-distribution histogram."""
    for rn, a in arrs.items():
        if len(a) == 0:
            continue
        pos, dep = np.arange(1, len(a) + 1), a
        safe = _safe_ref(rn)

        plt.figure(figsize=(14, 6))
        plt.plot(pos, dep, linewidth=0.5, alpha=0.7)
        plt.axhline(y=mind, color="r", linestyle="--", label="Min depth (%dx)" % mind)
        plt.xlabel("Position (bp)")
        plt.ylabel("Depth")
        plt.title("Coverage - %s" % rn)
        plt.legend()
        plt.tight_layout()
        plt.savefig(pre + "coverage_line_%s.png" % safe, dpi=300)
        plt.close()

        bsz = max(1, len(pos) // 100)
        binned = [np.mean(dep[i:i + bsz]) for i in range(0, len(dep), bsz)]
        cols = 10
        rows = (len(binned) + cols - 1) // cols
        hm = np.zeros((rows, cols))
        for i, d in enumerate(binned):
            if i // cols < rows:
                hm[i // cols, i % cols] = d
        fig, ax = plt.subplots(figsize=(12, max(6, rows * 0.5)))
        im = ax.imshow(hm, cmap="YlOrRd", aspect="auto", interpolation="nearest")
        fig.colorbar(im, ax=ax).set_label("Depth")
        ax.set_xticks(range(cols))
        ax.set_yticks(range(0, rows, max(1, rows // 20)))
        ax.set_title("Coverage Heatmap - %s" % rn)
        fig.tight_layout()
        fig.savefig(pre + "coverage_heatmap_%s.png" % safe, dpi=300)
        plt.close(fig)

        plt.figure(figsize=(10, 6))
        plt.hist(dep, bins=50, edgecolor="black", alpha=0.7)
        plt.axvline(x=mind, color="r", linestyle="--", label="Min depth (%dx)" % mind)
        plt.xlabel("Depth")
        plt.ylabel("Frequency")
        plt.title("Depth Distribution - %s" % rn)
        plt.legend()
        plt.yscale("log")
        plt.tight_layout()
        plt.savefig(pre + "depth_distribution_%s.png" % safe, dpi=300)
        plt.close()


def plot_variants_extra(pre, recs):
    """Overall variant position/quality scatter + quality histogram, then per-reference
    distribution, type pie chart and quality histogram. `recs` = read_vcf() records."""
    from collections import Counter
    if not recs:
        info("No variants found for plotting")
        return
    pos_all = [r["pos"] for r in recs]
    qual_all = [r["qual"] for r in recs]

    plt.figure(figsize=(12, 6))
    plt.scatter(pos_all, qual_all, alpha=0.6, edgecolors="black")
    plt.xlabel("Position")
    plt.ylabel("Quality")
    plt.title("Variant Positions and Quality")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(pre + "variant_positions.png", dpi=300)
    plt.close()

    plt.figure(figsize=(10, 6))
    plt.hist(qual_all, bins=30, edgecolor="black", alpha=0.7, color="coral")
    plt.xlabel("Quality Score")
    plt.ylabel("Frequency")
    plt.title("Variant Quality Distribution")
    plt.tight_layout()
    plt.savefig(pre + "variant_quality_dist.png", dpi=300)
    plt.close()

    tname = {"SNV": "SNP", "INS": "Insertion", "DEL": "Deletion"}
    by_ref = OrderedDict()
    for r in recs:
        d = by_ref.setdefault(r["chrom"], {"positions": [], "types": [], "qualities": []})
        d["positions"].append(r["pos"])
        d["types"].append(tname.get(r["type"], r["type"]))
        if r["qual"]:
            d["qualities"].append(r["qual"])
    for rn, data in by_ref.items():
        pos, types, quals = data["positions"], data["types"], data["qualities"]
        safe = _safe_ref(rn)

        plt.figure(figsize=(14, 6))
        for vt, col in [("SNP", "blue"), ("Insertion", "green"), ("Deletion", "red")]:
            vpos = [p for p, t in zip(pos, types) if t == vt]
            if vpos:
                plt.scatter(vpos, [1] * len(vpos), label=vt, alpha=0.6, s=20, color=col)
        plt.xlabel("Position (bp)")
        plt.ylabel("Variants")
        plt.title("Variant Distribution - %s" % rn)
        plt.legend()
        plt.yticks([])
        plt.tight_layout()
        plt.savefig(pre + "variant_distribution_%s.png" % safe, dpi=300)
        plt.close()

        tc = Counter(types).most_common()
        plt.figure(figsize=(8, 8))
        plt.pie([v for _, v in tc], labels=[k for k, _ in tc], autopct="%1.1f%%", startangle=90)
        plt.title("Variant Types - %s" % rn)
        plt.tight_layout()
        plt.savefig(pre + "variant_types_%s.png" % safe, dpi=300)
        plt.close()

        if quals:
            plt.figure(figsize=(10, 6))
            plt.hist(quals, bins=30, edgecolor="black", alpha=0.7)
            plt.xlabel("Quality")
            plt.ylabel("Frequency")
            plt.title("Variant Quality - %s" % rn)
            plt.tight_layout()
            plt.savefig(pre + "variant_quality_%s.png" % safe, dpi=300)
            plt.close()


def plot_ont_extra(sample, plots_dir, ref, fai, covtsv, arrs, calls, mind):
    """Run all the extra plots for one sample; never lets a plotting error fail the sample."""
    try:
        ref_info = build_ref_info(ref, fai)
        pre = os.path.join(plots_dir, sample + ".")
        plot_coverage_annotated(pre, covtsv, arrs, ref_info, mind)
        plot_coverage_per_reference(pre, arrs, mind)
        for caller, recs in calls.items():
            if recs:
                plot_variants_extra(os.path.join(plots_dir, "%s.%s." % (sample, caller)), recs)
    except Exception as e:  # noqa: BLE001
        info("[%s] WARNING: extra plots failed: %s" % (sample, e))


# --------------------------------------------------------------------------- #
#  Filter expressions (bcftools syntax, used for BOTH callers)
# --------------------------------------------------------------------------- #
def af_expr(caller, thr):
    if caller == "lofreq":
        return "INFO/AF>=%s" % thr
    return "(INFO/DP4[2]+INFO/DP4[3])/(INFO/DP4[0]+INFO/DP4[1]+INFO/DP4[2]+INFO/DP4[3])>=%s" % thr


def filter_expr(caller):
    return "QUAL>=%s && INFO/DP>=%s && %s" % (A.min_qual, A.min_depth, af_expr(caller, A.min_af))


# --------------------------------------------------------------------------- #
#  PHASE A  (step 1): QC -> trim -> QC -> MultiQC
# --------------------------------------------------------------------------- #
def phase_a(sample, fq):
    t = A.threads
    raw_html = os.path.join(DIRS["qc_raw"], sample + "_fastqc.html")
    if need(raw_html):
        info("[%s] FastQC (raw)" % sample)
        run("fastqc -t %d -q --outdir %s %s" % (t, q(DIRS["qc_raw"]), q(fq)), sample + ".fastqc_raw")
    trimmed = os.path.join(DIRS["trim"], sample + "_trimmed.fq.gz")
    if need(trimmed):
        info("[%s] Trim Galore" % sample)
        cmd = ["trim_galore", "--phred33", "--gzip", "--quality", A.trim_quality, "--length", A.min_length,
               "--stringency", A.stringency, "--cores", min(4, t), "--basename", sample, "-o", DIRS["trim"]]
        if A.max_length:
            cmd += ["--max_length", A.max_length]
        if A.clip5:
            cmd += ["--clip_R1", A.clip5]
        if A.clip3:
            cmd += ["--three_prime_clip_R1", A.clip3]
        if A.adapter:
            cmd += ["--adapter", A.adapter]
        cmd += [fq]
        run(" ".join(q(str(c)) for c in cmd), sample + ".trim_galore")
    if not os.path.exists(trimmed):
        raise RuntimeError("Trim Galore did not produce %s" % trimmed)
    trim_html = os.path.join(DIRS["qc_trim"], sample + "_trimmed_fastqc.html")
    if need(trim_html):
        info("[%s] FastQC (trimmed, validation)" % sample)
        run("fastqc -t %d -q --outdir %s %s" % (t, q(DIRS["qc_trim"]), q(trimmed)), sample + ".fastqc_trimmed")
    info("[%s] computing read statistics" % sample)
    return trimmed, fastq_stats(fq), fastq_stats(trimmed)


# --------------------------------------------------------------------------- #
#  PHASE B  (steps 2-4)
# --------------------------------------------------------------------------- #
def phase_b(sample, trimmed, st_raw, st_trim, ref, fai):
    t = A.threads
    log = lambda n: sample + "." + n  # noqa: E731
    D = DIRS
    sam = os.path.join(D["map"], sample + ".sam")
    sbam = os.path.join(D["map"], sample + ".sorted.bam")
    raw_fs = os.path.join(D["map"], sample + ".raw_alignments.flagstat.txt")
    sorted_fs = os.path.join(D["map"], sample + ".sorted.flagstat.txt")

    # ---------------- STEP 2: mapping ---------------- #
    if need(sam) and need(sbam):
        info("[%s] mapping with %s" % (sample, A.mapper))
        rg = "@RG\\tID:%s\\tSM:%s\\tPL:ONT" % (sample, sample)
        if A.mapper == "minimap2":
            cmd = "minimap2 -ax %s -t %d --secondary=no %s -R %s %s %s > %s" % (
                q(A.minimap2_preset), t, A.mapper_extra, q(rg), q(ref), q(trimmed), q(sam))
        else:
            cmd = "bwa-mem2 mem -x ont2d -t %d %s -R %s %s %s > %s" % (t, A.mapper_extra, q(rg), q(ref), q(trimmed), q(sam))
        run(cmd, log("mapping"))
    if need(raw_fs) and os.path.exists(sam):
        run("samtools flagstat -@ %d %s > %s" % (t, q(sam), q(raw_fs)), log("flagstat_raw"))
    if need(sbam):
        info("[%s] filter + sort + index BAM" % sample)
        run("samtools view -@ %d -u -F %d -q %d %s | samtools sort -@ %d -o %s -" % (
            t, A.exclude_flags, A.min_mapq, q(sam), t, q(sbam)), log("sort"))
    if need(sbam + ".bai"):
        run("samtools index -@ %d %s" % (t, q(sbam)), log("index"))
    if need(sorted_fs):
        run("samtools flagstat -@ %d %s > %s" % (t, q(sbam), q(sorted_fs)), log("flagstat_sorted"))
        run("samtools idxstats %s > %s" % (q(sbam), q(os.path.join(D["map"], sample + ".idxstats.txt"))), log("idxstats"))

    # ---------------- STEP 3a: duplicates ---------------- #
    fbam = os.path.join(D["dedup"], sample + ".final.bam")
    final_fs = os.path.join(D["dedup"], sample + ".final.flagstat.txt")
    if need(fbam):
        if A.dedup == "none":
            shutil.copyfile(sbam, fbam)
        else:
            info("[%s] samtools markdup (%s duplicates)" % (sample, A.dedup))
            rm = "-r " if A.dedup == "remove" else ""
            run("samtools sort -n -@ %d %s | samtools fixmate -m -@ %d - - | samtools sort -@ %d - | "
                "samtools markdup %s-@ %d -f %s - %s" % (
                    t, q(sbam), t, t, rm, t, q(os.path.join(D["dedup"], sample + ".markdup.stats.txt")), q(fbam)), log("markdup"))
    if need(fbam + ".bai"):
        run("samtools index -@ %d %s" % (t, q(fbam)), log("index_final"))
    if need(final_fs):
        run("samtools flagstat -@ %d %s > %s" % (t, q(fbam), q(final_fs)), log("flagstat_final"))
        run("samtools stats %s > %s" % (q(fbam), q(os.path.join(D["dedup"], sample + ".samtools.stats.txt"))), log("samtools_stats"))

    # ---------------- coverage / depth (samtools mpileup) ---------------- #
    depth = os.path.join(D["cov"], sample + ".depth.tsv")
    if need(depth):
        info("[%s] depth with samtools mpileup" % sample)
        run("samtools mpileup -aa -d 0 -Q 0 -q 0 -B -f %s %s | cut -f 1,2,4 > %s" % (q(ref), q(fbam), q(depth)), log("mpileup_depth"))
    covtsv = os.path.join(D["cov"], sample + ".coverage.tsv")
    if need(covtsv):
        run("samtools coverage %s > %s" % (q(fbam), q(covtsv)), log("coverage"))
    arrs = load_depth(depth, fai)
    mask = os.path.join(D["cov"], sample + ".lowdepth_mask.bed")
    n_mask = 0
    with open(mask, "w") as o:
        for c, a in arrs.items():
            low = (a < A.min_depth).astype(np.int8)
            if not low.any():
                continue
            d = np.diff(np.concatenate(([0], low, [0])))
            for s_, e_ in zip(np.where(d == 1)[0], np.where(d == -1)[0]):
                o.write("%s\t%d\t%d\n" % (c, s_, e_))
                n_mask += 1

    # ---------------- STEP 3b: SNP calling ---------------- #
    filtered = {}
    for caller in A.callers:
        raw_vcf = os.path.join(D["var"], "%s.%s.raw.vcf" % (sample, caller))
        filt_vcf = os.path.join(D["var"], "%s.%s.filtered.vcf" % (sample, caller))
        if need(raw_vcf):
            info("[%s] SNP calling: %s" % (sample, caller))
            if caller == "lofreq":
                src = fbam
                if A.indels:
                    src = os.path.join(D["var"], sample + ".indelqual.bam")
                    run("lofreq indelqual --dindel -f %s -o %s %s" % (q(ref), q(src), q(fbam)), log("lofreq_indelqual"))
                    run("samtools index %s" % q(src), log("lofreq_index"))
                cmd = "lofreq call-parallel --pp-threads %d -f %s -o %s -q %d -Q %d -m %d -C %d -B %s %s" % (
                    t, q(ref), q(raw_vcf), A.min_bq, A.min_bq, A.min_mapq, A.min_depth,
                    "--call-indels" if A.indels else "", q(src))
                run(cmd, log("lofreq_call"))
            else:
                skip = "" if A.indels else "-V indels"
                cmd = ("bcftools mpileup -f %s -B -Q %d -q %d -d 1000000 -L 1000000 -a FORMAT/AD,FORMAT/DP -Ou %s | "
                       "bcftools call -mv --ploidy 1 %s -Ov -o %s" % (q(ref), A.min_bq, A.min_mapq, q(fbam), skip, q(raw_vcf)))
                run(cmd, log("bcftools_call"))
        # reheader from the reference .fai so every contig is declared (lofreq omits ##contig lines)
        run("bcftools reheader -f %s %s | bcftools view -i %s -Ov -o %s" % (
            q(ref + ".fai"), q(raw_vcf), q(filter_expr(caller)), q(filt_vcf)), log(caller + "_filter"))
        run("bcftools stats %s > %s" % (q(filt_vcf), q(os.path.join(D["var"], "%s.%s.stats.txt" % (sample, caller)))), log(caller + "_stats"))
        filtered[caller] = filt_vcf

    # ---------------- STEP 4: consensus ---------------- #
    cc = A.consensus_from if A.consensus_from in filtered else list(filtered)[0]
    cons_vcf = os.path.join(D["var"], "%s.%s.consensus_sites.vcf.gz" % (sample, cc))
    cons_fa = os.path.join(D["cons"], sample + ".consensus.fasta")
    final_vcf = os.path.join(D["var"], sample + ".final.vcf")
    info("[%s] consensus with bcftools (variants from %s, AF>=%s, depth<%d masked)" % (sample, cc, A.consensus_af, A.min_depth))
    run("bcftools view -i %s %s | bcftools norm -f %s -c w -Oz -o %s" % (
        q(af_expr(cc, A.consensus_af)), q(filtered[cc]), q(ref), q(cons_vcf)), log("consensus_prep"))
    run("bcftools index -f %s" % q(cons_vcf), log("consensus_index"))
    run("bcftools view %s -Ov -o %s" % (q(cons_vcf), q(final_vcf)), log("final_vcf"))
    tmp_fa = cons_fa + ".tmp"
    # Apply variants WITHOUT bcftools' -m (it segfaults on some builds); mask low-depth sites ourselves.
    run("bcftools consensus -f %s -o %s %s" % (q(ref), q(tmp_fa), q(cons_vcf)), log("consensus"))
    seqs = read_fasta(tmp_fa)
    if n_mask:
        if all(c in seqs and len(seqs[c]) == l for c, l in fai.items()):
            # SNP-only consensus keeps reference coordinates -> mask directly with N
            for c, a in arrs.items():
                lst = list(seqs[c])
                for i in np.where(a < A.min_depth)[0]:
                    lst[i] = "N"
                seqs[c] = "".join(lst)
        else:
            # indels shifted coordinates: fall back to bcftools' own masking
            run("bcftools consensus -f %s -m %s -o %s %s" % (q(ref), q(mask), q(tmp_fa), q(cons_vcf)), log("consensus_mask"))
            seqs = read_fasta(tmp_fa)
    os.remove(tmp_fa)
    nN = tot = 0
    with open(cons_fa, "w") as o:
        for c, s in seqs.items():
            hdr = sample if len(seqs) == 1 else "%s|%s" % (sample, c)
            o.write(">%s\n" % hdr)
            for i in range(0, len(s), 80):
                o.write(s[i:i + 80] + "\n")
            nN += s.upper().count("N")
            tot += len(s)

    # ---------------- stats + plots (depth plot only, for the consensus) ---------------- #
    fr, ff, fs = parse_flagstat(raw_fs), parse_flagstat(final_fs), parse_flagstat(sorted_fs)
    reads_aligned = fs["mapped"]
    allv = np.concatenate(list(arrs.values())) if arrs else np.zeros(1)
    calls = OrderedDict((c, read_vcf(p)) for c, p in filtered.items())
    row = OrderedDict([
        ("sample", sample),
        ("reads_raw", st_raw["reads"]), ("bases_raw", st_raw["bases"]), ("mean_len_raw", round(st_raw["mean_len"], 1)),
        ("n50_raw", st_raw["n50"]), ("meanQ_raw", round(st_raw["mean_q"], 2)),
        ("reads_trimmed", st_trim["reads"]), ("bases_trimmed", st_trim["bases"]), ("mean_len_trimmed", round(st_trim["mean_len"], 1)),
        ("n50_trimmed", st_trim["n50"]), ("meanQ_trimmed", round(st_trim["mean_q"], 2)),
        ("pct_reads_retained", round(100.0 * st_trim["reads"] / st_raw["reads"], 2) if st_raw["reads"] else 0),
        ("reads_aligned_filtered", reads_aligned),
        ("pct_aligned", round(100.0 * fr["primary_mapped"] / fr["primary"], 2) if fr["primary"] else 0),
        ("reads_final_bam", ff["mapped"]), ("duplicates_flagged", fs["mapped"] - ff["mapped"] if A.dedup == "remove" else ff["dups"]),
        ("mean_depth", round(float(allv.mean()), 2)), ("median_depth", round(float(np.median(allv)), 2)),
        ("breadth_1x", round(100.0 * float((allv >= 1).mean()), 2)),
        ("breadth_min_depth", round(100.0 * float((allv >= A.min_depth).mean()), 2)),
        ("lofreq_variants", len(calls.get("lofreq", []))), ("bcftools_variants", len(calls.get("bcftools", []))),
        ("consensus_len", tot), ("consensus_N", nN), ("consensus_N_pct", round(100.0 * nN / tot, 2) if tot else 0),
    ])
    P = D["plots"]
    plot_qc(sample, st_raw, st_trim, OrderedDict([("raw", st_raw["reads"]), ("trimmed", st_trim["reads"]),
                                                  ("aligned", reads_aligned), ("final BAM", ff["mapped"])]),
            os.path.join(P, sample + ".qc.png"))
    plot_depth(sample, arrs, A.min_depth, os.path.join(P, sample + ".depth.png"))
    plot_vcf(sample, calls, fai, os.path.join(P, sample + ".vcf.png"))
    plot_ont_extra(sample, P, ref, fai, covtsv, arrs, calls, A.min_depth)
    return row


# --------------------------------------------------------------------------- #
#  main
# --------------------------------------------------------------------------- #
def parse_args():
    ap = argparse.ArgumentParser(
        prog="masoudi_pipeline.py",
        description="masoudi-pipeline: Oxford Nanopore reference-based QC -> mapping -> variants -> consensus",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""examples:
  python masoudi_pipeline.py -i reads.fastq.gz -r ref.fasta -o results -t 8
  python masoudi_pipeline.py -i s1.fastq.gz s2.fastq.gz s3.fastq.gz -r ref.fasta -o results
  python masoudi_pipeline.py -i fastq_pass/barcode01 fastq_pass/barcode02 -r ref.fasta -o results
  python masoudi_pipeline.py -i run/*.fastq.gz --merge-fastq --sample-name patient1 -r ref.fasta -o results
""")
    g = ap.add_argument_group("input / output")
    g.add_argument("-i", "--input", nargs="+", required=True, help="FASTQ file(s), glob(s) or barcode folder(s). Each file/folder = one sample")
    g.add_argument("-r", "--reference", required=True, help="Custom reference FASTA (one or more contigs)")
    g.add_argument("-o", "--outdir", default="masoudi_results", help="Output directory [masoudi_results]")
    g.add_argument("-t", "--threads", type=int, default=8, help="Threads [8]")
    g.add_argument("--merge-fastq", action="store_true", help="Concatenate ALL inputs into ONE sample")
    g.add_argument("--sample-name", default="merged", help="Sample name used with --merge-fastq [merged]")
    g.add_argument("--resume", action="store_true", help="Skip steps whose output files already exist")
    g = ap.add_argument_group("step 1: Trim Galore")
    g.add_argument("--trim-quality", type=int, default=10, help="Trim Galore --quality (Phred) [10]")
    g.add_argument("--min-length", type=int, default=200, help="Trim Galore --length: discard reads shorter than N bp [200]")
    g.add_argument("--max-length", type=int, default=0, help="Trim Galore --max_length: discard reads longer than N bp [off]")
    g.add_argument("--stringency", type=int, default=3, help="Trim Galore --stringency: min adapter overlap [3]")
    g.add_argument("--clip5", type=int, default=0, help="Trim N bp from the 5' end of every read [0]")
    g.add_argument("--clip3", type=int, default=0, help="Trim N bp from the 3' end of every read [0]")
    g.add_argument("--adapter", default="", help="Custom adapter sequence for Trim Galore [auto-detect]")
    g = ap.add_argument_group("step 2: mapping")
    g.add_argument("--mapper", choices=["minimap2", "bwa-mem2"], default="minimap2", help="Mapper [minimap2]")
    g.add_argument("--minimap2-preset", default="map-ont", help="minimap2 -x preset [map-ont]")
    g.add_argument("--mapper-extra", default="", help="Extra arguments passed verbatim to the mapper")
    g.add_argument("--min-mapq", type=int, default=10, help="samtools view -q: minimum mapping quality [10]")
    g.add_argument("--exclude-flags", type=int, default=2308, help="samtools view -F (2308 = unmapped+secondary+supplementary) [2308]")
    g.add_argument("--no-sam", action="store_true", help="Delete the (large) SAM file at the end")
    g = ap.add_argument_group("step 3: duplicates + variant calling")
    g.add_argument("--dedup", choices=["remove", "mark", "none"], default="remove", help="samtools markdup mode [remove]")
    g.add_argument("--callers", nargs="+", choices=["lofreq", "bcftools"], default=["lofreq", "bcftools"], help="Variant callers [lofreq bcftools]")
    g.add_argument("--min-bq", type=int, default=7, help="Minimum base quality for calling [7]")
    g.add_argument("--min-depth", type=int, default=10, help="Min depth: VCF filter + consensus N-masking [10]")
    g.add_argument("--min-qual", type=float, default=20, help="Min variant QUAL [20]")
    g.add_argument("--min-af", type=float, default=0.10, help="Min allele frequency for the filtered VCF [0.10]")
    g.add_argument("--indels", action="store_true", help="Also call indels (off by default: ONT indel errors in homopolymers)")
    g = ap.add_argument_group("step 4: consensus")
    g.add_argument("--consensus-from", choices=["lofreq", "bcftools"], default="lofreq", help="Caller used for consensus [lofreq]")
    g.add_argument("--consensus-af", type=float, default=0.5, help="Min AF of a variant to be written into the consensus [0.5]")
    g.add_argument("--version", action="version", version="masoudi-pipeline " + __version__)
    return ap.parse_args()


def main():
    global A, DIRS, LOGDIR, CMD_FILE, MAIN_LOG
    A = parse_args()
    print(banner())
    out = os.path.abspath(A.outdir)
    names = OrderedDict([
        ("ref", "00_reference"), ("input", "00_input"), ("qc_raw", "01_qc_raw"), ("trim", "02_trimmed"),
        ("qc_trim", "03_qc_trimmed"), ("mqc1", "04_multiqc_trimming"), ("map", "05_mapping"), ("dedup", "06_dedup"),
        ("var", "07_variants"), ("cons", "08_consensus"), ("cov", "09_coverage"), ("plots", "10_plots"),
        ("summary", "11_summary")])
    for k, v in names.items():
        DIRS[k] = os.path.join(out, v)
        os.makedirs(DIRS[k], exist_ok=True)
    LOGDIR = os.path.join(out, "logs")
    os.makedirs(LOGDIR, exist_ok=True)
    CMD_FILE = os.path.join(out, "commands.sh")
    MAIN_LOG = os.path.join(out, "masoudi_pipeline.log")
    if not os.path.exists(CMD_FILE):
        with open(CMD_FILE, "w") as fh:
            fh.write("#!/usr/bin/env bash\n# commands executed by masoudi-pipeline %s\nset -euo pipefail\n" % __version__)
    info("masoudi-pipeline %s | outdir: %s | threads: %d" % (__version__, out, A.threads))
    info("command line: " + " ".join(shlex.quote(a) for a in sys.argv))
    if not os.path.isfile(A.reference):
        die("Reference FASTA not found: %s" % A.reference)
    check_tools()

    samples = collect_inputs(A.input, A.merge_fastq, A.sample_name)
    info("%d sample(s): %s" % (len(samples), ", ".join(samples)))

    # reference + indexes
    ref = os.path.join(DIRS["ref"], "reference.fasta")
    if need(ref):
        shutil.copyfile(A.reference, ref)
    if need(ref + ".fai"):
        run("samtools faidx %s" % q(ref), "reference_index")
    if A.mapper == "bwa-mem2" and need(ref + ".bwt.2bit.64"):
        info("indexing reference for BWA-MEM2")
        run("bwa-mem2 index %s" % q(ref), "reference_bwa_index")
    fai = read_fai(ref + ".fai")
    info("reference: %d contig(s), %s bp" % (len(fai), "{:,}".format(sum(fai.values()))))

    # ---------- STEP 1 ----------
    info("=== STEP 1: FastQC -> Trim Galore -> FastQC -> MultiQC ===")
    prep, failed = OrderedDict(), []
    for s, files in samples.items():
        try:
            fq = prepare_fastq(s, files)
            prep[s] = phase_a(s, fq)
        except Exception as e:  # keep going with the other samples
            failed.append(s)
            info("[%s] FAILED in step 1: %s" % (s, e))
    if prep:
        run("multiqc -f -q -o %s -n masoudi_multiqc_qc_trimming %s %s %s" % (
            q(DIRS["mqc1"]), q(DIRS["qc_raw"]), q(DIRS["trim"]), q(DIRS["qc_trim"])), "multiqc_step1")
        info("MultiQC (step 1): %s" % os.path.join(DIRS["mqc1"], "masoudi_multiqc_qc_trimming.html"))

    # ---------- STEPS 2-4 ----------
    rows = []
    for s, (trimmed, st_raw, st_trim) in prep.items():
        info("=== STEPS 2-4 for sample %s ===" % s)
        try:
            if st_trim["reads"] == 0:
                raise RuntimeError("0 reads left after trimming - relax --trim-quality/--min-length")
            rows.append(phase_b(s, trimmed, st_raw, st_trim, ref, fai))
            info("[%s] done" % s)
        except Exception as e:
            failed.append(s)
            info("[%s] FAILED: %s" % (s, e))
            traceback.print_exc()

    # ---------- summary ----------
    if rows:
        tsv = os.path.join(DIRS["summary"], "qc_summary.tsv")
        with open(tsv, "w") as fh:
            fh.write("\t".join(rows[0].keys()) + "\n")
            for r in rows:
                fh.write("\t".join(str(v) for v in r.values()) + "\n")
        plot_overview(rows, A.min_depth, os.path.join(DIRS["summary"], "run_overview.png"))
        info("QC summary: %s" % tsv)
        try:
            run("multiqc -f -q -o %s -n masoudi_multiqc_full_report %s" % (
                q(DIRS["summary"]),
                " ".join(q(DIRS[k]) for k in ("qc_raw", "trim", "qc_trim", "map", "dedup", "var"))), "multiqc_final")
        except Exception as e:
            info("final MultiQC failed: %s" % e)
        if A.no_sam:
            for f in glob.glob(os.path.join(DIRS["map"], "*.sam")):
                os.remove(f)
    info("=== FINISHED: %d ok, %d failed %s ===" % (len(rows), len(failed), failed if failed else ""))
    print("\nResults in: %s\n  consensus : 08_consensus/\n  VCF       : 07_variants/*.final.vcf  (+ *.lofreq.filtered.vcf, *.bcftools.filtered.vcf)"
          "\n  BAM/BAI   : 06_dedup/*.final.bam(.bai)\n  plots     : 10_plots/\n  QC summary: 11_summary/qc_summary.tsv + MultiQC html\n" % out)
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()