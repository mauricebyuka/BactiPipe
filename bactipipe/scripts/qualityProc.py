#!/usr/bin/env python3
import subprocess
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from bactipipe.scripts.utils import time_print, simple_print, logger

import os, sys, shutil, bz2, lzma
from tqdm import tqdm

# Faster gzip if present
try:
    from isal import igzip as gzmod
except Exception:
    import gzip as gzmod

def fastp_version():
    vercmd = ["fastp", "--version"]
    try:
        process = subprocess.Popen(vercmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = process.communicate()
        fastp_ver = out.decode('utf-8', errors='replace').strip().split()[-1]
        fastp_ver = "v" + fastp_ver
    except Exception as e:
        fastp_ver = "Unknown"
    return fastp_ver

def _iter_fastq_with_progress(fastq_path, desc):
    total_bytes = os.path.getsize(fastq_path)
    raw = open(fastq_path, 'rb', buffering=1024 * 1024)

    if fastq_path.endswith('.gz'):
        fh = gzmod.GzipFile(fileobj=raw, mode='rb'); tell = raw.tell
    elif fastq_path.endswith('.bz2'):
        fh = bz2.BZ2File(raw); tell = raw.tell
    elif fastq_path.endswith(('.xz', '.lzma')):
        fh = lzma.LZMAFile(raw); tell = raw.tell
    else:
        fh = raw; tell = fh.tell

    cols = shutil.get_terminal_size(fallback=(100, 20)).columns
    # account for your desc length (expand tabs conservatively)
    desc_display = desc.expandtabs(2)
    # reserve ~40 chars for counters/ETA; cap bar width between 10 and 30
    bar_cols = max(10, min(30, cols - 40 - len(desc_display)))
    bar_format = (
        f"{{desc}} |{{bar:{bar_cols}}}| {{percentage:3.0f}}% "
        f"({{n_fmt}}/{{total_fmt}}) [{{elapsed}}<{{remaining}}]"
    )

    try:
        with tqdm(
            total=total_bytes,
            file=sys.stdout,
            ncols=cols,
            bar_format=bar_format,
            ascii=True,
            position=0,
            leave=False,
            mininterval=0.20,
            smoothing=0.1,
            disable=not sys.stdout.isatty()
        ) as pbar:
            prev = tell()
            r = fh.readline
            while True:
                title = r()
                if not title:
                    break
                seq = r(); plus = r(); qual = r()
                if not title.startswith(b'@') or not plus.startswith(b'+') or not seq.strip() or len(seq.rstrip(b'\r\n')) != len(qual.rstrip(b'\r\n')):
                    raise ValueError('Malformed or truncated FASTQ record')

                cur = tell(); delta = cur - prev
                if delta > 0:
                    pbar.update(delta); prev = cur

                yield qual
    finally:
        try: fh.close()
        finally:
            if fh is not raw:
                raw.close()

class ProcQuality:
    def __init__(self, fastq_files, genome_size=None, plot_path='output_plots', sample=None, logfile=None):
        """ Initialize the quality processing object """
        self.sample = sample
        self.fastq_files = fastq_files
        self.genome_size = genome_size
        self.plot_path = plot_path
        self.file_quality_counts = []  # Store quality counts for each file separately
        self.total_counts = []  # Store total counts of bases for each file separately
        self.combined_total_bases = 0
        self.combined_num_reads = 0
        self.combined_quality_sum = 0
        self.individual_stats = []  # Store individual stats for each file
        self.log = logfile

        # Process the FASTQ files and compute quality statistics
        self._process_fastq_files()

    def _process_fastq_files(self):
        """Process each fastq file to collect quality statistics separately
        Write important information to a log file
        """
        log = self.log
        if self.sample:
            time_print(f'QC assessment for sample : {self.sample}', "Header")
            logger(log, f'QC assessment for sample : {self.sample}', "Header")

        # Updated info messages (no pyfastx / no index)
        info = ('Reading quality data from FASTQ by streaming (xopen or gzip), '
                f'NumPy version: {np.__version__}')
        logger(log, info)
        time_print(info)

        outinfo = ('Loading quality and base count data into NumPy for fast vectorized operations '
                f'(NumPy {np.__version__})')
        logger(log, outinfo)
        time_print(outinfo)

        try:
            for i, fastq_file in enumerate(self.fastq_files):
                # Per-file accumulators — keep original variable names
                file_total_bases = 0
                file_num_reads = 0
                file_quality_sum = 0

                # Positional quality tallies (grow as needed)
                quality_counts_above_30 = np.zeros(0, dtype=np.int64)
                quality_counts_20_to_30 = np.zeros(0, dtype=np.int64)
                quality_counts_below_20 = np.zeros(0, dtype=np.int64)
                total_counts = np.zeros(0, dtype=np.int64)

                desc = f"\t---> Processing file {i+1}/{len(self.fastq_files)}"
                for qual in _iter_fastq_with_progress(fastq_file, desc=desc):
                    # Vectorized Phred+33 decode
                    # Strip newlines but keep as bytes for frombuffer
                    qv = np.frombuffer(qual.rstrip(b"\r\n"), dtype=np.uint8).astype(np.int32) - 33
                    if np.any(qv < 0) or np.any(qv > 93):
                        raise ValueError('Invalid Phred+33 quality character')
                    L = qv.size

                    # Grow positional arrays if this read is longer than what we've seen
                    if L > quality_counts_above_30.size:
                        pad = L - quality_counts_above_30.size
                        quality_counts_above_30 = np.pad(quality_counts_above_30, (0, pad))
                        quality_counts_20_to_30 = np.pad(quality_counts_20_to_30, (0, pad))
                        quality_counts_below_20 = np.pad(quality_counts_below_20, (0, pad))
                        total_counts = np.pad(total_counts, (0, pad))

                    # Update per-position categories (use original variable names)
                    quality_counts_above_30[:L] += (qv >= 30)
                    quality_counts_20_to_30[:L] += ((qv >= 20) & (qv < 30))
                    quality_counts_below_20[:L] += (qv < 20)
                    total_counts[:L] += 1

                    # Update per-file totals
                    file_total_bases += L
                    file_num_reads += 1
                    file_quality_sum += int(qv.sum())

                # Store quality counts and total counts for this file (original structure)
                self.file_quality_counts.append({
                    'above_30': quality_counts_above_30,
                    'between_20_and_30': quality_counts_20_to_30,
                    'below_20': quality_counts_below_20
                })
                self.total_counts.append(total_counts)

                # Accumulate combined statistics (original variables)
                self.combined_total_bases += file_total_bases
                self.combined_num_reads += file_num_reads
                self.combined_quality_sum += file_quality_sum

                # Save individual file statistics (original keys)
                avg_quality = file_quality_sum / file_total_bases if file_total_bases > 0 else 0
                perc_above_30 = (np.sum(quality_counts_above_30) / file_total_bases * 100) if file_total_bases > 0 else 0
                perc_20_to_30 = (np.sum(quality_counts_20_to_30) / file_total_bases * 100) if file_total_bases > 0 else 0
                perc_below_20 = (np.sum(quality_counts_below_20) / file_total_bases * 100) if file_total_bases > 0 else 0
                coverage = file_total_bases / self.genome_size if self.genome_size else "N/A"

                self.individual_stats.append({
                    'num_reads': file_num_reads,
                    'total_bases': file_total_bases,
                    'perc_above_30': perc_above_30,
                    'perc_20_to_30': perc_20_to_30,
                    'perc_below_20': perc_below_20,
                    'avg_quality': avg_quality,
                    'coverage': coverage
                })

        except Exception as e:
            error_info = f'Error processing FASTQ files: {e}'
            logger(log, error_info)
            time_print(error_info, "Fail")
            raise


    def plot_quality_distribution(self):
        """ Generate stacked bar plots for the proportion of bases in each quality category at each position for each file """
        log = self.log
        outinfo = f'Generating a quality distribution plot using the python module "matplotlib" (version: {matplotlib.__version__})'
        logger(log, outinfo)
        time_print(outinfo)

        q_plot = f"{self.plot_path}_quality_distribution.png"
        if not os.path.exists(q_plot):
            num_files = len(self.fastq_files)
            fig, axes = plt.subplots(1, num_files, figsize=(10 * num_files, 6))

            # If only one file, make axes iterable
            if num_files == 1:
                axes = [axes]

            for i, ax in enumerate(axes):
                quality_counts = self.file_quality_counts[i]
                total_counts = self.total_counts[i]

                # Calculate proportions
                proportion_below_20 = quality_counts['below_20'] / total_counts
                proportion_20_to_30 = quality_counts['between_20_and_30'] / total_counts
                proportion_above_30 = quality_counts['above_30'] / total_counts

                positions = np.arange(1, len(proportion_below_20) + 1)

                # Plot stacked bars
                ax.bar(positions, proportion_below_20, label='Q<20', color='red')
                ax.bar(positions, proportion_20_to_30, bottom=proportion_below_20, label='Q20-29', color='orange')
                ax.bar(positions, proportion_above_30, bottom=proportion_below_20 + proportion_20_to_30, label='Q30+', color='green')

                # Set custom x-axis labels at 10 intervals
                max_position = positions[-1]
                interval_positions = np.linspace(1, max_position, 10, dtype=int)
                ax.set_xticks(interval_positions)
                ax.set_xticklabels([str(pos) for pos in interval_positions])

                # Labels and title for each subplot
                ax.set_xlabel('Position in Read')
                ax.set_ylabel('Proportion of Bases')
                if num_files == 2:
                    ax.set_title(f"Pair {i+1} Base Quality Distribution", fontsize=14)
                else:
                    ax.set_title("Base Quality Distribution", fontsize=14)

            # Add a legend to the first plot only
            axes[0].legend()

            # Save the combined plot figure
            plt.tight_layout()
            plt.savefig(q_plot, dpi=300, bbox_inches='tight')
            plt.close()
        else:
            message1 = f'The quality distribution plot exists already. Skipping ploting!\n'

            time_print(message1)
        outinfo = f'Quality distribution plot: {q_plot}\n'
        time_print(outinfo, "Pass")
        logger(log, outinfo, "Pass")


    def plot_quality_metrics(self):
        log = self.log
        #Plot and save metrics for individual files and combined data
        outinfo = 'Generating a plot of quality scores using the python module "matplotlib"'
        logger(log, outinfo)
        time_print(outinfo)

        metrics_plot = f"{self.plot_path}_metrics.png"
        metrics_txt = f"{self.plot_path}_metrics.txt"
        if not os.path.exists(metrics_plot):
            num_files = len(self.fastq_files)
            fig2, metrics_ax = plt.subplots(figsize=(7, 2) if num_files == 2 else (1.5, 1))
            metrics_ax.axis('off')  # Hide axis
            with open(metrics_txt, 'a') as metrics_file:

                for i, stats in enumerate(self.individual_stats):
                    if num_files == 1 and self.genome_size:
                        stats_text = (f"Tot reads: {stats['num_reads']:,}\n"
                                    f"Tot bases: {stats['total_bases']:,}\n"
                                    f"Q30+: {stats['perc_above_30']:.2f}%\n"
                                    f"Q20-29: {stats['perc_20_to_30']:.2f}%\n"
                                    f"Q<20: {stats['perc_below_20']:.2f}%\n"
                                    f"Avg Quality: {stats['avg_quality']:.2f}\n"
                                    f"Genome cov: {stats['coverage']:.2f}X")
                        metrics_file.write(stats_text + "\n")
                    else:
                        stats_text = (f"Pair {i+1}:\n\n"
                                    f"Reads: {stats['num_reads']:,}\n"
                                    f"Bases: {stats['total_bases']:,}\n"
                                    f"Q30+: {stats['perc_above_30']:.2f}%\n"
                                    f"Q20-29: {stats['perc_20_to_30']:.2f}%\n"
                                    f"Q<20: {stats['perc_below_20']:.2f}%\n"
                                    f"Mean Quality: {stats['avg_quality']:.2f}\n")
                        metrics_file.write(stats_text + "\n")

                    if i == 0:
                        metrics_ax.text(0, 0.5, stats_text, fontsize=12, va='center', ha='left', 
                                        bbox=dict(facecolor='white', edgecolor='black', boxstyle="round,pad=0.3"))
                    else:
                        metrics_ax.text(0.70, 0.5, stats_text, fontsize=12, va='center', ha='left', 
                                        bbox=dict(facecolor='white', edgecolor='black', boxstyle="round,pad=0.3"))

                if num_files == 2:
                    combined_avg_quality = self.combined_quality_sum / self.combined_total_bases if self.combined_total_bases > 0 else 0
                    combined_coverage = self.combined_total_bases / self.genome_size if self.genome_size else "N/A"
                    comb_cov = f"{combined_coverage:.2f}X" if isinstance(combined_coverage, float) else combined_coverage

                    print(f"Pair 1 + Pair 2:\n\nAvg Quality: {combined_avg_quality}\nTotal reads: {self.combined_num_reads}\nTotal bases: {self.combined_total_bases}\nGenome cov: {comb_cov}")

                    combined_stats_text = (f"Pair 1 + Pair 2:\n\n"
                                        f"Avg Quality: {combined_avg_quality:.2f}\n"
                                        f"Total reads: {self.combined_num_reads:,}\n"
                                        f"Total bases: {self.combined_total_bases:,}\n"
                                        f"Genome cov: {comb_cov}")

                    metrics_file.write(combined_stats_text + "\n")

                    metrics_ax.text(0.31, 0.5, combined_stats_text, fontsize=12, va='center', ha='left', 
                                    bbox=dict(facecolor='white', edgecolor='black', boxstyle="round,pad=0.3"))

            plt.tight_layout()
            plt.savefig(metrics_plot, dpi=300, bbox_inches='tight')
        else:
            message2 = f'The metrics plot exists already. Skipping the ploting!\n'
            time_print(message2)
        outinfo = f'Quality metrics data: {metrics_plot}\n'
        time_print(outinfo, "Pass")
        logger(log, outinfo, "Pass")

    # Properties for retrieving statistics
    @property
    def reads(self):
        """Retrieve the number of reads for each file."""
        return [stat['num_reads'] for stat in self.individual_stats]

    @property
    def bases(self):
        """Retrieve the total number of bases for each file."""
        return [stat['total_bases'] for stat in self.individual_stats]

    @property
    def q30(self):
        """Retrieve the percentage of Q30+ bases for each file."""
        return [stat['perc_above_30'] for stat in self.individual_stats]

    @property
    def q20(self):
        """Retrieve the percentage of Q20-29 bases for each file."""
        return [stat['perc_20_to_30'] for stat in self.individual_stats]

    @property
    def q0(self):
        """Retrieve the percentage of Q<20 bases for each file."""
        return [stat['perc_below_20'] for stat in self.individual_stats]

    @property
    def tot_reads(self):
        """Retrieve the combined number of reads across all files."""
        return self.combined_num_reads

    @property
    def tot_bases(self):
        """Retrieve the combined number of bases across all files."""
        return self.combined_total_bases

    @property
    def coverage(self):
        """Calculate the combined genome coverage, if genome size is provided."""
        if self.genome_size:
            return self.combined_total_bases / self.genome_size
        else:
            return "N/A"

    @property
    def q_plot(self):
        """Retrieve the path to the quality distribution plot."""
        return f"{self.plot_path}_quality_distribution.png"

    @property
    def metrics_plot(self):
        """Retrieve the path to the metrics plot."""
        return f"{self.plot_path}_metrics.png"

    @property
    def file_average_quality(self):
        """Retrieve the average quality score for each individual file."""
        return [stat['avg_quality'] for stat in self.individual_stats]

    @property
    def average_quality(self):
        """Retrieve the combined average quality score across files."""
        if self.combined_total_bases > 0:
            return self.combined_quality_sum / self.combined_total_bases
        return 0
