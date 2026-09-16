from bactipipe.scripts.utils import time_print, simple_print, logger
from Bio import SeqIO
import subprocess
import os
import csv
import shlex

from bactipipe.scripts.pdf_renderer import render_pdf_from_tsv, LayoutConfig
from reportlab.lib.units import inch
from bactipipe.scripts.qc_plan import flye_size_arguments, overall_quality_verdict

def trimmer(sample, raw_reads, newReadsFolder, cpus=24, logfile=None):
    # with open (self.log, 'a') as log:
    message = '\n====== Trimming raw reads for quality control ======\n'
    time_print(message, "Header")
    if not os.path.exists(newReadsFolder):
        os.makedirs(newReadsFolder)

    print(f'\n\tTrimming reads for sample : {sample}\n')
    
    if len(raw_reads) == 2:
        fastq1 = raw_reads[0]
        fastq2 = raw_reads[1]
        newRead1 = os.path.join(newReadsFolder, sample, f'{sample}_R1.fastq.gz')
        newRead2 = os.path.join(newReadsFolder, sample, f'{sample}_R2.fastq.gz')
        cmd = ["fastp", "-i", fastq1, "-I", fastq2, "-o", newRead1, "-O", newRead2, "-w", str(cpus), "-j", "/dev/null", "-h", "/dev/null", "-q", "30", "-u", "10", "-l", "100", "-z", "4"]
        execute = None
        if not os.path.exists(newRead1) and not os.path.exists(newRead2):
            execute = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        outmessage = [
        '\n\tDone Trimming:\n',
        f'\t---> Trimmed read1: {newRead1}',
        f'\t---> Trimmed read2: {newRead2}\n'
        ]
      
    elif len(raw_reads) == 1:
        single_reads = raw_reads[0]
        newRead = os.path.join(newReadsFolder, sample, f'{sample}.fastq.gz')
        cmd = ["fastp", "-i", single_reads, "-o", newRead, "-w", str(cpus), "-j", "/dev/null", "-h", "/dev/null", "-q", "30", "-u", "10", "-l", "100", "-z", "4"]
        execute = None
        if not os.path.exists(newRead):
            execute = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        outmessage = [
        '\n\tDone Trimming:\n',
        f'\t---> Trimmed reads: {newRead}',
        ]

    if execute is not None and execute.returncode != 0:
        time_print("Failed to trimming reads for quality.", "Fail") # Check the log file for more details.")
        return
    else:
        message = '\n\t Trimmed reads exist already. Skiping QC!!!'
        # log.write(message +'\n')
        print(message)

    for message in outmessage:
        # log.write(message + '\n')
        time_print(message, "Pass")
# return newReadsFolder

def filter_genome(input_fasta, output_fasta, min_length=500):
    records = []
    for record in SeqIO.parse(input_fasta, 'fasta'):
        if len(record.seq) >= min_length:
            records.append(record)
    SeqIO.write(records, output_fasta, 'fasta')
    return output_fasta

def assemble(sample, reads, assembly_dir, assembler, sequencer, cpus=24, logfile=None, gsize=None, single=True, memory_gb=None):
    log = os.path.join(assembly_dir, sample, sample + ".log")
    try:
        if single:
            time_print(f'Assembling the genome for sample : {sample}', "Header")
        logger(log, f'Assembling the genome for sample : {sample}', "Header")
        
        assembler = assembler.strip().lower()
        if assembler == "unicycler":
            tool = "unicycler"
        elif assembler == "spades":
            tool = "spades.py"
        elif assembler == "flye":
            tool = "flye"
        elif assembler == "skesa":
            tool = "skesa"
        else:
            return_info = f"Assembler {assembler} is not supported."
            if single:
                time_print(return_info, "Fail")
            logger(log, return_info)
            return

        # Get tool version
        vercmd = [tool, "--version"]
        stdout = subprocess.run(vercmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out = (stdout.stdout or stdout.stderr).decode('utf-8', errors='replace')
        version = out.split()[-1] if out.split() else "Unknown"

        if single:
            simple_print(f'\t---> Assembler: {assembler} (version: {version})')
            simple_print(f'\t---> Number of fastq files: {len(reads)}')
            simple_print(f'\t---> Number of CPUs: {cpus}\n')
            
        logger(log, f'\t---> Assembler: {assembler}', mode="simple")
        logger(log, f'\t---> Number of fastq files: {len(reads)}', mode="simple")
        logger(log, f'\t---> Number of CPUs: {cpus}\n', mode="simple")

        if len(reads) == 2:
            read1 = reads[0]
            read2 = reads[1]
        elif len(reads) == 1:
            single_reads = reads[0]

        # assembly_dir = os.path.join(outDir, 'assemblies')   
        tempDir = os.path.join(assembly_dir, sample)
        genomesDir = os.path.join(assembly_dir, 'genomes')
        os.makedirs(genomesDir, exist_ok=True)

        finalAssembly = os.path.join(genomesDir, sample + '.fasta')

        if sequencer.lower() == "illumina":
            if assembler == "unicycler":
                if len(reads) == 2:
                    cmd1 = ["unicycler", "-1", read1, "-2", read2, "-o", tempDir, "--threads", str(cpus)]
                elif len(reads) == 1:
                    cmd1 = ["unicycler", "-s", single_reads, "-o", tempDir, "--threads", str(cpus)]
                draft_assembly = f'{tempDir}/assembly.fasta'
            elif assembler == "spades":
                if len(reads) == 2:
                    cmd1 = ["spades.py", "--isolate", "-1", read1, "-2", read2, "-o", tempDir, "--threads", str(cpus)]
                elif len(reads) == 1:
                    cmd1 = ["spades.py", "--isolate", "-s", single_reads, "-o", tempDir, "--threads", str(cpus)]
                if memory_gb is not None:
                    cmd1 += ['--memory', str(int(memory_gb))]
                draft_assembly = f'{tempDir}/scaffolds.fasta'
            elif assembler == "skesa":
                if len(reads) == 2:
                    cmd1 = ["skesa", "--reads", f"{read1},{read2}", "--contigs_out", f"{tempDir}/assembly.fasta", "--cores", str(cpus)]
                elif len(reads) == 1:
                    cmd1 = ["skesa", "--reads", single_reads, "--contigs_out", f"{tempDir}/assembly.fasta", "--cores", str(cpus)]
                if memory_gb is not None:
                    cmd1 += ['--memory', str(int(memory_gb))]
                draft_assembly = f'{tempDir}/assembly.fasta'
            else:
                return_info = f"Assembler {assembler} is not supported for Illumina reads. Please use Skesa, Spades or Unicycler."
                if single:
                    time_print(return_info, "Fail")
                logger(log, return_info)
                return
        
        elif sequencer.lower() == "nanopore":
            if assembler == "unicycler":
                cmd1 = ["unicycler", "-l", single_reads, "-o", tempDir, "--threads", str(cpus)]
                draft_assembly = f'{tempDir}/assembly.fasta'
            elif assembler == "flye":
                cmd1 = [
                    "flye",
                    "--nano-hq",
                    single_reads,
                    "--out-dir",
                    tempDir,
                    "--threads",
                    str(cpus),
                    *flye_size_arguments(gsize),
                    "--iterations",
                    "2",
                ]
                draft_assembly = f'{tempDir}/assembly.fasta'
            else:
                return_info = f"Assembler {assembler} is not supported for Nanopore reads. Please use Unicycler or Flye."
                if single:
                    time_print(return_info, "Fail")
                logger(log, return_info)
                return
        
        if not os.path.exists(finalAssembly):
            # message1 = f'\t---> Assembling the reads for sample {sample}'
            # print(message1)
            info2 = f"Running: {shlex.join(cmd1)}"
            if single:
                time_print(info2, s_type="command")
            logger(log, info2, s_type="command")
            # Sample QC workers deliberately cap numerical-library thread pools.
            # Do not let their OMP_NUM_THREADS=1 override SPAdes --threads (or
            # assemblers launched by Unicycler). Change only the child environment;
            # concurrent workers and Python/BLAS retain their existing limits.
            assembly_env = dict(os.environ, OMP_NUM_THREADS=str(max(1, int(cpus))))
            logger(log, f'Assembler OpenMP thread budget: {assembly_env["OMP_NUM_THREADS"]}', mode='simple')
            execute = subprocess.run(cmd1, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                     env=assembly_env)

            if execute.stdout:
                logger(log, execute.stdout.decode('utf-8', errors='replace'), mode='simple')

            if execute.returncode != 0:
                return_info = f"Assembly failed for {sample}: exit {execute.returncode}"
                time_print(return_info, "Fail")
                logger(log, return_info) #" Check the log file for more details.")
                return
            else:
                return_info = f"Command exit status: Success!"
                if single:
                    time_print(return_info, "Pass")
                logger(log, return_info)
            
            # Filter the genome to remove smaller contigs
            filter_genome(input_fasta=draft_assembly, output_fasta=finalAssembly, min_length=400)
        else:
            message2 = f'The assembly file for the sample {sample} exists already. Skipping the assembly!\n'
            # log.write(message2 + '\n')
            if single:
                time_print(message2)
            logger(log, message2)

        outinfo = f'Done with Assembly :: Assembly file: {finalAssembly}\n'
        if single:
            time_print(outinfo, "Pass")
        logger(log, outinfo)
    except Exception as e:
        time_print(f"Unexpected error in assemble for {sample}: {e}", "Fail")
        logger(log, f"Unexpected error in assemble for {sample}: {e}")
        return
                            
def checkM_stats(genomes_dir, outdir, cpus=24, logfile=None, exclude_samples=()):
    genomes_dir, outdir = os.fspath(genomes_dir), os.fspath(outdir)
    log = logfile
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    # Get checkM version
    vercmd = ["checkm", "-h"]
    stdout = subprocess.run(vercmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if stdout.returncode != 0:
        message = f"CheckM is not installed. Please install CheckM to proceed"
        logger(log, message)
        return {}, "Unknown"
    else:
        out = stdout.stdout.decode('utf-8')
        version = "Unknown"
        for line in out.split('\n'):
            if ":::" in line:
                version = line.split(":::")[1].strip().split()[-1]
                break
            
    if exclude_samples:
        from pathlib import Path
        from tempfile import TemporaryDirectory
        excluded = set(exclude_samples)
        # Read-only links leave final assemblies available for reporting/retries.
        with TemporaryDirectory(prefix="bactipipe-checkm-inputs-") as tmp:
            eligible = [path for path in Path(genomes_dir).glob("*.fasta") if path.stem not in excluded]
            logger(log, f"CheckM excludes Lambda controls: {', '.join(sorted(excluded))}")
            if not eligible:
                logger(log, "CheckM skipped: no bacterial assemblies to assess.")
                return {}, version
            for path in eligible:
                (Path(tmp) / path.name).symlink_to(path.resolve())
            return checkM_stats(tmp, outdir, cpus=cpus, logfile=log)

    checkm_intro = f'Assessing Genome Quality with CheckM (version: {version})'
    time_print(checkm_intro, "Header")
    logger(log, checkm_intro, "Header")
    
    checkM_out = os.path.join(outdir, 'checkM_stats.txt')
    pplacer_threads = max(1, int(int(cpus) / 4))
    cmd = ["checkm", "lineage_wf", "-t", str(cpus), "--pplacer_threads", str(pplacer_threads), "-q", "--tab_table", "-x", "fasta", genomes_dir, outdir]
    
    checkM_done = False
    if os.path.exists(checkM_out):
        with open(checkM_out, 'r') as file:
            if any(line.startswith("Bin Id\t") for line in file):
                message = f'The CheckM output file exists already. Skipping CheckM!!!\n'
                logger(log, message)
                time_print(message)
                checkM_done = True
    if not checkM_done:
        info3 = f"Running: {shlex.join(cmd)} > {checkM_out}"
        time_print(info3, s_type='command')
        logger(log, info3, s_type='command')
        with open(checkM_out, "w") as out_fh:
            execute = subprocess.run(cmd, stdout=out_fh, stderr=subprocess.STDOUT)
        if execute.returncode != 0:
            return_info = f"CheckM failed. For more detaisls, check the CheckM log file in the directory '{checkM_out}"
            time_print(return_info, "Fail")
            logger(log, return_info) #" Check the log file for more details.")
            return {}, version
        else:
            return_info = f"Command exit status: Success!"
            time_print(return_info, "Pass")
            logger(log, return_info)

    checkM_output = {}
    with open(checkM_out, 'r') as file:
        for line in file:
            fields = line.rstrip("\n").split('\t')
            if len(fields) < 3 or fields[0] == "Bin Id":
                continue
            try:
                float(fields[-3])
                float(fields[-2])
            except (TypeError, ValueError):
                continue
            checkM_output[fields[0]] = [fields[-3], fields[-2]] # Completeness and contamination

    return checkM_output, version

def make_summary(
    qc_summary,
    temp_qc_summary,
    header,
    checkm_out,
    logfile=None,
    min_completeness=90,
    max_contamination=10,
    defer_pdf=False,
):
    from bactipipe.lambda_control import is_lambda_control, CONTROL_IDENTIFICATION
    log = logfile
    with(open(qc_summary , 'w')) as qc_sum, open(temp_qc_summary, 'r') as temp_sum:
        writer = csv.writer(qc_sum, dialect='excel-tab')
        for line in header:
            writer.writerow([line])
        writer.writerow(["Quality Summary"])
        writer.writerow([""])
        writer.writerow(["Sample ID",  "Mean quality", "Expected organism", "Identified organism", "% Match", "Coverage Depth", "Required depth", "CheckM completeness", "CheckM contamination", "Overall Quality", "Isolate QC", "Control QC"])
        lines = list(temp_sum)
        # Stable ordering preserves the sample sheet order within each group.
        lines.sort(key=lambda line: not is_lambda_control(line.rstrip("\n").split("\t")[3]))
        for line in lines:
            sample, avqc, qc_verdict, organism, hit, match, coverage, min_depth, cov_verdict, tax_confirm = line.strip().split("\t")
            if sample == "Sample":
                continue
            control = is_lambda_control(organism)
            control_qc = tax_confirm.removeprefix("Control") if control and tax_confirm in {"ControlPass", "ControlFail", "ControlError"} else "Not assessed"
            if control and control_qc == "Not assessed":
                if not hit.lower().startswith("analysis failed"):
                    hit = CONTROL_IDENTIFICATION
                match = "N/A"
            if sample in checkm_out and not control:
                completeness, contamination = checkm_out[sample]
                from bactipipe.qc_status import threshold_verdict
                completeness_verdict = threshold_verdict(completeness, min_completeness)
                contamination_verdict = threshold_verdict(contamination, max_contamination, maximum=True)
            else:
                completeness, contamination = "N/A", "N/A"
                completeness_verdict = "N/A"
                contamination_verdict = "N/A"

            from bactipipe.qc_status import interpret
            interpretation = interpret(organism, hit, [qc_verdict, cov_verdict,
                                      completeness_verdict, contamination_verdict],
                                      isolate_qc_failure=tax_confirm == "Contamination", control_qc=control_qc)
            final_verdict = interpretation['overall']
            # Explicit taxonomy review flags may include independently detected
            # isolate concerns; do not turn those into a silent pass.
            if tax_confirm in {"Fail", "Review"} and final_verdict == "Pass":
                final_verdict = "Review"
            coverage = coverage.replace("X", "")
            cov_display = f'{float(coverage):.2f}X' if coverage != 'N/A' else 'N/A'
            match_display = (
                "N/A"
                if str(match).strip().upper() in {"", "N/A", "NA", "NONE", "NULL", "-"}
                else f"{float(match):.2f}"
            )

            writer.writerow([sample, f'{float(avqc):.2f}', organism, hit, match_display, cov_display, f"{min_depth}X", completeness, contamination, final_verdict, "Fail" if tax_confirm == "Contamination" else "Pass", control_qc if control else "N/A"])

            simple_print(f'{sample} :: QC {final_verdict}!\n', final_verdict)
            logger(log, f'{sample} :: QC {final_verdict}!', final_verdict, mode="simple")

            simple_print(f'---> Species confirmation: {tax_confirm} --- {hit}', tax_confirm)
            logger(log, f'---> Species confirmation: {tax_confirm} --- {hit}', tax_confirm, mode="simple")
            simple_print(f'---> Mean quality score: {qc_verdict} --- {float(avqc):.2f}', qc_verdict)
            logger(log, f'---> Mean quality score: {qc_verdict} --- {float(avqc):.2f}', qc_verdict, mode="simple")
            if coverage != 'N/A':
                simple_print(f'---> Genome coverage: {cov_verdict} --- {float(coverage):.2f}X', cov_verdict)
                logger(log, f'---> Genome coverage: {cov_verdict} --- {float(coverage):.2f}X', cov_verdict, mode="simple")
            else:
                simple_print(f'---> Genome coverage: N/A --- {coverage}')
                logger(log, f'---> Genome coverage: N/A --- {coverage}', mode="simple")
            
            simple_print(f'---> CheckM completeness: {completeness}', completeness_verdict)
            logger(log, f'---> CheckM completeness: {completeness}', completeness_verdict, mode="simple")
            simple_print(f'---> CheckM contamination: {contamination}\n', contamination_verdict)
            logger(log, f'---> CheckM contamination: {contamination}\n', contamination_verdict, mode="simple")
    # Write the pdf summary
    run_name = ""
    if os.path.exists(qc_summary) and not defer_pdf:
        with open(qc_summary, "r") as tsv:
            for line in tsv:
                if "Run name:" in line:
                    run_name = line.strip().split(":")[1] 
                    break
        
        title = f"Run {run_name}: QC Metrics"
        layout = LayoutConfig(
            reserved_cols=(2, 3),          
            reserved_width=1.7*inch,       
            last3_narrow_width=0.8*inch,   
            # left_align_cols=(0, 2, 3),
            # left_panel_ratio=0.65, 
            # right_panel_ratio=0.35
        )
        pdf_file = qc_summary.replace('.tsv', '.pdf')
        render_pdf_from_tsv(tsv_path=qc_summary, pdf_path=pdf_file, title=title, layout=layout, run_name=run_name)
        if os.path.exists(pdf_file):
            simple_print(f'PDF summary file generated: {pdf_file}', "Pass")
            logger(log, f'PDF summary file generated: {pdf_file}', mode="simple")
