__author__ = "Patrik Smeds"
__copyright__ = "Copyright 2023, Patrik Smeds"
__email__ = "patrik.smeds@gmail.com"
__license__ = "MIT"

import glob
import sys
from snakemake.shell import shell
from tempfile import TemporaryDirectory
from tempfile import mkdtemp

log = snakemake.log_fmt_shell(stdout=False, stderr=True)
extra = snakemake.params.get("extra", "")

#with TemporaryDirectory(delete=not args.keep_temp_files) as tmpdirname:
#with TemporaryDirectory(delete=False) as tmpdirname:
# with TemporaryDirectory() as tmpdirname:
tmpdirname = mkdtemp() 
    try:
        shell(
            """
            (paraphase --bam {snakemake.input.bam} \
            --reference {snakemake.input.fasta} \
            --out {tmpdirname} \
            {snakemake.params.genome} \
            {extra}) {log}
            """
        )
    except Exception as e:
        print(f"Error running paraphase: {e}")
        sys.exit(1)

    # Concatenating, reheadering, and sorting the zipped and indexed VCF files
    vcf_res = glob.glob(f"{tmpdirname}/**/*variants.vcf")
    
    if vcf_res:
        for vcf in vcf_res:
            bgzip_cmd = f"bgzip -c {vcf} > {vcf}.gz"
            shell(bgzip_cmd)
            index_cmd = f"bcftools index {vcf}.gz"
            shell(index_cmd)
            print(f"Compressed and indexed: {vcf}.gz")
        params_variant_files = " ".join([f"{vcf}.gz" for vcf in vcf_res]) 
        shell(
            f"bcftools concat -a -Oz {params_variant_files} | "
            f"bcftools annotate --header-lines {snakemake.input.vcf_header} | "
            f"bcftools sort -Oz -o {snakemake.output.merged_vcf} {log}" 
        )
        print(f"Merged, reheadered, and sorted VCF file created: {snakemake.output.merged_vcf}")
    else:
        print("No output VCF files were produced by paraphase, I hope this is what you were expecting, human?")
        shell(f"touch {snakemake.output.merged_vcf}")

 
