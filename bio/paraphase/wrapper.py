__author__ = "Patrik Smeds"
__copyright__ = "Copyright 2023, Patrik Smeds"
__email__ = "patrik.smeds@gmail.com"
__license__ = "MIT"

import glob
import sys
from snakemake.shell import shell
from tempfile import mkdtemp
import shutil

log = snakemake.log_fmt_shell(stdout=False, stderr=True)
extra = snakemake.params.get("extra", "")

# Using mkdtemp() to create a temporary directory that is not automatically deleted
tmpdirname = mkdtemp()
print("TMPDIR: ",tmpdirname)

try:
    shell(
        """
        (paraphase --bam {snakemake.input.bam} \
        --reference {snakemake.input.fasta} \
        --out {tmpdirname} \
        {snakemake.params.genome} \
        {extra}) {log}
        """,
        tmpdirname=tmpdirname  # Pass tmpdirname to the shell environment
    )

    # Concatenating, reheadering, and sorting the zipped and indexed VCF files
    vcf_res = glob.glob(f"/{tmpdirname}/**/*variants.vcf", recursive=True)
    print("RES: ", vcf_res )
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

except Exception as e:
    print(f"Error running paraphase: {e}")
    sys.exit(1)

finally:
    # Uncomment the line below if you decide you want to clean up the directory after debugging
    # shutil.rmtree(tmpdirname)
    pass 

 
