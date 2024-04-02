__author__ = "Patrik Smeds"
__copyright__ = "Copyright 2023, Patrik Smeds"
__email__ = "patrik.smeds@gmail.com"
__license__ = "MIT"

import glob
import sys
from snakemake.shell import shell
from tempfile import TemporaryDirectory


log = snakemake.log_fmt_shell(stdout=False, stderr=True)
extra = snakemake.params.get("extra", "")

# shell("touch {snakemake.output.merged_vcf} ")


with TemporaryDirectory() as tmpdirname:
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

    vcf_res = glob.glob(f"{tmpdirname}/**/*variants.vcf")
    vcf_res2 = glob.glob(f"{snakemake.params.outfolder}/**/*variants.vcf")
    print("VCFRES", tmpdirname, vcf_res, vcf_res2)
    
    if vcf_res:
        for vcf in vcf_res:
            bgzip_cmd = f"bgzip -c {vcf} > {vcf}.gz"
            shell(bgzip_cmd)
            index_cmd = f"bcftools index {vcf}.gz"
            shell(index_cmd)
            print(f"Compressed and indexed: {vcf}.gz")
        # Concatenating, reheadering, and sorting the zipped and indexed VCF files
        #concatenated_file = f"{tmpdirname}/concatenated.vcf.gz"
        params_variant_files = " ".join([f"{vcf}.gz" for vcf in vcf_res])  # Adjust if your file naming needs differ
        shell(
            f"bcftools concat -a -Oz {params_variant_files} | "
            f"bcftools annotate --header-lines {snakemake.input.vcf_header} | "
            f"bcftools sort -Oz -o {snakemake.output.merged_vcf} {log}" 
        )
        print(f"Merged, reheadered, and sorted VCF file created: {snakemake.output.merged_vcf}")
    else:
        print("No output VCF files were produced by paraphase, I hope this is what you were expecting, human?")
        shell(f"touch {snakemake.output.merged_vcf}")

 
#shell("find {snakemake.output.outres}/*_variants.vcf -type f -exec bgzip -f {{}} \\; ")
#shell("find {snakemake.output.outres}/*_variants.vcf.gz -type f -name '*_variants.vcf.gz' -exec bcftools index {{}} \\;")
#shell("bcftools concat -a -O v {params.variant_files}")


