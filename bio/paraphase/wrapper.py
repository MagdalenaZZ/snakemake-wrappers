import glob
import sys
from snakemake.shell import shell
from tempfile import TemporaryDirectory
import shutil

log = snakemake.log_fmt_shell(stdout=True, stderr=True)
extra = snakemake.params.get("extra", "")

try:
    with TemporaryDirectory() as tmpdirname:
        shell(
            f"""
            (paraphase --bam {snakemake.input.bam} \
            --reference {snakemake.input.fasta} \
            --out {tmpdirname} \
            {snakemake.params.genome} \
            {extra}) {log}
            """,
            tmpdirname=tmpdirname  # Pass tmpdirname to the shell environment
        )

        # Create a new VCF header
        shell(
            f"""
            python -c "import sys; fai_file=open('{snakemake.input.faidx}', 'r'); lines = fai_file.readlines(); fai_file.close(); output = open('{snakemake.output.vcf_header}', 'w'); [output.write(f'##contig=<ID={{line.split()[0]}},length={{line.split()[1]}}>\\n') for line in lines]; output.close()"
            """
        )

        # Concatenating, reheadering, and sorting the zipped and indexed VCF files
        vcf_res = glob.glob(f"{tmpdirname}/*_vcfs/*vcf")
        print("VCFRES: ", vcf_res)
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
                f"bcftools annotate --header-lines {snakemake.output.vcf_header} | "
                f"bcftools sort -Oz -o {snakemake.output.merged_vcf}"
            )
            print(f"Merged, reheadered, and sorted VCF file created: {snakemake.output.merged_vcf}")
            # Copy out bam and bai files
            bam_res = glob.glob(f"{tmpdirname}/*.bam")
            bai_res = glob.glob(f"{tmpdirname}/*.bai")
            print("BAM RES: ", bam_res, bai_res)
            shell(
                f"""
                cp -pr {' '.join(bam_res)} {snakemake.output.bam};
                cp -pr {' '.join(bai_res)} {snakemake.output.bai}
                """
            )
        else:
            print("No output VCF files were produced by paraphase, I hope this is what you were expecting, human?")
            shell(f"touch {snakemake.output.merged_vcf}")

except Exception as e:
    print(f"Error running paraphase: {e}")
    sys.exit(1)


