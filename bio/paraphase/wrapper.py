__author__ = "Patrik Smeds"
__copyright__ = "Copyright 2023, Patrik Smeds"
__email__ = "patrik.smeds@gmail.com"
__license__ = "MIT"

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=False, stderr=True)
extra = snakemake.params.get("extra", "")


#annotate = snakemake.input.get("annotate", "")
#if annotate:
#    annotate = f"--annotate {annotate}"

shell(
    """
    (TMPDIR={tmp_dir}; \
    head  \
    {snakemake.input.bam} \
    >  {snakemake.output.vcf} \
    ) {log}
    """
)




#with tempfile.TemporaryDirectory(dir=tmp_root) as tmp_dir:
#    shell(
#        """
#        (TMPDIR={tmp_dir}; \
#        pbmm2 align --num-threads {snakemake.threads} \
#            --preset {snakemake.params.preset} \
#            --sample {snakemake.params.sample} \
#            --log-level {snakemake.params.loglevel} \
#            {extra} \
#            {snakemake.input.reference} \
#            {snakemake.input.query} \
#            {snakemake.output.bam}) {log}
#        """
#    )
