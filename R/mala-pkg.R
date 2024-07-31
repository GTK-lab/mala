
options::define_option(
             "unibind-fpath",
             desc="resource path use for Unibind",
             default="https://unibind.uio.no/static/data/20220914/bulk_Robust/Homo_sapiens/damo_hg38_TFBS_per_TF.tar.gz",
             option_name="MalaUnibindSource",
             envvar_name="mala_unibind_source")

options::define_option(
             "unibind-rname",
             desc="resource name to use for Unibind",
             default="hg38-robust-damo-20220914",
             option_name="MalaUnibindSource",
             envvar_name="mala_unibind_source")
