comb_file="/PHShome_actual/j/jy1008/.local/share/mamba/envs/humann39_env/lib/python3.11/site-packages/humann/tests/data/demo.fastq"
METAPHLAN_V31_BOWTIE2_DB="/data/bwh-comppath-full/databases/metaphlan_v31" 
METAPHLAN_V31_INDEX="mpa_v31_CHOCOPhlAn_201901"
metaphlan_out_file="metaphlan31_demo_test_profile_v31.txt"

micromamba activate metaphlan31_env
metaphlan \
    "$comb_file" \
    --input_type fastq \
    --bowtie2db "$METAPHLAN_V31_BOWTIE2_DB" \
    --index "$METAPHLAN_V31_INDEX" \
    --nproc 8 \
    -t rel_ab_w_read_stats \
    -o "$metaphlan_out_file"

HUMANN_CHOCOPHLAN_DB="/data/bwh-comppath-full/databases/humann/chocophlan/chocophlan"
HUMANN_UNIREF_DB="/data/bwh-comppath-full/databases/humann/uniref/uniref"
sample_out_dir="humann39_demo_out"

micromamba activate humann39_env
humann \
    --input "$comb_file" \
    --output "$sample_out_dir" \
    --taxonomic-profile "$metaphlan_out_file" \
    --nucleotide-database "$HUMANN_CHOCOPHLAN_DB" \
    --protein-database "$HUMANN_UNIREF_DB" \
    --threads 8






# ===============
# Run Humann 3.9 with Metaphlan 4.1.1
# ===============
# find $(python -c "import humann; print(humann.__path__[0])") -iname "demo*"
# to grab data path
#
micromamba activate humann39_env

demo_data="/PHShome_actual/j/jy1008/.local/share/mamba/envs/humann39_env/lib/python3.11/site-packages/humann/tests/data/demo.fastq"
HUMANN_CHOCOPHLAN_DB="/data/bwh-comppath-full/databases/humann/chocophlan/chocophlan"
HUMANN_UNIREF_DB="/data/bwh-comppath-full/databases/humann/uniref/uniref"
humann --input "$demo_data" \
    --output demo_test_vJun23 \
    --threads 4 \
    --nucleotide-database "$HUMANN_CHOCOPHLAN_DB" \
    --protein-database "$HUMANN_UNIREF_DB" \
    --metaphlan-options "--bowtie2db /data/bwh-comppath-seq/databases/bowtie2 --index mpa_vJun23_CHOCOPhlAnSGB_202403"
