# GC-MS
python plot_dep2_results.py \
    --input /data/local/jy1008/SaMu/results/latest/proteomics_GC-MS/GC_MS_dep2_results.csv \
    --contrast Sarc_vs_NoSarc \
    --experiment-name GC-MS \
    --abundance-table /data/local/jy1008/SaMu/results/latest/proteomics_GC-MS/GC_MS_dep2_vsn_matrix.csv \
    --sample-metadata /data/local/jy1008/SaMu/results/latest/proteomics_GC-MS/GC_MS_dep2_vsn_imputed_metadata.csv \
    --output /data/local/jy1008/SaMu/results/latest/proteomics_GC-MS/GC_MS_dep2_volcano_top_features

# proteomics
python plot_dep2_results.py \
    --input /data/local/jy1008/SaMu/results/latest/proteomics/proteomics_dep2_results.csv \
    --contrast Sarc_vs_NoSarc \
    --experiment-name proteomics \
    --abundance-table /data/local/jy1008/SaMu/results/latest/proteomics/proteomics_dep2_vsn_matrix.csv \
    --sample-metadata /data/local/jy1008/SaMu/results/latest/proteomics/proteomics_dep2_vsn_imputed_metadata.csv \
    --output /data/local/jy1008/SaMu/results/latest/proteomics/proteomics_dep2_volcano_top_features

# NMR
python plot_nmr_quorum_results.py \
    --raw /data/local/jy1008/SaMu/results/latest/nmr/samu_nmr_FullSaMu.csv \
    --scaled /data/local/jy1008/SaMu/results/latest/nmr/samu_nmr_logrel_scaled.csv \
    --out /data/local/jy1008/SaMu/results/latest/nmr/nmr_top_features \
    --top-n 10

# QUORUM
python plot_nmr_quorum_results.py \
    --raw /data/local/jy1008/SaMu/results/latest/quorum/samu_quorum_raw_qsp.csv \
    --scaled /data/local/jy1008/SaMu/results/latest/quorum/samu_quorum_log_scaled_qsp.csv \
    --pseudocount 1e-6 \
    --value-label "Quant*Prob" \
    --feature-label "QSP" \
    --out /data/local/jy1008/SaMu/results/latest/quorum/quorum_qsp_fig \
    --stats-out /data/local/jy1008/SaMu/results/latest/quorum/quorum_qsp_stats.csv

python plot_nmr_quorum_results.py \
    --raw /data/local/jy1008/SaMu/results/latest/quorum/samu_quorum_raw_species.csv \
    --scaled /data/local/jy1008/SaMu/results/latest/quorum/samu_quorum_log_scaled_species.csv \
    --pseudocount 1e-6 \
    --value-label "Quant*Prob" \
    --feature-label "Species" \
    --out /data/local/jy1008/SaMu/results/latest/quorum/quorum_species_fig \
    --stats-out /data/local/jy1008/SaMu/results/latest/quorum/quorum_species_stats.csv

python plot_nmr_quorum_results.py \
    --raw /data/local/jy1008/SaMu/results/latest/quorum/samu_quorum_raw_microbial_target.csv \
    --scaled /data/local/jy1008/SaMu/results/latest/quorum/samu_quorum_log_scaled_microbial_target.csv \
    --pseudocount 1e-6 \
    --value-label "Quant*Prob" \
    --feature-label "Microbial target" \
    --out /data/local/jy1008/SaMu/results/latest/quorum/quorum_microbial_target_fig \
    --stats-out /data/local/jy1008/SaMu/results/latest/quorum/quorum_microbial_target_stats.csv