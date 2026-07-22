# =============================================================================
# matrix_utils.R
#
# Shared helpers for building sample x feature matrices and group
# annotations from a merged metadata + assay data.frame. Used by both
# 01a_nmr.R and 01b_quorum.R so this logic doesn't drift between them.
#
# Source after config.R and after dplyr/tibble are loaded (both callers do
# this already). Functions take their inputs as arguments rather than
# reading globals, so they're safe to source anywhere.
# =============================================================================

# ---------------------------------------------------------------------------
# build_matrix — coerce selected columns to numeric, drop all-NA rows,
# and return a matrix with record_id as row names.
#
# Args:
#   df   : data.frame containing record_id and the feature columns
#   cols : character vector of feature column names to include
# ---------------------------------------------------------------------------
build_matrix <- function(df, cols) {
  m <- df %>%
    mutate(across(all_of(cols), as.numeric)) %>%  # NAs-by-coercion expected
    select(all_of(cols)) %>%
    as.matrix()
  rownames(m) <- df$record_id
  m[rowSums(is.na(m)) < ncol(m), , drop = FALSE]
}

# ---------------------------------------------------------------------------
# group_annotation — ordered group annotation + the row gap between groups,
# aligned to a matrix's row order.
#
# Args:
#   mat                : sample x feature matrix (rownames = record_id)
#   df                 : data.frame containing record_id and group_var
#   group_var          : name of the grouping column (e.g. "sarc_status_bin")
#   group_levels       : c(reference, other) factor level order
#   cluster_within     : if TRUE, hierarchically cluster rows separately
#                        within each group level (using mat's values), then
#                        concatenate the blocks in group_levels order. If
#                        FALSE (default), rows within a group keep mat's
#                        original order. Either way the group split itself
#                        (and thus `gap`) is preserved, so this stays
#                        compatible with pheatmap's gaps_row + cluster_rows
#                        = FALSE — gaps_row requires no clustering be done
#                        by pheatmap itself, so this is how you get rows
#                        grouped by similarity within each split without
#                        losing the visual boundary between groups. The
#                        trade-off: no row dendrogram is drawn, since that
#                        requires handing pheatmap an actual clustering
#                        object for cluster_rows, which gaps_row can't
#                        accept.
#   clustering_method  : passed to hclust() when cluster_within = TRUE
#
# Returns list(anno, order, gap):
#   anno  : data.frame of group_var, ordered to match mat rows after grouping
#   order : integer index (into mat's original rows) to reorder mat to
#           match anno
#   gap   : row count of the first group level (for pheatmap's gaps_row)
# ---------------------------------------------------------------------------
group_annotation <- function(mat, df, group_var, group_levels,
                             cluster_within = FALSE,
                             clustering_method = "ward.D2") {
  anno <- df %>%
    filter(record_id %in% rownames(mat)) %>%
    arrange(match(record_id, rownames(mat))) %>%
    select(record_id, all_of(group_var)) %>%
    column_to_rownames("record_id")
  anno[[group_var]] <- factor(anno[[group_var]], levels = group_levels)

  if (cluster_within) {
    # Cluster each group's rows on complete-case columns only (NA columns
    # would otherwise make dist()/hclust() fail); this affects the row
    # ORDER only — the plotted matrix itself keeps every column.
    ordered_ids <- unlist(lapply(group_levels, function(lvl) {
      ids <- rownames(anno)[anno[[group_var]] == lvl]
      if (length(ids) < 2) return(ids)  # nothing to cluster
      sub_mat <- mat[ids, , drop = FALSE]
      sub_mat <- sub_mat[, colSums(is.na(sub_mat)) == 0, drop = FALSE]
      if (ncol(sub_mat) < 2) return(ids)  # not enough complete features
      hc <- hclust(dist(sub_mat), method = clustering_method)
      ids[hc$order]
    }))
    anno <- anno[ordered_ids, , drop = FALSE]
    ord  <- match(ordered_ids, rownames(mat))
  } else {
    ord  <- order(anno[[group_var]])
    anno <- anno[ord, , drop = FALSE]
  }

  list(anno = anno, order = ord,
       gap = sum(anno[[group_var]] == group_levels[1]))
}