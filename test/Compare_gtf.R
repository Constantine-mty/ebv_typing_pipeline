


# 我需要比较 B95-8/type 1(NC_007605.1)与 AG876/type 2(NC_009334.1)在关键基因上注释的差异

library(rtracklayer)
gtf1 <- rtracklayer::import("~/Project/0X.EBV_Genome/biomni_pipeline/ebv_typing_pipeline/resources/ebv1_genes.gtf")
gtf2 <- rtracklayer::import("~/Project/0X.EBV_Genome/biomni_pipeline/ebv_typing_pipeline/resources/ebv2_genes.gtf")

gtf1
# GRanges object with 94 ranges and 6 metadata columns:
#   seqnames        ranges strand |   source     type     score     phase              gene_id            gene_name
# <Rle>     <IRanges>  <Rle> | <factor> <factor> <numeric> <integer>          <character>          <character>
#   [1] NC_007605.1     1691-5856      + |   RefSeq     gene        NA      <NA> tegument protein G75 tegument protein G75
# [2] NC_007605.1    9631-10262      + |   RefSeq     gene        NA      <NA> interleukin-10 BCRF1 interleukin-10 BCRF1
# [3] NC_007605.1   11305-97654      + |   RefSeq     gene        NA      <NA>                EBNA1                EBNA1
# [4] NC_007605.1   11305-89482      + |   RefSeq     gene        NA      <NA>               EBNA3B               EBNA3B
# [5] NC_007605.1   11305-82962      + |   RefSeq     gene        NA      <NA>               EBNA3A               EBNA3A
# ...         ...           ...    ... .      ...      ...       ...       ...                  ...                  ...
# [90] NC_007605.1 166011-177679      + |   RefSeq     gene        NA      <NA>                LMP2A                LMP2A
# [91] NC_007605.1 166483-169088      - |   RefSeq     gene        NA      <NA>                 LMP1                 LMP1
# [92] NC_007605.1 166483-167067      - |   RefSeq     gene        NA      <NA>       protein BNLF2a       protein BNLF2a
# [93] NC_007605.1 166483-166836      - |   RefSeq     gene        NA      <NA>       protein BNLF2b       protein BNLF2b
# [94] NC_007605.1 169294-177679      + |   RefSeq     gene        NA      <NA>                LMP2B                LMP2B


# 我提取坐标信息，并且有什么好的可视化方法可以展示这几个基因在两个坐标下的差异

library(rtracklayer)
library(GenomicRanges)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(forcats)
library(readr)


focus.gene <- c(
  "EBNA1", "EBNA2", "EBNA3A", "EBNA3B", "EBNA3C",
  "LMP1", "LMP2A", "LMP2B",
  "BZLF1", "BZLF2", "BLLF1"
)

## =========================
## 1. GRanges -> data.frame
## =========================

gr_to_df <- function(gr, strain, type_label, genome_length = NA_integer_) {
  as.data.frame(gr) %>%
    mutate(
      strain = strain,
      EBV_type = type_label,
      genome_length = genome_length,
      gene_id = as.character(gene_id),
      gene_name = as.character(gene_name),
      feature_type = as.character(type),
      seqnames = as.character(seqnames),
      strand = as.character(strand)
    ) %>%
    select(
      strain, EBV_type, genome_length,
      seqnames, start, end, width, strand,
      source, feature_type, score, phase,
      gene_id, gene_name,
      everything()
    )
}

gtf1_df <- gr_to_df(
  gtf1,
  strain = "B95-8",
  type_label = "type 1",
  genome_length = 171823
)

gtf2_df <- gr_to_df(
  gtf2,
  strain = "AG876",
  type_label = "type 2",
  genome_length = 172764
)

gtf_all_df <- bind_rows(gtf1_df, gtf2_df)


table(gtf_all_df$strain, gtf_all_df$feature_type)

## =========================
## 2. focus gene 提取
## =========================

focus_gene_df <- gtf_all_df %>%
  filter(feature_type == "gene") %>%
  filter(gene_name %in% focus.gene | gene_id %in% focus.gene) %>%
  mutate(
    gene = case_when(
      gene_name %in% focus.gene ~ gene_name,
      gene_id %in% focus.gene ~ gene_id,
      TRUE ~ gene_name
    ),
    midpoint = (start + end) / 2,
    crosses_origin = !is.na(genome_length) & end > genome_length,
    start_mod = if_else(
      !is.na(genome_length),
      ((start - 1) %% genome_length) + 1,
      start
    ),
    end_mod = if_else(
      !is.na(genome_length),
      ((end - 1) %% genome_length) + 1,
      end
    )
  ) %>%
  arrange(gene, strain)

focus_gene_df



## =========================
## 3. alias 辅助匹配
## =========================

# alias_tbl <- tibble::tribble(
#   ~target_gene, ~pattern,
#   "EBNA1",  "EBNA1|BKRF1",
#   "EBNA2",  "EBNA2|BYRF1",
#   "EBNA3A", "EBNA3A|BERF1|BERF2",
#   "EBNA3B", "EBNA3B|BERF3|BERF4",
#   "EBNA3C", "EBNA3C|BERF5|BERF6",
#   "LMP1",   "LMP1|BNLF1",
#   "LMP2A",  "LMP2A",
#   "LMP2B",  "LMP2B",
#   "BZLF1",  "BZLF1|Zta|ZEBRA",
#   "BZLF2",  "BZLF2|gp42|glycoprotein gp42",
#   "BLLF1",  "BLLF1|gp350|gp220|major envelope glycoprotein"
# )
# 
# extract_focus_by_alias <- function(df, alias_tbl) {
#   map_df(seq_len(nrow(alias_tbl)), function(i) {
#     tg <- alias_tbl$target_gene[i]
#     pat <- alias_tbl$pattern[i]
#     
#     df %>%
#       filter(feature_type == "gene") %>%
#       filter(
#         str_detect(gene_name, regex(pat, ignore_case = TRUE)) |
#           str_detect(gene_id, regex(pat, ignore_case = TRUE))
#       ) %>%
#       mutate(gene = tg)
#   }) %>%
#     distinct(strain, EBV_type, seqnames, start, end, strand, gene_id, gene_name, gene, .keep_all = TRUE) %>%
#     mutate(
#       midpoint = (start + end) / 2,
#       crosses_origin = !is.na(genome_length) & end > genome_length,
#       start_mod = if_else(
#         !is.na(genome_length),
#         ((start - 1) %% genome_length) + 1,
#         start
#       ),
#       end_mod = if_else(
#         !is.na(genome_length),
#         ((end - 1) %% genome_length) + 1,
#         end
#       )
#     ) %>%
#     arrange(gene, strain, start)
# }
# 
# focus_gene_df <- extract_focus_by_alias(gtf_all_df, alias_tbl)
# 
# focus_gene_df %>%
#   select(
#     strain, EBV_type, gene, gene_id, gene_name,
#     seqnames, start, end, width, strand,
#     crosses_origin, start_mod, end_mod
#   )


focus_check <- focus_gene_df %>%
  dplyr::count(gene, strain) %>%
  complete(gene = focus.gene, strain = c("B95-8", "AG876"), fill = list(n = 0)) %>%
  arrange(gene, strain)

focus_check


## =========================
## 4. 宽表比较
## =========================

focus_gene_genelevel_df <- focus_gene_df %>%
  filter(feature_type == "gene") %>%
  group_by(strain, EBV_type, gene) %>%
  summarise(
    gene_id = paste(unique(gene_id), collapse = ";"),
    gene_name = paste(unique(gene_name), collapse = ";"),
    seqnames = paste(unique(seqnames), collapse = ";"),
    start = min(start, na.rm = TRUE),
    end = max(end, na.rm = TRUE),
    width = end - start + 1,
    strand = paste(unique(strand), collapse = ";"),
    genome_length = unique(genome_length)[1],
    crosses_origin = any(crosses_origin, na.rm = TRUE),
    start_mod = ((start - 1) %% genome_length) + 1,
    end_mod = ((end - 1) %% genome_length) + 1,
    .groups = "drop"
  ) %>%
  mutate(
    strain_key = case_when(
      strain == "B95-8" ~ "type1_B95_8",
      strain == "AG876" ~ "type2_AG876",
      TRUE ~ strain
    )
  )


focus_gene_genelevel_df %>%
  dplyr::count(gene, strain) %>%
  arrange(gene, strain)

coord_compare_df <- focus_gene_genelevel_df %>%
  select(
    gene, strain_key,
    gene_id, gene_name, seqnames,
    start, end, width, strand,
    crosses_origin, start_mod, end_mod
  ) %>%
  pivot_wider(
    id_cols = gene,
    names_from = strain_key,
    values_from = c(
      gene_id, gene_name, seqnames,
      start, end, width, strand,
      crosses_origin, start_mod, end_mod
    ),
    names_glue = "{.value}_{strain_key}"
  ) %>%
  mutate(
    across(
      c(
        start_type1_B95_8, end_type1_B95_8, width_type1_B95_8,
        start_mod_type1_B95_8, end_mod_type1_B95_8,
        start_type2_AG876, end_type2_AG876, width_type2_AG876,
        start_mod_type2_AG876, end_mod_type2_AG876
      ),
      as.numeric
    ),
    delta_start_AG876_minus_B95_8 = start_type2_AG876 - start_type1_B95_8,
    delta_end_AG876_minus_B95_8 = end_type2_AG876 - end_type1_B95_8,
    delta_width_AG876_minus_B95_8 = width_type2_AG876 - width_type1_B95_8,
    width_ratio_AG876_over_B95_8 = width_type2_AG876 / width_type1_B95_8
  ) %>%
  arrange(factor(gene, levels = focus.gene))

coord_compare_df


coord_compare_df %>%
  select(
    gene,
    start_type1_B95_8, end_type1_B95_8, width_type1_B95_8,
    start_type2_AG876, end_type2_AG876, width_type2_AG876,
    delta_start_AG876_minus_B95_8,
    delta_end_AG876_minus_B95_8,
    delta_width_AG876_minus_B95_8,
    width_ratio_AG876_over_B95_8
  )


## =========================
## Plot 1: genome-level gene location
## =========================

plot_df <- focus_gene_df %>%
  mutate(
    gene = factor(gene, levels = rev(focus.gene)),
    strain_label = paste0(strain, " / ", EBV_type),
    y = as.numeric(gene)
  )

p_gene_location <- ggplot(plot_df) +
  geom_segment(
    aes(
      x = start,
      xend = end,
      y = gene,
      yend = gene,
      color = gene
    ),
    linewidth = 5,
    lineend = "round"
  ) +
  geom_point(
    aes(x = start, y = gene),
    size = 1.8
  ) +
  geom_point(
    aes(x = end, y = gene),
    size = 1.8
  ) +
  facet_wrap(
    ~ strain_label,
    ncol = 1,
    scales = "free_x"
  ) +
  scale_x_continuous(
    labels = scales::comma,
    expand = expansion(mult = c(0.02, 0.04))
  ) +
  labs(
    x = "Genome coordinate",
    y = NULL,
    title = "Key EBV gene annotation positions in B95-8/type 1 and AG876/type 2",
    color = "Gene"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  )

p_gene_location

ggsave(
  "EBV_type1_type2_focus_gene_location.pdf",
  p_gene_location,
  width = 9,
  height = 6
)

ggsave(
  "EBV_type1_type2_focus_gene_location.png",
  p_gene_location,
  width = 9,
  height = 6,
  dpi = 300
)





