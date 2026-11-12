# Heatmap of viral_rpm.tsv: vOTUs in rows, samples in columns, reads per million in five bins.
#
#   Rscript plot.R results/viral_rpm.tsv [-s samples.txt] [-g column] [-o prefix]
#
# -s  text file with the samples to draw, one name per line, drawn in that order (default: all)
# -g  column of the table (for example condition) that colours and groups the samples; samples
#     with the same value next to each other form a group
# -o  output prefix, writes <prefix>.pdf and <prefix>.png (default: next to the table)
#
# vOTUs without reads in the drawn samples are left out. Rows are ordered by their highest RPM.

suppressPackageStartupMessages(library(ComplexHeatmap))
ht_opt(message = FALSE)
# without this ComplexHeatmap warns that it cannot find Arial
if (!"Arial" %in% names(pdfFonts())) pdfFonts(Arial = pdfFonts()$Helvetica)

args <- commandArgs(trailingOnly = TRUE)
flag <- function(f) if (f %in% args) args[match(f, args) + 1] else NULL
if (!length(args) || startsWith(args[1], "-")) stop("usage: Rscript plot.R viral_rpm.tsv [-s samples.txt] [-g column] [-o prefix]")
infile <- args[1]
out <- if (is.null(flag("-o"))) sub("\\.[^.]*$", "_heatmap", infile) else flag("-o")

BREAKS <- c(1, 10, 100, 1e3, 1e4)
BIN_COLS <- c("0" = "#E6E6E6", "1" = "#CAC5F3", "2" = "#A79FE1", "3" = "#8175CB", "4" = "#594CA0", "5" = "#312271")
BIN_LABELS <- c("< 1", "1–10", "10–100", "100–1,000", "1,000–10,000", "≥ 10,000")
GROUP_COLS <- c("#0072B2", "#D55E00", "#009E73", "#E69F00", "#56B4E9", "#CC79A7", "#F0E442")   # Okabe-Ito
CELL <- unit(3.2, "mm")
GAP <- unit(1.5, "mm")
gp <- function(...) gpar(fontsize = 6, fontfamily = "Arial", ...)

tab <- read.delim(infile, check.names = FALSE, colClasses = c(sample = "character"))
first <- grep("\\]$", names(tab))[1]
if (is.na(first)) stop("no vOTU column in ", infile)
if (!is.null(flag("-s"))) {
  samples <- trimws(readLines(flag("-s")))
  samples <- unique(samples[samples != "" & !startsWith(samples, "#")])
  missing <- setdiff(samples, tab$sample)
  if (length(missing)) stop("not in the table: ", paste(missing, collapse = ", "))
} else {
  samples <- tab$sample
}
if (!length(samples)) stop("no samples to draw")
drawn <- tab[match(samples, tab$sample), ]

votu <- names(tab)[first:ncol(tab)]
label <- sub(" \\[[^]]*\\]$", "", votu)
if (anyDuplicated(label)) label <- votu
rpm <- t(as.matrix(drawn[votu]))
top <- apply(rpm, 1, max)
keep <- which(top > 0)
if (!length(keep)) stop("no vOTU has reads in these samples")
keep <- keep[order(-top[keep])]
rpm <- rpm[keep, , drop = FALSE]
label <- label[keep]
bin <- matrix(as.character(findInterval(rpm, BREAKS)), nrow(rpm))
n_row <- nrow(bin)
n_col <- ncol(bin)

# samples with the same value next to each other form a group
group <- flag("-g")
if (!is.null(group)) {
  if (!group %in% names(tab)[seq_len(first - 1)])
    stop("no column ", group, " in the table; columns: ", paste(names(tab)[2:(first - 1)], collapse = ", "))
  value <- as.character(drawn[[group]])
  run <- cumsum(c(TRUE, value[-1] != value[-n_col]))
  every <- unique(as.character(tab[[group]]))          # colours from the whole table, so they do not change with -s
  colour <- setNames(GROUP_COLS[(seq_along(every) - 1) %% length(GROUP_COLS) + 1], every)
}

build <- function() {
  a <- list(bin, name = "rpm", col = BIN_COLS, show_heatmap_legend = FALSE,
            rect_gp = gpar(col = "white", lwd = 0.8),
            cluster_rows = FALSE, cluster_columns = FALSE, height = n_row * CELL,
            row_labels = label, row_names_side = "left", row_names_gp = gp(),
            row_names_max_width = max_text_width(label, gp = gp()),
            column_labels = samples, column_names_side = "bottom", column_names_gp = gp(),
            column_names_max_height = max_text_width(samples, gp = gp()))
  if (is.null(group)) {
    a$width <- n_col * CELL
  } else {
    title <- value[!duplicated(run)]
    # turn the titles when one is wider than its group
    wide <- convertWidth(max_text_width(title, gp = gp()), "mm", valueOnly = TRUE) > min(table(run)) * 3.2
    a <- c(a, list(width = n_col * CELL + (max(run) - 1) * GAP,
                   column_split = factor(run), column_gap = GAP,
                   column_title = title, column_title_gp = gp(), column_title_rot = if (wide) 90 else 0,
                   top_annotation = HeatmapAnnotation(group = value, col = list(group = colour),
                                                      simple_anno_size = unit(1.2, "mm"),
                                                      show_annotation_name = FALSE, show_legend = FALSE,
                                                      gp = gpar(col = NA))))
  }
  do.call(Heatmap, a)
}
legend <- function() Legend(title = "Reads per million", at = names(BIN_COLS), labels = BIN_LABELS,
                            legend_gp = gpar(fill = BIN_COLS), grid_width = CELL, grid_height = CELL,
                            border = "white", labels_gp = gp(), title_gp = gp(fontface = "bold"))
paint <- function() {
  lg <- legend()
  d <- draw(build(), heatmap_legend_list = list(lg), heatmap_legend_side = "right", padding = unit(rep(2, 4), "mm"))
  c(w = convertWidth(ComplexHeatmap:::width(d), "inch", valueOnly = TRUE),
    h = convertHeight(ComplexHeatmap:::height(d), "inch", valueOnly = TRUE),
    legend = convertHeight(ComplexHeatmap:::height(lg), "inch", valueOnly = TRUE))
}

# the first drawing is only to get the size of the figure
cairo_pdf(tempfile(fileext = ".pdf"), width = 20, height = 20, family = "Arial")
size <- paint()
invisible(dev.off())
w <- size[["w"]]
h <- max(size[["h"]], size[["legend"]] + 8 / 25.4)       # the legend can be higher than the heatmap
cairo_pdf(paste0(out, ".pdf"), width = w, height = h, family = "Arial")
invisible(paint())
invisible(dev.off())
png(paste0(out, ".png"), width = w, height = h, units = "in", res = 300, type = "cairo", family = "Arial")
invisible(paint())
invisible(dev.off())

message(n_row, " vOTUs x ", n_col, " samples -> ", out, ".pdf/.png (", round(w * 25.4), " mm wide)")
