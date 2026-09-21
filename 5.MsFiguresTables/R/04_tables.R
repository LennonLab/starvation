################################################################################
# 04_tables.R
#
# The manuscript's LaTeX tables, and an Overleaf upload staged beside them.
#
# The tables are not retyped. Each analysis project writes a markdown table
# with its caption; this converts that file, so a number in the manuscript
# cannot differ from the number the analysis produced. Re-run after re-running
# any project, then re-upload.
#
# The manuscript text itself is not in this repository -- it is drafted in
# Overleaf, where it changes constantly and a git history of it would be noise.
# What lives here is everything generated: the figures and these tables.
#
# Output: output/tables/*.tex   one file per table
#         overleaf/             figures/ and tables/, laid out the way
#                               manuscript.tex expects them. Not tracked;
#                               rebuild and re-upload when anything changes.
################################################################################

.setup <- Filter(file.exists, c("R/00_setup.R", "00_setup.R", "../R/00_setup.R"))
if (!length(.setup)) stop("Run this from the manuscript-figures project root (or R/).")
source(.setup[1])

TAB_DIR <- file.path(OUT_DIR, "tables")
dir.create(TAB_DIR, showWarnings = FALSE, recursive = TRUE)

## ---- markdown -> LaTeX -------------------------------------------------------

#' Escape the characters LaTeX treats specially. Backslash first, or the
#' escapes we add are themselves escaped.
tex_escape <- function(x) {
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  for (ch in c("&", "%", "$", "#", "_", "{", "}")) {
    x <- gsub(ch, paste0("\\", ch), x, fixed = TRUE)
  }
  x <- gsub("~", "\\textasciitilde{}", x, fixed = TRUE)
  x <- gsub("^", "\\textasciicircum{}", x, fixed = TRUE)
  x
}

# Symbols the projects write as UTF-8 and LaTeX wants as commands. Converting
# them keeps every .tex file plain ASCII, so the engine and the input encoding
# stop mattering -- which is what makes the Overleaf upload predictable.
SYMBOLS <- c(
  "σ" = "$\\sigma$",  "µ" = "$\\mu$",     "μ" = "$\\mu$",
  "θ" = "$\\theta$",  "τ" = "$\\tau$",    "β" = "$\\beta$",
  "Δ" = "$\\Delta$",  "±" = "$\\pm$",     "×" = "$\\times$",
  "≤" = "$\\le$",     "≥" = "$\\ge$",     "²" = "$^{2}$",
  "—" = "---",        "–" = "--",         "−" = "$-$",
  "‘" = "`",          "’" = "'",
  "“" = "``",         "”" = "''")

apply_symbols <- function(x) {
  for (k in names(SYMBOLS)) x <- gsub(k, SYMBOLS[[k]], x, fixed = TRUE)
  x
}

#' 3.10e+08 -> $3.10\times10^{8}$, and 10^6^ -> $10^{6}$.
#'
#' A leading < or > is pulled inside the maths. Left in text mode LaTeX sets
#' those characters as something else entirely.
sci_notation <- function(x) {
  x <- gsub("([0-9]+(?:\\.[0-9]+)?)[eE]\\+?(-?)0*([0-9]+)",
            "$\\1\\\\times10^{\\2\\3}$", x)
  x <- gsub("10\\^(-?[0-9]+)\\^", "$10^{\\1}$", x)
  # First: a sign butting against maths the rule above already made.
  x <- gsub("([<>])\\s*\\$", "$\\1", x)
  # Then: a sign against a plain number. The [^$] guard keeps this off the
  # ones just handled, whose sign now sits immediately after a $.
  gsub("(^|[^$])([<>])\\s*([0-9])", "\\1$\\2\\3$", x)
}

#' Markdown emphasis. Bold before italic, or ** is eaten as two * markers.
markup <- function(x) {
  x <- gsub("`([^`]+)`", "\\\\texttt{\\1}", x)
  x <- gsub("\\*\\*([^*]+)\\*\\*", "\\\\textbf{\\1}", x)
  gsub("\\*([^*]+)\\*", "\\\\textit{\\1}", x)
}

to_tex <- function(x) markup(sci_notation(apply_symbols(tex_escape(x))))

#' Column alignment from the markdown separator row.
md_align <- function(sep_cells) {
  vapply(sep_cells, function(s) {
    s <- trimws(s)
    if (grepl("^:.*:$", s)) "c" else if (grepl(":$", s)) "r" else "l"
  }, character(1), USE.NAMES = FALSE)
}

split_row <- function(line) {
  line <- sub("^\\s*\\|", "", sub("\\|\\s*$", "", trimws(line)))
  trimws(strsplit(line, "\\s*\\|\\s*")[[1]])
}

#' Convert one of the projects' table files.
#'
#' They are laid out as: "## Table N. Title", a blank line, one or more caption
#' paragraphs, then the pipe table. The title becomes the caption and the
#' paragraphs become the table note.
convert_table <- function(path, label, fullwidth = FALSE, fontsize = NULL,
                          italic_col = NULL, landscape = FALSE) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")

  title <- sub("^##\\s*", "", lines[grepl("^##\\s", lines)][1])
  title <- sub("^Table\\s*[0-9A-Za-z]+\\.\\s*", "", title)

  is_row  <- grepl("^\\s*\\|", lines)
  tbl     <- lines[is_row]
  caption <- lines[!is_row & !grepl("^#", lines) & nzchar(trimws(lines))]

  if (length(tbl) < 3) stop("No table found in ", path)

  header <- split_row(tbl[1])
  align  <- md_align(split_row(tbl[2]))
  body   <- lapply(tbl[-(1:2)], split_row)

  ncol <- length(header)
  body <- Filter(function(r) length(r) == ncol, body)

  # Gene names belong in italics in a manuscript. The projects' own tables set
  # them plain, because a markdown reader of those files does not need it.
  if (!is.null(italic_col)) {
    j <- match(italic_col, header)
    if (is.na(j)) stop("No column '", italic_col, "' in ", path)
    body <- lapply(body, function(r) {
      if (nzchar(r[j]) && !grepl("^[-–]$", r[j])) r[j] <- sprintf("*%s*", r[j])
      r
    })
  }

  # These tables now live in the supplement, which is a standard article, not
  # the ASM class -- so no \\begin{fullwidth} or \\begin{tablenotes}, both of
  # which the ASM class provides and a plain article does not. Wide tables are
  # shrunk to the line width by adjustbox rather than left to run off the page;
  # narrow ones are left alone, since max width only ever scales down. The
  # widest are set sideways instead, where shrinking would leave them unreadable.
  env <- if (landscape) "sidewaystable" else "table"
  c(sprintf("%% generated by R/04_tables.R from %s -- do not edit", basename(path)),
    sprintf("%% needs \\usepackage{booktabs,adjustbox%s} in the preamble",
            if (landscape) ",rotating" else ""),
    sprintf("\\begin{%s}%s", env, if (landscape) "" else "[htbp]"),
    "\\centering",
    sprintf("\\caption{%s}%%", to_tex(title)),
    sprintf("\\label{%s}", label),
    if (!is.null(fontsize)) sprintf("\\%s", fontsize) else NULL,
    "\\begin{adjustbox}{max width=\\linewidth}",
    sprintf("\\begin{tabular}{%s}", paste(align, collapse = "")),
    "\\toprule",
    paste0(paste(sprintf("\\textbf{%s}", to_tex(header)), collapse = " & "), " \\\\"),
    "\\midrule",
    vapply(body, function(r) paste0(paste(to_tex(r), collapse = " & "), " \\\\"),
           character(1)),
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{adjustbox}",
    if (length(caption)) c("\\par\\smallskip",
                           sprintf("{\\footnotesize\\raggedright %s\\par}",
                                   to_tex(paste(caption, collapse = " ")))) else NULL,
    sprintf("\\end{%s}", env))
}

## ---- which tables ------------------------------------------------------------

TABLES <- list(
  list(out = "tab1_population_descriptors", project = "pop",
       file = "table2_descriptors.md", label = "tab:popdesc",
       fullwidth = TRUE, fontsize = "small"),
  list(out = "tab2_mutations_by_fraction", project = "mut",
       file = "table1_mutations_by_fraction.md", label = "tab:mutfrac"),
  list(out = "tab3_repeated_genes", project = "mut",
       file = "table2_repeated_genes.md", label = "tab:parallel",
       fullwidth = TRUE, fontsize = "small", italic_col = "Gene"),
  list(out = "tab4_growth_group_models", project = "growth",
       file = "table1_group_models.md", label = "tab:growthmodels",
       fontsize = "small", landscape = TRUE),   # ten columns
  list(out = "tab5_biofilm_group_models", project = "biofilm",
       file = "table1_group_models.md", label = "tab:biofilmmodels",
       fullwidth = TRUE, fontsize = "small"),
  list(out = "tab6_biofilm_lineage_groups", project = "biofilm",
       file = "table4_lineage_groups.md", label = "tab:biofilmlineage")
)

built <- 0L
for (t in TABLES) {
  src <- file.path(project_dir(t$project), "output", t$file)
  if (!file.exists(src)) {
    warning(sprintf("missing: %s -- run %s first", t$file, PROJECTS[[t$project]]),
            call. = FALSE)
    next
  }
  writeLines(convert_table(src, t$label, fullwidth = isTRUE(t$fullwidth),
                           fontsize = t$fontsize, italic_col = t$italic_col,
                           landscape = isTRUE(t$landscape)),
             file.path(TAB_DIR, paste0(t$out, ".tex")))
  built <- built + 1L
  cat(sprintf("  tables/%s.tex  <- %s/%s\n", t$out, PROJECTS[[t$project]], t$file))
}

## ---- stage the Overleaf upload ----------------------------------------------

# Laid out the way manuscript.tex refers to them: \input{tables/...} and
# \includegraphics{figures/...}. Drag both folders into the Overleaf project.
OVERLEAF <- file.path(MS_ROOT, "overleaf")
unlink(OVERLEAF, recursive = TRUE)
dir.create(file.path(OVERLEAF, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OVERLEAF, "tables"),  recursive = TRUE, showWarnings = FALSE)

# Stage the figures the documents include, not everything in output/. Staging
# every figure*.pdf uploaded the Figure 2 variants and the superseded four-panel
# version alongside the real one, which is how a stale Figure 2 went unnoticed.
# The documents live outside this repository; if they cannot be found, fall
# back to staging everything and say so.
TEXT_DIR <- path.expand(Sys.getenv("STARVATION_TEXT", "~/Desktop/6.Manuscript"))
docs <- file.path(TEXT_DIR, c("main.tex", "supplementary.tex"))
docs <- docs[file.exists(docs)]

wanted <- unique(unlist(lapply(docs, function(d) {
  x <- readLines(d, warn = FALSE)
  x <- sub("(?<!\\\\)%.*", "", x, perl = TRUE)
  m <- regmatches(x, gregexpr("\\\\includegraphics(\\[[^]]*\\])?\\{[^}]+\\}", x))
  basename(sub(".*\\{(.*)\\}$", "\\1", unlist(m)))
})))

if (length(wanted)) {
  have <- basename(list.files(OUT_DIR, pattern = "\\.pdf$"))
  missing <- setdiff(wanted, have)
  if (length(missing)) {
    warning("The documents include figure(s) that do not exist in output/: ",
            paste(missing, collapse = ", "), call. = FALSE)
  }
  figs <- file.path(OUT_DIR, intersect(wanted, have))
  cat(sprintf("  staging the %d figure(s) included by %s\n", length(figs),
              paste(basename(docs), collapse = " and ")))
} else {
  figs <- list.files(OUT_DIR, pattern = "^figure.*\\.pdf$", full.names = TRUE)
  cat("  documents not found under", TEXT_DIR, "-- staging every figure*.pdf\n")
}
file.copy(figs, file.path(OVERLEAF, "figures"), overwrite = TRUE)
file.copy(list.files(TAB_DIR, pattern = "\\.tex$", full.names = TRUE),
          file.path(OVERLEAF, "tables"), overwrite = TRUE)

cat(sprintf("\nWrote %d table(s), and staged overleaf/ with %d figure(s).\n",
            built, length(figs)))
cat("Upload overleaf/figures and overleaf/tables to the Overleaf project.\n")
