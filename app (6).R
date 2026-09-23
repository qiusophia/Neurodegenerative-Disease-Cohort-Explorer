library(shiny)
library(ggplot2)
library(dplyr)
library(stringr)
library(ggrepel)

read_tab <- function(path)
  read.delim(path, sep = "\t", header = TRUE, check.names = FALSE,
             stringsAsFactors = FALSE, quote = "", comment.char = "")

df_raw <- read_tab("lab.txt")

parse_tx <- function(df, idx) {
  df %>%
    transmute(
      canon     = str_split_fixed(RefSeq_AA_change, ",", 2)[, idx],
      canon_pos = as.integer(str_match(
        str_split_fixed(str_split_fixed(RefSeq_AA_change, ",", 2)[, 1], ":", 5)[, 5],
        "[A-Z](\\d+)")[, 2]),
      clinvar = ClinVar_clinical_significance,
      af      = coalesce(suppressWarnings(as.numeric(AllOfUs_AF_all)), 0),
      cadd    = coalesce(suppressWarnings(as.numeric(CADD_phred)), 0),
      revel   = coalesce(suppressWarnings(as.numeric(REVEL_score)), 0)
    )
}

ALZFORUM_CURATION <- c(
  "M84V" = "Pathogenic", "V103_S104delinsG" = "Pathogenic", "T116I" = "Pathogenic",
  "P117A" = "Pathogenic", "T119I" = "Pathogenic", "M139V" = "Pathogenic",
  "I143T" = "Pathogenic", "M146L" = "Pathogenic", "H163R" = "Pathogenic",
  "L173F" = "Pathogenic", "G209R" = "Pathogenic", "Q223K" = "Pathogenic",
  "P264L" = "Pathogenic", "E280A" = "Pathogenic", "P284L" = "Pathogenic",
  "I416T" = "Pathogenic",
  "I162S" = "Likely pathogenic", "I227L" = "Likely pathogenic", "A431V" = "Likely pathogenic",
  "V94M"  = "Uncertain significance", "T354I" = "Uncertain significance",
  "I439V" = "Uncertain significance", "I408T" = "Uncertain significance",
  "E318G" = "Benign", "I427V" = "Benign", "T291A" = "Benign", "D333G" = "Benign",
  "P355S" = "Likely Benign",
  "V412I" = "Not curated"
)

variants <- bind_rows(parse_tx(df_raw, 1), parse_tx(df_raw, 2)) %>%
  filter(canon != "" & !is.na(canon)) %>%
  mutate(
    gene       = str_split_fixed(canon, ":", 5)[, 1],
    transcript = str_split_fixed(canon, ":", 5)[, 2],
    exon       = as.integer(str_remove(str_split_fixed(canon, ":", 5)[, 3], "exon")),
    aa_change  = str_split_fixed(canon, ":", 5)[, 5],
    p          = str_remove(aa_change, "^p\\."),
    ref_aa     = str_match(p, "^([A-Z])")[, 2],
    position   = as.integer(str_match(p, "^[A-Z](\\d+)")[, 2]),
    alt_aa     = str_match(p, "^[A-Z]\\d+([A-Z*])$")[, 2],
    label      = if_else(str_detect(p, "delins|del|ins|dup"),
                         p, paste0(ref_aa, position, alt_aa)),
    txid       = paste0(gene, ": ", transcript),
    alzforum   = coalesce(unname(ALZFORUM_CURATION[label]), "Not reported"),
    clin_cat = case_when(
      clinvar %in% c("Pathogenic", "Pathogenic/Likely_pathogenic")     ~ "Pathogenic",
      clinvar == "Likely_pathogenic"                                   ~ "Likely pathogenic",
      clinvar == "Uncertain_significance"                              ~ "Uncertain significance",
      clinvar == "Conflicting_classifications_of_pathogenicity"        ~ "Conflicting",
      clinvar %in% c("Benign", "Likely_benign", "Benign/Likely_benign")~ "Benign / Likely benign",
      TRUE                                                              ~ "Not reported"
    )
  ) %>%
  filter(transcript == "NM_000021")

CASE_GROUP    <- "HALB"
CONTROL_GROUP <- "SL"

QC_KEEP_VALUE <- "Keep"

DIAG_LEVELS <- c("Alzheimer's Disease", "Frontotemporal Dementia",
                 "Amyotrophic Lateral Sclerosis", "Early Onset Dementia",
                 "Healthy at evaluation", "Other", "Missing")
DIAG_COLORS <- c("Alzheimer's Disease"           = "#f2994a",
                 "Frontotemporal Dementia"       = "#bb6bd9",
                 "Amyotrophic Lateral Sclerosis" = "#eb5757",
                 "Early Onset Dementia"          = "#f2c94c",
                 "Healthy at evaluation"         = "#6fcf97",
                 Other = "#9aa0a6", Missing = "#6f6f6f")

meta_raw <- read_tab("filter_individuals.txt")
names(meta_raw)[1:6] <- c("vcf_ID", "phenotype", "diagnosis", "aoo", "aoe", "cohort")
meta <- meta_raw %>%
  transmute(
    vcf_ID    = as.character(vcf_ID),
    qc        = QC_KEEP_VALUE,
    diagnosis = dplyr::recode(as.character(diagnosis),
                              "AD"      = "Alzheimer's Disease",
                              "FTD"     = "Frontotemporal Dementia",
                              "ALS"     = "Amyotrophic Lateral Sclerosis",
                              "EOD"     = "Early Onset Dementia",
                              "AFM"     = "Healthy at evaluation",
                              "CN"      = "Healthy at evaluation",
                              "CONTROL" = "Healthy at evaluation",
                              "Control" = "Healthy at evaluation"),
    cohort    = as.character(cohort),
    aoo       = suppressWarnings(as.numeric(aoo)),
    aoe       = suppressWarnings(as.numeric(aoe))
  ) %>%
  filter(!is.na(vcf_ID))

nm21_pos <- function(s) {
  if (is.na(s)) return(NA_integer_)
  parts <- str_split(s, ",")[[1]]
  hit   <- parts[str_detect(parts, "NM_000021")]
  if (!length(hit)) return(NA_integer_)
  p <- str_split_fixed(hit[1], ":", 5)[, 5]
  as.integer(str_match(p, "[A-Z](\\d+)")[, 2])
}
variant_pos <- vapply(df_raw$RefSeq_AA_change, nm21_pos, integer(1))

sample_cols <- intersect(names(df_raw), meta$vcf_ID)

GTm         <- as.matrix(df_raw[, sample_cols])
storage.mode(GTm) <- "character"
carrier_ind <- matrix(str_detect(GTm, "[1-9]"),
                      nrow = nrow(GTm), dimnames = list(NULL, sample_cols))
carrier_ind[is.na(carrier_ind)] <- FALSE

ci <- which(carrier_ind, arr.ind = TRUE)
carrier_df <- data.frame(
  position = variant_pos[ci[, 1]],
  vcf_ID   = sample_cols[ci[, 2]],
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(position)) %>%
  left_join(meta, by = "vcf_ID") %>%
  mutate(prefix = case_when(
    str_starts(vcf_ID, CASE_GROUP)    ~ "Cases",
    str_starts(vcf_ID, CONTROL_GROUP) ~ "Controls",
    TRUE                              ~ "Other prefix"
  ))

DIAG_CHOICES <- intersect(DIAG_LEVELS, unique(meta$diagnosis))
if (length(DIAG_CHOICES) == 0) DIAG_CHOICES <- sort(unique(meta$diagnosis))
DIAG_DEFAULT <- DIAG_CHOICES

COHORT_CHOICES <- sort(unique(meta$cohort))
COHORT_DEFAULT <- COHORT_CHOICES
cohort_idx     <- setNames(seq_along(COHORT_CHOICES), COHORT_CHOICES)
COHORT_NAMES   <- lapply(COHORT_CHOICES, function(x)
  HTML(paste0(x, "<sup>", cohort_idx[[x]], "</sup>")))

fin <- function(x, d) { x <- x[is.finite(x)]; if (length(x)) x else d }
AOO_MIN <- floor(min(fin(meta$aoo, 0)));   AOO_MAX <- ceiling(max(fin(meta$aoo, 100)))
AOE_MIN <- floor(min(fin(meta$aoe, 0)));   AOE_MAX <- ceiling(max(fin(meta$aoe, 100)))

PROTEIN_LENGTHS <- c(
  "PSEN1: NM_000021" = 467
)

af_pos    <- variants$af[variants$af > 0]
AF_MIN    <- floor(log10(min(af_pos)))
AF_MAX    <- ceiling(log10(max(af_pos)))
CADD_MIN  <- 0
CADD_MAX  <- 1

base_colors <- c("#ff9999", "#ffcc99", "#ffff99", "#ccff99", "#99ff99",
                 "#99ffff", "#99ccff", "#9999ff", "#e6b3ff", "#ffb3e6")

CLIN_CHOICES <- c("Pathogenic", "Likely pathogenic", "Uncertain significance",
                  "Conflicting", "Benign / Likely benign", "Not reported")
ALZ_CHOICES  <- c("Pathogenic", "Likely pathogenic", "Uncertain significance",
                  "Benign", "Likely Benign", "Not curated", "Not reported")
TX_CHOICES   <- sort(unique(variants$txid))

ui <- fluidPage(
  tags$head(tags$style(HTML("
    @import url('https://fonts.googleapis.com/css2?family=Playfair+Display:wght@600;700&family=Inter:wght@400;500;600&display=swap');

    body {
      background: radial-gradient(circle at top, #1c1a14 0%, #080808 70%);
      color: #e8e0c5;
      font-family: 'Inter', system-ui, -apple-system, sans-serif;
    }
    .container-fluid { padding-top: 10px; }

    /* ---- Title ---- */
    .app-title {
      font-family: 'Playfair Display', Georgia, serif;
      font-size: 42px;
      font-weight: 700;
      letter-spacing: 1px;
      color: #d4af37;
      text-align: center;
      margin: 18px 0 6px 0;
      text-shadow: 0 2px 12px rgba(212,175,55,0.30);
    }
    .app-title:after {
      content: '';
      display: block;
      width: 140px; height: 2px;
      margin: 12px auto 0 auto;
      background: linear-gradient(90deg, transparent, #d4af37, transparent);
    }

    /* ---- Control panel ---- */
    .well {
      background: linear-gradient(180deg, #1a1813 0%, #121212 100%);
      border: 1px solid #d4af37;
      border-radius: 12px;
      box-shadow: 0 6px 28px rgba(0,0,0,0.65), inset 0 0 0 1px rgba(212,175,55,0.12);
      padding: 20px 22px;
    }

    label, .control-label {
      color: #d4af37 !important;
      font-weight: 600;
      letter-spacing: 0.3px;
    }

    select.form-control, .selectize-input, select, .selectize-dropdown {
      background-color: #0e0e0e !important;
      color: #f0e6c8 !important;
      border: 1px solid #6b5a1e !important;
      border-radius: 6px !important;
    }
    .selectize-dropdown .active { background-color: #2a2410 !important; color: #f0e6c8 !important; }

    /* ---- Checkboxes ---- */
    input[type=checkbox] { accent-color: #d4af37; width: 15px; height: 15px; }
    .checkbox label, .shiny-options-group label { color: #e8e0c5 !important; font-weight: 400; }

    /* ---- Sliders (ionRangeSlider) ---- */
    .irs-bar, .irs-bar-edge, .irs--shiny .irs-bar, .irs--shiny .irs-bar-edge {
      background: linear-gradient(90deg, #b8860b, #d4af37) !important;
      border-color: #d4af37 !important;
    }
    .irs-line, .irs--shiny .irs-line { background: #2a2a2a !important; border-color: #2a2a2a !important; }
    .irs-handle, .irs--shiny .irs-handle { border: 1px solid #d4af37 !important; }
    .irs--shiny .irs-handle > i:first-child { background: #d4af37 !important; }
    .irs-from, .irs-to, .irs-single,
    .irs--shiny .irs-from, .irs--shiny .irs-to, .irs--shiny .irs-single {
      background: #d4af37 !important; color: #0a0a0a !important; font-weight: 600;
    }
    .irs-min, .irs-max, .irs-grid-text,
    .irs--shiny .irs-min, .irs--shiny .irs-max, .irs--shiny .irs-grid-text {
      color: #9a8c5c !important; background: transparent !important;
    }

    #afthresh { color: #c9b676; font-style: italic; margin-top: 10px; }

    .panel-heading-lux {
      font-family:'Playfair Display',Georgia,serif; font-weight:700;
      color:#d4af37; letter-spacing:0.5px; margin:0 0 16px 0;
      text-shadow:0 2px 8px rgba(212,175,55,0.25);
    }

    .shiny-plot-output { border-radius: 12px; overflow: hidden; margin-top: 14px; }

    /* ---- Footnotes ---- */
    .page-footnotes { color:#9a8c5c; font-size:12.5px; line-height:1.55; margin:8px 6px 34px 6px; }
    .page-footnotes p { margin:0 0 10px 0; }
    .page-footnotes sup { color:#d4af37; font-weight:700; }
    .page-footnotes b { color:#c9b676; }
    .page-footnotes a { color:#d4af37; text-decoration:none; word-break:break-word; }
    .page-footnotes a:hover { text-decoration:underline; }
  "))),
  titlePanel(div(class = "app-title", "Variant Visualizer"),
             windowTitle = "Variant Visualizer"),

  wellPanel(
    selectInput("transcript", "Gene",
                choices = TX_CHOICES, selected = TX_CHOICES[1])
  ),

  uiOutput("lolli_ui"),

  wellPanel(
    tags$h4("Filter Individuals", class = "panel-heading-lux"),
    fluidRow(
      column(6, checkboxGroupInput("ind_diag", "Diagnosis",
                                   choices = DIAG_CHOICES, selected = DIAG_DEFAULT, inline = TRUE)),
      column(6, checkboxGroupInput("ind_cohort", "Cohort",
                                   choiceNames = COHORT_NAMES, choiceValues = as.list(COHORT_CHOICES),
                                   selected = COHORT_DEFAULT, inline = TRUE))
    ),
    fluidRow(
      column(6,
             sliderInput("aoo_rng", "Age at onset (applies to patients with reported data)",
                         min = AOO_MIN, max = AOO_MAX, value = c(AOO_MIN, AOO_MAX), step = 1),
             checkboxInput("aoo_na", "Include individuals with unknown age at onset", TRUE)),
      column(6,
             sliderInput("aoe_rng", "Age at evaluation",
                         min = AOE_MIN, max = AOE_MAX, value = c(AOE_MIN, AOE_MAX), step = 1),
             checkboxInput("aoe_na", "Include individuals with unknown age at evaluation", TRUE))
    )
  ),

  wellPanel(
    tags$h4("Filter Variants", class = "panel-heading-lux"),
    checkboxGroupInput("clin", "ClinVar curation",
                       choices = CLIN_CHOICES, selected = CLIN_CHOICES, inline = TRUE),
    checkboxGroupInput("alz", "Alzforum curation",
                       choices = ALZ_CHOICES, selected = ALZ_CHOICES, inline = TRUE),
    sliderInput("af", "Max allele frequency (log10, AllOfUs all)",
                min = AF_MIN, max = AF_MAX, value = AF_MAX, step = 0.1),
    checkboxInput("incl_absent",
                  "Include variants not observed in AllOfUs (frequency = 0)", value = TRUE),
    fluidRow(
      column(6, sliderInput("cadd", "CADD (0\u20131)",
                            min = CADD_MIN, max = CADD_MAX, value = c(CADD_MIN, CADD_MAX), step = 0.1)),
      column(6, sliderInput("revel", "REVEL",
                            min = 0, max = 1, value = c(0, 1), step = 0.05))
    ),
    textOutput("afthresh")
  ),

  tags$div(
    class = "page-footnotes",
    tags$hr(style = "border-color:#3a3320; margin-top:22px;"),
    do.call(tagList, lapply(COHORT_CHOICES, function(x)
      tags$p(HTML(paste0(
        "<sup>", cohort_idx[[x]], "</sup> <b>", x,
        "</b> &mdash; synthetic example cohort included with the demo dataset."
      )))))
  )
)

server <- function(input, output, session) {

  tx_variants <- reactive({
    variants %>% filter(txid == input$transcript)
  })

  protein_len <- reactive({
    tx  <- tx_variants()
    key <- input$transcript
    if (!is.null(key) && key %in% names(PROTEIN_LENGTHS)) {
      PROTEIN_LENGTHS[[key]]
    } else {
      fallback <- suppressWarnings(max(tx$position, na.rm = TRUE))
      if (!is.finite(fallback)) 500 else fallback
    }
  })

  exon_map <- reactive({
    plen <- protein_len()
    em <- tx_variants() %>%
      filter(!is.na(exon) & !is.na(position)) %>%
      group_by(exon) %>%
      summarise(min_pos = min(position), max_pos = max(position), .groups = "drop") %>%
      arrange(exon)

    if (nrow(em) > 1) {
      em$x0 <- em$min_pos
      em$x1 <- em$max_pos
      for (i in 1:(nrow(em) - 1)) {
        mid <- floor((em$max_pos[i] + em$min_pos[i + 1]) / 2)
        em$x1[i]     <- mid
        em$x0[i + 1] <- mid
      }
      em$x0[1] <- 1
      em$x1[nrow(em)] <- plen
    } else if (nrow(em) == 1) {
      em$x0 <- 1
      em$x1 <- plen
    }
    if (nrow(em) > 0) em$color <- rep(base_colors, length.out = nrow(em))
    em
  })

  output$afthresh <- renderText({
    sprintf("Showing variants with allele frequency \u2264 %.2e", 10^input$af)
  })

  inc_carriers <- reactive({
    d <- carrier_df
    d <- d[d$qc == QC_KEEP_VALUE &
             d$diagnosis %in% input$ind_diag &
             d$cohort %in% input$ind_cohort, , drop = FALSE]
    keep_aoo <- ifelse(is.na(d$aoo), isTRUE(input$aoo_na),
                       d$aoo >= input$aoo_rng[1] & d$aoo <= input$aoo_rng[2])
    keep_aoe <- ifelse(is.na(d$aoe), isTRUE(input$aoe_na),
                       d$aoe >= input$aoe_rng[1] & d$aoe <= input$aoe_rng[2])
    d[keep_aoo & keep_aoe, , drop = FALSE]
  })

  pdata <- reactive({
    em   <- exon_map()
    plen <- protein_len()

    df <- tx_variants() %>%
      filter(!is.na(position), clin_cat %in% input$clin, alzforum %in% input$alz) %>%
      filter(af <= 10^input$af) %>%
      filter(input$incl_absent | af > 0) %>%
      filter(cadd  >= input$cadd[1]  & cadd  <= input$cadd[2]) %>%
      filter(revel >= input$revel[1] & revel <= input$revel[2]) %>%
      group_by(position) %>%
      summarise(label = paste(unique(label), collapse = " | "), .groups = "drop") %>%
      arrange(position)

    if (nrow(df) >= 1) {
      gap <- 8.2
      xd  <- df$position
      n   <- nrow(df)
      if (n >= 2) {
        for (i in 2:n) if (xd[i] - xd[i - 1] < gap) xd[i] <- xd[i - 1] + gap
        if (xd[n] > plen) {
          xd <- xd - (xd[n] - plen)
          for (i in (n - 1):1) if (xd[i + 1] - xd[i] < gap) xd[i] <- xd[i + 1] - gap
        }
      }
      df$xd <- xd
    }

    inc <- inc_carriers()
    posset <- df$position
    inc_disp <- inc[inc$position %in% posset, , drop = FALSE]

    present   <- unique(inc_disp$diagnosis)
    diag_rows <- c(intersect(DIAG_LEVELS, present), setdiff(present, c(DIAG_LEVELS, NA)))
    diag_rows <- diag_rows[!is.na(diag_rows)]
    tracks    <- diag_rows

    mk <- function(sub, trackname) {
      cnt <- as.data.frame(table(position = sub$position), stringsAsFactors = FALSE)
      if (nrow(cnt)) cnt$position <- as.integer(cnt$position)
      out <- data.frame(position = posset, track = trackname, stringsAsFactors = FALSE)
      out$n <- cnt$Freq[match(out$position, cnt$position)]
      out$n[is.na(out$n)] <- 0L
      out
    }
    parts <- list()
    for (d in diag_rows)
      parts[[length(parts) + 1]] <- mk(inc_disp[inc_disp$diagnosis == d, , drop = FALSE], d)
    counts <- bind_rows(parts)

    if (nrow(df) > 0 && nrow(counts) > 0)
      counts <- left_join(counts, df[, c("position", "xd")], by = "position")

    list(df = df, em = em, plen = plen, tracks = tracks, counts = counts)
  })

  output$lolli_ui <- renderUI({
    n_tr <- length(pdata()$tracks)
    h    <- 470 + n_tr * 72
    tags$div(style = "overflow-x: auto;",
             plotOutput("lolli", width = "1500px", height = paste0(h, "px")))
  })

  output$lolli <- renderPlot({
    pd     <- pdata()
    em     <- pd$em; plen <- pd$plen
    df     <- pd$df; counts <- pd$counts
    tracks <- pd$tracks
    n_tr   <- length(tracks)

    Y_HEAD <- 1; Y_LINK <- 1.40; Y_NAME <- 1.50
    Y_T0   <- 2.50; Y_STEP <- 0.30
    track_y <- setNames(Y_T0 + (seq_len(n_tr) - 1) * Y_STEP, tracks)
    Y_TOP   <- Y_T0 + (n_tr - 1) * Y_STEP + 0.40

    track_col <- c(Cases = "#f2c94c", Controls = "#7ec8f0", DIAG_COLORS)
    col_for   <- function(tr) if (tr %in% names(track_col)) unname(track_col[[tr]]) else "#cfcfcf"

    lab_x <- -plen * 0.03

    p <- ggplot()

    if (nrow(em) > 0) {
      p <- p +
        geom_segment(data = em,
                     aes(x = x0, xend = x1, y = 0, yend = 0, color = as.factor(exon)),
                     linewidth = 8, inherit.aes = FALSE) +
        scale_color_manual(values = setNames(em$color, em$exon)) +
        geom_text(data = em,
                  aes(x = (x0 + x1) / 2, y = 0, label = paste("Exon", exon)),
                  size = 3.4, color = "#2d3436", fontface = "bold")
    }

    if (nrow(df) > 0) {
      p <- p +
        geom_segment(data = df, aes(x = position, xend = position, y = 0, yend = Y_HEAD),
                     color = "#8a8f96", linewidth = 0.8) +
        geom_point(data = df, aes(x = position, y = Y_HEAD),
                   size = 2.8, color = "#d4af37") +
        geom_segment(data = df, aes(x = position, xend = xd, y = Y_HEAD + 0.03, yend = Y_LINK),
                     color = "#6f6a52", linewidth = 0.4) +
        geom_text(data = df, aes(x = xd, y = Y_NAME, label = label),
                  angle = 90, hjust = 0, color = "#e8e0c5", size = 3.7)

      for (tr in tracks) {
        cdat <- counts[counts$track == tr & !is.na(counts$xd) & counts$n > 0, , drop = FALSE]
        col  <- col_for(tr)
        if (nrow(cdat) > 0) {
          p <- p + geom_label(data = cdat, aes(x = xd, label = n), y = track_y[[tr]],
                              size = 4.5, fontface = "bold", color = col, fill = "#141414",
                              label.size = 0.4, label.padding = unit(0.14, "lines"))
        }
        p <- p + annotate("text", x = lab_x, y = track_y[[tr]], label = tr,
                          hjust = 1, fontface = "bold", color = col, size = 5.0)
      }
    }

    p +
      labs(x = input$transcript) +
      scale_x_continuous(breaks = seq(0, plen, 50)) +
      coord_cartesian(xlim = c(0, plen + 5), ylim = c(-0.2, Y_TOP), clip = "off") +
      theme_minimal(base_size = 14) +
      theme(
        legend.position  = "none",
        plot.background  = element_rect(fill = "#111111", color = NA),
        panel.background = element_rect(fill = "#111111", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y     = element_blank(),
        axis.text.y      = element_blank(),
        axis.ticks.y     = element_blank(),
        axis.text.x      = element_text(color = "#c9b676"),
        axis.ticks.x     = element_line(color = "#7a6f4e"),
        axis.title.x     = element_text(face = "bold", color = "#d4af37", margin = margin(t = 15)),
        plot.margin      = margin(t = 10, r = 12, b = 10, l = 260)
      )
  })
}

shinyApp(ui, server)
