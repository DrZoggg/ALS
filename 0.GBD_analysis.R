# Integrate files from multiple folders
library(dplyr)

# Set working directory
setwd("GBD_analysis")

# Helper: list all CSVs under a folder and row-bind them
read_and_bind <- function(dir_path) {
  csv_files <- list.files(path = dir_path, pattern = "\\.csv$", full.names = TRUE)
  lapply(csv_files, read.csv) %>% bind_rows()
}

# -----------------------------
# Cause of death or injury/*
# -----------------------------
DALYs <- read_and_bind("Cause of death or injury/DALYs/")
str(DALYs)

Deaths <- read_and_bind("Cause of death or injury/Deaths/")
str(Deaths)

Incidence <- read_and_bind("Cause of death or injury/Incidence/")
str(Incidence)

Prevalence <- read_and_bind("Cause of death or injury/Prevalence/")
str(Prevalence)

YLDs <- read_and_bind("Cause of death or injury/YLDs/")
str(YLDs)

YLLs <- read_and_bind("Cause of death or injury/YLLs/")
str(YLLs)

# ----------------------------------------
# Probability of death, by cause/*
# ----------------------------------------
Probability <- read_and_bind("Probability of death, by cause/")
str(Probability)

# ---------------------------------------------------------------
# Casuse of death or injury(total percentage change)/*
# (Note: keeping folder names exactly as provided; no changes)
# ---------------------------------------------------------------
DALYsPercentChange <- read_and_bind("Casuse of death or injury(total percentage change)/DALYs/")
str(DALYsPercentChange)

DeathsPercentChange <- read_and_bind("Casuse of death or injury(total percentage change)/Deaths/")
str(DeathsPercentChange)

IncidencePercentChange <- read_and_bind("Casuse of death or injury(total percentage change)/Incidence/")
str(IncidencePercentChange)

MMRPercentChange <- read_and_bind("Casuse of death or injury(total percentage change)/MMR/")
str(MMRPercentChange)

PrevalencePercentChange <- read_and_bind("Casuse of death or injury(total percentage change)/Prevalence/")
str(PrevalencePercentChange)

# Intentionally repeat this block to preserve the original script's behavior
PrevalencePercentChange <- read_and_bind("Casuse of death or injury(total percentage change)/Prevalence/")
str(PrevalencePercentChange)

YLDsPercentChange <- read_and_bind("Casuse of death or injury(total percentage change)/YLDs/")
str(YLDsPercentChange)

YLLsPercentChange <- read_and_bind("Casuse of death or injury(total percentage change)/YLLs/")
str(YLLsPercentChange)

# Combine and save files
MND <- do.call(rbind, list(
  DALYs, Deaths, Incidence, Prevalence,
  Probability, YLDs, YLLs
))
str(MND)

MND_PercentChange <- do.call(rbind, list(
  DALYsPercentChange,
  DeathsPercentChange,
  IncidencePercentChange,
  MMRPercentChange,
  PrevalencePercentChange,
  YLDsPercentChange,
  YLLsPercentChange
))
str(MND_PercentChange)

# Use qs for fast compressed saving
library(qs)
# Save objects with compression
qsave(MND, "MND.qs", nthreads = 16, compress_level = 9)               
qsave(MND_PercentChange, "MND_PercentChange.qs", nthreads = 16, compress_level = 9)



# -------------------------------
# Calculate EAPC values
# -------------------------------
library(tidyverse)
library(qs)

MND <- qread("MND.qs")
MND_PercentChange <- qread("MND_PercentChange.qs")


# Check available measure names
unique(MND$measure_name)
unique(MND_PercentChange$measure_name)

# Select measures of interest
selected_measure_name <- c(
  "DALYs (Disability-Adjusted Life Years)",
  "Deaths",
  "Incidence",
  "Prevalence",
  "YLDs (Years Lived with Disability)",
  "YLLs (Years of Life Lost)"
)
selected_measure_name

# Inspect column names
colnames(MND)
colnames(MND_PercentChange)


# --------------------------------
# Extract data for EAPC calculation
# --------------------------------
df1 <- MND %>%
  filter(measure_name %in% selected_measure_name) %>%
  filter(age_name == "Age-standardized") %>%
  filter(!metric_name == "Percent") %>%
  mutate(
    id = paste0(
      measure_id, measure_name, location_id, location_name,
      sex_id, sex_name, age_id, age_name,
      cause_id, cause_name, metric_id, metric_name
    )
  )

# Calculate EAPC using linear regression
library(broom)
decx <- 3

df_eapc <- df1 %>%
  group_by(id) %>%
  do(tidy(lm(log(val) ~ year, data = .), conf.int = TRUE)) %>%
  mutate(
    EAPC    = 100 * (exp(estimate) - 1),
    lower_CI = 100 * (exp(conf.low) - 1),
    upper_CI = 100 * (exp(conf.high) - 1)
  ) %>%
  filter(term %in% "year") %>%
  mutate(
    EAPCs = sprintf(
      paste0("%.", decx, "f(%.", decx, "f,%.", decx, "f)"),
      EAPC, lower_CI, upper_CI
    )
  )


# ----------------------------------------
# Extract 1990 and 2021 values for display
# ----------------------------------------
df1 <- MND %>%
  filter(age_name == "Age-standardized") %>%
  filter(year %in% c(1990, 2021)) %>%
  filter(!metric_name == "Percent") %>%
  filter(measure_name %in% selected_measure_name) %>%
  mutate(
    num = sprintf(
      paste0("%.", 3, "f(%.", 3, "f-%.", 3, "f)"),
      val, lower, upper
    )
  ) %>%
  select(contains("name"), contains("id"), year, num)

dfa <- df1 %>%
  pivot_wider(names_from = year, values_from = num) %>%
  mutate(
    id = paste0(
      measure_id, measure_name, location_id, location_name,
      sex_id, sex_name, age_id, age_name,
      cause_id, cause_name, metric_id, metric_name
    )
  ) %>%
  select(id, contains("name"), contains("id"), "1990", "2021")


# ---------------------------
# Extract Percent Change data
# ---------------------------
df2 <- MND_PercentChange %>%
  filter(age_name == "Age-standardized") %>%
  filter(year_start %in% "1990") %>%
  filter(year_end %in% "2021") %>%
  filter(!metric_name == "Percent") %>%
  filter(measure_name %in% selected_measure_name) %>%
  mutate(
    PC = sprintf(
      paste0("%.", 3, "f(%.", 3, "f-%.", 3, "f)"),
      val, lower, upper
    )
  ) %>%
  select(contains("name"), contains("id"), PC) %>%
  mutate(
    id = paste0(
      measure_id, measure_name, location_id, location_name,
      sex_id, sex_name, age_id, age_name,
      cause_id, cause_name, metric_id, metric_name
    )
  )

dfb <- df2 %>% select(id, PC)


# --------------------------
# Merge data frames together
# --------------------------
dfall <- left_join(dfa, dfb, multiple = "first")
head(dfall, 3)

df_all <- left_join(dfall, df_eapc) %>%
  select(-id)
head(df_all)

saveRDS(df_all, file = "MND_all_age.standarded.rds")

# Save as Excel
library(openxlsx)
write.xlsx(df_all, file = "Age-standardized_EAPC_PC_1990_2021_MND.xlsx")


# -------- Trend plots (GBD-based), faceted by SDI groups ----------
library(tidyverse)
library(scales)

# Source data
df <- MND

# Locations to include (SDI groups + Global = 1)
sdi_ids <- c(44634, 44635, 44636, 44637, 44639, 1)

# A small helper: build a readable y-axis label
ys_label <- function(measure) {
  paste0("Age-standardized ", measure, " rate per 100,000 population")
}

# A small helper: filter data for a given measure under your fixed conditions
filter_measure <- function(measure) {
  df %>%
    filter(
      measure_name == measure,
      location_id %in% sdi_ids,
      cause_name == "Motor neuron disease",
      age_id == 27,                # Age-standardized in GBD
      metric_name == "Rate",
      sex_name == "Both"
    )
}

# Ensure an output directory exists
out_dir <- "Trend_Plot"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# --- generic plotting function with title/subtitle/caption ---
plot_trend_faceted <- function(data, ylab_text,
                               title = NULL,
                               subtitle = NULL,
                               caption = NULL) {
  ggplot(data, aes(x = year, y = val, fill = location_name)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.20, show.legend = TRUE) +
    geom_line(aes(color = location_name, linetype = location_name), linewidth = 0.7) +
    facet_wrap(~ location_name, scales = "free") +
    theme_classic(base_size = 12) +
    scale_y_continuous(labels = scales::label_comma()) +
    scale_x_continuous(
      breaks = seq(min(data$year, na.rm = TRUE), max(data$year, na.rm = TRUE), by = 5),
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      caption = caption,
      x = "Year",
      y = ylab_text,
      color = "", linetype = "", fill = ""
    ) +
    guides(
      color = guide_legend(title = "", nrow = 1),
      linetype = guide_legend(title = "", nrow = 1),
      fill = "none"
    ) +
    theme(
      plot.title.position = "plot",
      plot.title   = element_text(face = "bold", size = 14, hjust = 0),
      plot.subtitle= element_text(size = 12, colour = "grey30"),
      plot.caption = element_text(size = 10, colour = "grey40"),
      legend.position = "top",
      legend.title = element_text(size = 12),
      legend.text  = element_text(size = 12),
      strip.background = element_blank(),
      strip.text = element_text(size = 12),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      axis.text.x  = element_text(size = 12, angle = 60, hjust = 1, vjust = 1,
                                  margin = margin(t = 2)),
      axis.text.y  = element_text(size = 12),
      plot.margin  = margin(t = 8, r = 8, b = 14, l = 8)
    )
}

ys_label <- function(measure) {
  paste0("Age-standardized ", measure, " rate per 100,000 population")
}

# Prevalence
df_prev <- filter_measure("Prevalence")
p_prev <- plot_trend_faceted(
  df_prev,
  ylab_text = ys_label("Prevalence"),
  title     = "Age-standardized prevalence rate of motor neuron disease by SDI, 1990–2021 (GBD)",
  subtitle  = "Lines show mean estimates; ribbons show 95% uncertainty intervals. Sex: Both; Metric: rate per 100,000."
)
ggsave(file.path(out_dir, "trend_prevalence_age-standardized_rate.pdf"),
       p_prev, width = 9, height = 6, units = "in")

# DALYs
df_dalys <- filter_measure("DALYs (Disability-Adjusted Life Years)")
p_dalys <- plot_trend_faceted(
  df_dalys,
  ylab_text = ys_label("DALYs"),
  title     = "Age-standardized DALY rate of motor neuron disease by SDI, 1990–2021 (GBD)",
  subtitle  = "Lines: mean; ribbons: 95% UI. Sex: Both; Metric: rate per 100,000."
)
ggsave(file.path(out_dir, "trend_dalys_age-standardized_rate.pdf"),
       p_dalys, width = 9, height = 6, units = "in")

# Deaths
df_death <- filter_measure("Deaths")
p_death <- plot_trend_faceted(
  df_death,
  ylab_text = ys_label("Deaths"),
  title     = "Age-standardized death rate of motor neuron disease by SDI, 1990–2021 (GBD)",
  subtitle  = "Lines: mean; ribbons: 95% UI. Sex: Both; Metric: rate per 100,000."
)
ggsave(file.path(out_dir, "trend_deaths_age-standardized_rate.pdf"),
       p_death, width = 9, height = 6, units = "in")

# Incidence
df_incid <- filter_measure("Incidence")
p_incid <- plot_trend_faceted(
  df_incid,
  ylab_text = ys_label("Incidence"),
  title     = "Age-standardized incidence rate of motor neuron disease by SDI, 1990–2021 (GBD)",
  subtitle  = "Lines: mean; ribbons: 95% UI. Sex: Both; Metric: rate per 100,000."
)
ggsave(file.path(out_dir, "trend_incidence_age-standardized_rate.pdf"),
       p_incid, width = 9, height = 6, units = "in")

# ---- export the plotting data to Excel (one file per figure) ----
library(openxlsx)
write.xlsx(df_prev,  file.path(out_dir, "data_prevalence_age-standardized_rate.xlsx"),  overwrite = TRUE)
write.xlsx(df_dalys, file.path(out_dir, "data_dalys_age-standardized_rate.xlsx"),       overwrite = TRUE)
write.xlsx(df_death, file.path(out_dir, "data_deaths_age-standardized_rate.xlsx"),      overwrite = TRUE)
write.xlsx(df_incid, file.path(out_dir, "data_incidence_age-standardized_rate.xlsx"),   overwrite = TRUE)
getwd()




####################################
# MAP SECTION
####################################
setwd('GBD_analysis')
library(tidyverse)
library(readxl)
library(tidyr)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(mapdata)
library(openxlsx)

# Load helper objects and GBD-related data (e.g., df_world, namex, etc.)
load("function&GBDdata/GBD.Rdata")

# Load main data using qs (compressed format)
library(qs)
MND <- qread("MND.qs")
MND_PercentChange <- qread("MND_PercentChange.qs")


# ================================
# Percent Change (PC) WORLD MAPS
# ================================
# Note: In this section, PC values are multiplied by 100 (later expressed as percentage).
#       Rounding is kept to 1 decimal place to match the original behavior.

setwd("Percent Change across 204 countries")

for (dn in setdiff(unique(MND_PercentChange[["measure_name"]]), "Maternal mortality ratio")) {  # Exclude MMR for this analysis
  labelx <- dn
  
  dfx <- MND_PercentChange %>%
    filter(year_start == "1990") %>%
    filter(year_end   == "2021") %>%
    filter(measure_name == labelx) %>%
    filter(sex_name   == "Both") %>%
    filter(age_name   == "Age-standardized") %>%
    filter(metric_name == "Rate") %>%
    mutate(
      val   = round(val   * 100, 1),   # multiply by 100 and keep 1 decimal place
      upper = round(upper * 100, 1),
      lower = round(lower * 100, 1)
    )
  
  # ---- Data for map drawing and CSV export ----
  dfplot <- dfx %>%
    filter(location_id %in% namex$location_id) %>%
    select(location_id, val) %>%
    left_join(., namex)
  
  dfplot_csv <- dfx %>%
    filter(location_id %in% namex$location_id) %>%
    select(location_id, val, upper, lower, year_start, year_end, measure_name) %>%
    left_join(., namex)
  
  write.xlsx(dfplot_csv, file = paste0(dn, "_204_PC_detail.xlsx"))
  
  # ---- Color palette (10 bins) ----
  color1 <- c("#67001f", "#b2182b", "#d6604d", "#f4a582", "#fddbc7",
              "#FFF8DC", "#d1e5f0", "#BBFFFF", "#4393c3", "#2166ac")
  
  # ---- Join with world geometry ----
  ASR <- dfplot %>%
    select(location_id, location, val) %>%
    mutate(val = val / 1)            # keep original scale (no change)
  
  df_asr <- left_join(df_world, ASR)
  
  # ---- Define bin breaks and labels (quantile-based), keep 1 decimal place ----
  xmin <- floor(min(na.omit(df_asr$val)) * 10) / 10
  xmax <- ceiling(max(na.omit(df_asr$val)) * 10) / 10
  
  breaks <- quantile(df_asr$val, probs = seq(0, 1, length.out = 11), na.rm = TRUE)
  breaks[1] <- xmin
  breaks[length(breaks)] <- xmax
  breaks <- round(breaks, 1)
  
  formatted_breaks <- vector("character", length(breaks) - 1)
  for (i in 1:(length(breaks) - 1)) {
    formatted_breaks[i] <- paste0(
      format(breaks[i], nsmall = 1), " - ",
      format(breaks[i + 1], nsmall = 1)
    )
  }
  
  # Safety check: ensure breaks and labels match
  if ((length(breaks) - 1) != length(formatted_breaks)) {
    stop("breaks and formatted_breaks length mismatch!")
  }
  
  # ---- Bin values and drop NA bins ----
  df_asr <- df_asr %>%
    mutate(asr = replace_na(val, -99)) %>%
    mutate(asr_cut = cut(asr, breaks = breaks, labels = formatted_breaks, include.lowest = TRUE)) %>%
    filter(!is.na(asr_cut))
  
  # ---- Main world map ----
  p <- ggplot(df_asr %>% na.omit()) +
    geom_sf(aes(geometry = geometry, fill = asr_cut), size = 0.1) +
    geom_sf(
      data  = filter(df_asr, asr == -99),
      aes(geometry = geometry),
      fill  = NA, color = "black", linetype = "dashed", size = 0.2
    ) +
    scale_fill_manual(values = rev(color1), guide = guide_legend(reverse = TRUE)) +
    guides(fill = guide_legend(ncol = 2, title = labelx))
  
  p1 <- p + theme(
    axis.text.x       = element_blank(),
    axis.ticks.x      = element_blank(),
    legend.position   = c(0.13, 0.29),
    legend.background = element_blank(),
    legend.key        = element_blank(),
    legend.title      = element_text(size = 8),
    legend.text       = element_text(size = 8),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    panel.border      = element_blank(),
    panel.background  = element_blank()
  )
  
  # ---- Sub-maps and layout (patchwork) ----
  library(patchwork)
  
  theme_map_sub <- theme_void() + labs(x = "", y = "") + theme_bw() +
    theme(
      text              = element_text(size = 9),
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank(),
      legend.position   = "none",
      axis.text         = element_blank(),
      axis.ticks        = element_blank(),
      plot.background   = element_rect(fill = "transparent"),
      plot.title        = element_text(vjust = 0.01, hjust = 0.5)
    )
  
  sub1 <- p1 + coord_sf(xlim = c(-92,  -59),  ylim = c(7,   28),  expand = FALSE) +
    ggtitle("Caribbean and central America") + theme(legend.position = "none") + theme_map_sub
  sub2 <- p1 + coord_sf(xlim = c(45,   55.8), ylim = c(21,  31.5), expand = FALSE) +
    ggtitle("Persian Gulf") + theme(legend.position = "none") + theme_map_sub
  sub3 <- p1 + coord_sf(xlim = c(12.5, 32),   ylim = c(34.5, 50),  expand = FALSE) +
    ggtitle("Balkan Peninsula") + theme(legend.position = "none") + theme_map_sub
  sub4 <- p1 + coord_sf(xlim = c(94.9, 119.1), ylim = c(-9.2, 9),  expand = FALSE) +
    ggtitle("Southeast Asia") + theme(legend.position = "none") + theme_map_sub
  sub5 <- p1 + coord_sf(xlim = c(-17.8, -7),  ylim = c(6.5, 15.8), expand = FALSE) +
    ggtitle("West Africa") + theme(legend.position = "none") + theme_map_sub
  sub6 <- p1 + coord_sf(xlim = c(30.5, 38.5), ylim = c(28.4, 35.9), expand = FALSE) +
    ggtitle("Eastern\n Mediterranean") + theme(legend.position = "none") + theme_map_sub
  sub7 <- p1 + coord_sf(xlim = c(4.7,  27.5), ylim = c(48,  59),   expand = FALSE) +
    ggtitle("Northern Europe") + theme(legend.position = "none") + theme_map_sub
  
  plot1 <- (sub1 + sub2 + sub3 + sub4) + plot_layout(nrow = 1)
  plot2 <- (sub5 | sub6) / sub7 + plot_layout(heights = c(2.8, 2.5))
  plot3 <- plot1 | plot2 + plot_layout(widths = c(1, 15))
  px    <- p1 / plot3 + plot_layout(heights = c(2, 1), widths = c(2, 1.2))
  
  ggsave(px, file = paste0(dn, "_PC_world.map.pdf"), height = 8, width = 14)
}


# ================================
# 2021 RATE (LEVEL) WORLD MAPS
# ================================
# Note: Here we map 2021 age-standardized rates (not PC). Values are typically small,
#       so keep original scale/precision (no rounding to 1 decimal place here).

setwd("GBD_analysis/Change 204 Countries")

for (dn in setdiff(unique(MND_PercentChange[["measure_name"]]), "Maternal mortality ratio")) {  # Exclude MMR
  labelx <- dn
  
  dfx <- MND %>%
    filter(year == "2021") %>%
    filter(measure_name == labelx) %>%
    filter(sex_name     == "Both") %>%
    filter(age_name     == "Age-standardized") %>%
    filter(metric_name  == "Rate")
  
  # ---- Data for map drawing and CSV export ----
  dfplot <- dfx %>%
    filter(location_id %in% namex$location_id) %>%
    select(location_id, val) %>%
    left_join(., namex)
  
  dfplot_csv <- dfx %>%
    filter(location_id %in% namex$location_id) %>%
    select(location_id, val, upper, lower, year, measure_name) %>%
    left_join(., namex)
  
  write.xlsx(dfplot_csv, file = paste0(dn, "_204_detail.xlsx"))
  
  # ---- Color palette (10 bins) ----
  color1 <- c("#67001f", "#b2182b", "#d6604d", "#f4a582", "#fddbc7",
              "#FFF8DC", "#d1e5f0", "#BBFFFF", "#4393c3", "#2166ac")
  
  # ---- Join with world geometry ----
  ASR <- dfplot %>%
    select(location_id, location, val) %>%
    mutate(val = val / 1)            # keep original scale (no change)
  
  df_asr <- left_join(df_world, ASR)
  
  # ---- Define bin breaks and labels (quantile-based), expand edges slightly ----
  xmin <- min(na.omit(df_asr$val)) - 0.0001
  xmax <- max(na.omit(df_asr$val)) + 0.0001
  
  breaks <- quantile(ASR$val, probs = seq(0, 1, length.out = 11), na.rm = TRUE)
  breaks[1] <- xmin
  breaks[length(breaks)] <- xmax
  
  formatted_breaks <- vector("character", length(breaks) - 1)
  for (i in 1:(length(breaks) - 1)) {
    formatted_breaks[i] <- paste(
      format(breaks[i],     digits = 4, nsmall = 3), "-",
      format(breaks[i + 1], digits = 4, nsmall = 3)
    )
  }
  
  # Safety check: ensure breaks and labels match
  if ((length(breaks) - 1) != length(formatted_breaks)) {
    stop("breaks and formatted_breaks length mismatch!")
  }
  
  # ---- Bin values and drop NA bins ----
  df_asr <- df_asr %>%
    mutate(asr = replace_na(val, -99)) %>%
    mutate(asr_cut = cut(asr, breaks = breaks, labels = formatted_breaks, include.lowest = TRUE)) %>%
    filter(!is.na(asr_cut))
  
  # ---- Main world map ----
  p <- ggplot(df_asr %>% na.omit()) +
    geom_sf(aes(geometry = geometry, fill = asr_cut), size = 0.1) +
    geom_sf(
      data  = filter(df_asr, asr == -99),
      aes(geometry = geometry),
      fill  = NA, color = "black", linetype = "dashed", size = 0.2
    ) +
    scale_fill_manual(values = rev(color1), guide = guide_legend(reverse = TRUE)) +
    guides(fill = guide_legend(ncol = 2, title = labelx))
  
  p1 <- p + theme(
    axis.text.x       = element_blank(),
    axis.ticks.x      = element_blank(),
    legend.position   = c(0.13, 0.29),
    legend.background = element_blank(),
    legend.key        = element_blank(),
    legend.title      = element_text(size = 8),
    legend.text       = element_text(size = 8),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    panel.border      = element_blank(),
    panel.background  = element_blank()
  )
  
  # ---- Sub-maps and layout (patchwork) ----
  library(patchwork)
  
  theme_map_sub <- theme_void() + labs(x = "", y = "") + theme_bw() +
    theme(
      text              = element_text(size = 9),
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank(),
      legend.position   = "none",
      axis.text         = element_blank(),
      axis.ticks        = element_blank(),
      plot.background   = element_rect(fill = "transparent"),
      plot.title        = element_text(vjust = 0.01, hjust = 0.5)
    )
  
  sub1 <- p1 + coord_sf(xlim = c(-92,  -59),  ylim = c(7,   28),  expand = FALSE) +
    ggtitle("Caribbean and central America") + theme(legend.position = "none") + theme_map_sub
  sub2 <- p1 + coord_sf(xlim = c(45,   55.8), ylim = c(21,  31.5), expand = FALSE) +
    ggtitle("Persian Gulf") + theme(legend.position = "none") + theme_map_sub
  sub3 <- p1 + coord_sf(xlim = c(12.5, 32),   ylim = c(34.5, 50),  expand = FALSE) +
    ggtitle("Balkan Peninsula") + theme(legend.position = "none") + theme_map_sub
  sub4 <- p1 + coord_sf(xlim = c(94.9, 119.1), ylim = c(-9.2, 9),  expand = FALSE) +
    ggtitle("Southeast Asia") + theme(legend.position = "none") + theme_map_sub
  sub5 <- p1 + coord_sf(xlim = c(-17.8, -7),  ylim = c(6.5, 15.8), expand = FALSE) +
    ggtitle("West Africa") + theme(legend.position = "none") + theme_map_sub
  sub6 <- p1 + coord_sf(xlim = c(30.5, 38.5), ylim = c(28.4, 35.9), expand = FALSE) +
    ggtitle("Eastern\n Mediterranean") + theme(legend.position = "none") + theme_map_sub
  sub7 <- p1 + coord_sf(xlim = c(4.7,  27.5), ylim = c(48,  59),   expand = FALSE) +
    ggtitle("Northern Europe") + theme(legend.position = "none") + theme_map_sub
  
  plot1 <- (sub1 + sub2 + sub3 + sub4) + plot_layout(nrow = 1)
  plot2 <- (sub5 | sub6) / sub7 + plot_layout(heights = c(2.8, 2.5))
  plot3 <- plot1 | plot2 + plot_layout(widths = c(1, 15))
  px    <- p1 / plot3 + plot_layout(heights = c(2, 1), widths = c(2, 1.2))
  
  ggsave(px, file = paste0(dn, "_world.map.pdf"), height = 8, width = 14)
}


# ================================
# EAPC WORLD MAPS
# ================================
# Note: EAPC values come from pre-computed results (MND_all_age.standarded.rds).
#       Here we only map EAPC (point estimate) by country using the same layout.

MND_eapce <- readRDS("MND_all_age.standarded.rds")

setwd("EAPC Change across 204 countries")

for (dn in unique(MND_eapce[["measure_name"]])) {  # Exclude MMR already handled earlier
  labelx <- dn
  
  dfx <- MND_eapce %>%
    filter(measure_name == labelx) %>%
    filter(sex_name     == "Both") %>%
    filter(age_name     == "Age-standardized") %>%
    filter(metric_name  == "Rate") %>%
    mutate(val = EAPC)   # map the EAPC point estimate
  
  # ---- Data for map drawing and CSV export ----
  dfplot <- dfx %>%
    filter(location_id %in% namex$location_id) %>%
    select(location_id, val) %>%
    left_join(., namex)
  
  dfplot_csv <- dfx %>%
    filter(location_id %in% namex$location_id) %>%
    select(location_id, EAPC, lower_CI, upper_CI, measure_name) %>%
    left_join(., namex)
  
  write.xlsx(dfplot_csv, file = paste0(dn, "_204_eapc.detail.xlsx"))
  
  # ---- Color palette (10 bins) ----
  color1 <- c("#67001f", "#b2182b", "#d6604d", "#f4a582", "#fddbc7",
              "#FFF8DC", "#d1e5f0", "#BBFFFF", "#4393c3", "#2166ac")
  
  # ---- Join with world geometry ----
  ASR <- dfplot %>%
    select(location_id, location, val) %>%
    mutate(val = val / 1)            # keep original scale (no change)
  
  df_asr <- left_join(df_world, ASR)
  
  # ---- Define bin breaks and labels (quantile-based), expand edges slightly ----
  xmin <- min(na.omit(df_asr$val)) - 0.0001
  xmax <- max(na.omit(df_asr$val)) + 0.0001
  
  breaks <- quantile(ASR$val, probs = seq(0, 1, length.out = 11), na.rm = TRUE)
  breaks[1] <- xmin
  breaks[length(breaks)] <- xmax
  
  formatted_breaks <- vector("character", length(breaks) - 1)
  for (i in 1:(length(breaks) - 1)) {
    formatted_breaks[i] <- paste(
      format(breaks[i],     digits = 4, nsmall = 3), "-",
      format(breaks[i + 1], digits = 4, nsmall = 3)
    )
  }
  
  # Safety check: ensure breaks and labels match
  if ((length(breaks) - 1) != length(formatted_breaks)) {
    stop("breaks and formatted_breaks length mismatch!")
  }
  
  # ---- Bin values and drop NA bins ----
  df_asr <- df_asr %>%
    mutate(asr = replace_na(val, -99)) %>%
    mutate(asr_cut = cut(asr, breaks = breaks, labels = formatted_breaks, include.lowest = TRUE)) %>%
    filter(!is.na(asr_cut))
  
  # ---- Main world map ----
  p <- ggplot(df_asr %>% na.omit()) +
    geom_sf(aes(geometry = geometry, fill = asr_cut), size = 0.1) +
    geom_sf(
      data  = filter(df_asr, asr == -99),
      aes(geometry = geometry),
      fill  = NA, color = "black", linetype = "dashed", size = 0.2
    ) +
    scale_fill_manual(values = rev(color1), guide = guide_legend(reverse = TRUE)) +
    guides(fill = guide_legend(ncol = 2, title = labelx))
  
  p1 <- p + theme(
    axis.text.x       = element_blank(),
    axis.ticks.x      = element_blank(),
    legend.position   = c(0.13, 0.29),
    legend.background = element_blank(),
    legend.key        = element_blank(),
    legend.title      = element_text(size = 8),
    legend.text       = element_text(size = 8),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    panel.border      = element_blank(),
    panel.background  = element_blank()
  )
  
  # ---- Sub-maps and layout (patchwork) ----
  library(patchwork)
  
  theme_map_sub <- theme_void() + labs(x = "", y = "") + theme_bw() +
    theme(
      text              = element_text(size = 9),
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank(),
      legend.position   = "none",
      axis.text         = element_blank(),
      axis.ticks        = element_blank(),
      plot.background   = element_rect(fill = "transparent"),
      plot.title        = element_text(vjust = 0.01, hjust = 0.5)
    )
  
  sub1 <- p1 + coord_sf(xlim = c(-92,  -59),  ylim = c(7,   28),  expand = FALSE) +
    ggtitle("Caribbean and central America") + theme(legend.position = "none") + theme_map_sub
  sub2 <- p1 + coord_sf(xlim = c(45,   55.8), ylim = c(21,  31.5), expand = FALSE) +
    ggtitle("Persian Gulf") + theme(legend.position = "none") + theme_map_sub
  sub3 <- p1 + coord_sf(xlim = c(12.5, 32),   ylim = c(34.5, 50),  expand = FALSE) +
    ggtitle("Balkan Peninsula") + theme(legend.position = "none") + theme_map_sub
  sub4 <- p1 + coord_sf(xlim = c(94.9, 119.1), ylim = c(-9.2, 9),  expand = FALSE) +
    ggtitle("Southeast Asia") + theme(legend.position = "none") + theme_map_sub
  sub5 <- p1 + coord_sf(xlim = c(-17.8, -7),  ylim = c(6.5, 15.8), expand = FALSE) +
    ggtitle("West Africa") + theme(legend.position = "none") + theme_map_sub
  sub6 <- p1 + coord_sf(xlim = c(30.5, 38.5), ylim = c(28.4, 35.9), expand = FALSE) +
    ggtitle("Eastern\n Mediterranean") + theme(legend.position = "none") + theme_map_sub
  sub7 <- p1 + coord_sf(xlim = c(4.7,  27.5), ylim = c(48,  59),   expand = FALSE) +
    ggtitle("Northern Europe") + theme(legend.position = "none") + theme_map_sub
  
  plot1 <- (sub1 + sub2 + sub3 + sub4) + plot_layout(nrow = 1)
  plot2 <- (sub5 | sub6) / sub7 + plot_layout(heights = c(2.8, 2.5))
  plot3 <- plot1 | plot2 + plot_layout(widths = c(1, 15))
  px    <- p1 / plot3 + plot_layout(heights = c(2, 1), widths = c(2, 1.2))
  
  ggsave(px, file = paste0(dn, "_world_eapc.map.pdf"), height = 8, width = 14)
}


#########################################
# Global Trends by Age and Sex (GBD-based)
#########################################

setwd("GBD_analysis")

library(tidyverse)
library(patchwork)
library(scales)

# Load helper objects (e.g., dfage with age ordering, etc.) and GBD-related data
load("function&GBDdata/GBD.Rdata")

# Source data
df <- MND
head(df, 3)

# Output directory for figures/tables
setwd("Global Trends by Age and Sex")

# Loop over measures, excluding "Probability of death"
for (i in setdiff(unique(df[["measure_name"]]), "Probability of death")) {
  
  # -----------------------------
  # 1) Filter to Global, 2021, sex-specific (exclude Both), MND, non-Percent, selected age groups
  # -----------------------------
  dfx <- df %>%
    filter(measure_name == i) %>%
    filter(location_id %in% c(1)) %>%                 # Global
    filter(year %in% c(2021)) %>%                     # Year 2021
    filter(!sex_name == "Both") %>%                   # Sex-specific
    filter(
      cause_name  == "Motor neuron disease",
      !metric_name == "Percent"
    ) %>%
    filter(age_id %in% c(1, 5:20, 30, 31, 32, 235))   # Selected age groups
  
  # Standardize age labels for 80+ (append "years" to the label)
  dfx <- dfx %>%
    mutate(age_name = if_else(age_name %in% c("80-84", "85-89", "90-94"),
                              paste(age_name, "years"), age_name))
  
  # If dfage exists from the loaded .Rdata (age ordering/id), keep consistent suffix
  dfage <- dfage %>%
    mutate(age_name = if_else(age_name %in% c("80-84", "85-89", "90-94"),
                              paste(age_name, "years"), age_name))
  
  # Merge age id-order mapping and ensure factor order of age_name
  dfid <- dfx %>%
    select(age_id, age_name) %>%
    distinct(age_id, .keep_all = TRUE) %>%
    left_join(., dfage) %>%
    arrange(id)
  
  dfx <- left_join(dfx, dfid) %>%
    filter(!is.na(id))
  
  dfx$age_name <- factor(dfx$age_name, levels = unique(dfid$age_name[order(dfid$id)]))
  
  df1 <- dfx
  
  # -----------------------------
  # 2) Split into Number and Rate; compute scaling so both can share a primary axis
  # -----------------------------
  df1_number <- df1 %>% filter(metric_name == "Number")
  df1_rate   <- df1 %>% filter(metric_name == "Rate")
  
  max_number     <- max(df1_number$val)
  max_rate       <- max(df1_rate$val)
  scaling_factor <- max_number / max_rate
  
  # -----------------------------
  # 3) Plot: bars = Number; lines + ribbons = Rate (scaled to primary axis)
  # -----------------------------
  ggplot() +
    # Bars: Number
    geom_bar(
      data = df1_number,
      aes(x = age_name, y = val, fill = sex_name),
      stat = "identity",
      position = position_dodge(),
      color = "black"
    ) +
    # Error bars: Number
    geom_errorbar(
      data = df1_number,
      aes(x = age_name, ymin = lower, ymax = upper, group = sex_name),
      position = position_dodge(0.9),
      width = 0.25
    ) +
    # Lines: Rate (scaled)
    geom_line(
      data = df1_rate,
      aes(x = age_name, y = val * scaling_factor, group = sex_name, color = sex_name),
      position = position_dodge(0.9)
    ) +
    # Ribbons: Rate 95% UI (scaled, semi-transparent)
    geom_ribbon(
      data = df1_rate,
      aes(x = age_name, ymin = lower * scaling_factor, ymax = upper * scaling_factor,
          group = sex_name, fill = sex_name),
      alpha = 0.5,
      position = position_dodge(width = 0.9)
    ) +
    # Labels
    labs(x = "", y = i) +
    # Primary y-axis for Number; secondary axis back-transforms to Rate scale
    scale_y_continuous(
      labels = label_number(unit = "K"),
      sec.axis = sec_axis(~ . / scaling_factor, name = "Rate per 100,000 population")
    ) +
    # Colors for bars/lines
    scale_fill_manual(values = c("Male" = "skyblue", "Female" = "lightpink"), name = "") +
    scale_color_manual(values = c("Male" = "blue", "Female" = "red"), name = "") +
    theme_classic() +
    theme(
      legend.position = "right",
      axis.text.x = element_text(angle = 65, hjust = 1)
    )
  
  # -----------------------------
  # 4) Save figure and export plotting data
  # -----------------------------
  ggsave(paste0(i, "_change.with.age_sex.numerber_rate.pdf"), height = 6, width = 8)
  
  df_001 <- rbind(df1_number, df1_rate)
  df_001 <- df_001 %>% select(-id)
  
  write.csv(df_001, file = paste0(i, "_change.with.age_sex.numerber_rate.csv"), row.names = FALSE)
}



#############################################
# Correlation between SDI and MND Metrics
# (21 GBD regions + Global, Age-standardized Rate, SDI-2019)
#############################################

setwd("GBD_analysis")

library(tidyverse)
library(patchwork)
library(ggsci)
library(scales)
library(openxlsx)

load("function&GBDdata/GBD.Rdata")  # provides: MND, SDI2019, name_21, etc.

# Load main data using qs (compressed format)
library(qs)
MND <- qread("MND.qs")
df <- MND
head(df,3)
setwd("Correlation between SDI and MND Metrics")

for (i in setdiff(unique(df$'measure_name'),'Probability of death')) {
  dfx = df %>% filter(measure_name == i) %>% 
    filter(location_id %in% c(1, name_21)) %>% 
    filter(sex_name == "Both") %>% 
    filter(metric_name == 'Rate') %>%
    filter(age_id %in% c(27))
  head(dfx)
  
  # Define the specific order for the location names
  location_order <- c(
    "Global", "High-income Asia Pacific", "High-income North America", "Western Europe",
    "Australasia", "Andean Latin America", "Tropical Latin America", "Central Latin America",
    "Southern Latin America", "Caribbean", "Central Europe", "Eastern Europe", "Central Asia",
    "North Africa and Middle East", "South Asia", "Southeast Asia", "East Asia", "Oceania",
    "Western Sub-Saharan Africa", "Eastern Sub-Saharan Africa", "Central Sub-Saharan Africa",
    "Southern Sub-Saharan Africa"
  )
  
  library(ggsci)
  library(openxlsx)  # <-- added for Excel export
  colors <- c(
    pal_npg("nrc", alpha = 0.7)(9),
    pal_aaas("default", alpha = 0.7)(9),
    pal_nejm("default", alpha = 0.7)(8),
    pal_jama("default", alpha = 0.7)(7)
  )
  
  ######### 1.1  Incidence SDI 2019 in 21
  ############### ############### ############### ############### ############### 
  # bin data
  DALY_2017 = dfx %>% select(location_id, val, year)
  df11 = left_join(DALY_2017, SDI2019)
  # get df
  df3 = df11 %>% filter(sdi > 0.2)
  
  # plot
  # Filter out rows with NA in 'sdi' or 'val'
  df11_filtered <- df11 %>% 
    filter(!is.na(sdi) & !is.na(val)) %>%
    mutate(location_name = factor(location_name, levels = location_order))
  
  # --- NEW: export the raw plotting data to Excel ---
  write.xlsx(
    df11_filtered %>% mutate(location_name = as.character(location_name)),
    file = paste0(i, "_21_sdi_raw.xlsx"),
    overwrite = TRUE
  )
  
  # Calculate Spearman's correlation
  spearman_cor <- cor.test(df11_filtered$sdi, df11_filtered$val, method = "spearman")
  # Format r and p values for display
  r <- spearman_cor$estimate
  p <- as.numeric(format(spearman_cor$p.value, scientific = TRUE))
  spearx <- sprintf("r=%.4f, p=%.4e", r, p)
  
  p1 <- ggplot(df11_filtered, aes(x = sdi, y = val, color = location_name, shape = location_name)) +
    geom_point() +  # Plot points
    geom_smooth(method = "loess", se = TRUE, aes(group = 1), color = "#92A8D1") +  # Add a smoothing line
    scale_shape_manual(values = 1:22, breaks = location_order, labels = location_order) +  # Use manually defined shapes
    scale_color_manual(values = colors, breaks = location_order, labels = location_order) +
    labs(x = "SDI", 
         shape = "", color = "", x = paste0("SDI"),
         y = paste0("Age-standardised ", i, " Rate (per 100,000 population)")) +
    theme_bw() +
    theme(
      legend.key.size = unit(0.03, "line"),
      panel.border = element_blank(), 
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(), 
      legend.background = element_blank(),
      legend.title = element_text(size = 12),
      axis.line = element_line(colour = "black"),
      legend.text = element_text(size = 8)
    ) +
    guides(shape = guide_legend(nrow = 22)) +
    # 在图表左上角添加 r 和 p 值
    annotate("text", x = 0.2, y = max(df11_filtered$val), label = spearx, hjust = 0, size = 4, color = "black")
  
  ggsave(p1, file = paste0(i, "_21.sdi.pdf"), height = 9, width = 12)
}





