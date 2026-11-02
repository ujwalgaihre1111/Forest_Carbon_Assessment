inventory <- read.csv("./data/bagmati/bagmati_raw_inventory.csv")

1. INPUT AND OUTPUT ------------------------------------------
file_path <- "data/bagmati/bagmati_raw_inventory.csv"
if (!file.exists(file_path)) file_path <- "bagmati_raw_inventory.csv"
if (!file.exists(file_path)) stop("CSV not found. Check the file location!")
dir.create("outputs", showWarnings = FALSE, recursive = TRUE)

inventory <- read.csv(file_path, stringsAsFactors = FALSE,
                      na.strings = c("", "NA", "N/A", "NULL"))
required <- c("Plot_ID", "Stand_ID", "Tree_ID", "Species", "DBH_cm",
              "Height_m", "Wood_density", "Biomass_kg", "Plot_area_ha",
              "Elevation_m", "Latitude", "Longitude",
              "Measurement_Date", "Census_Year")
if (!all(required %in% names(inventory))) {
  stop("Missing columns: ", paste(setdiff(required, names(inventory)), collapse = ", "))
}
if (nrow(inventory) == 0) stop("The dataset is empty.")
inventory$CSV_Row <- seq_len(nrow(inventory)) + 1L

# 2. DATA EXPLORATION ------------------------------------------
cat("\n=== DATA EXPLORATION ===\n")
cat("Rows:", nrow(inventory), " Columns:", length(required), "\n")
cat("Stands:", length(unique(inventory$Stand_ID)), "\n")
cat("Unique stand-plot pairs:",
    nrow(unique(inventory[c("Stand_ID", "Plot_ID")])), "\n")

str(inventory)
print(head(inventory, 6))
print(summary(inventory))
missing <- sapply(inventory[required], function(x) sum(is.na(x)))
print(missing)

write.csv(data.frame(Variable = names(missing), Missing = as.integer(missing)),
          "outputs/missing_values.csv", row.names = FALSE)
writeLines(capture.output(str(inventory), summary(inventory),
                          table(inventory$Census_Year, useNA = "ifany")),
           "outputs/data_exploration.txt")
# Convert numeric variables safely; nonnumeric values become NA and get flagged.
numeric_vars <- c("DBH_cm", "Height_m", "Wood_density", "Biomass_kg",
                  "Plot_area_ha", "Elevation_m", "Latitude", "Longitude",
                  "Census_Year")
for (v in numeric_vars) inventory[[v]] <- suppressWarnings(as.numeric(inventory[[v]]))


# 3. QUALITY CONTROL -------------------------------------------
# One record may have multiple flags; the audit file keeps each reason.
flags <- list()
add_flag <- function(test, severity, issue) {
  rows <- which(!is.na(test) & test)
  if (length(rows) > 0) {
    flags[[length(flags) + 1L]] <<- data.frame(
      CSV_Row = inventory$CSV_Row[rows],
      Stand_ID = inventory$Stand_ID[rows],
      Plot_ID = inventory$Plot_ID[rows],
      Tree_ID = inventory$Tree_ID[rows],
      Severity = severity, Issue = issue
    )
  }
}

invalid_text <- function(x) is.na(x) | trimws(as.character(x)) == ""
bad_id <- invalid_text(inventory$Stand_ID) |
          invalid_text(inventory$Plot_ID) |
          invalid_text(inventory$Tree_ID)

bad_biomass <- !is.finite(inventory$Biomass_kg) | inventory$Biomass_kg <= 0
bad_area <- !is.finite(inventory$Plot_area_ha) | inventory$Plot_area_ha <= 0
add_flag(bad_id, "CRITICAL", "Missing stand, plot, or tree ID: exclude")
add_flag(bad_biomass, "CRITICAL", "Biomass missing/nonpositive: exclude")
add_flag(bad_area, "CRITICAL", "Plot area missing/nonpositive: exclude")

# Screening thresholds are REVIEW triggers, not universal biological limits.
add_flag(!is.finite(inventory$DBH_cm) | inventory$DBH_cm <= 0,
         "REVIEW", "DBH missing, zero or negative; verify field sheet")
add_flag(is.finite(inventory$DBH_cm) & inventory$DBH_cm > 150,
         "REVIEW", "DBH > 150 cm; verify field sheet")
add_flag(!is.finite(inventory$Height_m) | inventory$Height_m <= 0,
         "REVIEW", "Height missing or nonpositive; verify field sheet")
add_flag(!is.finite(inventory$Wood_density) |
         inventory$Wood_density < 0.2 | inventory$Wood_density > 1.2,
         "REVIEW", "Wood density outside 0.2-1.2 g/cm3")
add_flag(is.finite(inventory$Biomass_kg) & inventory$Biomass_kg > 10000,
         "REVIEW", "Biomass > 10000 kg; check tree and biomass model")
add_flag(invalid_text(inventory$Species), "REVIEW", "Missing species name")

# Tree identity: same tree number can occur in DIFFERENT plots.
key <- paste(inventory$Stand_ID, inventory$Plot_ID,
             inventory$Census_Year, inventory$Tree_ID, sep = "|")
exact_duplicate <- duplicated(inventory[required])
key_groups <- split(seq_len(nrow(inventory)), key)
conflict <- rep(FALSE, nrow(inventory))
for (idx in key_groups) {
  if (length(idx) > 1 &&
      nrow(unique(inventory[idx, required, drop = FALSE])) > 1) {
    conflict[idx] <- TRUE
  }
}
add_flag(exact_duplicate, "CRITICAL", "Exact duplicate row: keep first")
add_flag(conflict, "CRITICAL", "Conflicting entries for same tree: exclude all")


# A single plot must have a consistent area within each census.
plot_key <- paste(inventory$Stand_ID, inventory$Plot_ID,
                  inventory$Census_Year, sep = "|")
area_conflict <- rep(FALSE, nrow(inventory))
for (idx in split(seq_len(nrow(inventory)), plot_key)) {
  a <- unique(inventory$Plot_area_ha[idx])
  a <- a[is.finite(a) & a > 0]
  if (length(a) > 1) area_conflict[idx] <- TRUE
}
add_flag(area_conflict, "CRITICAL", "Conflicting area within plot: exclude plot")


# Do not mix partial remeasurement (2026) with complete baseline (2025).
# Select automatically the year containing the most distinct sampled plots.
years <- sort(unique(inventory$Census_Year[is.finite(inventory$Census_Year)]))
if (!length(years)) stop("No valid census year!")
plot_coverage <- sapply(years, function(y) {
  d <- inventory[inventory$Census_Year == y &
                 !is.na(inventory$Census_Year), , drop = FALSE]
  nrow(unique(d[c("Stand_ID", "Plot_ID")]))
})
selected_year <- years[which.max(plot_coverage)]
other_year <- is.na(inventory$Census_Year) |
              inventory$Census_Year != selected_year
add_flag(other_year, "EXCLUDE", "Different/unknown census year: analyze separately")

keep <- !(bad_id | bad_biomass | bad_area | exact_duplicate |
          conflict | area_conflict | other_year)
keep[is.na(keep)] <- FALSE
clean <- inventory[keep, , drop = FALSE]
if (!nrow(clean)) stop("No observations remain after QC.")

# Species-name normalization for richness only; preserve original Species.
clean$Species_clean <- trimws(clean$Species)
aliases <- c("s. robusta" = "Shorea robusta",
             "shorea rubusta" = "Shorea robusta",
             "schima wallichi" = "Schima wallichii",
             "c. indica" = "Castanopsis indica")
found <- match(tolower(clean$Species_clean), names(aliases))
clean$Species_clean[!is.na(found)] <- aliases[found[!is.na(found)]]
clean$Species_clean[invalid_text(clean$Species_clean)] <- NA_character_

qc_log <- if (length(flags)) do.call(rbind, flags) else data.frame()
write.csv(qc_log, "outputs/quality_control_flags.csv", row.names = FALSE)
write.csv(data.frame(CSV_Row = inventory$CSV_Row, Included = keep,
                     Stand_ID = inventory$Stand_ID, Plot_ID = inventory$Plot_ID,
                     Tree_ID = inventory$Tree_ID),
          "outputs/row_inclusion_audit.csv", row.names = FALSE)


# 4. CUSTOM FUNCTION: STAND-LEVEL CARBON -----------------------
assess_stand <- function(d, carbon_fraction = 0.47,
                         min_trees = 2, min_plots = 2) {
  needed <- c("Stand_ID", "Plot_ID", "Tree_ID", "Biomass_kg",
              "Plot_area_ha", "Species_clean")
  if (!is.data.frame(d) || !all(needed %in% names(d)))
    stop("Invalid stand data frame or missing variables.")
  if (length(unique(d$Stand_ID)) != 1) stop("Input must contain one stand.")
  if (!is.numeric(carbon_fraction) || length(carbon_fraction) != 1 ||
      !is.finite(carbon_fraction) || carbon_fraction <= 0 ||
      carbon_fraction > 1) stop("Invalid carbon fraction.")
  if (nrow(d) < min_trees) stop("Insufficient trees in stand.")
  if (any(!is.finite(d$Biomass_kg) | d$Biomass_kg <= 0))
    stop("Missing/invalid tree biomass: run QC first.")
  if (any(!is.finite(d$Plot_area_ha) | d$Plot_area_ha <= 0))
    stop("Invalid plot area: run QC first.")
  review <- (!is.finite(d$DBH_cm) | d$DBH_cm <= 0 | d$DBH_cm > 150 |
             !is.finite(d$Height_m) | d$Height_m <= 0 |
             !is.finite(d$Wood_density) | d$Wood_density < 0.2 |
             d$Wood_density > 1.2)
  if (any(review)) message("Stand ", d$Stand_ID[1], ": ", sum(review),
                           " trees need biological measurement review.")

  plots <- split(d, d$Plot_ID)
  if (length(plots) < min_plots) stop("Insufficient plots in stand.")
  areas <- sapply(plots, function(p) {
    a <- unique(p$Plot_area_ha)
    if (length(a) != 1 || !is.finite(a) || a <= 0)
      stop("Inconsistent plot area: ", p$Plot_ID[1])
    a
  })
  sampled_area <- sum(areas)   # IMPORTANT: count each plot ONCE
  d$Carbon_kg <- d$Biomass_kg * carbon_fraction
  d$Carbon_Mg <- d$Carbon_kg / 1000
  total_carbon <- sum(d$Carbon_Mg)

  plot_result <- do.call(rbind, lapply(plots, function(p) {
    c_mg <- sum(p$Biomass_kg * carbon_fraction / 1000)
    data.frame(Stand_ID = p$Stand_ID[1], Plot_ID = p$Plot_ID[1],
               Trees = nrow(p), Plot_area_ha = unique(p$Plot_area_ha),
               Carbon_Mg = c_mg,
               Carbon_density_Mg_C_ha = c_mg / unique(p$Plot_area_ha))
  }))
  species <- unique(na.omit(d$Species_clean))
  stand_result <- data.frame(
    Stand_ID = d$Stand_ID[1], Total_trees = nrow(d),
    Number_of_plots = length(plots), Species_richness = length(species),
    Max_tree_biomass_kg = max(d$Biomass_kg),
    Total_carbon_Mg_C = total_carbon,
    Total_sampled_area_ha = sampled_area,
    Carbon_density_Mg_C_ha = total_carbon / sampled_area
  )
  list(stand = stand_result, trees = d, plots = plot_result)
}

# 5. DYNAMIC PARTITIONING: split(), sapply(), by() --------------
stand_list <- split(clean, clean$Stand_ID)
cat("\nTREES BY STAND (sapply):\n")
print(sapply(stand_list, nrow))
results <- lapply(stand_list, assess_stand)
stand_summary <- do.call(rbind, lapply(results, function(x) x$stand))
tree_results <- do.call(rbind, lapply(results, function(x) x$trees))
plot_results <- do.call(rbind, lapply(results, function(x) x$plots))
rownames(stand_summary) <- NULL

# Independent validation with by().
biomass_by_stand <- by(clean$Biomass_kg, clean$Stand_ID, sum)
check_carbon <- as.numeric(unlist(biomass_by_stand)[stand_summary$Stand_ID]) * 0.47 / 1000
stopifnot(isTRUE(all.equal(check_carbon, stand_summary$Total_carbon_Mg_C,
                           check.attributes = FALSE)))


# 6. FINAL MATRIX: STANDS = COLUMNS; METRICS = ROWS ------------
metrics <- c("Total_trees", "Species_richness", "Max_tree_biomass_kg",
             "Total_carbon_Mg_C", "Total_sampled_area_ha",
             "Carbon_density_Mg_C_ha")
summary_matrix <- t(as.matrix(stand_summary[metrics]))
colnames(summary_matrix) <- stand_summary$Stand_ID
rownames(summary_matrix) <- c("Total trees", "Species richness",
                              "Maximum tree biomass (kg)",
                              "Total carbon (Mg C)", "Sampled area (ha)",
                              "Carbon density (Mg C/ha)")
print(round(summary_matrix, 3))
write.csv(summary_matrix, "outputs/stand_summary_matrix.csv")
write.csv(stand_summary, "outputs/stand_summary.csv", row.names = FALSE)
write.csv(tree_results, "outputs/tree_level_carbon.csv", row.names = FALSE)
write.csv(plot_results, "outputs/plot_level_carbon.csv", row.names = FALSE)


# 7. QC REPORT -------------------------------------------------
report <- c(
  "BAGMATI FOREST CARBON: QUALITY CONTROL REPORT",
  paste("Original records:", nrow(inventory)),
  paste("Selected census year:", selected_year),
  paste("Records included:", sum(keep)),
  paste("Records excluded:", sum(!keep)),
  paste("Exact duplicates:", sum(exact_duplicate)),
  paste("Conflicting tree identities:", sum(conflict)),
  paste("Conflicting plot areas:", sum(area_conflict)),
  paste("Invalid biomass:", sum(bad_biomass)),
  paste("Invalid area:", sum(bad_area)),
  paste("Other-year records:", sum(other_year)),
  "CRITICAL: exclude missing IDs, invalid biomass/area, duplicated or",
  "conflicting records. REVIEW: keep valid supplied biomass but verify",
  "unusual DBH, height, density, and biomass against field records.",
  "Species aliases standardized for richness; raw names preserved.",
  "Plot area summed once per unique plot (never once per tree).",
  "See quality_control_flags.csv and row_inclusion_audit.csv."
)
writeLines(report, "outputs/quality_control_report.txt")


# 8. SCIENTIFIC INTERPRETATION (about 330 words) ---------------
ranked <- stand_summary[order(-stand_summary$Carbon_density_Mg_C_ha), ]
highest <- ranked[1, ]
lowest <- ranked[nrow(ranked), ]
ranking <- paste(sprintf("%s (%.2f Mg C/ha)", ranked$Stand_ID,
                         ranked$Carbon_density_Mg_C_ha), collapse = "; ")
interpretation <- c(
  "# Scientific interpretation: Forest above-ground carbon in Bagmati",
  "",
  sprintf(paste0("The inventory included %d tree observations across %d ",
    "forest stands and %d sampled plots after quality-control screening. ",
    "The most spatially complete census year (%d) was used, while records ",
    "from other years were kept separate to prevent double counting. ",
    "Above-ground carbon was estimated by multiplying recorded biomass ",
    "by 0.47 and converting kilograms to megagrams. Stand-level carbon ",
    "density was obtained by dividing summed stand carbon by the total ",
    "sampled plot area, with each plot counted only once."),
    nrow(clean), nrow(ranked), sum(ranked$Number_of_plots), selected_year),
  "",
  sprintf(paste0("Carbon density ranked from highest to lowest as follows: ",
    "%s. %s had the highest recorded density (%.2f Mg C/ha), whereas ",
    "%s had the lowest (%.2f Mg C/ha). The difference between them was ",
    "%.2f Mg C/ha. This comparison describes the carbon contained in ",
    "the measured above-ground biomass within the sampled plots, rather ",
    "than the total forest carbon stock across the entire landscape."),
    ranking, highest$Stand_ID, highest$Carbon_density_Mg_C_ha,
    lowest$Stand_ID, lowest$Carbon_density_Mg_C_ha,
    highest$Carbon_density_Mg_C_ha - lowest$Carbon_density_Mg_C_ha),
  "",
  paste0("Several ecological factors could contribute to variation in ",
    "carbon density, including differences in tree diameter, height, ",
    "species composition, stand development, and the concentration ",
    "of large-biomass trees. However, this descriptive analysis does ",
    "not test those factors statistically and therefore cannot assign ",
    "causal explanations to the observed differences. Species names ",
    "were standardized before richness was calculated, while the ",
    "original species labels were retained for traceability."),
  "",
  sprintf(paste0("Quality-control screening began with %d original records ",
    "and excluded %d from the selected single-year assessment. Exact ",
    "duplicate records were removed, and observations from a separate ",
    "census were not pooled with the baseline. Missing or suspicious ",
    "diameters, biologically unusual measurements, and other warnings ",
    "were recorded for field verification. Because above-ground biomass ",
    "was supplied directly in the dataset, questionable diameter or ",
    "height values did not automatically invalidate an otherwise ",
    "positive biomass value. Every inclusion decision can be traced ",
    "using the accompanying quality-control audit files."),
    nrow(inventory), sum(!keep)),
  "",
  paste0("The reported estimates cover living-tree above-ground carbon ",
    "only. Soil carbon, below-ground roots, deadwood, and litter are not ",
    "included. The accuracy of the estimates depends on the reliability ",
    "of the supplied biomass calculations and the representativeness ",
    "of the field plots. Future assessments should verify unusual ",
    "records against original field sheets, document the biomass ",
    "equations used, and quantify plot-level sampling uncertainty ",
    "before extrapolating carbon density to broader forest areas.")
)
writeLines(interpretation, "outputs/scientific_interpretation.md")
cat("\nDONE! Files saved in the outputs folder.\n")
