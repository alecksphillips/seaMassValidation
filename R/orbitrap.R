#' Orbitrap
#'
#' Blah.
#'
#' @import data.table
#' @export
orbitrap <- function() {
  library(MSqRob)
  library(MSnbase)
  library(seaMass)
  library(ggplot2)
  library(data.table)
  setwd("~/gonzales1/seaMass/unwin_orbitrap/")
  sm.path <- "~/bc4/Scratch/seaMass/"

  ### Split out 4v4
  assays.rm <- setdiff(
    c("A1-5", "A1-6", "A1-7", "A2-4", "A3-5", "A3-6", "A4-5", "A4-6", "A5-3", "A6-8", "B1-3", "B2-4", "B3-5", "B3-6", "B4-5", "B4-6", "B5-3", "B6-3", "C1-3", "C2-3", "C3-3", "C4-3"),
    c("A1-6", "A3-6", "A5-3", "A6-8", "B2-4", "B3-5", "B4-5", "B5-3", "C1-3", "C2-3", "C3-3", "C4-3")
  )

  DT.raw <- fread("input/peptides.txt")
  for (assay in assays.rm) DT.raw[, (grep(paste0("^.* ", assay), colnames(DT.raw))) := NULL]
  fwrite(DT.raw, "input4v4/peptides.txt", sep = "\t", quote = F)

  DT.raw <- fread("input/proteinGroups.txt")
  for (assay in assays.rm) DT.raw[, (grep(paste0("^.* ", assay), colnames(DT.raw))) := NULL]
  fwrite(DT.raw, "input4v4/proteinGroups.txt", sep = "\t", quote = F)


  ### Import MaxQuant data with MSqRob's own protein grouping ######################################################

  # Based on CPTAC part of: https://htmlpreview.github.io/?https://github.com/statOmics/MSqRob/blob/master/vignettes/MSqRob.html

  # function to import how MSqRob likes it, but with our own ground truth median normalisation
  import_with_MSqRob_protein_groups <- function(minIdentified = 2, normalise = T) {
    peptides <- preprocess_MaxQuant(
      import2MSnSet("input4v4/peptides.txt", filetype = "MaxQuant", remove_pattern = T),
      exp_annotation = data.frame(
        run = c("A1-6", "A3-6", "A5-3", "A6-8", "B2-4", "B3-5", "B4-5", "B5-3"),
        condition = c("A", "A", "A", "A", "B", "B", "B", "B")
      ),
      useful_properties = c("Proteins", "Sequence", "PEP", "Evidence.IDs"), # add Evidence IDs for seaMass three-level model
      normalisation = "none",
      remove_only_site = T,
      file_proteinGroups = "input4v4/proteinGroups.txt",
      minIdentified = minIdentified
    )

    # mean-centred median normalisation on ground truth rat proteins only
    data <- exprs(peptides)[grep("_RAT", fData(peptides)$Proteins),]
    data <- data[complete.cases(data),]
    expos <- apply(data - rowMeans(data, na.rm = T), 2, median)
    exprs(peptides) <- sweep(exprs(peptides), 2, expos)

    return(peptides)
  }

  # input data for MSqRob and seaMass
  peptides2 <- import_with_MSqRob_protein_groups() # MSqRob default
  peptides1 <- import_with_MSqRob_protein_groups(minIdentified = 1) # Run MSqRob with 'minIdentified = 1'
  peptides <- import_with_MSqRob_protein_groups(minIdentified = 1, normalise = F) # For seaMass three-level model don't need to normalise and use 'minIdentified = 1'

  ### MSQROB ########################################################################

  # function to execute MSqRob
  execute_MSqRob <- function(peptides, name, shrinkage = TRUE) {
    if (!file.exists(paste(name, "rds", sep = "."))) {
      # Convert data to a protdata object
      proteins <- MSnSet2protdata(peptides, accession = "Proteins")

      # Fit MSqRob model
      if (shrinkage) {
        models <- fit.model(
          proteins,
          response = "quant_value",
          fixed = "condition",
          random = c("Sequence","run")
        )
      } else {
        models <- fit.model(
          proteins,
          response = "quant_value",
          fixed = "condition",
          shrinkage.fixed = c(0, 0),
          random = c("Sequence","run")
        )
      }

      # Test the appropriate research hypotheses
      results <- test.contrast_adjust(models, makeContrast(
        contrasts = c("conditionB - conditionA"),
        levels = c("conditionA", "conditionB")
      ))

      results <- cbind(Group = rownames(results), results)
      rownames(results) <- NULL
      colnames(results) <- c("Group", "m", "s", "df", "tvalue", "pvalue", "qvalue", "signif")
      results <- add_seaMass_spikein_truth(results)

      saveRDS(results, paste(name, "rds", sep = "."))
    } else {
      results <- readRDS(paste(name, "rds", sep = "."))
    }

    g <- plot_volcano(results, x.col = "m", y.col = "qvalue", stdev.col = "s")
    ggsave(paste0(name, ".pdf"))

    return(g)
  }

  # Run MSqRob with default settings
  execute_MSqRob(peptides2, "QR_uo4_d2_nMG__qr_sT")
  execute_MSqRob(peptides2, "QR_uo4_d2_nMG__qr_sF", shrinkage = F)

  # Run MSqRob with 'minIdentified = 1'
  execute_MSqRob(peptides1, "QR_uo4_d1_nMG__qr_sT")
  execute_MSqRob(peptides1, "QR_uo4_d1_nMG__qr_sF", shrinkage = F)


  plot_pr(list(
    "MSqRob d=2 n=mt -> MSqRob s=T" = readRDS("QR_uo4_d2_nMG__qr_sT.rds"),
    "MSqRob d=2 n=mt -> MSqRob s=F" = readRDS("QR_uo4_d2_nMG__qr_sF.rds"),
    "MSqRob d=1 n=mt -> MSqRob s=T" = readRDS("QR_uo4_d1_nMG__qr_sT.rds"),
    "MSqRob d=1 n=mt -> MSqRob s=F" = readRDS("QR_uo4_d1_nMG__qr_sF.rds")
  ))

  plot_pr(list(
    "MSqRob d=2 n=mt -> MSqRob s=T" = readRDS("QR_uo4_d2_nMG__qr_sT.rds"),
  ))


  plot_pr(list(
    "MSqRob d=2 n=mt -> MSqRob s=T" = readRDS("uo__qr_d2_nMG__qr_sT.rds")[["conditionB - conditionA"]],
    "MSqRob d=2 n=mt -> MSqRob s=F" = readRDS("uo__qr_d2_nMG__qr_sF.rds")[["conditionB - conditionA"]],
    #"MSqRob d=1 n=mt -> MSqRob s=T" = readRDS("uo__qr_d1_nMG__qr_sF.rds")[["conditionB - conditionA"]],
    "MSqRob d=1 n=mt -> MSqRob s=F" = readRDS("uo__qr_d1_nMG__qr_sF.rds")[["conditionB - conditionA"]]
  ), plot.fdr = F, legend.nrow = 1)
  ggplot2::ggsave("msqrob__msqrob.pdf", width = 12, height = 6)
















  ### SEAMASS ########################################################################

  # function to execute seaMass
  execute_seaMass <- function(data, name, error.models, missingness.models, two.level = FALSE) {
    data.design <- new_assay_design(data)
    data.design$Condition <- substr(data.design$Sample, 1, 1)

    for (i in 1:length(error.models)) {
      #message(i)
      for (j in 1:length(missingness.models)) {
        #message(j)
        path <- paste0(sm.path, name, "__sm", ifelse(two.level, "2", "3"), "_e", names(error.models)[i], "_m", names(missingness.models)[j])

        schedule <- schedule_slurm(
          submit.prefix = "ssh bc4login.acrc.bris.ac.uk /mnt/storage/home/ad16243/Scratch/seaMass/",
          cpus_per_task = 28,
          mem = "128000m", # unfortunately this dataset takes more than 64Gb!
          mail_user = "ad16243@bristol.ac.uk"
        )

        if (two.level) {
          sigma <- seaMass_sigma(data, data.design, path = path, run = F, control = sigma_control(
            keep = "raw.group.quants", error.model = error.models[i], missingness.model = missingness.models[j], schedule = schedule, nthread = 7,
            component.model = NULL, component.eb.min = 1, assay.model = "measurement"
          ))
        } else {
          sigma <- seaMass_sigma(data, data.design, path = path, run = F, control = sigma_control(
            keep = "raw.group.quants", error.model = error.models[i], missingness.model = missingness.models[j], schedule = schedule, nthread = 7
          ))
        }


        norm.all <- c("M" = "median", "Q" = "quantile", "T" = "theta")
        norm.groups.all <- c("A" = ".*", "G" = "_RAT")
        for (k in 1:length(norm.all)) {
          for (l in 1:length(norm.groups.all)) {

            norm <- paste0("n", names(norm.all)[k], names(norm.groups.all)[l])

            #seaMass_delta(
            #  sigma,
            #  name = paste0(norm, "_tW"),
            #  norm.groups = norm.groups.all[l],
            #  control = delta_control(
            #    keep = "group.fdr",
            #    norm.model = norm.all[k]
            #  )
            #)

            seaMass_delta(
              sigma,
              name = paste0(norm, "_tS"),
              norm.groups = norm.groups.all[l],
              control = delta_control(
                keep = "group.fdr",
                norm.model = norm.all[k]
              ),
              random = ~ Sample,
              prior = list(
                G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 25^2)),
                R = list(V = 1, nu = 2e-4)
              )
            )

          }
        }

        run(sigma)
      }
    }
  }

  # Run seaMass three-level model with 'minIdentified = 1'
  data <- import_MSqRob(fData(peptides), exprs(peptides), "input/evidence.txt")
  execute_seaMass(data, "uo__qr_d1_n0", c("L" = "lognormal", "P" = "poisson"), c("0" = "", "r" = "rm", "1" = "one", "M" = "minimum", "C3" = "censored3", "C" = "censored", "C5" = "censored5"))

  # Run seaMass two-level model
  data <- import_MSqRob(fData(peptides2), exprs(peptides2))
  execute_seaMass(data, "uo__qr_d2_n0", c("L" = "lognormal", "P" = "poisson"), c("0" = "", "r" = "rm", "1" = "one", "M" = "minimum", "C" = "censored"), two.level = T)

  # Run seaMass two-level model with 'minIdentified = 1'
  data <- import_MSqRob(fData(peptides1), exprs(peptides1))
  execute_seaMass(data, "uo__qr_d1_n0", c("L" = "lognormal", "P" = "poisson"), c("0" = "", "r" = "rm", "1" = "one", "M" = "minimum", "C" = "censored"), two.level = T)

  # Run seaMass three-level model with MaxQuant's protein groups
  data <- import_MaxQuant("input/proteinGroups.txt", "input/evidence.txt")
  execute_seaMass(data, "uo__mq", c("L" = "lognormal"), c("0" = "", "r" = "rm", "1" = "one", "M" = "minimum", "C3" = "censored3", "C" = "censored", "C5" = "censored5"))

  # Run seaMass three-level model with MaxQuant's protein groups
  data <- import_MaxQuant("input_nomatching/proteinGroups.txt", "input_nomatching/evidence.txt")
  levels(data$Run)[10] <- "A1-6"
  execute_seaMass(data, "uo__mq0", c("L" = "lognormal"), c("0" = "", "r" = "rm", "1" = "one", "M" = "minimum", "C3" = "censored3", "C" = "censored", "C5" = "censored5"))



  open_seaMass <- function(path, name, model, effect, data.truth) {
    data <- group_fdr(open_seaMass_delta(open_seaMass_sigma(path), name))
    data <- data[data$Model == model & data$Effect == effect,]
    data <- add_truth(data, data.truth)
    return(data)
  }

  plot_pr(list(
    "qr d2 nMT -> sm3 eL m0" = open_seaMass(file.path(sm.path, "uo__mq__sm3_eL_m0.seaMass"), "nMG_tS", "contrast B - A", "Condition", data.truth)
  ), plot.fdr = F, legend.nrow = 1)



  data.uo__mq0__sm3_eL_m0.nMG_tS <- add_truth(group_fdr(open_seaMass_delta(open_seaMass_sigma(file.path(sm.path, "uo__mq0__sm3_eL_m0")), "nMG_tS")), data.truth)
  data.uo__mq0__sm3_eL_mR.nMG_tS <- add_truth(group_fdr(open_seaMass_delta(open_seaMass_sigma(file.path(sm.path, "uo__mq0__sm3_eL_mr")), "nMG_tS")), data.truth)
  data.uo__mq0__sm3_eL_m1.nMG_tS <- add_truth(group_fdr(open_seaMass_delta(open_seaMass_sigma(file.path(sm.path, "uo__mq0__sm3_eL_m1")), "nMG_tS")), data.truth)
  data.uo__mq0__sm3_eL_mC3.nMG_tS <- add_truth(group_fdr(open_seaMass_delta(open_seaMass_sigma(file.path(sm.path, "uo__mq0__sm3_eL_mC3")), "nMG_tS")), data.truth)
  data.uo__mq0__sm3_eL_mC.nMG_tS <- add_truth(group_fdr(open_seaMass_delta(open_seaMass_sigma(file.path(sm.path, "uo__mq0__sm3_eL_mC")), "nMG_tS")), data.truth)
  data.uo__mq0__sm3_eL_mC5.nMG_tS <- add_truth(group_fdr(open_seaMass_delta(open_seaMass_sigma(file.path(sm.path, "uo__mq0__sm3_eL_mC5")), "nMG_tS")), data.truth)

  data.uo__mq__sm3_eL_m0.nMG_tS <- add_truth(group_fdr(open_seaMass_delta(open_seaMass_sigma(file.path(sm.path, "uo__mq__sm3_eL_m0")), "nMG_tS")), data.truth)
  data.uo__mq__sm3_eL_mR.nMG_tS <- add_truth(group_fdr(open_seaMass_delta(open_seaMass_sigma(file.path(sm.path, "uo__mq__sm3_eL_mr")), "nMG_tS")), data.truth)
  data.uo__mq__sm3_eL_m1.nMG_tS <- add_truth(group_fdr(open_seaMass_delta(open_seaMass_sigma(file.path(sm.path, "uo__mq__sm3_eL_m1")), "nMG_tS")), data.truth)
  data.uo__mq__sm3_eL_mC3.nMG_tS <- add_truth(group_fdr(open_seaMass_delta(open_seaMass_sigma(file.path(sm.path, "uo__mq__sm3_eL_mC3")), "nMG_tS")), data.truth)
  data.uo__mq__sm3_eL_mC.nMG_tS <- add_truth(group_fdr(open_seaMass_delta(open_seaMass_sigma(file.path(sm.path, "uo__mq__sm3_eL_mC")), "nMG_tS")), data.truth)
  data.uo__mq__sm3_eL_mC5.nMG_tS <- add_truth(group_fdr(open_seaMass_delta(open_seaMass_sigma(file.path(sm.path, "uo__mq__sm3_eL_mC5")), "nMG_tS")), data.truth)


  plot_pr(list(
    "uo__mq__sm3_eL_m0.nMG_tS" = data.uo__mq__sm3_eL_m0.nMG_tS[data.uo__mq__sm3_eL_m0.nMG_tS$Batch == "Condition.contrast B - A",],
    "uo__mq__sm3_eL_mR.nMG_tS" = data.uo__mq__sm3_eL_mR.nMG_tS[data.uo__mq__sm3_eL_mR.nMG_tS$Batch == "Condition.contrast B - A",],
    "uo__mq__sm3_eL_m1.nMG_tS" = data.uo__mq__sm3_eL_m1.nMG_tS[data.uo__mq__sm3_eL_m1.nMG_tS$Batch == "Condition.contrast B - A",],
    "uo__mq__sm3_eL_mC3.nMG_tS" = data.uo__mq__sm3_eL_mC3.nMG_tS[data.uo__mq__sm3_eL_mC3.nMG_tS$Batch == "Condition.contrast B - A",],
    "uo__mq__sm3_eL_mC.nMG_tS" = data.uo__mq__sm3_eL_mC.nMG_tS[data.uo__mq__sm3_eL_mC.nMG_tS$Batch == "Condition.contrast B - A",],
    "uo__mq__sm3_eL_mC5.nMG_tS" = data.uo__mq__sm3_eL_mC5.nMG_tS[data.uo__mq__sm3_eL_mC5.nMG_tS$Batch == "Condition.contrast B - A",],
    "uo__mq0__sm3_eL_mC.nMG_tS" = data.uo__mq0__sm3_eL_mC.nMG_tS[data.uo__mq0__sm3_eL_mC.nMG_tS$Batch == "Condition.contrast B - A",],
    "uo__mq0__sm3_eL_mC5.nMG_tS" = data.uo__mq0__sm3_eL_mC5.nMG_tS[data.uo__mq0__sm3_eL_mC5.nMG_tS$Batch == "Condition.contrast B - A",]
  ), plot.fdr = F, legend.nrow = 2)








  data.uo__qr_d1_n0__sm3_eP_m0.nMG_tS <- add_truth(group_fdr(open_seaMass_delta(open_seaMass_sigma(file.path(sm.path, "uo__qr_d1_n0__sm3_eP_m0")), "nMG_tS")), data.truth)
  data.uo__qr_d1_n0__sm3_eP_m1.nTA_tS <- add_truth(group_fdr(open_seaMass_delta(open_seaMass_sigma(file.path(sm.path, "uo__qr_d1_n0__sm3_eP_m1")), "nTA_tS")), data.truth)
  data.uo__qr_d1_n0__sm3_eP_mM.nTA_tS <- add_truth(group_fdr(open_seaMass_delta(open_seaMass_sigma(file.path(sm.path, "uo__qr_d1_n0__sm3_eP_mM")), "nTA_tS")), data.truth)
  data.uo__qr_d1_n0__sm3_eP_mC.nTA_tS <- add_truth(group_fdr(open_seaMass_delta(open_seaMass_sigma(file.path(sm.path, "uo__qr_d1_n0__sm3_eP_mC")), "nTA_tS")), data.truth)
  data.uo__qr_d1_n0__sm3_eL_m0.nTA_tS <- add_truth(group_fdr(open_seaMass_delta(open_seaMass_sigma(file.path(sm.path, "uo__qr_d1_n0__sm3_eL_m0")), "nTA_tS")), data.truth)
  data.uo__qr_d1_n0__sm3_eL_m1.nTA_tS <- add_truth(group_fdr(open_seaMass_delta(open_seaMass_sigma(file.path(sm.path, "uo__qr_d1_n0__sm3_eL_m1")), "nTA_tS")), data.truth)
  data.uo__qr_d1_n0__sm3_eL_mM.nTA_tS <- add_truth(group_fdr(open_seaMass_delta(open_seaMass_sigma(file.path(sm.path, "uo__qr_d1_n0__sm3_eL_mM")), "nTA_tS")), data.truth)
  data.uo__qr_d1_n0__sm3_eL_mC.nTA_tS <- add_truth(group_fdr(open_seaMass_delta(open_seaMass_sigma(file.path(sm.path, "uo__qr_d1_n0__sm3_eL_mC")), "nTA_tS")), data.truth)

  plot_pr(list(
    "uo__qr_d1_n0__sm3_eP_m0.nTA_tS" = data.uo__qr_d1_n0__sm3_eP_m0.nTA_tS[data.uo__qr_d1_n0__sm3_eP_m0.nTA_tS$Batch == "Condition.contrast B - A",],
    "uo__qr_d1_n0__sm3_eL_m0.nTA_tS" = data.uo__qr_d1_n0__sm3_eL_m0.nTA_tS[data.uo__qr_d1_n0__sm3_eL_m0.nTA_tS$Batch == "Condition.contrast B - A",],
    "uo__qr_d1_n0__sm3_eP_m1.nTA_tS" = data.uo__qr_d1_n0__sm3_eP_m1.nTA_tS[data.uo__qr_d1_n0__sm3_eP_m1.nTA_tS$Batch == "Condition.contrast B - A",],
    "uo__qr_d1_n0__sm3_eL_m1.nTA_tS" = data.uo__qr_d1_n0__sm3_eL_m1.nTA_tS[data.uo__qr_d1_n0__sm3_eL_m1.nTA_tS$Batch == "Condition.contrast B - A",],
    "uo__qr_d1_n0__sm3_eP_mM.nTA_tS" = data.uo__qr_d1_n0__sm3_eP_mM.nTA_tS[data.uo__qr_d1_n0__sm3_eP_mM.nTA_tS$Batch == "Condition.contrast B - A",],
    "uo__qr_d1_n0__sm3_eL_mM.nTA_tS" = data.uo__qr_d1_n0__sm3_eL_mM.nTA_tS[data.uo__qr_d1_n0__sm3_eL_mM.nTA_tS$Batch == "Condition.contrast B - A",],
    "uo__qr_d1_n0__sm3_eP_mC.nTA_tS" = data.uo__qr_d1_n0__sm3_eP_mC.nTA_tS[data.uo__qr_d1_n0__sm3_eP_mC.nTA_tS$Batch == "Condition.contrast B - A",],
    "uo__qr_d1_n0__sm3_eL_mC.nTA_tS" = data.uo__qr_d1_n0__sm3_eL_mC.nTA_tS[data.uo__qr_d1_n0__sm3_eL_mC.nTA_tS$Batch == "Condition.contrast B - A",]
  ), plot.fdr = F, legend.nrow = 2)

  plot_pr(list(
    "uo__qr_d1_n0__sm3_eL_m0.nTA_tS" = data.uo__qr_d1_n0__sm3_eL_m0.nTA_tS[data.uo__qr_d1_n0__sm3_eL_m0.nTA_tS$Batch == "Condition.contrast B - A",],
    "uo__qr_d1_n0__sm3_eL_m1.nTA_tS" = data.uo__qr_d1_n0__sm3_eL_m1.nTA_tS[data.uo__qr_d1_n0__sm3_eL_m1.nTA_tS$Batch == "Condition.contrast B - A",],
    "uo__qr_d1_n0__sm3_eL_mM.nTA_tS" = data.uo__qr_d1_n0__sm3_eL_mM.nTA_tS[data.uo__qr_d1_n0__sm3_eL_mM.nTA_tS$Batch == "Condition.contrast B - A",],
    "uo__qr_d1_n0__sm3_eL_mC.nTA_tS" = data.uo__qr_d1_n0__sm3_eL_mC.nTA_tS[data.uo__qr_d1_n0__sm3_eL_mC.nTA_tS$Batch == "Condition.contrast B - A",]
  ), plot.fdr = F, legend.nrow = 2)


  data.uo__qr_d1_n0__sm3_eL_mC.nMA_tS <- add_truth(group_fdr(open_seaMass_delta(open_seaMass_sigma(file.path(sm.path, "uo__qr_d1_n0__sm3_eL_mC")), "nMA_tS")), data.truth)
  data.uo__qr_d1_n0__sm3_eL_mC.nQA_tS <- add_truth(group_fdr(open_seaMass_delta(open_seaMass_sigma(file.path(sm.path, "uo__qr_d1_n0__sm3_eL_mC")), "nQA_tS")), data.truth)
  data.uo__qr_d1_n0__sm3_eL_mC.nTA_tS <- add_truth(group_fdr(open_seaMass_delta(open_seaMass_sigma(file.path(sm.path, "uo__qr_d1_n0__sm3_eL_mC")), "nTA_tS")), data.truth)
  data.uo__qr_d1_n0__sm3_eL_mC.nTA_tW <- add_truth(group_fdr(open_seaMass_delta(open_seaMass_sigma(file.path(sm.path, "uo__qr_d1_n0__sm3_eL_mC")), "nTA_tW")), data.truth)

  plot_pr(list(
    "uo__qr_d1_n0__sm3_eL_mC.nMA_tS" = data.uo__qr_d1_n0__sm3_eL_mC.nMA_tS[data.uo__qr_d1_n0__sm3_eL_mC.nMA_tS$Batch == "Condition.contrast B - A",],
    "uo__qr_d1_n0__sm3_eL_mC.nQA_tS" = data.uo__qr_d1_n0__sm3_eL_mC.nQA_tS[data.uo__qr_d1_n0__sm3_eL_mC.nQA_tS$Batch == "Condition.contrast B - A",],
    "uo__qr_d1_n0__sm3_eL_mC.nTA_tS" = data.uo__qr_d1_n0__sm3_eL_mC.nTA_tS[data.uo__qr_d1_n0__sm3_eL_mC.nTA_tS$Batch == "Condition.contrast B - A",],
    "uo__qr_d1_n0__sm3_eL_mC.nTA_tW" = data.uo__qr_d1_n0__sm3_eL_mC.nTA_tW[data.uo__qr_d1_n0__sm3_eL_mC.nTA_tW$Batch == "Condition.contrast B - A",]
  ), plot.fdr = F, legend.nrow = 2)

  plot_pr(list(
    "uo__qr_d1_n0__sm3_eL_mC.nMA_tS" = data.uo__qr_d1_n0__sm3_eL_mC.nMA_tS[data.uo__qr_d1_n0__sm3_eL_mC.nMA_tS$Batch == "Condition.contrast C - A",],
    "uo__qr_d1_n0__sm3_eL_mC.nQA_tS" = data.uo__qr_d1_n0__sm3_eL_mC.nQA_tS[data.uo__qr_d1_n0__sm3_eL_mC.nQA_tS$Batch == "Condition.contrast C - A",],
    "uo__qr_d1_n0__sm3_eL_mC.nTA_tS" = data.uo__qr_d1_n0__sm3_eL_mC.nTA_tS[data.uo__qr_d1_n0__sm3_eL_mC.nTA_tS$Batch == "Condition.contrast C - A",],
    "uo__qr_d1_n0__sm3_eL_mC.nTA_tW" = data.uo__qr_d1_n0__sm3_eL_mC.nTA_tW[data.uo__qr_d1_n0__sm3_eL_mC.nTA_tW$Batch == "Condition.contrast C - A",]
  ), plot.fdr = F, legend.nrow = 2)





  plot_pr(list(
    "g=qr d=1 n=mt -> seaMass3 e=l m=c1 n=mt" = add_truth(group_fdr(open_delta_fit("uo__mq__sm3_eL_m0")), data.truth),
    "g=qr d=1 n=mt -> seaMass3 e=l m=c2 n=mt" = add_truth(group_fdr(open_delta_fit("msqrob__d_1__n_NA.seamass3__e_l__m_c2__n_mt")), data.truth),
    "g=qr d=1 n=mt -> seaMass3 e=l m=c3 n=mt" = add_truth(group_fdr(open_delta_fit("msqrob__d_1__n_NA.seamass3__e_l__m_c3__n_mt")), data.truth),
    #"g=qr d=1 n=mt -> seaMass3 e=l m=c4 n=mt" = add_truth(group_fdr(open_delta_fit("msqrob__d_1__n_NA.seamass3__e_l__m_c4__n_mt")), data.truth),
    "g=qr d=1 n=mt -> seaMass3 e=l m=c2_1 n=mt" = add_truth(group_fdr(open_delta_fit("msqrob__d_1__n_NA.seamass3__e_l__m_c2_1__n_mt")), data.truth)
    #"g=qr d=1 n=mt -> seaMass3 e=l m=c2_2 n=mt" = add_truth(group_fdr(open_delta_fit("msqrob__d_1__n_NA.seamass3__e_l__m_c2_2__n_mt")), data.truth)
  ), plot.fdr = F, ymax = 0.2, legend.nrow = 3)
  ggplot2::ggsave("msqrob__seamass3.pdf", width = 12, height = 6)

}



