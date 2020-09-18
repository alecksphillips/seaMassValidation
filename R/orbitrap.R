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
  library(MSstats)
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

  DT.raw <- fread("input/evidence.txt")
  for (assay in assays.rm) DT.raw <- DT.raw[Experiment != assay,]
  fwrite(DT.raw, "input4v4/evidence.txt", sep = "\t", quote = F)
  

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

  ### MSstats ########################################################################
  # Function to run msstats on maxquant's output
  execute_MSstats <- function(name, normalize = T, rm.missing = 1.0, useLeadingRazorProteins = F) {
    
    if (!file.exists(paste(name, "rds", sep = "."))) {
      
      evidence <- read.csv("input4v4/evidence.txt", sep ='\t')
      #annotation <- data.frame(
      #  Raw.file = c("VH_210319_A1-5","VH_210319_A1-6", "VH_210319_A1-7", "VH_210319_A2-4", "VH_210319_A3-5", "VH_210319_A3-6", "VH_210319_A4-5", "VH_210319_A4-6","VH_210319_A5-3", "VH_210319_A6-8", "VH_210319_B1-3", "VH_210319_B2-4", "VH_210319_B3-5", "VH_210319_B3-6", "VH_210319_B4-5", "VH_210319_B4-6", "VH_210319_B5-3", "VH_210319_B6-3", "VH_210319_C1-3", "VH_210319_C2-3", "VH_210319_C3-3", "VH_210319_C4-3"),
      #  Condition = c("A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "B", "B", "B", "B", "B", "B", "B", "B", "C", "C", "C", "C"),
      #  BioReplicate = c("A1", "A1", "A1", "A2", "A3", "A3", "A4", "A4", "A5", "A6", "B1", "B2", "B3", "B3", "B4", "B4", "B5", "B6", "C1", "C2", "C3", "C4"),
      #  IsotopeLabelType = rep_len(c("L"),22),
      #  Run = c("A1-5" ,"A1-6" ,"A1-7" ,"A2-4" ,"A3-5" ,"A3-6" ,"A4-5" ,"A4-6" ,"A5-3" ,"A6-8" ,"B1-3" ,"B2-4" ,"B3-5" ,"B3-6" ,"B4-5" ,"B4-6" ,"B5-3" ,"B6-3" ,"C1-3" ,"C2-3" ,"C3-3" ,"C4-3")
      # )
      annotation <- data.frame(
        Raw.file = c("VH_210319_A1-6", "VH_210319_A3-6", "VH_210319_A5-3", "VH_210319_A6-8", "VH_210319_B2-4", "VH_210319_B3-5", "VH_210319_B4-5", "VH_210319_B5-3", "VH_210319_C1-3", "VH_210319_C2-3", "VH_210319_C3-3", "VH_210319_C4-3"),
        Condition = c("A", "A", "A", "A", "B", "B", "B", "B", "C", "C", "C", "C"),
        BioReplicate = c("A1","A3", "A5", "A6", "B2", "B3", "B4", "B5", "C1", "C2", "C3", "C4"),
        IsotopeLabelType = rep_len(c("L"),12),
        Run = c("A1-6", "A3-6", "A5-3" ,"A6-8" ,"B2-4" ,"B3-5", "B4-5", "B5-3", "C1-3" ,"C2-3" ,"C3-3" ,"C4-3")
      )
      proteinGroups <- read.csv("input4v4/proteinGroups.txt", sep ='\t')
    
      
      if (useLeadingRazorProteins) {
        inputData <- MaxQtoMSstatsFormat(
          evidence = evidence, annotation = annotation, proteinGroups = proteinGroups,
          proteinID = 'Leading.razor.protein'
        )
        
      } else {
        inputData <- MaxQtoMSstatsFormat(
          evidence = evidence, annotation = annotation, proteinGroups = proteinGroups
        )
      }
    
      QuantData <- dataProcess(inputData)
    
      exposure_correction_MSstats <- function(data, normalisation.proteins = levels(data$Protein)) {
        DT.ProcessedData <- as.data.table(data$ProcessedData)
        DT.RunlevelData <- as.data.table(data$RunlevelData)
      
        # compute exposures
        DT.exposure <- DT.ProcessedData[, .(exposure = median(ABUNDANCE[PROTEIN %in% normalisation.proteins], na.rm = T)), keyby = RUN]
        DT.exposure[,exposure := exposure - mean(exposure, na.rm = T)]
      
        DT.exposure.intensity <- DT.ProcessedData[, .(exposure = median(log2(INTENSITY[PROTEIN %in% normalisation.proteins]), na.rm = T)), keyby = RUN]
        DT.exposure.intensity[,exposure := exposure - mean(exposure, na.rm = T)]
      
        DT.Runlevel.exposure <- DT.RunlevelData[, .(exposure = median(LogIntensities[Protein %in% normalisation.proteins], na.rm = T)), keyby = RUN]
        DT.Runlevel.exposure[,exposure := exposure - mean(exposure, na.rm = T)]
      
        # apply exposure
        DT.ProcessedData[DT.exposure, ABUNDANCE := ABUNDANCE - exposure, on = "RUN"]
        DT.ProcessedData[DT.exposure, INTENSITY := INTENSITY * 2^(-exposure), on = "RUN"]
        DT.RunlevelData[DT.Runlevel.exposure, LogIntensities := LogIntensities - exposure, on = "RUN"]
      
        data$ProcessedData <- setDF(DT.ProcessedData)
        data$RunlevelData <- setDF(DT.RunlevelData)
        return(data)
      }
    
      de_MSstats <- function(data, rm.missing = 1.0) {
        DT <- as.data.table(data)
      
        #DT[, .(Group = Protein, issue, MissingPercentage, ImputationPercentage, m = log2FC, s = SE, df = DF, t = Tvalue, pvalue, qvalue = adj.pvalue)]
      
        DT <- DT[, .(
          Group = Protein,
          m = log2FC,
          s = SE,
          df = DF,
          tvalue = Tvalue,
          pvalue = pvalue,
          qvalue = adj.pvalue,
          MissingPercentage
        )]
        #DT$pvalue[is.na(DT$pvalue)] <- 1
        #DT$p.value[is.nan(DT$FDR)] <- 1
        #DT$qvalue[DT$pvalue == 0] <- 0
        #DT$qvalue[is.nan(DT$qvalue)] <- 1
      
        DT[, abs.m := abs(m)]
        setorder(DT, pvalue, -abs.m, MissingPercentage, na.last = T)
        DT[, abs.m := NULL]
      
        DT <- DT[MissingPercentage <= rm.missing,]
        
        DT[, qvalue := p.adjust(pvalue, method = "BH")]
        setorder(DT, qvalue)
      
        data <- setDF(DT)
        return(data)
      }
    
      #Each row corresponds to a comparison being made, columns are conditions
      comparison <- matrix(c(-1,1,0), nrow=1)
      row.names(comparison) <- "B-A"
      
      normProteins <- levels(QuantData$ProcessedData$PROTEIN)
      normProteins <- normProteins[grepl("_RAT$",normProteins)]
      
      if (normalize){
        QuantData <- exposure_correction_MSstats(QuantData, normalisation.proteins = normProteins)
      }
      
      # Tests for differentially abundant proteins with models:
      testResultOneComparison <- groupComparison(contrast.matrix=comparison, data=QuantData)
      
      results <- de_MSstats(testResultOneComparison$ComparisonResult, rm.missing = rm.missing)
      results <- seaMass::add_seaMass_spikein_truth(results)
      
      saveRDS(results, paste(name, "rds", sep = "."))
    }
    else {
      #Read previous results
      results <- readRDS(paste(name, "rds", sep = "."))
    }
    
    g <- seaMass::plot_volcano(results, x.col = "m", y.col = "qvalue", stdev.col = "s")
    ggsave(paste0(name, ".pdf"))
  
    return(g)
  
  }

  execute_MSstats("MSstats_truthNorm", normalize = T, rm.missing = 1.0)
  execute_MSstats("MSstats_msStatNorm", normalize = F, rm.missing = 1.0)
  execute_MSstats("MSstats_truthNorm_rmMissing", normalize = T, rm.missing = 0.0)
  execute_MSstats("MSstats_msStatNorm_rmMissing", normalize = F, rm.missing = 0.0)
  execute_MSstats("MSstats_truthNorm_rmMissing0.5", normalize = T, rm.missing = 0.5)
  execute_MSstats("MSstats_msStatNorm_rmMissing0.5", normalize = F, rm.missing = 0.5)
  
  execute_MSstats("MSstats_truthNorm_razor", normalize = T, rm.missing = 1.0, useLeadingRazorProteins = T)
  execute_MSstats("MSstats_msStatNorm_razor", normalize = F, rm.missing = 1.0, useLeadingRazorProteins = T)
  execute_MSstats("MSstats_truthNorm_rmMissing_razor", normalize = T, rm.missing = 0.0, useLeadingRazorProteins = T)
  execute_MSstats("MSstats_msStatNorm_rmMissing_razor", normalize = F, rm.missing = 0.0, useLeadingRazorProteins = T)
  execute_MSstats("MSstats_truthNorm_rmMissing0.5_razor", normalize = T, rm.missing = 0.5, useLeadingRazorProteins = T)
  execute_MSstats("MSstats_msStatNorm_rmMissing0.5_razor", normalize = F, rm.missing = 0.5, useLeadingRazorProteins = T)
  
  plot_pr(list(
    "MSstats n=mt rm=F pid=proteins" = readRDS("MSstats_truthNorm.rds"),
    "MSstats n=ms rm=F pid=proteins" = readRDS("MSstats_msStatNorm.rds"),
    #"MSstats n=mt rm=T pid=proteins" = readRDS("MSstats_truthNorm_rmMissing.rds"),
    #"MSstats n=ms rm=T pid=proteins" = readRDS("MSstats_msStatNorm_rmMissing.rds"),
    #"MSstats n=mt rm=0.5 pid=proteins" = readRDS("MSstats_truthNorm_rmMissing0.5.rds")
    #"MSstats n=ms rm=0.5 pid=proteins" = readRDS("MSstats_msStatNorm_rmMissing0.5.rds")
    "MSstats n=mt rm=F pid=razor" = readRDS("MSstats_truthNorm_razor.rds"),
    "MSstats n=ms rm=F pid=razor" = readRDS("MSstats_msStatNorm_razor.rds")
    #"MSstats n=mt rm=T pid=razor" = readRDS("MSstats_truthNorm_rmMissing_razor.rds"),
    #"MSstats n=ms rm=T pid=razor" = readRDS("MSstats_msStatNorm_rmMissing_razor.rds"),
    #"MSstats n=mt rm=0.5 pid=razor" = readRDS("MSstats_truthNorm_rmMissing0.5_razor.rds"),
    #"MSstats n=ms rm=0.5 pid=razor" = readRDS("MSstats_msStatNorm_rmMissing0.5_razor.rds")
  ), plot.fdr = T, legend.nrow = 2)
  ggplot2::ggsave("msstats_comparison.pdf", width = 12, height = 6)














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



