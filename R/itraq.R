#' iTraq
#'
#' Blah.
#'
#' @import data.table
#' @export
itraq <- function() {
  data <- rbind(
    seaMass::import_ProteinPilot(system.file("UItraqOABC.PP/PeptideSummary1.txt.bz2", package = "seaMassValidation")),
    seaMass::import_ProteinPilot(system.file("UItraqOABC.PP/PeptideSummary2.txt.bz2", package = "seaMassValidation"))
  )


  #FullpooledRunA

  setwd("/mnt/storage/home/ad16243/bc3/Scratch/seaMass/")

  library(seaMass)

  filename <- "~/Repositories/bayesprotValidation/inst/data/UnwinItraqSpikeinFullpooledRunA.ProteinPilot/20181107_AT_8_PLEX_1_full_pooled_PeptideSummary.txt.bz2"
  dds <- list(
    default = import_ProteinPilot(filename),
    dpt = import_ProteinPilot(filename, use.discordant.peptide.type = T),
    ni = import_ProteinPilot(filename, use.no.itraq = T),
    ws = import_ProteinPilot(filename, use.weak.signal = T),
    sp = import_ProteinPilot(filename, use.shared.peptides = T),
    dpt_sp = import_ProteinPilot(filename, use.discordant.peptide.type = T, use.shared.peptides = T),
    dpt_sp = import_ProteinPilot(filename, use.discordant.peptide.type = T, use.shared.peptides = T),
    dpt_sp_ni = import_ProteinPilot(filename, use.discordant.peptide.type = T, use.shared.peptides = T, use.no.itraq = T),
    dpt_sp_ni_ws = import_ProteinPilot(filename, use.discordant.peptide.type = T, use.shared.peptides = T, use.no.itraq = T, use.weak.signal = T)
  )

  data <- dds[[1]]
  data.design <- new_assay_design(data)
  data.design$Assay <- factor(c("A1", "A2", "B1", "B2", "A3", "A4", "B3", "B4"))
  data.design$Condition <- factor(c("A", "A", "B", "B", "A", "A", "B", "B"))
  data.design$Sample <- factor(c("A1", "A2", "B1", "B2", "A3", "A4", "B3", "B4"))

  # fit.sigma <- seaMass_sigma(
  #   data,
  #   data.design,
  #   path = "FullpooledRunA3new5",
  #   norm.groups = "_RAT",
  #   control = sigma_control(norm.model = "median")
  # )
  #
  # fit.delta <- seaMass_delta(fit.sigma)

  for (i in 1:length(dds)) {
    data <- dds[[i]]

    # fit.sigma <- seaMass_sigma(
    #   data, data.design, run = F, path = file.path("/mnt/storage/home/ad16243/bc4/Scratch/seaMass/", paste0("FP.", names(dds)[i])),
    #   norm.groups = "_RAT",
    #   control = sigma_control(
    #     norm.model = "median",
    #     schedule = schedule_slurm(
    #       submit.prefix = "ssh bc4login.acrc.bris.ac.uk /mnt/storage/home/ad16243/Scratch/seaMass/",
    #       cpus_per_task = 14,
    #       mem = "64000m",
    #       mail_user = "ad16243@bristol.ac.uk"
    #     )
    #   )
    # )

    fit.sigma <- seaMass_sigma(
      data, data.design, run = F, path = file.path("/mnt/storage/home/ad16243/bc3/Scratch/seaMass", paste0("FP.", names(dds)[i])),
      norm.groups = "_RAT",
      control = sigma_control(
        norm.model = "median",
        schedule = schedule_pbs(
          submit.prefix = "ssh bluecrystalp3.acrc.bris.ac.uk /newhome/ad16243/Scratch/seaMass/",
          q = "veryshort",
          walltime = "12:00:00",
          ppn = 16,
          mem = "64000m",
          M = "ad16243@bristol.ac.uk"
        )
      )
    )

    fit.delta <- seaMass_delta(fit.sigma)
    run(fit.sigma)
  }


  results <- c("default_gtm", "udpt_gtm", "uni_gtm", "uws_gtm", "usp_gtm")
  data.fdr <- lapply(results, function(name) add_seaMass_spikein_truth(standardised_group_deviations_fdr(open_delta(open_sigma(name)))))
  names(data.fdr) <- results

  # add ground truth and plot precision-recall curve
  plot_pr(data.fdr, y.max=0.2)


}



