#' iTraq
#'
#' Blah.
#'
#' @import data.table
#' @export
itraq <- function() {

  #FullpooledRunA

  setwd("/mnt/storage/home/ad16243/bc3/Scratch/seaMass/")

  library(seaMass)

  input <- c("UItraqSFA.PP", "UItraqOA.PP")
  files <- file.path(system.file(input, package = "seaMassValidation"), "PeptideSummary.txt.bz2")

  funcs <- list(
    default = function(file) import_ProteinPilot(file),
    dpt = function(file) import_ProteinPilot(file, use.discordant.peptide.type = T),
    ni = function(file) import_ProteinPilot(file, use.no.itraq = T),
    ws = function(file) import_ProteinPilot(file, use.weak.signal = T),
    sp = function(file) import_ProteinPilot(file, use.shared.peptides = T),
    dpt_sp = function(file) import_ProteinPilot(file, use.discordant.peptide.type = T, use.shared.peptides = T),
    dpt_sp_ni = function(file) import_ProteinPilot(file, use.discordant.peptide.type = T, use.shared.peptides = T, use.no.itraq = T),
    dpt_sp_ni_ws = function(file) import_ProteinPilot(file, use.discordant.peptide.type = T, use.shared.peptides = T, use.no.itraq = T, use.weak.signal = T)
  )

  for (file in files) {
    for (i in 1:length(funcs)) {
      data <- funcs[[i]](file)

      for (j in 1:3) {
        message(names(funcs)[i])

        data.design <- new_assay_design(data)
        data.design$Assay <- factor(c("A1", "A2", "B1", "B2", "A3", "A4", "B3", "B4"))
        data.design$Condition <- factor(c("A", "A", "B", "B", "A", "A", "B", "B"))
        data.design$Sample <- factor(c("A1", "A2", "B1", "B2", "A3", "A4", "B3", "B4"))

        fit.sigma <- seaMass_sigma(
          data, data.design, run = F, path = file.path("/mnt/storage/home/ad16243/bc3/Scratch/seaMass", paste0(basename(dirname(file)), ".", names(funcs)[i], ".", j)),
          norm.groups = "_RAT",
          control = sigma_control(
            random.seed = j,
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
    }
  }


  results <- c("UItraqFPA.PP.default.1", "UItraqFPA.PP.default.2", "UItraqFPA.PP.default.3", "UItraqFPA.PP.dpt.1", "UItraqFPA.PP.dpt.2", "UItraqFPA.PP.dpt.3")
  results <- c("UItraqFPA.PP.default.1", "UItraqFPA.PP.dpt.1", "UItraqFPA.PP.ni.2", "UItraqFPA.PP.sp.1", "UItraqFPA.PP.ws.1")
  results <- c("UItraqFPA.PP.dpt.1", "UItraqFPA.PP.dpt.2", "UItraqFPA.PP.dpt.3", "UItraqFPA.PP.dpt_sp.2", "UItraqFPA.PP.dpt_sp.3")
  results <- c("UItraqFPA.PP.dpt_sp.2", "UItraqFPA.PP.dpt_sp.3", "UItraqFPA.PP.dpt_sp_ni.2", "UItraqFPA.PP.dpt_sp_ni.3")
  results <- c("UItraqFPA.PP.dpt_sp.2", "UItraqFPA.PP.dpt_sp.3", "UItraqFPA.PP.dpt_sp_ni_.2", "UItraqFPA.PP.dpt_sp_ni.3")
  data.fdr <- lapply(results, function(name) add_seaMass_spikein_truth(standardised_group_deviations_fdr(open_delta(open_sigma(name)))))
  names(data.fdr) <- results

  # add ground truth and plot precision-recall curve
  plot_pr(data.fdr, y.max=0.2)


}



