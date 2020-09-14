library(seaMass)

schedule <- schedule_pbs(
  submit.prefix = "ssh bluecrystalp3.acrc.bris.ac.uk /newhome/ad16243/Scratch/seaMass/",
  q = "veryshort",
  walltime = "12:00:00",
  ppn = 16,
  mem = "64000m",
  M = "ad16243@bristol.ac.uk"
)

schedule <- schedule_slurm(
  submit.prefix = "ssh bc4login.acrc.bris.ac.uk /mnt/storage/home/ad16243/Scratch/seaMass/",
  cpus_per_task = 14,
  mem = "64000m",
  mail_user = "ad16243@bristol.ac.uk"
)

itraq("/mnt/storage/home/ad16243/bc4/Scratch/seaMass", schedule)

# results <- c("UItraqFPA.PP.default.1", "UItraqFPA.PP.default.2", "UItraqFPA.PP.default.3", "UItraqFPA.PP.dpt.1", "UItraqFPA.PP.dpt.2", "UItraqFPA.PP.dpt.3")
# results <- c("UItraqFPA.PP.default.1", "UItraqFPA.PP.dpt.1", "UItraqFPA.PP.ni.2", "UItraqFPA.PP.sp.1", "UItraqFPA.PP.ws.1")
# results <- c("UItraqFPA.PP.dpt.1", "UItraqFPA.PP.dpt.2", "UItraqFPA.PP.dpt.3", "UItraqFPA.PP.dpt_sp.2", "UItraqFPA.PP.dpt_sp.3")
# results <- c("UItraqFPA.PP.dpt_sp.2", "UItraqFPA.PP.dpt_sp.3", "UItraqFPA.PP.dpt_sp_ni.2", "UItraqFPA.PP.dpt_sp_ni.3")
# results <- c("UItraqFPA.PP.dpt_sp.2", "UItraqFPA.PP.dpt_sp.3", "UItraqFPA.PP.dpt_sp_ni_.2", "UItraqFPA.PP.dpt_sp_ni.3")
# data.fdr <- lapply(results, function(name) add_seaMass_spikein_truth(standardised_group_deviations_fdr(open_delta(open_sigma(name)))))
# names(data.fdr) <- results

# add ground truth and plot precision-recall curve
# plot_pr(data.fdr, y.max=0.2)
