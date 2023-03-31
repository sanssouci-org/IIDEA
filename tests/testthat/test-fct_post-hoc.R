test_that("posthoc_bound2", {
  p.values <- runif(100)
  thr <- sort(runif(100))
  posthoc_bound2(p.values,
    S = seq_along(p.values), thr = thr, lab = NULL,
    what = c("TP", "FDP"), all = FALSE
  )

  obj <- SansSouciSim(m = 502, rho = 0.5, n = 100, pi0 = 0.8, SNR = 3, prob = 0.5)
  res <- fit(obj, B = 100, alpha = 0.1)
  # post hoc bound on the set of all hypotheses
  posthoc_bound2(pValues(res), thr = thresholds(res))
})
