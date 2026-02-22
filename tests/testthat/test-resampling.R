library(testthat)

test_that("BootData returns same number of rows.", {
  d <- data.frame(time = 1:5, status = c(1,1,0,2,1))
  withr::local_seed(20)
  b <- BootData(d)
  expect_equal(nrow(b), nrow(d))
  expect_equal(names(b), names(d))
})

test_that("StratBoot preserves strata and total n.", {
  d <- data.frame(time = 1:6, status = 1, arm = c(0,0,0,1,1,1), strata = c(1,1,1,2,2,2))
  withr::local_seed(21)
  b <- StratBoot(d)
  expect_equal(nrow(b), 6)
  expect_equal(sort(unique(b$strata)), c(1, 2))
})

test_that("PermData permutes arm labels without replacement.", {
  d <- data.frame(time = 1:5, status = 1, arm = c(0, 0, 1, 1, 1))
  withr::local_seed(22)
  p <- PermData(d)
  expect_equal(sort(table(p$arm)), sort(table(d$arm)))
  expect_equal(length(unique(p$arm)), 2)
  # With same seed, permutation is deterministic; arm counts preserved
  withr::local_seed(22)
  p2 <- PermData(d)
  expect_equal(table(p$arm), table(p2$arm))
})

test_that("StratPerm preserves stratum sizes and permutes within stratum.", {
  d <- data.frame(
    time = 1:6,
    status = 1,
    arm = c(0, 1, 0, 1, 0, 1),
    strata = c(1, 1, 1, 2, 2, 2)
  )
  withr::local_seed(23)
  p <- StratPerm(d)
  expect_equal(nrow(p), 6)
  expect_equal(as.vector(table(d$strata)), as.vector(table(p$strata)))
  expect_equal(as.vector(table(d$arm)), as.vector(table(p$arm)))
})

test_that("CalcP returns p-value in [0, 1].", {
  expect_equal(CalcP(rep(0, 100)), 2 * 1 / 101)
  expect_equal(CalcP(rep(1, 100)), 1)
  expect_true(CalcP(c(0, 0, 1)) >= 0 && CalcP(c(0, 0, 1)) <= 1)
})

test_that("BootCI returns lower and upper.", {
  withr::local_seed(99)
  x <- rnorm(1000)
  ci <- BootCI(x, alpha = 0.05)
  expect_equal(names(ci), c("lower", "upper"))
  expect_lt(ci[["lower"]], ci[["upper"]])
  expect_lt(ci[["lower"]], mean(x))
  expect_gt(ci[["upper"]], mean(x))
})

test_that("BootSim runs and returns expected structure.", {
  d <- GenTwoSampleData(n1 = 25, n0 = 25, tau = 5)
  d$strata <- 1
  obs <- SumStats(d, sum_stat = "AUC", param = 5, return_strata = FALSE)
  withr::local_seed(24)
  boot <- BootSim(data = d, obs_stats = obs, sum_stat = "AUC", param = 5, alpha = 0.05, reps = 50)
  expect_true(all(c("pval", "cis", "sim") %in% names(boot)))
  expect_equal(nrow(boot$sim), 50)
  expect_true(boot$pval >= 0 && boot$pval <= 1)
  expect_true(all(c("lower", "upper") %in% names(boot$cis)))
})

test_that("PermSim runs and returns expected structure.", {
  d <- GenTwoSampleData(n1 = 25, n0 = 25, tau = 5)
  d$strata <- 1
  obs <- SumStats(d, sum_stat = "Rate", param = 3, return_strata = FALSE)
  withr::local_seed(25)
  perm <- PermSim(data = d, obs_stats = obs, sum_stat = "Rate", param = 3, alpha = 0.05, reps = 50)
  expect_true(all(c("pvals", "sim") %in% names(perm)))
  expect_equal(nrow(perm$sim), 50)
  expect_true(all(perm$pvals$perm_p >= 0 & perm$pvals$perm_p <= 1))
})
