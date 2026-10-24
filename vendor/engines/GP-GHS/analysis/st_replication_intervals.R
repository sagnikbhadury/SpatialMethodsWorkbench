## ================================================================
## Bootstrap intervals for the pooled ratios of Table 3 (2026-10-22).
##
## The ratio of a comparison type is sum(common edges) / sum(expected
## common edges) over the zone-pair comparisons of that type. The
## interval resamples the comparisons of the type with replacement
## (B = 10,000) and takes the 2.5% and 97.5% quantiles of the pooled
## ratio; it treats the comparisons of a type as exchangeable, although
## comparisons that share a zone are not independent. Also computes the
## same interval for the pooled overlap.
##
## Usage (R 4.4.1): Rscript st_replication_intervals.R
## Input : st_results/replication/unit_pairs.csv
## Output: st_results/replication/table3_intervals.csv
## ================================================================
set.seed(20261022)
up <- read.csv("st_results/replication/unit_pairs.csv", stringsAsFactors = FALSE)
cats <- c("replicate sections, same zone", "replicate sections, different zone",
          "same slide, different zone", "different patients, same zone",
          "different patients, different zone")
B <- 10000
out <- do.call(rbind, lapply(cats, function(ct) {
  d <- up[up$category == ct, ]
  n <- nrow(d)
  ratio <- sum(d$both) / sum(d$expected)
  smaller <- pmin(d$edges_a, d$edges_b)
  overlap <- sum(d$both) / sum(smaller)
  bs <- replicate(B, { i <- sample.int(n, n, replace = TRUE)
    c(sum(d$both[i]) / sum(d$expected[i]), sum(d$both[i]) / sum(smaller[i])) })
  data.frame(category = ct, n = n, common = sum(d$both), expected = round(sum(d$expected), 1),
             ratio = ratio, ratio_lo = quantile(bs[1, ], 0.025), ratio_hi = quantile(bs[1, ], 0.975),
             overlap = overlap, overlap_lo = quantile(bs[2, ], 0.025), overlap_hi = quantile(bs[2, ], 0.975),
             min_ratio = min(d$fold), max_ratio = max(d$fold))
}))
rownames(out) <- NULL
write.csv(out, "st_results/replication/table3_intervals.csv", row.names = FALSE)
print(out, digits = 3)
