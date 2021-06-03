rm(list = ls())
gc()
options(scipen = 999)
# # # # #
library(data.table)
library(rstan)

# PRINT AND PLOT RESULTS

# load Bird's Eye models for reference
pred.names = c(
  "intercept",
  "coledu",
  "age",
  "age2",
  "married",
  "student",
  "lowstat",
  "coledu_lowstat"
)
load(file = 'data/analysis/full_fit_all.RData')

extracted_all = extract(full_fit_all , 
                        pars = c('beta',
                                 'phi',
                                 'psi',
                                 'eta','tau_eta',
                                 'lambda',
                                 'tau',
                                 'gamma',
                                 'rho','rho_star'), 
                        permuted = TRUE, 
                        inc_warmup = FALSE,
                        include = TRUE)

#load estimation and data files for Tunisia
load(file = 'data/analysis/full_fit_tunisia.RData')
load(file = 'data/analysis/cleaned.tunisia.short.RData')

# results
print(
  full_fit_tunisia,
  pars = c(
    'beta',
    'phi',
    'psi',
    'eta',
    'tau_eta',
    'lambda',
    'tau',
    'gamma',
    'rho',
    'rho_star'
  )
)

print(full_fit_tunisia, pars = c('beta'))
print(full_fit_tunisia, pars = c('phi'))
print(full_fit_tunisia, pars = c('psi'))
print(full_fit_tunisia, pars = c('eta'))
print(full_fit_tunisia, pars = c('lambda'))
print(full_fit_tunisia, pars = c('tau'))
print(full_fit_tunisia, pars = c('gamma'))
print(full_fit_tunisia, pars = c('rho'))
print(full_fit_tunisia, pars = c('rho_star'))
print(full_fit_tunisia, pars = c('lp__'))

summary_fit_tunisia = summary(full_fit_tunisia)

################# MAIN PLOTS

# relative deprivation effects
pred.names.tunisia = c(
  "intercept",
  "coledu",
  "age",
  "age2",
  "married",
  "student",
  "lowstat",
  "coledu_lowstat",
  "pop_10plus_2014",
  "population_density",
  "pct_agri_2014",
  "pct_higher_edu_2014",
  "unemp_rate_2014",
  "dip_unemp_rate_2014",
  "sqrt_distance_to_libya",
  "sqrt_post_rev_protest_events",
  "pct_nahdha_2014",
  "pct_nahdha_2011"
)

plot.names.tunisia =  c(
  "Intercept",
  "College edu.",
  "Age",
  "Age^2",
  "Married",
  "Student",
  "Low status",
  "College edu.*Low status",
  "Population",
  "Population density",
  "% Agriculture",
  "% Higher edu.",
  "Unemp. rate",
  "Grad. unemp. rate",
  "Distance to Libya (sqrt.)",
  "Post-rev. protest (sqrt.)",
  "% Ennahdha 2014",
  "% Ennahdha 2011"
)

summary_fit_tunisia = summary(full_fit_tunisia)

extracted_tunisia = extract(
  full_fit_tunisia ,
  pars = c(
    'beta',
    'phi',
    'psi',
    'eta',
    'tau_eta',
    'lambda',
    'tau',
    'gamma',
    'rho',
    'rho_star'
  ),
  permuted = TRUE,
  inc_warmup = FALSE,
  include = TRUE
)

rel.dep.effects = data.table(expand.grid(coledu = c(0, 1),
                                         lowstat = c(0, 1)))

rel.dep.effects$coledu_lowstat = rel.dep.effects$coledu * rel.dep.effects$lowstat

rel.dep.effects.center.tunisia = c(
  mean(dt_tunisia_short$coledu),
  mean(dt_tunisia_short$lowstat),
  mean(dt_tunisia_short$coledu_lowstat)
)
rel.dep.effects.scale.tunisia = c(
  sd(dt_tunisia_short$coledu),
  sd(dt_tunisia_short$lowstat),
  sd(dt_tunisia_short$coledu_lowstat)
)

rel.dep.effects.std.tunisia = (rel.dep.effects - t(array(
  rel.dep.effects.center.tunisia, rev(dim(rel.dep.effects))
))) / t(array(rel.dep.effects.scale.tunisia, rev(dim(rel.dep.effects))))

logit.total.effect.tunisia =
  sapply(1:dim(extracted_tunisia$beta)[1], function(i) {
    as.matrix(rel.dep.effects.std.tunisia) %*%
      extracted_tunisia$beta[i, match(colnames(rel.dep.effects), pred.names.tunisia)]
  })

row.names(logit.total.effect.tunisia) =
  c(
    "no college education & high status",
    "college educated & high status",
    "no college education & low status",
    "college educated & low status"
  )

pdf(file = 'plots/tunisia/total.effect_relative.deprivation_tunisia.pdf',
    height = 5,
    width = 7.5)
plot(
  density(logit.total.effect.tunisia[1,]),
  xlim = c(
    min(logit.total.effect.tunisia),
    max(logit.total.effect.tunisia)
  ),
  main = "",
  col = NA,
  xlab = 'logit probability of recruitment',
  zero.line = FALSE
)
for (i in 1:dim(logit.total.effect.tunisia)[1]) {
  d <- density(logit.total.effect.tunisia[i,])
  polygon(
    d,
    main = '',
    xaxt = "n",
    yaxt = "n",
    ylab = "",
    xlab = "",
    border = NA,
    col = adjustcolor(c(
      "skyblue", "lightgreen", "lightcoral", "grey"
    )[i], 0.75)
  )
}
abline(
  v = rowMeans(logit.total.effect.tunisia),
  lty = 2,
  col = adjustcolor(c("blue", "darkgreen", "red", "black"), 0.5),
  lwd = 2
)
text(
  x = rowMeans(logit.total.effect.tunisia),
  y = 2,
  labels = rownames(logit.total.effect.tunisia),
  srt = 60,
  cex = 0.85,
  pos = 1,
  col = c("blue", "darkgreen", "red", "black")
)
dev.off()

# regression coefficients

pdf(file = 'plots/tunisia/regression.coefficients_ind_tunisia.pdf',
    height = 8,
    width = 12)
par(mfrow = c(3, 3))
for (i in 1:dim(extracted_all$beta)[2]) {
  plot(
    density(extracted_tunisia$beta[, i]),
    xlim = c(ifelse(
      min(extracted_tunisia$beta[, i]) < 0,
      min(extracted_tunisia$beta[, i]) - 3,-3
    ),
    ifelse(
      max(extracted_tunisia$beta[, i]) > 0,
      max(extracted_tunisia$beta[, i]) + 3,
      3
    )),
    zero.line = FALSE,
    main = plot.names.tunisia[i],
    col = NA,
    xlab = 'logit probability of recruitment'
  )
  d <- density(extracted_tunisia$beta[, i])
  polygon(
    d,
    main = '',
    xaxt = "n",
    yaxt = "n",
    ylab = "",
    xlab = "",
    border = NA,
    col = adjustcolor("orange", 0.75)
  )
  abline(
    v = colMeans(extracted_tunisia$beta)[i],
    lty = 2,
    col = adjustcolor("purple", 0.75),
    lwd = 2
  )
  abline(v = 0, lty = 1, lwd = 1)
  if (i == 1) {
    legend.var = c(
      paste("E[beta_", i, "] = ", round(mean(
        extracted_tunisia$beta[, i]
      ), 3), sep = ""),
      paste("sd(beta_", i, ") = ", round(sd(
        extracted_tunisia$beta[, i]
      ), 3), sep = ""),
      ifelse(
        mean(extracted_tunisia$beta[, i]) > 0,
        paste("Pr(beta_", i, ">0) = ", round(
          sum(extracted_tunisia$beta[, i] > 0) / length(extracted_tunisia$beta[, i]),
          3
        ), sep = ""),
        paste("Pr(beta_", i, "<0) = ", round(
          sum(extracted_tunisia$beta[, i] < 0) / length(extracted_tunisia$beta[, i]),
          3
        ), sep = "")
      )
    )
  } else{
    legend.var = c(
      paste("E[beta_", i, "] = ", round(mean(
        extracted_tunisia$beta[, i]
      ), 3), sep = ""),
      paste("sd(beta_", i, ") = ", round(sd(
        extracted_tunisia$beta[, i]
      ), 3), sep = ""),
      ifelse(
        mean(extracted_tunisia$beta[, i]) > 0,
        paste("Pr(beta_", i, ">0) = ", round(
          sum(extracted_tunisia$beta[, i] > 0) / length(extracted_tunisia$beta[, i]),
          3
        ), sep = ""),
        paste("Pr(beta_", i, "<0) = ", round(
          sum(extracted_tunisia$beta[, i] < 0) / length(extracted_tunisia$beta[, i]),
          3
        ), sep = "")
      ),
      paste("mean(X) = ", round(mean(
        unlist(dt_tunisia_short[, pred.names[i], with = FALSE])
      ), 2)),
      paste("sd(X) = ", round(sd(
        unlist(dt_tunisia_short[, pred.names[i], with = FALSE])
      ), 2))
    )
  }
  legend(
    ifelse(colMeans(extracted_tunisia$beta)[i] > 0, "topleft", "topright"),
    legend = legend.var,
    cex = 0.9,
    bg = "white"
  )
}
dev.off()

# lambda density

pdf(file = 'plots/tunisia/lambda_tunisia.pdf',
    height = 5,
    width = 5)
plot(
  density(extracted_tunisia$lambda),
  xlim = c(0, 1),
  zero.line = FALSE,
  main = "",
  col = NA,
  xlab = 'logit probability of recruitment'
)
polygon(
  density(extracted_tunisia$lambda),
  main = '',
  xaxt = "n",
  yaxt = "n",
  ylab = "",
  xlab = "",
  border = NA,
  col = adjustcolor("orange", 0.75)
)
abline(
  v = mean(extracted_tunisia$lambda),
  lty = 2,
  col = adjustcolor("purple", 0.75),
  lwd = 2
)
legend(
  "topright",
  legend = paste("MCMC mean:", round(mean(
    extracted_tunisia$lambda
  ), 3)),
  lty = 2,
  lwd = 2,
  col = 'purple'
)
dev.off()