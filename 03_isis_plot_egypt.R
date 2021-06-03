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

#load estimation and data files for Egypt
load(file = 'data/analysis/full_fit_egypt.RData')
load(file = 'data/analysis/cleaned.egypt.short.RData')

# source additional utility functions
# source("utils/nb_data_funs.R")

print(
  full_fit_egypt,
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

print(full_fit_egypt, pars = c('beta'))
print(full_fit_egypt, pars = c('phi'))
print(full_fit_egypt, pars = c('psi'))
print(full_fit_egypt, pars = c('eta'))
print(full_fit_egypt, pars = c('lambda'))
print(full_fit_egypt, pars = c('tau'))
print(full_fit_egypt, pars = c('gamma'))
print(full_fit_egypt, pars = c('rho'))
print(full_fit_egypt, pars = c('rho_star'))
print(full_fit_egypt, pars = c('lp__'))

summary_fit_egypt = summary(full_fit_egypt)

################# MAIN PLOTS

# relative deprivation effects
pred.names.egypt = c(
  "intercept",
  "coledu",
  "age",
  "age2",
  "married",
  "student",
  "lowstat",
  "coledu_lowstat",
  "population_density",
  "total_population_2006",
  "christian_2006_pct",
  "university_2006_pct",
  "agriculture_2006_pct",
  "mursi_vote_2012_pct",
  "sqrt_killed_at_rabaa",
  "unemployment_2013q4_pct",
  "protest_post_Mubarak"
)

plot.names.egypt = c(
  "Intercept",
  "College edu.",
  "Age",
  "Age^2",
  "Married",
  "Student",
  "Low status",
  "College edu.*Low status",
  "Population density",
  "Population",
  "% Christian",
  "% College edu.",
  "% Agriculture",
  "% Mursi 2012",
  "Killed at Rabaa (sqrt.)",
  "Unemp. rate",
  "Post-rev. protest (sqrt.)"
)

summary_fit_egypt = summary(full_fit_egypt)

extracted_egypt = extract(
  full_fit_egypt ,
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
rel.dep.effects.center.egypt = c(
  mean(dt_egypt_short$coledu),
  mean(dt_egypt_short$lowstat),
  mean(dt_egypt_short$coledu_lowstat)
)
rel.dep.effects.scale.egypt = c(
  sd(dt_egypt_short$coledu),
  sd(dt_egypt_short$lowstat),
  sd(dt_egypt_short$coledu_lowstat)
)

rel.dep.effects.std.egypt = (rel.dep.effects - t(array(rel.dep.effects.center.egypt, rev(
  dim(rel.dep.effects)
)))) / t(array(rel.dep.effects.scale.egypt, rev(dim(rel.dep.effects))))

logit.total.effect.egypt =
  sapply(1:dim(extracted_egypt$beta)[1], function(i) {
    as.matrix(rel.dep.effects.std.egypt) %*%
      extracted_egypt$beta[i, match(colnames(rel.dep.effects), pred.names.egypt)]
  })

row.names(logit.total.effect.egypt) =
  c(
    "no college education & high status",
    "college educated & high status",
    "no college education & low status",
    "college educated & low status"
  )

pdf(file = 'plots/egypt/total.effect_relative.deprivation_egypt.pdf',
    height = 5,
    width = 7.5)
plot(
  density(logit.total.effect.egypt[1, ]),
  xlim = c(
    min(logit.total.effect.egypt),
    max(logit.total.effect.egypt)
  ),
  main = "",
  col = NA,
  xlab = 'logit probability of recruitment',
  zero.line = FALSE
)
for (i in 1:dim(logit.total.effect.egypt)[1]) {
  d <- density(logit.total.effect.egypt[i, ])
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
  v = rowMeans(logit.total.effect.egypt),
  lty = 2,
  col = adjustcolor(c("blue", "darkgreen", "red", "black"), 0.5),
  lwd = 2
)
text(
  x = rowMeans(logit.total.effect.egypt),
  y = 0.8,
  labels = rownames(logit.total.effect.egypt),
  srt = 60,
  cex = 0.85,
  pos = 1,
  col = c("blue", "darkgreen", "red", "black")
)
dev.off()

# regression coefficients
pdf(file = 'plots/egypt/regression.coefficients_ind_egypt.pdf',
    height = 8,
    width = 12)
par(mfrow = c(3, 3))
for (i in 1:dim(extracted_all$beta)[2]) {
  plot(
    density(extracted_egypt$beta[, i]),
    xlim = c(ifelse(
      min(extracted_egypt$beta[, i]) < 0,
      min(extracted_egypt$beta[, i]) - 3,
      -3
    ),
    ifelse(
      max(extracted_egypt$beta[, i]) > 0,
      max(extracted_egypt$beta[, i]) + 3,
      3
    )),
    zero.line = FALSE,
    main = plot.names.egypt[i],
    col = NA,
    xlab = 'logit probability of recruitment'
  )
  d <- density(extracted_egypt$beta[, i])
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
    v = colMeans(extracted_egypt$beta)[i],
    lty = 2,
    col = adjustcolor("purple", 0.75),
    lwd = 2
  )
  abline(v = 0, lty = 1, lwd = 1)
  if (i == 1) {
    legend.var = c(
      paste("E[beta_", i, "] = ", round(mean(
        extracted_egypt$beta[, i]
      ), 3), sep = ""),
      paste("sd(beta_", i, ") = ", round(sd(
        extracted_egypt$beta[, i]
      ), 3), sep = ""),
      ifelse(
        mean(extracted_egypt$beta[, i]) > 0,
        paste("Pr(beta_", i, ">0) = ", round(
          sum(extracted_egypt$beta[, i] > 0) / length(extracted_egypt$beta[, i]),
          3
        ), sep = ""),
        paste("Pr(beta_", i, "<0) = ", round(
          sum(extracted_egypt$beta[, i] < 0) / length(extracted_egypt$beta[, i]),
          3
        ), sep = "")
      )
    )
  } else{
    legend.var = c(
      paste("E[beta_", i, "] = ", round(mean(
        extracted_egypt$beta[, i]
      ), 3), sep = ""),
      paste("sd(beta_", i, ") = ", round(sd(
        extracted_egypt$beta[, i]
      ), 3), sep = ""),
      ifelse(
        mean(extracted_egypt$beta[, i]) > 0,
        paste("Pr(beta_", i, ">0) = ", round(
          sum(extracted_egypt$beta[, i] > 0) / length(extracted_egypt$beta[, i]),
          3
        ), sep = ""),
        paste("Pr(beta_", i, "<0) = ", round(
          sum(extracted_egypt$beta[, i] < 0) / length(extracted_egypt$beta[, i]),
          3
        ), sep = "")
      ),
      paste("mean(X) = ", round(mean(
        unlist(dt_egypt_short[, pred.names[i], with = FALSE])
      ), 2)),
      paste("sd(X) = ", round(sd(
        unlist(dt_egypt_short[, pred.names[i], with = FALSE])
      ), 2))
    )
  }
  legend(
    ifelse(colMeans(extracted_egypt$beta)[i] > 0, "topleft", "topright"),
    legend = legend.var,
    cex = .9,
    bg = "white"
  )
}
dev.off()

# lambda density
pdf(file = 'plots/egypt/lambda_egypt.pdf',
    height = 5,
    width = 5)
plot(
  density(extracted_egypt$lambda),
  xlim = c(0, 1),
  zero.line = FALSE,
  main = "",
  col = NA,
  xlab = 'logit probability of recruitment'
)
polygon(
  density(extracted_egypt$lambda),
  main = '',
  xaxt = "n",
  yaxt = "n",
  ylab = "",
  xlab = "",
  border = NA,
  col = adjustcolor("orange", 0.75)
)
abline(
  v = mean(extracted_egypt$lambda),
  lty = 2,
  col = adjustcolor("purple", 0.75),
  lwd = 2
)
legend(
  "topright",
  legend = paste("MCMC mean:", round(mean(
    extracted_egypt$lambda
  ), 3)),
  lty = 2,
  lwd = 2,
  col = 'purple'
)
dev.off()