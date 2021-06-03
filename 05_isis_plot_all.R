rm(list = ls())
gc()
options(scipen = 999)
# # # # #
library(data.table)
library(rstan)
library(xtable)

# PRINT AND PLOT RESULTS

# load estimation and data files for all countries
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

load(file = 'data/analysis/cleaned.all.short.RData')

summary_fit_all = summary(full_fit_all)
# results
print(
  full_fit_all,
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

print(full_fit_all, pars = c('beta'))
print(full_fit_all, pars = c('phi'))
print(full_fit_all, pars = c('psi'))
print(full_fit_all, pars = c('eta'))
print(full_fit_all, pars = c('lambda'))
print(full_fit_all, pars = c('tau'))
print(full_fit_all, pars = c('gamma'))
print(full_fit_all, pars = c('rho'))
print(full_fit_all, pars = c('rho_star'))
print(full_fit_all, pars = c('lp__'))

summary_fit_all = summary(full_fit_all)

################# MAIN PLOTS

# relative deprivation effects
extracted_all = extract(
  full_fit_all ,
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
rel.dep.effects.center = c(0.2, 0.23, 0.02)
rel.dep.effects.scale = c(0.4, 0.42, 0.16)
rel.dep.effects.std = (rel.dep.effects - t(array(rel.dep.effects.center, rev(
  dim(rel.dep.effects)
)))) / t(array(rel.dep.effects.scale, rev(dim(rel.dep.effects))))

logit.total.effect =
  sapply(1:dim(extracted_all$beta)[1], function(i) {
    as.matrix(rel.dep.effects.std) %*%
      extracted_all$beta[i, match(colnames(rel.dep.effects), pred.names)]
  })
row.names(logit.total.effect) =
  c(
    "no college education & high status",
    "college educated & high status",
    "no college education & low status",
    "college educated & low status"
  )

pdf(file = 'plots/all/total.effect_relative.deprivation_all.pdf',
    height = 5,
    width = 7.5)
plot(
  density(logit.total.effect[1,]),
  xlim = c(min(logit.total.effect), max(logit.total.effect)),
  main = "",
  col = NA,
  xlab = 'logit probability of recruitment',
  zero.line = FALSE
)
for (i in 1:dim(logit.total.effect)[1]) {
  d <- density(logit.total.effect[i,])
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
  v = rowMeans(logit.total.effect),
  lty = 2,
  col = adjustcolor(c("blue", "darkgreen", "red", "black"), 0.5),
  lwd = 2
)
text(
  x = rowMeans(logit.total.effect),
  y = 7.5,
  labels = rownames(logit.total.effect),
  srt = 60,
  cex = 0.85,
  pos = 1,
  col = c("blue", "darkgreen", "red", "black")
)
dev.off()

# b) regression coefficients

plot.names = c(
  "Intercept",
  "College edu.",
  "Age",
  "Age^2",
  "Married",
  "Student",
  "Low status",
  "College edu.*Low status"
)

pdf(file = 'plots/all/regression.coefficients_ind_all.pdf',
    height = 8,
    width = 12)
par(mfrow = c(3, 3))
for (i in 1:dim(extracted_all$beta)[2]) {
  plot(
    density(extracted_all$beta[, i]),
    xlim = c(ifelse(
      min(extracted_all$beta[, i]) < 0, min(extracted_all$beta[, i]) - 3,-3
    ),
    ifelse(
      max(extracted_all$beta[, i]) > 0, max(extracted_all$beta[, i]) + 3, 3
    )),
    zero.line = FALSE,
    main = plot.names[i],
    col = NA,
    xlab = 'logit probability of recruitment'
  )
  d <- density(extracted_all$beta[, i])
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
    v = colMeans(extracted_all$beta)[i],
    lty = 2,
    col = adjustcolor("purple", 0.75),
    lwd = 2
  )
  abline(v = 0, lty = 1, lwd = 1)
  if (i == 1) {
    legend.var = c(
      paste("E[beta_", i, "] = ", round(mean(
        extracted_all$beta[, i]
      ), 3), sep = ""),
      paste("sd(beta_", i, ") = ", round(sd(
        extracted_all$beta[, i]
      ), 3), sep = ""),
      ifelse(
        mean(extracted_all$beta[, i]) > 0,
        paste("Pr(beta_", i, ">0) = ", round(
          sum(extracted_all$beta[, i] > 0) / length(extracted_all$beta[, i]),
          3
        ), sep = ""),
        paste("Pr(beta_", i, "<0) = ", round(
          sum(extracted_all$beta[, i] < 0) / length(extracted_all$beta[, i]),
          3
        ), sep = "")
      )
    )
  } else{
    legend.var = c(
      paste("E[beta_", i, "] = ", round(mean(
        extracted_all$beta[, i]
      ), 3), sep = ""),
      paste("sd(beta_", i, ") = ", round(sd(
        extracted_all$beta[, i]
      ), 3), sep = ""),
      ifelse(
        mean(extracted_all$beta[, i]) > 0,
        paste("Pr(beta_", i, ">0) = ", round(
          sum(extracted_all$beta[, i] > 0) / length(extracted_all$beta[, i]),
          3
        ), sep = ""),
        paste("Pr(beta_", i, "<0) = ", round(
          sum(extracted_all$beta[, i] < 0) / length(extracted_all$beta[, i]),
          3
        ), sep = "")
      ),
      paste("mean(X) = ", round(mean(
        unlist(dt_all_short[, pred.names[i], with = FALSE])
      ), 2)),
      paste("sd(X) = ", round(sd(
        unlist(dt_all_short[, pred.names[i], with = FALSE])
      ), 2))
    )
  }
  legend(
    ifelse(colMeans(extracted_all$beta)[i] > 0, "topleft", "topright"),
    legend = legend.var,
    cex = 1,
    bg = "white"
  )
}
dev.off()

# c) lambda density
pdf(file = 'plots/all/lambda_all.pdf',
    height = 5,
    width = 5)
plot(
  density(extracted_all$lambda),
  xlim = c(0, 1),
  zero.line = FALSE,
  main = "",
  col = NA,
  xlab = 'logit probability of recruitment'
)
polygon(
  density(extracted_all$lambda),
  main = '',
  xaxt = "n",
  yaxt = "n",
  ylab = "",
  xlab = "",
  border = NA,
  col = adjustcolor("orange", 0.75)
)
abline(
  v = mean(extracted_all$lambda),
  lty = 2,
  col = adjustcolor("purple", 0.75),
  lwd = 2
)
legend(
  "topleft",
  legend = paste("MCMC mean:", round(mean(
    extracted_all$lambda
  ), 3)),
  lty = 2,
  lwd = 2,
  col = 'purple'
)
dev.off()



# d) predicted probabilities of all possible profiles on a distirbution
profiles =
  expand.grid(
    coledu = unique(dt_all_short$coledu),
    age = quantile(dt_all_short$age, probs = c(0.1, 0.25, 0.5, 0.75, 0.9)),
    married = unique(dt_all_short$married),
    student = unique(dt_all_short$student),
    lowstat = unique(dt_all_short$lowstat),
    governorate = unique(paste(dt_all_short$governorate))
  )
profiles$age2 = profiles$age ^ 2
profiles$coledu_lowstat = profiles$coledu * profiles$lowstat
profiles$country = unlist(sapply(profiles$governorate, function(x) {
  strsplit(as.character(unlist(x)), " - ")[[1]][1]
}))
profiles$governorate = unlist(sapply(profiles$governorate, function(x) {
  strsplit(as.character(unlist(x)), " - ")[[1]][2]
}))

# calculate total effect

temp = cbind(
  temp.0 = mean(extracted_all$beta[, which(pred.names == "intercept")]),
  temp.1 =  mean(extracted_all$beta[, which(pred.names == "coledu")]) *
    (profiles$coledu - mean(dt_all_short$coledu)) / sd(dt_all_short$coledu),
  temp.2 =  mean(extracted_all$beta[, which(pred.names == "age")]) * (profiles$age -
                                                                        mean(dt_all_short$age)) / sd(dt_all_short$age),
  temp.3 =  mean(extracted_all$beta[, which(pred.names == "married")]) *
    (profiles$married - mean(dt_all_short$married)) / sd(dt_all_short$married),
  temp.4 =  mean(extracted_all$beta[, which(pred.names == "student")]) *
    (profiles$student - mean(dt_all_short$student)) / sd(dt_all_short$student),
  temp.5 =  mean(extracted_all$beta[, which(pred.names == "lowstat")]) *
    (profiles$lowstat - mean(dt_all_short$lowstat)) / sd(dt_all_short$lowstat),
  temp.6 =  mean(extracted_all$beta[, which(pred.names == "age2")]) * (profiles$age ^
                                                                         2 - mean(dt_all_short$age ^ 2)) / sd(dt_all_short$age ^ 2),
  temp.7 =  mean(extracted_all$beta[, which(pred.names == "coledu_lowstat")]) *
    (
      profiles$coledu_lowstat - mean(dt_all_short$coledu_lowstat)
    ) / sd(dt_all_short$coledu_lowstat),
  temp.8 =  mean(sqrt(1 / extracted_all$tau)) * colMeans(extracted_all$gamma[, match(
    paste(profiles$country, profiles$governorate, sep = " - "),
    levels(dt_all_short$governorate)
  )]),
  temp.9 =  mean(sqrt(1 / extracted_all$tau_eta)) * colMeans(extracted_all$eta[, match(profiles$country, levels(dt_all_short$country))])
)
profiles$total_effect = rowSums(temp)

pdf(file = 'plots/all/predicted.probs_all.pdf',
    height = 5,
    width = 7.5)
d = density(profiles$total_effect)
plot(
  d,
  ylim = c(0, max(d$y) + 0.1),
  zero.line = FALSE,
  xlim = c(min(d$x), max(d$x)),
  main = "",
  col = NA,
  xlab = 'logit probability of recruitment'
)
polygon(
  d,
  xaxt = "n",
  yaxt = "n",
  ylab = "",
  xlab = "",
  border = "orange",
  col = adjustcolor("orange", 0.75)
)
profiles = profiles[order(profiles$total_effect),]
quantiles = c(0, 0.1, 0.5, 0.9, 1)
rank = quantile(order(profiles$total_effect), probs = quantiles)
temp = profiles[, c("coledu",
                    "age",
                    "married",
                    "student",
                    "lowstat",
                    "country",
                    "governorate")]
temp$coledu = ifelse(temp$coledu == 1, "college edu.", "no college edu.")
temp$age = paste("age:", temp$age)
temp$married = ifelse(temp$married == 1, "married", "not married")
temp$student = ifelse(temp$student == 1, "student", "not a student")
temp$lowstat = ifelse(temp$lowstat == 1, "low status", "high status")
for (r in 1:length(rank)) {
  points(
    x = profiles$total_effect[rank[r]],
    y = d$y[which(abs(d$x - profiles$total_effect[rank[r]]) == min(abs(d$x -
                                                                         profiles$total_effect[rank[r]])))],
    col = 'purple',
    pch = c(0:4)[r],
    lwd = 2
  )
  text(
    x = profiles$total_effect[rank[r]],
    y = 0.01 + d$y[which(abs(d$x - profiles$total_effect[rank[r]]) ==
                           min(abs(d$x - profiles$total_effect[rank[r]])))],
    col = 'black',
    cex = 0.7,
    pos = 3,
    labels = apply(temp, 1, paste0, collapse = "\n")[rank[r]]
  )
}
segments(
  x0 = median(profiles$total_effect),
  x1 = median(profiles$total_effect),
  y0 = 0,
  y1 = d$y[which(abs(d$x - profiles$total_effect[rank[which(quantiles ==
                                                              0.5)]]) == min(abs(d$x - profiles$total_effect[rank[which(quantiles == 0.5)]])))],
  col = 'purple',
  lty = 2,
  lwd = 2
)
legend(
  "topleft",
  legend = c(
    "quantile: 0",
    "quantile: 0.1",
    "quantile: 0.5",
    "quantile: 0.9",
    "quantile: 1"
  ),
  pch = c(0:4),
  col = 'purple',
  lwd = 2,
  lty = NA
)
dev.off()

print(xtable(cbind(quantile = rev(quantiles),
                   profiles[rev(rank), c("coledu",
                                         "age",
                                         "married",
                                         "student",
                                         "lowstat",
                                         "country",
                                         "governorate")])),
      include.rownames = FALSE)
print(xtable(cbind(profiles[rev(order(profiles$total_effect))[1:10],
                            c("coledu",
                              "age",
                              "married",
                              "student",
                              "lowstat",
                              "country",
                              "governorate")],
                   pred.prob = rev(
                     exp(profiles$total_effect) / (1 + exp(profiles$total_effect))
                   )[1:10]), digits = 5), include.rownames = TRUE)