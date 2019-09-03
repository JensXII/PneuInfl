################
### Analyses ###
################

#' GLM backward selection.
#'
#' @param regression an glm regression object.
#' @param min minimum significance (default = 0.1)
#' @param verbose show stepwise summaries (default = FALSE)
#' @return Return regression object with selected variables.
#' @examples
#'
bw.glm <- function(regression, min, verbose = FALSE) {
  res <- regression
  covars <- names(summary(res)$coeff[-1,4])
  while ((length(covars) >= 0) & (max(summary(res)$coeff[-1,4]) >= min)) {
    toselect <- summary(res)$coeff[-1,4] < max(summary(res)$coeff[-1,4])
    # select sig. covars
    covars <- names(toselect)[toselect == TRUE]
    if (length(covars) == 0) break
    # formula with only sig covars
    formula <- as.formula(paste(res$formula[[2]], "~", paste(c("1", covars), collapse = " + ")))
    res <- update(res, formula)
    if (verbose) print(summary(res))
  }
  return(res)
}

### Data ###

load("data/InflData.RData")
load("data/PneuData.RData")
PneuData <- subset(PneuData, Pneu == 1)
PneuInfl.data <- CreateData(1, 15)

PneuInfl.data$sin52 <- sin((2*pi/(365.25)) * as.numeric(PneuInfl.data$date))
PneuInfl.data$cos52 <- cos((2*pi/(365.25)) * as.numeric(PneuInfl.data$date))
PneuInfl.data$sin26 <- sin((4*pi/(365.25)) * as.numeric(PneuInfl.data$date))
PneuInfl.data$cos26 <- cos((4*pi/(365.25)) * as.numeric(PneuInfl.data$date))
PneuInfl.data$season <- as.numeric(format(PneuInfl.data$date, "%Y")) - (as.numeric(format(PneuInfl.data$date, "%m")) <= 9)


### P(Infl | Pneu.prev) = NInflPneu.prev / NPneu.prev
plot(PneuInfl.data$NPneu.prev, PneuInfl.data$NInflPneu.prev)

# Rate
aggregate(NInflPneu.prev ~ 1 , data = PneuInfl.data, sum) / aggregate(NPneu.prev ~ 1 , data = PneuInfl.data, sum)

Infl.res <- bw.glm(glm(NInflPneu.prev ~ NPneu.prev + NInfl + NInfl.prev + temp, offset = log(NPneu.prev), poisson, data = subset(PneuInfl.data, NPneu.prev > 0)), 0.1, verbose = TRUE)
summary(Infl.res)
exp(cbind(coefficients(Infl.res), confint(Infl.res)))

Infl.res <- bw.glm(glm(NInflPneu.prev ~ sin52 + cos52 + sin26 + cos26, offset = log(NPneu.prev), poisson, data = subset(PneuInfl.data, NPneu.prev > 0)), 0.1)
Infl.res <- glm(NInflPneu.prev ~ sin52 + cos52 + sin26 + cos26, offset = log(NPneu.prev), poisson, data = subset(PneuInfl.data, NPneu.prev > 0))
summary(Infl.res)
exp(cbind(coefficients(Infl.res), confint(Infl.res)))


### Over calendar time

Infl <- glm(NInfl ~ NInfl.prev + factor(season) + sin52 + cos52 + sin26 + cos26, offset = log(N), poisson, data = PneuInfl.data)
Infl <- glm(NInfl ~ splines::bs(date, df = 10) + sin52 + cos52 + sin26 + cos26, offset = log(N), poisson, data = PneuInfl.data)
summary(Infl)
PneuInfl.data$EnoPneu <- exp(predict.glm(Infl, se.fit=TRUE)$fit)

Infl <- glm(NInfl ~ NPneu.prev + NInfl.prev + factor(season) + sin52 + cos52 + sin26 + cos26, offset = log(N), quasipoisson, data = PneuInfl.data)
# Infl <- glm(NInfl ~ NPneu + NInfl.prev + splines::bs(date, df = 10) + sin52 + cos52 + sin26 + cos26, offset = log(N), poisson, data = PneuInfl.data)
summary(Infl)
exp(cbind(coefficients(Infl), confint(Infl)))

PneuInfl.data$Epred <- exp(predict.glm(Infl, se.fit=TRUE)$fit)
PneuInfl.data.B <- PneuInfl.data
PneuInfl.data.B[, c("NInfl.prev", "NPneu.prev")] <- 0
PneuInfl.data$EB <- exp(predict.glm(Infl, newdata=PneuInfl.data.B, se.fit=TRUE)$fit)
rm(PneuInfl.data.B)

library(ggplot2)
# png(filename = paste0("output/Infl.png"), width = 1280, height = 768, units = "px", pointsize = 12)

ggplot(PneuInfl.data, aes(x = date)) +
  geom_point(aes(date, NInfl, colour = "NInfl"), size = 0.2) +
  geom_line(aes(date, EB, colour = "EB")) +
  geom_line(aes(date, Epred, colour = "Epred")) +
  geom_line(aes(date, EnoPneu, colour = "EnoPneu")) +
  ggtitle("Influenza") + theme(plot.title = element_text(hjust = 0.5)) +
  xlab("date") + ylab("Number of cases per day") +
  scale_colour_manual("",
                      breaks = c("NInfl", "EnoPneu", "EB", "Epred"),
                      values = c("NInfl" = "black", "EnoPneu" = "green", "EB" = "red", "Epred" = "blue"),
                      guide = guide_legend(override.aes = list(linetype = c("blank", rep("solid", 3)),
                                                               shape = c(16, rep(NA, 3))))) +
  theme(legend.position = "bottom")

# dev.off()

### Pneu

### P(Pneu | Infl.prev) = NPneuInfl.prev / NInfl.prev
plot(PneuInfl.data$NInfl.prev, PneuInfl.data$NPneuInfl.prev)

# Rate
aggregate(NPneuInfl.prev ~ 1 , data = PneuInfl.data, sum) / aggregate(NInfl.prev ~ 1 , data = PneuInfl.data, sum)

Pneu.res <- bw.glm(glm(NPneuInfl.prev ~ NPneu + NPneu.prev + NInfl + NInfl.prev, offset = log(NInfl.prev), poisson, data = subset(PneuInfl.data, NInfl.prev > 0)), 0.1)
summary(Pneu.res)
exp(cbind(coefficients(Pneu.res), confint(Pneu.res)))

### Over calendat time

Pneu <- glm(NPneu ~ NPneu.prev + factor(season) + sin52 + cos52 + sin26 + cos26, offset = log(N), poisson, data = PneuInfl.data)
# Pneu <- glm(NPneu ~ splines::bs(date, df = 10) + sin52 + cos52 + sin26 + cos26, offset = log(N), poisson, data = PneuInfl.data)
PneuInfl.data$EnoInfl <- exp(predict.glm(Pneu, se.fit=TRUE)$fit)

Pneu <- glm(NPneu ~ NPneu.prev + NInfl.prev + factor(season) + sin52 + cos52 + sin26 + cos26, offset = log(N), quasipoisson, data = PneuInfl.data)
# Pneu <- glm(NPneu ~ NPneu.prev + NInfl.prev + splines::bs(date, df = 10) + sin52 + cos52 + sin26 + cos26, offset = log(N), poisson, data = PneuInfl.data)
summary(Pneu)
exp(cbind(coefficients(Pneu), confint(Pneu)))

PneuInfl.data$Epred <- exp(predict.glm(Pneu, se.fit=TRUE)$fit)
PneuInfl.data.B <- PneuInfl.data
PneuInfl.data.B$NInfl.prev <- 0
# PneuInfl.data.B[, c("NInfl.prev", "NPneu.prev")] <- 0
PneuInfl.data$EB <- exp(predict.glm(Pneu, newdata=PneuInfl.data.B, se.fit=TRUE)$fit)
rm(PneuInfl.data.B)

library(ggplot2)
# png(filename = paste0("output/Pneu.png"), width = 1280, height = 768, units = "px", pointsize = 12)

ggplot(PneuInfl.data, aes(x = date)) +
  geom_point(aes(date, NPneu, colour = "NPneu"), size = 0.2) +
  geom_line(aes(date, EB, colour = "EB")) +
  geom_line(aes(date, Epred, colour = "Epred")) +
  geom_line(aes(date, EnoInfl, colour = "EnoInfl")) +
  ggtitle("Pneumonia") + theme(plot.title = element_text(hjust = 0.5)) +
  xlab("date") + ylab("Number of cases per day") +
  scale_colour_manual("",
                      breaks = c("NPneu", "EnoInfl", "EB", "Epred"),
                      values = c("NPneu" = "black", "EnoInfl" = "green", "EB" = "red", "Epred" = "blue"),
                      guide = guide_legend(override.aes = list(linetype = c("blank", rep("solid", 3)),
                                                               shape = c(16, rep(NA, 3))))) +
  theme(legend.position = "bottom")

# dev.off()




#summary(Pneu.res)$dispersion
#sum(Pneu.res[["weights"]]*Pneu.res[["residuals"]]^2)/df.residual(Pneu.res)
#sum(residuals(Pneu.res, type = "deviance")^2)/df.residual(Pneu.res)
