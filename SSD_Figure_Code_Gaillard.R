

## This code produces Figures 1, 2 and 3 of the paper by Gaillard et al. in
## response to Tombak et al. (2024)


## ------------------------------------
##
## Load libraries, functions and data
##
## ------------------------------------


library(ggplot2)
library(ggrepel)
library(MASS)
library(scales)
library(mgcv)
library(gtools)

## Updated data set corrected for errors
load("ssd_ngy.RData")

logit <- function(x) {log(x / (1 - x))}
ilogit <- function(x) {exp(x) / (1 + exp(x))}


## ---------------------
##
## Code to plot Fig. 1
##
## ---------------------

## Fit allometric and isometric models to the standard deviation in MALE body
## size across mammals

lm1 <- lm(log(SDmassM) ~ log(massM), data = ssd9)
summary(lm1)
lm2 <- lm(log(SDmassM) ~ offset(log(massM)), data = ssd9)

## Compute predict values, 95% CI and format the results for convienient plotting
new <- data.frame(
    massM = seq(min(ssd$massM, na.rm = T), max(ssd$massM, na.rm = T), le = 1000)
)
new <- cbind(new, as.data.frame(predict(lm1, new, se.fit = T)))
new$lo <- exp(new$fit - 1.96 * new$se.fit)
new$up <- exp(new$fit + 1.96 * new$se.fit)
new$fit <- exp(new$fit)
new$off <- exp(predict(lm2, new))

## Plot observed points and fitted allometric relationships for MALES
fig_comb <-
ggplot(ssd9, aes(x =  massM, y =  SDmassM)) +
    geom_point(alpha = 0.25, size =  2, colour = "black") +
    geom_line(data =  new, aes(y =  fit, x = massM)) +
    geom_line(data =  new, aes(y =  off, x = massM), linetype =  2) +
    geom_ribbon(data =  new, aes(ymin = lo, ymax = up, x = massM),
        fill = "black", alpha = 0.3, inherit.aes = FALSE
    ) +
    theme_bw(base_size =  20) +
    scale_x_log10(
        breaks = trans_breaks("log10", function(x) 10^x),
        labels = trans_format("log10", math_format(10^.x))
    ) +
    scale_y_log10(
        breaks = trans_breaks("log10", function(x) 10^x),
        labels = trans_format("log10", math_format(10^.x))
    ) +
    annotation_logticks(sides="b") +
    xlab("Species Body Mass (g)") +
    ylab("SD(Body Mass)") +
    theme(legend.position="none")

## Fit allometric and isometric models to the standard deviation in FEMALE body
## size across mammals

lm1 <- lm(log(SDmassF) ~ log(massF), data = ssd9)
summary(lm1)
lm2 <- lm(log(SDmassF) ~ offset(log(massF)), data = ssd9)

## Compute predicted values and associated 95% CI, and format the results for
## convienient plotting
new <- data.frame(
    massF = seq(min(ssd$massF, na.rm = T), max(ssd$massF, na.rm = T), le = 1000)
)
new <- cbind(new, as.data.frame(predict(lm1, new, se.fit = T)))
new$lo <- exp(new$fit - 1.96 * new$se.fit)
new$up <- exp(new$fit + 1.96 * new$se.fit)
new$fit <- exp(new$fit)
new$off <- exp(predict(lm2, new))

## Combine MALE plot with observed points and fitted allometric relationships
## for FEMALES
fig_comb <- fig_comb +
    geom_point(data = ssd9, aes(x =  massF, y =  SDmassF),
        alpha = 0.25, size =  2, colour = "#E69F00") +
    geom_line(data =  new, aes(y =  fit, x = massF), colour = "#E69F00") +
    geom_line(data =  new, aes(y =  off, x = massF), linetype =  2) +
    geom_ribbon(data =  new, aes(ymin = lo, ymax = up, x = massF), fill = "#E69F00",
        alpha = 0.3, inherit.aes = FALSE)

fig_comb

## Not run
##
## dev.copy2pdf(file =  "Fig_combined.pdf")
## dev.off()


## ---------------------
##
## Code to plot Fig. 2
##
## ---------------------

## Convert the male dimorphic factor varible into 0/1 numbers
ssd9$is.Mdimo <- as.numeric(ssd9$massDimorphism == "Male-Biased Dimorphic")
## Compute log10 of male body mass
ssd9$log10_massM <- log10(ssd9$massM)

## Run logistic regression to describe the link between male body size and the probability for the species to show a signitifcant male-biased sexual dimorphism
glm1 <- glm(is.Mdimo ~ log10_massM, data =  ssd9, family = "binomial")
summary(glm1)

## Compute predicted values and associated 95% CI, and format the results for
## convienient plotting
new <- data.frame(log10_massM = seq(min(ssd9$log10_massM, na.rm = T),
    max(ssd9$log10_massM, na.rm = T), le = 500))
new <- cbind.data.frame(new, predict(glm1, new = new, se.fit = TRUE))
new$up <-  ilogit(new$fit + 1.96 * new$se.fit)
new$lo <-  ilogit(new$fit - 1.96 * new$se.fit)
new$fit <- ilogit(new$fit)
new$massM <- 10^new$log10_massM
head(new)
dim(new)

## Check raw values for the proportion of male dimorphic species by deciles of
## body mass distribution
moy <- data.frame(
    obs = tapply(ssd9$is.Mdimo, quantcut(ssd9$massM, 20), mean, na.rm =  T),
    mass =  tapply(ssd9$massM, quantcut(ssd9$massM, 20), mean, na.rm =  T)
)

## Plot observed data and predicted values from logistic model
p <- ggplot(ssd9, aes(massM, is.Mdimo)) +
    geom_point(position = position_jitter(height = 0, width=0.03), alpha =  0.1, size = 3) +
    xlab("Species Body Mass (g)") + ylab("Pr(Male dimorphic)") +
    theme_bw(base_size = 20) +
    scale_x_log10(
        breaks = trans_breaks("log10", function(x) 10^x),
        labels = trans_format("log10", math_format(10^.x))
    ) +
    annotation_logticks(sides="b") +
    geom_ribbon(data = new, aes(ymin = lo, ymax = up, x = massM), fill = "#E69F00",
        alpha = 0.3, inherit.aes = FALSE) +
    geom_line(data = new, aes(x = massM, y =  fit), colour = "#E69F00", size =  2) +
    geom_line(data = new, aes(x = massM, y =  fit)) +
    geom_hline(
        yintercept = 0.5,
        linetype="dotted"
    )
p

## Not run
##
## dev.copy2pdf(file =  "Fig_prob_all.pdf")
## dev.off()


## ---------------------
##
## Code to plot Fig. 3
##
## ---------------------


## Calculate the proportion of male-dimorphic species by Order
dimo.table <- data.frame(
    order = levels(ssd9$Order),
    N   = as.numeric(table(ssd9$Order)),
    fit = as.numeric(table(ssd9$is.Mdimo, ssd9$Order)[1, ] / table(ssd9$Order)),
    lo = NA,
    up = NA,
    loc = 1:17
)
for(i in 1:dim(dimo.table)[1]) {
    dimo.table[i, 4:5] <- binom.test(
        table(ssd9$is.Mdimo, ssd9$Order)[1, i],
        table(ssd9$Order)[i])$conf.int
    ## print(binom.test(
    ##     table(ssd9$is.Mdimo, ssd9$Order)[1, i],
    ##     table(ssd9$Order)[i])$conf.int)
}
dimo.table
o <- order(dimo.table$fit)
## Not run
##
## Export table in CSV format
## write.csv2(dimo.table[o, ], file =  "Table1.csv")

dimo.table <- dimo.table[o, ]
dimo.table$loc <- 1:17

## Plot observed values and associated 95% confidence intervals
g <- ggplot(data = dimo.table, aes(loc, fit))+
    geom_errorbar(aes(ymax =  up, ymin =  lo), width =  0, col = "grey42") +
    geom_point(
        color='white',
        size =  7,
        aes(fill = fit)
    ) +
    geom_point(
        color='black',
        shape = 21,
        size =  4,
        aes(fill = fit)
    ) +
    theme_bw(base_size =  20) +
    theme(
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        aspect.ratio=1
    ) +
    geom_text_repel(
        aes(label = order),
        size = 5,
        box.padding   = unit(0.1, "lines"),
        point.padding = unit(15, "lines"),
        force = 50,
        segment.size = 0.2,
        segment.color = "grey75"
    ) +
    xlab("") + ylab("Pr(Male dimorphic)") +
    theme(legend.position = "none")
g

## Not run
##
## dev.copy2pdf(file =  "Fig_prob_all.pdf")
## dev.off()

