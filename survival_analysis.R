library(data.table)
library(survminer)
library(survival)
library(ggplot2)

# read in data from authors
df <- fread("~/Downloads/Hakimi_data.txt")

# rename TTF months and TTF censor column
names(df)[26:27] <- c("TTF_months", "TTF")

# create additional column with definition from previous studies
df$pbrm1_new <- "other"
df[which(df$pbrm1_status == "lof"),]$pbrm1_new <- "lof"
df[which(df$pbrm1_mutation_types == "splice"),]$pbrm1_new <- "lof"

# create additional column with drug class assigned as authors did for their Cox PH models
df$drug_class_new <- df$drug_class
df[which(df$drug_class_new == "IO-IO"),]$drug_class_new <- "IO"

# recreate Figure 4B from authors paper
fit <- survfit(Surv(TTF_months, TTF) ~ pbrm1_status, data = df)
surv_pvalue(fit)
surv.plot <- ggsurvplot(fit,
                        conf.int = FALSE,
                        pval = TRUE,
                        risk.table = "absolute") +
  ggtitle("PBRM1 in ccRCC: TTF")
print(surv.plot)

# create figure using definition from previous studies
fit <- survfit(Surv(TTF_months, TTF) ~ pbrm1_new, data = df)
surv_pvalue(fit)
surv.plot <- ggsurvplot(fit, 
                        conf.int = FALSE,
                        palette = c("#552583", "#FDB927"),
                        pval = TRUE,
                        risk.table = "absolute") +
  ggtitle("PBRM1 in ccRCC: TTF")
print(surv.plot)

# recreate authors Cox PH results
surv_object <- Surv(df$TTF_months, df$TTF)
df$pbrm1_new <- factor(df$pbrm1_new, levels = c("other", "lof"))
df$pbrm1_status <- factor(df$pbrm1_status, levels = c("wt", "non_lof", "lof"))

# PBRM1 Status
fit.coxph <- coxph(surv_object ~ pbrm1_status, data = df)
fit.coxph

# TMB
fit.coxph <- coxph(surv_object ~ tmb, data = df)
fit.coxph

# Drug Class
fit.coxph <- coxph(surv_object ~ drug_class_new, data = df)
fit.coxph

# Multivariate
fit.coxph <- coxph(surv_object ~ pbrm1_status + tmb + drug_class_new, data = df)
fit.coxph


# Cox PH with definitions from previous studies
fit.coxph <- coxph(surv_object ~ pbrm1_new, data = df)
fit.coxph

fit.coxph <- coxph(surv_object ~ pbrm1_status + tmb + drug_class_new, data = df)
fit.coxph


# log TMB no longer significant predictor
# with padding to account for TMB = 0 cases
df$log_tmb <- log((df$tmb + min(df[which(df$tmb > 0), ]$tmb)))
surv_object <- Surv(df[which(df$log_tmb > 0), ]$TTF_months, df[which(df$log_tmb > 0), ]$TTF)
fit.coxph <- coxph(surv_object ~ log_tmb, data = df[which(df$log_tmb > 0), ])
fit.coxph

# and even when just removing those TMB = 0 samples
df$log_tmb <- log(df$tmb)
surv_object <- Surv(df[which(df$log_tmb > 0), ]$TTF_months, df[which(df$log_tmb > 0), ]$TTF)
fit.coxph <- coxph(surv_object ~ log_tmb, data = df[which(df$log_tmb > 0), ])
fit.coxph

# TTF based on TMB tertiles 
tertiles <- quantile(df$tmb, c(0.33, 0.66))
df$tmb_tertile <- 1
df[which(df$tmb > tertiles[1] & df$tmb <= tertiles[2]), ]$tmb_tertile <- 2
df[which(df$tmb > tertiles[2]), ]$tmb_tertile <- 3

# plot TTF curves for TMB tertiles
# (middle tertile performs worst suggesting continuous TMB is not a thing in ccRCC)
fit <- survfit(Surv(TTF_months, TTF) ~ tmb_tertile, data = df)
surv_pvalue(fit)
surv.plot <- ggsurvplot(fit, 
                        conf.int = FALSE,
                        palette = c("#552583", "#FDB927", "black"),
                        pval = TRUE,
                        risk.table = "absolute") +
  ggtitle("PBRM1 in ccRCC: TTF")
print(surv.plot)


# by removing PBRM1 from the TMB equation, there is no survival effect based on TMB
df$tmb_tertile_pbrm1 <- paste0(df$pbrm1_new, " - ", df$tmb_tertile)

fit <- survfit(Surv(TTF_months, TTF) ~ tmb_tertile_pbrm1, 
               data = df[which(df$tmb_tertile_pbrm1 %in% c("other - 1", "other - 3")), ])
surv_pvalue(fit)
surv.plot <- ggsurvplot(fit, 
                        conf.int = FALSE,
                        palette = c("#552583", "#FDB927"),
                        pval = TRUE,
                        risk.table = "absolute") +
  ggtitle("PBRM1 in ccRCC: TTF")
print(surv.plot)


# same applies for quartiles
quartiles <- quantile(df$tmb, c(0.25, 0.5, 0.75))
df$tmb_quartile <- 1
df[which(df$tmb > quartiles[1] & df$tmb <= quartiles[2]), ]$tmb_quartile<- 2
df[which(df$tmb > quartiles[2] & df$tmb <= quartiles[3]), ]$tmb_quartile<- 3
df[which(df$tmb > quartiles[3]), ]$tmb_quartile <- 4
df$tmb_quartile_pbrm1 <- paste0(df$pbrm1_new, " - ", df$tmb_quartile)

fit <- survfit(Surv(TTF_months, TTF) ~ tmb_quartile_pbrm1, 
               data = df[which(df$tmb_quartile_pbrm1 %in% c("other - 1", "other - 4")), ])
surv_pvalue(fit)
surv.plot <- ggsurvplot(fit, 
                        conf.int = FALSE,
                        palette = c("#552583", "#FDB927"),
                        pval = TRUE,
                        risk.table = "absolute") +
  ggtitle("PBRM1 in ccRCC: TTF")
print(surv.plot)

