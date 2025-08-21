library(readr)
library(dplyr)
library(tibble)
library(metafor)
library(car)
library(performance)
library(ggplot2)

# 1.Convert SE in my original data into SD and prepare data
data <- read_csv("data.csv")
data <- data %>%
  filter(!is.na(i_bd_error_value), !is.na(data$c_bd_error_value))
data$i_bd_error_value <- data$i_bd_error_value * sqrt(data$i_sample_size)
data$c_bd_error_value <- data$c_bd_error_value * sqrt(data$c_sample_size)
write_csv(data, "data_sd.csv")

# Combine project member's data with my own data and save as combined_data.csv and import
data <- read_csv("combined_data.csv")

# 2.Calculate LRR (yi) and sampling variance (vi)
#Ensure numeric data are in numeric form
data <- data %>%
  mutate(
    i_bd_mean        = as.numeric(i_bd_mean),
    i_bd_error_value = as.numeric(i_bd_error_value),
    i_sample_size    = as.numeric(i_sample_size),
    c_bd_mean        = as.numeric(c_bd_mean),
    c_bd_error_value = as.numeric(c_bd_error_value),
    c_sample_size    = as.numeric(c_sample_size)
  )

#Epsilon: calculate 5th percentile of positive means for each metric type
eps_tbl <- data %>%
  group_by(bd_metric_type) %>%
  summarise(
    eps_i = quantile(i_bd_mean[i_bd_mean > 0], probs = 0.05, na.rm = TRUE), 
    eps_c = quantile(c_bd_mean[c_bd_mean > 0], probs = 0.05, na.rm = TRUE))

#Replace with lower bound if mean is too small
data <- data %>%
  left_join(eps_tbl, by = "bd_metric_type") %>%
  mutate(
    eps_i = ifelse(is.na(eps_i) | eps_i <= 0, 0.01, eps_i),
    eps_c = ifelse(is.na(eps_c) | eps_c <= 0, 0.01, eps_c),
    i_bd_mean_adj = ifelse(is.na(i_bd_mean) | i_bd_mean <= 0, eps_i, pmax(i_bd_mean, eps_i)),
    c_bd_mean_adj = ifelse(is.na(c_bd_mean) | c_bd_mean <= 0, eps_c, pmax(c_bd_mean, eps_c))
  )

#m1i = mean of intervention,sd1i = sd of intervention, n1i = sample size of intervention
es <- escalc(
  measure = "ROM",
  m1i = i_bd_mean_adj, sd1i = i_bd_error_value, n1i = i_sample_size,
  m2i = c_bd_mean_adj, sd2i = c_bd_error_value, n2i = c_sample_size,
  data = data)

#add calculated yi(LRR) and vi(sampling variance), and filter NA and keep only target functional groups
data <- bind_cols(data, es[, c("yi", "vi")]) %>%
  filter(
    !is.na(yi), !is.na(vi),
    !is.na(bd_metric_type),
    !is.na(functional_group),
    !is.na(scale_of_intervention)) %>%
  filter(functional_group %in% c("natural enemies","pests","pollinators")) %>%
  filter(bd_metric_type %in% c("abundance", "richness")) %>%
  filter(scale_of_intervention %in% c("small (2)", "medium (3)", "large (4)"))

write.csv(data, "effect_size.csv", row.names = FALSE)

#Factorise variables
data <- data %>%
  mutate(
    functional_group      = factor(functional_group,
                                   levels = c("pests", "natural enemies", "pollinators")),
    bd_metric_type        = factor(bd_metric_type),
    scale_of_intervention = factor(scale_of_intervention,
                                   levels = c("small (2)", "medium (3)", "large (4)")),
    study_ID              = factor(study_ID),
    effect_size_ID        = factor(effect_size_ID)
  )

data <- data %>% 
  filter(vi < 5000)

# 3. Baseline Model
m0 <- rma.mv(
  yi = yi, V = vi,
  random = ~ 1 | study_ID/effect_size_ID,
  method = "REML",
  data = data
)

summary(m0)

#Overall effect and prediction interval (back‐transformed to % change)
overall <- predict(m0, transf = exp)
pct_change  <- (overall$pred - 1) * 100
pct_ci_lb   <- (overall$ci.lb - 1) * 100
pct_ci_ub   <- (overall$ci.ub - 1) * 100
pct_pi_lb   <- (overall$pi.lb - 1) * 100
pct_pi_ub   <- (overall$pi.ub - 1) * 100
cat(sprintf("\n[Overall effect] %.1f%% (95%%CI: %.1f%%–%.1f%%; PI: %.1f%%–%.1f%%)\n",
            pct_change, pct_ci_lb, pct_ci_ub, pct_pi_lb, pct_pi_ub))

#I2 of baseline model
i2_decomp <- function(fit, mean_vi){
  sig2 <- fit$sigma2
  names(sig2) <- c("between_studies","within_study")
  tot <- sum(sig2) + mean_vi
  tibble(
    component = c("Level-3 (between studies)","Level-2 (within study)","Total I²"),
    I2 = 100 * c(sig2[1]/tot, sig2[2]/tot, sum(sig2)/tot)
  )
}

mean_vi <- mean(data$vi, na.rm = TRUE)
i2_m0   <- i2_decomp(m0, mean_vi)
print(i2_m0)
write_csv(i2_m0, "i2_components_m0.csv")

# 4. Check collinearity among moderators
mm <- lm(yi ~ functional_group + bd_metric_type + scale_of_intervention, data = data)
vif(mm)
check_collinearity(mm)

# 5.Three-level random-effects model with 2 moderators
m_main <- rma.mv(
  yi = yi, V = vi,
  mods = ~ functional_group + bd_metric_type,
  random = ~ 1 | study_ID/effect_size_ID,
  method = "REML",
  data = data
)

summary(m_main)
anova(m_main)

# Prediction interval
pred_main <- predict(m_main, transf = exp)
pred_tbl_main <- transform(pred_main,
                           pct = (pred-1)*100,
                           pct_ci.lb = (ci.lb-1)*100,
                           pct_ci.ub = (ci.ub-1)*100,
                           pct_pi.lb = (pi.lb-1)*100,
                           pct_pi.ub = (pi.ub-1)*100)

write_csv(as_tibble(pred_tbl_main), "prediction_interval_m_main.csv")

#Coefficient table
coef_main <- as.data.frame(coef(summary(m_main)))
coef_main$term <- rownames(coef_main); rownames(coef_main) <- NULL
coef_main$ci.lb  <- coef_main$estimate - 1.96*coef_main$se
coef_main$ci.ub  <- coef_main$estimate + 1.96*coef_main$se
coef_main$RR       <- exp(coef_main$estimate)
coef_main$RR_ci.lb <- exp(coef_main$ci.lb)
coef_main$RR_ci.ub <- exp(coef_main$ci.ub)
coef_main$pct      <- (coef_main$RR - 1)*100
coef_main$pct_lb   <- (coef_main$RR_ci.lb - 1)*100
coef_main$pct_ub   <- (coef_main$RR_ci.ub - 1)*100
coef_main <- coef_main[, c("term","estimate","se","zval","pval","ci.lb","ci.ub",
                           "RR","RR_ci.lb","RR_ci.ub","pct","pct_lb","pct_ub")]

write_csv(as_tibble(coef_main), "coefficients_m_main.csv")
print(coef_main)

# 6.forest plot
ord <- order(data$yi)
slabs <- with(data, paste0(study_ID, ":", effect_size_ID))[ord]
forest(m_main, order = ord, slab = slabs,
       xlab = "Log Response Ratio (LRR)", cex = 0.5,
       ylim = c(-1, length(slabs) + 3))

# 7.Publication bias
#Funnel plot
funnel(m_main)

#PET test (Precision Effect Test)
PET <- rma.mv(yi = yi, V = vi,
              mods = ~ I(sqrt(vi)),
              random = ~ 1 | study_ID/effect_size_ID,
              method = "REML", data = data)
summary(PET)

# PEESE test (Precision Effect Estimate with SE as explanatory variable)
PEESE <- rma.mv(yi = yi, V = vi,
                mods = ~ vi,
                random = ~ 1 | study_ID/effect_size_ID,
                method = "REML", data = data)
summary(PEESE)

# 8. Sensitivity analysis (leave-one-out)
studies <- unique(data$study_ID)
loo_results <- tibble(study = studies, intercept = NA_real_)
for (i in seq_along(studies)) {
  dat_sub <- subset(data, study_ID != studies[i])
  m_sub <- rma.mv(
    yi = yi, V = vi,
    mods = ~ functional_group + bd_metric_type,
    random = ~ 1 | study_ID/effect_size_ID,
    method = "REML", data = dat_sub
  )
  loo_results$intercept[i] <- coef(m_sub)[1]
}

write_csv(loo_results, "loo_main_model.csv")
print(loo_results)

#plot
ggplot(loo_results, aes(x = study, y = intercept, group = 1)) +
  geom_line() + geom_point() +
  labs(x = "Excluded study",
       y = "Overall log RR (leave-one-out)",
       title = "Leave-one-out sensitivity (main model)") +
  theme_bw()

# A1.confounding check by cross-tabs
xt_scale_x_fg     <- as.data.frame(table(data$scale_of_intervention, data$functional_group))
xt_scale_x_metric <- as.data.frame(table(data$scale_of_intervention, data$bd_metric_type))

write_csv(xt_scale_x_fg,     "appendix_A1_crosstab_scale_x_functional_group.csv")
write_csv(xt_scale_x_metric, "appendix_A1_crosstab_scale_x_bd_metric_type.csv")

# A2: overall effects by scale
small <- filter(data, scale_of_intervention == "small (2)")
medlarge <- filter(data, scale_of_intervention %in% c("medium (3)", "large (4)"))

fit_overall_by_group <- function(df) {
  fit  <- rma.mv(yi = yi, V = vi,
                 random = ~ 1 | study_ID/effect_size_ID,
                 method = "REML", data = df)
  pred <- predict(fit, transf = exp)
  tibble(
    k = nrow(df),
    RR = as.numeric(pred$pred),
    RR_ci.lb = as.numeric(pred$ci.lb), RR_ci.ub = as.numeric(pred$ci.ub),
    RR_pi.lb = as.numeric(pred$pi.lb), RR_pi.ub = as.numeric(pred$pi.ub),
    pct = (RR - 1) * 100,
    pct_ci.lb = (RR_ci.lb - 1) * 100, pct_ci.ub = (RR_ci.ub - 1) * 100,
    pct_pi.lb = (RR_pi.lb - 1) * 100, pct_pi.ub = (RR_pi.ub - 1) * 100)}

A2_small <- fit_overall_by_group(small) %>% 
  mutate(scale_group = "small (2)")
A2_medlarge <- fit_overall_by_group(medlarge) %>% 
  mutate(scale_group = "medium (3) + large (4)")

A2 <- bind_rows(A2_small, A2_medlarge) %>%
  select(scale_group, k, RR, RR_ci.lb, RR_ci.ub, RR_pi.lb, RR_pi.ub,
         pct, pct_ci.lb, pct_ci.ub, pct_pi.lb, pct_pi.ub)

write_csv(A2, "appendix_A2_scale_overall_effects.csv")
print(A2)
