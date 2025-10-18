
#install.packages("fixest")
#install.packages("data.table")
library(fixest)
library(data.table)
set.seed(20200403)

# 100,000 people with differing levels of covid symptoms
N_people = 100000
df = data.table(person = 1:N_people)



# Potential outcomes
# Y(O): life-span if no vent
# Y(0) - normal distribution with a mean of 9.4 and standard deviation of 4;
#        replace negative values with 0.

df[, y0 := rnorm(N_people, 9.4, 4)]
df[y0 < 0, y0 := 0]



# Y(1): life span if assigned to vents
#       follows a normal distribution with a mean of 10 and standard deviation of 4;
#       replace negative values with 0.

df[, y1 := rnorm(N_people, 10, 4)]
df[y1 < 0, y1 := 0]

# Define individual treatment effect
df[, delta := y1 - y0]

# Perfect doctor assigns vents (the treatment) only to those who benefit
df[, vents := (delta > 0)]

# Calculate all aggregate Causal Parameters (ATE, ATT, ATU)
ate = df[, mean(delta)]
att = df[vents == TRUE, mean(delta)]
atu = df[vents == FALSE, mean(delta)]

cat(sprintf("ATE = %.03f\n", ate))
cat(sprintf("ATT = %.03f\n", att))
cat(sprintf("ATU = %.03f\n", atu))

# Use the switching equation to select realised outcomes from potential outcomes based on treatment assignment given by the Perfect Doctor
df[, y := vents * y1 + (1 - vents) * y0]

# Calculate E[Y(0)] for vent group and no vent group so we can calculate selection bias
ey01 = df[vents == TRUE, mean(y0)]
ey00 = df[vents == FALSE, mean(y0)]

# Calculate selection bias based on the previous conditional expectations
selection_bias = (ey01 - ey00)

cat(sprintf("Selection Bias = %.03f - %.03f = %.03f \n",ey01,ey00,selection_bias))

# Calculate the share of units treated with vents (pi)
pi = mean(df$vents)

# Manually calculate the simple difference in mean health outcomes between the vent and non-vent group
ey1 = df[vents == TRUE, mean(y)]
ey0 = df[vents == FALSE, mean(y)]
sdo = ey1 - ey0
cat(sprintf("Simple Difference-in-Outcomes = %.03f - %.03f = %.03f \n",ey1,ey0,sdo))


# Calculate the simple difference in mean health outcomes between the vent and non-vent group using an OLS specification
reg = feols(y ~ vents,data = df,vcov = "hc1")
reg



# Fill out table with all this information

# Were you able to estimate the ATE, ATT or the ATU using the SDO?  Why/why not?
sdo_check = ate + selection_bias + (1 - pi) * (att - atu)

# --- Visualization: Simple density plot of observed outcomes -----------------
library(ggplot2)

burgundy <- "#800020"

# (Optional) subsample for speed if needed (uncomment next 2 lines)
# set.seed(20200403)
# df_plot <- df[sample(.N, 5000)] else:
df_plot <- df

# Group means (already computed above as ey1, ey0, but recompute to be safe)
mean_treated  <- df_plot[vents == TRUE,  mean(y)]
mean_control  <- df_plot[vents == FALSE, mean(y)]
SDO_value     <- mean_treated - mean_control

p_obs <- ggplot(df_plot, aes(x = y, fill = factor(vents), color = factor(vents))) +
  geom_density(alpha = 0.25, linewidth = 0.9) +
  scale_fill_manual(values = c("FALSE" = "gray70", "TRUE" = burgundy),
                    labels = c("FALSE" = "0 = No Vent (Control)",
                               "TRUE"  = "1 = Vent (Treated)")) +
  scale_color_manual(values = c("FALSE" = "gray50", "TRUE" = burgundy),
                     labels = c("FALSE" = "0 = No Vent (Control)",
                                "TRUE"  = "1 = Vent (Treated)")) +
  geom_vline(xintercept = mean_treated, color = burgundy, linetype = "dashed") +
  geom_vline(xintercept = mean_control, color = "gray40", linetype = "dashed") +
  labs(
    title = "Observed Outcomes by Treatment (Perfect Doctor Assignment)",
    subtitle = sprintf("Treated mean = %.2f | Control mean = %.2f | SDO = %.2f",
                       mean_treated, mean_control, SDO_value),
    x = "Observed Outcome (y)", y = "Density",
    fill = "Group", color = "Group"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", color = burgundy, hjust = 0.5),
    plot.subtitle = element_text(color = "gray30", hjust = 0.5),
    legend.position = "top"
  )

print(p_obs)


# Perfect Doctor calculations


# Causal parameters
ATE <- df[, mean(delta)]
ATT <- df[vents == TRUE, mean(delta)]
ATU <- df[vents == FALSE, mean(delta)]

# Selection bias terms
EY0_treated  <- df[vents == TRUE,  mean(y0)]
EY0_untreat  <- df[vents == FALSE, mean(y0)]
selection_bias <- EY0_treated - EY0_untreat

# Pi (share treated)
pi <- mean(df$vents)

# Simple differences
EY_treated   <- df[vents == TRUE,  mean(y)]
EY_untreat   <- df[vents == FALSE, mean(y)]
SDO_manual   <- EY_treated - EY_untreat

# OLS SDO check
reg <- feols(y ~ vents, data = df, vcov = "hc1")
SDO_OLS <- coef(reg)["ventsTRUE"]

# SDO decomposition check
SDO_decomp <- ATE + selection_bias + (1 - pi) * (ATT - ATU)

Obs <- nrow(df)

# Print
cat("\n=============================\nPerfect Doctor Summary\n=============================\n")
cat(sprintf("ATE                 = %.3f\n", ATE))
cat(sprintf("ATT                 = %.3f\n", ATT))
cat(sprintf("ATU                 = %.3f\n\n", ATU))
cat(sprintf("E[Y0|D=1]           = %.3f\n", EY0_treated))
cat(sprintf("E[Y0|D=0]           = %.3f\n", EY0_untreat))
cat(sprintf("Selection Bias      = %.3f\n\n", selection_bias))
cat(sprintf("Pi (share treated)  = %.3f\n", pi))
cat(sprintf("SDO (manual)        = %.3f\n", SDO_manual))
cat(sprintf("SDO (OLS)           = %.3f\n", SDO_OLS))
cat(sprintf("SDO (decomposition) = %.3f\n", SDO_decomp))
cat(sprintf("Observations        = %d\n", Obs))



