# -----------------------------
# FONDUE ALCOHOL EVAPORATION MODEL
# WITH UNIFORM CONSUMPTION
# -----------------------------

# Load libraries
library(ggplot2)

# -----------------------------
# USER PARAMETERS (modifiable)
# -----------------------------

# --- Cooking (pre-serving) ---
cook_time <- 7   # minutes of simmering with boiling before serving (typical)
boiling_evaporation_constant <- -log(0.40) / 15  
# According to a study by the USDA, about 40% of the alcohol remains after 15 minutes of simmering while boiling

# --- Under-burner (serving/eating) ---
eat_time <- 20   # time while eating (minutes)
burner_evaporation_constant <- -log(0.60) / 40  
# Assuming 40% loss (i.e., 60% remaining) after 40 minutes under the burner
# NOTE: keep in mind the numeric here is the k that satisfies exp(-k * 40) = 0.60

# --- Alcohol content of ingredients ---
wine_vol <- 240    # ml of wine
wine_abv <- 0.12   # 12% ABV
kirsch_vol <- 30   # ml of kirsch
kirsch_abv <- 0.40 # 40% ABV

# --- Beer comparison ---
beer_vol <- 330    # ml, standard European beer volume
beer_abv <- 0.05   # 5% ABV (typical)
ethanol_density <- 0.789  # g/ml

# -----------------------------
# CALCULATIONS
# -----------------------------

# 1. Compute initial ethanol
ethanol_init_ml <- wine_vol * wine_abv + kirsch_vol * kirsch_abv
ethanol_init_g <- ethanol_init_ml * ethanol_density

# 2. Exponential decay during cooking and eating
k_cook <- boiling_evaporation_constant
k_burner <- burner_evaporation_constant

# After cooking (simmering)
fraction_after_cook <- exp(-k_cook * cook_time)
ethanol_after_cook_g <- ethanol_init_g * fraction_after_cook

# -----------------------------
# UNIFORM CONSUMPTION MODEL
# -----------------------------
# We assume uniform consumption throughout the eating period.
# The fondue is consumed at a constant rate of (1/eat_time) per minute.
# After eat_time minutes, 100% of the fondue has been consumed.
#
# At time s (where 0 ≤ s ≤ eat_time), the ethanol concentration is:
#   C(s) = ethanol_after_cook_g * exp(-k_burner * s)
#
# The total ethanol consumed is the integral:
#   ∫[0 to eat_time] (1/eat_time) * C(s) ds
#   = (1/eat_time) * ethanol_after_cook_g * ∫[0 to eat_time] exp(-k_burner * s) ds
#   = (1/eat_time) * ethanol_after_cook_g * [-1/k_burner * exp(-k_burner * s)]|[0 to eat_time]
#   = (ethanol_after_cook_g / (k_burner * eat_time)) * (1 - exp(-k_burner * eat_time))

# Calculate total ethanol consumed with uniform consumption
if (k_burner > 0) {
  ethanol_consumed_g <- (ethanol_after_cook_g / (k_burner * eat_time)) * 
                        (1 - exp(-k_burner * eat_time))
} else {
  # If no evaporation during eating (k_burner = 0), all ethanol is consumed
  ethanol_consumed_g <- ethanol_after_cook_g
}

# 3. Beer comparison
beer_ethanol_ml <- beer_vol * beer_abv
beer_ethanol_g <- beer_ethanol_ml * ethanol_density
beer_equiv <- ethanol_consumed_g / beer_ethanol_g

# -----------------------------
# MODEL CURVE FOR PLOT
# -----------------------------
time_cook <- seq(0, cook_time, by = 0.1)
time_eat <- seq(0, eat_time, by = 0.1)

# Stage 1 (cooking)
ethanol_cook <- ethanol_init_g * exp(-k_cook * time_cook)

# Stage 2 (eating) - showing the ethanol in the pot as it evaporates
ethanol_eat <- tail(ethanol_cook, 1) * exp(-k_burner * time_eat)

# Combine full timeline
time_total <- c(time_cook, cook_time + time_eat)
ethanol_total <- c(ethanol_cook, ethanol_eat)

# Assemble dataframe with additional, more relatable metrics
df <- data.frame(Time = time_total, Ethanol_g = ethanol_total)
df$Ethanol_ml <- df$Ethanol_g / ethanol_density
# ml of beer that contain the same amount of ethanol (i.e. beer-equivalent volume)
df$Beer_ml_equiv <- df$Ethanol_ml / beer_abv
df$Beers_equiv <- df$Beer_ml_equiv / beer_vol
df$Pct_remain <- df$Ethanol_g / ethanol_init_g * 100
df$Stage <- ifelse(df$Time <= cook_time, "Cooking", "Serving")

# -----------------------------
# PLOT
# -----------------------------
decay_eqs <- paste0(
  "Cooking: exp(-", round(k_cook,4), " * t)\n",
  "Serving: exp(-", round(k_burner,4), " * t)")

subtitle_text <- paste0("With uniform consumption during eating: equivalent to ", 
                        round(beer_equiv,2), " beers (approx)")

# Improved main plot: Beer-equivalent (ml) over time with a secondary axis showing beers
library(scales)
main_plot <- ggplot(df, aes(x = Time)) +
  geom_line(aes(y = Beer_ml_equiv, color = Stage), size = 1.4) +
  geom_area(aes(y = Beer_ml_equiv, fill = Stage), alpha = 0.12) +
  geom_vline(xintercept = cook_time, linetype = "dashed", color = "gray40") +
  geom_point(data = df[which.max(df$Time), ], aes(x = Time, y = Beer_ml_equiv), 
             color = "black", size = 2) +
  scale_color_manual(values = c("Cooking" = "#D73027", "Serving" = "#4575B4")) +
  scale_fill_manual(values = c("Cooking" = "#D73027", "Serving" = "#4575B4")) +
  scale_y_continuous(
    name = "Beer-equivalent (ml)",
    sec.axis = sec_axis(~ . / beer_vol, name = paste0("Beers (", beer_vol, " ml)"))
  ) +
  labs(title = "Beer-equivalent Alcohol Remaining in Cheese Fondue",
       subtitle = subtitle_text,
       x = "Time (minutes)") +
  annotate("text", x = cook_time, y = max(df$Beer_ml_equiv)*0.95,
           label = "Serving starts →", hjust = -0.1, size = 3.6, color = "gray20") +
  annotate("text", x = max(df$Time)*0.72, y = max(df$Beer_ml_equiv)*0.25,
           label = decay_eqs, hjust = 0, size = 3.3, color = "black") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top",
        legend.title = element_blank())

# Secondary plot: percent remaining (makes exponential decay shape clearer)
pct_plot <- ggplot(df, aes(x = Time, y = Pct_remain)) +
  geom_line(color = "darkgreen", size = 1.3) +
  geom_hline(yintercept = c(100, 50, 25), linetype = "dotted", color = "gray85") +
  labs(x = "Time (minutes)", y = "Percent of initial ethanol (%)",
       title = "Percent of Initial Ethanol Remaining in Pot") +
  scale_y_continuous(labels = function(x) paste0(round(x,0), "%")) +
  theme_minimal(base_size = 12)

# Print both plots (user can arrange/ggsave as desired)
print(main_plot)
#print(pct_plot) # Uncomment to view percent remaining plot

# -----------------------------
# PRINT SUMMARY
# -----------------------------
cat("Initial ethanol:", round(ethanol_init_g,2), "g\n")
cat("After cooking (", cook_time, "min):", round(fraction_after_cook*100,1), "% remains\n")
cat("Ethanol consumed (uniform consumption over", eat_time, "min):", 
    round(ethanol_consumed_g,2), "g\n")
cat("Beer equivalent (approx):", round(beer_equiv,2), " beers\n")


# -----------------------------
# RESULTS BELOW HERE
# -----------------------------

####################################################################################################################
### With UNIFORM CONSUMPTION throughout eating (100% consumed by end of eat_time):
### Half a serving per person (with default parameters) is about 0.78 beers equivalent.
###
### This corresponds to a blood alcohol concentration (BAC) of approximately 0.0125% for an average adult, 
### which is described as a slight "buzz", or a feeling of warmth and relaxation, but generally not impairing.
###
### To get "legally drunk" (0.08% BAC), one would require approximately 2.6 servings of Fondue.
####################################################################################################################

# -----------------------------
# COMMENTS on METHODS USED
# -----------------------------
# I used the Widmark approach to estimate how much pure ethanol (grams) would
# be needed to reach a target BAC (blood alcohol concentration).
#
# Standard Widmark form (percent BAC units):
#   BAC_percent = (A_g / (r * weight_kg)) * 0.1  -  beta * t_hours
# where:
#   - A_g is the mass of ethanol consumed (grams)
#   - r is the Widmark distribution factor (approx. 0.70 for average male,
#     0.60 for average female)
#   - weight_kg is body weight in kg
#   - beta is elimination rate in % BAC per hour (optional if using a
#     metabolized-grams correction)
#   - t_hours is time over which alcohol was consumed (hours)
#
# Note on the exponential-decay equation used in the model:
#   - I modeled alcohol loss as first-order (exponential): fraction(t) = exp(-k * t),
#     where k (1/min) is the rate constant and t is time in minutes.
#   - The half-life is t1/2 = ln(2) / k; this assumes the loss rate is proportional
#     to the amount remaining (a common simple mass-transfer approximation).
#
# UNIFORM CONSUMPTION MODEL:
#   - The model assumes fondue is consumed uniformly throughout the eating period,
#     with a constant consumption rate = 1/eat_time per minute.
#   - By the end of eat_time, 100% of the fondue has been consumed.
#   - The total ethanol consumed is calculated by integrating the product of:
#     (consumption rate) × (ethanol concentration at time t) over [0, eat_time].
#   - This gives: ∫[0,T] (1/T) * E₀ * exp(-k*t) dt = (E₀/(k*T)) * (1 - exp(-k*T))
#     where E₀ is ethanol after cooking, k is burner evaporation rate, T is eat_time.
#   - This integral represents the time-weighted average of the ethanol content.
#   - Since the fondue is consumed continuously while evaporation is ongoing, the
#     total ethanol intake accounts for both the higher concentration at the beginning
#     and the lower concentration toward the end of the eating period.