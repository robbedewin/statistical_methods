# Flu Shot Analysis - Logistic Regression
# Biostatistics Course Project
# Group 7

# Load required packages --------------------------------------------------
library(MASS)    # for stepAIC
library(ggplot2) # for visualization

# Set working directory 
setwd("C:/Users/robbe/Documents/KULeuven/Bioinformatics_2024-2025/Statistical_Methods/Assignment 2025")


# 1. Data Preparation -----------------------------------------------------
# Load data from assignment_data.txt (TSV format)
flu_data <- read.table("assignment_data.txt", header = TRUE, sep = "\t")

# Convert gender to factor
flu_data$gender <- factor(flu_data$gender, levels = c(0, 1), labels = c("Female", "Male"))

# Convert Y (flu shot status) to a labeled factor
flu_data$flu_shot_status <- factor(flu_data$Y, 
                                   levels = c(0, 1), 
                                   labels = c("No Flu Shot", "Flu Shot"))

# Check the first few rows of the data
head(flu_data)

# 2. Exploratory Analysis -------------------------------------------------
cat("\n=== Descriptive Statistics ===\n")
summary(flu_data)

cat("\n=== Flu Shot Distribution ===\n")
table(flu_data$flu_shot_status)

cat("\n=== Gender Distribution ===\n")
table(flu_data$gender)


# 3. Plot: Age Distribution by Flu Shot Status ----------------------------
library(ggplot2)
ggplot(flu_data, aes(x = age, fill = flu_shot_status)) +
  geom_histogram(position = "dodge", binwidth = 1, color = "black") +
  labs(
    title = "Age Distribution by Flu Shot Status",
    x = "Age",
    y = "Count",
    fill = "Flu Shot Status"
  ) +
  scale_y_continuous(breaks = scales::breaks_extended(n = 10)) + # Force integer breaks
  theme_minimal()

# 4. Plot: Health Awareness by Flu Shot Status ----------------------------
ggplot(flu_data, aes(x = flu_shot_status, y = HA, fill = flu_shot_status)) +
  geom_boxplot() +
  labs(
    title = "Health Awareness by Flu Shot Status",
    x = "Flu Shot Status",
    y = "Health Awareness Index (HA)",
    fill = "Flu Shot Status"
  ) +
  scale_fill_manual(values = c("#F8766D", "#00BFC4")) +
  theme_minimal() +
  theme(legend.position = "none")


ggplot(flu_data, aes(x = flu_shot_status, y = HA, fill = flu_shot_status)) +
  geom_boxplot() +
  scale_fill_manual(values = c("Not Taken" = "orange", "Taken" = "green")) +
  labs(
    title = "Health Awareness by Flu Shot Status",
    x = "Flu Shot",
    y = "Health Awareness Index"
  ) +
  scale_x_discrete(labels = c("No Flu Shot", "Flu Shot")) +
  theme_minimal() +
  theme(legend.position = "none")  # Remove color legend

# Health Awareness Exploration Plot
ggplot(flu_data, aes(x = age, y = HA)) +
  geom_point(aes(color = flu_shot_status, shape = gender), 
             size = 3, alpha = 0.8, position = position_jitter(width = 0.3)) +
  geom_smooth(aes(color = flu_shot_status), 
              method = "loess", se = FALSE, linewidth = 1.2) +
  facet_wrap(~gender) +
  labs(
    title = "Health Awareness Patterns",
    subtitle = "Relationship between age, HA score, gender, and flu shot status",
    x = "Age",
    y = "Health Awareness Index (HA)",
    color = "Flu Shot Status",
    shape = "Gender"
  ) +
  scale_color_manual(values = c("Taken" = "#2c7bb6", "Not Taken" = "#d7191c")) +
  scale_shape_manual(values = c("Female" = 16, "Male" = 17)) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.spacing = unit(1.5, "lines"),
    strip.text = element_text(face = "bold", size = 12)
  )



# 3. Logistic Regression Analysis ------------------------------------------
# Initial model
full_model <- glm(Y ~ age + HA + gender, 
                  family = binomial(link = "logit"), 
                  data = flu_data)

cat("\n=== Initial Model Summary ===\n")
summary(full_model)

# Model selection using AIC
step_model <- stepAIC(full_model, direction = "both", trace = 0)

cat("\n=== Final Model after Stepwise Selection ===\n")
summary(step_model)

# 4. Model Diagnostics ----------------------------------------------------
# Odds Ratios and CI
cat("\n=== Odds Ratios with 95% CI ===\n")
exp(cbind(OR = coef(step_model), confint(step_model)))

# Model comparison
cat("\n=== Model Comparison ===\n")
anova(step_model, test = "Chisq")

# 5. Visualization --------------------------------------------------------
# Coefficient plot
coef_plot <- data.frame(
  Predictor = names(coef(step_model)),
  OR = exp(coef(step_model))
)

ggplot(coef_plot[-1,], aes(x = Predictor, y = OR)) +
  geom_point() +
  geom_errorbar(aes(ymin = exp(confint(step_model)[-1,1]),
                    ymax = exp(confint(step_model)[-1,2])),
                width = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(title = "Odds Ratios with Confidence Intervals") +
  coord_flip()

# 6. Prediction Analysis --------------------------------------------------
# Create prediction plot for main effects
new_data <- expand.grid(
  age = seq(min(flu_data$age), max(flu_data$age), length = 50),
  HA = mean(flu_data$HA),
  gender = levels(flu_data$gender)
)

new_data$pred <- predict(step_model, newdata = new_data, type = "response")

ggplot(new_data, aes(x = age, y = pred, color = gender)) +
  geom_line() +
  labs(title = "Predicted Probability of Flu Shot by Age and Gender",
       y = "Predicted Probability") +
  theme_minimal()

