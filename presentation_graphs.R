library(tidyverse)
library(reldist)

# Load the saved census data and perform statistical analysis

load(file = "ct_demographics.Rda")

# Create a tible without Charles Parker's household, to show how much his wealth accumulation affected the Gini coefficient for the 50-59 age cohort in 1860
ct_1860_hh_excl <- ct_1860_hh %>%
  dplyr::filter(family != 33643101)

age_cats <- c("20 - 29", "30 - 39", "40 - 49", "50 - 59", "60 and over")

# Define a function to calculate Gini coefficient by age category for a given set of household data
calc_gini_by_age <- function(data, town, wealth_col) {
  result <- c()
  
  for (age in age_cats) {
    gini_val <- data %>%
      dplyr::filter(town == town, AGE_CAT == age) %>%
      summarise(
        n = n(),
        gini = ifelse(n >= 5, reldist::gini({{ wealth_col }}), NA)
      ) %>%
      pull(gini)
    
    result <- c(result, gini_val)
  }
  
  return(tibble(gini = result))
}

# Calculate all four series
gini_1850_realprop <- calc_gini_by_age(ct_1850_hh, 'Meriden', FAMILY_REALPROP)
gini_1860_realprop <- calc_gini_by_age(ct_1860_hh, 'Meriden', FAMILY_REALPROP)
gini_1860_wealth <- calc_gini_by_age(ct_1860_hh, 'Meriden',FAMILY_WEALTH)
gini_1860_wealth_excl <- calc_gini_by_age(ct_1860_hh_excl, 'Meriden', FAMILY_WEALTH)

# Combine
gini_all <- tibble(
  age_cat = age_cats,
  realprop_1850 = gini_1850_realprop$gini,
  realprop_1860 = gini_1860_realprop$gini,
  wealth_1860 = gini_1860_wealth$gini,
  wealth_1860_excl = gini_1860_wealth_excl$gini
)

write.csv(gini_all, "gini.csv", row.names = FALSE)

w1 <- ggplot(gini_all, aes(x = age_cat)) +
  geom_line(aes(y = realprop_1850, color = "Real Property 1850"), group = 1) +
  geom_line(aes(y = realprop_1860, color = "Real Property 1860"), group = 1) +
  geom_line(aes(y = wealth_1860, color = "Total Wealth 1860"), group = 1) +
  ylim(0.65, 1) +
  labs(title = "Meriden Wealth Inequality",
       x = "Age Cohort",
       y = "Gini Coefficient") +
  scale_color_manual(values = c("blue", "red", "green")) +
  theme_minimal() +
  theme(legend.title = element_blank(), legend.position = "bottom" ) +  # Larger base
  theme(
    plot.title = element_text(hjust = 0.5, size = 24, face = "bold", margin = margin(b = 10)),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 12),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5)
  )
ggsave("meriden_wealth_inequality.png", plot=w1, width = 10, height = 6)

w2 <- ggplot(gini_all, aes(x = age_cat)) +
  geom_line(aes(y = wealth_1860, color = "Total Wealth 1860"), group = 1) +
  ylim(0.65, 1) +
  labs(title = "Meriden Wealth Inequality",
       x = "Age Cohort",
       y = "Gini Coefficient") +
  scale_color_manual(values = c("green")) +
  theme_minimal() +
  theme(legend.title = element_blank(), legend.position = "bottom") +  # Larger base
  theme(
    plot.title = element_text(hjust = 0.5, size = 24, face = "bold", margin = margin(b = 10)),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 12),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5)
  )
ggsave("1860_meriden_wealth_inequality.png", plot=w2, width = 10, height = 6)

w3<- ggplot(gini_all, aes(x = age_cat)) +
  geom_line(aes(y = wealth_1860_excl, color = "Total Wealth 1860, excluding Charles Parker"), group = 1) +
  ylim(0.65, 1) +
  labs(title = "Meriden Wealth Inequality",
       x = "Age Cohort",
       y = "Gini Coefficient") +
  scale_color_manual(values = c("orange")) +
  theme_minimal() +
  theme(legend.title = element_blank(), legend.position = "bottom") +  # Larger base
  theme(
    plot.title = element_text(hjust = 0.5, size = 24, face = "bold", margin = margin(b = 10)),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 12),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5)
  )
ggsave("1860_meriden_wealth_inequality_excluding.png", plot=w3, width = 10, height = 6)

# Meriden has a disproportionate number of native-born, young-adult males.
# Stonington and New London, other Know Nothing hotbeds, show the same pattern,
# as do the cities of Hartford and New Haven.

# Define a function to identify maxima
find_peaks <- function(x, m = 10) {
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- sapply(which(shape < 0), FUN = function(i) {
    z <- i - m + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + m + 1
    w <- ifelse(w < length(x), w, length(x))
    if (all(x[c(z:i, (i + 2):w)] <= x[i + 1])) {
      return(i + 1)
    } else {
      return(numeric(0))
    }
  })
  pks <- unlist(pks)
  pks
}

ct_native_male_1850 <- ct_1850 %>%
  dplyr::filter(BIRTH == "native" & SEX == 1)
ct_native_male_1860 <- ct_1860 %>%
  dplyr::filter(BIRTH == "native" & SEX == 1)
ct_age_distribution <- ct_native_male_1850 %>%
  select(AGE) %>%
  arrange(AGE) %>%
  group_by(AGE) %>%
  summarize(n = n())

# Create the smoothed spline data as a data frame for ggplot
ct_native_male.ss <- smooth.spline(ct_age_distribution)
spline_data <- data.frame(
  AGE = ct_native_male.ss$x,
  n = ct_native_male.ss$y
)

# Plot with ggplot2
p1 <- ggplot(ct_age_distribution, aes(x = AGE, y = n)) +
  geom_point(color = "black", size = 2) +
  geom_line(data = spline_data, aes(x = AGE, y = n), color = "blue", linewidth = 1) +
  labs(
    title = "Statewide",
    x = "Age",
    y = "Count"
  ) +
  theme_minimal(base_size = 16) +  # Larger base
  theme(
    plot.title = element_text(hjust = 0.5, size = 24, face = "bold", margin = margin(b = 10)),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 12),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5)
  )
ggsave("connecticut_age_distribution.png", plot=p1, width = 10, height = 6)

for (t in distinct(ct_1850 %>% select(town))[[1]]) {
  native_male_1850 <- ct_native_male_1850 %>%
    dplyr::filter(town == t) %>%
    select(AGE) %>%
    arrange(AGE) %>%
    group_by(AGE) %>%
    summarize(n = n())
  
  native_male_1850.ss <- smooth.spline(native_male_1850)
  
  native_male_1860 <- ct_native_male_1860 %>%
    dplyr::filter(town == t) %>%
    select(AGE) %>%
    arrange(AGE) %>%
    group_by(AGE) %>%
    summarize(n = n())
  
  native_male_1860.ss <- smooth.spline(native_male_1860)
  
  # Find the maxima in the smoothed age-distribution data
  peaks <- find_peaks(native_male_1850.ss$y, m = 5)
  
  if (length(peaks) < 3) {
    for (peak in peaks) {
      if (t == 'Meriden' & native_male_1850.ss$y[peak] - native_male_1850.ss$y[10] > 10) {
        max_age <- max(max(native_male_1850$n), max(native_male_1850.ss$y), max(native_male_1860.ss$y))
        
        # Prepare spline data frames
        spline_1850 <- data.frame(AGE = native_male_1850.ss$x, n = native_male_1850.ss$y)
        spline_1860 <- data.frame(AGE = native_male_1860.ss$x, n = native_male_1860.ss$y)
        label_data <- bind_rows(
          spline_1850 %>% filter(AGE == max(AGE)) %>% mutate(label = "1850 smoothed", color_group = "blue"),
          spline_1860 %>% filter(AGE == max(AGE)) %>% mutate(label = "1860 smoothed", color_group = "red")
        )
        
        # Create ggplot
        p <- ggplot() +
          geom_point(data = native_male_1850, aes(x = AGE, y = n, color = "1850 ages"), size = 2) +
          geom_line(data = spline_1850, aes(x = AGE, y = n, color = "1850 smoothed ages"), linewidth = 1) +
          geom_line(data = spline_1860, aes(x = AGE, y = n, color = "1860 smoothed ages"), linewidth = 1, linetype = "dashed") +
          
          coord_cartesian(ylim = c(0, max_age)) +
          labs(
            title = t,
            x = "Age",
            y = "Count",
            color = ""
          ) +
          scale_color_manual(values = c("1850 ages" = "black", 
                                        "1850 smoothed ages" = "blue", 
                                        "1860 smoothed ages" = "red")) +
          theme_minimal(base_size = 16) +  # Larger base
          theme(
            plot.title = element_text(hjust = 0.5, size = 24, face = "bold", margin = margin(b = 10)),
            axis.title = element_text(size = 16),
            axis.text = element_text(size = 12),
            plot.margin = margin(5.5, 5.5, 5.5, 5.5)
          )
        
        ggsave("meriden_age_distribution.png", plot = p, width = 10, height = 6)
        break
      }
    }
  }
}