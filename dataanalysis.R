require(tidyverse)
library(tools) # For file_path_sans_ext

file_names <- list.files(pattern = "*.csv")
df_list <- map(file_names, ~ read_csv(.x))
max_length <- max(sapply(df_list, nrow))
first_col <- seq(0, (max_length - 1) * 0.1, by = 0.1)

padded_columns <- map(df_list, ~ {
  column <- .x[[2]] # Extract the second column
  length(column) <- max_length # Pad with NA to max_length
  column
})

# Create a new dataframe by aligning the second columns
data_combined <- bind_cols(padded_columns) %>%
  set_names(file_path_sans_ext(file_names)) %>%
  mutate("time" = first_col) 


# Create plots for each combination of Mh and log age

all_plots = list()

plot_dir <- "plots"
dir.create(plot_dir, showWarnings = FALSE) # Create directory for plots

value_columns <- colnames(data_combined)[1:6]

# Create a list of ggplot objects for each column
plot_list <- map(value_columns, ~ {
  # Create the plot
  p <- ggplot(data_combined %>% filter(time >= 20 & time <= 40), aes(time, !!sym(.x))) +
    geom_line() +
    labs(
      title = paste("Time-Angular Velocity Plot for", .x),
      x = "Time (s)",
      y = "Angular Velocity (rad/s)"
    )
  
  # Save the plot to the 'plots' directory
  file_name <- file.path("plots", paste0(.x, ".png"))
  ggsave(file_name, plot = p, width = 6, height = 4, dpi = 300)
  
  # Return the plot
  p
})


multiplot <- gridExtra::arrangeGrob(grobs=plot_list)
ggsave(filename = file.path(plot_dir, "gridplot.png"), plot = multiplot, width = 15, height = 20)

L1 = 0.189
L2 = 0.325
M = 0.0301
g = 9.81
I = 0.01359375
calculate_frequency <- function(signal) {
  # Perform a Fourier transform to find the frequency
  signal <- na.omit(signal)
  fft_result <- fft(signal)
  
  # Compute the frequency axis
  n <- length(signal)
  sample_rate <- 10  # Assuming evenly spaced data
  freq_axis <- seq(0, sample_rate / 2, length.out = n / 2)
  
  # Take the absolute value of the fft and identify the peak frequency
  magnitude <- Mod(fft_result)
  peak_freq <- freq_axis[which.max(magnitude[1:(n/2)])]  # Ignore the negative frequencies
  return(peak_freq)
}

final_data <- data_combined %>%
  filter(time >= 20 & time <= 40) %>%
  pivot_longer(!time, names_to = "run", values_to = "ang_vel") %>%
  group_by(run) %>%
  summarise(
    ang_vel_mean = abs(mean(ang_vel, na.rm=TRUE)),
    ang_vel_err = sd(ang_vel, na.rm=TRUE),
    nut_freq = calculate_frequency(ang_vel)
  ) %>%
  mutate(prec_period = (2 * pi) / ang_vel_mean, prec_period_err = prec_period * (ang_vel_err / ang_vel_mean)) %>%
  mutate(length = if_else(str_detect(run, "L1"), L1, L2)) %>%
  mutate(spin_period = c(0.16,0.50,0.42,0.45,0.18,0.38)) %>%
  mutate(calc_period = (4*pi^2*I)/(M*g*length*spin_period))

m = lm(calc_period ~ prec_period, final_data)
slope = m$coefficients[2]

ggplot(final_data) + 
  geom_point(aes(x=prec_period, y=calc_period)) +
  geom_line(aes(x=prec_period, y=m$fitted), col="red") +
  labs(title="Plot of Calculated and Measured Precession Period.", x="Measured Precession Period (s)", y="Calculated Precession Period (s)")
ggsave(file.path("plots", "CalcVsMeas.png"))
t.test(final_data$calc_period, final_data$prec_period, paired=TRUE)

