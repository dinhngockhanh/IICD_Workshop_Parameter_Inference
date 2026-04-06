library(IICDWorkshopParameterInference)

setwd( #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    "/Users/kndinh/RESEARCH AND EVERYTHING/Projects/GITHUB/IICD_Workshop_Parameter_Inference/R"
) #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
files_sources <- list.files(pattern = "\\.[rR]$") #<<<<<<<<<<<<<<<<<<<<<
sapply(files_sources, source) #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
setwd( #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    "/Users/kndinh/RESEARCH AND EVERYTHING/Projects/GITHUB/IICD_Workshop_Parameter_Inference/vignettes"
) #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
library(ggplot2)
data <- read.csv( #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    "/Users/kndinh/RESEARCH AND EVERYTHING/Projects/GITHUB/IICD_Workshop_Parameter_Inference/data/data.csv"
) #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
initial_state <- list(
    rabbit = 700,
    fox = 300
)
parameters <- list(
    alpha = 2.0,
    beta = 4.0,
    gamma = 0.8,
    delta = 1.7
)
predicted_populations <- Lotka_Volterra_ODE_solver(
    parameters = parameters,
    times = seq(0, 14, by = 0.1),
    initial_state = initial_state
)
color_rabbit <- "#EFC000"
color_fox <- "#BC3C29"
obs_df <- rbind(
    data.frame(time = data$time, species = "Rabbit", count = data$rabbit),
    data.frame(time = data$time, species = "Fox", count = data$fox)
)
obs_df$species <- factor(obs_df$species, levels = c("Rabbit", "Fox"))
fit_df <- rbind(
    data.frame(time = predicted_populations$time, species = "Rabbit", count = predicted_populations$rabbit),
    data.frame(time = predicted_populations$time, species = "Fox", count = predicted_populations$fox)
)
fit_df$species <- factor(fit_df$species, levels = c("Rabbit", "Fox"))
p_fit <- ggplot() +
    geom_line(
        data = fit_df,
        aes(x = time, y = count, colour = species, group = species),
        linewidth = 1.2
    ) +
    geom_point(
        data = obs_df,
        aes(x = time, y = count, colour = species),
        size = 7,
        stroke = 0.3
    ) +
    scale_colour_manual(
        values = c(Rabbit = color_rabbit, Fox = color_fox),
        name = NULL
    ) +
    labs(x = "Time", y = "Count") +
    theme_bw(base_size = 20) +
    theme(
        plot.title = element_blank(),
        panel.grid = element_blank(),
        legend.position = c(0.98, 0.98),
        legend.justification = c(1, 1),
        legend.background = element_rect(
            fill = "white",
            colour = "grey35",
            linewidth = 0.5
        ),
        legend.margin = margin(6, 8, 6, 8)
    )
ggsave(
    filename = "Interactive_session_1.png",
    plot = p_fit,
    width = 10,
    height = 6,
    dpi = 300,
    bg = "white"
)
