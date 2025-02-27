.libPaths("C:/Users/24976197/OneDrive - UTS/Desktop/R/Packages")
input.dir <- "C:/Users/24976197/OneDrive - UTS/Desktop/PhD UTS Data folder/004 Bioinformatics/Z. Other people"
output.dir <- "C:/Users/24976197/OneDrive - UTS/Desktop/PhD UTS Data folder/004 Bioinformatics/Z. Other people"

##############################################################################
################## Load used libraries within script #########################
##############################################################################
library(ggplot2)
library(readxl)

##############################################################################
#################### Load in raw data for analysis ###########################
##############################################################################
data <- read_excel(
    path = file.path(input.dir, "Samara Bridge - E7G2 naresh sample.xlsx"),
    sheet = 'tethaQuick', # Change to sheet name in Excel
    skip = 40, # Amount of rows to skip until big table
    col_names = TRUE
) # If there is header present

##############################################################################
#################### Determine period of measurement #########################
##############################################################################
start.point <- 392 # Time:s value
end.point <- 4399 # Time:s value

# Helps further down to select the period you are interested in
select.period <- which(data[, "Time:s"] == start.point):which(
    data[, "Time:s"] == end.point
)

# To select correct samples for the correct variable
conductance.samples <- grep(pattern = "Gm", x = colnames(data))
capacitance.samples <- grep(pattern = "Cm", x = colnames(data))

# Assign the time points and data points to their own variable
time.points <- data[select.period, "Time:s"]

################################################################################
# Conductance samples #
data.points <- as.matrix(data[select.period, conductance.samples])

# Normalizing datapoints to the starting point
data.points <- as.data.frame(sweep(
    x = data.points,
    MARGIN = 2,
    STATS = data.points[1, ],
    FUN = "/"
))

# Combining the time and normalized data
conductance <- cbind(time.points, data.points)
colnames(conductance) <- gsub(
    pattern = ":.*",
    replacement = "",
    colnames(conductance)
)

# Setting the start point as Time == 0
conductance$Time <- conductance$Time - conductance$Time[1]

# Transform table for graphing purposes
conductance <- reshape2::melt(conductance, id = "Time")

################################################################################
# Capacitance samples #
data.points <- as.matrix(data[select.period, capacitance.samples])

# Normalizing datapoints to the starting point
data.points <- as.data.frame(sweep(
    x = data.points,
    MARGIN = 2,
    STATS = data.points[1, ],
    FUN = "/"
))

# Combining the time and normalized data
capacitance <- cbind(time.points, data.points)
colnames(capacitance) <- gsub(
    pattern = ":.*",
    replacement = "",
    colnames(capacitance)
)

# Setting the start point as Time == 0
capacitance$Time <- capacitance$Time - capacitance$Time[1]

# Transform table for graphing purposes
capacitance <- reshape2::melt(capacitance, id = "Time")


##############################################################################
#################### Creating the conductance/capacitance figure #############
##############################################################################
ggplot(conductance, aes(x = Time, y = value)) +

    # Graphing the individual samples
    stat_summary(
        aes(colour = variable),
        fun.y = "mean",
        geom = "point",
        alpha = 0.2
    ) +
    stat_summary(
        aes(colour = variable),
        fun.data = mean_cl_normal,
        geom = "line",
        alpha = 0.2
    ) +
    # Colours of the individual samples
    scale_color_manual(
        values = c("grey10", "grey26", "grey42", "grey58", "grey74", "grey90")
    ) +

    # Graphing the mean, standard error of the mean and the line through the mean
    stat_summary(
        fun.y = "mean",
        geom = "line",
        size = 1,
        position = position_dodge(0.2),
        col = "firebrick1"
    ) +
    stat_summary(
        fun.data = "mean_se",
        geom = "errorbar",
        position = position_dodge(0.2),
        col = "firebrick1"
    ) +
    stat_summary(
        fun.y = "mean",
        geom = "point",
        size = 1,
        position = position_dodge(0.2),
        col = "firebrick1"
    ) +

    # Graphing the mean and the 95% confidence interval around the line
    stat_summary(geom = "line", fun.y = mean, size = 1, col = "deepskyblue") +
    stat_summary(
        geom = "ribbon",
        fun.data = mean_cl_normal,
        alpha = 0.1,
        col = "deepskyblue",
        fill = "deepskyblue"
    ) +

    # Labels
    labs(x = "Time (s)", y = "Normalized conductance [μS]/ capacitance [nF]") + # Change the axis titles
    guides(color = guide_legend(title = "Sample")) + # Change the title of the colour legend

    theme_bw()

# Saving the figure
ggsave(
    filename = "conductance.tiff",
    path = file.path(output.dir),
    dpi = 300,
    width = 18,
    height = 18,
    units = "cm"
)
