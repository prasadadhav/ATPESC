# This is the file to post process the Results from XDEM+OF Direct coupled simulations
#!/usr/bin/Rscript

# Have the xdem reader directory along with the R script file location
# devtools::load_all("xdemreader-0.0.1.0000/")
devtools::load_all("xdemreader/")
library(xdemreader)

library(tidyverse)
library(ggplot2)
library(plotly)
library(esquisse)

library(extrafont)
font_import()
loadfonts(device="win")   #Register fonts for Windows bitmap output

fonts()                       #vector of font family names
##  [1] "Andale Mono"                  "AppleMyungjo"                
##  [3] "Arial Black"                  "Arial"                       
##  [5] "Arial Narrow"                 "Arial Rounded MT Bold" 

library(viridis)
# ----------------------------------------------------------------------------------------------
#                 Load experimental data
# ----------------------------------------------------------------------------------------------

data_filename_1 <- "pi_HPC_100000000.csv"
df1 <- read_csv(data_filename_1, col_names = TRUE)

df1 <- df1 %>%
  mutate(`PI value` = as.double(`PI value`)) %>%
  mutate(`Time [s]` = as.double(`Time [s]`)) %>%
  mutate(`OMP Threads` = as.integer(`OMP Threads`)) %>%
  mutate(`Iteration No.` = as.factor(`Iteration No.`)) 



# # Group by 'OMP Threads' and calculate the average time for each group
# df_avg <- df1 %>%
#   group_by(`OMP Threads`) %>%
#   summarise(Avg_Time = mean(`Time [s]`))

# Calculate the average time for each combination of OMP Threads, Iteration No., and Label
avg_time_df <- aggregate(`Time [s]` ~ `OMP Threads`  + `Iteration No.` + Label, data = df1, FUN = mean)


# Create the plot using ggplot2
ggplot(avg_time_df, aes(x = `OMP Threads`, y = `Time [s]`, color = Label)) +
  # geom_line() +
  geom_point() +
  labs(title = "Time vs OMP Threads",
       x = "OMP Threads",
       y = "Time [s]",
       color = "Label") +
  scale_color_viridis(discrete = TRUE) +  # Use Viridis color scale
  theme_minimal()


ggplot(df1, aes(x = `OMP Threads`, y = `Time [s]`, color = Label)) +
  geom_point() +
  # geom_smooth(formula = y ~ poly(x, 4)) +
  labs(title = "Time vs OMP Threads",
       x = "OMP Threads",
       y = "Time [s]",
       color = "Label") +
  scale_color_viridis(discrete = TRUE) +  # Use Viridis color scale
  theme_linedraw() +
  theme(legend.position = c(0.75, 0.8), legend.title= element_blank(), legend.key.width = unit(2.5,"line"), text = element_text(family="Times New Roman", size = 20), legend.text = element_text(size = 12)) 
ggsave('./2024_03_03_pi_HPC_100000000.png', width=12, height=6)























