# Clear workspace
rm(list = ls())

# Load necessary libraries
library(ggplot2)
library(sf)
library(dplyr)
library(maps)
library(ggthemes)
library(readr)
library(kableExtra)

# Source external functions
source("functions.R")

# Load US map data and remove Alaska
states <- map_data("state") %>%
  filter(region != "alaska")

# Set seed for reproducibility
set.seed(1)

# Load and preprocess data
df <- read_csv("57.to.2008.raw.data.Aug4.2014.csv") %>%
  filter(year > 1957, state != 2) %>%
  arrange(year, state)

# Select variable based on index
ind <- 5
Y <- switch(ind,
            `1` = df$`inc tax pi`,
            `2` = df$nonresrevpi,
            `3` = df$totexppi,
            `4` = df$`edu exp pi`,
            `5` = df$`savings pi`)

X <- df$resrevpi

# Define constants
N <- 49
T <- 51

# Create ID and time vectors
id <- numeric(N * T)
for (i in 1:N) {
  for (t in 1:T) {
    id[(t - 1) * N + i] <- i
  }
}

time <- numeric(N*T)
for (i in (1:N)){
  for (t in (1:T)){
    time[((t-1)*N+1):(t*N)]<-t
  }
}

# Standardize Y and X
data <- data.frame(
  Y = c(as.matrix(Y)),
  X = c(as.matrix(X)),
  id = id,
  time = time
)

data <- data %>%
  mutate(
    Y = (Y - mean(Y)) / sd(Y),
    X = (X - mean(X)) / sd(X)
  )

# Define regression formula
formula <- Y ~ X - 1
formula_vars <- all.vars(formula)
dependent_var <- formula_vars[1]
independent_vars <- formula_vars[-1]

# Reshape data
Y_matrix <- reshape_to_matrix(data, "id", "time", dependent_var)
X_list <- lapply(independent_vars, function(var) reshape_to_matrix(data, "id", "time", var))
X_combined <- do.call(cbind, X_list)

# Perform clustering
set.seed(1)
cluster_result <- cluster_general(Y_matrix, X_list, N, T, type = "long")
G <- cluster_result$clusters
klong <- cluster_result$res

# Create cluster data and merge with map data
cluster_data <- data.frame(
  region = unique(states$region),
  cluster = klong$cluster
)

map_data_clustered <- states %>%
  left_join(cluster_data, by = "region")

# Plot the map
p <- ggplot(map_data_clustered, aes(long, lat, group = group, fill = as.factor(cluster))) +
  geom_polygon(color = "white") +
  coord_fixed(1.3) +
  theme_void() +
  scale_fill_brewer(palette = "Set3", name = "Cluster") +
  theme(plot.margin = margin(20, 20, 20, 20))

# Save the plot
ggsave("map.pdf", plot = p, width = 10, height = 6)

# Perform tall clustering
clustert <- cluster_general(Y_matrix, X_list, N, T, type = "tall")
C <- clustert$clusters
ktall <- clustert$res

# Prepare cluster centers for LaTeX output
dclusi <- data.frame(t(klong$centers)) %>%
  mutate(across(everything(), ~ round(., 3)))
dclusi %>%
  kable("latex", align = "c", linesep = "", escape = FALSE)

# Prepare cluster assignments for LaTeX output
clust <- data.frame(Year = 1958:2008, Cluster = ktall$cluster)
clust %>%
  kable("latex", align = "c", linesep = "", escape = FALSE, row.names = FALSE)

# Prepare cluster centers for LaTeX output
dclust <- data.frame(Cluster = 1:9, Center = ktall$centers) %>%
  mutate(across(everything(), ~ round(., 3)))
dclust %>%
  kable("latex", align = "c", linesep = "", escape = FALSE, row.names = FALSE)
