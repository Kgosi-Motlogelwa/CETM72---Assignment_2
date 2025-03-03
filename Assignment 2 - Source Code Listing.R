# Full Data Pre-processing
# ==============================
# ==============================

#Data Loading
# ==============================
library(missForest)
library(Hmisc) # Imputation
library(VIM) # Visualization

wisconsin <- read.csv(
  "N:/My Drive/MISSION 30K/MSc Computer Science - Sunderland/CETM72 - Data Science Principles/wisconsin.csv",
  header = TRUE, stringsAsFactors = FALSE
)


# Visualizing Missing Data
# ==============================
png("visualize_missing_data.png", width = 800, height = 600)
aggr(wisconsin, col = c('navyblue', 'yellow'), 
     numbers = TRUE, sortVars = TRUE,
     labels = names(wisconsin), cex.axis = .7, 
     gap = 3, ylab = c("Missing data", "Pattern"))


# Missing Data Imputation
# ==============================
wisconsin$Bare.nuclei <- factor(wisconsin$Bare.nuclei, ordered = TRUE)
# Collapse level 6 into level 5
levels(wisconsin$Bare.nuclei)[levels(wisconsin$Bare.nuclei) == "6"] <- "5"

impute_arg <- aregImpute(
  Bare.nuclei ~ Cl.thickness + Cell.size + Cell.shape + 
    Marg.adhesion + Epith.c.size + Bl.cromatin + 
    Normal.nucleoli + Mitoses + Class, 
  data = wisconsin, 
  n.impute = 5,
  nk = 0 # Disables splines to avoid knot issues
)

# View imputed values
impute_arg$imputed$Bare.nuclei

# Convert 'Bare.nuclei' to an ordered factor
wisconsin$Bare.nuclei <- factor(wisconsin$Bare.nuclei, ordered = TRUE)

# Collapse level 6 into level 5
levels(wisconsin$Bare.nuclei)[levels(wisconsin$Bare.nuclei) == "6"] <- "5"

# Perform multiple imputations on 'Bare.nuclei'
impute_arg <- aregImpute(
  Bare.nuclei ~ Cl.thickness + Cell.size + Cell.shape + 
    Marg.adhesion + Epith.c.size + Bl.cromatin + 
    Normal.nucleoli + Mitoses + Class, 
  data = wisconsin, 
  n.impute = 5,
  nk = 0  # Disables splines to avoid knot issues
)

# Compute mode for each imputed value
final_imputed_values <- apply(impute_arg$imputed$Bare.nuclei, 1, function(x) {
  as.numeric(names(sort(table(x), decreasing = TRUE))[1])  # Mode (most frequent value)
})

# Assign imputed values
wisconsin$Bare.nuclei[is.na(wisconsin$Bare.nuclei)] <- final_imputed_values

# Convert 'Class' column to numeric (keeping 2 and 4)
wisconsin$Class <- ifelse(wisconsin$Class == "benign", 2, 4)


# Converting to Numeric for Transformation & PCA
# ===============================================
wisconsin_numeric <- wisconsin

wisconsin_numeric[] <- lapply(wisconsin_numeric, function(x) {
  if (is.factor(x) || is.character(x)) {
    return(as.numeric(as.character(x)))  # Convert factor/character to numeric
  } else {
    return(x)  # Keep numeric columns unchanged
  }
})


# DATA TRANSFORMATION
# ================================
# Load dataset
# save(wisconsin, file = "wisconsin.RData")
# load("wisconsin.RData")

# Create a directory to store the plots
dir.create("transformation_plots_log", showWarnings = FALSE)

# Loop through all numeric columns (excluding 'Class')
for (col in colnames(wisconsin_log)[-ncol(wisconsin_log)]) {  
  # Define the filename
  filename <- paste0("transformation_plots_log/", col, "_hist_qq.png")
  
  # Open a PNG device
  png(filename, width = 1200, height = 800)
  
  # Arrange plots in a 2x2 grid
  par(mfrow = c(1, 2))  
  
  # Histogram
  hist(wisconsin_log[[col]], 
       main = paste("Histogram of", col), 
       freq = FALSE, 
       col = "lightblue", 
       ylim=c(0, .4))
  
  # Q-Q Plot
  qqnorm(wisconsin_log[[col]], 
         main = paste("Q-Q Plot of", col))
  qqline(wisconsin_log[[col]], col = "red")  # Add normal reference line
  
  # Close the PNG device
  dev.off()
}

# Reset plot layout
par(mfrow = c(1, 1))


# Transformation using logmarithmic scale
wisconsin_log <- wisconsin_numeric

wisconsin_log[] <- lapply(wisconsin_log, function(x) {
  if (is.numeric(x)) {
    x[x == 0] <- NA  # Avoid log(0) errors
    return(log(x))
  } else {
    return(x)  # Keep non-numeric columns unchanged
  }
})

# Check transformed data
summary(wisconsin_log)

# Normalizing or scaling a variable. Z-score normalization
wisconsin_scaled <- as.data.frame(scale(wisconsin_numeric))
summary(wisconsin_scaled)

# DATA REDUCTION | PCA
# ================================
# Convert all non-numeric columns to numeric (except 'Class')
wisconsin_scaled_numeric <- wisconsin_scaled[, !names(wisconsin_scaled) %in% "Class"]

# Ensure all columns are numeric
wisconsin_scaled_numeric[] <- lapply(wisconsin_scaled_numeric, as.numeric)

# Perform PCA on numeric dataset (excluding 'Class')
pca_wisconsin <- prcomp(wisconsin_scaled_numeric[, !names(wisconsin_scaled_numeric) %in% "Class"], scale = TRUE)

# Summary of PCA
summary(pca_wisconsin)

pca_wisconsin$sdev^2

# Create the "PCA" folder if it doesn't exist
if (!dir.exists("PCA")) {
  dir.create("PCA")
}

# Save biplot as a PNG file
png("PCA/biplot_pca_wisconsin.png", width = 1200, height = 800)
biplot(pca_wisconsin, cex = 0.8)
abline(h = 0, v = 0, lty = 2, col = 8)
dev.off()  # Close the PNG device

# Save scree plot as a PNG file
png("PCA/screeplot_pca_wisconsin.png", width = 1200, height = 800)
screeplot(pca_wisconsin, type = "lines", col = 3)
dev.off()  # Close the PNG device

cor(wisconsin_numeric)
cor(wisconsin_log)

png("cor_wisconsin_log.png", width = 1200, height = 800)
plot(wisconsin_log)
png("cor_wisconsin_numeric.png", width = 1200, height = 800)
plot(wisconsin_numeric)

dev.off()  # Close the PNG device

# Assuming 'Class' is your target variable and wisconsin_numeric contains numeric columns

# Create an empty list to store the results
regression_results <- list()

# Loop through all numeric columns (excluding 'Class')
for (col in colnames(wisconsin_numeric)[-which(names(wisconsin_numeric) == "Class")]) {
  
  # Build the linear model for each variable
  formula <- as.formula(paste("Class ~", col))
  lmmodel <- lm(formula, data = wisconsin_numeric)
  
  # Store the summary in the results list
  regression_results[[col]] <- summary(lmmodel)
  
  # Print the summary for each model
  print(paste("Linear Regression between Class and", col))
  print(summary(lmmodel))
}

# Optionally, you can view the stored results after the loop
# To view the results of the first comparison (e.g., between Class and studytime):
print(regression_results$studytime)

# Create a directory to save plots
dir.create("scatterplots_with_regression_lines", showWarnings = FALSE)

# Loop through all numeric columns (excluding 'Class')
for (col in colnames(wisconsin_numeric)[-which(names(wisconsin_numeric) == "Class")]) {
  
  # Define the filename for saving the plot
  filename <- paste0("scatterplots_with_regression_lines/", col, "_scatterplot.png")
  
  # Open a PNG device to save the plot
  png(filename, width = 1200, height = 800)
  
  # Plot the scatter plot for 'Class' vs current column (col)
  plot(wisconsin_numeric[[col]], wisconsin_numeric$Class,
       main = paste("Scatterplot of Class vs", col),
       xlab = col, ylab = "Class", pch = 19)
  
  # Add the regression line (using `Class ~ .` to refer to the dependent variable)
  model <- lm(Class ~ wisconsin_numeric[[col]], data = wisconsin_numeric)
  abline(model, col = "red")
  
  # Close the PNG device
  dev.off()
}


# Simple Linear Regression
# ===========================
# Create a directory to store the plots
dir.create("regression_plots", showWarnings = FALSE)

# Loop through each column in the dataset (excluding 'Class')
for (col in colnames(wisconsin_numeric)) {
  if (col != "Class") {
    # Fit the linear model for each predictor against Class
    lmmodel <- lm(Class ~ get(col), data = wisconsin_numeric)
    
    # Create a file name for the plot
    filename <- paste0("regression_plots/plot_", col, ".png")
    
    # Open a PNG device to save the plot
    png(filename, width = 800, height = 600)
    
    # Plotting the regression for each predictor
    plot(wisconsin_numeric[[col]], as.numeric(wisconsin_numeric$Class),
         main=paste("Simple Linear Regression: Class ~", col),
         xlab=col, ylab="Class", pch=19)
    
    # Add the regression line
    abline(lmmodel, col="red")
    
    # Close the PNG device (saving the plot)
    dev.off()
    
    # Print model summary
    print(paste("Model for", col))
    print(summary(lmmodel))
  }
}


# Multiple Linear Regression
# ===========================
# Predict Class using all available variables
attach(wisconsin_numeric)

# Fit the Multiple Linear Regression model using all predictors
mrmodel <- lm(Class ~ Cl.thickness + Cell.size + Cell.shape + Marg.adhesion + 
                Epith.c.size + Bare.nuclei + Bl.cromatin + Normal.nucleoli + Mitoses, 
              data = wisconsin_numeric)

# View the summary of the multiple linear regression model
summary(mrmodel)

# Perform diagnostics to check multicollinearity (optional)
library(mctest)
imcdiag(mrmodel)

# Create a data frame for new predictions
newdata <- data.frame(Cl.thickness = 5, 
                      Cell.size = 2, 
                      Cell.shape = 3, 
                      Marg.adhesion = 1, 
                      Epith.c.size = 4, 
                      Bare.nuclei = 6, 
                      Bl.cromatin = 2, 
                      Normal.nucleoli = 1, 
                      Mitoses = 3)

# Predicting the outcome for new data with confidence intervals
predict(mrmodel, newdata, interval = "confidence")

# Plotting the relationship between Cl.thickness and Class (one predictor for simplicity)
plot(Cl.thickness, Class, main = "Scatterplot: Cl.thickness vs Class", 
     xlab = "Cl.thickness", ylab = "Class", pch = 19)

# Add regression line to the plot
abline(lm(Class ~ Cl.thickness), col = "red")

# Optional: Plot diagnostic plots for the multiple regression model
par(mfrow = c(2, 2))  # Set up a 2x2 grid for diagnostic plots
plot(mrmodel)

# Check diagnostics (multicollinearity) for the model
imcdiag(mrmodel)


# Diagnostic Plot
# =======================
# Fit the multiple linear regression model
fit_all <- lm(Class ~ Cl.thickness + Cell.size + Cell.shape + Marg.adhesion + 
                Epith.c.size + Bare.nuclei + Bl.cromatin + Normal.nucleoli + Mitoses, 
              data = wisconsin_numeric)

# View the summary of the model
summary(fit_all)

# Set up a 2x2 grid for diagnostic plots
par(mfrow = c(2, 2))

# Plot diagnostic plots for the regression model
plot(fit_all)

# Perform the diagnostic tests for multicollinearity
library(mctest)
imcdiag(fit_all)


# Preparing Data
# ===============

# Load necessary libraries
library(rpart)
library(rpart.plot)
library(RColorBrewer)
library(rattle)
library(caret)  # For confusionMatrix function
library(ggplot2)
library(graphics)
library(cluster)

# Load your dataset
# Ensure 'Class' is a factor in your dataset
wisconsin$Class <- as.factor(wisconsin$Class)

# Set seed for reproducibility
set.seed(42)

# Split the dataset into training and testing sets (70% train, 30% test)
trainIndex <- createDataPartition(wisconsin$Class, p = 0.7, list = FALSE)
trainData <- wisconsin[trainIndex, ]
testData <- wisconsin[-trainIndex, ]


# Decision Tree Visualization
# ============================

# Build the decision tree model using rpart
treeModel <- rpart(Class ~ ., data = trainData, method = "class",
                   control = rpart.control(minsplit = 1, minbucket = 1, cp = 0.001))

# Visualize the Decision Tree
# Basic plot
rpart.plot(treeModel, main = "Decision Tree for Wisconsin Dataset")

# Enhanced plot
fancyRpartPlot(treeModel, main = "Enhanced Decision Tree for Wisconsin Dataset")


# Model Evaluation: Confusion Matrix
# ===================================
# Ensure Class is a factor before splitting
wisconsin_numeric$Class <- as.factor(wisconsin_numeric$Class)

# Split the dataset into training (70%) and testing (30%) sets
set.seed(123)  # For reproducibility
index <- sample(nrow(wisconsin_numeric), 0.7 * nrow(wisconsin_numeric)) 
wisconsin_train <- wisconsin_numeric[index, ]   # Training set
wisconsin_test <- wisconsin_numeric[-index, ]   # Testing set

# Build the decision tree model
treeModel <- rpart(Class ~ ., data = wisconsin_train, method = "class",
                   control = rpart.control(minsplit = 1, minbucket = 1, cp = 0.001))

# Predict the classes for the test dataset
predictions <- predict(treeModel, newdata = wisconsin_test, type = "class")

# Ensure predictions and actual values have the same factor levels
predictions <- factor(predictions, levels = levels(wisconsin_test$Class))

# Evaluate the model using confusionMatrix from caret package
conf_mat <- confusionMatrix(predictions, wisconsin_test$Class)
print(conf_mat)

# Save the decision tree plot with a larger size
png("fancy_decision_tree_plot_large.png", width = 1600, height = 1200)
fancyRpartPlot(treeModel)
dev.off()

# Create confusion matrix as a table
confusion_matrix <- table(Predicted = predictions, Actual = wisconsin_test$Class)

# Convert the confusion matrix into a data frame for ggplot2
conf_matrix_df <- as.data.frame(confusion_matrix)
names(conf_matrix_df) <- c("Predicted", "Actual", "Count")

# Save the confusion matrix as an image using ggplot2
png("confusion_matrix_large.png", width = 1600, height = 1200)

ggplot(conf_matrix_df, aes(x = Actual, y = Predicted, fill = Count)) +
  geom_tile() +
  geom_text(aes(label = Count), color = "white", size = 6) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Confusion Matrix", x = "Actual", y = "Predicted") +
  theme(plot.title = element_text(hjust = 0.5))

dev.off()


# Classification Plots
# ====================
# Define the output folder
output_folder <- "classification_plots"
dir.create(output_folder, showWarnings = FALSE)  # Create folder if it doesn't exist

# Assign colors to classes (Benign = Blue, Malignant = Red)
class_colors <- c("blue", "red")[as.factor(wisconsin_numeric$Class)]

# Function to generate and save scatter plots of features vs Class
classification_plots_1 <- function(data, colors, folder) {
  feature_names <- colnames(data)[colnames(data) != "Class"]  # Exclude Class
  num_features <- length(feature_names)
  
  for (i in 1:num_features) {
    # Define the file name
    file_name <- paste0(folder, "/", feature_names[i], "_vs_Class.png")
    
    # Open a PNG device
    png(file=file_name, width=800, height=600)
    
    # Create the scatter plot of the feature vs Class
    plot(data[[feature_names[i]]], as.factor(data$Class), 
         pch=19, col=colors, 
         xlab=feature_names[i], ylab="Class", 
         main=paste(feature_names[i], "vs Class"))
    
    # Add a legend
    legend("topright", legend=c("Benign", "Malignant"), 
           col=c("blue", "red"), pch=19, bty="n")
    
    # Close the PNG device
    dev.off()
  }
}

# Call the function to generate and save plots
classification_plots_1(wisconsin_numeric, class_colors, output_folder)

# Print completion message
cat("Plots saved in:", output_folder)

# Assuming `wisconsin_scaled_numeric` is your dataset
d <- wisconsin_scaled_numeric
plot(d)

# K-means clustering (4 clusters)
# ===============================
set.seed(42)  # For reproducibility
cl <- kmeans(d, 5, nstart = 25)

# Plot the K-means clustering results
png("kmeans_plot.png", width=1600, height=1200)
plot(d, col = cl$cluster, pch = 19, main="K-means Clustering")
points(cl$centers, col = 1:5, pch = 15, cex = 2)  # Add cluster centers
dev.off()

# Silhouette plot (to evaluate the clustering quality)
# ====================================================
sk <- silhouette(cl$cluster, dist(d))

# Save silhouette plot
png("silhouette_plot.png", width=1600, height=1200)
plot(sk, main="Silhouette Plot for K-means Clustering")
dev.off()


# Hierarchical clustering
# ==========================
hca <- hclust(dist(d))  # Compute hierarchical clustering using the distance matrix

# Save hierarchical clustering plot
png("scaled_hca_plot.png", width=3000, height=2400)  # Further scale up the width and height
plot(hca, main="Hierarchical Clustering")
rect.hclust(hca, k=4, border="red")  # Highlight 4 clusters with red borders
dev.off()