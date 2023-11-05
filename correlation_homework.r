
library(ppcor)
library(corrplot)


# Question 3

file_path <- "/home/sam/Documents/my_files/technical_biology/technical_biology_2/scientific_methodology/homework/correlation/hypervent.csv"

data <- read.table(file_path, header = TRUE, sep = ";")


x <- c(data$Normal)
y <- c(data$Hypervent)

x_norm <- (x - min(x)) / (max(x) - min(x))
y_norm <- (y - min(y)) / (max(y) - min(y))


plot(x_norm, y_norm, xlab = "Normalized Normal", ylab = "Normalized Hypervent",
     main = "Scatter Plot of Normalized Data", col = "blue", pch = 19)


pearson_cor1 <- cor(x, y, method = "pearson")
spearman_cor1 <- cor(x, y, method = "spearman")

print("Q 3")
cat("Pearson Correlation Coefficient:", pearson_cor1, "\n")
cat("Spearman Correlation Coefficient:", spearman_cor1, "\n")


# The reason one correlation coefficient is greater than the other
# lies in the type of relationship between the data and the 
# characteristics of the dataset. The Pearson correlation coefficient 
# measures linear relationships, and when the relationship between the 
# times is approximately linear, the Pearson coefficient tends to be larger. 
# On the other hand, the Spearman correlation coefficient evaluates
# monotonic relationships, regardless of the specific form of the relationship, 
# and it is less sensitive to non-linear relationships and outliers. If the 
# relationship between the times is not strictly linear and/or if outliers are 
# present, the Spearman coefficient tends to be smaller than the Pearson 
# coefficient, highlighting the distinct applications of these correlation
# measures based on data characteristics.





# Question 4

path <- "/home/sam/Documents/my_files/technical_biology/technical_biology_2/scientific_methodology/homework/correlation/faith.csv"

faith <- read.table(path, header = TRUE, sep = ",")

erup <- c(faith$eruptions)
wait <- c(faith$waiting)

erup_norm <- (erup - min(erup)) / (max(erup) - min(erup))
wait_norm <- (wait - min(wait)) / (max(wait) - min(wait))



plot(erup_norm, wait_norm, xlab = "Normalized Eruptions",
     ylab = "Normalized Waiting", main = "Scatter Plot of Normalized Data",
     col = "blue", pch = 19)


pearson_cor2 <- cor.test(erup, wait, method = "pearson")
print("Q 4")
print(pearson_cor2)

# The data is binomial ("Daten-Inseln" können Korrelation vortäuschen)






# Question 5

ice_path <- "/home/sam/Documents/my_files/technical_biology/technical_biology_2/scientific_methodology/homework/correlation/ice.csv"

ice_data <- read.table(ice_path, header = TRUE, sep = ";")


time <- ice_data$T
num_ice <- ice_data$Ice


plot(time, num_ice, xlab = "Time", ylab = "Ice",
     main = "Time vs. Ice", col = "blue", pch = 19)


pearson_cor3 <- cor(time, num_ice, method = "pearson")
spearman_cor3 <- cor(time, num_ice, method = "spearman")
kendall_cor3 <- cor(time, num_ice, method = "kendall")

print("Q 5")
cat("Pearson Correlation Coefficient:", pearson_cor3, "\n")
cat("Spearman Correlation Coefficient:", spearman_cor3, "\n")
cat("Kendall Correlation Coefficient:", kendall_cor3, "\n")




# Question 6

ant_path <- "/home/sam/Documents/my_files/technical_biology/technical_biology_2/scientific_methodology/homework/correlation/antilope.csv"

antilope <- read.csv(ant_path, header = TRUE, sep = ";")

pearson_kitze <- cor(antilope$Kitze, antilope$Population)
pearson_niederschlag <- cor(antilope$Niederschlag, antilope$Population)

print("Q 6")

print("Pearson correlation between Kitze and Population:")
print(pearson_kitze)
print("Pearson correlation between Niederschlag and Population:")
print(pearson_niederschlag)


pearson_test <- cor.test(antilope$Population,
                         antilope$Kitze, method = "pearson")

print("Pearson correlation test for Population and Kitze:")
print(pearson_test)

pearson_test_niederschlag <- cor.test(antilope$Population,
                                      antilope$Niederschlag, method = "pearson")

print("Pearson correlation test for Population and Niederschlag:")
print(pearson_test_niederschlag)


cor_matrix <- cor(antilope[, c("Kitze", "Niederschlag", "Population")],
                  method = "spearman")

print("Correlation matrix:")
print(cor_matrix)



# Question 7

blood_path <- "/home/sam/Documents/my_files/technical_biology/technical_biology_2/scientific_methodology/homework/correlation/blood.csv"


blood <- read.csv(blood_path, header = TRUE, sep = ",")

model <- lm(sys ~ age + weight, data = blood)

print("Q 7")

summary(model)