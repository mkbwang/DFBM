
library(vegan)


rm(list=ls())
data <- readRDS("experiment/COHRA2/phyasv_visit24.rds")
taxonomy_table <- tax_table(data) |> data.frame()
species_names <- paste(taxonomy_table$Genus,
                       taxonomy_table$Species, sep=" ") |> as.character()
names(species_names) <- rownames(taxonomy_table)


metadata <- read.csv("experiment/COHRA2/metadata_yr2.csv", row.names=1)

diagnoses <-metadata$CaseEver == "Case"

raw_counts <- read.csv("experiment/COHRA2/16S_counts_yr2.csv", row.names=1)
raw_relabd <- t(apply(raw_counts, 1, function(row) {
  row_sum <- sum(row)
  row / row_sum
}))

F_statistics <- rep(0, 9)
for (j in 1:9){
  denoised_counts <- read.csv(sprintf("experiment/COHRA2/denoise/16S_yr2_%d.csv", j))
  denoised_relabd <- t(apply(denoised_counts, 1, function(row) {
    row_sum <- sum(row)
    row / row_sum
  }))
  denoised_bray <- vegdist(denoised_relabd, method="bray")
  denoised_permanova <- adonis2(formula=denoised_bray ~ diagnoses, permutations=999)
  F_statistics[j] <- denoised_permanova$F[1]
}
denoised_counts <- read.csv(sprintf("experiment/COHRA2/denoise/16S_yr2_%d.csv", 4))
denoised_relabd <- t(apply(denoised_counts, 1, function(row) {
  row_sum <- sum(row)
  row / row_sum
}))
denoised_bray <- vegdist(denoised_relabd, method="bray")
denoised_permanova <- adonis2(formula=denoised_bray ~ diagnoses, permutations=999)

taxa_names <- colnames(raw_counts)

library(ggplot2)

# MDS plot based on bray curtis distance
raw_bray <- vegdist(raw_relabd, method="bray")
raw_permanova <- adonis2(formula=raw_bray ~ diagnoses, permutations=999)
raw_coordinates <- cmdscale(raw_bray, k=2)
raw_df <- cbind(raw_coordinates, diagnoses) |> as.data.frame()
raw_df$diagnoses <- as.factor(raw_df$diagnoses)
raw_df$Type <- "Original"
ggplot(raw_df, aes(x=V1, y=V2, color=diagnoses)) +
  geom_point()
library(cluster)
library(compositions)
raw_counts[raw_counts == 0] <- 0.5
clr_raw_counts <- matrix(0, nrow=nrow(raw_counts), ncol=ncol(raw_counts))
for (j in 1:nrow(clr_raw_counts)){
  clr_raw_counts[j, ] <- clr(raw_counts[j, ])
}



denoised_coordinates <- cmdscale(denoised_bray, k=2)
denoised_silhouette <- silhouette(diagnoses, denoised_bray)
denoised_df <- cbind(denoised_coordinates, diagnoses) |> as.data.frame()
denoised_df$diagnoses <- as.factor(denoised_df$diagnoses)
denoised_df$Type <- "Denoised"


combined_coordinates <- rbind(raw_df, denoised_df)
combined_coordinates$Type <- factor(combined_coordinates$Type,
                                    levels=c("Original", "Denoised"))


cmap_vals <- c("#006622", "#ff4d4d")
ggplot(combined_coordinates, aes(x=V1, y=V2, color=diagnoses)) +
  geom_point(alpha=0.7) + labs(x="PcoA1", y="PcoA2") +
  scale_color_manual(values=cmap_vals)+
  facet_grid(cols=vars(Type))


clr_denoised_counts <- matrix(0, nrow=nrow(denoised_counts),
                              ncol=ncol(denoised_counts))
for (j in 1:nrow(clr_denoised_counts)){
  clr_denoised_counts[j, ] <- clr(denoised_counts[j, ])
}


diagnoses <- as.integer(diagnoses)


# logistic lasso
fit_logistic_lasso <- function(input, output, seed=1){

  set.seed(seed)
  index_positive <- which(output == TRUE)
  split_mask <- rbinom(n=length(index_positive), 1, prob=0.8)
  train_positive <- index_positive[split_mask == 1]
  test_positive <- index_positive[split_mask == 0]

  index_negative <- which(output == FALSE)
  split_mask <- rbinom(n=length(index_negative), 1, prob=0.8)
  train_negative <- index_negative[split_mask == 1]
  test_negative <- index_negative[split_mask == 0]

  train_labels <- c(output[train_positive], output[train_negative]) |> as.integer()
  test_labels <- c(output[test_positive], output[test_negative]) |> as.integer()

  train_data <- input[c(train_positive, train_negative), ]
  test_data <- input[c(test_positive, test_negative), ]

  cv_fit <- cv.glmnet(x=train_data, y=train_labels, family="binomial",
                      alpha=1, type.measure = "auc", nfolds=5)

  optimal_coefs <- coef(cv_fit, s = "lambda.min") |> as.vector()

  predicted_probabilities_train <-
    predict(cv_fit, s = "lambda.min", newx = train_data, type = "response") |> as.vector()

  roc_obj <- roc(train_labels, predicted_probabilities_train)
  train_auc_value <- auc(roc_obj) |> as.numeric()

  predicted_probabilities_test <-
    predict(cv_fit, s = "lambda.min", newx = test_data, type = "response") |> as.vector()

  roc_obj <- roc(test_labels, predicted_probabilities_test)
  test_auc_value <- auc(roc_obj) |> as.numeric()

  output <- list(coefs=optimal_coefs, train_auc=train_auc_value, test_auc=test_auc_value)
  return(output)

}

library(glmnet)
library(doParallel)
library(foreach)
numCores <- detectCores() - 1
cl <- makeCluster(numCores)

registerDoParallel(cl)

logitlasso_raw <- foreach(j=1:100, .packages=c("glmnet", "pROC")) %dopar% {
  result <- fit_logistic_lasso(input=clr_raw_counts,
                               output=diagnoses, seed=j)
  result
}


logitlasso_denoised <- foreach(j=1:100, .packages=c("glmnet", "pROC")) %dopar% {
  result <- fit_logistic_lasso(input=clr_denoised_counts,
                               output=diagnoses, seed=j)
  result
}



raw_aucs<- do.call(c, lapply(logitlasso_raw, function(result) result$test_auc))
raw_coefs <- do.call(rbind, lapply(logitlasso_raw, function(result) result$coefs))

num_selected_variables_raw <- rowSums(raw_coefs[, -1] != 0)
selected_frequency_raw <- colSums(raw_coefs[, -1] != 0)
top_features <- sort(selected_frequency_raw, index.return=TRUE, decreasing=TRUE)
top_feature_names <- taxa_names[top_features$ix[1:10]]
top_feature_frequencies <- top_features$x[1:10]
original_top_features <- data.frame(Name=top_feature_names,
                                    Frequency=top_feature_frequencies,
                                    Taxonomy = species_names[top_feature_names])

write.csv(original_top_features, "experiment/COHRA2/plots/original_top_features_yr2.csv",
          row.names=FALSE)



denoised_aucs <- do.call(c, lapply(logitlasso_denoised, function(result) result$test_auc))
denoised_coefs <- do.call(rbind, lapply(logitlasso_denoised, function(result) result$coefs))
num_selected_variables_denoised <- rowSums(denoised_coefs[, -1] != 0)
selected_frequency_denoised <- colSums(denoised_coefs[, -1] != 0)
top_features <- sort(selected_frequency_denoised, index.return=TRUE, decreasing=TRUE)
top_feature_names <- taxa_names[top_features$ix[1:10]]
top_feature_frequencies <- top_features$x[1:10]
denoised_top_features <- data.frame(Name=top_feature_names,
                                    Frequency=top_feature_frequencies,
                                    Taxonomy = species_names[top_feature_names])


write.csv(denoised_top_features, "experiment/COHRA2/plots/denoised_top_features_yr2.csv",
          row.names=FALSE)





auc_df <- data.frame(AUC=c(raw_aucs, denoised_aucs),
                    Type=rep(c("Original", "Denoised"), each=100))
auc_df$Type <- factor(auc_df$Type, levels=c("Original", "Denoised"))
ggplot(auc_df, aes(x=Type, y=AUC))+
  geom_boxplot() + labs(y="Test AUC", x=NULL)+
  scale_y_continuous(breaks = seq(0.5, 0.95, by = 0.1), limits=c(0.5, 0.95))



selected_vars <- data.frame(Count=c(num_selected_variables_raw, num_selected_variables_denoised),
                            Type=rep(c("Original", "Denoised"), each=100))

selected_vars$Type <- factor(selected_vars$Type, levels=c("Original", "Denoised"))
ggplot(selected_vars, aes(x=Type, y=Count))+
  geom_boxplot() + labs(y="Number of Selected Taxa", x=NULL)+
  scale_y_continuous(breaks = seq(0, 25, by = 5), limits=c(0, 25))


# selected_frequency_denoised <- colSums(denoised_coefs[, -1] != 0)


# try using random forest

library(pROC)
library(randomForest)
library(caret)


raw_auc <- rep(0, 25)
denoised_auc <- rep(0, 25)




for (j in 1:25){

  print(j)
  train_ids <- sample(nrow(denoised_relabd), nrow(denoised_relabd) * 0.8)
  test_ids <- setdiff(1:nrow(denoised_relabd), train_ids)

  train_raw <- cbind(diagnoses[train_ids], raw_relabd[train_ids, ]) |>
    as.data.frame()
  colnames(train_raw)[1] <- "diagnoses"
  train_raw$diagnoses <- as.factor(train_raw$diagnoses)

  test_raw <- cbind(diagnoses[test_ids], raw_relabd[test_ids, ]) |>
    as.data.frame()
  colnames(test_raw)[1] <- "diagnoses"
  test_raw$diagnoses <- as.factor(test_raw$diagnoses)

  train_denoised <- cbind(diagnoses[train_ids], denoised_relabd[train_ids, ]) |>
    as.data.frame()
  colnames(train_denoised)[1] <- "diagnoses"
  train_denoised$diagnoses <- as.factor(train_denoised$diagnoses)

  test_denoised <- cbind(diagnoses[test_ids], denoised_relabd[test_ids, ]) |>
    as.data.frame()
  colnames(test_denoised)[1] <- "diagnoses"
  test_denoised$diagnoses <- as.factor(test_denoised$diagnoses)

  train_control <- trainControl(method = "cv", number = 5)  # 5-fold CV

  tune_grid <- expand.grid(mtry = c(2, 5, 10, 15))

  rf_raw <- train(diagnoses ~ ., data = train_raw,
                  method = "rf",
                  trControl = train_control,
                  tuneGrid = tune_grid,
                  metric = "Accuracy")

  predictions_raw <- predict(rf_raw, newdata = test_raw,
                             type="prob")
  roc_raw <- roc(test_raw$diagnoses, predictions_raw[, 2])
  raw_auc[j] <- roc_raw$auc

  rf_denoised <- train(diagnoses ~ ., data = train_denoised,
                       method = "rf",
                       trControl = train_control,
                       tuneGrid = tune_grid,
                       metric = "Accuracy")

  predictions_denoised <- predict(rf_denoised, newdata = test_denoised,
                                  type="prob")
  roc_denoised <- roc(test_denoised$diagnoses, predictions_denoised[, 2])
  denoised_auc[j] <- roc_denoised$auc

}








