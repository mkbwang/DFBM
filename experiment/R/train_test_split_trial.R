

# load all the compositions
flocyt_composition <- read.csv("inst/extdata/flocyt_proj_composition_leaf_merged.csv")


# split into train and test
for (j in 1:20){
  print(j)
  set.seed(j)
  train_ids <- sample(seq(1:nrow(flocyt_composition)),
                      0.7*nrow(flocyt_composition))
  test_ids <- setdiff(seq(1:nrow(flocyt_composition)), train_ids)

  train_compositions <- flocyt_composition[train_ids, ]
  test_compositions <- flocyt_composition[test_ids, ]


  begin <- proc.time()
  train_dfbm_result <- dfbm(A=train_compositions, alpha=0.55)
  end <- proc.time()
  train_dfbm_denoised <- train_dfbm_result$denoised |> as.data.frame()


  test_dfbm_result <- predict.dfbm(train_dfbm_result, newdata=test_compositions)
  test_dfbm_denoised <- test_dfbm_result$denoised |> as.data.frame()


  denoise_result <- list(train_ids=train_ids, test_ids=test_ids,
                         train_original=train_compositions,
                         test_original=test_compositions,
                         train_denoised=train_dfbm_denoised,
                         test_denoised=test_dfbm_denoised)
  saveRDS(denoise_result, sprintf("experiment/data/denoise_flocyt_merged_%d.rds", j))

}




flocyt_composition <- read.csv("inst/extdata/flocyt_proj_composition_leaf.csv")
for (j in 1:20){
  print(j)
  set.seed(j)
  train_ids <- sample(seq(1:nrow(flocyt_composition)),
                      0.7*nrow(flocyt_composition))
  test_ids <- setdiff(seq(1:nrow(flocyt_composition)), train_ids)

  train_compositions <- flocyt_composition[train_ids, ]
  test_compositions <- flocyt_composition[test_ids, ]


  begin <- proc.time()
  train_dfbm_result <- dfbm(A=train_compositions, alpha=0.55)
  end <- proc.time()
  train_dfbm_denoised <- train_dfbm_result$denoised |> as.data.frame()


  test_dfbm_result <- predict.dfbm(train_dfbm_result, newdata=test_compositions)
  test_dfbm_denoised <- test_dfbm_result$denoised |> as.data.frame()


  denoise_result <- list(train_ids=train_ids, test_ids=test_ids,
                         train_original=train_compositions,
                         test_original=test_compositions,
                         train_denoised=train_dfbm_denoised,
                         test_denoised=test_dfbm_denoised)
  saveRDS(denoise_result, sprintf("experiment/data/denoise_flocyt_unmerged_%d.rds", j))

}



