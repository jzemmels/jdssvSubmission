## This script fits various models to the 9 visual diagnostic features. These
## models are explored in used in Chapter 3 of
## <https://github.com/jzemmels/cartridgeCaseLitReview/blob/main/docs/thesis.pdf>

## Fit binary classifier models to visual diagnostic data, consider
## upsampling/downsampling to balance matches/non-matches
library(tidyverse)
library(caret)

# collect the visual diagnostic features
visualDiagFeatures <- map_dfr(list.files("trainingFeatures/",full.names = TRUE),
        ~ {

          return(readRDS(.))

        }) %>%
  distinct() %>%
  group_by(comparisonName,type,cellBased_eps,cellBased_minPts) %>%
  summarize(across(tidyselect::where(is.numeric),~ mean(.))) %>%
  mutate(cellBased_clusterInd = !is.na(cellBased_clusterSize)) %>%
  ungroup() %>%
  distinct() %>%
  select(-c(contains("fullScan_neighborhoodSizeAve_sd"),contains("fullScan_neighborhoodSizeSD_sd"),
            contains("fullScan_differenceCor_sd"),contains("fullScan_filteredRatio_sd"),
            contains("cellBased_neighborhoodSizeAve_sd"),contains("cellBased_neighborhoodSizeSD_sd"),
            contains("cellBased_differenceCor_sd"),
            contains("fullScan_ccfMean"),contains("cellBased_ccfMean"),contains("cellBased_ccfSD"))) %>%
  select(comparisonName,type,
         contains("differenceCor"),
         contains("filteredRatio"),
         contains("neighborhoodSize")) %>%
  distinct()

saveRDS(visualDiagFeatures,file = "~/topMatchData_210scanTrain/visualDiagFeatures.rds")

# impute missing values for the visual diagnostic features
visualDiag_imputed <- visualDiagFeatures %>%
  mutate(
    across(everything(),
           .fns = ~ ifelse(is.na(.x) | is.infinite(.x) | is.nan(.x),
                           median(.x,na.rm = TRUE),.x))
  )

# fit 27 models (3 feaure groups x 3 classifier models x 3 sampling schemes)
visualDiag_fittedModels <- map2_dfr(list(visualDiag_imputed %>%
                                           select(type,contains("fullScan")),
                                         visualDiag_imputed %>%
                                           select(type,contains("cellBased")),
                                         visualDiag_imputed),
                                    c("fullScan","cellBased","all"),
                                    function(caretDat,featureGroup){

                                      lr_upsample <- caret::train(form = type ~ .,
                                                                  data = caretDat,
                                                                  trControl = trainControl(method = "repeatedcv",
                                                                                           number = 10,
                                                                                           repeats = 3,
                                                                                           classProbs = TRUE,
                                                                                           summaryFunction = twoClassSummary,
                                                                                           sampling = 'up',
                                                                                           seeds = NULL,
                                                                                           savePredictions = TRUE),
                                                                  method = "glm",
                                                                  family = "binomial")

                                      rf_upsample <- caret::train(form = type ~ .,
                                                                  data = caretDat,
                                                                  trControl = trainControl(method = "repeatedcv",
                                                                                           number = 10,
                                                                                           repeats = 3,
                                                                                           classProbs = TRUE,
                                                                                           summaryFunction = twoClassSummary,
                                                                                           sampling = 'up',
                                                                                           savePredictions = TRUE),
                                                                  tuneGrid = data.frame(mtry = c(2,
                                                                                                 floor(sqrt(ncol(caretDat) - 1)),
                                                                                                 ceiling(sqrt(ncol(caretDat) - 1)),
                                                                                                 (ncol(caretDat) - 1))),
                                                                  method = "rf")

                                      cart_upsample <- caret::train(form = type ~ .,
                                                                    data = caretDat,
                                                                    trControl = trainControl(method = "repeatedcv",
                                                                                             number = 10,
                                                                                             repeats = 3,
                                                                                             classProbs = TRUE,
                                                                                             summaryFunction = twoClassSummary,
                                                                                             sampling = 'up',
                                                                                             savePredictions = TRUE),
                                                                    method = "rpart")

                                      lr_downsample <- caret::train(form = type ~ .,
                                                                    data = caretDat,
                                                                    trControl = trainControl(method = "repeatedcv",
                                                                                             number = 10,
                                                                                             repeats = 3,
                                                                                             classProbs = TRUE,
                                                                                             summaryFunction = twoClassSummary,
                                                                                             sampling = 'down',
                                                                                             seeds = NULL,
                                                                                             savePredictions = TRUE),
                                                                    method = "glm",
                                                                    family = "binomial")

                                      rf_downsample <- caret::train(form = type ~ .,
                                                                    data = caretDat,
                                                                    trControl = trainControl(method = "repeatedcv",
                                                                                             number = 10,
                                                                                             repeats = 3,
                                                                                             classProbs = TRUE,
                                                                                             summaryFunction = twoClassSummary,
                                                                                             sampling = 'down',
                                                                                             savePredictions = TRUE),
                                                                    tuneGrid = data.frame(mtry = c(2,
                                                                                                   floor(sqrt(ncol(caretDat) - 1)),
                                                                                                   ceiling(sqrt(ncol(caretDat) - 1)),
                                                                                                   (ncol(caretDat) - 1))),
                                                                    method = "rf")

                                      cart_downsample <- caret::train(form = type ~ .,
                                                                      data = caretDat,
                                                                      trControl = trainControl(method = "repeatedcv",
                                                                                               number = 10,
                                                                                               repeats = 3,
                                                                                               classProbs = TRUE,
                                                                                               summaryFunction = twoClassSummary,
                                                                                               sampling = 'down',
                                                                                               savePredictions = TRUE),
                                                                      method = "rpart")

                                      lr_nosample <- caret::train(form = type ~ .,
                                                                  data = caretDat,
                                                                  trControl = trainControl(method = "repeatedcv",
                                                                                           number = 10,
                                                                                           repeats = 3,
                                                                                           classProbs = TRUE,
                                                                                           summaryFunction = twoClassSummary,
                                                                                           # sampling = 'none',
                                                                                           seeds = NULL,
                                                                                           savePredictions = TRUE),
                                                                  method = "glm",
                                                                  family = "binomial")

                                      rf_nosample <- caret::train(form = type ~ .,
                                                                  data = caretDat,
                                                                  trControl = trainControl(method = "repeatedcv",
                                                                                           number = 10,
                                                                                           repeats = 3,
                                                                                           classProbs = TRUE,
                                                                                           summaryFunction = twoClassSummary,
                                                                                           savePredictions = TRUE),
                                                                  tuneGrid = data.frame(mtry = c(2,
                                                                                                 floor(sqrt(ncol(caretDat) - 1)),
                                                                                                 ceiling(sqrt(ncol(caretDat) - 1)),
                                                                                                 (ncol(caretDat) - 1))),
                                                                  method = "rf")

                                      cart_nosample <- caret::train(form = type ~ .,
                                                                    data = caretDat,
                                                                    trControl = trainControl(method = "repeatedcv",
                                                                                             number = 10,
                                                                                             repeats = 3,
                                                                                             classProbs = TRUE,
                                                                                             summaryFunction = twoClassSummary,
                                                                                             savePredictions = TRUE),
                                                                    method = "rpart")

                                      tibble(featureGroup = featureGroup,
                                             modelType = rep(c("lr","rf","cart"),each = 3),
                                             sampling = rep(c("none","up","down"),times = 3),
                                             fittedModel =
                                               list(lr_nosample,
                                                    lr_upsample,
                                                    lr_downsample,
                                                    rf_nosample,
                                                    rf_upsample,
                                                    rf_downsample,
                                                    cart_nosample,
                                                    cart_upsample,
                                                    cart_downsample))

                                    })

saveRDS(visualDiag_fittedModels,file = "~/topMatchData_210scanTrain/visualDiag_fittedModels.rds")

visualDiag_fittedModelsResults <- visualDiag_fittedModels %>%
  pmap_dfr(~ {

    data.frame(featureGroup = ..1,
               modelType = ..2,
               sampling = ..3,
               ROC = max(..4$results$ROC))

  })

visualDiag_fittedModelsResults %>%
  arrange(desc(ROC))
