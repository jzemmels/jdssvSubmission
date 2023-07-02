# this script fits various models to the features computed by running
# comparing_script.R -- we assume that both preProcesssing_script.R and
# comparing_script.R have already been run.

# Recall from comparing_script.R that we computed the training features across a
# grid of epsilon and minPts values. We will want to fit a model for each of
# these grid combinations. For the sake of space (especifally when we
# parallelize), we split the training features based on the eps/minPts
# combinations and save each data set to a distinct file.

library(tidyverse)
# install.packages("caret")
library(caret)

# first read the full "trainingFeatures" data set into memory
trainingFeatures <- map_dfr(list.files("trainingFeatures/",full.names = TRUE),
                         ~ {

                           return(readRDS(.))

                         }) %>%
  distinct() %>%
  group_by(comparisonName,type,cellBased_eps,cellBased_minPts) %>%
  summarize(across(tidyselect::where(is.numeric),~ mean(.))) %>%
  mutate(cellBased_clusterInd = !is.na(cellBased_clusterSize)) %>%
  ungroup() %>%
  distinct()

# then, split the trainingFeatures data set based on the value of epsilon and
# minPts. We remove some features from the data set that we're not interesting
# in using to fit models
walk(trainingFeatures %>%
       group_by(cellBased_eps,cellBased_minPts) %>%
       group_split(),
     function(dat){

       ret <- dat %>%
         select(-c(contains("fullScan_neighborhoodSizeAve_sd"),contains("fullScan_neighborhoodSizeSD_sd"),
                   contains("fullScan_differenceCor_sd"),contains("fullScan_filteredRatio_sd"),
                   contains("cellBased_neighborhoodSizeAve_sd"),contains("cellBased_neighborhoodSizeSD_sd"),
                   contains("cellBased_differenceCor_sd"),
                   contains("fullScan_ccfMean"),contains("cellBased_ccfMean"),contains("cellBased_ccfSD")))

       saveRDS(ret,file =
                 paste0("trainingFeatures_split/trainingFeatures_split_eps",
                        unique(ret$cellBased_eps),"_minPts",
                        unique(ret$cellBased_minPts),".rds"))

     })

# first, we fit logistic regression models to each training data set using all
# 19 ACES features.
future:::ClusterRegistry("stop")

future::plan(future::multisession(workers = round(future::availableCores()/2)))

set.seed(442023)

lr_allACES <- furrr::future_map_dfr(
  list.files("trainingFeatures_split/",full.names = TRUE),
  function(fileName){

    dat <- readRDS(fileName)

    caretDat <- dat %>%
      mutate(
        across(everything(),
               .fns = ~ ifelse(is.na(.x) | is.infinite(.x) | is.nan(.x),
                               median(.x,na.rm = TRUE),.x))
      ) %>%
      select(-c(comparisonName,cellBased_eps,cellBased_minPts))

    lr <- caret::train(form = type ~ .,
                       data = caretDat,
                       trControl = trainControl(method = "repeatedcv",
                                                number = 10,
                                                repeats = 3,
                                                classProbs = TRUE,
                                                summaryFunction = twoClassSummary,
                                                seeds = NULL,
                                                savePredictions = TRUE),
                       method = "glm",
                       family = "binomial")

    tibble(cellBased_eps = unique(dat$cellBased_eps),
           cellBased_minPts = unique(dat$cellBased_minPts),
           lr = list(lr))

  })

if(!dir.exists("models")){
  dir.create("models")
}
saveRDS(lr_allACES,file = "models/lr_allACES.rds")

# repeat the logistic regression fit, but now only use the C_0 &
# registration-based features

future:::ClusterRegistry("stop")
future::plan(future::multisession(workers = future::availableCores() - 5))

set.seed(442023)

lr_baselineRegistration <-
  furrr::future_map_dfr(
    list.files("trainingFeatures_split/",full.names = TRUE),
    function(fileName){

      dat <- readRDS(fileName)

      caretDat <- dat %>%
        select(comparisonName,type,cellBased_eps,cellBased_minPts,
               cellBased_clusterInd,
               cellBased_xTransSD,cellBased_yTransSD,cellBased_thetaRotSD,
               cellBased_pairwiseCompCorAve,cellBased_pairwiseCompCorSD,
               fullScan_pairwiseCompCorAve) %>%
        mutate(
          across(everything(),
                 .fns = ~ ifelse(is.na(.x) | is.infinite(.x) | is.nan(.x),
                                 median(.x,na.rm = TRUE),.x))
        ) %>%
        select(-c(comparisonName,cellBased_eps,cellBased_minPts))

      lr <- caret::train(form = type ~ .,
                         data = caretDat,
                         trControl = trainControl(method = "repeatedcv",
                                                  number = 10,
                                                  repeats = 3,
                                                  classProbs = TRUE,
                                                  summaryFunction = twoClassSummary,
                                                  seeds = NULL,
                                                  savePredictions = TRUE),
                         method = "glm",
                         family = "binomial")

      tibble(cellBased_eps = unique(dat$cellBased_eps),
             cellBased_minPts = unique(dat$cellBased_minPts),
             lr = list(lr))

    })

saveRDS(lr_baselineRegistration,file = "models/lr_baselineRegistration.rds")

# next, we fit random forest models using the full ACES and C_0 + Registration
# data sets

future:::ClusterRegistry("stop")

future::plan(future::multisession(workers = future::availableCores() - 5))

set.seed(442023)

rf_allACES <-
  furrr::future_map_dfr(
    list.files("trainingFeatures_split/",full.names = TRUE),
    function(fileName){

      dat <- readRDS(fileName)

      caretDat <- dat %>%
        mutate(
          across(everything(),
                 .fns = ~ ifelse(is.na(.x) | is.infinite(.x) | is.nan(.x),
                                 median(.x,na.rm = TRUE),.x))
        ) %>%
        select(-c(comparisonName,cellBased_eps,cellBased_minPts))

      rf <- caret::train(form = type ~ .,
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

      tibble(cellBased_eps = unique(dat$cellBased_eps),
             cellBased_minPts = unique(dat$cellBased_minPts),
             rf = list(rf))

    })

saveRDS(rf_allACES,file = "models/rf_allACES.rds")

rf_baselineRegistration <-
  furrr::future_map_dfr(
    list.files("trainingFeatures_split/",full.names = TRUE),
    function(fileName){

      dat <- readRDS(fileName)

      caretDat <- dat %>%
        select(comparisonName,type,cellBased_eps,cellBased_minPts,
               cellBased_clusterInd,
               cellBased_xTransSD,cellBased_yTransSD,cellBased_thetaRotSD,
               cellBased_pairwiseCompCorAve,cellBased_pairwiseCompCorSD,
               fullScan_pairwiseCompCorAve) %>%
        mutate(
          across(everything(),
                 .fns = ~ ifelse(is.na(.x) | is.infinite(.x) | is.nan(.x),
                                 median(.x,na.rm = TRUE),.x))
        ) %>%
        select(-c(comparisonName,cellBased_eps,cellBased_minPts))

      rf <- caret::train(form = type ~ .,
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

      tibble(cellBased_eps = unique(dat$cellBased_eps),
             cellBased_minPts = unique(dat$cellBased_minPts),
             rf = list(rf))

    })



saveRDS(rf_baselineRegistration,file = "models/rf_baselineRegistration.rds")
