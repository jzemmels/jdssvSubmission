# This script assumes preProcessing_script.R has been run and the scans are
# pre-processed and available in the .

library(tidyverse)
library(cmcR)
# devtools::install_github("jzemmels/scored")
library(scored)

# get the names of the scans saved in the preprocessed folder

allScans <- data.frame(file = list.files("DFSC-scans-preprocessed",recursive = TRUE,pattern = "\\.x3p")) %>%
  dplyr::mutate(scanName = file %>%
                  stringr::str_remove("G[0-9]{3}/") %>%
                  stringr::str_remove("\\.x3p") %>%
                  stringr::str_remove("-[0-9]{1,}"))

# this function computes the ACES features for a particular comparison of two
# scans
calculateFeatures <- function(compName,
                              x3p_df,
                              registFolder = "trainingRegistrations",
                              featureFolder = "trainingFeatures",
                              numCells = c(4,4),
                              eps = 3:15,
                              minPts = 3:10,
                              threshold = function(x3p1,x3p2){sd(abs(x3p1$surface.matrix - x3p2$surface.matrix),na.rm = TRUE)}){

  scanNames <- compName %>%
    str_remove("\\.x3p") %>%
    str_split("_vs_") %>%
    .[[1]]

  refName <- scanNames[1]
  targName <- scanNames[2]

  refScan <- x3p_df %>%
    filter(scanName == refName) %>%
    pull(processedScan) %>%
    .[[1]]

  targScan <- x3p_df %>%
    filter(scanName == targName) %>%
    pull(processedScan) %>%
    .[[1]]

  if(!isTRUE(all.equal(refScan$header.info$incrementX,targScan$header.info$incrementX))){

    if(refScan$header.info$incrementX > targScan$header.info$incrementX){

      targScan <- x3ptools::x3p_interpolate(targScan,resx = refScan$header.info$incrementX)

    }
    else{

      refScan <- x3ptools::x3p_interpolate(refScan,resx = targScan$header.info$incrementX)

    }
  }

  fullScanRegistrations <- scored::comparison_fullScan(refScan,targScan) %>%
    group_by(direction) %>%
    filter(fft_ccf == max(fft_ccf)) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(comparisonName = paste0(refName,"_vs_",targName)) %>%
    arrange(direction) %>%
    mutate(cellHeightValues = map(cellHeightValues,
                                  function(dat){

                                    if(!all(is.na(dat))){
                                      dat$surface.matrix <- dat$surface.matrix*dat$cmcR.info$scaleByVal # convert to micron scale
                                    }

                                    return(dat)
                                  }),
           alignedTargetCell = map(alignedTargetCell,
                                   function(dat){

                                     if(!all(is.na(dat))){
                                       dat$surface.matrix <- dat$surface.matrix*dat$cmcR.info$scaleByVal # convert to micron scale
                                     }

                                     return(dat)
                                   }))

  # compute registration and visual diagnostic features from full scans
  fullScanFeatures <- fullScanRegistrations %>%
    group_by(comparisonName,direction) %>%
    feature_aLaCarte(threshold = threshold) %>%
    group_by(comparisonName) %>%
    summarize(across(tidyselect::where(is.numeric),~ mean(.,na.rm = TRUE))) %>%
    set_names(paste0("fullScan_",names(.))) %>%
    rename(comparisonName = fullScan_comparisonName)

  # compute cell-based registrations
  cellBasedRegistrations <-
    bind_rows(scored::comparison_cellBased(reference = fullScanRegistrations$cellHeightValues[[1]],
                                           target = fullScanRegistrations$alignedTargetCell[[1]],
                                           direction = "one",
                                           numCells = numCells,
                                           maxMissingProp = .99,
                                           sideLengthMultiplier = 1.1,
                                           returnX3Ps = TRUE,
                                           thetas = -2:2) %>%
                mutate(direction = "reference_vs_target"),
              scored::comparison_cellBased(reference = fullScanRegistrations$cellHeightValues[[2]],
                                           target = fullScanRegistrations$alignedTargetCell[[2]],
                                           direction = "one",
                                           numCells = numCells,
                                           maxMissingProp = .99,
                                           sideLengthMultiplier = 1.1,
                                           returnX3Ps = TRUE,
                                           thetas = -2:2) %>%
                mutate(direction = "target_vs_reference")) %>%
    mutate(comparisonName = paste0(refName,"_vs_",targName)) %>%
    arrange(direction,theta) %>%
    group_by(direction,cellIndex) %>%
    # to save on space, only keep x3ps associated with the max CCF for each cell index
    mutate(cellHeightValues = ifelse(fft_ccf == max(fft_ccf),cellHeightValues,NA),
           alignedTargetCell = ifelse(fft_ccf == max(fft_ccf),alignedTargetCell,NA)) %>%
    ungroup()

  # compute registration and visual diagnsotics
  cellBasedFeatures_registrationVisual <- cellBasedRegistrations %>%
    filter(!is.na(cellHeightValues)) %>%
    mutate(cellHeightValues = map(cellHeightValues,
                                  function(dat){

                                    if(!all(is.na(dat))){
                                      dat$surface.matrix <- dat$surface.matrix*dat$cmcR.info$scaleByVal # convert to micron scale
                                    }

                                    return(dat)
                                  }),
           alignedTargetCell = map(alignedTargetCell,
                                   function(dat){

                                     if(!all(is.na(dat))){
                                       dat$surface.matrix <- dat$surface.matrix*dat$cmcR.info$scaleByVal # convert to micron scale
                                     }

                                     return(dat)
                                   })) %>%
    group_by(comparisonName,direction) %>%
    scored::feature_aLaCarte(features = c("registration","visual"),quiet = TRUE,
                             threshold = threshold) %>%
    group_by(comparisonName) %>%
    summarize(across(tidyselect::where(is.numeric),~ mean(.)))

  # compute the densty-based features over a grid of eps and minPts values
  cellBasedFeatures_density <-
    pmap_dfr(expand_grid("eps" = eps,
                         "minPts" = eps),
             ~ {

               cellBasedRegistrations %>%
                 select(-c(cellHeightValues,alignedTargetCell)) %>%
                 group_by(comparisonName,direction) %>%
                 feature_aLaCarte(features = c("density"),quiet = TRUE,
                                  eps = ..1,minPts = ..2) %>%
                 ungroup() %>%
                 mutate(eps = ..1,minPts = ..2,
                        thetaDiff = as.numeric(thetaDiff),
                        translationDiff = as.numeric(translationDiff),
                        clusterSize = as.numeric(clusterSize))

             }) %>%
    group_by(eps,minPts) %>%
    summarize(across(tidyselect::where(is.numeric),~ mean(.))) %>%
    mutate(clusterInd = !is.na(clusterSize))

  # combine the full scan and cell-based registrations into one data frame
  allRegistrations <-
    bind_rows(fullScanRegistrations %>%
                mutate(comparisonType = "fullScan"),
              cellBasedRegistrations %>%
                mutate(comparisonType = "cellBased"))

  # combine the cell-based visual diagnostic, registration, and density based
  # features
  cellBasedFeatures <- bind_cols(cellBasedFeatures_registrationVisual,
                                 cellBasedFeatures_density) %>%
    set_names(paste0("cellBased_",names(.))) %>%
    select(-cellBased_comparisonName)

  # combine the full scan and cell-based features into one data frame and add
  # the ground-truth to the comparison
  allFeatures <-
    bind_cols(fullScanFeatures,cellBasedFeatures)%>%
    distinct() %>%
    tidyr::separate(col = comparisonName,into = c("reference","target"),sep = "_vs_",remove = FALSE) %>%
    mutate(refBarrel = str_sub(reference,-2,-2),
           targBarrel = str_sub(target,-2,-2),
           type = ifelse(refBarrel == targBarrel,"match","non.match")) %>%
    select(-c(refBarrel,targBarrel,reference,target)) %>%
    select(comparisonName,type,
           cellBased_eps_4x4,cellBased_minPts_4x4,
           contains("fullScan_"),contains("cellBased_"))

  # create the specified registration/feature folders if they don't currently
  # exist
  if(!dir.exists(registFolder)){
    dir.create(registFolder)
  }
  if(!dir.exists(featureFolder)){
    dir.create(featureFolder)
  }

  saveRDS(allRegistrations,file = paste0(registFolder,"/",refName,"_vs_",targName,".RData"))
  saveRDS(allFeatures,file = paste0(featureFolder,"/",refName,"_vs_",targName,".RData"))

}

# we selected a subset of 210 scans for training the models, which are
# identified here
trainingScans <- c("K002eG1", "K002eG2", "K002eG3", "K002yF1", "K002yF2", "K002yF3",
                   "K009gG1", "K009gG2", "K009gG3", "K009mF1", "K009mF2", "K009mF3",
                   "K011dU1", "K011dU2", "K011dU3", "K011sR1", "K011sR2", "K011sR3",
                   "K012eN1", "K012eN2", "K012eN3", "K013mB1", "K013mB2", "K013mB3",
                   "K013pC1", "K013pC2", "K013pC3", "K013sA1", "K013sA2", "K013sA3",
                   "K014fB1", "K014fB2", "K014fB3", "K014pC1", "K014pC2", "K014pC3",
                   "K014vA1", "K014vA2", "K014vA3", "K015iZ1", "K015iZ2", "K015iZ3",
                   "K015rW1", "K015rW2", "K015rW3", "K016sF1", "K016sF2", "K016sF3",
                   "K016xG1", "K016xG2", "K016xG3", "K017dW1", "K017dW2", "K017dW3",
                   "K017xZ1", "K017xZ2", "K017xZ3", "K018eR1", "K018eR2", "K018eR3",
                   "K018jU1", "K018jU2", "K018jU3", "K019rF1", "K019rF2", "K019rF3",
                   "K019vG1", "K019vG2", "K019vG3", "K027gA1", "K027gA2", "K027gA3",
                   "K027uB1", "K027uB2", "K027uB3", "K027xC1", "K027xC2", "K027xC3",
                   "K028yF1", "K028yF2", "K028yF3", "K028zG1", "K028zG2", "K028zG3",
                   "K031fZ1", "K031fZ2", "K031fZ3", "K031zW1", "K031zW2", "K031zW3",
                   "K032mU1", "K032mU2", "K032mU3", "K032yR1", "K032yR2", "K032yR3",
                   "K036vN1", "K036vN2", "K036vN3", "K038gR1", "K038gR2", "K038gR3",
                   "K038jU1", "K038jU2", "K038jU3", "K043jR1", "K043jR2", "K043jR3",
                   "K043yU1", "K043yU2", "K043yU3", "K044gR1", "K044gR2", "K044gR3",
                   "K044uU1", "K044uU2", "K044uU3", "K045uW1", "K045uW2", "K045uW3",
                   "K045yZ1", "K045yZ2", "K045yZ3", "K046fU1", "K046fU2", "K046fU3",
                   "K046uR1", "K046uR2", "K046uR3", "K047pA1", "K047pA2", "K047pA3",
                   "K047rB1", "K047rB2", "K047rB3", "K047sC1", "K047sC2", "K047sC3",
                   "K048fZ1", "K048fZ2", "K048fZ3", "K048vW1", "K048vW2", "K048vW3",
                   "K049jC1", "K049jC2", "K049jC3", "K049pA1", "K049pA2", "K049pA3",
                   "K049uB1", "K049uB2", "K049uB3", "K050mR1", "K050mR2", "K050mR3",
                   "K050zU1", "K050zU2", "K050zU3", "K051iG1", "K051iG2", "K051iG3",
                   "K051vF1", "K051vF2", "K051vF3", "K054jW1", "K054jW2", "K054jW3",
                   "K054rZ1", "K054rZ2", "K054rZ3", "K055pR1", "K055pR2", "K055pR3",
                   "K055zU1", "K055zU2", "K055zU3", "K060eA1", "K060eA2", "K060eA3",
                   "K060fC1", "K060fC2", "K060fC3", "K060xB1", "K060xB2", "K060xB3",
                   "K061jW1", "K061jW2", "K061jW3", "K061zZ1", "K061zZ2", "K061zZ3",
                   "K065jN1", "K065jN2", "K065jN3", "K066yN1", "K066yN2", "K066yN3",
                   "K067dG1", "K067dG2", "K067dG3", "K067uF1", "K067uF2", "K067uF3",
                   "K155dG1", "K155dG2", "K155dG3", "K155sF1", "K155sF2", "K155sF3"
)

trainingScan_df <- allScans %>%
  dplyr::filter(scanName %in% trainingScans) %>%
  mutate(processedScan = map(file,~ {

    return(x3ptools::x3p_read(paste0("DFSC-scans-preprocessed/G",str_sub(file,2,4),"/",.)))

  }))

# create a data frame of all pairwise comparisons between the training scans
comparisonNames_train <-
  expand_grid(source = trainingScan_df$scanName,
              target = trainingScan_df$scanName) %>%
  left_join(trainingScan_df %>%
              select(-c(processedScan,file)) %>%
              mutate(scanInd = 1:nrow(.)),
            by = c("source" = "scanName")) %>%
  rename(sourceInd = scanInd) %>%
  left_join(trainingScan_df %>%
              select(-c(processedScan,file)) %>%
              mutate(scanInd = 1:nrow(.)),
            by = c("target" = "scanName")) %>%
  rename(targetInd = scanInd) %>%
  filter(sourceInd < targetInd) %>%
  select(-c(sourceInd,targetInd)) %>%
  mutate(comparisonName = paste0(source,"_vs_",target))


# future:::ClusterRegistry("stop") # this function will kill any currently running parallel jobs
future::plan(future::multisession(workers = round(future::availableCores()/2)))

comparisonNames_train %>%
  pull(comparisonName) %>%
  furrr::future_walk(~ calculateFeatures(compName = .,
                                         x3p_df = trainingScan_df,
                                         registFolder = "trainingRegistrations",
                                         featureFolder = "trainingFeatures"))

# now repeat the process, but this time for the 300 test scans.

testScans <- c("K002fJ1", "K002fJ2", "K002fJ3", "K002mH1", "K002mH2", "K002mH3",
               "K002uK1", "K002uK2", "K002uK3", "K009dK1", "K009dK2", "K009dK3",
               "K009eH1", "K009eH2", "K009eH3", "K009sJ1", "K009sJ2", "K009sJ3",
               "K011uV1", "K011uV2", "K011uV3", "K011yT1", "K011yT2", "K011yT3",
               "K011zS1", "K011zS2", "K011zS3", "K012dM1", "K012dM2", "K012dM3",
               "K012sQ1", "K012sQ2", "K012sQ3", "K012uP1", "K012uP2", "K012uP3",
               "K012zL1", "K012zL2", "K012zL3", "K013uD1", "K013uD2", "K013uD3",
               "K013yE1", "K013yE2", "K013yE3", "K014rD1", "K014rD2", "K014rD3",
               "K014sE1", "K014sE2", "K014sE3", "K015d%1", "K015d%2", "K015d%3",
               "K015eY1", "K015eY2", "K015eY3", "K015jX1", "K015jX2", "K015jX3",
               "K016dJ1", "K016dJ2", "K016dJ3", "K016gH1", "K016gH2", "K016gH3",
               "K016uK1", "K016uK2", "K016uK3", "K017g%1", "K017g%2", "K017g%3",
               "K017uY1", "K017uY2", "K017uY3", "K017vX1", "K017vX2", "K017vX3",
               "K018fT1", "K018fT2", "K018fT3", "K018mV1", "K018mV2", "K018mV3",
               "K018vS1", "K018vS2", "K018vS3", "K019gK1", "K019gK2", "K019gK3",
               "K019sH1", "K019sH2", "K019sH3", "K019xJ1", "K019xJ2", "K019xJ3",
               "K027jD1", "K027jD2", "K027jD3", "K027yE1", "K027yE2", "K027yE3",
               "K028fJ1", "K028fJ2", "K028fJ3", "K028jK1", "K028jK2", "K028jK3",
               "K028pH1", "K028pH2", "K028pH3", "K031d%1", "K031d%2", "K031d%3",
               "K031eY1", "K031eY2", "K031eY3", "K031xX1", "K031xX2", "K031xX3",
               "K032dV1", "K032dV2", "K032dV3", "K032jT1", "K032jT2", "K032jT3",
               "K032zS1", "K032zS2", "K032zS3", "K036gP1", "K036gP2", "K036gP3",
               "K036jL1", "K036jL2", "K036jL3", "K036uM1", "K036uM2", "K036uM3",
               "K036xQ1", "K036xQ2", "K036xQ3", "K038mV1", "K038mV2", "K038mV3",
               "K038yS1", "K038yS2", "K038yS3", "K038zT1", "K038zT2", "K038zT3",
               "K043dV1", "K043dV2", "K043dV3", "K043sS1", "K043sS2", "K043sS3",
               "K043zT1", "K043zT2", "K043zT3", "K044iV1", "K044iV2", "K044iV3",
               "K044pS1", "K044pS2", "K044pS3", "K044sT1", "K044sT2", "K044sT3",
               "K045fY1", "K045fY2", "K045fY3", "K045pX1", "K045pX2", "K045pX3",
               "K045r%1", "K045r%2", "K045r%3", "K046dV1", "K046dV2", "K046dV3",
               "K046iS1", "K046iS2", "K046iS3", "K046rT1", "K046rT2", "K046rT3",
               "K047dE1", "K047dE2", "K047dE3", "K047eD1", "K047eD2", "K047eD3",
               "K048pY1", "K048pY2", "K048pY3", "K048uX1", "K048uX2", "K048uX3",
               "K048y%1", "K048y%2", "K048y%3", "K049iD1", "K049iD2", "K049iD3",
               "K049sE1", "K049sE2", "K049sE3", "K050eT1", "K050eT2", "K050eT3",
               "K050xV1", "K050xV2", "K050xV3", "K050yS1", "K050yS2", "K050yS3",
               "K051eK1", "K051eK2", "K051eK3", "K051gJ1", "K051gJ2", "K051gJ3",
               "K051mH1", "K051mH2", "K051mH3", "K054fY1", "K054fY2", "K054fY3",
               "K054vX1", "K054vX2", "K054vX3", "K054y%1", "K054y%2", "K054y%3",
               "K055gS1", "K055gS2", "K055gS3", "K055jV1", "K055jV2", "K055jV3",
               "K055yT1", "K055yT2", "K055yT3", "K060rE1", "K060rE2", "K060rE3",
               "K060zD1", "K060zD2", "K060zD3", "K061gY1", "K061gY2", "K061gY3",
               "K061pX1", "K061pX2", "K061pX3", "K061x%1", "K061x%2", "K061x%3",
               "K065dQ1", "K065dQ2", "K065dQ3", "K065iM1", "K065iM2", "K065iM3",
               "K065mL1", "K065mL2", "K065mL3", "K065zP1", "K065zP2", "K065zP3",
               "K066fP1", "K066fP2", "K066fP3", "K066iL1", "K066iL2", "K066iL3",
               "K066rQ1", "K066rQ2", "K066rQ3", "K066zM1", "K066zM2", "K066zM3",
               "K067gH1", "K067gH2", "K067gH3", "K067rK1", "K067rK2", "K067rK3",
               "K067sJ1", "K067sJ2", "K067sJ3", "K155gK1", "K155gK2", "K155gK3",
               "K155jH1", "K155jH2", "K155jH3", "K155xJ1", "K155xJ2", "K155xJ3"
)

comparisonNames_test <-
  expand_grid(source = testScans,
              target = testScans) %>%
  left_join(data.frame(scanName = testScans,
                       scanInd = 1:length(testScans)),
            by = c("source" = "scanName")) %>%
  rename(sourceInd = scanInd) %>%
  left_join(data.frame(scanName = testScans,
                       scanInd = 1:length(testScans)),
            by = c("target" = "scanName")) %>%
  rename(targetInd = scanInd) %>%
  filter(sourceInd < targetInd) %>%
  select(-c(sourceInd,targetInd)) %>%
  mutate(comparisonName = paste0(source,"_vs_",target))

testScan_df <- allScans %>%
  dplyr::filter(scanName %in% testScans) %>%
  mutate(processedScan = map(file,~ {

    return(x3ptools::x3p_read(paste0("DFSC-scans-preprocessed/G",str_sub(file,2,4),"/",.)))

  }))


# future:::ClusterRegistry("stop") # this function will kill any currently running parallel jobs
future::plan(future::multisession(workers = round(future::availableCores()/2)))

comparisonNames_test %>%
  pull(comparisonName) %>%
  furrr::future_walk(~ calculateFeatures(compName = .,
                                         x3p_df = testScan_df,
                                         eps = 3,
                                         minPts = 3,
                                         numCells = c(4,4),
                                         registFolder = "testRegistrations",
                                         featureFolder = "testsFeatures"))
