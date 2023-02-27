library(tidymodels)
library(sf)
library(mapview)
library(raster)
library(here)
library(fasterize)
library(randomForest)
library(doParallel)
library(parallel)
library(CAST)
library(caret)
library(exactextractr)
library(plotly)
library(stringr)
library(svMisc)

areas <- read_sf(here("data","FieldData","buffered.geojson"))  %>% st_transform(4326)  %>% 
  dplyr::select(-lon, -lat) %>% mutate(province = case_when(area %in% c("Chokwe", "Xai-Xai") ~ "Gaza",TRUE ~ "Manica"))

trainingPoly <- st_read(here("data","FieldData", "newTD_20210727.geojson")) %>%
  dplyr::select(-area) %>% mutate(landcover_level2 = case_when(landcover_level2 == "Cropland irrigated 1" ~ "Cropland irrigated",
                                                               landcover_level2 == "Cropland irrigated 2" ~ "Cropland irrigated",
                                                               landcover_level2 == "Plantation" ~ "Cropland irrigated",
                                                               TRUE ~ landcover_level2),
                                  landcover_code = case_when( landcover_code == 0 ~ "1",
                                                              landcover_code == 7 ~ "2",
                                                              TRUE ~landcover_code)) %>%
  mutate(
    code_level1 = as.character(as.numeric(as.factor(landcover_level1))),
    code_level2 = as.character(as.numeric(as.factor(landcover_level2))),
    PolygonID = seq(1, nrow(.),1),
    DEA_id = 1:nrow(.)) %>% 
  filter(!st_is_empty(.)) %>%
  mutate(point_type = ifelse(is.na(point_type), "Manually drawn", point_type)) 

trainingPoly <- st_join(trainingPoly, areas) 
unique_values <- trainingPoly %>% st_drop_geometry() %>% summarise(landcover_level2  = unique(landcover_level2 ),
                                                                   code = unique(code_level2)) %>% arrange(code)
unique_values

names_S2 <- c('red', 'green', 'blue', 'nir', 'swir_1', 'swir_2', 'red_edge_1', 'red_edge_2', "smad", "emad", "bcmad", "NDVI", "BSI", "MNDWI", "CIRE")
names_S1 <- c("VV", "VH", 'RVI')

extract_pixel_function <- function(location){
  no_cores <- detectCores() - 2
  cl <- makeCluster(no_cores)
  registerDoParallel(cl)
  
  
  
  mo_6_Q41 <- brick(here("data","DEA_round2", paste0(location,"_6mo_2019Q4_2020Q1.tif"))) %>% dropLayer(12) #layer 12 is clear pixel count
  cire_6_Q41 <- (mo_6_Q41[[4]]/mo_6_Q41[[8]]) - 1
  mo_6_Q41 <- addLayer(mo_6_Q41, cire_6_Q41 ) 
  names(mo_6_Q41) <- paste0(names_S2, "_6m_Q41")
  print("mo_6_Q41 ready")
  
  mo_6_Q23 <- brick(here("data","DEA_round2", paste0(location,"_6mo_2020Q2_2020Q3.tif"))) %>% dropLayer(12) #layer 12 is clear pixel count
  cire_6_Q23 <- (mo_6_Q23[[4]]/mo_6_Q23[[8]]) - 1
  mo_6_Q23 <- addLayer(mo_6_Q23, cire_6_Q23 ) 
  names(mo_6_Q23) <- paste0(names_S2, "_6m_Q23")
  print("mo_6_Q23 ready")
  
  
  mo_6_Q41_SAR <- brick(here("data","DEA_round2",paste0(location,"_6mo_2019Q4_2020Q1_SAR.tif")))
  names(mo_6_Q41_SAR) <- paste0(names_S1, "_6m_Q41")
  print("mo_6_Q41_SAR ready")
  
  mo_6_Q23_SAR <- brick(here("data","DEA_round2",paste0(location,"_6mo_2020Q2_2020Q3_SAR.tif")))
  names(mo_6_Q23_SAR) <- paste0(names_S1, "_6m_Q23")
  print("mo_6_Q23_SAR ready")
  
  raster_ready_6m <- addLayer(mo_6_Q41,     mo_6_Q23, 
                              mo_6_Q41_SAR, mo_6_Q23_SAR)
  
  extracted_df <- exact_extract(raster_ready_6m , trainingPoly, force_df = TRUE)
  print("extraction ready")
  
  extracted_binded <- bind_rows(extracted_df, .id = "PolygonID")%>% 
    dplyr::select(-coverage_fraction) %>% na.omit()
  
  TD_df <- merge(extracted_binded, trainingPoly %>% st_drop_geometry(), by = "PolygonID")  %>% 
    select(all_of(names(raster_ready_6m)), "code_level2", "PolygonID")
  
  return(TD_df)

}

catandica_vals <- extract_pixel_function("Catandica")
manica_vals <- extract_pixel_function("Manica")
xai_vals <- extract_pixel_function("Xai-Xai")
chokwe_vals <- extract_pixel_function("Chokwe")

extracted_pixel_vals_all <- rbind(manica_vals, catandica_vals) # comment one of the two
extracted_pixel_vals_all <- rbind(xai_vals, chokwe_vals)

nrow(extracted_pixel_vals_all)

table(extracted_pixel_vals_all$code_level2)
set.seed(100)
polys_split <- initial_split(extracted_pixel_vals_all, prop = .8, code_level2) #prop defines the amount of split

TD_df <- training(polys_split)
nrow(TD_df)
nrow(TD_df)/nrow(extracted_pixel_vals_all)*100

pixels_per_class <- TD_df%>% 
  group_by(code_level2 ) %>% 
  summarise(total = n())
pixels_per_class
min(round(pixels_per_class$total,-2))

predictors <- paste0(names_S1, names_S2)
response <- "code_level2"

model_train <- function(trainDat, composite_lengths, algorithm){

  indices <- CreateSpacetimeFolds(trainDat,
                                  spacevar = "PolygonID",
                                  k=3,
                                  class="code_level2")
  trainDat <- trainDat %>% select(-PolygonID)
  trainDat <- mutate(trainDat, code_level2 = as.factor(code_level2))

  # set.seed = 100
  if(algorithm == "knn"){
    ctrl <- trainControl(method="cv", 
                         index = indices$index,
                         savePredictions = TRUE,
                         allowParallel= F,
                         number = 5, 
                         verboseIter = TRUE)
  }else{
    ctrl <- trainControl(method="cv", 
                         index = indices$index,
                         savePredictions = TRUE,
                         allowParallel= TRUE,
                         number = 5, 
                         verboseIter = TRUE)}
  
 if(algorithm == "nnet"){
    print("nnet model")
    model_ffs <- train(
      code_level2 ~ .,
      data = trainDat,
      method=algorithm, 
      metric="Accuracy",
      trControl=ctrl,
      importance=TRUE,
      withinSE = TRUE,
      tuneLength = 5,
      na.rm = TRUE ,
      preProcess = c("center", "scale"))
  }else if(algorithm == "svmRadial"){
    print("svm model")
    model_ffs <- train(
      code_level2 ~ .,
      data = trainDat,
      method=algorithm,
      metric="Accuracy",
      trControl=ctrl,
      importance=TRUE,
      withinSE = TRUE,
      tuneLength = 5,
      na.rm = TRUE ,
      preProcess = c("center", "scale"))
  }    else if(algorithm == "rf"){
    print( "rf model")
    model_ffs <- train(
      code_level2 ~ .,
      data = trainDat,
      method=algorithm, 
      metric="Accuracy",
      trControl=ctrl,
      importance=TRUE,
      withinSE = TRUE,
      tuneLength = 5,
      na.rm = TRUE )
  }
  return(model_ffs)
  
}

set.seed(100)
seed_vals <- sort(sample.int(1000, 25))

accu_vect <- data.frame(overall_accuracy = double(), 
                        overall_accuracy_model = double(), 
                        fraction = double(),
                        pixels = double(),
                        seed = double(), 
                        algo = character(), 
                        best_model = double(),
                        data_set = character())

accu_vect_irri <- data.frame(irri_UA = double(), 
                             irri_PA = double(),
                             irri_UA_model = double(), 
                             irri_PA_model = double(),
                             fraction = double(),
                             pixels_irri = double(),
                             seed = double(), 
                             algo = character(), 
                             best_model = double(),
                             data_set = character())

# scenario 1
fracs <- c(0.01, 0.05, 0.1,0.2,0.4,0.6,0.8,1)

# scenario 2
min_pixs <- min(round(pixels_per_class$total,-1))
max_pixels <- round((min_pixs*(length(unique(pixels_per_class$code_level2))-1))*100/99) #1% irri, 99% equally divided other classes
balanced_vals <- seq(50, min(round(pixels_per_class$total,-2)), round((min(round(pixels_per_class$total,-2))-50)/6))

# scenario 3
forced_imblace_irri <- round(c(0.01*max_pixels,0.05*max_pixels,0.1*max_pixels, 0.2*max_pixels, 0.5*max_pixels,0.8*max_pixels,0.9*max_pixels, 0.95*max_pixels, 0.99*max_pixels))
forced_imblace_rest <- round((max_pixels-forced_imblace_irri)/(length(unique(pixels_per_class$code_level2))-1))
forced_imblace_irri+(forced_imblace_rest*(length(unique(pixels_per_class$code_level2))-1)) #test if all sets have equal total numbers

# scenario 4
mislabelled_num <- round(c(nrow(TD_df)*0.01, nrow(TD_df)*0.05, nrow(TD_df)*0.1, nrow(TD_df)*0.2, nrow(TD_df)*0.4 ))
class_vals <- c("2", "3", "6")


repeat_classification_function <- function(data_set, algo, location){
  
  no_cores <- detectCores() - 2
  cl <- makeCluster(no_cores)
  registerDoParallel(cl)
  
  if(data_set == "balanced_vals") {
    data_set0 = balanced_vals
    print("balanced_vals")
  } else if(data_set == "fracs") {
    data_set0 = fracs
    print("fracs")
  } else if(data_set == "forced_imbalance") {
    data_set0 = forced_imblace_irri
    print("forced_imbalance")
  } else if(data_set == "misclassified"){
    data_set0 = mislabelled_num
    print("misclassified")
  }else{
    print("Wrong option, choose from `fracs`, `balanced_vals`,  `forced_imbalance` or 'misclassified'")}
  
  for(i in seq_along(data_set0)){

    for (j in seq_along(seed_vals)){
      set.seed(seed_vals[j])
      
      if(data_set == "balanced_vals") {
        trainDat <- TD_df %>% group_by(code_level2 ) %>% dplyr::sample_n(balanced_vals[i]) %>% as.data.frame()
      } else if(data_set == "fracs") {
        trainDat <- TD_df %>% group_by(code_level2 ) %>% dplyr::sample_frac(fracs[i]) %>% as.data.frame()
      } else if(data_set == "forced_imbalance") {
        a <- TD_df %>% filter(code_level2 ==2)  %>% dplyr::sample_n(forced_imblace_irri[i],replace = TRUE) %>% as.data.frame()
        b <- TD_df %>% filter(code_level2 !=2) %>% group_by(code_level2 )  %>% dplyr::sample_n(forced_imblace_rest[i], replace = TRUE) %>% as.data.frame()
        trainDat <- rbind(a,b) 
      }else if(data_set == "misclassified") {
        replace_vals <- sample(class_vals,mislabelled_num[i],replace = TRUE)
        trainDat <- TD_df 
        trainDat$code_level2 <- replace(trainDat$code_level2, sample(nrow(trainDat),mislabelled_num[i]), replace_vals)
      } 
      
      
      model_norm <- model_train(trainDat,"6m" , algo)
      
      if(data_set == "balanced_vals") {
        data_set0 = balanced_vals
      } else if(data_set == "fracs") {
        data_set0 = fracs
      } else if(data_set == "forced_imbalance") {
        data_set0 = forced_imblace_irri
      }else if(data_set == "misclassified") {
        data_set0 = mislabelled_num}
      
      testing_data <- testing(polys_split) %>% select(-code_level2, - PolygonID)
      testing_code <- testing(polys_split) %>% select(code_level2) %>% mutate(code_level2 = as.factor(code_level2))
      predictions <- predict(model_norm, testing_data) %>% as.data.frame()%>% mutate(prediction = as.factor(.)) %>%
        select(prediction)
      levels(predictions$prediction) <- levels(testing_code$code_level2)
      
      cm <- confusionMatrix((testing_code$code_level2), (predictions$prediction))
      
      
      if(algo == "rf") {
        cvPredictions <- model_norm$pred[model_norm$pred$mtry==model_norm$bestTune$mtry,]
      } else if(algo == "svmRadial"){
        cvPredictions <- model_norm$pred[model_norm$pred$C==model_norm$bestTune$C,]
      } else if(algo == "nnet"){
        cvPredictions <- model_norm$pred[model_norm$pred$size==model_norm$bestTune$size &
                                           model_norm$pred$decay==model_norm$bestTune$decay  ,]
      } else if(algo == "knn"){
        cvPredictions <- model_norm$pred[model_norm$pred$k==model_norm$bestTune$k,]
      }
      cm_model <- confusionMatrix(cvPredictions$pred,cvPredictions$obs)
      
      # overall accuracy
      overall_accu <- cm$overall['Accuracy'] %>% round(2)
      overall_accu_model <- cm_model$overall['Accuracy'] %>% round(2)
      
      accu_vect <- rbind(accu_vect, 
                         data.frame(overall_accuracy = overall_accu[[1]],
                                    overall_accuracy_model = overall_accu_model[[1]], 
                                    fraction = data_set0[i],
                                    pixels = nrow(trainDat),
                                    seed =  j,
                                    algo = algo,
                                    best_model = model_norm$bestTune,
                                    data_set=data_set))     
      
      #irrigated accuracy on tested 
      UP_accuracies <- data.frame(UA = round(cm$byClass[,3],3), PA = round(cm$byClass[,1],3))
      UP_accuracies_model <- data.frame(UA_model = round(cm_model$byClass[,3],3), PA_model = round(cm_model$byClass[,1],3))
      
      UP_accuracies <- rownames_to_column(UP_accuracies, "Class")
      UP_accuracies_model <- rownames_to_column(UP_accuracies_model, "Class")
      
      irri_accu <- UP_accuracies %>% filter(Class == "Class: 2")
      irri_accu_model <- UP_accuracies_model %>% filter(Class == "Class: 2")
      
      accu_vect_irri <- rbind(accu_vect_irri, 
                              data.frame(irri_UA = irri_accu[2], 
                                         irri_PA = irri_accu[3],
                                         irri_UA_model = irri_accu_model[2], 
                                         irri_PA_model = irri_accu_model[3],
                                         fraction = data_set0[i],
                                         pixels_irri = nrow(trainDat %>% filter(code_level2 == 2)), 
                                         seed =  j,
                                         algo = algo,
                                         best_model = model_norm$bestTune,
                                         data_set=data_set))
      
      if(data_set == "balanced_vals") {
        out <-   paste("model", i, "of", length(balanced_vals), 
                       "- seed", j, 'of', length(seed_vals),
                       "- Pixels per class:", balanced_vals[i], 
                       "- Overall accuracy:",overall_accu,
                       "- Algo:", algo)
      } else if(data_set == "fracs") {
        out <-   paste("model", i, "of", length(fracs), 
                       "- seed", j, 'of', length(seed_vals),
                       "- Fraction:", fracs[i], 
                       "- Overall accuracy:",overall_accu)
      } else if(data_set == "forced_imbalance") {
        out <-   paste("model", i, "of", length(forced_imblace_irri), 
                       "- seed", j, 'of', length(seed_vals),
                       "- Irri pixels:", forced_imblace_irri[i], 
                       "- Overall accuracy:",overall_accu)
      } else if(data_set == "misclassified") {
        out <-   paste("model", i, "of", length(mislabelled_num), 
                       "- seed", j, 'of', length(seed_vals),
                       "- Pixels mislabelled:", mislabelled_num[i], 
                       "- Overall accuracy:",overall_accu)
      }
      
      
      print(out) 
      
    }
    
  }
  
  a <- list(accu_vect, accu_vect_irri)
  
  saveRDS(accu_vect, here("output", "TD", "overall_accus", paste0(location, "_accu_",data_set,  algo,".rds")))
  saveRDS(accu_vect_irri, here("output", "TD", "overall_accus", paste0(location, "_accu_",data_set,  algo, "_IRRI.rds")))
  
  return(a)
}

model_a <- repeat_classification_function("forced_imbalance", "svmRadial", "Gaza_province") #data_set = fracs or balanced_vals forced_imbalance misclassified
model_b <- repeat_classification_function("forced_imbalance", "rf", "Gaza_province") #provinces = Gaza_province or Manica_province
model_c <- repeat_classification_function("forced_imbalance", "nnet", "Gaza_province")

