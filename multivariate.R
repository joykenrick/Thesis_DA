### multivariate prediction model ###
library(tidymodels)

cancer_levels <- c('AML','CLL','LYMPH','MYEL','CRC','LUNGC','GLIOM','BRC', 'CVX','ENDC','OVC','PRC','HCC','PAN')
run_tidy_multiclass_lasso <- function(train_data,test_data){
  set.seed(222)
  # should I test with different values for step_corr?
  # Define recipe
  cancer_rec <-
    recipe(Disease  ~ ., data = train_data) %>% 
    update_role(DAid,new_role = 'ID') %>% 
    #added step_corr on 4/18
    step_corr(all_predictors()) %>%
    step_nzv(all_predictors()) %>% 
    #tested out bottom line on 4/14 to see if could adjust based on sex 
    #step_dummy(all_nominal_predictors()) %>% 
    #step_dummy('sex') --> maybe include for bivariate? 
    step_impute_knn(all_numeric())
  
  # Define model 
  # try tuning the mixture value as well? 
  logistic_mod <-
    multinom_reg(mixture = 1,penalty = tune()) %>% 
    set_mode('classification') %>% 
    set_engine('glmnet')
  
  # Define grid
  logistic_grid <- 
    grid_regular(extract_parameter_set_dials(logistic_mod),levels = 10)
  
  # Define tuning process
  model_control <- control_grid(save_pred = T)
  #model_metrics <- metric_set(accuracy, sens, spec, mn_log_loss, roc_auc)
  model_metrics <- metric_set(roc_auc)
  
  # Tune model
  logistic_res <- tune_grid(
    logistic_mod,
    cancer_rec,
    grid = logistic_grid,
    control = model_control,
    metrics = model_metrics,
    resamples = k_folds_data
  )
  
  final_param <- logistic_res %>% show_best("roc_auc") %>% slice(1) %>% select(penalty)
  
  'logistic_res %>% 
    collect_predictions() %>% 
    inner_join(final_param) %>% 
    conf_mat(truth = Disease, estimate = .pred_class) '
  
  #Finalize Model
  final_model <- workflow() %>% 
    add_model(logistic_mod) %>% 
    add_recipe(cancer_rec)
  
  #update logistic res model with the identified tuning parameter (penalty)
  final_model <- finalize_workflow(final_model, final_param)
  
  # estimate parameters for model given set of data
  final_fit <- fit(final_model,train_data)
  
  cancer_res <- augment(final_fit,test_data)
  
  variable_importance <-
    final_fit %>% 
    extract_fit_parsnip() %>% 
    vi()
  
  cancer_res <- augment(final_fit,test_data)
  
  roc_curve <- cancer_res %>% 
    roc_curve(Disease,.pred_AML:.pred_PRC) %>% autoplot()
  
  roc_res <- cancer_res %>% 
    roc_curve(Disease,.pred_AML:.pred_PRC)
  
  ## extract features per class
  coefficients_imp_features <- tidy(final_fit$fit$fit,exponentiate=T) %>%
    filter(term!='(Intercept)') %>% 
    split(.$class) %>% 
    map(~select(.x,term,estimate) %>% 
          arrange(-abs(estimate)))
  
  return(list(multiclass_res = cancer_res,
              variable_importance = coefficients_imp_features,
              roc_res = roc_res,
              roc_curve = roc_curve))
}

run_tidy_multiclass_mixture <- function(train_data,test_data){
  set.seed(222)
  
  # Define recipe
  cancer_rec <-
    recipe(Disease  ~ ., data = train_data) %>% 
    update_role(DAid,new_role = 'ID') %>% 
    step_corr(all_predictors()) %>%
    step_nzv(all_predictors()) %>%
    step_impute_knn(all_numeric())
  
  # Define model 
  # try tuning the mixture value as well? 
  logistic_mod <-
    multinom_reg(penalty = tune()) %>% 
    set_mode('classification') %>% 
    set_engine('glmnet')
  
  # Define grid
  logistic_grid <- 
    grid_regular(extract_parameter_set_dials(logistic_mod),levels = 10)
  
  # Define tuning process
  model_control <- control_grid(save_pred = T)
  #model_metrics <- metric_set(accuracy, sens, spec, mn_log_loss, roc_auc)
  model_metrics <- metric_set(roc_auc)
  
  # Tune model
  logistic_res <- tune_grid(
    logistic_mod,
    cancer_rec,
    grid = logistic_grid,
    control = model_control,
    metrics = model_metrics,
    resamples = k_folds_data
  )
  
  final_param <- logistic_res %>% show_best("roc_auc") %>% slice(1) %>% select(penalty)
  
  'logistic_res %>% 
    collect_predictions() %>% 
    inner_join(final_param) %>% 
    conf_mat(truth = Disease, estimate = .pred_class) '
  
  #Finalize Model
  final_model <- workflow() %>% 
    add_model(logistic_mod) %>% 
    add_recipe(cancer_rec)
  
  #update logistic res model with the identified tuning parameter (penalty)
  final_model <- finalize_workflow(final_model, final_param)
  
  # estimate parameters for model given set of data
  final_fit <- fit(final_model,train_data)
  
  cancer_res <- augment(final_fit,test_data)
  
  variable_importance <-
    final_fit %>% 
    extract_fit_parsnip() %>% 
    vi()
  
  cancer_res <- augment(final_fit,test_data)
  
  roc_curve <- cancer_res %>% 
    roc_curve(Disease,.pred_AML:.pred_PRC) %>% autoplot()
  
  roc_res <- cancer_res %>% 
    roc_curve(Disease,.pred_AML:.pred_PRC)
  
  ## extract features per class
  coefficients_imp_features <- tidy(final_fit$fit$fit,exponentiate=T) %>%
    filter(term!='(Intercept)') %>% 
    split(.$class) %>% 
    map(~select(.x,term,estimate) %>% 
          arrange(-abs(estimate)))
  
  return(list(multiclass_res = cancer_res,
              variable_importance = coefficients_imp_features,
              roc_res = roc_res,
              roc_curve = roc_curve))
}

run_tidy_multiclass_no_healthy <- function(train_data,test_data){
  set.seed(222)
  
  # Define recipe
  cancer_rec <-
    recipe(Disease  ~ ., data = train_data) %>% 
    step_corr(all_predictors()) %>%
    step_nzv(all_predictors()) %>% 
    update_role(DAid, new_role = 'ID') %>% 
    step_impute_knn(all_numeric())
  
  # Define model 
  logistic_mod <-
    multinom_reg(mixture = 1,penalty = tune()) %>% 
    set_mode('classification') %>% 
    set_engine('glmnet')
  
  # Define grid
  logistic_grid <- 
    grid_regular(extract_parameter_set_dials(logistic_mod),levels = 10)
  
  # Define tuning process
  model_control <- control_grid(save_pred = T)
  model_metrics <- metric_set(accuracy, sens, spec, mn_log_loss, roc_auc)
  
  # Tune model
  logistic_res <- tune_grid(
    logistic_mod,
    cancer_rec,
    grid = logistic_grid,
    control = model_control,
    metrics = model_metrics,
    resamples = k_folds_data
  )
  
  final_param <- logistic_res %>% show_best("roc_auc") %>% slice(1) %>% select(penalty)
  
  logistic_res %>% 
    collect_predictions() %>% 
    inner_join(final_param) %>% 
    conf_mat(truth = Disease, estimate = .pred_class) 
  
  #Finalize Model
  final_model <- workflow() %>% 
    add_model(logistic_mod) %>% 
    add_recipe(cancer_rec)
  
  #update logistic res model with the identified tuning parameter (penalty)
  final_model <- finalize_workflow(final_model, final_param)
  
  # estimate parameters for model given set of data
  final_fit <- fit(final_model,train_data)
  
  cancer_res <- augment(final_fit,test_data)
  
  variable_importance <-
    final_fit %>% 
    extract_fit_parsnip() %>% 
    vi()
  
  cancer_res <- augment(final_fit,test_data)
  
  roc_curve <- cancer_res %>% 
    roc_curve(Disease,.pred_AML:.pred_PRC) %>% autoplot()
  
  roc_res <- cancer_res %>% 
    roc_curve(Disease,.pred_AML:.pred_PRC)
  
  ## extract features per class
  coefficients_imp_features <- tidy(final_fit$fit$fit,exponentiate=T) %>%
    filter(term!='(Intercept)') %>% 
    split(.$class) %>% 
    map(~select(.x,term,estimate) %>% 
          arrange(-abs(estimate)))
  
  return(list(multiclass_res = cancer_res,
              variable_importance = coefficients_imp_features,
              roc_res = roc_res,
              roc_curve = roc_curve))
}


library(yardstick)

# Obtain predicted probabilities for each class
cancer_probs <- predict(final_fit, test_data, type = "prob")

# Pivot the data frame to long format
cancer_probs_long <- cancer_probs %>%
  pivot_longer(cols = starts_with(".pred_"), 
               names_to = "class", 
               values_to = "predicted_prob",
               names_pattern = "\\.pred_(.*)")

# Compute the AUC for each cancer
cancer_aucs <- cancer_probs_long %>%
  group_by(class) %>%
  summarize(auc = roc_auc(truth = test_data$Disease, 
                          .pred = predicted_prob))

cancer_aucs <- cancer_probs_long %>%
  mutate(predicted_prob = as.numeric(as.character(predicted_prob))) %>%
  group_by(class) %>%
  summarize(auc = roc_auc(truth = test_data$Disease, .pred = predicted_prob))

cancer_aucs <- cancer_probs_long %>%
  mutate(predicted_prob = as.numeric(as.character(predicted_prob))) %>%
  mutate(class = factor(class)) %>%
  group_by(class) %>%
  summarize(auc = roc_auc(truth = as.numeric(as.character(test_data$Disease)), 
                          .pred = as.numeric(predicted_prob)))
cancer_aucs <- cancer_probs_long %>%
  mutate(predicted_prob = as.numeric(as.character(predicted_prob))) %>%
  mutate(class = factor(class)) %>%
  group_by(class) %>%
  summarize(auc = roc_auc(truth = as.numeric(as.character(test_data$Disease)), 
                          .pred = as.numeric(predicted_prob)))
cancer_aucs <- cancer_probs_long %>%
  mutate(predicted_prob = as.numeric(as.character(predicted_prob))) %>%
  mutate(class = factor(class)) %>%
  group_by(class) %>%
  summarize(auc = roc_auc(truth = factor(test_data$Disease), 
                          .pred = as.numeric(predicted_prob)))


cancer_res %>%
  #filter(Resample == "Fold01") %>%
  select(Disease, .pred_AML:.pred_PRC) %>% 
  pivot_longer(-Disease, names_to = 'class', values_to = 'prob') %>% 
  group_nest(class) %>% 
  mutate(
    auc = map2_dbl(class, data, function(x, y) {
      obs <- factor(y$obs == x)
      est <- y$prob
      roc_auc_vec(obs, est, event_level = 'second')
    })
  )



# View the AUCs for each cancer
cancer_aucs



set.seed(222)

#convert cancer and healthy data to wide format
#filt_cancer_dat <- select(filt_cancer,DAid,Assay,NPX,Disease,Age,Sex) %>% 
  #filt_cancer_dat <- select(filt_IGT,DAid,Assay,NPX,Disease,Age) %>% 
  filt_cancer_dat <- select(filt_cancer,DAid,Assay,NPX,Disease) %>% 
    spread(Assay,NPX,-1)

#health_dat <- select(filt_IGT,DAid,Assay,NPX,Disease,Age,Sex) %>% 
 # health_dat <- select(filt_IGT,DAid,Assay,NPX,Disease,Age) %>% 
health_dat <- select(filt_IGT,DAid,Assay,NPX,Disease) %>% 
  spread(Assay,NPX,-1) %>% 
  drop_na() 

random_samps <- sample(health_dat$DAid,250,replace = F)
rand_healthy <-
  health_dat %>%
  filter(DAid %in% random_samps)

# combine healthy and cancer data
withIGT_combined_data <- bind_rows(filt_cancer_dat,rand_healthy)  
# convert Disease to factor
withIGT_combined_data <- mutate(combined_data, Disease = factor(Disease))

filt_cancer_dat <- mutate(filt_cancer_dat, Disease = factor(Disease))

withIGT_data_split <- initial_split(withIGT_combined_data, prop = .70,strata = Disease)
filt_cancer_split <- initial_split(filt_cancer_dat, prop = .70,strata = Disease)


# Create data frames for the two sets:
withIGT_train_data <- training(withIGT_data_split)
withIGT_test_data  <- testing(withIGT_data_split)

cancer_train_data <- training(filt_cancer_split)
cancer_test_data <- testing(filt_cancer_split)

withIGT_k_folds_data <- vfold_cv(withIGT_train_data,v = 5)
cancer_kfolds_data <- vfold_cv(cancer_train_data,v=5)
k_folds_data <- cancer_kfolds_data
## Run multiclass
multiclass_res <- march_version_run_tidy_multiclass(train_data,test_data)

IGT_multiclass_res <- run_tidy_multiclass(train_data,test_data)
cancer_multiclass_res <- run_tidy_multiclass_lasso(cancer_train_data,cancer_test_data)


mix_IGT_multiclass_res <- run_tidy_multiclass_mixture(train_data,test_data)
mix_cancer_multiclass_res <- run_tidy_multiclass_mixture(cancer_train_data,cancer_test_data)



#saveRDS(multiclass_res,'multiclass_res_1403.rds')

### generate top prots list based on limma and prediction model ### 
# for loop to extract top variable proteins for each cancer in prediction model res
for (cancer in cancer_levels){
  res <- IGT_multiclass_res$variable_importance[[as.character(cancer)]] %>% head(30)
  res <- filter(res,!estimate == 0)
  nam <- paste('IGT_multiclass_top_prots_', cancer, sep = "")
  assign(nam,res)
}

for (cancer in cancer_levels){
  res <- cancer_multiclass_res$variable_importance[[as.character(cancer)]] %>% head(150) %>% 
    mutate(Cancer = cancer) %>% 
    mutate(absval_estimate = abs(estimate)) %>% 
    mutate(scaled = rescale(absval_estimate,to=c(0,100)))
  res <- filter(res,!estimate == 0)
  nam <- paste('cancer_multiclass_top_prots_', cancer, sep = "")
  assign(nam,res)
}


cancer_multiclass_top_prots <- 
  rbind(cancer_multiclass_top_prots_AML,
        cancer_multiclass_top_prots_BRC,
        cancer_multiclass_top_prots_CLL,
        cancer_multiclass_top_prots_CRC,
        cancer_multiclass_top_prots_CVX,
        cancer_multiclass_top_prots_ENDC,
        cancer_multiclass_top_prots_GLIOM,
        cancer_multiclass_top_prots_HCC,
        cancer_multiclass_top_prots_LUNGC,
        cancer_multiclass_top_prots_LYMPH,
        cancer_multiclass_top_prots_MYEL,
        cancer_multiclass_top_prots_OVC,
        cancer_multiclass_top_prots_PAN,
        cancer_multiclass_top_prots_PRC)
cancer_multiclass_top_prots<- 
  mutate(cancer_multiclass_top_prots,absval_estimate = abs(cancer_multiclass_top_prots$estimate)) %>% 
  mutate(cancer_multiclass_top_prots,scaled = rescale(cancer_multiclass_top_prots$absval_estimate,to=c(0,100)))
