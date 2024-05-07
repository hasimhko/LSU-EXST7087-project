library(ranger) # random forest
library(kernlab) # svm
library(caret) # training pipeline
library(MLmetrics) # precision-recall curve and AUC
library(glmnet) # penalized Log reg
library(dplyr) # data wrangling
library(ggplot2) # plots
library(patchwork) # arranging multiple plots
library(ggsci) # color palette
library(broom) # to create a logistic regression model summary table

###############################################################################
############################### data preparation ##############################
###############################################################################

ml_df <- read.csv("Project/synth_data.csv", row.names = 1)

# cutoffs for categorical variables
ml_df <- ml_df %>%
  mutate(group_label=ifelse(group>0.55, 1, 0)) %>%
  mutate(sb_label=ifelse(sb>0.5, 1, 0)) %>%
  mutate(afb_label=ifelse(afb>0.55, 1, 0)) %>%
  select(-group)

ml_df$group_label <- as.factor(ifelse(ml_df$group_label==1, "VR", "VS"))
ml_df$sb_label <- as.factor(ifelse(ml_df$sb_label==1, "pos", "neg"))
ml_df$afb_label <- as.factor(ifelse(ml_df$afb_label==1, "pos", "neg"))

ml_df <- ml_df[, -c(2,3)]

###############################################################################
###### lists for storing model objects, accuracy and confusion matrices #######
###############################################################################

# lists for storing model objects per repetition
rf_models <- vector("list", 50)
svm_models <- vector("list", 50)
logreg_models <- vector("list", 50)

# lists for storing model test accuracy per repetition
rf_acc <- vector("list", 50)
svm_acc <- vector("list", 50)
logreg_acc <- vector("list", 50)

# lists for storing model objects per repetition
rf_confMat <- vector("list", 50)
svm_confMat <- vector("list", 50)
logreg_confMat <- vector("list", 50)

###############################################################################
############## model training, testing, and performance metrics ###############
###############################################################################

tunegrid_rf <- expand.grid(mtry = (3:6),
                           min.node.size = c(2,6,10),
                           splitrule = "gini") 
tunegrid_svm <- expand.grid(sigma = 2^(seq(-2, 2, by = 1)),
                            C = 2^(seq(-2, 2, by = 1))) 

for(i in 1:50){
  
  # important for reproducing CV splits and weights
  set.seed(i)
  ## for random forest cv splits and hyperparameters
  seeds_rf <- vector("list", length = 6)
  for(j in 1:5) seeds_rf[[j]] <- sample.int(n=1000, 12)
  seeds_rf[[6]] <- sample.int(1000, 1) # for the last RF model
  ## for support vector machine cv splits and hyperparameters
  seeds_svm <- vector("list", length = 6)
  for(j in 1:5) seeds_svm[[j]]<- sample.int(n=1000, 25)
  seeds_svm[[6]] <- sample.int(1000, 1) # for the last SVM model
  
  # train-test split
  ind <- createDataPartition(ml_df$group_label, p = 0.7, list = FALSE)
  train <- ml_df[ ind,] %>% mutate_if(is.numeric, scale)
  test  <- ml_df[-ind,] %>% mutate_if(is.numeric, scale)
  
  # training parameters of random forest
  tunecontrol_rf <- trainControl(method = "cv",
                                 number = 5,
                                 search = "grid",
                                 classProbs = FALSE,
                                 seeds = seeds_rf)
  
  # training parameters of support vector machine
  tunecontrol_svm <- trainControl(method = "cv",
                                  number = 5,
                                  search = "grid",
                                  classProbs = FALSE,
                                  seeds = seeds_svm)
  
  # training random forests
  rf <- train(group_label ~ .,
              data = train,
              method = "ranger",
              num.trees = 1000,
              importance = "impurity",
              trControl = tunecontrol_rf,
              tuneGrid = tunegrid_rf,
              metric = "Accuracy")
  
  # training support vectors
  svm <- train(group_label ~ .,
               data = train,
               method = "svmRadial",
               trControl = tunecontrol_svm,
               tuneGrid = tunegrid_svm,
               metric = "Accuracy")
  
  # training logistic regression
  logreg <- train(group_label ~ . + dwv:kv + dwv:vdv1 + kv:vdv1, 
                  data = train, 
                  method = "glm", 
                  family = "binomial",
                  metric = "Accuracy")
  
  # saving model instances for each repetition
  svm_models[[i]] <- svm
  rf_models[[i]] <- rf
  logreg_models[[i]] <- logreg
  
  # predicting and saving test set accuracy and confusion matrix for each repetition
  ## random forest
  rf_pred_temp <- predict(rf, test)
  rf_acc[[i]] <- postResample(rf_pred_temp, test$group_label)[1]
  rf_confMat[[i]] <- confusionMatrix(rf_pred_temp, test$group_label, mode = "prec_recall")
  ## support vecctor machine
  svm_pred_temp <- predict(svm, test)
  svm_acc[[i]] <- postResample(svm_pred_temp, test$group_label)[1]
  svm_confMat[[i]] <- confusionMatrix(svm_pred_temp, test$group, mode = "prec_recall")
  ## logistic regression
  logreg_pred_temp <- predict(logreg, test)
  logreg_acc[[i]] <- unname(postResample(logreg_pred_temp, test$group_label)[1])
  logreg_confMat[[i]] <- confusionMatrix(logreg_pred_temp, test$group_label, mode = "prec_recall")
  
}

###############################################################################
############## boxplots of model accuracy across all repetitions ##############
###############################################################################

acc_df <- data.frame(acc = c(unlist(rf_acc), unlist(svm_acc), unlist(logreg_acc)),
                     model = c(rep("Random Forest", 50), 
                               rep("Support Vector Machine", 50), 
                               rep("Logistic Regression", 50)))

acc_bp <- ggplot(acc_df, aes(x = model, y = acc, color = model)) +
  geom_boxplot() +
  geom_jitter() +
  xlab("") +
  ylab("Accuracy") +
  scale_color_npg() +
  theme_bw() +
  theme(legend.position="none")

###############################################################################
############# extracting best model and plotting confusion matrix #############
###############################################################################

best_rf_idx <- unname(which(unlist(rf_acc) == max(unlist(rf_acc))))[1]
best_svm_idx <- unname(which(unlist(svm_acc) == max(unlist(svm_acc))))[1]
best_logreg_idx <- unname(which(unlist(logreg_acc) == max(unlist(logreg_acc))))[1]

best_rf_cm <- rf_confMat[[best_rf_idx]]
best_svm_cm <- svm_confMat[[best_svm_idx]]
best_logreg_cm <- logreg_confMat[[best_logreg_idx]]

# combining data from all confusion matrices
confMat_all <- data.frame(pred_plot_label=factor(c("VR", "VR", "VS", "VS")),
                          ref_plot_label=factor(c("VR", "VS", "VR", "VS")),
                          rf=c(best_rf_cm[["table"]]),
                          svm=c(best_svm_cm[["table"]]),
                          logreg=c(best_logreg_cm[["table"]]))

# plotting the confusion matrix for the best-performing random forest
rf_cm_plot <- ggplot(confMat_all, mapping = aes(x = pred_plot_label, y = ref_plot_label)) +
  geom_tile(aes(fill = rf)) +
  geom_text(aes(label = sprintf("%1.0f", rf)), vjust = 1) +
  scale_fill_gradient(low="#E64B35B2", high="#4DBBD5B2") +
  xlab("Predicition") +
  ylab("Reference") +
  ggtitle("Best-performing random forest") +
  theme_bw() + 
  theme(legend.position="none", panel.grid=element_blank())

# plotting the confusion matrix for the best-performing support vector machine
svm_cm_plot <- ggplot(confMat_all, mapping = aes(x = pred_plot_label, y = ref_plot_label)) +
  geom_tile(aes(fill = svm)) +
  geom_text(aes(label = sprintf("%1.0f", svm)), vjust = 1) +
  scale_fill_gradient(low="#E64B35B2", high="#4DBBD5B2") +
  xlab("Predicition") +
  ylab("Reference") +
  ggtitle("Best-performing support vector machine") +
  theme_bw() + 
  theme(legend.position="none", panel.grid=element_blank())

# plotting the confusion matrix for the best-performing logistic regression
logreg_cm_plot <- ggplot(confMat_all, mapping = aes(x = pred_plot_label, y = ref_plot_label)) +
  geom_tile(aes(fill = logreg)) +
  geom_text(aes(label = sprintf("%1.0f", logreg)), vjust = 1) +
  scale_fill_gradient(low="#E64B35B2", high="#4DBBD5B2") +
  xlab("Predicition") +
  ylab("Reference") +
  ggtitle("Best-performing logistic regression") +
  theme_bw() + 
  theme(legend.position="none", panel.grid=element_blank())

# combining all plots in one and saving it
png("Project/model_test_results.png", width = 3000, height = 2000, res = 300)
(acc_bp + rf_cm_plot) / (svm_cm_plot + logreg_cm_plot) + plot_annotation(tag_levels = "A")
dev.off()

###############################################################################
############# extracting variable importance and coeff. estimates #############
###############################################################################

rf_varImp <- varImp(rf_models[[best_rf_idx]], scale = FALSE)[[1]]
rf_varImp$var <- c("Varroa", "Nosema", "ABPV", "BQCV", "DWV", "KV", "VDV1", "SB", "AFB")

svm_varImp <- varImp(svm_models[[best_svm_idx]], scale = FALSE)[[1]]
svm_varImp$var <- c("Varroa", "Nosema", "ABPV", "BQCV", "DWV", "KV", "VDV1", "SB", "AFB")

rf_varImp_plot <- rf_varImp %>%
  mutate(name = forcats::fct_reorder(var, Overall)) %>%
  ggplot(aes(x=name, y=Overall)) +
  geom_segment( aes(xend=name, yend=0)) +
  geom_point( size=4, color="#E64B35FF") +
  coord_flip() +
  ylab("Variable importance") +
  xlab("") +
  ggtitle("Best-performing random forest") +
  theme_bw()

svm_varImp_plot <- svm_varImp %>%
  mutate(name = forcats::fct_reorder(var, VS)) %>%
  ggplot(aes(x=name, y=VS)) +
  geom_segment( aes(xend=name, yend=0)) +
  geom_point( size=4, color="#E64B35FF") +
  coord_flip() +
  ylab("Variable importance") +
  xlab("") +
  ggtitle("Best-performing support vector machine") +
  theme_bw()

add_sig_sym <- function(pvalues){
  ifelse(pvalues < 0.001, "***", 
         ifelse(pvalues < 0.01, "**",
                ifelse(pvalues < 0.05, "*", NA)))
}

logreg_fit_term <- tidy(logreg_models[[sample(best_logreg_idx, 1)]][["finalModel"]], conf.int = TRUE) %>%
  filter(term !="(Intercept)") %>%
  mutate(p.stars = add_sig_sym(p.value)) %>%
  mutate(group = ifelse(estimate > 0, "pos", "neg"))

# write.csv(logreg_fit_term, "glm_fit_term.csv")

y_label <- rev(c("vdv1" = "VDV1", 
                 "varroa_count" = "Varroa", 
                 "sb_labelpos" = "SB[+]", 
                 "nosema" = "Nosema", 
                 "kv" = "KV",
                 "dwv" = "DWV", 
                 "bqcv" = "BQCV", 
                 "afb_labelpos" = "AFB[+]", 
                 "abpv" = "ABPV",
                 "`kv:vdv1`" = "KV \u00D7 VDV1",
                 "`dwv:vdv1`" = "DWV \u00D7 VDV1",
                 "`dwv:kv`" = "DWV \u00D7 KV"))

logreg_coeff_plot <- logreg_fit_term %>%
  mutate(name = forcats::fct_reorder(term, estimate)) %>%
  ggplot(aes(x=estimate, y=name, color=group, label=p.stars)) +
  geom_point(size=0.7) + ylab("") + xlab("Estimate") +
  geom_text(vjust=0, nudge_y=0.01, size=2.5, color=1) +
  geom_errorbar(aes(xmin=conf.low, xmax=conf.high, width=0.2)) +
  geom_vline(xintercept=0, linetype="dashed", color = "black", linewidth=0.3) + 
  scale_y_discrete(labels=y_label) +
  scale_x_continuous(breaks = c(-0.5, 0, 0.5, 1),
                     labels = c(paste0("\u2212","0.5"), "0", "0.5", "1")) +
  scale_color_npg(guide="none") + 
  ggtitle("Best-performing logistic regression") +
  theme_bw()

png("Project/model_variables.png", width = 3000, height = 2500, res = 300)
(rf_varImp_plot + svm_varImp_plot) / logreg_coeff_plot + plot_annotation(tag_levels = "A")
dev.off()