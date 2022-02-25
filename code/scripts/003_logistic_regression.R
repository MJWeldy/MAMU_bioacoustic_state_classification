logistic_accuracy <- list()
logistic_AUC <- list()
logistic_models <- list()

# Number of calls
df_1 <- data.frame(state = ifelse(site_states$response < 3, 0, 1), 
                   total_calls = rowSums(HMM_df[,2:ncol(HMM_df)], na.rm = TRUE))
logistic_models[[1]] <- glm(state ~ total_calls, data = df_1[!is.na(df_1$state),], family = "binomial")
model_glm_pred <- ifelse(predict(logistic_models[[1]], type = "link") > 0, 1, 0)
CM <- confusionMatrix(factor(df_1$state[!is.na(df_1$state)], levels = c(0,1)), 
                factor(model_glm_pred, levels = c(0,1)),
                mode = "prec_recall")
logistic_accuracy[[1]] <- CM$overall["Accuracy"] 
test_prob <- predict(logistic_models[[1]], newdata = df_1, type = "response")
logistic_AUC[[1]] = roc(df_1$state ~ test_prob, plot = FALSE, print.auc = TRUE)

# Additive effects of week and number of calls
x <- HMM_df
names(x) <- c("site","x13","x14","x15","x16","x17","x18","x19","x20","x21","x22","x23",  
              "x24","x25","x26","x27","x28","x29","x30","x31","x32","x33","x34","x35",  
              "x36","x37")
df_2 <- data.frame(state = ifelse(site_states$response < 3, 0, 1),x) %>% 
  pivot_longer(cols = starts_with("x"), names_to = "week", values_to = "num_calls")
df_2$week <- as.numeric(str_sub(df_2$week, 2)) 
logistic_models[[2]] <- glm(state ~ num_calls + week, data = df_2, family = "binomial")
model_glm_pred <- ifelse(predict(logistic_models[[2]], type = "link") > 0, 1, 0)
CM <- confusionMatrix(factor(df_2$state[!is.na(df_2$state)], levels = c(0,1)), 
                      factor(model_glm_pred[!is.na(df_2$state)], levels = c(0,1)),
                      mode = "prec_recall")
logistic_accuracy[[2]] <- CM$overall["Accuracy"] 
test_prob <- predict(logistic_models[[2]], newdata = df_2, type = "response")
logistic_AUC[[2]] = roc(df_2$state ~ test_prob, plot = FALSE, print.auc = TRUE)

# Interactive effect of week and number of calls

logistic_models[[3]] <- glm(state ~ num_calls * week, data = df_2, family = "binomial")
summary(logistic_models[[3]])

model_glm_pred <- ifelse(predict(logistic_models[[3]], type = "link") > 0, 1, 0)
CM <- confusionMatrix(factor(df_2$state[!is.na(df_2$state)], levels = c(0,1)), 
                      factor(model_glm_pred[!is.na(df_2$state)], levels = c(0,1)),
                      mode = "prec_recall")

logistic_accuracy[[3]] <- CM$overall["Accuracy"] 
test_prob <- predict(logistic_models[[3]], newdata = df_2, type = "response")
logistic_AUC[[3]] = roc(df_2$state ~ test_prob, plot = FALSE, print.auc = TRUE)