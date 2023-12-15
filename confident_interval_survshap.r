library(survex)
library(ggplot2)
library(survival)
library(ranger) 
# read data
Final_patient_after_AI<-read.csv("AI_After_21genus.csv",header=TRUE)
Final_patient_before_AI<-read.csv("AI_Before_21genus.csv",header=TRUE)
#Final_patient_info_nonAI<-read.csv("Final_patient_info_cox0.0001_After_nonAI.csv",header=TRUE)
run_model = function(Final_patient_info_AI){
num_time_points <- 14
cd_auc_values_matrix <- matrix(NA, nrow = 100, ncol = num_time_points)
# round in 100 times

for (i in 1:100) {
    mymodel <- ranger(Surv(Time, Liver_metastasis) ~ ., data = Final_patient_info_AI)

    explainer <- explain(mymodel,
                         data = Final_patient_info_AI[,-c(1,2)],
                         y = Surv(Final_patient_info_AI$Time, Final_patient_info_AI$Liver_metastasis))
    

    perf <- model_performance(explainer)

    cd_auc_values_matrix[i, ] <- perf$result$`C/D AUC`
  }


mean_cd_auc <- apply(cd_auc_values_matrix, 2, mean, na.rm = TRUE)
std_cd_auc <- apply(cd_auc_values_matrix, 2, sd, na.rm = TRUE)


print("Mean C/D AUC at each time point:")
print(mean_cd_auc)
print("Standard Deviation of C/D AUC at each time point:")
print(std_cd_auc)

#create plot data frame

num_time_points <- length(mean_cd_auc)
num_actual_points <- length(mean_cd_auc)
time_points <- seq(20, 70, length.out = num_actual_points)
data <- data.frame(Time = time_points, Mean = mean_cd_auc, Std = std_cd_auc)

data$CI_lower <- data$Mean - 1.96 * data$Std
data$CI_upper <- data$Mean + 1.96 * data$Std
return(data)
}
data2 = run_model(Final_patient_before_AI)
data1 = run_model(Final_patient_after_AI)

gg1 = ggplot() +
  geom_line(data = data1, aes(x = Time, y = Mean, color = "Batch Effect Removed Data")) +
  geom_ribbon(data = data1, aes(x = Time, ymin = CI_lower, ymax = CI_upper), alpha = 0.2, fill = "purple") +
  geom_line(data = data2, aes(x = Time, y = Mean, color = "Original Data")) +
  geom_ribbon(data = data2, aes(x = Time, ymin = CI_lower, ymax = CI_upper), alpha = 0.2, fill = "orange") +
  scale_colour_manual(values = c("Batch Effect Removed Data" = "purple", "Original Data" = "orange")) +
  theme_minimal() +
  labs(title = "Mean C/D AUC Over Time (SurvSHAP(t))",
       x = "Time Point", y = "Mean C/D AUC",
       color = "Data Set")+theme(legend.text = element_text(size = 9.5))+theme(legend.position = "bottom")
print(gg1)

ggsave("survshap_auc_plot.pdf", plot = gg1, width = 8, height = 8, units = "in")