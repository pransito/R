install.packages('cvTools')
library(cvTools)
fit = lmrob(abstinenztage_vor_baseline ~ abstinenztage_baseline_bis_t2,data=cur_data)
cvFit(fit, data = cur_data, y = cur_data$abstinenztage_vor_baseline, cost = rtmspe, 
      K = 10, R = 1000, costArgs = list(trim = 0.1), seed = 1234)
