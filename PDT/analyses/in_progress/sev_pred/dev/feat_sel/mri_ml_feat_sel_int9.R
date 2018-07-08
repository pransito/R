# training data to fool around with
cur_dat = feature_clusters_cv[[kk]][[ii]][pred_ind]
cur_lab = feature_clusters_cv[[kk]][[ii]]$HCPG

# further preprocessing
correlationMatrix = cor(cur_dat)
highlyCorrelated  = findCorrelation(correlationMatrix, cutoff=0.90)
cur_dat           = cur_dat[-highlyCorrelated]

# parallel
cl    = makeCluster(detectCores()-1)
registerDoSNOW(cl)
# random forrest
control = rfeControl(functions=rfFuncs, method="cv", number=10,allowParallel = T,verbose=T)
# run the RFE algorithm
results <- rfe(cur_dat, cur_lab, sizes=c(1:length(cur_dat)), rfeControl=control)
# stop cluster
stopCluster(cl)
# summarize the results
print(results)
# list the chosen features
predictors(results)
# plot the results
plot(results, type=c("g", "o"))

ctrl <- rfeControl(functions = lmFuncs,
                   method = "repeatedcv",
                   repeats = 5,
                   verbose = TRUE,
                   allowParallel = F)

lmProfile <- rfe(cur_dat, cur_lab,
                 sizes = c(1:length(cur_dat)),
                 rfeControl = ctrl)

lmProfile
