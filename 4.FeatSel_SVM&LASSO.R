setwd("ML_features_select/SVM&LASSO_features_filter")
library(tidyverse)
library(glmnet)
source("../../msvmRFE.R")
library(VennDiagram)
library(sigFeature)
library(e1071)
library(caret)
library(randomForest)
### 上下调模块基因
load("../../DEA&EA/up&down_regulated_module_DEGs.rda")
### 加载WGCNA输入矩阵（只是过滤了一个离群sample）
load("../../WGCNA/Step01-WGCNA_input.Rda")
### 临床信息
load("../../Data_preprocessing/disc_cohort.rda")
all_dif_module <- c(down_dif_express_module,up_dif_express_module)
diagnosis <- disc_cohort_clin[rownames(datExpr),]$diagnosis
datExpr <- as.data.frame(cbind(diagnosis,datExpr))
datExpr <- datExpr[,c('diagnosis',all_dif_module)]
str(datExpr)

# 设置第一列的列名为'group'
colnames(datExpr)[1] <- 'group'
datExpr$group <- as.factor(datExpr$group)
# 将除了第一列,表达数据设置为数值型
datExpr[, -1] <- lapply(datExpr[, -1], function(x) as.numeric(as.character(x)))
# 查看数据结构确认修改
str(datExpr)

#后面svmRFE函数要求group的类型为factor
train <- datExpr
train[1:4,1:4]

### 找特征
# 转为lasso需要的格式：
# train行是样本名，列是基因名
# 第一列是样本的分类；后面的列就是基因表达数据
x <- as.matrix(train[,-1]) %>%  scale()
y <- ifelse(train$group == "CON", 0,1) #把分组信息换成01

fit = glmnet(x, y, family = "binomial", alpha = 1, lambda = NULL)

options(repr.plot.width=5, repr.plot.height=5)

plot(fit, xvar = "dev", label = F)

pdf("A_lasso.pdf", width = 5, height = 5)
plot(fit, xvar = "dev", label = F)
dev.off()

set.seed(123)

cvfit = cv.glmnet(x, y, 
                  nfold=10, #例文描述：10-fold cross-validation；将队列分为训练集，验证集，并进行10次交叉验证
                  family = "binomial", type.measure = "class")
pdf(file='LASSO.Coefficient_Paths_Plot.pdf',height=5,width=5)
plot(cvfit)
dev.off()
plot(cvfit)
cvfit$lambda.min #查看最佳lambda

# 获取LASSO选出来的特征
myCoefs <- coef(cvfit, s="lambda.1se")#lambda.min或lambda.1se#lambda.min是最优；lambda.1se变量最少
lasso_fea <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
(lasso_fea <- lasso_fea[-1])
# 把lasso找到的特征保存到文件
library(openxlsx)
write.xlsx(as.data.frame(lasso_fea),"feature_lasso.xlsx")

predict <- predict(cvfit, newx = x[1:nrow(x),], s = "lambda.1se", type = "class")
table(predict,y)



################################
#############SVM-features-filtering
input <- train
input[,-1] <- scale(input[,-1])

#采用5折交叉验证 (k-fold crossValidation）
nfold = 5
nrows = nrow(input)
folds = rep(1:nfold, len=nrows)[sample(nrows)]
folds = lapply(1:nfold, function(x) which(folds == x))
results = lapply(folds, svmRFE.wrap, input, k=5, halve.above=100) #特征选择

library(openxlsx)
top.features = WriteFeatures(results, input, save=F) #查看主要变量
head(top.features)
#把SVM-REF找到的特征保存到文件
write.xlsx(top.features,"feature_svm.xlsx")

library(e1071)
library(kernlab)
library(caret)
library(parallel)
set.seed(123)

Profile=rfe(x=input[,-1], 
            y=as.numeric(as.factor(input[,1])),
            sizes = c(2, 4, 6, 8, seq(10, 449, by = 3)),
            rfeControl = rfeControl(functions = caretFuncs, method = "cv"),
            method="svmRadial")
save(Profile,file = "Profile.RData")

# 假设 x 和 y 是从 Profile$results 中提取的数据
x = Profile$results$Variables
y = Profile$results$RMSE

# 找到最低 RMSE 对应的特征数
wmin = which.min(y)
wmin.x = x[wmin]
wmin.y = y[wmin]
pdf(file='SVM最佳features折线图.pdf',height=5,width=5)
# 绘制折线图，显示 RMSE 随特征数的变化
plot(x, y, type="l", col="black", 
     xlab="Number of Features", ylab="Root Mean Squared Error (RMSE)", 
     main="RMSE vs Number of Features for SVM Model")
# 绘制最低 RMSE 点
points(wmin.x, wmin.y, col="black", pch=16)
text(wmin.x, wmin.y, labels=paste0('Optimal n=', wmin.x), pos=3, col="red")

# 添加一个水平线，突出显示最低 RMSE 的位置
abline(h=wmin.y, col="red", lty=2)
dev.off()
plot(x, y, type="l", col="black", 
     xlab="Number of Features", ylab="Root Mean Squared Error (RMSE)", 
     main="RMSE vs Number of Features for SVM Model")
# 绘制最低 RMSE 点
points(wmin.x, wmin.y, col="black", pch=16)
text(wmin.x, wmin.y, labels=paste0('Optimal n=', wmin.x), pos=3, col="red")

# 添加一个水平线，突出显示最低 RMSE 的位置
abline(h=wmin.y, col="red", lty=2)
save.image(file='All.RData')

saveRDS(lasso_fea,file='LASSO_features.rds')#保留LASSO筛选的基因集
saveRDS(top.features[1:379,"FeatureName"],file='SVM_features.rds')#保留SVM筛选的基因集
