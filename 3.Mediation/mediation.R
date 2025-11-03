### refer to code from Jiarui Chen
library(mediation)

####INput demon data
taxon <- read.csv("02_Taxa_Species.csv", row.names = 1)
metadata <- read.csv("01_Metadata.csv", row.names = 1)

taxon <- as.data.frame(t(taxon))
taxon <- taxon[ order(row.names(taxon)), ]
metadata <- metadata[ order(row.names(metadata)), ]

taxon <-taxon[, colSums(taxon > 0) >= 0.1*nrow(taxon)] #Filter species with prevalence lower than 10% ###modified

key_species <- "Eubacterium_eligens"

#####investigate whether the key species will attenuate the progression of Metabolic syndrome via the mediation of BMI
data_for_mediation <- as.data.frame(cbind(taxon$Eubacterium_eligens,metadata$Group,metadata$BMI))
rownames(data_for_mediation) <- rownames(taxon)
colnames(data_for_mediation) <- c("Key_Species", "Group", "BMI")
str(data_for_mediation)

data_for_mediation$Key_Species <- as.numeric(data_for_mediation$Key_Species)
data_for_mediation$BMI <- as.numeric(data_for_mediation$BMI)
data_for_mediation$Group <- factor(data_for_mediation$Group)
str(data_for_mediation)

####Key species是否会通过BMI来影响代谢综合征进展？
model_1 <- lm(BMI ~ Key_Species, data=data_for_mediation) ####自变量对中介变量的影响
summary(model_1)


model_2 <- glm(Group ~ Key_Species + BMI, data=data_for_mediation, family=binomial(link="probit")) ####自变量对中介变量的影响
summary(model_2)



model_mediation <- mediation::mediate(model_1, model_2, sims=1000, treat="Key_Species", mediator= "BMI")
summary(model_mediation) #中介效应存在: P-ACME < 0.05 (Average Causal Mediation Effect); Prop. Mediated平均中介比例：46.1%


####Key species是否会通过代谢综合征进展来影响BMI？
model_1 <- glm(Group ~ Key_Species, data=data_for_mediation, family=binomial(link="probit")) ####自变量对中介变量的影响
summary(model_1)
model_2 <- lm(BMI ~ Key_Species + Group, data=data_for_mediation) ####自变量对中介变量的影响
summary(model_2)
model_mediation <- mediation::mediate(model_1, model_2, sims=1000, treat="Key_Species", mediator= "BMI")
summary(model_mediation) #中介效应存在: P-ACME > 0.05 (Average Causal Mediation Effect); Prop. Mediated平均中介比例：0%
