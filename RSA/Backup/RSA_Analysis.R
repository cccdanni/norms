library("yarrr")
library("afex")

folder = "C:\\Users\\Psychology\\Desktop\\Research\\social-norms-learning\\Exp1_Norms_Learning\\rawData\\RSA\\data"
setwd(folder)

library(readxl)
pre_alltrials_full_trans <- read_excel("pre_alltrials_full_trans.xlsx")
post_alltrials_full_trans <- read_excel("post_alltrials_full_trans.xlsx")

pre_alltrials_full_trans$Compare[pre_alltrials_full_trans$Compare==4] = 3
post_alltrials_full_trans$Compare[post_alltrials_full_trans$Compare==4] = 3

pre_alltrials_full_trans$part = "pre"
post_alltrials_full_trans$part = "post"

pre_alltrials_full_trans$Compare[pre_alltrials_full_trans$Compare==1] = "higher"
post_alltrials_full_trans$Compare[post_alltrials_full_trans$Compare==1] = "higher"
pre_alltrials_full_trans$Compare[pre_alltrials_full_trans$Compare==2] = "lower"
post_alltrials_full_trans$Compare[post_alltrials_full_trans$Compare==2] = "lower"
pre_alltrials_full_trans$Compare[pre_alltrials_full_trans$Compare==3] = "consistent"
post_alltrials_full_trans$Compare[post_alltrials_full_trans$Compare==3] = "consistent"

pre_alltrials_full_trans$Group[pre_alltrials_full_trans$Group==1] = "ingroup"
post_alltrials_full_trans$Group[post_alltrials_full_trans$Group==1] = "ingroup"
pre_alltrials_full_trans$Group[pre_alltrials_full_trans$Group==2] = "outgroup"
post_alltrials_full_trans$Group[post_alltrials_full_trans$Group==2] = "outgroup"

alltrials_full_trans = rbind(pre_alltrials_full_trans, post_alltrials_full_trans)

# Neural disimilarity -----------
# 2 by 2 by 3
afex::aov_ez(id = "SubjectID", dv = "meanNeuralRDM2Morph", data = alltrials_full_trans,  within = c("Group", "Compare", "part"))
yarrr::pirateplot(meanNeuralRDM2Morph ~ Group+Compare+part, data = alltrials_full_trans)

# 2 by 3
post_alltrials_full_trans$meanNeuralRDM2Morph_diff = post_alltrials_full_trans$meanNeuralRDM2Morph - pre_alltrials_full_trans$meanNeuralRDM2Morph
afex::aov_ez(id = "SubjectID", dv = "meanNeuralRDM2Morph_diff", data = post_alltrials_full_trans,  within = c("Group", "Compare"))
yarrr::pirateplot(meanNeuralRDM2Morph ~ Group+Compare, data = alltrials_full_trans)

# post: 2 by 2 
this.data = subset(alltrials_full_trans, Compare!="consistent")
afex::aov_ez(id = "SubjectID", dv = "meanNeuralRDM2Morph", data = this.data,  
             within = c("Compare", "part","Group"))
pirateplot(meanNeuralRDM2Morph~Compare+part+Group, 
           data = this.data)
afex::aov_ez(id = "SubjectID", dv = "meanNeuralRDM2Morph", data = this.data[this.data$Group=="ingroup",],  
             within = c("Compare", "part"))
afex::aov_ez(id = "SubjectID", dv = "meanNeuralRDM2Morph", data = this.data[this.data$Group=="outgroup",],  
             within = c("Compare", "part"))

this.data = subset(alltrials_full_trans, Compare!="consistent" & Group=="ingroup")
this.data$part = factor(this.data$part, levels = c("pre","post"))
this.data.sum = summarySE(data = this.data, measurevar = "meanNeuralRDM2Morph", groupvars = c("Compare", "part"), na.rm = TRUE)
this.data$subjectID = factor(this.data$SubjectID)
ggplot(data = this.data.sum, aes(x=part, y=meanNeuralRDM2Morph, group=Compare, fill=Compare))+
  geom_bar(stat = "identity",position="dodge",width = 0.7)+
  geom_point(data = this.data, aes(x=part, y=meanNeuralRDM2Morph, group=Compare),size =2,position=position_jitterdodge(dodge.width = .7, jitter.width = .2),alpha=.2)+
  geom_errorbar(aes(ymin=meanNeuralRDM2Morph-ci, ymax=meanNeuralRDM2Morph+ci, color = Compare),size=1.2, width=0.1, 
               stat="identity",position=position_dodge(.7))+
  geom_errorbar(aes(y=meanNeuralRDM2Morph, ymax=meanNeuralRDM2Morph, ymin=meanNeuralRDM2Morph, color = Compare), width = 0.7, size = 1.1,position=position_dodge(.7))+
  theme_bw()+
  labs(x="",y="Neural Disimilarity")+
  scale_fill_manual(values=c(in_color,color1_in),labels=c("higher","lower"))+
  scale_color_manual(values=c(in_color,color1_in),labels=c("higher","lower"))+
  mytheme+
  geom_blank(aes(y=1.5))+
  geom_blank(aes(y=1.5))+
  scale_y_continuous(breaks=seq(0, 1.6, 0.2))+
  ggsignif::geom_signif(y_position = 1.08, xmin = 1, xmax = 2, annotation = "+", tip_length = 0.0)+
  theme(legend.title = element_blank())+
  # geom_hline(aes(yintercept = 0), size = 1, linetype = "dashed")+
  ylim (0.95, 1.1)

this.data = subset(alltrials_full_trans, Compare!="consistent" & Group=="outgroup")
this.data$part = factor(this.data$part, levels = c("pre","post"))
this.data.sum = summarySE(data = this.data, measurevar = "meanNeuralRDM2Morph", groupvars = c("Compare", "part"), na.rm = TRUE)
this.data$subjectID = factor(this.data$SubjectID)
ggplot(data = this.data.sum, aes(x=part, y=meanNeuralRDM2Morph, group=Compare, fill=Compare))+
  geom_bar(stat = "identity",position="dodge",width = 0.7)+
  geom_point(data = this.data, aes(x=part, y=meanNeuralRDM2Morph, group=Compare),size =2,position=position_jitterdodge(dodge.width = .7, jitter.width = .2),alpha=.2)+
  geom_errorbar(aes(ymin=meanNeuralRDM2Morph-ci, ymax=meanNeuralRDM2Morph+ci, color = Compare),size=1.2, width=0.1, 
                stat="identity",position=position_dodge(.7))+
  geom_errorbar(aes(y=meanNeuralRDM2Morph, ymax=meanNeuralRDM2Morph, ymin=meanNeuralRDM2Morph, color = Compare), width = 0.7, size = 1.1,position=position_dodge(.7))+
  theme_bw()+
  labs(x="",y="Neural Disimilarity")+
  scale_fill_manual(values=c(out_color,color2_put),labels=c("higher","lower"))+
  scale_color_manual(values=c(out_color,color2_put),labels=c("higher","lower"))+
  mytheme+
  geom_blank(aes(y=1.5))+
  geom_blank(aes(y=1.5))+
  scale_y_continuous(breaks=seq(0, 1.6, 0.2))+
  ggsignif::geom_signif(y_position = 1.08, xmin = 1, xmax = 2, annotation = "n.s.", tip_length = 0.0)+
  theme(legend.title = element_blank())+
  # geom_hline(aes(yintercept = 0), size = 1, linetype = "dashed")+
  ylim (0.95, 1.1)



# Beha disimilarity --------------
# post: 2 by 2 
this.data = subset(alltrials_full_trans, Compare!="consistent")
afex::aov_ez(id = "SubjectID", dv = "meanBehaRDM2Morph", data = this.data,  
             within = c("Compare", "part","Group"))
pirateplot(meanBehaRDM2Morph~Compare+part+Group, 
           data = this.data)
afex::aov_ez(id = "SubjectID", dv = "meanBehaRDM2Morph", data = this.data[this.data$Group=="ingroup",],  
             within = c("Compare", "part"))
afex::aov_ez(id = "SubjectID", dv = "meanBehaRDM2Morph", data = this.data[this.data$Group=="outgroup",],  
             within = c("Compare", "part"))

this.data = subset(alltrials_full_trans, Compare!="consistent" & Group=="ingroup")
this.data$part = factor(this.data$part, levels = c("pre","post"))
this.data.sum = summarySE(data = this.data, measurevar = "meanBehaRDM2Morph", groupvars = c("Compare", "part"), na.rm = TRUE)
this.data$subjectID = factor(this.data$SubjectID)
ggplot(data = this.data.sum, aes(x=part, y=meanBehaRDM2Morph, group=Compare, fill=Compare))+
  geom_bar(stat = "identity",position="dodge",width = 0.7)+
  geom_point(data = this.data, aes(x=part, y=meanBehaRDM2Morph, group=Compare),size =2,position=position_jitterdodge(dodge.width = .7, jitter.width = .2),alpha=.2)+
  geom_errorbar(aes(ymin=meanBehaRDM2Morph-ci, ymax=meanBehaRDM2Morph+ci, color = Compare),size=1.2, width=0.1, 
                stat="identity",position=position_dodge(.7))+
  geom_errorbar(aes(y=meanBehaRDM2Morph, ymax=meanBehaRDM2Morph, ymin=meanBehaRDM2Morph, color = Compare), width = 0.7, size = 1.1,position=position_dodge(.7))+
  theme_bw()+
  labs(x="",y="Attractiveness Rating Disimilarity")+
  scale_fill_manual(values=c(in_color,color1_in),labels=c("higher","lower"))+
  scale_color_manual(values=c(in_color,color1_in),labels=c("higher","lower"))+
  mytheme+
  geom_blank(aes(y=1.5))+
  geom_blank(aes(y=1.5))+
  scale_y_continuous(breaks=seq(1, 6.5, 0.5))+
  ggsignif::geom_signif(y_position = 6, xmin = 1, xmax = 2, annotation = "***", tip_length = 0.0)+
  theme(legend.title = element_blank())+
  # geom_hline(aes(yintercept = 0), size = 1, linetype = "dashed")+
  ylim (1, 6.5)

this.data = subset(alltrials_full_trans, Compare!="consistent" & Group=="outgroup")
this.data$part = factor(this.data$part, levels = c("pre","post"))
this.data.sum = summarySE(data = this.data, measurevar = "meanBehaRDM2Morph", groupvars = c("Compare", "part"), na.rm = TRUE)
this.data$subjectID = factor(this.data$SubjectID)
ggplot(data = this.data.sum, aes(x=part, y=meanBehaRDM2Morph, group=Compare, fill=Compare))+
  geom_bar(stat = "identity",position="dodge",width = 0.7)+
  geom_point(data = this.data, aes(x=part, y=meanBehaRDM2Morph, group=Compare),size =2,position=position_jitterdodge(dodge.width = .7, jitter.width = .2),alpha=.2)+
  geom_errorbar(aes(ymin=meanBehaRDM2Morph-ci, ymax=meanBehaRDM2Morph+ci, color = Compare),size=1.2, width=0.1, 
                stat="identity",position=position_dodge(.7))+
  geom_errorbar(aes(y=meanBehaRDM2Morph, ymax=meanBehaRDM2Morph, ymin=meanBehaRDM2Morph, color = Compare), width = 0.7, size = 1.1,position=position_dodge(.7))+
  theme_bw()+
  labs(x="",y="Attractiveness Rating Disimilarity")+
  scale_fill_manual(values=c(out_color,color2_put),labels=c("higher","lower"))+
  scale_color_manual(values=c(out_color,color2_put),labels=c("higher","lower"))+
  mytheme+
  geom_blank(aes(y=1.5))+
  geom_blank(aes(y=1.5))+
  scale_y_continuous(breaks=seq(1, 6.5, 0.5))+
  ggsignif::geom_signif(y_position = 6, xmin = 1, xmax = 2, annotation = "*", tip_length = 0.0)+
  theme(legend.title = element_blank())+
  # geom_hline(aes(yintercept = 0), size = 1, linetype = "dashed")+
  ylim (1, 6.5)


# Correlation ---------------
this.data = subset(alltrials_full_trans, Compare!="consistent")
afex::aov_ez(id = "SubjectID", dv = "NeuralBehaCorr", data = this.data,  
             within = c("Compare", "part","Group"))
pirateplot(NeuralBehaCorr~Compare+part+Group, 
           data = this.data)
afex::aov_ez(id = "SubjectID", dv = "NeuralBehaCorr", data = this.data[this.data$Group=="ingroup",],  
             within = c("Compare", "part"))
afex::aov_ez(id = "SubjectID", dv = "NeuralBehaCorr", data = this.data[this.data$Group=="outgroup",],  
             within = c("Compare", "part"))

this.data = subset(alltrials_full_trans, Compare!="consistent" & Group=="ingroup")
this.data$part = factor(this.data$part, levels = c("pre","post"))
this.data.sum = summarySE(data = this.data, measurevar = "NeuralBehaCorr", groupvars = c("Compare", "part"), na.rm = TRUE)
this.data$subjectID = factor(this.data$SubjectID)
ggplot(data = this.data.sum, aes(x=part, y=NeuralBehaCorr, group=Compare, fill=Compare))+
  geom_bar(stat = "identity",position="dodge",width = 0.7)+
  geom_point(data = this.data, aes(x=part, y=NeuralBehaCorr, group=Compare),size =2,position=position_jitterdodge(dodge.width = .7, jitter.width = .2),alpha=.2)+
  geom_errorbar(aes(ymin=NeuralBehaCorr-ci, ymax=NeuralBehaCorr+ci, color = Compare),size=1.2, width=0.1, 
                stat="identity",position=position_dodge(.7))+
  geom_errorbar(aes(y=NeuralBehaCorr, ymax=NeuralBehaCorr, ymin=NeuralBehaCorr, color = Compare), width = 0.7, size = 1.1,position=position_dodge(.7))+
  theme_bw()+
  labs(x="",y="r (Neural-Attractiveness rating)")+
  scale_fill_manual(values=c(in_color,color1_in),labels=c("higher","lower"))+
  scale_color_manual(values=c(in_color,color1_in),labels=c("higher","lower"))+
  mytheme+
  geom_blank(aes(y=1.5))+
  geom_blank(aes(y=1.5))+
  scale_y_continuous(breaks=seq(-1, 6.5, 0.5))+
  theme(legend.title = element_blank())+
  geom_hline(aes(yintercept = 0), size = 1, linetype = "dashed")+
  ylim (-.25, .25)

this.data = subset(alltrials_full_trans, Compare!="consistent" & Group=="outgroup")
this.data$part = factor(this.data$part, levels = c("pre","post"))
this.data.sum = summarySE(data = this.data, measurevar = "NeuralBehaCorr", groupvars = c("Compare", "part"), na.rm = TRUE)
this.data$subjectID = factor(this.data$SubjectID)
ggplot(data = this.data.sum, aes(x=part, y=NeuralBehaCorr, group=Compare, fill=Compare))+
  geom_bar(stat = "identity",position="dodge",width = 0.7)+
  geom_point(data = this.data, aes(x=part, y=NeuralBehaCorr, group=Compare),size =2,position=position_jitterdodge(dodge.width = .7, jitter.width = .2),alpha=.2)+
  geom_errorbar(aes(ymin=NeuralBehaCorr-ci, ymax=NeuralBehaCorr+ci, color = Compare),size=1.2, width=0.1, 
                stat="identity",position=position_dodge(.7))+
  geom_errorbar(aes(y=NeuralBehaCorr, ymax=NeuralBehaCorr, ymin=NeuralBehaCorr, color = Compare), width = 0.7, size = 1.1,position=position_dodge(.7))+
  theme_bw()+
  labs(x="",y="r (Neural-Attractiveness rating")+
  scale_fill_manual(values=c(out_color,color2_put),labels=c("higher","lower"))+
  scale_color_manual(values=c(out_color,color2_put),labels=c("higher","lower"))+
  mytheme+
  geom_blank(aes(y=1.5))+
  geom_blank(aes(y=1.5))+
  scale_y_continuous(breaks=seq(-1, 6.5, 0.5))+
  #ggsignif::geom_signif(y_position = 6, xmin = 1, xmax = 2, annotation = "*", tip_length = 0.0)+
  theme(legend.title = element_blank())+
  geom_hline(aes(yintercept = 0), size = 1, linetype = "dashed")+
  ylim (-0.25, 0.25)







t.test(meanNeuralRDM2Morph~Compare, data = alltrials_full_trans[alltrials_full_trans$Compare!="consistent" & alltrials_full_trans$Group=="ingroup" & alltrials_full_trans$part=="post",])
pirateplot(meanNeuralRDM2Morph~Compare, data = alltrials_full_trans[alltrials_full_trans$Compare!="consistent" & alltrials_full_trans$Group=="ingroup" & alltrials_full_trans$part=="post",])

afex::aov_ez(id = "SubjectID", dv = "meanNeuralRDM2Morph", data = alltrials_full_trans[alltrials_full_trans$Compare!=3 & alltrials_full_trans$Group==1,],  within = c("Compare", "part"))
yarrr::pirateplot(meanNeuralRDM2Morph ~ Compare+part, data = alltrials_full_trans[alltrials_full_trans$Compare!=3 & alltrials_full_trans$Group==1,])


afex::aov_ez(id = "SubjectID", dv = "meanBehaRDM2Morph", data = alltrials_full_trans,  within = c("Group", "Compare", "part"))
yarrr::pirateplot(meanBehaRDM2Morph ~ Group+Compare+part, data = alltrials_full_trans)

afex::aov_ez(id = "SubjectID", dv = "NeuralBehaCorr", data = alltrials_full_trans,  within = c("Group", "Compare", "part"))
yarrr::pirateplot(NeuralBehaCorr ~ Group+Compare+part, data = alltrials_full_trans)
t.test(NeuralBehaCorr~part, data = alltrials_full_trans[alltrials_full_trans$Compare!=3,], paired = TRUE)
