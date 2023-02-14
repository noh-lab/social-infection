library(car)
library(lme4)
library(parameters)
library(effectsize)
library(ggplot2)

manuscript_theme = ggplot2::theme_bw() + ggplot2::theme(axis.line.x = element_line(color="black", size = .7),axis.line.y = element_line(color="black", size = .7),text=element_text(size=11))


#### host fitness #### 
spore <- read.table("host_fitness.clean.20220711.txt", h=T)

spore$symbiont <- relevel(spore$symbiont, "bb859")
spore$host <- factor(spore$host, levels=c("qs17", "qs18", "qs4", "qs859", "qs395", "qs433"))
spore$MOI <- factor(spore$MOI, levels=c("control","0.05%_0.3","0.1%_0.6","0.5%_3","2.5%_15"))
spore$date <- as.factor(spore$date)
str(spore)


# overall patterns in the two variables measured during the experiment
ggplot(subset(spore, MOI!="control"), aes(x=MOI, y=flow_infected, fill=type)) + geom_boxplot()
ggplot(subset(spore, MOI!="control"), aes(x=symbiont, y=percent_spore, fill=type)) +geom_boxplot()


# what affects host fitness, including host type, symbiont
fitness.1 <- lmer(data=spore, percent_spore ~ flow_infected*symbiont*type + (1|date) + (1|host)) 
Anova(fitness.1) 

fitness.2 <- lmer(data=spore, percent_spore ~ flow_infected*symbiont + type + (1|date) + (1|host))
Anova(fitness.2) 

anova(fitness.1, fitness.2) # not significant; type is never significant

hist(residuals(fitness.2))
plot(fitted(fitness.2), residuals(fitness.2))    

#model_parameters(fitness.2, effects="fixed")
#model_parameters(anova(fitness.2))

epsilon_squared(fitness.2)


# visualize effect of symbiont, and lack of effect of host type
# full factorial design
ggplot(spore, aes(x=flow_infected, y=percent_spore, fill=symbiont, shape=type)) +geom_point(size=2) + scale_shape_manual(values=c(21, 22, 23)) + geom_smooth(method = "lm", se=FALSE) + facet_grid(host~symbiont) + ylim(60,110)

# grouped by type - first one is easiest to see
ggplot(spore, aes(x=flow_infected, y=percent_spore, fill=symbiont, shape=type)) +geom_point(size=1.5) + scale_shape_manual(values=c(21, 22, 23)) + geom_smooth(method = "lm", se=FALSE, aes(group=symbiont, color=symbiont)) + facet_grid(type~symbiont) + ylab("Host fitness") + xlab("Infection prevalence") + manuscript_theme

ggplot(spore, aes(x=flow_infected, y=percent_spore, fill=symbiont, shape=type)) +geom_point(size=2) + scale_shape_manual(values=c(21, 22, 23)) + geom_smooth(method = "lm", se=FALSE, aes(color=symbiont)) + facet_grid(type~.) + ylim(60,110) + manuscript_theme


#### symbiont transmission ####
horiz <- read.table("symbiont_transmission.clean.20220711.txt", h=T)

horiz$symbiont <- relevel(horiz$symbiont, "bb859")
horiz$host <- factor(horiz$host, levels=c("qs17", "qs18", "qs4", "qs859", "qs395", "qs433"))
horiz$MOI <- factor(horiz$MOI, levels=c("0.05%_0.3","0.1%_0.6","0.2%_1.2","0.5%_3","0.7%_4.2","1%_6","1.5%_9","2%_12","2.5%_15","5%_30"))
horiz$date <- as.factor(horiz$date)
str(horiz)


# overall patterns in the two variables measured during the experiment
ggplot(horiz, aes(x=symbiont, y=flow_infected, fill=type)) +geom_boxplot() 
ggplot(horiz, aes(x=symbiont, y=transmission, fill=type)) +geom_boxplot() 


# what affects transmission
transmit.1 <- lmer(data=horiz, transmission ~ flow_infected*symbiont*type + (1|date) + (1|host)) 
Anova(transmit.1) 

transmit.2 <- lmer(data=horiz, transmission ~ flow_infected*symbiont + flow_infected*type + symbiont*type + (1|date) + (1|host))
Anova(transmit.2, type=2) 

transmit.3 <- lmer(data=horiz, transmission ~ flow_infected*symbiont + type + (1|date) + (1|host))
Anova(transmit.3, type=2) 

anova(transmit.1, transmit.2) # not significant
anova(transmit.2, transmit.3) # significant

hist(residuals(transmit.2))
plot(fitted(transmit.2), residuals(transmit.2))    

#model_parameters(transmit.2, effects="fixed")
#model_parameters(anova(transmit.2))

epsilon_squared(transmit.2)


# visualize effect of symbiont, and lack of effect of host type
# full factorial design
ggplot(horiz, aes(x=flow_infected, y=transmission, fill=symbiont, shape=type)) +geom_point(size=2) + scale_shape_manual(values=c(21, 22, 23)) + geom_smooth(method = "lm", se=FALSE, aes(group=symbiont, color=symbiont)) + facet_grid(host~symbiont) 

# grouped by type - first one is easiest to see
ggplot(horiz, aes(x=flow_infected, y=transmission, fill=symbiont, shape=type)) +geom_point(size=1.5) + scale_shape_manual(values=c(21, 22, 23)) + geom_smooth(method = "lm", se=FALSE, aes(group=symbiont, color=symbiont)) + facet_grid(type~symbiont) + ylab("Horizontal transmission") + xlab("Infection prevalence") + manuscript_theme

ggplot(horiz, aes(x=flow_infected, y=transmission, fill=symbiont, shape=type)) +geom_point(size=2) + scale_shape_manual(values=c(21, 22, 23)) + geom_smooth(method = "lm", se=FALSE, aes(group=symbiont, color=symbiont)) + facet_grid(type~.) + manuscript_theme


#### analyses repeated with strict definition of native host ####
spore$strict <- spore$type
spore[row.names(subset(spore, host=="qs395" & symbiont!="bb395")),9] <- c("nonnative")
spore[row.names(subset(spore, host=="qs433" & symbiont!="bb433")),9] <- c("nonnative")
spore[row.names(subset(spore, host=="qs859" & symbiont!="bb859")),9] <- c("nonnative")

horiz$strict <- horiz$type
horiz[row.names(subset(horiz, host=="qs395" & symbiont!="bb395")),8] <- c("nonnative")
horiz[row.names(subset(horiz, host=="qs433" & symbiont!="bb433")),8] <- c("nonnative")
horiz[row.names(subset(horiz, host=="qs859" & symbiont!="bb859")),8] <- c("nonnative")

# overall patterns with the two variables measured during the experiment
ggplot(subset(spore, MOI!="control"), aes(x=MOI, y=flow_infected, fill=strict)) + geom_boxplot()
ggplot(subset(spore, MOI!="control"), aes(x=symbiont, y=percent_spore, fill=strict)) +geom_boxplot()

ggplot(horiz, aes(x=symbiont, y=flow_infected, fill=strict)) +geom_boxplot() 
ggplot(horiz, aes(x=symbiont, y=transmission, fill=strict)) +geom_boxplot() 


# what affects fitness
fitness.1s <- lmer(data=spore, percent_spore ~ flow_infected*symbiont*strict + (1|date) + (1|host)) 
Anova(fitness.1s) 

fitness.2s <- lmer(data=spore, percent_spore ~ flow_infected*symbiont + strict + (1|date) + (1|host))
Anova(fitness.2s) 
# similar result, strict type doesn't matter

anova(fitness.1s, fitness.2s) # not sig
epsilon_squared(fitness.2s)

ggplot(spore, aes(x=flow_infected, y=percent_spore, fill=symbiont, shape=strict)) +geom_point(size=1.5) + scale_shape_manual(values=c(21, 22, 23)) + geom_smooth(method = "lm", se=FALSE, aes(group=symbiont, color=symbiont)) + facet_grid(strict~symbiont) + ylab("Host fitness") + xlab("Infection prevalence") + manuscript_theme


# what affects transmission
transmit.1s <- lmer(data=horiz, transmission ~ flow_infected*symbiont*strict + (1|date) + (1|host)) 
Anova(transmit.1s) 

transmit.2s <- lmer(data=horiz, transmission ~ flow_infected*symbiont + flow_infected*strict + strict + (1|date) + (1|host))
Anova(transmit.2s, type=2) 

transmit.3s <- lmer(data=horiz, transmission ~ flow_infected*symbiont + strict + (1|date) + (1|host))
Anova(transmit.3s, type=2) 
# similar result, strict type doesn't matter

anova(transmit.1s, transmit.2s) # nonsig
anova(transmit.2s, transmit.3s) # nonsig

epsilon_squared(transmit.3s)

ggplot(horiz, aes(x=flow_infected, y=transmission, fill=symbiont, shape=strict)) +geom_point(size=1.5) + scale_shape_manual(values=c(21, 22, 23)) + geom_smooth(method = "lm", se=FALSE, aes(group=symbiont, color=symbiont)) + facet_grid(strict~symbiont) + ylab("Horizontal transmission") + xlab("Infection prevalence") + manuscript_theme


#### difference between two experimental rounds ####
library(ggbeeswarm)
library(gridExtra)

r1 <- ggplot(subset(spore, MOI!="control"), aes(x=symbiont, y=flow_infected, fill=symbiont, shape=type)) + geom_beeswarm(priority="density", cex=2) + scale_shape_manual(values=c(21, 22, 23)) + theme(legend.position="none") + ylab("Infection prevalence") + xlab("Symbiont") + ggtitle("(a) Host fitness exp.") + theme(plot.title = element_text(size=11)) + ylim(0,100) + manuscript_theme

r2 <- ggplot(horiz, aes(x=symbiont, y=flow_infected, fill=symbiont, shape=type)) + geom_beeswarm(priority = "density", cex=2) + scale_shape_manual(values=c(21, 22, 23)) + theme(legend.position="none") + ylab("") + xlab("Symbiont") + ggtitle("(b) Symbiont transmission exp.") + theme(plot.title = element_text(size=11)) + ylim(0,100) + manuscript_theme

grid.arrange(r1, r2, ncol=2)

r1.data <- subset(spore, MOI!="control")[,c(1,2,3,4,6,8)]
r2.data <- horiz[c(1,2,3,4,5,7)]
r1.data$stage <- c("one")
r2.data$stage <- c("two")
social <- rbind(r1.data, r2.data)
rm(r1.data, r2.data)

round <- lmer(flow_infected ~ stage*symbiont + (1|host), data=social)
Anova(round) # all significant
hist(residuals(round))
plot(fitted(round), residuals(round))    

epsilon_squared(round)


postscript(file=paste("fig3",format(Sys.time(),"%Y%m%d"),"eps",sep="."),onefile=F, width=6, height=4)
ggplot(spore, aes(x=flow_infected, y=percent_spore, fill=symbiont, shape=type)) +geom_point(size=1.5) + scale_shape_manual(values=c(21, 22, 23)) + geom_smooth(method = "lm", se=FALSE, aes(group=symbiont, color=symbiont)) + facet_grid(type~symbiont) + ylab("Host fitness") + xlab("Infection prevalence") + manuscript_theme
dev.off()


postscript(file=paste("fig4",format(Sys.time(),"%Y%m%d"),"eps",sep="."),onefile=F, width=6, height=4)
ggplot(horiz, aes(x=flow_infected, y=transmission, fill=symbiont, shape=type)) +geom_point(size=1.5) + scale_shape_manual(values=c(21, 22, 23)) + geom_smooth(method = "lm", se=FALSE, aes(group=symbiont, color=symbiont)) + facet_grid(type~symbiont) + ylab("Horizontal transmission") + xlab("Infection prevalence") + manuscript_theme
dev.off()

postscript(file=paste("fig5",format(Sys.time(),"%Y%m%d"),"eps",sep="."),onefile=F, width=8, height=3)
grid.arrange(r1, r2, ncol=2)
dev.off()


