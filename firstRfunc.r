# Baby's first R function
# I don't even remember what this was for...
mem2 <-
function(x,output) {
out <- lme(x~groups,random=~1|chips,na.action=na.exclude)
pval <- anova(out)$'p-value'[2]
#tval <- summary(model)$tTable[,'t-value'][2]
cat(pval,file=output,append=T,fill=T)}

apply(myData,1,mem2,output=‘output’)
