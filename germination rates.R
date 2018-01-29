germination<-read.csv("germination rates.csv")
head(germination)

#load plyr and summarize sample size and total germination
library(plyr)
sumgermination<- ddply(germination, "Alpha", summarise,
                 N = length(Alpha),
                 sum = sum(Total.growth..)
)

#check output
head(sumgermination)

write.csv(sumgermination, file="sumgermination.csv")
