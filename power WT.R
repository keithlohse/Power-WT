#By Keith Lohse PhD, Auburn University, 2015
#Descriptive statistics for these power calculations come from the NHANES dataset
#Heymsfield et al. (2014) Scaling of adult bodyweight to height across sex and race/ethnic groups: Relevant to BMI.
#American Journal of Clinical Nutrition, 100, 1455-1461.

########
#Weight#
########

#set parameters
sigma1 = 26.5 #approximate SD for adult, white, non-hispanic men in the US
sigma2 = 26.0 #approximate SD for adult, white, non-hispanic women in the US

mu1 = 89.4 #approximate mean WT in kg for men in the US
mu2 = 74.7 #approximate mean WT in kg for women in the US

n1 = 55 #size of sample 1
n2 = 55 #size of sample 2

pop1 <- rnorm(100000, mean=mu1, sd=sigma1) #set the mean, sigma, and number of datapoints in the parent distribution 
pop2 <- rnorm(100000, mean=mu2, sd=sigma2) #set the mean, sigma, and number of datapoints in the parent distribution


#set the length of your index to control the number of experiments you will run (default is 25)
index<-c(1:10000)

DATA<-data.frame(index) #We will create an empty dataframe to store our output in...
DATA$sp<-NULL
DATA$ES <- NULL
DATA$SES <- NULL

#This is the meat of the script, this for loop will create random samples (based on the index you specified above)
#These samples come from two different parent populations. These populations have the means and standard deviations that 
#you specified above.

for (i in 1:length(DATA$index)) {
  s1 <- sample(pop1, n1, replace = FALSE, prob = NULL)
  s2 <- sample(pop2, n2, replace = FALSE, prob = NULL)
  
  m1<-mean(s1)
  m2<-mean(s2)
  
  sd1<-sd(s1)
  sd2<-sd(s1)
  
  DATA$sp[i] = sqrt(((n1-1)*(sd1*sd1)+(n2-1)*(sd2*sd2))/(n1+n2-2))
  DATA$ES[i]<-(m1-m2)
  DATA$SES[i]<-(m1-m2)/DATA$sp[i]
  
}

mean(DATA$ES) #The mean ES should be approximatley 15(closer to 15 with more samples)
mean(DATA$SES) #The means SES should be approximately (mu1-mu2)/sigma

#Calculate the Standard Error in order to get MOEs
DATA$se<-sqrt((DATA$sp*DATA$sp/n1)+(DATA$sp*DATA$sp)/n2)

#Calculate the t_critical value for all samples
tcrit<-abs(qt(.025,(n1+n2-2)))
DATA$MOE<-tcrit*DATA$se #We can use the critical t-value to create the margin of error (MOE) for each test

#Converting the observed effect sizes into t-values
DATA$t_obs<-DATA$ES/DATA$se

#Calculating p-values for each t_obs
p_value<-2*(1-pt(DATA$t_obs,df=(n1+n2-2))) #p-value for a two-tailed test
DATA$p_value<-p_value #We can then add each of these p-values into our dataset

#We can also look at our p-value as a dichotomous yes/no decision
DATA$sig<-as.numeric(DATA$t_obs>abs(qt(0.025,(n1+n2-2))))
sum(DATA$sig)

tail(DATA) #We keep adding things to our dataframe, let's look at it again to make sure things are written properly...

#Let's visualize some of these p-values
plot(p_value, cex=2, lwd=2, font=2, font.lab=2, main='d~0.58, n/group=XX', ylab='Observed P-Value', ylim=c(0,1))
abline(h = 0.05, col = "red", lty = 3, lwd=3) #draws a red line at p = 0.05
abline(h = 0.1, col = "blue", lty = 3, lwd=3) #draws a blue line at p = 0.10 

#We can also get a sense of power by looking at the proportion of studies that find a significant effect
frac05<-sum(p_value>0.05)
frac05 #This is the number of tests that have p-values > 0.05 (essentially the number of Type 2 errors)

frac10<-sum(p_value>0.10)
frac10 #This is the number of tests that have p-values > 0.10

#Outputting the data frame as a text file
#Create a name for the file that makes sense to you
write.csv(DATA, file="power WT n55.csv")



