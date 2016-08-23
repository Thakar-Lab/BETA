source(Bayes_functions.R)
P <- 0.13#No.of TF targets in the background set/ No. of genes in the background set
N <- 188#Number of differentially expressed (DE) genes
X <- 59#Number of TF targets in the DE genes

bayes = calculate_bayes(x=X,N=N,p=P,max_sigma=max_sigma,length_sigma=length_sigma)    
plot(bayes,type="l")
computeSigmaP(bayes)#gives p-value
calcBayesOutputInfo(bayes)##outputs mean and CI interval
