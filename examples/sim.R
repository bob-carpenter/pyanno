library("MCMCpack")

rdisc <- function(phi) {
  K <- length(phi);
  u <- runif(1,0,1);
  sum <- 0;
  for (k in 1:K) {
    sum <- sum + phi[k];
    if (u < sum) 
      return(k);
  }
  return(K);
}

K <- 3;   # categories
J <- 20;  # annotators
I <- 1000; # items

theta <- array(NA,dim=c(J,K,K));
for (k in 1:K) {
  alpha <- rep(1,K);
  alpha[k] <- 3 * K;
  for (j in 1:J)
    theta[j,k,1:K] <- rdirichlet(1,alpha);
 }

pi <- rdirichlet(1,3*(1:K))

z <- rep(0,I);
for (i in 1:I)
  z[i] <- rdisc(pi);

y_full <- array(NA,c(I,J));
for (i in 1:I) 
  for (j in 1:J)
    y_full[i,j] <- rdisc(theta[j,z[i],1:K]);

missing <- 0.7;
keep <- array(rbinom(I*J, 1, 1 - missing), c(I,J));

N <- sum(keep);
ii <- rep(NA,N);
jj <- rep(NA,N);
y <- rep(NA,N);

n <- 1;
for (i in 1:I) {
  for (j in 1:J) {
    if (keep[i,j]) {
      ii[n] <- i;
      jj[n] <- j;
      y[n] <- y_full[i,j];
      n <- n + 1;
    }
  }
 }

write.table(cbind(ii,jj,y), file="sim.tsv", sep='\t', col.names=FALSE, row.names=FALSE);
