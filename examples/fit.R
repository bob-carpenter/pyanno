# K: number of cats
# J: number of annotators
# I: number of annotations
# N: number of labels
# ii[N] in 1:I:  item for label n
# jj[N] in 1:J:  annotator for label n
# y[N] in 1:K:   category for label n


### READ DATA
data <- read.csv("test.csv");
ii <- data[[1]];
jj <- data[[2]];
y <- data[[3]];

N <- length(ii);
K <- max(y);
J <- max(jj);
I <- max(ii);

print(paste("K=",K, "; N=",N, "; J=",J, ";I=",I));

### INITS
theta_hat <- array(NA,c(J,K,K));
for (j in 1:J)
  for (k in 1:K)
    for (k2 in 1:K)
      theta_hat[j,k,k2] <- ifelse(k==k2, 0.7, 0.3/K);

pi_hat <- array(1/K,K);

### EM ITERATIONS
epoch <- 1;
min_relative_diff <- 1E-8;
last_log_likelihood = - Inf;
E_z <- array(1/K, c(I,K));
for (epoch in 1:20) {
  ### E step 
  for (i in 1:I)
    E_z[i,] <- pi_hat;
  for (n in 1:N)
    for (k in 1:K)
      E_z[ii[n],k] <- E_z[ii[n],k] * theta_hat[jj[n],k,y[n]];
  for (i in 1:I)
    E_z[i,] <- E_z[i,] / sum(E_z[i,]);

  ### M step
  pi_hat <- rep(0,K);
  for (i in 1:I)
    pi_hat <- pi_hat + E_z[i,];
  pi_hat <- pi_hat / sum(pi_hat);

  count <- array(0,c(J,K,K));
  for (n in 1:N)
    for (k in 1:K)
      count[jj[n],k,y[n]] <- count[jj[n],k,y[n]] + E_z[ii[n],k];
  for (j in 1:J)
    for (k in 1:K)
      theta_hat[j,k,] <- count[j,k,] / sum(count[j,k,]);

  p <- array(0,c(I,K));
  for (i in 1:I)
    p[i,] <- pi_hat;
  for (n in 1:N)
    for (k in 1:K)
      p[ii[n],k] <- p[ii[n],k] * theta_hat[jj[n],k,y[n]];
  log_likelihood <- 0.0;
  for (i in 1:I)
    log_likelihood <- log_likelihood + log(sum(p[i,]));
  if (epoch == 1)
    print(paste("epoch=",epoch," log likelihood=", log_likelihood));
  if (epoch > 1) {
    diff <- log_likelihood - last_log_likelihood;
    relative_diff <- abs(diff / last_log_likelihood);
    print(paste("epoch=",epoch," log likelihood=", log_likelihood," relative_diff=",relative_diff));
    if (relative_diff < min_relative_diff) {
      print("FINISHED.");
      break;
    }
  }
  last_log_likelihood <- log_likelihood;
}


# estimands
# pi_hat[k]:  estimate of pi[k]
# theta_hat[j,k,k']:  estimate of theta[j,k,k']
# E_z[i,k]:   probabilistic estimate of z[i], Pr[z[i]==k | data] 

pi_out <- array(0,dim=c(K,2),row.names=FALSE,dimnames=list(NULL,c("category","prob")));
pos <- 1;
for (k in 1:K) {
  pi_out[pos,] <- c(k,pi_hat[k]);
  pos <- pos + 1;
}
write.table(pi_out,sep='\t',row.names=FALSE,file="pi_hat.tsv",quote=FALSE);

theta_out <- array(0,dim=c(J*K*K,4),dimnames=list(NULL,c("annotator","reference","response","prob")));
pos <- 1;
for (j in 1:J) {
  for (ref in 1:K) {
    for (resp in 1:K) {
      theta_out[pos,] <- c(j,ref,resp,theta_hat[j,ref,resp]);
      pos <- pos + 1;
    }
  }
}
write.table(theta_out,sep='\t',row.names=FALSE,file="theta_hat.tsv",quote=FALSE);

z_out <- array(0,dim=c(I*K,3),dimnames=list(NULL,c("item","category","prob")));
pos <- 1;
for (i in 1:I) {
  for (k in 1:K) {
    z_out[pos,] = c(i,k,E_z[i,k]);
    pos <- pos + 1;
  }
}
write.table(z_out,sep='\t',row.names=FALSE,file="z_hat.tsv",quote=FALSE);

