X0 = X_mix
X0 = X0[, -which(apply(X0, 2, var)<1e-5)]

pv= c()
y=  log1p(apply( Y_cor,1,sum) )
for ( k in 1: ncol(X0)){

  pv =c (pv, summary ( lm(y~X0[,k]))$coefficients[2,4])
}

plot( -log10(pv))
