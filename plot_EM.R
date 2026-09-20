prior_history <- read_csv("/project2/mstephens/wdenault/susie_mix/results_em/prior_history.csv")
prior_history
View(prior_history)

rec= prior_history[which(prior_history$iteration== max(prior_history$iteration)) ,]
rec[order(rec$pi_rec, decreasing = TRUE),]

rec[order(rec$pi_add ),]

rec$ncs=rec$pip_add+rec$pip_rec+rec$pip_dom
hist(, nclas=100)

plot(rec$pi_rec, rec$ncs)

plot(rec$pi_dom, rec$ncs)
library(ggplot2)

tt= prior_history[ which(prior_history$iteration==1),]

tt$iteration=0
tt$pi_add=0.3
tt$pi_rec=0.3
tt$pi_dom=0.3
prior_hist=rbind(tt,prior_history  )




ggplot(prior_hist, aes(x = iteration, color = tissue)) +
  geom_line(aes(y = pi_add)) +
  geom_point(aes(y = pi_add)) +
  geom_line(aes(y = pi_rec)) +
  geom_point(aes(y = pi_rec)) +
  geom_line(aes(y = pi_dom)) +
  geom_point(aes(y = pi_dom)) +
  ylab("Estimated prior inclusion probability") +
  ylim(0, 1)+
  geom_hline(yintercept = .8)+
  geom_hline(yintercept = .85)+

  geom_hline(yintercept = .9)


