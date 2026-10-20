k=1

out [[k]]$fit_slide_prior$sets
 k=k+1


 k=k-1
 out [[k]]$fit_slide_prior$delta[,out [[k]]$fit_slide_prior$sets$cs [[1]][1]]


 CS_delta=  sum(out [[k]]$fit_slide_prior$delta[1, ]*out [[k]]$fit_slide_prior$alpha  [1,])
 CS_delta
 out [[k]]$fit_slide_prior$set
 out [[k]]$fit_slide_prior$delta[,out [[k]]$fit_slide_prior$sets$cs [[1]][1]]


 path="/project2/mstephens/wdenault/susie_mix/results_slide_prior_em/iteration_000/results/"
 lf=list.files(path)
 delta_pos=c()
 delta_lead=c()
for ( k in 1:length( lf) ){

   out <- readRDS(paste0(path, lf[k]))
if ( length(out)>2){
  for (o in 1: length(out )){

    if(!is.null(out [[o]]$fit_slide_prior$sets$cs)){
      CS_delta=  sum(out [[o]]$fit_slide_prior$delta[1, ]*out [[o]]$fit_slide_prior$alpha  [1,])
      CS_delta
      delta_pos=c(delta_pos,CS_delta)

      delta_lead=c(delta_lead, out [[o]]$fit_slide_prior$delta[1,out [[o]]$fit_slide_prior$sets$cs [[1]][1]])
    }
  }
}


print(k)
}


hist(delta_lead)
hist(delta_pos)




path <- "/project2/mstephens/wdenault/susie_mix/results_slide_prior_em/iteration_000/results/"

lf <- list.files(
  path, pattern = "\\.rds$", full.names = TRUE, ignore.case = TRUE
)
length(lf)
# Summaries conditional on the selected SNP, or on membership in the CS.
summarize_delta <- function(weights, grid, prefix) {
  prob <- weights / sum(weights)

  ans <- c(
    delta_mean       = sum(prob * grid),
    p_recessive_like = sum(prob[grid < 0]),
    p_recessive      = sum(prob[grid == -1]),
    p_additive       = sum(prob[grid == 0]),
    p_dominant_like  = sum(prob[grid > 0])
  )

  setNames(ans, paste0(prefix, names(ans)))
}

rows <- list()

for (k in seq_along(lf)) {
  message(k, "/", length(lf), ": ", basename(lf[k]))
  out <- readRDS(lf[k])

  # Preserve your original file-level filter.
  if (length(out) <= 2L) next

  for (o in seq_along(out)) {
    fit <- out[[o]]$fit_slide_prior
    if (is.null(fit$sets$cs) || !length(fit$sets$cs)) next

    tissue <- if (is.null(names(out))) as.character(o) else names(out)[o]

    for (s in seq_along(fit$sets$cs)) {
      cs  <- fit$sets$cs[[s]]
      ell <- fit$sets$cs_index[s]

      # Lead = largest component-specific alpha within this CS.
      alpha_cs <- fit$alpha[ell, cs]
      lead <- cs[which.max(alpha_cs)]

      # Sum joint SNP–delta probabilities over CS members.
      cs_weights <- apply(
        fit$alpha_delta[ell, cs, , drop = FALSE],
        3, sum
      )

      cs_summary <- summarize_delta(
        cs_weights, fit$delta_grid, "cs_"
      )
      lead_summary <- summarize_delta(
        fit$alpha_delta[ell, lead, ], fit$delta_grid, "lead_"
      )

      rows[[length(rows) + 1L]] <- data.frame(
        file        = basename(lf[k]),
        tissue      = tissue,
        cs          = names(fit$sets$cs)[s],
        component   = ell,
        n_snps      = length(cs),
        cs_coverage = sum(alpha_cs),

        lead_index        = lead,
        lead_alpha        = unname(fit$alpha[ell, lead]),
        lead_delta_forced = unname(fit$delta_forced[lead]),
        lead_n0           = unname(fit$genotype_counts[lead, "0"]),
        lead_n1           = unname(fit$genotype_counts[lead, "1"]),
        lead_n2           = unname(fit$genotype_counts[lead, "2"]),
        min_obs           = fit$min_obs,

        cs_forced_mass = sum(
          alpha_cs * fit$delta_forced[cs]
        ) / sum(alpha_cs),

        as.list(cs_summary),
        as.list(lead_summary),
        row.names = NULL
      )
    }
  }
}

results <- if (length(rows)) do.call(rbind, rows) else data.frame()

# Summarize forced status and plot only unforced leads.
if (nrow(results)) {
  print(table(results$lead_delta_forced))

  hist(
    results$lead_delta_mean[!results$lead_delta_forced],
   # breaks = seq(-1, 1, length.out = 101),
    main = "Lead SNPs without forced additive coding",
    xlab = "Posterior mean delta", nclass=200
  )
}


library(ggplot2)
df= data.frame(x=results$lead_delta_mean[!results$lead_delta_forced],
               tissue=as.factor(results$tissue[!results$lead_delta_forced]))




P1= ggplot(df, aes(x = x)) +
  geom_density() +
  geom_density(data = data.frame(x = df$x), colour = "red") +
  facet_wrap(~ tissue) +
  geom_vline(xintercept = 0)+
  ggtitle("lead CS SNP delta")

P1

df_cs= data.frame(x=results$cs_delta_mean[!results$lead_delta_forced],
                  tissue=as.factor(results$tissue[!results$lead_delta_forced]))


P2= ggplot(df_cs, aes(x = x)) +
  geom_density() +
  geom_density(data = data.frame(x = df_cs$x), colour = "red") +
  facet_wrap(~ tissue) +
  geom_vline(xintercept = 0)+
  ggtitle(" CS delta")

P2

gridExtra::grid.arrange(P1,P2, ncol=2)

results <- if (length(rows)) do.call(rbind, rows) else data.frame()


par(mfrow=c(3,1))
hist(results$cs_p_additive[-which(results$cs_delta_mean==0)],
     main= "posterior probability of a CS having delta = 0,  additive",
     nclass = 100,
     xlab = 'posterior probability of a CS having delta = 0')
hist(results$cs_p_recessive_like[-which(results$cs_delta_mean==0)],

     main= "posterior probability of a CS having delta <0, recessive like",
     nclass = 100 ,
     xlab = 'posterior probability of a CS having delta <0')
hist(results$cs_p_dominant_like[-which(results$cs_delta_mean==0)],

     main= "posterior probability of a CS haveing delta >0, dominant like",
     nclass = 100  ,
     xlab = 'posterior probability of a CS having delta >0')



par(mfrow=c(1,1))


#table(results$cs_p_dominant_like +results$cs_p_additive+results$cs_p_recessive_like)
hist(results$lead_p_additive , nclass = 100)
hist(results$lead_p_dominant[-which(results$cs_delta_mean==0)] , nclass = 100)
hist(results$lead_p_recessive[-which(results$cs_delta_mean==0)], nclass = 100)


hist(results$lead_delta_mean[-which(results$lead_delta_mean==0)] , nclass = 150)

quantile(results$lead_delta_mean[-which(results$lead_delta_mean==0)] )

summary(lm(results$lead_delta_mean[-which(results$lead_delta_mean==0)] ~1))

hist(results$cs_delta_mean[-which(results$cs_delta_mean==0)]  , nclass = 100, main= "posterior mean of delta per CS",
    xlab="posterior mean of delta" )

# Your original vectors, now covering every reported CS.
delta_pos  <- as.numeric(results$cs_delta_mean)
delta_lead <- as.numeric(results$lead_delta_mean)



head(results)


