## Results: avg. parameter, bias, variance, mean^2 error
true_mu <- mu[1]
Ba_hat <- B_hat[,2]
Ba_hat.avg <- mean(Ba_hat)
Ba_hat.bias <- Ba_hat.avg - true_mu
Ba_hat.pcbias <- (Ba_hat.bias / true_mu) * 100
Ba_hat.var <- var(Ba_hat)
Ba_hat.mean2 <- Ba_hat.bias^2 + Ba_hat.var

## Tabulated results
results.summary <- data.frame(mu=true_mu,
                              mean=Ba_hat.avg,
                              bias=Ba_hat.bias,
                              pc_bias=Ba_hat.pcbias,
                              var=Ba_hat.var,
                              mean2_err=Ba_hat.mean2)

## Make histogram
## Ba_hat_hist_full has all n values
## Ba_hat_hist is cleaned out outliers
## Ba_hat_outliers contains the outliers
Ba_hat_hist_full <- data.frame(Ba_hat=Ba_hat)
Ba_hat_hist <- data.frame(Ba_hat=Ba_hat[Ba_hat > Ba_hat.avg - 3*sd(Ba_hat) & Ba_hat < Ba_hat.avg + 3*sd(Ba_hat)])
Ba_hat_outliers <- Ba_hat[Ba_hat < Ba_hat.avg - 3*sd(Ba_hat) | Ba_hat > Ba_hat.avg + 3*sd(Ba_hat)]

## Plot histogram WITHOUT outliers
ggplot(Ba_hat_hist, aes(x=Ba_hat)) +
  geom_histogram(bins = 20) +
  geom_vline(aes(xintercept=mean(Ba_hat)),
             color="red", linetype="dashed", size=1, alpha=.5) +
  ggtitle("Histogram of Ba_hat (no outliers)")

## Plot histogram WITH outliers
ggplot(Ba_hat_hist_full, aes(x=Ba_hat)) +
  geom_histogram(bins = 20) +
  geom_vline(aes(xintercept=mean(Ba_hat)),
             color="red", linetype="dashed", size=1, alpha=.5) +
  ggtitle("Histogram of Ba_hat (all values)")