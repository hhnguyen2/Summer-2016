# Input: Ba_hat, a vector of Ba_hat estimates
# Output: 2 histograms -- one with outliers and one without. 
gen_hist <- function(Ba_hat){
  Ba_hat_hist_full <- data.frame(Ba_hat=Ba_hat)
  ggplot(Ba_hat_hist_full, aes(x=Ba_hat)) +
    geom_histogram(bins = 20) +
    geom_vline(aes(xintercept=mean(Ba_hat)),
               color="red", linetype="dashed", size=1, alpha=.5) +
    ggtitle("Histogram of Ba_hat (all values)")
  
  Ba_hat_hist <- data.frame(Ba_hat=Ba_hat[Ba_hat > mean(Ba_hat) - 3*sd(Ba_hat) & Ba_hat < mean(Ba_hat) + 3*sd(Ba_hat)])
  Ba_hat_outliers <- Ba_hat[Ba_hat < mean(Ba_hat) - 3*sd(Ba_hat) | Ba_hat > mean(Ba_hat) + 3*sd(Ba_hat)]
  
  ## Plot histogram WITHOUT outliers
  ggplot(Ba_hat_hist, aes(x=Ba_hat)) +
    geom_histogram(bins = 20) +
    geom_vline(aes(xintercept=mean(Ba_hat)),
               color="red", linetype="dashed", size=1, alpha=.5) +
    ggtitle("Histogram of Ba_hat (no outliers)")
}