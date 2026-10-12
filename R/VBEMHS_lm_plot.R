### Functions for visualizing results of VBEMHS_int and VBEHMS

# Making tidy data for trace plot (requires tidyverse packages)
VBEMHS_lm.trace.tidy <- function(res){
  wide <- dplyr::bind_cols(tau2=res$tau2, df_s2=res$df_s2, tibble::as_tibble(res$beta))
  beta_long <- wide |>
    tidyr:: pivot_longer(cols=!c("tau2", "df_s2"), names_to = "Variable", values_to = "beta")
  wide <- dplyr::bind_cols(tau2=res$tau2, tibble::as_tibble(res$se_beta))
  se_beta_long <- wide |> 
    tidyr:: pivot_longer(cols=!tau2, names_to = "Variable", values_to = "se_beta")

  long <- beta_long |> 
    dplyr::full_join(se_beta_long,by=c("tau2","Variable")) |>
    dplyr::mutate(t_stat = ifelse(se_beta>0, beta/se_beta, 0)) |>
    dplyr::mutate(pd = pt(abs(t_stat), df_s2)) 
  long$Variable <- forcats::fct_inorder(long$Variable)

  wide <- dplyr::bind_cols(tau2=res$tau2,
                           sigma2=res$sigma2,
                           lb=res$lb)  
  return(tibble::lst(beta=long, trace=wide))
}

# plot beta.mean with CI (requires ggplot2)
plot.VBEMHS_lm.beta <- function(res, alpha=0.05, limits=FALSE, breaks=FALSE,
                             col.select=c("#D55E00","black")) # Okabe-Ito vermillion
{
  res.sum <- VBEMHS_lm.summary(res, CI=TRUE, alpha=alpha) |> rownames_to_column(var="Variable")
  res.sum$Variable <- forcats::fct_inorder(res.sum$Variable)
  res.sum$Selected <- cut(res.sum$pd, breaks=c(0,1-alpha/2,1), labels=c("No","Yes"))

  if (limits) var.limits <- res.sum$Variable[res.sum$Selected=="Yes"]
    else var.limits <- NULL
  if (breaks) var.breaks <- res.sum$Variable[res.sum$Selected=="Yes"]
    else var.breaks <- waiver()
  
  p <- ggplot(res.sum, aes(x = Variable, y = beta.mean, colour = Selected, shape = Selected)) + 
    scale_x_discrete(limits=var.limits,breaks=var.breaks) + 
    geom_point(size=2, na.rm=TRUE) + theme(legend.position = "none") + 
    geom_errorbar(aes(ymin = beta.lower, ymax = beta.upper), width=0.5) + ylab("beta") +
    scale_color_manual(values = c("Yes" = col.select[1], "No" = col.select[2])) + 
    scale_shape_manual(values=c(21,16)) + 
    geom_hline(yintercept=0,col="darkgray") + 
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5), axis.title.x=element_blank())
  p
}

# Trace plot (requires ggplot2 and scales package)
plot.VBEMHS_lm.trace.beta <- function(res.tidy, alpha=0.05, tau2_opt=NULL,
                                      col.tau2.vline="black", legend=TRUE)
{
  res <- res.tidy$beta
  res$Selected <- cut(abs(res$pd), breaks=c(0,1-alpha/2,1), labels=c("No","Yes"))
  p <- ggplot(res, aes(x = tau2, y = beta, group = Variable, shape = Selected, colour = Variable)) + 
    scale_x_log10(labels=scales::trans_format(log10, scales::math_format(10^.x))) + 
    geom_point(size=2) + scale_shape_manual(values=c(21,16)) + 
    geom_line(linetype="dotted") + scale_colour_hue(l=50)
  if (!legend) p <- p + theme(legend.position = "none")
  if (!is.null(tau2_opt)) 
    p <- p + geom_vline(xintercept=tau2_opt,
                        linetype="dashed",colour=col.tau2.vline,show.legend = FALSE)
  p
}

plot.VBEMHS_lm.trace.tstat <- function(res.tidy, alpha=0.05, tau2_opt=NULL,
                                       col.tau2.vline="black", legend=TRUE)
{
  res <- res.tidy$beta
  res$Selected <- cut(abs(res$pd), breaks=c(0,1-alpha/2,1), labels=c("No","Yes"))
  p <- ggplot(res, aes(x = tau2, y = t_stat, group = Variable, shape = Selected, colour = Variable)) + 
    scale_x_log10(labels=scales::trans_format(log10, scales::math_format(10^.x))) + 
    geom_point(size=2) + scale_shape_manual(values=c(21,16)) + 
    geom_line(linetype="dotted") + scale_colour_hue(l=50)
  if (!legend) p <- p + theme(legend.position = "none")
  if (!is.null(tau2_opt)) 
    p <- p + geom_vline(xintercept=tau2_opt,
                        linetype="dashed",colour=col.tau2.vline,show.legend = FALSE)
  p
}

plot.VBEMHS_lm.trace.sigma2 <- function(res.tidy,tau2_opt=NULL,col.tau2.vline="black")
{
  res <- res.tidy$trace
  p <- ggplot(res, aes(x = tau2, y = sigma2)) + 
    scale_x_log10(labels=scales::trans_format(log10, scales::math_format(10^.x))) + 
    geom_point(size=2) + geom_line(linetype="dotted") + ylim(0,max(res$sigma2))
  if (!is.null(tau2_opt)) 
    p <- p + geom_vline(xintercept=tau2_opt,
                        linetype="dashed",colour=col.tau2.vline,show.legend = FALSE)
  p
}

plot.VBEMHS_lm.trace.lb <- function(res.tidy,tau2_opt=NULL,col.tau2.vline=NULL)
{
  res <- res.tidy$trace
  p <- ggplot(res, aes(x = tau2, y = lb)) + 
    scale_x_log10(labels=scales::trans_format(log10, scales::math_format(10^.x))) + 
    geom_point(size=2) + geom_line(linetype="dotted")
  if (!is.null(tau2_opt)) 
    p <- p + geom_vline(xintercept=tau2_opt,
                        linetype="dashed",colour=col.tau2.vline,show.legend = FALSE)
  p
}

