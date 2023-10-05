library(deSolve)
library(tidyverse)

theme_set(theme_bw())


model <- function(t, y, pars) {
  with(as.list(c(y, pars)), {
    r  <- beta / (1 - risk) * risk / ppv
    
    dy <- c(
      beta * (xf + xt) - r * xn,
      (1 - ppv) * r * xn - ppv * r * xf - beta * xf,
           ppv  * r * xn + ppv * r * xf - beta * xt
    )
    
    return(list(dy, c(Prop = xt / (xt + xf))))
  })
}


p <- c(beta=0.02, r=0.001, ppv=0.61)



fn <- function(risk) {
  ys <- ode(c(xn = 1, xf = 0, xt = 0), seq(0, 1000, 0.5), model, c(beta=1 / 70, risk=risk, ppv=0.61))
  ys[nrow(ys), "Prop"]
}


x <- seq(0.01, 0.5, 0.01)
prs <- sapply(x, fn)



res <- tibble(
  ltr = x,
  prop = prs,
  ppv = p['ppv'],
  prtreated = 0.1 * p['ppv'] / prs
) 

res %>% 
  ggplot() + 
  geom_line(aes(x = ltr, y = prop)) +
  geom_hline(yintercept = 0.61, linetype = 2) +
  scale_x_continuous("Life-time of TB infection, %", labels = scales::percent) +
  scale_y_continuous("True TB treatment history, %", labels = scales::percent) +
  expand_limits(x= 0, y = c(0, 1))


res %>% 
  ggplot() + 
  geom_line(aes(x = ltr, y = prtreated)) +
  scale_x_continuous("Life-time of TB infection, %", labels = scales::percent) +
  scale_y_continuous("True TB treatment history, %", labels = scales::percent) +
  expand_limits(x= 0, y = 0)











