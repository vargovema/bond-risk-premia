if (!require("tidyverse")) install.packages("tidyverse"); library("tidyverse")
if (!require("lmtest")) install.packages("lmtest"); library("lmtest")
if (!require("sandwich")) install.packages("sandwich"); library("sandwich")
if (!require("knitr")) install.packages("knitr"); library("knitr")
if (!require("kableExtra")) install.packages("kableExtra"); library("kableExtra")

load("data/data.RData")
start_date <- (12*11-4) #index
end_date <- start_date + 468 + 11 +12 #index

excess.returns <- function(FBlog, FBm, start_date, end_date, t, tp1) {
  # adding y_1, f_2, ..., f_5
  FBlog <- FBlog %>%
    mutate(
      ym = FBm$ym[start_date:end_date],
      y1 = -X1,
      f2 = X1 - X2,
      f3 = X2 - X3,
      f4 = X3 - X4,
      f5 = X4 - X5
    ) %>%
    rename(
      "p1" = "X1",
      "p2" = "X2",
      "p3" = "X3",
      "p4" = "X4",
      "p5" = "X5"
    )
  
  # Excess returns
  
  FBlog$er1 = rep(0,nrow(FBlog))
  FBlog$er2 = c(
    rep(NA,12),
    (FBlog$p1[tp1]-FBlog$p2[t])-FBlog$y1[t]
  )
  FBlog$er3 = c(
    rep(NA,12),
    (FBlog$p2[tp1]-FBlog$p3[t])-FBlog$y1[t]
  )
  FBlog$er4 = c(
    rep(NA,12),
    (FBlog$p3[tp1]-FBlog$p4[t])-FBlog$y1[t]
  )
  FBlog$er5 = c(
    rep(NA,12),
    (FBlog$p4[tp1]-FBlog$p5[t])-FBlog$y1[t]
  )
  
  FBlog <- FBlog %>%
    mutate(avg_er = (er2+er3+er4+er5)/4) %>%
    mutate(
      fs2 = f2 - y1,
      fs3 = f3 - y1,
      fs4 = f4 - y1,
      fs5 = f5 - y1
    )
}


# data frame of log prices
FBlog <- FBm[start_date:end_date,]
FBlog <- data.frame(
  matrix(c(-(FBlog$fb.y01), -2*(FBlog$fb.y02), -3*(FBlog$fb.y03), -4*(FBlog$fb.y04), -5*(FBlog$fb.y05)), 
         ncol = 5))

t <- 1:(nrow(FBlog)-12)
tp1 <- 13:nrow(FBlog)
FBlog <- excess.returns(FBlog, FBm, start_date, end_date, t, tp1)


# regression of excess returns with maturity 2-5 on forward rates
reg_er2 <- lm(er2[tp1] ~ y1[t] + f2[t] + f3[t] + f4[t] + f5[t], data = FBlog)
reg_er3 <- lm(er3[tp1] ~ y1[t] + f2[t] + f3[t] + f4[t] + f5[t], data = FBlog)
reg_er4 <- lm(er4[tp1] ~ y1[t] + f2[t] + f3[t] + f4[t] + f5[t], data = FBlog)
reg_er5 <- lm(er5[tp1] ~ y1[t] + f2[t] + f3[t] + f4[t] + f5[t], data = FBlog)


# regression of average (across maturities) returns on forward rates
reg_er_avg <- lm(avg_er[tp1] ~ y1[t] + f2[t] + f3[t] + f4[t] + f5[t], data = FBlog)
g <- reg_er_avg$coefficients

# 4 regressions of excess returns (n=3,...,5) on gamma'*f_t
FBlog <- FBlog %>%
  mutate(gf = g[1]+g[2]*y1+g[3]*f2+g[4]*f3+g[5]*f4+g[6]*f5)

reg2_er2 <- lm(er2[tp1] ~ gf[t] - 1, data = FBlog)
reg2_er3 <- lm(er3[tp1] ~ gf[t] - 1, data = FBlog)
reg2_er4 <- lm(er4[tp1] ~ gf[t] - 1, data = FBlog)
reg2_er5 <- lm(er5[tp1] ~ gf[t] - 1, data = FBlog)


# NW,18
reg_er_avg.NeweyWest <- coeftest(reg_er_avg, NeweyWest(reg_er_avg, lag=6))
reg2_er2.NeweyWest <- coeftest(reg2_er2, NeweyWest(reg2_er2, lag=6))
reg2_er3.NeweyWest <- coeftest(reg2_er3, NeweyWest(reg2_er3, lag=6))
reg2_er4.NeweyWest <- coeftest(reg2_er4, NeweyWest(reg2_er4, lag=6))
reg2_er5.NeweyWest <- coeftest(reg2_er5, NeweyWest(reg2_er5, lag=6))

chi_squared_5 <- sum(reg_er_avg.NeweyWest[,3]^2)

png(file="out/reg_coeffs_1.png", width=8, height=4, units="in", res=600)
par(mfrow=c(1,2), mar=c(2.5,2,1.5,0.1), cex=0.6, mgp=c(1.5,0.5,0))
plot(reg_er5$coefficients[2:6], xlab="Forward rates (at time t)", ylab="", main="Unrestricted")
lines(reg_er5$coefficients[2:6])
points(reg_er4$coefficients[2:6],pch=6)
lines(reg_er4$coefficients[2:6],lty=2)
points(reg_er3$coefficients[2:6],pch=0)
lines(reg_er3$coefficients[2:6],lty=3)
points(reg_er2$coefficients[2:6],pch=5)
lines(reg_er2$coefficients[2:6],lty=4)
legend(x="topright", legend=c("5","4","3","2"), lty=c(1,2,3,4), pch=c(1,6,0,5))

plot(reg_er_avg$coefficients[2:6]*reg2_er5$coefficients, xlab="Forward rates (at time t)", ylab="", main="Restricted")
lines(reg_er_avg$coefficients[2:6]*reg2_er5$coefficients)
points(reg_er_avg$coefficients[2:6]*reg2_er4$coefficients,pch=6)
lines(reg_er_avg$coefficients[2:6]*reg2_er4$coefficients,lty=2)
points(reg_er_avg$coefficients[2:6]*reg2_er3$coefficients,pch=0)
lines(reg_er_avg$coefficients[2:6]*reg2_er3$coefficients,lty=3)
points(reg_er_avg$coefficients[2:6]*reg2_er2$coefficients,pch=5)
lines(reg_er_avg$coefficients[2:6]*reg2_er2$coefficients,lty=4)
legend(x="topright", legend=c("5","4","3","2"), lty=c(1,2,3,4), pch=c(1,6,0,5))
dev.off()

table1.a.1 <- rbind(reg_er_avg.NeweyWest[,1], reg_er_avg.NeweyWest[,2])
table1.a.1 <- cbind(table1.a.1, c(NA,summary(reg_er_avg)$r.squared), c(NA,chi_squared_5))
rownames(table1.a.1) <- c("OLS estimates","NW, 18 lags")
opts <- options(knitr.kable.NA = "")
table1.a.1.html <- kbl(table1.a.1, digits = 2,booktabs = TRUE,escape = FALSE, caption="Estimates of the single-factor model using the Fama-Bliss yields.",
    col.names = c("$\\gamma_0$", "$\\gamma_1$", "$\\gamma_2$", "$\\gamma_3$", 
                  "$\\gamma_4$", "$\\gamma_5$","$R^2$","$\\chi^2(5)$")) %>%
  kable_styling(position = "left", latex_options = c("hold_position"), full_width = T, font_size = 10) %>% 
  column_spec(1, width = "4cm")  %>%
  add_header_above(c("A. Estimates of the return-forecasting factor"=9), align="l", bold=T, line=F, line_sep = 5)
save_kable(table1.a.1.html, file = "out/table1.a.1.html")
print(table1.a.1.html)

# Panel B
b_n <- c(reg2_er2$coefficients[1], reg2_er3$coefficients[1],
         reg2_er4$coefficients[1], reg2_er5$coefficients[1])
#large T goes here

r_sq_r <- c(summary(reg2_er2)$r.squared, summary(reg2_er3)$r.squared,
            summary(reg2_er4)$r.squared, summary(reg2_er5)$r.squared)

r_sq_u <- c(summary(reg_er2)$r.squared, summary(reg_er3)$r.squared,
            summary(reg_er4)$r.squared, summary(reg_er5)$r.squared)

large_t <-c(summary(reg2_er2)$coefficients[2], summary(reg2_er3)$coefficients[2], 
            summary(reg2_er4)$coefficients[2], summary(reg2_er5)$coefficients[2])


table1.b.1 <- cbind(2:5, b_n, large_t, r_sq_r, r_sq_u)
rownames(table1.b.1) <- NULL
table1.b.1.html <- kbl(table1.b.1, digits = 2, booktabs = TRUE,escape = FALSE,
    col.names = c("$n$", "$b_n$", "Large $T$","$R^2$","$R^2$")) %>%
  kable_styling(position = "left", latex_options = c("hold_position"), full_width = F, font_size = 10) %>%
  add_header_above(c("Restricted" = 4, "Unrestricted" = 1)) %>%
  add_header_above(c("B. Individual-bond regressions"=5), bold=T, line=F, line_sep = 5, align="l")
save_kable(table1.b.1.html, file = "out/table1.b.1.html")
print(table1.b.1.html)

# Table 2
reg_fb_er2 <- lm(er2[tp1] ~ fs2[t], data = FBlog)
reg_fb_er3 <- lm(er3[tp1] ~ fs3[t], data = FBlog)
reg_fb_er4 <- lm(er4[tp1] ~ fs4[t], data = FBlog)
reg_fb_er5 <- lm(er5[tp1] ~ fs5[t], data = FBlog)

beta_fb <- c(reg_fb_er2$coefficients[2], reg_fb_er3$coefficients[2],
             reg_fb_er4$coefficients[2], reg_fb_er5$coefficients[2])

r_sq_fb <- c(summary(reg_fb_er2)$r.squared, summary(reg_fb_er3)$r.squared,
             summary(reg_fb_er4)$r.squared, summary(reg_fb_er5)$r.squared)

reg_fb_er2.NeweyWest <- coeftest(reg_fb_er2,NeweyWest(reg_fb_er2,lag=6))
reg_fb_er3.NeweyWest <- coeftest(reg_fb_er3,NeweyWest(reg_fb_er3,lag=6))
reg_fb_er4.NeweyWest <- coeftest(reg_fb_er4,NeweyWest(reg_fb_er4,lag=6))
reg_fb_er5.NeweyWest <- coeftest(reg_fb_er5,NeweyWest(reg_fb_er5,lag=6))

reg_fb_er.NeweyWest <- c(reg_fb_er2.NeweyWest[4], reg_fb_er3.NeweyWest[4], 
                         reg_fb_er4.NeweyWest[4], reg_fb_er5.NeweyWest[4])


table2.1 <- cbind(2:5 ,beta_fb,reg_fb_er.NeweyWest, r_sq_fb)
rownames(table2.1) <- NULL
table2.1.html <- kbl(table2.1, digits = 2, booktabs = TRUE, escape = FALSE, caption="Famma-Bliss excess return regressions using the Famma-Bliss yields.",
    col.names = c("Maturity $n$", "$\\beta$", "NW, 18 lags","$R^2$")) %>%
  kable_styling(latex_options = c("hold_position"), full_width = F, font_size = 10)
save_kable(table2.1.html, file = "out/table2.1.html")
print(table2.1.html)

# data frame of log prices
GSWlog <- GSWm[start_date:end_date,]
GSWlog <- data.frame(
  matrix(c(-(GSWlog$GSW.y01), -2*(GSWlog$GSW.y02), -3*(GSWlog$GSW.y03), -4*(GSWlog$GSW.y04), -5*(GSWlog$GSW.y05)),
         ncol = 5))

FBlog <- excess.returns(GSWlog, GSWm, start_date, end_date, t, tp1)


# regression of excess returns with maturity 2-5 on forward rates
reg_er2 <- lm(er2[tp1] ~ y1[t] + f2[t] + f3[t] + f4[t] + f5[t], data = FBlog)
reg_er3 <- lm(er3[tp1] ~ y1[t] + f2[t] + f3[t] + f4[t] + f5[t], data = FBlog)
reg_er4 <- lm(er4[tp1] ~ y1[t] + f2[t] + f3[t] + f4[t] + f5[t], data = FBlog)
reg_er5 <- lm(er5[tp1] ~ y1[t] + f2[t] + f3[t] + f4[t] + f5[t], data = FBlog)


# regression of average (across maturities) returns on forward rates
reg_er_avg <- lm(avg_er[tp1] ~ y1[t] + f2[t] + f3[t] + f4[t] + f5[t], data = FBlog)

# 4 regressions of excess returns (n=3,...,5) on gamma'*f_t
FBlog <- FBlog %>%
  mutate(
    gf = g[1]+g[2]*y1+g[3]*f2+g[4]*f3+g[5]*f4+g[6]*f5
  )

reg2_er2 <- lm(er2[tp1] ~ gf[t] - 1, data = FBlog)
reg2_er3 <- lm(er3[tp1] ~ gf[t] - 1, data = FBlog)
reg2_er4 <- lm(er4[tp1] ~ gf[t] - 1, data = FBlog)
reg2_er5 <- lm(er5[tp1] ~ gf[t] - 1, data = FBlog)


# NW,18
reg_er_avg.NeweyWest <- coeftest(reg_er_avg, NeweyWest(reg_er_avg, lag=6))
reg2_er2.NeweyWest <- coeftest(reg2_er2, NeweyWest(reg2_er2, lag=6))
reg2_er3.NeweyWest <- coeftest(reg2_er3, NeweyWest(reg2_er3, lag=6))
reg2_er4.NeweyWest <- coeftest(reg2_er4, NeweyWest(reg2_er4, lag=6))
reg2_er5.NeweyWest <- coeftest(reg2_er5, NeweyWest(reg2_er5, lag=6))

chi_squared_5 <- sum(reg_er_avg.NeweyWest[,3]^2)

png(file="out/reg_coeffs_2.png", width=8, height=4, units="in", res=600)
par(mfrow=c(1,2), mar=c(2.5,2,1.5,0.1), cex=0.6, mgp=c(1.5,0.5,0))
plot(reg_er5$coefficients[2:6], xlab="Forward rates (at time t)", ylab="", main="Unrestricted")
lines(reg_er5$coefficients[2:6])
points(reg_er4$coefficients[2:6],pch=6)
lines(reg_er4$coefficients[2:6],lty=2)
points(reg_er3$coefficients[2:6],pch=0)
lines(reg_er3$coefficients[2:6],lty=3)
points(reg_er2$coefficients[2:6],pch=5)
lines(reg_er2$coefficients[2:6],lty=4)
legend(x="bottomright", legend=c("5","4","3","2"), lty=c(1,2,3,4), pch=c(1,6,0,5))

plot(reg_er_avg$coefficients[2:6]*reg2_er5$coefficients, xlab="Forward rates (at time t)", ylab="", main="Restricted")
lines(reg_er_avg$coefficients[2:6]*reg2_er5$coefficients)
points(reg_er_avg$coefficients[2:6]*reg2_er4$coefficients,pch=6)
lines(reg_er_avg$coefficients[2:6]*reg2_er4$coefficients,lty=2)
points(reg_er_avg$coefficients[2:6]*reg2_er3$coefficients,pch=0)
lines(reg_er_avg$coefficients[2:6]*reg2_er3$coefficients,lty=3)
points(reg_er_avg$coefficients[2:6]*reg2_er2$coefficients,pch=5)
lines(reg_er_avg$coefficients[2:6]*reg2_er2$coefficients,lty=4)
legend(x="bottomright", legend=c("5","4","3","2"), lty=c(1,2,3,4), pch=c(1,6,0,5))
dev.off()

table1.a.2 <- rbind(reg_er_avg.NeweyWest[,1], reg_er_avg.NeweyWest[,2])
table1.a.2 <- cbind(table1.a.2, c(NA,summary(reg_er_avg)$r.squared), c(NA,chi_squared_5))
rownames(table1.a.2) <- c("OLS estimates","NW, 18 lags")
opts <- options(knitr.kable.NA = "")
table1.a.2.html <- kbl(table1.a.2, digits = 2,booktabs = TRUE,escape = FALSE, caption="Estimates of the single-factor model using GSW yields.",
    col.names = c("$\\gamma_0$", "$\\gamma_1$", "$\\gamma_2$", "$\\gamma_3$", 
                  "$\\gamma_4$", "$\\gamma_5$","$R^2$","$\\chi^2(5)$")) %>%
  kable_styling(position = "left", latex_options = c("hold_position"), full_width = T, font_size = 10) %>% 
  column_spec(1, width = "4cm")  %>%
  add_header_above(c("A. Estimates of the return-forecasting factor"=9), align="l", bold=T, line=F, line_sep = 5)
save_kable(table1.a.2.html, file = "out/table1.a.2.html")
print(table1.a.2.html)

# Panel B
b_n <- c(reg2_er2$coefficients[1], reg2_er3$coefficients[1],
         reg2_er4$coefficients[1], reg2_er5$coefficients[1])

r_sq_r <- c(summary(reg2_er2)$r.squared, summary(reg2_er3)$r.squared,
            summary(reg2_er4)$r.squared, summary(reg2_er5)$r.squared)

r_sq_u <- c(summary(reg_er2)$r.squared, summary(reg_er3)$r.squared, 
            summary(reg_er4)$r.squared, summary(reg_er5)$r.squared)

large_t <-c(summary(reg2_er2)$coefficients[2], summary(reg2_er3)$coefficients[2], 
            summary(reg2_er4)$coefficients[2], summary(reg2_er5)$coefficients[2])

table1.b.2 <- cbind(2:5, b_n, large_t, r_sq_r, r_sq_u)
rownames(table1.b.2) <- NULL
table1.b.2.html <- kbl(table1.b.2, digits = 2, booktabs = TRUE,escape = FALSE,
    col.names = c("$n$", "$b_n$", "Large $T$","$R^2$","$R^2$")) %>%
  kable_styling(position = "left", latex_options = c("hold_position"), full_width = F, font_size = 10) %>%
  add_header_above(c("Restricted" = 4, "Unrestricted" = 1)) %>%
  add_header_above(c("B. Individual Bond regressions"=5), bold=T, line=F, line_sep = 5, align="l")
save_kable(table1.b.2.html, file = "out/table1.b.2.html")
print(table1.b.2.html)

# Table 2
reg_fb_er2 <- lm(er2[tp1] ~ fs2[t], data = FBlog)
reg_fb_er3 <- lm(er3[tp1] ~ fs3[t], data = FBlog)
reg_fb_er4 <- lm(er4[tp1] ~ fs4[t], data = FBlog)
reg_fb_er5 <- lm(er5[tp1] ~ fs5[t], data = FBlog)

beta_fb <- c(reg_fb_er2$coefficients[2], reg_fb_er3$coefficients[2], 
             reg_fb_er4$coefficients[2], reg_fb_er5$coefficients[2])

r_sq_fb <- c(summary(reg_fb_er2)$r.squared, summary(reg_fb_er3)$r.squared, 
             summary(reg_fb_er4)$r.squared, summary(reg_fb_er5)$r.squared)

reg_fb_er2.NeweyWest <- coeftest(reg_fb_er2,NeweyWest(reg_fb_er2,lag=6))
reg_fb_er3.NeweyWest <- coeftest(reg_fb_er3,NeweyWest(reg_fb_er3,lag=6))
reg_fb_er4.NeweyWest <- coeftest(reg_fb_er4,NeweyWest(reg_fb_er4,lag=6))
reg_fb_er5.NeweyWest <- coeftest(reg_fb_er5,NeweyWest(reg_fb_er5,lag=6))

reg_fb_er.NeweyWest <- c(reg_fb_er2.NeweyWest[4], reg_fb_er3.NeweyWest[4], 
                         reg_fb_er4.NeweyWest[4], reg_fb_er5.NeweyWest[4])

table2.2 <- cbind(2:5 ,beta_fb,reg_fb_er.NeweyWest, r_sq_fb)
rownames(table2.2) <- NULL
table2.2.html <- kbl(table2.2, digits = 2, booktabs = TRUE, escape = FALSE, caption="Famma-Bliss excess return regressions using GSW yields.",
    col.names = c("Maturity $n$", "$\\beta$", "NW, 18 lags","$R^2$")) %>%
  kable_styling(latex_options = c("hold_position"), full_width = F, font_size = 10)
save_kable(table2.2.html, file = "out/table2.2.html")
print(table2.2.html)

start_date <- (12*11-4) #index
end_date <- start_date + 719 #index
#FBm$ym[c(start_date,end_date)] 

# data frame of log prices
FBlog <- FBm[start_date:end_date,]
FBlog <- data.frame(
  matrix(c(-(FBlog$fb.y01), -2*(FBlog$fb.y02), -3*(FBlog$fb.y03), -4*(FBlog$fb.y04), -5*(FBlog$fb.y05)),
         ncol = 5))

t <- 1:(nrow(FBlog)-12)
tp1 <- 13:nrow(FBlog)
FBlog <- excess.returns(FBlog, FBm, start_date, end_date, t, tp1)

# regression of excess returns with maturity 2-5 on forward rates
reg_er2 <- lm(er2[tp1] ~ y1[t] + f2[t] + f3[t] + f4[t] + f5[t], data = FBlog)
reg_er3 <- lm(er3[tp1] ~ y1[t] + f2[t] + f3[t] + f4[t] + f5[t], data = FBlog)
reg_er4 <- lm(er4[tp1] ~ y1[t] + f2[t] + f3[t] + f4[t] + f5[t], data = FBlog)
reg_er5 <- lm(er5[tp1] ~ y1[t] + f2[t] + f3[t] + f4[t] + f5[t], data = FBlog)

# regression of average (across maturities) returns on forward rates
reg_er_avg <- lm(avg_er[tp1] ~ y1[t] + f2[t] + f3[t] + f4[t] + f5[t], data = FBlog)
g <- reg_er_avg$coefficients

# 4 regressions of excess returns (n=3,...,5) on gamma'*f_t
FBlog <- FBlog %>%
  mutate(gf = g[1]+g[2]*y1+g[3]*f2+g[4]*f3+g[5]*f4+g[6]*f5)

reg2_er2 <- lm(er2[tp1] ~ gf[t] - 1, data = FBlog)
reg2_er3 <- lm(er3[tp1] ~ gf[t] - 1, data = FBlog)
reg2_er4 <- lm(er4[tp1] ~ gf[t] - 1, data = FBlog)
reg2_er5 <- lm(er5[tp1] ~ gf[t] - 1, data = FBlog)


# NW,18
reg_er_avg.NeweyWest <- coeftest(reg_er_avg, NeweyWest(reg_er_avg, lag=6))
reg2_er2.NeweyWest <- coeftest(reg2_er2, NeweyWest(reg2_er2, lag=6))
reg2_er3.NeweyWest <- coeftest(reg2_er3, NeweyWest(reg2_er3, lag=6))
reg2_er4.NeweyWest <- coeftest(reg2_er4, NeweyWest(reg2_er4, lag=6))
reg2_er5.NeweyWest <- coeftest(reg2_er5, NeweyWest(reg2_er5, lag=6))

chi_squared_5 <- sum(reg_er_avg.NeweyWest[,3]^2)

png(file="out/reg_coeffs_3.png", width=8, height=4, units="in", res=600)
par(mfrow=c(1,2), mar=c(2.5,2,1.5,0.1), cex=0.6, mgp=c(1.5,0.5,0))
plot(reg_er5$coefficients[2:6], xlab="Forward rates (at time t)", ylab="", main="Unrestricted")
lines(reg_er5$coefficients[2:6])
points(reg_er4$coefficients[2:6],pch=6)
lines(reg_er4$coefficients[2:6],lty=2)
points(reg_er3$coefficients[2:6],pch=0)
lines(reg_er3$coefficients[2:6],lty=3)
points(reg_er2$coefficients[2:6],pch=5)
lines(reg_er2$coefficients[2:6],lty=4)
legend(x="topright", legend=c("5","4","3","2"), lty=c(1,2,3,4), pch=c(1,6,0,5))

plot(reg_er_avg$coefficients[2:6]*reg2_er5$coefficients, xlab="Forward rates (at time t)", ylab="", main="Restricted")
lines(reg_er_avg$coefficients[2:6]*reg2_er5$coefficients)
points(reg_er_avg$coefficients[2:6]*reg2_er4$coefficients,pch=6)
lines(reg_er_avg$coefficients[2:6]*reg2_er4$coefficients,lty=2)
points(reg_er_avg$coefficients[2:6]*reg2_er3$coefficients,pch=0)
lines(reg_er_avg$coefficients[2:6]*reg2_er3$coefficients,lty=3)
points(reg_er_avg$coefficients[2:6]*reg2_er2$coefficients,pch=5)
lines(reg_er_avg$coefficients[2:6]*reg2_er2$coefficients,lty=4)
legend(x="topright", legend=c("5","4","3","2"), lty=c(1,2,3,4), pch=c(1,6,0,5))
dev.off()

table1.a.3 <- rbind(reg_er_avg.NeweyWest[,1], reg_er_avg.NeweyWest[,2])
table1.a.3 <- cbind(table1.a.3, c(NA,summary(reg_er_avg)$r.squared), c(NA,chi_squared_5))
rownames(table1.a.3) <- c("OLS estimates","NW, 18 lags")
opts <- options(knitr.kable.NA = "")
table1.a.3.html <- kbl(table1.a.3, digits = 2,booktabs = TRUE,escape = FALSE, caption="Estimates of the single-factor model the recent Fama-Bliss yields.",
    col.names = c("$\\gamma_0$", "$\\gamma_1$", "$\\gamma_2$", "$\\gamma_3$", 
                  "$\\gamma_4$", "$\\gamma_5$","$R^2$","$\\chi^2(5)$")) %>%
  kable_styling(position = "left", latex_options = c("hold_position"), full_width = T, font_size = 10) %>% 
  column_spec(1, width = "4cm")  %>%
  add_header_above(c("A. Estimates of the return-forecasting factor"=9), align="l", bold=T, line=F, line_sep = 5)
save_kable(table1.a.3.html, file = "out/table1.a.3.html")
print(table1.a.3.html)

# Panel B
b_n <- c(reg2_er2$coefficients[1], reg2_er3$coefficients[1],
         reg2_er4$coefficients[1], reg2_er5$coefficients[1])

r_sq_r <- c(summary(reg2_er2)$r.squared, summary(reg2_er3)$r.squared,
            summary(reg2_er4)$r.squared, summary(reg2_er5)$r.squared)

r_sq_u <- c(summary(reg_er2)$r.squared, summary(reg_er3)$r.squared,
            summary(reg_er4)$r.squared, summary(reg_er5)$r.squared)

large_t <-c(summary(reg2_er2)$coefficients[2], summary(reg2_er3)$coefficients[2], 
            summary(reg2_er4)$coefficients[2], summary(reg2_er5)$coefficients[2])

table1.b.3 <- cbind(2:5, b_n, large_t, r_sq_r, r_sq_u)
rownames(table1.b.3) <- NULL
table1.b.3.html <- kbl(table1.b.3, digits = 2, booktabs = TRUE,escape = FALSE,
    col.names = c("$n$", "$b_n$", "Large $T$","$R^2$","$R^2$")) %>%
  kable_styling(position = "left", latex_options = c("hold_position"), full_width = F, font_size = 10) %>%
  add_header_above(c("Restricted" = 4, "Unrestricted" = 1)) %>%
  add_header_above(c("B. Individual-bond regressions"=5), bold=T, line=F, line_sep = 5, align="l")
save_kable(table1.b.3.html, file = "out/table1.b.3.html")
print(table1.b.3.html)

# Table 2
reg_fb_er2 <- lm(er2[tp1] ~ fs2[t], data = FBlog)
reg_fb_er3 <- lm(er3[tp1] ~ fs3[t], data = FBlog)
reg_fb_er4 <- lm(er4[tp1] ~ fs4[t], data = FBlog)
reg_fb_er5 <- lm(er5[tp1] ~ fs5[t], data = FBlog)

beta_fb <- c(reg_fb_er2$coefficients[2], reg_fb_er3$coefficients[2],
             reg_fb_er4$coefficients[2], reg_fb_er5$coefficients[2])

r_sq_fb <- c(summary(reg_fb_er2)$r.squared, summary(reg_fb_er3)$r.squared,
             summary(reg_fb_er4)$r.squared, summary(reg_fb_er5)$r.squared)

reg_fb_er2.NeweyWest <- coeftest(reg_fb_er2,NeweyWest(reg_fb_er2,lag=6))
reg_fb_er3.NeweyWest <- coeftest(reg_fb_er3,NeweyWest(reg_fb_er3,lag=6))
reg_fb_er4.NeweyWest <- coeftest(reg_fb_er4,NeweyWest(reg_fb_er4,lag=6))
reg_fb_er5.NeweyWest <- coeftest(reg_fb_er5,NeweyWest(reg_fb_er5,lag=6))

reg_fb_er.NeweyWest <- c(reg_fb_er2.NeweyWest[4], reg_fb_er3.NeweyWest[4], 
                         reg_fb_er4.NeweyWest[4], reg_fb_er5.NeweyWest[4])

table2.3 <- cbind(2:5 ,beta_fb,reg_fb_er.NeweyWest, r_sq_fb)
rownames(table2.3) <- NULL
table2.3.html <- kbl(table2.3, digits = 2, booktabs = TRUE, escape = FALSE, caption="Famma-Bliss excess return regressions using the recent Famma-Bliss yields.",
    col.names = c("Maturity $n$", "$\\beta$", "NW, 18 lags","$R^2$")) %>%
  kable_styling(latex_options = c("hold_position"), full_width = F, font_size = 10)
save_kable(table2.3.html, file = "out/table2.3.html")
print(table2.3.html)


FBlog$er2.pred_4f <- NA
FBlog$er3.pred_4f <- NA
FBlog$er4.pred_4f <- NA
FBlog$er5.pred_4f <- NA

FBlog$avg.pred <- NA

FBlog$er2.pred_1f <- NA
FBlog$er3.pred_1f <- NA
FBlog$er4.pred_1f <- NA
FBlog$er5.pred_1f <- NA

FBlog$er2.pred_1fe <- NA
FBlog$er3.pred_1fe <- NA
FBlog$er4.pred_1fe <- NA
FBlog$er5.pred_1fe <- NA

FBlog$avg_er.hist <- NA
FBlog$er2.hist <- NA
FBlog$er3.hist <- NA
FBlog$er4.hist <- NA
FBlog$er5.hist <- NA

t0 <- which(FBlog$ym == "Jan 2004")

for(i in t0:nrow(FBlog)) {
  t <- 1:(nrow(FBlog[1:(i-1),])-12)
  tp1 <- 13:nrow(FBlog[1:(i-1),])
  
  # regression of excess returns with maturity 2-5 on forward rates
  reg_er2 <- lm(er2[tp1] ~ y1[t] + f2[t] + f3[t] + f4[t] + f5[t], data = FBlog[1:(i-1),])
  reg_er3 <- lm(er3[tp1] ~ y1[t] + f2[t] + f3[t] + f4[t] + f5[t], data = FBlog[1:(i-1),])
  reg_er4 <- lm(er4[tp1] ~ y1[t] + f2[t] + f3[t] + f4[t] + f5[t], data = FBlog[1:(i-1),])
  reg_er5 <- lm(er5[tp1] ~ y1[t] + f2[t] + f3[t] + f4[t] + f5[t], data = FBlog[1:(i-1),])
  
  reg_er_avg <- lm(avg_er[tp1] ~ y1[t] + f2[t] + f3[t] + f4[t] + f5[t], data = FBlog[1:(i-1),])
  g <- reg_er_avg$coefficients
  
  FBlog <- FBlog %>% mutate(gf = g[1]+g[2]*y1+g[3]*f2+g[4]*f3+g[5]*f4+g[6]*f5)
  
  reg2_er2 <- lm(er2[tp1] ~ gf[t] - 1, data = FBlog[1:(i-1),])
  reg2_er3 <- lm(er3[tp1] ~ gf[t] - 1, data = FBlog[1:(i-1),])
  reg2_er4 <- lm(er4[tp1] ~ gf[t] - 1, data = FBlog[1:(i-1),])
  reg2_er5 <- lm(er5[tp1] ~ gf[t] - 1, data = FBlog[1:(i-1),])
  
  reg_fb_er2 <- lm(er2[tp1] ~ fs2[t], data = FBlog[1:(i-1),])
  reg_fb_er3 <- lm(er3[tp1] ~ fs3[t], data = FBlog[1:(i-1),])
  reg_fb_er4 <- lm(er4[tp1] ~ fs4[t], data = FBlog[1:(i-1),])
  reg_fb_er5 <- lm(er5[tp1] ~ fs5[t], data = FBlog[1:(i-1),])
  
  FBlog$er2.pred_4f[i] <- predict(reg_er2, newdata = FBlog[i,])
  FBlog$er3.pred_4f[i] <- predict(reg_er3, newdata = FBlog[i,])
  FBlog$er4.pred_4f[i] <- predict(reg_er4, newdata = FBlog[i,])
  FBlog$er5.pred_4f[i] <- predict(reg_er5, newdata = FBlog[i,])
  
  FBlog$avg.pred[i] <- predict(reg_er_avg, newdata = FBlog[i,])
  
  FBlog$er2.pred_1f[i] <- predict(reg2_er2, newdata = FBlog[i,])
  FBlog$er3.pred_1f[i] <- predict(reg2_er3, newdata = FBlog[i,])
  FBlog$er4.pred_1f[i] <- predict(reg2_er4, newdata = FBlog[i,])
  FBlog$er5.pred_1f[i] <- predict(reg2_er5, newdata = FBlog[i,])
  
  FBlog$er2.pred_1fe[i] <- predict(reg_fb_er2, newdata = FBlog[i,])
  FBlog$er3.pred_1fe[i] <- predict(reg_fb_er3, newdata = FBlog[i,])
  FBlog$er4.pred_1fe[i] <- predict(reg_fb_er4, newdata = FBlog[i,])
  FBlog$er5.pred_1fe[i] <- predict(reg_fb_er5, newdata = FBlog[i,])
  
  FBlog$avg_er.hist[i] <- mean(FBlog$avg_er[1:(i-1)], na.rm=T)
  
  FBlog$er2.hist[i] <- mean(FBlog$er2[1:(i-1)], na.rm=T)
  FBlog$er3.hist[i] <- mean(FBlog$er3[1:(i-1)], na.rm=T)
  FBlog$er4.hist[i] <- mean(FBlog$er4[1:(i-1)], na.rm=T)
  FBlog$er5.hist[i] <- mean(FBlog$er5[1:(i-1)], na.rm=T)
  
}

r.squared.os.er2.4f <- 1 - sum((FBlog$er2[t0:nrow(FBlog)] - FBlog$er2.pred_4f[t0:nrow(FBlog)])^2, na.rm = TRUE) /
  sum((FBlog$er2[t0:nrow(FBlog)] - FBlog$er2.hist[t0:nrow(FBlog)])^2, na.rm = TRUE) 

r.squared.os.er3.4f <- 1 - sum((FBlog$er3[t0:nrow(FBlog)] - FBlog$er3.pred_4f[t0:nrow(FBlog)])^2, na.rm = TRUE) /
  sum((FBlog$er3[t0:nrow(FBlog)] - FBlog$er3.hist[t0:nrow(FBlog)])^2, na.rm = TRUE) 

r.squared.os.er4.4f <- 1 - sum((FBlog$er4[t0:nrow(FBlog)] - FBlog$er4.pred_4f[t0:nrow(FBlog)])^2, na.rm = TRUE) /
  sum((FBlog$er4[t0:nrow(FBlog)] - FBlog$er4.hist[t0:nrow(FBlog)])^2, na.rm = TRUE) 

r.squared.os.er5.4f <- 1 - sum((FBlog$er5[t0:nrow(FBlog)] - FBlog$er5.pred_4f[t0:nrow(FBlog)])^2, na.rm = TRUE) /
  sum((FBlog$er5[t0:nrow(FBlog)]-FBlog$er5.hist[t0:nrow(FBlog)])^2, na.rm = TRUE) 

r.squared.os.4f <- as.data.frame(t(c(r.squared.os.er2.4f, r.squared.os.er3.4f, 
                                     r.squared.os.er4.4f, r.squared.os.er5.4f)))

r.squared.os.avg <- 1 - sum((FBlog$avg_er[t0:nrow(FBlog)] - FBlog$avg.pred[t0:nrow(FBlog)])^2, na.rm = TRUE) /
  sum((FBlog$avg_er[t0:nrow(FBlog)]-FBlog$avg_er.hist[t0:nrow(FBlog)])^2, na.rm = TRUE) 




r.squared.os.er2.1f <- 1 - sum((FBlog$er2[t0:nrow(FBlog)] - FBlog$er2.pred_1f[t0:nrow(FBlog)])^2, na.rm = TRUE) /
  sum((FBlog$er2[t0:nrow(FBlog)] - FBlog$er2.hist[t0:nrow(FBlog)])^2, na.rm = TRUE) 

r.squared.os.er3.1f <- 1 - sum((FBlog$er3[t0:nrow(FBlog)] - FBlog$er3.pred_1f[t0:nrow(FBlog)])^2, na.rm = TRUE) /
  sum((FBlog$er3[t0:nrow(FBlog)] - FBlog$er3.hist[t0:nrow(FBlog)])^2, na.rm = TRUE) 

r.squared.os.er4.1f <- 1 - sum((FBlog$er4[t0:nrow(FBlog)] - FBlog$er4.pred_1f[t0:nrow(FBlog)])^2, na.rm = TRUE) /
  sum((FBlog$er4[t0:nrow(FBlog)] - FBlog$er4.hist[t0:nrow(FBlog)])^2, na.rm = TRUE) 

r.squared.os.er5.1f <- 1 - sum((FBlog$er5[t0:nrow(FBlog)] - FBlog$er5.pred_1f[t0:nrow(FBlog)])^2, na.rm = TRUE) /
  sum((FBlog$er5[t0:nrow(FBlog)]-FBlog$er5.hist[t0:nrow(FBlog)])^2, na.rm = TRUE) 

r.squared.os.1f <- as.data.frame(t(c(r.squared.os.er2.1f, r.squared.os.er3.1f, 
                                     r.squared.os.er4.1f, r.squared.os.er5.1f)))




r.squared.os.er2.1fe <- 1 - sum((FBlog$er2[t0:nrow(FBlog)] - FBlog$er2.pred_1fe[t0:nrow(FBlog)])^2, na.rm = TRUE) /
  sum((FBlog$er2[t0:nrow(FBlog)] - FBlog$er2.hist[t0:nrow(FBlog)])^2, na.rm = TRUE) 

r.squared.os.er3.1fe <- 1 - sum((FBlog$er3[t0:nrow(FBlog)] - FBlog$er3.pred_1fe[t0:nrow(FBlog)])^2, na.rm = TRUE) /
  sum((FBlog$er3[t0:nrow(FBlog)] - FBlog$er3.hist[t0:nrow(FBlog)])^2, na.rm = TRUE) 

r.squared.os.er4.1fe <- 1 - sum((FBlog$er4[t0:nrow(FBlog)] - FBlog$er4.pred_1fe[t0:nrow(FBlog)])^2, na.rm = TRUE) /
  sum((FBlog$er4[t0:nrow(FBlog)] - FBlog$er4.hist[t0:nrow(FBlog)])^2, na.rm = TRUE) 

r.squared.os.er5.1fe <- 1 - sum((FBlog$er5[t0:nrow(FBlog)] - FBlog$er5.pred_1fe[t0:nrow(FBlog)])^2, na.rm = TRUE) /
  sum((FBlog$er5[t0:nrow(FBlog)]-FBlog$er5.hist[t0:nrow(FBlog)])^2, na.rm = TRUE) 

r.squared.os.1fe <- as.data.frame(t(c(r.squared.os.er2.1fe, r.squared.os.er3.1fe, 
                                      r.squared.os.er4.1fe, r.squared.os.er5.1fe)))


r.squared.os <- cbind(r.squared.os.4f, r.squared.os.avg, r.squared.os.1f, r.squared.os.1fe)
rownames(r.squared.os) <- NULL
r.squared.os.html <- kbl(r.squared.os, digits = 2, booktabs = TRUE, escape = FALSE, caption="Out-of-sampe performance evaluation using $R^2_{OS,r_2}$.", align="c", col.names = c("n=2", "n=3", "n=4","n=5","$\\bar {rx}$","n=2", "n=3", "n=4","n=5","n=2", "n=3", "n=4","n=5")) %>%
  kable_styling(latex_options = c("hold_position"), full_width = F, font_size = 9) %>%
  add_header_above(linebreak(c("Unrestricted \n model"=4,"Return-forecasting \n factor"=1,"Restricted \n model"=4,"Excess return \n regressions"=4), linebreaker = "\n", align="c"), bold=T) %>%
  add_header_above(c("Cochrane & Piazzesi (2005) models"=9,"Fama & Bliss (1987) model"=4))
save_kable(r.squared.os.html, file = "out/r.squared.os.html")
print(r.squared.os.html)
