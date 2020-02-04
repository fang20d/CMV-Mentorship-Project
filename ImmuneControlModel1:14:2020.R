M46346 <- read_excel("ARV01 Plasma VL_RG11202018_in tabs.xlsx", sheet = 1)
M46346 <- M46346[-c(1:4),c(3,7)]
M46346 <- M46346[rowSums(is.na(M46346)) == 0,]
M46346 <- M46346[-c(10:27),]
M46346 <- cbind(M46346, sd = 0.05)
colnames(M46346) <- c("time", "logV", "sd")
M46346$logV <- log(M46346$logV)

DataLogV <- as.matrix(datas)

dv0 <- as.numeric(datas[2,2])-as.numeric(datas[1,2])

HIV <- function (pars, V_0 = as.numeric(datas[1,2]), dV_0 = dv0, T_0 = 5, E_0 = 3.301) 
{
  derivs <- function(time, y, pars) 
  {
    with (as.list(c(pars, y)), 
          {
            dT <- lam - rho * T - bet * T * V
            dI <- bet * T * V - k * I * E
            dV <- p * I - c * V
            dE <- alph * I * E - delt * E
            return(list(c(dT, dI, dV, dE), logV = V))
          })
  }
  
  # initial conditions
  I_0 <- with(as.list(pars), (dV_0 + c * V_0) / (p))
  y <- c(T = T_0, I = I_0, V = V_0, E = E_0)
  
  times <- c(seq(0, 8, 0.01))
  out <- ode(y = y, parms = pars, times = times, func = derivs)
  
  as.data.frame(out)
}

pars <- c(bet = 1.4467	, rho = 0.190546, delt = 2.585710, c = 0.5173521, lam = 354.81, p = 1.443912e+04, alph = 1.2, k = 0.000001)

out <- HIV(pars=pars)

DataT <- as.matrix(cbind(out[,c(1,2)], sd=0.45))

HIVcost <- function (pars) 
{
  out <- HIV(pars)
  cost <- modCost(model = out, obs = DataLogV, err = "sd")
  return(modCost(model = out, obs = DataT, err = "sd", cost = cost))
}

HIVcost2 <- function(lpars)
  HIVcost(c(exp(lpars)))

Pars <- log(pars[1:8] * 2)
Fit <- modFit(f = HIVcost2, p = Pars)
coef(Fit)

final <- HIV(pars = c(coef(Fit)))

a <- ggplot(final) + geom_line(aes(x = time, y = logV), color = "royalblue2") + geom_point(data = M46346, aes(x = time, y = logV), color = "black")
print(a+ggtitle("logV vs. Time for Model and Data")+ theme(plot.title = element_text(hjust = 0.5))
)
summary(Fit)
# par(mfrow = c(1, 1))
# plot(M46346$time, M46346$logV, xlab = "time", ylab = "logV",log = "y", main = "Data Points vs. Model Points", col = "red")
# lines(final$time, final$logV)
# legend("bottomright", c("data", "fitted"),
#        lty = c(NA,1), pch = c(1, NA), col = c("red", "black"))

