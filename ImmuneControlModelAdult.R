A42870 <- read_excel("MIV05 Plasma VL_RG11202018 .xlsx")
A42870 <- A42870[-c(1:2),c(1,9)]
A42870 <- A42870[rowSums(is.na(A42870)) == 0,]
A42870 <- A42870[-c(8:25),]
A42870 <- cbind(A42870, sd = 0.05)
colnames(A42870) <- c("time", "logV", "sd")
A42870$logV <- as.numeric((A42870$logV))
A42870$logV <- log(A42870$logV)

DataLogV <- as.numeric(A42870)

dv0 <- as.numeric(A42870[2,2])-as.numeric(A42870[1,2])

HIV <- function (pars, V_0 = 15, dV_0 = dv0, T_0 = 3000, E_0 = 20) 
{
  derivs <- function(time, y, pars) 
  {
    with (as.list(c(pars, y)), 
          {
            dT <- lam - rho * T - bet * T * V
            dI <- bet * T * V - k * I * E
            dV <- n * I - c * V
            dE <- alph * I * E - delt * E
            return(list(c(dT, dI, dV, dE), logV = log(V)))
          })
  }
  
  # initial conditions
  I_0 <- with(as.list(pars), (dV_0 + c * V_0) / (n))
  y <- c(T = T_0, I = I_0, V = V_0, E = E_0)
  
  times <- c(seq(0, 8, 0.1))
  out <- ode(y = y, parms = pars, times = times, func = derivs)
  
  as.data.frame(out)
}

pars <- c(bet = 0.000012467	, rho = 0.190546, delt = 2.585710, c = 0.5173521, lam = 354.81, n = 1.443912e+04, alph = 1.2, k = 0.000001)

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
exp(coef(Fit))

final <- HIV(pars = c(exp(coef(Fit))))

a <- ggplot(final) + geom_line(aes(x = time, y = logV), color = "royalblue2") + geom_point(data = A42870, aes(x = as.numeric(time), y = as.numeric(logV)), color = "black")
print(a+ggtitle("logV vs. Time for Model and Data")+ theme(plot.title = element_text(hjust = 0.5))
)
# par(mfrow = c(1, 1))
# plot(A42870$time, A42870$logV, xlab = "time", ylab = "logV",log = "y", main = "Data Points vs. Model Points", col = "red")
# lines(final$time, final$logV)
# legend("bottomright", c("data", "fitted"),
#        lty = c(NA,1), pch = c(1, NA), col = c("red", "black"))

