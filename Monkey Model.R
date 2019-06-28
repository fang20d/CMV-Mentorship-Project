my_data <- read_excel("PVL Data.xlsx")

RQc19 <- read_excel("PVL Data.xlsx", range = cell_cols("A:B"))

RQc19 <- RQc19[-c(4,16,19),]

RQc19 <- cbind(RQc19, sd = 0.45)

colnames(RQc19) <- c("time", "logV", "sd")

DataLogV <- as.matrix(RQc19)

dv0 <- as.numeric(RQc19[2,2])-as.numeric(RQc19[1,2])

HIV <- function (pars, V_0 = 60, dV_0 = dv0, T_0 = 100) 
{
  I_0 <- with(as.list(pars), (dV_0 + c * V_0) / (n * delt))
  y <- c(T = T_0, I = I_0, V = V_0)
  times <- c(0,1,2,4,6,8,9,10,12,14,16,18,20,24,28,32)
  out <- ode(y = y, parms = pars, times = times, func = "derivshiv",
             initfunc = "inithiv", nout = 1, outnames = "logV", dllname = "FME")
  
  as.data.frame(out)
}

pars <- c(bet = exp(-5.3), rho = exp(-0.96), delt = exp(0.55), c = exp(3.8), lam = exp(2.47), n = 900)

out <- HIV(pars=pars)

DataT <- as.matrix(cbind(out[,c(1,2)], sd=0.45))

HIVcost <- function (pars) 
{
  out <- HIV(pars)
  cost <- modCost(model = out, obs = DataLogV, err = "sd")
  return(modCost(model = out, obs = DataT, err = "sd", cost = cost))
}

HIVcost2 <- function(lpars)
  HIVcost(c(exp(lpars), n = 900))

Pars <- log(pars[1:5] * 2)
Fit <- modFit(f = HIVcost2, p = Pars)
exp(coef(Fit))

final <- HIV(pars = c(exp(coef(Fit)), n = 900))

plot(RQc19$time, RQc19$logV, xlab = "time", ylab = "logV",log = "y")
lines(final$WeeksPostInfection, final$RQc19)