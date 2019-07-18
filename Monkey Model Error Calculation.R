my_data <- read_excel("PVL Data.xlsx")

Sim <- read_excel("PVL Data.xlsx")

Sets <- read_excel("Monkey data fitting_DOE_no results (1).xlsx")

Sim <- Sim[,c(1,8)]

Sim <- Sim[-c(4,16,19),]

Sim <- cbind(Sim, sd = 0.45)

colnames(Sim) <- c("time", "logV", "sd")

Sim <- Sim[c(1:6),]

Sim$logV <- log(Sim$logV)

DataLogV <- as.matrix(Sim)

dv0 <- as.numeric(Sim[2,2])-as.numeric(Sim[1,2])

for (row in 1:nrow(Sets))
{
  HIV <- function (pars, V_0 = 60, dV_0 = dv0, T_0 = 100) 
  {
    derivs <- function(time, y, pars) 
    {
      with (as.list(c(pars, y)), 
            {
              dT <- lam - rho * T - bet * T * V
              dI <- bet * T * V - delt * I
              dV <- n * delt * I - c * V
              return(list(c(dT, dI, dV), logV = log(V)))
            })
    }
    
    # initial conditions
    I_0 <- with(as.list(pars), (dV_0 + c * V_0) / (n * delt))
    y <- c(T = T_0, I = I_0, V = V_0)
    
    times <- c(seq(0, 8, 1))
    out <- ode(y = y, parms = pars, times = times, func = derivs)
    
    as.data.frame(out)
  }
  
  pars <- c(bet = as.numeric(Sets[row, 2]), rho = as.numeric(Sets[row, 3]), delt = as.numeric(Sets[row, 4]), c = as.numeric(Sets[row, 5]), lam = as.numeric(Sets[row, 6]), n = 1500)
  
  out <- HIV(pars=pars)
  
  DataT <- as.matrix(cbind(out[,c(1,2)], sd=0.45))
  
  HIVcost <- function (pars) 
  {
    out <- HIV(pars)
    cost <- modCost(model = out, obs = DataLogV, err = "sd")
    return(modCost(model = out, obs = DataT, err = "sd", cost = cost))
  }
  
  HIVcost2 <- function(lpars)
    HIVcost(c(exp(lpars), n = 1500))
  
  Pars <- log(pars[1:5] * 2)
  Fit <- modFit(f = HIVcost2, p = Pars)
  exp(coef(Fit))
  
  final <- HIV(pars = c(exp(coef(Fit)), n = 1500))
  
  model <- final[-c(3, 5, 7),]
  
  for (point in c(2,4,6))
  {
    for (error in c(7,8,9))
    {
      Sets[row, error] <- Sim[point, 2]-model[point, 5]
    }
  }
}

par(mfrow=c(1,2))
plot(Sim$time, Sim$logV, xlab = "time", ylab = "logV",log = "y", main = "Data Points")
lines(final$time, final$logV)
plot(out$time, out$logV, xlab = "time", ylab = "logV", log = "y", main = "Model Points")
