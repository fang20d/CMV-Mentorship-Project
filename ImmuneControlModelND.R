M46346 <- read_excel("ARV01 Plasma VL_RG11202018_in tabs.xlsx", sheet = 6)
M46346 <- M46346[-c(1:4),c(3,7)]
M46346 <- M46346[rowSums(is.na(M46346)) == 0,]
M46346 <- M46346[-c(10:27),]
M46346 <- cbind(M46346, sd = 0.45)
colnames(M46346) <- c("time", "logV", "sd")
s = 1
M46346$logV <- M46346$logV/s
M46346$logV <- log(M46346$logV)

DataLogV <- as.matrix(M46346)

dv0 <- as.numeric(M46346[2,2])-as.numeric(M46346[1,2])

HIV <- function (pars, V_0 = 10, dV_0 = dv0, T_0 = 100, E_0 = 20) 
{
  derivs <- function(time, y, pars) 
  {
    with (as.list(c(pars, y)), 
          {
            dT <- 1 -  T -  T * V
            dI <- a * T * V - I * E
            dV <- b * I - c * V
            dE <- I * E - d * E
            return(list(c(dT, dI, dV, dE), logV = log(V)))
          })
  }
  
  # initial conditions
  I_0 <- with(as.list(pars), (dV_0 + c * V_0) / (b))
  y <- c(T = T_0, I = I_0, V = V_0, E = E_0)
  
  times <- c(seq(0, 8, 0.1))
  out <- ode(y = y, parms = pars, times = times, func = derivs)
  
  as.data.frame(out)
}

# pars <- c(a = 30724.85, b = 0.3978971939, c = 12.3847, d = 0.837465)
# 2: pars <- c(a = 3074.85, b = 0.3978971939, c = 2.3847, d = 0.837465)
# 3: pars <- c(a = 6149.7000000, b = 0.7957944, c = 4.7694000, d = 1.6749300)
pars <- c(a = 1.230390e+05, b = 3.567974e-04, c = 4.284548e-04, d = 1.227347e+02 )
# pars <- c(a = 22324.551422, b = 3.203226, c = 19.419577, d = 6.617040)
# pars <- c(a = 22960.94514, b = 4017.72432, c = 72.56745, d = 50.72646)
# 7: pars <- c(a = 79693.94995, b = 12.44663, c = 85.38133, d = 25.62006)
# 8: pars <- c(a = 148852.70684, b = 26.96342, c = 188.52497, d = 45.16352 )
# 9: pars <- c(a = 3.920012e+4, b = 1.730743e+13, c = 2000, d = 2.255056e+04 )

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

Pars <- log(pars[1:4] * 2)
Fit <- modFit(f = HIVcost2, p = Pars)
exp(coef(Fit))

final <- HIV(pars = c(exp(coef(Fit))))

par(mfrow = c(1, 1))
plot(M46346$time, M46346$logV, xlab = "time", ylab = "logV",log = "y", main = "Data Points vs. Model Points", col = "red")
lines(final$time, final$logV)
legend("bottomright", c("data", "fitted"),
       lty = c(NA,1), pch = c(1, NA), col = c("red", "black"))

