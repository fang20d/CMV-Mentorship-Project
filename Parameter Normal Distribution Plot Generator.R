Ilam_mean <- 0.4211822/7
Alam_mean <- exp(-0.68)
Ilam_sd <- 0.06
Alam_sd <- 0.07
Ilower_bound <- Ilam_mean - Ilam_sd
Iupper_bound <- Ilam_mean + Ilam_sd
Alower_bound <- Alam_mean - Alam_sd
Aupper_bound <- Alam_mean + Alam_sd

x <- seq(-4, 4, length = 1000) * Ilam_sd + Ilam_mean
y <- dnorm(x, Ilam_mean, Ilam_sd)

g <- seq(-4, 4, length = 1000) * Alam_sd + Alam_mean
h <- dnorm(g, Alam_mean, Alam_sd)

plot(x, y, type="n", xlab = "1/day", ylab = "", main = "Distribution of Rho for Infants vs. Adults", xlim = c(0,1))
lines(x, y, col = "red")
points(g,h, type = "n")
lines(g, h, col = "blue")
legend("topright", c("infant", "adult"),
       lty = c(1,1), pch = c(NA, NA), col = c("red", "blue"))
