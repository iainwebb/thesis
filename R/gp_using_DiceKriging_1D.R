# ------------
# a 1D example
# ------------

library(DiceKriging)

# Number of observations
n <- 5

# Specify the inputs
control_inputs <- matrix(seq(0,2*pi,len=n+2), ncol = 1)[2:(n+1),]

# Real-life function
real_life_process <- function(x){
  sin(x)
}

# Observations
f <- real_life_process(control_inputs)
f <- as.matrix(f, ncol=1)

nf_data_points <- data.frame(x=control_inputs,
                             y=f)
colnames(nf_data_points) <- c("control_inputs", "f")


x <- control_inputs
y <- real_life_process(control_inputs)

formula <- y~0

theta <- 1
sigma <- 1 
# trend <- c(1.2,-0.4)

model1 <- km(~0, design=data.frame(x=x), response=data.frame(y=y), 
            covtype="gauss", 
            # coef.trend=trend,
            coef.cov=theta,
            coef.var=sigma^2
            )
model1

# model2 <- km(formula=formula, design=data.frame(x=x), response=data.frame(y=y),
#             covtype="gauss")
# 
# model2

model <- model1

# Bounds
expand_quantity <- 3
x_upper <- 2*pi + expand_quantity
x_lower <- 0 - expand_quantity

tmin <- x_lower; tmax <- x_upper
t <- seq(from=tmin, to=tmax, by=0.005)
color <- list(SK="black", UK="blue")

# Results with Universal Kriging formulae (mean and and 95% intervals)
p.UK <- predict(model, newdata=data.frame(x=t), type="UK")
plot(t, p.UK$mean, type="l", ylim=c(min(p.UK$lower95),max(p.UK$upper95)),
     xlab="x", ylab="y")
lines(t, p.UK$trend, col="violet", lty=2)
lines(t, p.UK$lower95, col=color$UK, lty=2)
lines(t, p.UK$upper95, col=color$UK, lty=2)
points(x, y, col="red", pch=19)
abline(h=0)

# Results with Simple Kriging (SK) formula. The difference between the width of
# SK and UK intervals are due to the estimation error of the trend parameters
# (but not to the range parameters, not taken into account in the UK formulae).
p.SK <- predict(model, newdata=data.frame(x=t), type="SK")
lines(t, p.SK$mean, type="l", ylim=c(-7,7), xlab="x", ylab="y")
lines(t, p.UK$trend, col="violet", lty=2) # Same as for UK
lines(t, p.SK$lower95, col=color$SK, lty=2)
lines(t, p.SK$upper95, col=color$SK, lty=2)
points(x, y, col="red", pch=19)
abline(h=0)
legend.text <- c("Universal Kriging (UK)", "Simple Kriging (SK)")
legend(x=tmin, y=max(p.UK$upper), legend=legend.text,
       text.col=c(color$UK, color$SK), col=c(color$UK, color$SK),
       lty=3, bg="white")