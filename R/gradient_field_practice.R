### Attempt 1

```{r}
df <- expand.grid(x=seq(-2, 2, .1), y=seq(-2, 2, .1)); k=2; b=4
df$z <- with(df, (y^2)-(x^2) )
library(raster)
library(rasterVis)
r <- rasterFromXYZ(df)
projection(r) <- CRS("+proj=longlat +datum=WGS84")
vectorplot(r, par.settings=RdBuTheme)
```

```{r}
df <- expand.grid(x=seq(-2, 2, .1), y=seq(-2, 2, .1)); k=2; b=4
df$z <- with(df, (y^2)-(x^2) )
library(raster)
library(rasterVis)
r <- rasterFromXYZ(df)
projection(r) <- CRS("+proj=longlat +datum=WGS84")
vectorplot(r, par.settings=RdBuTheme)
```

```{r}
# var 1 to 3 are given as inputs and do not change.
var1 = 3
var2 = 1
var3 = 0.02

# x and y are additional variables, but they are plotted within a range 
# on the chart, xlim = c(0, 24), ylim = c(0, 0.05)

# Given var 1 to 3, and the ranges x and y, plot the result (z) as a heat map and/or contour
result <- function(var1, var2, var3,x, y){
  var2 * (1 + (x/12) * var3) * (1 + -var1 * y)
}

# sample output at point x = 12, y = 0.02, and equals 0.9588 in this example.
result(var1 = 3, var2 = 1,var3 = 0.02,
       x=12, 
       y=0.02)

```

```{r}
result <- function(var1, var2, var3,x, y,d){
  ydata=numeric(length(seq(y[1],y[2],d)))
  mat=matrix(NA,length(seq(y[1],y[2],d)),length(seq(x[1],x[2],d)))
  yy=seq(y[1],y[2],d)
  xx=seq(x[1],x[2],d)
  for(i in 1:length(xx)){
    for(n in 1:length(yy)){
      ydata[n]=var2 * (1 + (xx[i]/12) * var3) * (1 + -var1 * yy[n])
    }
    mat[,i]=ydata
  }
  return(mat)}
```

```{r}
rotate <- function(x) t(apply(x, 2, rev))
```

```{r}
mydata=result(var1 = 3, var2 = 1,var3 = 0.02,x=c(-1200,1200),y=c(-3005,2008),d=10)
image((mydata))
par(new=TRUE)
contour(mydata, add = TRUE)
```

\newpage

Let the posterior mean function for a particular GP, over a 2D input space

```{r, include=F}
x=c(-1200,1200)
```

$$x\in[-1200,1200]$$
  
  ```{r, eval=F}
var1 = 3; var2 = 1; var3 = 0.02
y=c(-3005,2008)
d=10

ydata=numeric(length(seq(y[1],y[2],d)))
mat=matrix(NA,length(seq(y[1],y[2],d)),length(seq(x[1],x[2],d)))

for(i in 1:length(xx)){
  for(n in 1:length(yy)){
    ydata[n]=var2 * (1 + (xx[i]/12) * var3) * (1 + -var1 * yy[n])
  }
  mat[,i]=ydata
}
```

be represented by:
  
  ```{r}
result <- function(var1, var2, var3,x, y,d){
  ydata=numeric(length(seq(y[1],y[2],d)))
  mat=matrix(NA,length(seq(y[1],y[2],d)),length(seq(x[1],x[2],d)))
  yy=seq(y[1],y[2],d)
  xx=seq(x[1],x[2],d)
  for(i in 1:length(xx)){
    for(n in 1:length(yy)){
      ydata[n]=var2 * (1 + (xx[i]/12) * var3) * (1 + -var1 * yy[n])
    }
    mat[,i]=ydata
  }
  return(mat)}
```



```{r}
mydata=result(var1 = 3, var2 = 1,var3 = 0.02,x=c(-1200,1200),y=c(-3005,2008),d=10)
dim(mydata)
image((mydata))
par(new=TRUE)
contour(mydata)

yy=seq(y[1],y[2],10*d)
xx=seq(x[1],x[2],10*d)
df2 <- data.frame(rep(xx,length(yy)), sort(rep(yy,length(xx))), rep(NA, length(xx)*length(yy)))
colnames(df2) <- c("x1", "x2", "y")
for (j in 1:length(yy)){
  for (i in 1:length(xx)){
    df2[(i-1)*length(yy)+j,3] <- var2 * (1 + (xx[i]/12) * var3) * (1 + -var1 * yy[j])
  }
}
df2
tail(df2)
r <- rasterFromXYZ(df2)
projection(r) <- CRS("+proj=longlat +datum=WGS84")
vectorplot(r, par.settings=RdBuTheme)

df <- expand.grid(x=seq(-2, 2, .1), y=seq(-2, 2, .1)); k=2; b=4
df$z <- with(df, (y^2)-(x^2) )
library(raster)
library(rasterVis)
r <- rasterFromXYZ(df)
projection(r) <- CRS("+proj=longlat +datum=WGS84")
vectorplot(r, par.settings=RdBuTheme)
```

```{r}
mf2 <- expression(
  var2 * (1 + (x1/12) * var3) * (1 + -var1 * x2)
)
D(mf2, 'x1')
D(mf2, 'x2')
xgrid <- seq
```

### Attempt 2

df <- expand.grid(x1=seq(-2, 2, .1), x2=seq(-2, 2, .1)); k=2; b=4
df$z <- with(df, (x2^2)-(x1^2) )
library(raster)
library(rasterVis)
r <- rasterFromXYZ(df)
projection(r) <- CRS("+datum=WGS84")
vectorplot(r, par.settings=RdBuTheme)

mf2 <- expression(
  (x2^2)-(x1^2)
)
D(mf2, "x1")
D(mf2, "x2")
Dmf2_x1 <- function(input){
                    -(2 * input)     
}
Dx1 <- lapply(x1, Dmf2_x1)
Dmf2_x2 <- function(input){
  2 * input
}
Dx2 <- lapply(x2, Dmf2_x2)

library(gcookbook) # For the isabel data set       
isabel
# Need to load grid for arrow() function
library(grid)

# Make the plot with the subset, and use an arrowhead 0.1 cm long
ggplot(isabel, aes(x = x, y = y)) +
  geom_segment(aes(xend = x+vx/50, yend = y+vy/50),
               arrow = arrow(length = unit(0.1, "cm")), size = 0.25)

df22 <- expand.grid(x1=seq(-2, 2, 1), x2=seq(-2, 2, 1))
df22$z <- with(df22, )
df$z <- with(df, (x2^2)-(x1^2) )

DF <- data.frame(x1 = (rep(NA,41)),
                 x2 = (rep(NA,41)),
                 Dx1 = (rep(NA,41)),
                 Dx2 = (rep(NA,41)))
DF$x1 <- x1
DF <- data.frame(x1 = x1, x2 = x2, Dx1 = Dx1, Dx2 = Dx2)
head(DF)

ggplot(isabel, aes(x = x1, y = x2)) +
  geom_segment(aes(xend = x1+Dx1/50, yend = x2+Dx2/50))

#########################

mf2 <- expression(
  x1^2*x2
)
D(mf2, "x1")
D(mf2, "x2")
Dmf2_x1 <- function(x1, x2){
  2 * x1 * x2    
}
Dmf2_x2 <- function(x1, x2){
  x1^2   
}

df44 <- expand.grid(x1=seq(-5, 5, 1), x2=seq(-5, 5, 1))
df44$vx1 <- with(df44, Dmf2_x1(x1, x2))
df44$vx2 <- with(df44, Dmf2_x2(x1, x2))
# df33$vx1 <- c(2, 2, 2, 1, 1, 1, 0.5, 0.5, 0.5)
# df33$vx2 <- c(0, 0, 0, 1, 1, 1, 1.5, 1.5, 1.5)
library(grid)
arrow_scale <- 50
ggplot(df44, aes(x = x1, y = x2)) +
  geom_segment(aes(xend = x1+vx1/arrow_scale, yend = x2+vx2/arrow_scale),
               arrow = arrow(length = unit(0.1, "cm")), size = 0.25)


Dx1 <- lapply(x1, Dmf2_x1)
Dmf2_x2 <- function(input){
  2 * input
}
Dx2 <- lapply(x2, Dmf2_x2)


df33 <- expand.grid(x1=seq(-1, 1, 1), x2=seq(-1, 1, 1))
# df22$z <- with(df22, )
df33$vx1 <- c(2, 2, 2, 1, 1, 1, 0.5, 0.5, 0.5)
df33$vx2 <- c(0, 0, 0, 1, 1, 1, 1.5, 1.5, 1.5)
library(grid)
arrow_scale <- 5
ggplot(df33, aes(x = x1, y = x2)) +
  geom_segment(aes(xend = x1+vx1/arrow_scale, yend = x2+vx2/arrow_scale),
                   arrow = arrow(length = unit(0.1, "cm")), size = 0.25)
