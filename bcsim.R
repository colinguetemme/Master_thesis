library(PhaseTypeGenetics)

lambda = function(alpha, x){
  gamma(2)/(gamma(2-alpha)-gamma(alpha))*(x^(1-alpha))*((1-x)^alpha)
}

fx = function(x){
  x^(k-2) * (1-x)^(m-k) * (gamma(2)/(gamma(2-alpha)-gamma(alpha))*(x^(1-alpha))*((1-x)^alpha))
}

fx = function(m,k,alpha=1){
  lmk=beta(k-alpha,m-k+alpha)/beta(alpha,2-alpha)
  return(lmk)
  }
  
}
f2x = function(x){
  x^(k-2) * (1-x)^(m-k)
}

n=5
alpha = 1.2
Vecrate = diag(rep(0,n))
prob = 0

for (m in n:2){
  for (k in 2:m){
    prob = choose(m,k)*integrate(fx,0,1)$value[1]
    Vecrate[(n-m+1), (n-m+k)] = prob
  }
}

n=5
alpha = 1.2
Vecrate = diag(rep(0,n))
prob = 0

for (m in n:2){
  for (k in 2:m){
    prob = fx(m,k,alpha)
    Vecrate[(n-m+1), (n-m+k)] = prob
  }
}

diag(Vecrate)=-apply(Vecrate,1,sum)
print(Vecrate)

Exrate = Vecrate[,n]
Vecrate = Vecrate[-n,-n]

a = contphasetype(c(1,0,0,0),Vecrate)
p=seq(0,8,0.01)
b=dphasetype(a,p)
plot(b)

plot(NULL,xlim=c(0,100), ylim=c(0,0.1))
y=NULL
m=10
for (k in 2:((m/2)+1)){
  for (i in seq(0,100)){
    y[i] = f2x(i/100)
  }
  lines(y, col=k)
}

plot(NULL,xlim=c(0,100), ylim=c(0,10))
y=NULL
m=10
for (alpha in seq(1,1.99,0.1)){
  for (i in seq(0,100)){
    y[i] = lambda(alpha,i/100)
  }
  lines(y, col=(alpha-1)*10)
  Sys.sleep(0.5)
}

plot(NULL,xlim=c(0,100), ylim=c(0,0.1))
y=NULL
m=10
for (alpha in seq(1.9999999,1,length.out = 10)){
  for (i in seq(0,1,0.01)){
    y[i/0.01] = dbeta(i,2-alpha,alpha)/100
  }
  lines(y, col=(alpha-1)*10)
  Sys.sleep(0.5)
}

