library(PhaseTypeGenetics)

lambda = function(alpha, x){
  gamma(2)/(gamma(2-alpha)-gamma(alpha))*(x^(1-alpha))*((1-x)^alpha)
}

fx = function(x){
  x^(k-2) * (1-x)^(m-k) * (gamma(2)/(gamma(2-alpha)-gamma(alpha))*(x^(1-alpha))*((1-x)^alpha))
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
diag(Vecrate)=-apply(Vecrate,1,sum)
print(Vecrate)

Exrate = Vecrate[,n]
Vecrate = Vecrate[-n,-n]

a = contphasetype(c(1,0,0,0),Vecrate)
p=seq(0,8,0.01)
b=dphasetype(a,p)
plot(b)

y=NULL
m=8
for (k in 2:m){
  for (i in seq(0,1000)){
    y[i] = f2x(i/1000)
  }
  plot(y)
}
