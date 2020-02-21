using Plots
dt =0.01 #pas de temps
N =10
dx = 0.01 #pas d'espace
Nx = convert(Int, 1/dx)+1
T = [i*dt for i=1:N]
#parametres
delta = 50
d = 100
a = 0.2
b = 1.3
#concentration = (U, V)
#U[n][x] concentration temps n, en x
U = zeros(N, Nx)
V = zeros(N, Nx)
#conditions initiales
for j=1:N
  U[1,j] = 2
  V[1,j] = 2
end

function fu(u,v)
  return (a-u+v*u^2)
end

function fv(u, v)
  return (b-v*u^2)
end

for n=1:(N-1)
  for i=2:(Nx-1)
    U[n+1,i] = U[n,i] + (dt/dx^2)*(U[n,i+1]-2*U[n,i]+U[n,i-1])+dt*delta*fu(U[n,i], V[n,i])
    V[n+1,i] = V[n,i] + d*(dt/dx^2)*(V[n,i+1]-2*V[n,i]+V[n,i-1])+dt*delta*fv(U[n,i],V[n,i])
  end
end
gr()
for x=1:Nx
  ux = zeros(N)
  for t=1:N
    ux[t] = U[t,x]
  end
  plot(T, ux)
end
