using Plots
using PyPlot

dt = 0.001 #pas de temps
N = 100
dx = 0.001 #pas d'espace
Nx = convert(Int, 1 / dx) + 1
println("Nx = ", Nx)
T = [i * dt for i = 1:N]
println("t_0 = ", T[1])
println("t_n = ", T[N])

#parametres
delta = 50
d = 100
a = 0.2
b = 1.3

# ---- RESOLUTION ----

#concentration = (U, V)
#U[n][x] concentration temps n, en x
U = zeros(N, Nx)
V = zeros(N, Nx)
#conditions initiales
for j = 1:N
  U[1, j] = 2
  V[1, j] = 2
end

function fu(u, v)
  return (a - u + v * u^2)
end

function fv(u, v)
  return (b - v * u^2)
end

for n = 1:(N-1)
  for i = 2:(Nx-1)
    U[n+1, i] = U[n, i] + (dt / dx^2) * (U[n, i+1] - 2 * U[n, i] + U[n, i-1]) +
                dt * delta * fu(U[n, i], V[n, i])
    V[n+1, i] = V[n, i] +
                d * (dt / dx^2) * (V[n, i+1] - 2 * V[n, i] + V[n, i-1]) +
                dt * delta * fv(U[n, i], V[n, i])
  end
end



# ----- AFFICHAGE + SAUVEGARDE FIGURES -----

# ----- EVOLUTION TEMPORELLE ------

gr()

X = zeros(Int, 4)
X[1] = 1
X[2] = 2
X[3] = 3
X[4] = 4


ux1 = zeros(N)
ux2 = zeros(N)
ux3 = zeros(N)
ux4 = zeros(N)

for t = 1:N
  ux1[t] = U[t, X[1]]
end

for t = 1:N
  ux2[t] = U[t, X[2]]
end

for t = 1:N
  ux3[t] = U[t, X[3]]
end

for t = 1:N
  ux4[t] = U[t, X[4]]
end

p1 = Plots.plot(T, ux1, label = string("x = ", X[1]))
p2 = Plots.plot(T, ux2, label = string("x = ", X[2]))
p3 = Plots.plot(T, ux3, label = string("x = ", X[3]))
p4 = Plots.plot(T, ux4, label = string("x = ", X[4]))

Plots.plot(p1, p2, p3, p4, layout = (2, 2))


#sauvegarde figures
"""
ux = zeros(N)

for x = 995:1000
  for t = 1:N
    ux[t] = U[t,x]
  end
  plt.figure()
  plt.plot(T,ux,label = string("x = ", x))
  xlabel("temps")
  ylabel("concentration")
  plt.savefig(string("solution_concentration_x=",x), dpi = 300)

end
"""

# ----- EVOLUTION SPATIALE -----

X = [x for x = 20:Nx-10]
t = 20
println("t = ", T[t])
ut = zeros(Nx-29)
for i = 1:(Nx-29)
  ut[i] = U[t,X[i]]
end
Plots.plot(X,ut)
