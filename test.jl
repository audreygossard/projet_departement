using PyPlot

dt = 0.001 #pas de temps
N = 100
dx = 0.001 #pas d'espace
Nx = convert(Int, 1 / dx) + 1

#T = [i * dt for i = 1:N]
#X = [i * dx for i = 1:Nx]
#U = zeros(N, Nx)

T = [1,2,3]
X = [1,2,3]
U = zeros((3,3))
U[1,1] = 1
U[1,2] = 2
U[2,2] = 2

fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")

ax.plot_surface(X, T, U, rstride=10, cstride=10)

plt.savefig("solution 3D", dpi = 300)




# ----- AFFICHAGE 3D -----

fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")

ax.plot_surface(T, X, U, rstride=10, cstride=10)
#ou alors avec plot_wireframe c'est peut etre mieux ?
plt.savefig("solution 3D", dpi = 300)
