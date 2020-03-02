using PyPlot

# Données du problème
a = 0.2
b = 1.3
delta = 175
d = 35

Nx = 30
dt = 0.001
Nt = 5000

stock = zeros(Nt, 2 * Nx)  # tableau pour stocker tous les c

c = zeros(2 * Nx)    #variable pour implicite

r = randn(2*Nx)
o = ones(Nx)
c[1:Nx] = (a+b)*o +(r[1:Nx]*2-o)*10^-4
c[Nx+1:end] = b/((a+b)^2)*o +(r[Nx+1:end]*2-o)*10^-4
stock[1,1:Nx] = c[1:Nx]
stock[1,1+Nx:end] = c[Nx+1:end]

c0 = copy(c)

D = zeros(2 * Nx, 2 * Nx)  #matrice de diffusion
A = zeros(2 * Nx, 2 * Nx)  #matrice de l'opérateur laplacien
I = zeros(2 * Nx, 2 * Nx)   #matrice identié

#Remplissage des matrices A et D

for i = 1:2*Nx
    D[i, i] = 1
    I[i,i] = 1
    if (i >= Nx)
        D[i, i] = d
    end
end

for i = 1:2*Nx
    A[i, i] = -2
    if (i % Nx != 0)
        A[i, i+1] = 1
        A[i+1, i] = 1
    end
end

#conditions cycliques
A[1,Nx]=1
A[Nx,1]=1
A[Nx+1,2*Nx]=1
A[2*Nx,Nx+1]=1

A = Nx^2 * A

M=D*A

# comment faire pour le carre? c[1:Nx]*c[1:Nx]?????
# function f(c)
#     L = zeros(2*Nx)
#     L[1:Nx] = a*o - c[1:Nx] + c[1+Nx:end] * c[1:Nx]*c[1:Nx]
#     L[1+Nx:end] = b*o - c[1+Nx:end] * c[1:Nx]*c[1:Nx]
#     return L
# end

function f(c)
    L = zeros(2 * Nx)
    for x = 1:Nx
        L[x] = a - c[x] + c[x+Nx] * c[x]^2
        L[x+Nx] = b - c[x+Nx] * c[x]^2
    end
    return L
end


B = inv(I-dt*M)
for t = 1:Nt-1   #Incrémentation de c
    E = B*(c+dt*delta*f(c))
    stock[t+1,:] = E[:]
    c[:] = stock[t,:]
end


# ---- Affichage 3D -----

#PyPlot.use("tkagg")
figure()

dx = 1/Nx
X = [i*dx for i = 1:Nx]
T = [i*dt for i = 1:Nt]


u = copy(stock[:,1:Nx])
v = copy(stock[:,Nx+1:end])


surf(T, X, u', rstride=10, cstride=10)
surf(T, X, v', rstride=10, cstride=10)
xlabel("temps")
ylabel("espace")
zlabel("concentrations")
plt.show()
plt.savefig("affichage_3D_u", dpi = 300)
