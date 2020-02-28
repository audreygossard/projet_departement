using Plots

# Données du problème
a = 0.2
b = 1.3
delta = 5
d = 1

Nx = 30
dt = 0.05
T = 5
Nt = 50 # T/dt

stock = zeros(Nt, 2 * Nx)  # tableau pour stocker tous les c

c1 = zeros(2 * Nx)    #variable pour explicite
X=0:1/Nx:1-1/Nx

for x = 1:Nx  # initialisation de c
    c1[x] = a+b+rand()*10^-4
    c1[x+Nx] = b/((a+b)^2)+rand()*10^-4
end

c2=c1         #variable pour implicite
cmod1 = c1    #variable auxiliaire pour l'incrémentation
comd2 = c2

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

function f(c)
    L = zeros(2 * Nx)
    for x = 1:Nx
        L[x] = a - c[x] + c[x+Nx] * c[x]^2
        L[x+Nx] = b - c[x+Nx] * c[x]^2
    end
    return L
end

for t = 1:Nt   #Incrémentation de c
    cmod1 = c1 + dt * (M * c1 + delta*f(c1)) #explicite
    cmod2 = inv(I-dt*M)*(c2+dt*f(c2)) #implicite
    for x =1:2*Nx
        c1[x]=cmod1[x]
        c2[x]=cmod2[x]
    end
    #println(c1[1:Nx])
end

#ploting
u=c2[1:Nx]
v=c1[1:Nx]

plotly()

p1=plot(X,u,ylims=(0.,2.))
p2=plot(X,v,ylims=(0.,2.))

plot(p1,p2,layout=(2,1))

print("end")
