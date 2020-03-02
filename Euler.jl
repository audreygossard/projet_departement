using PyPlot

# Données du problème
a = 0.2
b = 1.3
delta = 175
d = 35

Nx = 30
dt = 0.001
Nt = 5000

stock1 = zeros(Nt, 2 * Nx)  # tableau pour stocker tous les c
stock2 = zeros(Nt, 2 * Nx)

c1 = zeros(2 * Nx)    #variable pour explicite
c2 = zeros(2 * Nx)
X=0:1/Nx:1-1/Nx


for x = 1:Nx  # initialisation de c
    c1[x] = a+b +(rand()*2-1)*10^-4
    c1[x+Nx] = b/((a+b)^2) +(rand()*2-1)*10^-4
    c2[x]=c1[x]
    c2[x+Nx]=c1[x+Nx]
    stock1[1,x]=c1[x]
    stock1[1,x+Nx]=c1[x+Nx]
    stock2[1,x]=c1[x]
    stock2[1,x+Nx]=c1[x+Nx]
end

u01 = copy(c1)
u02 = copy(c2)

#c2=c1         #variable pour implicite
#cmod1 = c1    #variable auxiliaire pour l'incrémentation
#comd2 = c2

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

for t = 1:Nt-1   #Incrémentation de c
    E = inv(I-dt*M)*(c2+dt*delta*f(c2))
    for x = 1:2*Nx
        #stock1[t+1,x] = c1[x] + dt * ((M * c1)[x] + delta*f(c1)[x]) #explicite
        stock1[t+1,x]=0
        stock2[t+1,x] = E[x] #implicite
        c1[x]=stock1[t,x]
        c2[x]=stock2[t,x]
    end
end

#ploting
u1=c1[1:Nx]
u2=c2[1:Nx]
v1=c1[Nx+1:2*Nx]
v2=c2[Nx+1:2*Nx]

print("end")


#p1=Plots.plot(X,u1,label="u schéma explicite")
#p2=Plots.plot(X,u2,label="u schéma implicite")
#p3=Plots.plot(X,v1,label="v schéma explicite")
#p4=Plots.plot(X,v2,label="v schéma implicite")

#Plots.plot(p1,p2,p3,p4,layout=(2,2))



# ---- Affichage 3D -----

#PyPlot.use("tkagg")
figure()

dx = 1/Nx
X = [i*dx for i = 1:Nx]
T = [i*dt for i = 1:Nt]


u = zeros((Nt,Nx))
v = zeros((Nt,Nx))
for t = 1:Nt
    for x = 1:Nx
        u[t,x] = stock2[t,x]
        v[t,x] = stock2[t,x+Nx]
    end
end


#Plots.plot(T,X,u) # fonctionne pas??

surf(T, X, u', rstride=10, cstride=10)
surf(T, X, v', rstride=10, cstride=10)
plt.show()
plt.savefig("affichage_3D_u", dpi = 300)

#ax_v.plot_wireframe(T, X, u, rstride=10, cstride=10)
#plt.savefig("affichage_3D_v", dpi = 300)
