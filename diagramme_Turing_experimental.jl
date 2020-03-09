using PyPlot

# Données du problème
a = 0.2
b = 1.3
Nx = 30
dx = 1/Nx
Nt = 3000
dt = 0.001

X = [i*dx for i = 1:Nx]
T = [i*dt for i = 1:Nt]

A = zeros(2 * Nx, 2 * Nx)  #matrice de l'opérateur laplacien
I = zeros(2 * Nx, 2 * Nx)   #matrice identié

#Remplissage de la matrice A + conditions cycliques
for i = 1:2*Nx
    A[i, i] = -2
    if (i % Nx != 0)
        A[i, i+1] = 1
        A[i+1, i] = 1
    end
end

A[1,Nx]=1
A[Nx,1]=1
A[Nx+1,2*Nx]=1
A[2*Nx,Nx+1]=1

A = Nx^2 * A


function f(c)
    L = zeros(2 * Nx)
    for x = 1:Nx
        L[x] = a - c[x] + c[x+Nx] * c[x]^2
        L[x+Nx] = b - c[x+Nx] * c[x]^2
    end
    return L
end


function conditions_initiales_alea() #pour ne pas avoir a en regenerer a chaque fois
    #renvoie c et stock avec conditions initiales aux premiers indices
    stock = zeros(Nt, 2 * Nx)
    c = zeros(2 * Nx)
    r = randn(2*Nx)
    o = ones(Nx)
    c[1:Nx] = (a+b)*o +(r[1:Nx]*2-o)*10^-4
    c[Nx+1:end] = b/((a+b)^2)*o +(r[Nx+1:end]*2-o)*10^-4
    stock[1,1:Nx] = c[1:Nx]
    stock[1,1+Nx:end] = c[Nx+1:end]
    return (c,stock)
end


function stabilite(u_final, v_final, seuil) #renvoie un booleen qui indique si la solution est stable ou non
    stable = true
    x = 1
    while (stable && x<=Nx)
        stable = (abs(u_final[x]-a-b)<seuil && abs(v_final[x]-b/(a+b)^2)<seuil)
        x += 1
    end
    return stable
end


function resout_euler(delta, d, c0, stock0)
    c = copy(c0)
    stock = copy(stock0)
    D = zeros(2 * Nx, 2 * Nx)  #matrice de diffusion
    #Remplissage de la matrice de diffusion
    for i = 1:2*Nx
        D[i, i] = 1
        I[i,i] = 1
        if (i >= Nx)
            D[i, i] = d
        end
    end

    M=D*A
    B = inv(I-dt*M)

    for t = 1:Nt-1   #Incrémentation de c
        E = B*(c+dt*delta*f(c))
        stock[t+1,:] = E[:]
        c[:] = stock[t,:]
    end

    u = copy(stock[:,1:Nx])
    v = copy(stock[:,Nx+1:end])
    return (u,v)
end

c0,stock0 = conditions_initiales_alea()


tab_d = zeros(100)
tab_delta = zeros(100)
d_stables = []
delta_stables = []

for i=1:100
    tab_d[i] = i
    tab_delta[i] = i*4
end

for d in tab_d
    for delta in tab_delta
        u,v = resout_euler(delta, d, c0, stock0)
        if stabilite(u[Nt,:], v[Nt,:], 0.05)
            append!(d_stables, d)
            append!(delta_stables, delta)
        end
    end
end

plt.figure()
plt.plot(d_stables,delta_stables, marker = "s", linestyle = "None", label = "solution stable")
xlabel("d")
ylabel("delta")
legend()
plt.savefig("diagrammes_turing_exp.png", dpi = 300)
