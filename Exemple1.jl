using Plots

h=0.01
N=10000
L=h*N

mu  = 5
mu1 = 5
mu2 = 5
mu3 = 5
# pour l'oscillateur
#mu->se comporte un facteur d'amortissement

T=zeros(N)
U1=zeros(N)
U2=zeros(N)
U3=zeros(N)
U4=zeros(N)

X=zeros(N)
DX=zeros(N)

X1=zeros(N)
DX1=zeros(N)

X2=zeros(N)
DX2=zeros(N)

X3=zeros(N)
DX3=zeros(N)

U1[1]=0.01    #solution pour mu = -1 donc stable
U2[1]=0.01   #solution pour mu = +2 donc instable
U3[1]=1.01   #solution pour mu = +2 donc stable
U4[1]=-1.01  #solution pour mu = -(+ pour f2)1 donc instable

X[1]=1.
DX[1]=1.

X1[1]=-1.
DX1[1]=1.

X2[1]=-1.
DX2[1]=-1.

X3[1]=1.
DX3[1]=-1.

T[1]=0.


function f1(u,mu)
    return u*(mu-u)
end

function f2(u,mu)
    return u*(mu-u^2)
end

function f3(x,dx,mu)
    return (dx,(mu-x^2)*dx-x)
end

for k = 1:(N-1)
    T[k+1]=T[k]+h
    # U1[k+1]=U1[k]+h*f1(U1[k],mu1)
    # U2[k+1]=U2[k]+h*f1(U2[k],mu2)
    # U3[k+1]=U3[k]+h*f1(U3[k],mu2)
    # U4[k+1]=U4[k]+h*f1(U4[k],mu1)

    X[k+1]=X[k]+h*f3(X[k],DX[k],mu)[1]
    DX[k+1]=DX[k]+h*f3(X[k+1],DX[k],mu)[2]

    X1[k+1]=X1[k]+h*f3(X1[k],DX1[k],mu1)[1]
    DX1[k+1]=DX1[k]+h*f3(X1[k+1],DX1[k],mu1)[2]

    X2[k+1]=X2[k]+h*f3(X2[k],DX2[k],mu2)[1]
    DX2[k+1]=DX2[k]+h*f3(X2[k+1],DX2[k],mu2)[2]

    X3[k+1]=X3[k]+h*f3(X3[k],DX3[k],mu3)[1]
    DX3[k+1]=DX3[k]+h*f3(X3[k+1],DX3[k],mu3)[2]

end

gr()

# p1=plot(T,U1,label = "u1")
# p2=plot(T,U2,label = "u2")
# p3=plot(T,U3,label = "u3")
# p4=plot(T,U4,label = "u4")

p1=plot(X,DX,label = "dx/dt(x)")
p2=plot(X1,DX1,label = "dx/dt(x)")
p3=plot(X2,DX2,label = "dx/dt(x)")
p4=plot(X3,DX3,label = "dx/dt(x)")
plot(p1,p2,p3,p4, layout =(2,2))
