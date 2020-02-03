using Plots

h=0.01
N=1000
L=h*N
mu1 = -1.
mu2 = 4.

T=zeros(N)
U1=zeros(N)
U2=zeros(N)
U3=zeros(N)
U4=zeros(N)

X=zeros(N)
DX=zeros(N)

#ui ::int8 = -1.
U1[1]=0.01    #solution pour mu = -1 donc stable
U2[1]=0.01   #solution pour mu = +4 donc instable
U3[1]=2.01   #solution pour mu = +4 donc stable
U4[1]=-2.01  #solution pour mu = +4 donc stable

T[1]=0.


function f(u,mu)
    return u*(mu-u^2)
end


for k = 1:(N-1)
    T[k+1]=T[k]+h
    U1[k+1]=U1[k]+h*f(U1[k],mu1)
    U2[k+1]=U2[k]+h*f(U2[k],mu2)
    U3[k+1]=U3[k]+h*f(U3[k],mu2)
    U4[k+1]=U4[k]+h*f(U4[k],mu2)


end

gr()

p1=plot(T,U1,label = "u1")
p2=plot(T,U2,label = "u2")
p3=plot(T,U3,label = "u3")
p4=plot(T,U4,label = "u4")
plot(p1,p2,p3,p4, layout =(2,2))
