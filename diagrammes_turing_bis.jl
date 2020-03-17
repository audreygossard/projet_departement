using PyPlot
using LinearAlgebra

a = 0.2
b = 1.3
d_crit = ((1+sqrt(2*b*(a+b)))*(a+b)/(a-b))^2

function stabilite(n,d,delta)
    B = zeros(2,2)
    pi2 = pi*pi
    n2 = n*n
    B[1,1] = -n2*pi2+delta*(b-a)/(a+b)
    B[1,2] = delta*(a+b)*(a+b)
    B[2,1] = -2*b*delta/(a+b)
    B[2,2] = -d*n2*pi2-delta*(a+b)*(a+b)
    return (det(B)>0)
end

couleurs = ["b","g","r","c","m","y"]

function diagramme_turing(d1,d2,d_nb,delta1,delta2,delta_nb,n1,n2)
    #fig,axs = subplots(1, n2-n1+1, sharey = true)
    figure()
    list_d = range(d1,stop=d2,length = d_nb)
    list_delta = range(delta1,stop=delta2, length = delta_nb)
    compt_label = n1 #pour savoir si on a deja mis le label pour n
    for n = n1:n2
        d_prec = -20
        for (j,delta) in enumerate(list_delta)
            i = 1
            stable = true
            while (stable && i<=d_nb)
                d = list_d[i]
                stable = stabilite(n,d,delta)
                if (!stable)
                    if (d_prec!=-20)
                        if (compt_label != n)
                            plt.plot([d_prec,d], [list_delta[j-1],delta], color = couleurs[n-n1+1], label = string("n = ", n))
                            compt_label += 1
                        else
                            plt.plot([d_prec,d], [list_delta[j-1],delta], color = couleurs[n-n1+1])
                        end
                    end
                    d_prec = d
                end
                i += 1
            end
        end
    end
    xlabel("d")
    ylabel("delta")
    plt.plot()
    legend()
    plt.savefig("diagrammes_turing_bis.png", dpi = 300)
    #show() #ne fonctionne pas! :'(
end


diagramme_turing(d_crit/100,d_crit*2,1000,0.00005,500,1000,0,5)
