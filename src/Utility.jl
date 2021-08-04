function write(mystr, FREQ, RR)

    for (c,t) = enumerate([10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90])

        len = length(FREQ[c])
    
        T = zeros(Float64, len, 2)

        T[:,1] = FREQ[c]
        T[:,2] = RR[c]
        
        f = writedlm(mystr*"_$t.csv", T)

    end

    
end


using Plots

function myplot(dat, poff = 1.0, n = 1:14)

    offset = 0.0

    plot()
    for (c,t) = zip(n, [10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90][n])
        freq = dat[c, :, 1]
        vals = dat[c, :, 2]

        a = findall(x->x < 1e-5 , freq )[1] - 1

        

        plot!(freq[1:a], (vals[1:a]) .+ offset, label=string(t))
        plot!(freq[1:a], (zeros(a)) .+ offset, label="", color="gray", line=:dot)

        offset += poff
    end
    xlims!(60, 180)
    ylims!(-3.0, 110.0/8.0+1.0)
    display(plot!())
end


function fano(x, E_res, G, q)

    return 1.0 .- (q * G / 2 .+ (x .- E_res)	 ).^2 ./ ( (G/2)^2 .+ (x .- E_res).^2)

end

function fano_pos(x, E_res, G, q)

    return max.(1.0 .- (q * G / 2 .+ (x .- E_res)	 ).^2 ./ ( (G/2)^2 .+ (x .- E_res).^2) , 0.0)

end

function lorentzian(x, E_res, gamma)

    return (gamma/2)^ 2 ./ ((gamma/2)^ 2 .+ (x .- E_res).^2)
    
end



function fano_other(x, E_res, G, q)

    eps = (E_res .- x)/G
    t = (q .+ eps).^2 ./ (1.0 .+ eps.^2) .* (1/(q^2))
    dx = x[2] - x[1]

    int = dx*sum(t)
    return  t / int
    
end

    
