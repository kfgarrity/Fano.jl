function ham(freq, npts,  coupling, E1, gamma1, E2, gamma2, E3, gamma3, gamma_mag, fano_param) 

    N = npts + 3
    H = zeros( N, N)

    fmin = minimum(freq)
    fmax = maximum(freq)

    dvals = collect(fmin: (fmax-fmin)/(npts-1) : (fmax+1e-10))

#    println(typeof(dvals))
#    println(typeof(E1))
    v = [E1; E2; E3; dvals]
#    println(typeof(v))


    H[:,:] = diagm(v)

    #    f = fano_pos(dvals, fano_param[1], fano_param[2], fano_param[3])
    f = fano_other(dvals, fano_param[1], fano_param[2], fano_param[3])


#    H[1,4:end] = coupling * sqrt.(f)
#    H[4:end, 1] = coupling * sqrt.(f)

#    H[2,4:end] = coupling * sqrt.(f)
#    H[4:end,2] = coupling * sqrt.(f)

#    H[3,4:end] = coupling * sqrt.(f)
#    H[4:end,3] = coupling * sqrt.(f)

    H[1,4:end] = coupling * (f).^0.5 / sqrt(npts)
    H[4:end, 1] = coupling * (f).^0.5 / sqrt(npts)

    H[2,4:end] = coupling * (f).^0.5 / sqrt(npts)
    H[4:end,2] = coupling * (f).^0.5 / sqrt(npts)

    H[3,4:end] = coupling * (f).^0.5 / sqrt(npts)
    H[4:end,3] = coupling * (f).^0.5 / sqrt(npts)
    
    
    vals, vects = eigen(H)
    
    t1 = real(vects[1,:] .* conj(vects[1,:]))
    t2 = real(vects[2,:] .* conj(vects[2,:]))
    t3 = real(vects[3,:] .* conj(vects[3,:]))

    vc = zeros(N-3, N)
    for i = 1:N
        vc[:,i] = real((vects[4:end,i] .* conj(vects[4:end,i]))[:] .* f[:] )
    end

    vct = sum(vc, dims=1) * 200 / npts
    
    Alocal1 = zeros(length(freq))
    Alocal2 = zeros(length(freq))
    Alocal3 = zeros(length(freq))
    Adisp = zeros(length(freq))


    for (c,w) in enumerate(freq)
        for i = 1:N
            Alocal1[c] += gamma1^2 * t1[i] / ((vals[i] - w)^2 + gamma1^2)
            Alocal2[c] += gamma2^2 * t2[i] / ((vals[i] - w)^2 + gamma2^2)
            Alocal3[c] += gamma3^2 * t3[i] / ((vals[i] - w)^2 + gamma3^2)
            
            Adisp[c] += gamma_mag^2 * vct[i] / ((vals[i] - w)^2 + gamma_mag^2)
        end
    end

    return Alocal1, Alocal2, Alocal3, Adisp
    
end
    
