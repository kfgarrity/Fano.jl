function load()

    VH = zeros(Float64, 14, 1200, 2)
    
    for (c,t) = enumerate([10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90])
        println("VH $t")
        f = readdlm("MnPSe3_Tdep_VH/MnPSe3_cryo_$(t).0K_515nm_400uW_slitsat10_60sexp_10acc_VH_12cm-1Limit.txt")
        a,b = size(f)
        VH[c,1:a,:] = Float64.(f)
    end

    VV = zeros(Float64, 14, 1200, 2)
    
    for (c,t) = enumerate([10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90])
        println("VV $t")
        f=missing
        try
            f = readdlm("MnPSe3_Tdep_VV/MnPSe3_cryo_$(t).0K_515nm_400uW_slitsat10_60sexp_10acc_VV_12cm-1Limit.txt")
        catch

                    
            f = readdlm("MnPSe3_Tdep_VV/MnPSe3_cryo_$(t).0K_515nm_400uW_slitsat10_60sexp_30acc_VV_12cm-1Limit.txt")
        end
        a,b = size(f)
        VV[c,1:a,:] = Float64.(f)
    end

    return VH, VV

end
