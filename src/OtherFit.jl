function fit_fano(data,  E, G, q, m, background; t1=missing, t2=missing, doplot=false, ftype=:fano)

    freq = data[:,1]

    if !ismissing(t1)
        a1 = findall(x->x > t1 , freq )[1] - 1
        a1 = max(a1, 1)
    else
        a1 = 1
    end

    if !ismissing(t2)
        a2 = findall(x->x > t2 , freq )[1]  + 1
        a2 = min(a2, length(freq))
    else
        a2 = length(freq)
    end


    
    freq = data[a1:a2,1]
    val = data[a1:a2,2]

    
    function f(x)
        d2 = x[4] * fano(freq, x[1], x[2], x[3]) .+ x[5]
        return sum((d2 - val).^2)
    end

    function f_lor(x)
        d2 = x[3] * lorentzian(freq, x[1], x[2]) .+ x[4]
        return sum((d2 - val).^2)
    end

    

#    res = optimize(f)
    if ftype == :fano
        res = optimize(f, [E,G,q,m, background], LBFGS(); autodiff = :forward)

    else
        res = optimize(f_lor, [E,G,m, background], LBFGS(); autodiff = :forward)
    end

    x = Optim.minimizer(res)

    println(x)
    if doplot

        plot(freq, val, color="blue", line=:solid)
        if ftype == :fano
            display(plot!(freq, x[4].* fano(freq, x[1], x[2], x[3]) .+ x[5] , color="orange", line=:dash))
        else
            display(plot!(freq, x[3].* lorentzian(freq, x[1], x[2]) .+ x[4] .+ offset, color="orange", line=:dash, label=""))
        end
    end    

    return x, freq, val

end

function manyfit(data,  E; G=1.5, q=0.002, m=2.7, background=0.05, t1=missing, t2=missing, doplot=false, plot_offset = 1.0, ftype=:fano)

    N = size(data)[1]

    T = zeros(Float64, N, 5)

    offset = 0.0

    if doplot
        plot()
    end
    for n = 1:N

        x, freq, vals = fit_fano(data[n,:,:],  E, G, q, m, background; t1=t1, t2=t2, doplot=false, ftype=ftype)

        if ftype == :fano
            E, G, q, m, background  = x
            T[n,:] = [E, G, q, m, background]
        else
            E, G, m, background  = x
            T[n,:] = [E, G, 0.0, m, background]
        end
        
        if doplot
            plot!(freq, vals .+ offset, color="orange", line=:solid, label="")
            x = T[n,:]

            if ftype == :fano
                plot!(freq, x[4].* fano(freq, x[1], x[2], x[3]) .+ x[5] .+ offset, color="blue", line=:dash, label="")
            else
                plot!(freq, x[4].* lorentzian(freq, x[1], x[2]) .+ x[5] .+ offset, color="blue", line=:dash, label="")
            end
            
            offset += plot_offset
            
        end
    end
    if doplot
        display(plot!())
    end
        
    return T
        
end







#

function fit_ham(data, npts;  coupling=0.1 , E1=83, E2=108, E3=147.5, gamma_phonon=0.5, gamma_phonon2 = 0.5, gamma_mag = 0.5, f1=130, f2 = 5.0, f3 = 0.1, A1 = 5.0, A2 = 10.0, A3 = 100.0, Amag = 1.0, background = 0.0, t1=missing, t2=missing, doplot=false)

    println("npts ", npts)
    
    freq = data[:,1]

    if !ismissing(t1)
        a1 = findall(x->x > t1 , freq )[1] - 1
        a1 = max(a1, 1)
    else
        a1 = 1
    end

    if !ismissing(t2)
        a2 = findall(x->x > t2 , freq )[1]  + 1
        a2 = min(a2, length(freq))
    else
        a2 = length(freq)
    end

    freq = data[a1:a2,1]
    val = data[a1:a2,2]

    weight = ones(length(freq))

    if t2 > 140.0
        wcut = findall(x->x > 140 , freq )[1] - 1
        weight[wcut:end] .= 0.01
    end
    

    

    
    function f(x)

        coupling_, E1_, E2_ , E3_, gamma_phonon_, gamma_phonon2_, f1_,f2_,f3_,A1_,A2_, A3_, Amag_, background_ = x
        
        R1,R2,R3,Rm = ham(freq, npts, coupling_, E1_, gamma_phonon_, E2_, gamma_phonon_, E3_, gamma_phonon2_, gamma_mag, [f1_,f2_,f3_])

        R = A1_ * R1 + A2_ * R2 + A3_ * R3 + Amag_ * Rm .+ background_

        ret = sum((R - val).^2 .* weight)
#        println("ret ", ret)
        
        return ret
    end

    x0 = [coupling, E1, E2, E3, gamma_phonon, gamma_phonon2, f1,f2,f3,A1,A2,A3,Amag,background]

    if false
        res = optimize(f, x0 )
        x = Optim.minimizer(res)
    else
        x = x0
    end    

    coupling, E1, E2, E3, gamma_phonon, gamma_phonon2, f1,f2,f3,A1,A2, A3, Amag, background = x

    
    
    R1,R2,R3,Rm = ham(freq, npts, coupling, E1, gamma_phonon, E2, gamma_phonon, E3, gamma_phonon2, gamma_mag, [f1,f2,f3])
    R = A1 * R1 + A2 * R2 + A3 * R3 + Amag * Rm .+ background

    println(x)
    if doplot

        plot(freq, val, color="blue", line=:solid)

        
        display(plot!(freq, R, color="orange", line=:dash, label=""))
    end    

    return x, freq, val, R

end


function run_fit_ham(data ; npts=200, start=missing, t1=70, t2=151, doplot=true, poff = 2.0)

    if ismissing(start)
        #        x0 = [   0.09625539025761233;  82.469471880995; 107.98245601636306; 147.25455660996678;   0.6065022978970545;   130.47301375134032;   5.076080899489914;   0.13520606738509242;   4.034556458625623;   5.732477100823; 122.86862282522425;   0.8836092733259274;   0.008834640553286704]

        x0 = [0.411813;  82.4566;  108.101;  147.086;  0.610807; 0.610807; 130.79;   5.11908;  0.203162;  3.10525/10;  5.91651/10;  123.594/10;  1.17198;   -0.027727]
    else
        x0 = deepcopy(start)
    end

    gamma_mag = 0.5
    coupling, E1, E2, E3, gamma_phonon,gamma_phonon2,  f1,f2,f3,A1,A2, A3, Amag, background = x0

    N = size(data)[1]
    
    T = zeros(Float64, N, 14)

    RR = []
    
    if doplot
        plot()
    end

    offset = 0.0
    
    for n = 1:N

        println("$n xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")
        println()

        
        x, freq, val, R = fit_ham(data[n,:,:], npts;  coupling=coupling , E1=E1, E2=E2, E3=E3, gamma_phonon=gamma_phonon,gamma_phonon2=gamma_phonon2, gamma_mag = gamma_mag, f1=f1, f2 = f2, f3 = f3, A1 = A1, A2 = A2, A3 = A3, Amag = Amag, background = background, t1=t1, t2=t2, doplot=false)

        coupling, E1, E2, E3, gamma_phonon, gamma_phonon2, f1,f2,f3,A1,A2, A3, Amag, background = x

        T[n,:] = x

        push!(RR, R)
        
        if doplot
            plot!(         freq, val .+ offset, color="blue", line=:solid, label="")

            display(plot!(freq, R .+ offset , color="orange", line=:dash, label=""))

            offset += poff
        end
        
    end

    if doplot
        display(plot!())
    end

    return T, RR
    
end




function run_fit_ham_big(data ; npts=200, start=missing, t1=70, t2=151, doplot=true, poff = 2.0)

    if ismissing(start)
        x0 = [0.411813;  82.4566;  108.101;  147.086;  0.610807; 0.610807; 130.79;   5.11908;  0.203162;  3.10525/10;  5.91651/10;  123.594/10;  1.17198;   -0.027727]
    else
        x0 = deepcopy(start)
    end

    gamma_mag = 0.5
    coupling, E1, E2, E3, gamma_phonon,gamma_phonon2,  f1,f2,f3,A1,A2, A3, Amag, background = x0

    N = size(data)[1]
    
    T = zeros(Float64, N, 14)

    RR = []
    
    if doplot
        plot()
    end

    offset = 0.0
    
    for n = 1:N

        println("$n xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")
        println()

        
        x, freq, val, R = fit_ham(data[n,:,:], npts;  coupling=coupling , E1=E1, E2=E2, E3=E3, gamma_phonon=gamma_phonon,gamma_phonon2=gamma_phonon2, gamma_mag = gamma_mag, f1=f1, f2 = f2, f3 = f3, A1 = A1, A2 = A2, A3 = A3, Amag = Amag, background = background, t1=t1, t2=t2, doplot=false)

        coupling, E1, E2, E3, gamma_phonon, gamma_phonon2, f1,f2,f3,A1,A2, A3, Amag, background = x

        T[n,:] = x

        push!(RR, R)
        
        if doplot
            plot!(         freq, val .+ offset, color="blue", line=:solid, label="")

            display(plot!(freq, R .+ offset , color="orange", line=:dash, label=""))

            offset += poff
        end
        
    end

    if doplot
        display(plot!())
    end

    return T, RR
    
end



function muliti_fit_ham_big(temps, data ; npts=200, start=missing, t1=90, t2=140, doplot=true, poff = 1.0)
    #t1=70, t2=151
    
    if ismissing(start)
#        x0 = [0.411813;  82.4566;  108.101;  147.086;  0.610807; 0.610807; 130.79;   5.11908;  0.203162;  3.10525/10;  5.91651/10;  123.594/10;  1.17198;   -0.027727]

        #        x0 = [1.5, 108.101, 0.61, 130.79,5.11908,0.203162, 5.91651, 1.17198, 0.0, 80.0 , 20.0, -0.1 , -0.01 , 0.01 ]
        #        x0 = [1.0, 110, 0.61, 133 ,3.0,0.203162, 10.0, 1.17198, 0.0, 50.0 , 20.0, -0.02 , -0.01 , 0.2 ]
        #        x0 =  [0.32691685135591336, 108.19782810323653, 0.7439879555299554, 131.4081151636247, 1.0165476080215046, 0.06976731936297491, 8.589680201839695, 1.358464630502035, 0.1457806316447262, 24.916386025010222, 22.849795316808144, -0.057471067762370875, -0.0032093509243590152, 0.3043259239353536]

        x0 =  [0.90, 109.5, 0.7439879555299554, 131.4081151636247, 1.0165476080215046, 0.06976731936297491, 8.589680201839695, 1.358464630502035, 0.1457806316447262, 60, 20.0, -0.057471067762370875, -0.007, 0.35, 0.05]        
        x0 = [0.36993440791288007, 108.21274899696148, 0.6998429730008063, 131.61102812049847, 1.5322385494518524, 0.0789652195672034, 8.369405953413189, 1.346883847308227, 0.21058041952323242, 26.359924303660428, 22.401161118150615, -0.06019156812355979, -0.003301713917499625, 0.2555723977206325]

        x0 = [0.6654661973552438, 108.79721774543371, 0.5223324717900983, 132.06323925706238, 1.2298904245214135, 0.11745564217670998, 6.959911797673964, 1.7613326400519393, 0.21833814515833969, 67.20859801126461, 17.32717510700243, -0.05574040870386513, -0.004322890346137954, 0.2469631657914773, 0.01320180157759427]

#        x0 = [0.7088438188464762, 110.3448472765819, 0.4413130829327534, 132.2534249303064, 2.096783885176483, -0.0006560379556130206, 7.71675260619097, 1.817481724353533, 0.24685010285324835, 64.67746792389698, 14.414248645819448, -0.05583106319124072, -0.0075043358331730306, 0.28821145329799375, 0.007378107053313115]

#        x0 =  [0.6246505604498853, 108.90058083951655, 0.6016992015001379, 131.4166419486252, 1.3393316135328646, 0.02717953383302059, 6.355551320409767, 1.8099163361189838, 0.23148625539104592, 45.874897307803955, 16.605404692344113, -0.05360903815100381, -0.004417250899865743, 0.25992340313108064, 0.02477246311559988]

        x0 = [0.6264046139951935, 109.25998552429667, 0.5061721064506493, 131.58000359617756, 0.6009733372042408, 0.028544627167931515, 6.302766107482095, 2.283586509011992, 0.23757885738601006, 58.189189287786995, 14.792012079494532, -0.054850402144352085, -0.004823250188262926, 0.2856147955495232, 0.01086930912139069]

        #0.75
        x0 = [0.6456983703031358, 109.17020432397409, 0.6845573662069631, 131.4748812442583, -0.04350314227492411, 0.01783417892346154, 4.658684836512608, 2.3741764625628754, 0.2363534080680801, 61.75562539824748, 14.523554284962055, -0.0593046086640189, -0.004720282610935692,  -0.00003, 0.32229306336177016, 0.002491387467372571]

        x0 = [0.6557774312337368, 109.35736791068133, 0.7888368652296203, 131.00746623971943, -0.04115044943606008, 0.010050191201780053, 4.108482360370584, 2.2637555080750653, 0.23588220281221467, 59.75537648993504, 14.640172749458806, -0.05991141363080115, -0.0032261658922605346, -2.25854249331533e-5, 0.33397136114567016, 0.004876140336456884]

        x0 = [1.016277045511919, 109, 0.6958379991827707, 129.67509070124063, -0.0444641864058478, 0.02894238313933098, 4.514551844013133, 2.2494866377969878, 0.2499344638496038, 105.17211077548912, 15.650493538171116, -0.037648714584561004, 0.0049808990968535655, -0.0002083247229616496, 0.28385456193341024, -0.014915457270877874, 0.02]

        x0 = [0.6216582925495537, 108.16285407584833, 0.5896062622398496, 129.22579138549008, -0.04266621362839789, 0.005801072612091236, 5.8263446119177775, 2.0987905036845698, 0.24951884910656158, 84.5071536626195, 15.639718576114305, -0.03175101920965793, 0.004624812720576591, -0.0001920727171778543, 0.3116682898942314, -0.015163655871586346, 0.03543930692135847]

#        x0 = [0.7677683472836454, 107.9819739172338, 0.4228595745297955, 129.9270663662193, -0.04158458808015177, 0.03554476208094479, 4.493472832266291, 2.6952201536377425, 0.2417592821255878, 80.09936915858367, 13.274412971664153, -0.02675321149846141, 0.0025737425840624086, -0.00016470492466507628, 0.17770126951605913, -0.017552559980420726, 0.0405870772856814]

#        x0[1] = 0.0
        
    else
        x0 = deepcopy(start)
    end

    gamma_mag = 0.5
    coupling, E1, gamma_phonon, f1,f2,f3,A1, Amag, background, D1, Dmag, f1_lin, f1_quad, f1_cube, f2_lin, A1_lin, E1_lin  = x0


    N = size(data)[1]


    RR = []
    FREQ = []
    VAL = []

    c = 0

    temp_weights = ones(14)
    temp_weights[12:14] .= 5.0
    
    
    function sum_plot()

        plot()
 
        offset = 0.0
        
        for n = 1:N

#        println("$n xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")
#        println()

            freq = FREQ[n]
            val = VAL[n]
            R = RR[n]        
            
            if doplot
                plot!( freq, val .+ offset, color="blue", line=:solid, label="", lw=2)
                
                display(plot!(freq, R .+ offset , color="orange", line=:solid, label=""))
                
                offset += poff
            end
            
        end
        display(plot!())
        
    end

    
    function f(x)

        coupling_, E1_, gamma_phonon_, f1_,f2_,f3_,A1_, Amag_, background_, D1_, Dmag_, f1_lin_, f1_quad_, f1_cube_, f2_lin_, A1_lin_, E1_lin_ = x

#=
    println("coupling $coupling_")
    println("E1 $E1_")
    println("gamma_phonon $gamma_phonon_")
    println("f1 $f1_")
    println("f2 $f2_")
    println("f3 $f3_")
    println("A1 $A1_")
    println("Amag $Amag_")
    println("background $background_")
    println("D1 $D1_")
    println("Dmag $Dmag_")
    println("f1_lin $f1_lin_")
    println("f1_quad $f1_quad_")
    println("f1_cube $f1_cube_")
    println("f2_lin $f2_lin_")
    println("A1_lin $A1_lin_")
    println("E1_lin $E1_lin_")
=#
        
        err = 0.0

        RR = []
        FREQ = []
        VAL = []
        
        for n = 1:N

            freq = data[n,:,1]

            if !ismissing(t1)
                a1 = findall(x->x > t1 , freq )[1] - 1
                a1 = max(a1, 1)
            else
                a1 = 1
            end

            if !ismissing(t2)
                a2 = findall(x->x > t2 , freq )[1]  + 1
                a2 = min(a2, length(freq))
            else
                a2 = length(freq)
            end

            freq = data[n,a1:a2,1]
            val = data[n, a1:a2,2]

            push!(FREQ, freq)
            push!(VAL, val)
            
            weight = ones(length(freq))

            if t2 > 140.0
                wcut = findall(x->x > 140 , freq )[1] - 1
                weight[wcut:end] .= 0.01
            end

            
            temp = temps[n]
            
            A1t_ = (A1_ + A1_lin_ * temp)  * exp(- temp / D1_)
            Amagt_ = Amag_ * exp(- temp / Dmag_)

            f1t_ = f1_ + f1_lin_ * temp + f1_quad_ * temp^2 + f1_cube_ * temp^3
            f2t_ = f2_ + f2_lin_ * temp           #+ f1_lin * temp^2            


            E1t_ = E1_ + E1_lin_ * temp
            
            
            R1,R2,R3,Rm = ham(freq, npts, coupling_, E1t_, gamma_phonon_, -100.0, gamma_phonon_, +200.0 , 0.5, gamma_mag, [f1t_,f2t_,f3_])

            R = A1t_ * R1  + Amagt_ * Rm .+ background_

            push!(RR, R)
            
            err += sum((log.(R) - log.(val)).^2 .* weight) * temp_weights[n] + 0.1 * sum((R - val).^2 .* weight) * temp_weights[n]

            
        end

        println("err $err")

        c += 1
        if c > 100
            c = 0
            println("x $x")
            if doplot
                sum_plot()
            end
        end
        
        return err

    end


#    x0 = [coupling, E1, E2, E3, gamma_phonon, gamma_phonon2, f1,f2,f3,A1,A2,A3,Amag,background]
    x0 = [coupling, E1, gamma_phonon, f1,f2,f3,A1, Amag, background, D1, Dmag, f1_lin, f1_quad, f1_cube, f2_lin, A1_lin, E1_lin ]
    err = f(x0)
    if doplot
        sum_plot()
    end
        
    
    if true
        res = optimize(f, x0 )
        x = Optim.minimizer(res)
    else
        x = x0
        err = f(x0)
    end    
    
    println("final x $x")
    coupling, E1, gamma_phonon, f1,f2,f3,A1, Amag, background, D1, Dmag, f1_lin, f1_quad, f1_cube, f2_lin, A1_lin, E1_lin  = x

    
    println("coupling $coupling")
    println("E1 $E1")
    println("gamma_phonon $gamma_phonon")
    println("f1 $f1")
    println("f2 $f2")
    println("f3 $f3")
    println("A1 $A1")
    println("Amag $Amag")
    println("background $background")
    println("D1 $D1")
    println("Dmag $Dmag")
    println("f1_lin $f1_lin")
    println("f1_quad $f1_quad")
    println("f1_cube $f1_cube")
    println("f2_lin $f2_lin")
    println("A1_lin $A1_lin")
    println("E1_lin $E1_lin")
    
    
    if doplot
        sum_plot()
    end

    return x
    
end


function muliti_fit_ham_big_more(temps, data ; npts=200, start=missing, t1=70, t2=151 , doplot=true, poff = 1.0)
    #t1=70, t2=151
    
    if ismissing(start)
#        x0 = [0.411813;  82.4566;  108.101;  147.086;  0.610807; 0.610807; 130.79;   5.11908;  0.203162;  3.10525/10;  5.91651/10;  123.594/10;  1.17198;   -0.027727]

        #        x0 = [1.5, 108.101, 0.61, 130.79,5.11908,0.203162, 5.91651, 1.17198, 0.0, 80.0 , 20.0, -0.1 , -0.01 , 0.01 ]
        #        x0 = [1.0, 110, 0.61, 133 ,3.0,0.203162, 10.0, 1.17198, 0.0, 50.0 , 20.0, -0.02 , -0.01 , 0.2 ]
        #        x0 =  [0.32691685135591336, 108.19782810323653, 0.7439879555299554, 131.4081151636247, 1.0165476080215046, 0.06976731936297491, 8.589680201839695, 1.358464630502035, 0.1457806316447262, 24.916386025010222, 22.849795316808144, -0.057471067762370875, -0.0032093509243590152, 0.3043259239353536]

        x0 =  [0.90, 109.5, 0.7439879555299554, 131.4081151636247, 1.0165476080215046, 0.06976731936297491, 8.589680201839695, 1.358464630502035, 0.1457806316447262, 60, 20.0, -0.057471067762370875, -0.007, 0.35, 0.05]        
        x0 = [0.36993440791288007, 108.21274899696148, 0.6998429730008063, 131.61102812049847, 1.5322385494518524, 0.0789652195672034, 8.369405953413189, 1.346883847308227, 0.21058041952323242, 26.359924303660428, 22.401161118150615, -0.06019156812355979, -0.003301713917499625, 0.2555723977206325]

        x0 = [0.6654661973552438, 108.79721774543371, 0.5223324717900983, 132.06323925706238, 1.2298904245214135, 0.11745564217670998, 6.959911797673964, 1.7613326400519393, 0.21833814515833969, 67.20859801126461, 17.32717510700243, -0.05574040870386513, -0.004322890346137954, 0.2469631657914773, 0.01320180157759427]

#        x0 = [0.7088438188464762, 110.3448472765819, 0.4413130829327534, 132.2534249303064, 2.096783885176483, -0.0006560379556130206, 7.71675260619097, 1.817481724353533, 0.24685010285324835, 64.67746792389698, 14.414248645819448, -0.05583106319124072, -0.0075043358331730306, 0.28821145329799375, 0.007378107053313115]

#        x0 =  [0.6246505604498853, 108.90058083951655, 0.6016992015001379, 131.4166419486252, 1.3393316135328646, 0.02717953383302059, 6.355551320409767, 1.8099163361189838, 0.23148625539104592, 45.874897307803955, 16.605404692344113, -0.05360903815100381, -0.004417250899865743, 0.25992340313108064, 0.02477246311559988]

        x0 = [0.6264046139951935, 109.25998552429667, 0.5061721064506493, 131.58000359617756, 0.6009733372042408, 0.028544627167931515, 6.302766107482095, 2.283586509011992, 0.23757885738601006, 58.189189287786995, 14.792012079494532, -0.054850402144352085, -0.004823250188262926, 0.2856147955495232, 0.01086930912139069]

        #0.75
        x0 = [0.6456983703031358, 109.17020432397409, 0.6845573662069631, 131.4748812442583, -0.04350314227492411, 0.01783417892346154, 4.658684836512608, 2.3741764625628754, 0.2363534080680801, 61.75562539824748, 14.523554284962055, -0.0593046086640189, -0.004720282610935692,  -0.00003, 0.32229306336177016, 0.002491387467372571]

        x0 = [0.6557774312337368, 109.35736791068133, 0.7888368652296203, 131.00746623971943, -0.04115044943606008, 0.010050191201780053, 4.108482360370584, 2.2637555080750653, 0.23588220281221467, 59.75537648993504, 14.640172749458806, -0.05991141363080115, -0.0032261658922605346, -2.25854249331533e-5, 0.33397136114567016, 0.004876140336456884]

        x0 = [1.016277045511919, 109, 0.6958379991827707, 129.67509070124063, -0.0444641864058478, 0.02894238313933098, 4.514551844013133, 2.2494866377969878, 0.2499344638496038, 105.17211077548912, 15.650493538171116, -0.037648714584561004, 0.0049808990968535655, -0.0002083247229616496, 0.28385456193341024, -0.014915457270877874, 0.02]

        #VV good
        x0 = [0.6216582925495537, 108.16285407584833, 0.5896062622398496, 129.22579138549008, -0.04266621362839789, 0.005801072612091236, 5.8263446119177775, 2.0987905036845698, 0.24951884910656158, 84.5071536626195, 15.639718576114305, -0.03175101920965793, 0.004624812720576591, -0.0001920727171778543, 0.3116682898942314, -0.015163655871586346, 0.03543930692135847]

#        x0[1] = 0.0
        
        #VH good
#        x0 = [0.7677683472836454, 107.9819739172338, 0.4228595745297955, 129.9270663662193, -0.04158458808015177, 0.03554476208094479, 4.493472832266291, 2.6952201536377425, 0.2417592821255878, 80.09936915858367, 13.274412971664153, -0.02675321149846141, 0.0025737425840624086, -0.00016470492466507628, 0.17770126951605913, -0.017552559980420726, 0.0405870772856814]

#        x0[1] = 0.0

#        x2 = [82, 2 , 40, 0.01, 146.5, 50, 40, 0.01] 
#        x2 = [82.6357850826726, 2.722351711340347, 56.21106714674642, -0.016055817832190934, 146.8002684966559, 193.1921338129491, 29.4694679440709, -0.019504891214054146, 0.6]

        #VV good
        x2 = [82.5960653258821, 3.6264588233480333, 58.06820522697674, -0.02712998956684886, 146.76047033603248, 193.05459925680287, 34.0187847575256, -0.01061495594554369, 0.6111736131584395, -0.007963257884389823]
#        x2 = [82.63627749546839, 2.7506724738482538, 56.0787898143584, -0.016321766735781874, 146.73409643804197, 215.7322107918633, 32.82865497171554, -0.017594549660580713, 0.4439883702390163, 0.0]

        #        x2 = [82.63627749546839, 2.7506724738482538, 56.0787898143584, -0.016321766735781874,  147.11602605586614, 22.095587247807313, 31.447256515347114, 0.0009163238952651408, 0.4279867822606822]
        #x2 = [82.63882014600242, 3.6097901240438524, 53.886582704530035, -0.04646390797790811,  147.11602605586614, 22.095587247807313, 31.447256515347114, 0.0009163238952651408, 0.4279867822606822]
        #        x2 = [82.64021658937791, 3.3482780934752814, 40.0, 0.0, 147.09407344179206, 27.951045530822537, 25.99918834391115, 0.0009316161180310912, 0.6057374026562841]

#        x2 = [82.64132104245196, 4.584127998853803, 20.48656989223468, 0.001, 147.09056989938293, 27.93109258661703, 25.765110554243282, 0.010082858696903164, 0.6074814979087447, 0.01]

        #VH good
#        x2 = [82.66107308408675, 4.169018590107433, 28.600194476889136, 0.001, 147.0788883829565, 26.8494209952009, 26.052318744478097, 0.044334027384477045, 0.6410064009526978, -0.10996015580351662]
    else
        x0 = deepcopy(start)
    end

    gamma_mag = 0.5
    coupling, E1, gamma_phonon, f1,f2,f3,A1, Amag, background, D1, Dmag, f1_lin, f1_quad, f1_cube, f2_lin, A1_lin, E1_lin  = x0

    E2, A2, D2, A2_lin, E3, A3, D3, A3_lin, gamma_phonon2, back2 = x2

    N = size(data)[1]


    RR = []
    FREQ = []
    VAL = []

    RR1 = []
    RR2 = []
    RR3 = []
    RRm = []

    
    c = 0

    temp_weights = ones(14)
    temp_weights[12:14] .= 5.0
    
    
    function sum_plot()

        plot()
 
        offset = 0.0
        
        for n = 1:N

#        println("$n xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")
#        println()

            freq = FREQ[n]
            val = VAL[n]
            R = RR[n]        

            R1 = RR1[n]
            R2 = RR2[n]
            R3 = RR3[n]
            Rm = RRm[n]
            
            
            if doplot
#                plot!( freq, val .+ offset, color="blue", line=:solid, label="", lw=2)
#                display(plot!(freq, R .+ offset , color="orange", line=:solid, label=""))

                if n == 1

#                    plot!( freq, val .+ offset, color="blue", line=:solid, label="Experiment", lw=2)
#                    display(plot!(freq, R .+ offset , color="orange", line=:solid, label="Model"))

                    display(plot!(freq, R3 .+ offset , color="yellow", line=:solid, label="phonon3"))
                    display(plot!(freq, R1 .+ offset , color="magenta", line=:solid, label="phonon2"))
                    display(plot!(freq, Rm .+ offset , color="cyan", line=:solid, label="2mag"))
                    display(plot!(freq, R2 .+ offset , color="green", line=:solid, label="phonon1"))
                else
                    display(plot!(freq, R3 .+ offset , color="yellow", line=:solid, label=""))
                    display(plot!(freq, R1 .+ offset , color="magenta", line=:solid, label=""))
                    display(plot!(freq, Rm .+ offset , color="cyan", line=:solid, label=""))
                    display(plot!(freq, R2 .+ offset , color="green", line=:solid, label=""))

#                    plot!( freq, val .+ offset, color="blue", line=:solid, label="", lw=2)
#                    display(plot!(freq, R .+ offset , color="orange", line=:solid, label=""))

                end                
                
                offset += poff
            end
            
        end
        ylims!(0, 25)
        xlabel!("Freq (cm-1)")
        ylabel!("Raman Intensity")
        
        display(plot!(framestyle=:box))
        
    end

    
    function f(x)

        coupling_, E1_, gamma_phonon_, f1_,f2_,f3_,A1_, Amag_, background_, D1_, Dmag_, f1_lin_, f1_quad_, f1_cube_, f2_lin_, A1_lin_, E1_lin_ = x0


        E2_, A2_, D2_, A2_lin_, E3_, A3_, D3_, A3_lin_, gamma_phonon2_, back2_ = x

#        println([E2_, A2_, D2_, A2_lin_, E3_, A3_, D3_, A3_lin_, gamma_phonon2_])
        
#=
    println("coupling $coupling_")
    println("E1 $E1_")
    println("gamma_phonon $gamma_phonon_")
    println("f1 $f1_")
    println("f2 $f2_")
    println("f3 $f3_")
    println("A1 $A1_")
    println("Amag $Amag_")
    println("background $background_")
    println("D1 $D1_")
    println("Dmag $Dmag_")
    println("f1_lin $f1_lin_")
    println("f1_quad $f1_quad_")
    println("f1_cube $f1_cube_")
    println("f2_lin $f2_lin_")
    println("A1_lin $A1_lin_")
    println("E1_lin $E1_lin_")
=#
        
        err = 0.0

        RR = []
        FREQ = []
        VAL = []

        RR1 = []
        RR2 = []
        RR3 = []
        RRm = []
        
        
        for n = 1:N

            freq = data[n,:,1]

            if !ismissing(t1)
                a1 = findall(x->x > t1 , freq )[1] - 1
                a1 = max(a1, 1)
            else
                a1 = 1
            end

            if !ismissing(t2)
                a2 = findall(x->x > t2 , freq )[1]  + 1
                a2 = min(a2, length(freq))
            else
                a2 = length(freq)
            end

            freq = data[n,a1:a2,1]
            val = data[n, a1:a2,2]

            push!(FREQ, freq)
            push!(VAL, val)
            
            weight = ones(length(freq))

            if t2 > 140.0
                wcut = findall(x->x > 140 , freq )[1] - 1
                weight[wcut:end] .= 0.01
            end

            
            temp = temps[n]
            
            A1t_ = (A1_ + A1_lin_ * temp)  * exp(- temp / D1_)
            Amagt_ = Amag_ * exp(- temp / Dmag_)

            f1t_ = f1_ + f1_lin_ * temp + f1_quad_ * temp^2 + f1_cube_ * temp^3
            f2t_ = f2_ + f2_lin_ * temp           #+ f1_lin * temp^2            


            E1t_ = E1_ + E1_lin_ * temp

            #            A2t_ = (A2_ + A2_lin_ * temp)  * exp(- temp / D2_)
            A2t_ = (A2_ +  A2_lin_ * temp)  * exp(- temp / D2_)
            A3t_ = (A3_ + A3_lin_ * temp)  * exp(- temp / D3_)

            
            
            R1,R2,R3,Rm = ham(freq, npts, coupling_, E1t_, gamma_phonon_, E2_, gamma_phonon_, E3_ , gamma_phonon2_ , gamma_mag, [f1t_,f2t_,f3_])

            R = A1t_ * R1  + Amagt_ * Rm .+ background_ .+ R2 * A2t_ .+ R3 * A3t_ .+ back2_

            Ra1 =A1t_ * R1
            Ra2 =A2t_ * R2
            Ra3 =A3t_ * R3
            Ram = Amagt_ * Rm

#            Ra1 =  R1
#            Ra2 =  R2
#            Ra3 = R3
#            Ram = Rm / 5
            

            
            push!(RR, deepcopy(R))

            push!(RR1, Ra1)
            push!(RR2, Ra2)
            push!(RR3, Ra3)
            push!(RRm, Ram)
            

            a10 = findall(X->X > 100.0 , freq )[1] - 1
            a10 = max(a10, 1)

            a20 = findall(X->X > 140.0 , freq )[1] - 1
            a20 = min(a20, length(val))


            try
#                err += sum( (log.( R[1:a10]) - log.(val[1:a10])).^2 .* weight[1:a10] )* temp_weights[n]
#                err += sum( (log.(R[a20:end]) - log.(val[a20:end])).^2 .* weight[a20:end] ) * temp_weights[n]

                err += sum( (( R[1:a10]) - (val[1:a10])).^2 .* weight[1:a10] )* temp_weights[n]
                err += sum( ((R[a20:end]) - (val[a20:end])).^2 .* weight[a20:end] ) * temp_weights[n]
            catch
                err = 1000.0
            end
                
#            println(a10, " a10 ", freq[a10])
#            println(a20, " a20 ", freq[a20])
            
            #            err += sum((log.(R) - log.(val)).^2 .* weight) * temp_weights[n] + 0.1 * sum((R - val).^2 .* weight) * temp_weights[n]

#            plot(freq[1:a10], R[1:a10] - val[1:a10])
#            plot!(freq[1:a10], R[1:a10])
#            display(plot!(freq[1:a10],  val[1:a10]))

#            plot!(freq[a20:end], R[a20:end] - val[a20:end])
#            plot!(freq[a20:end], R[a20:end])
#            display(plot!(freq[a20:end],  val[a20:end]))
            
#            sleep(5.0)
            
#            err += sum(( log.( max.(R[1:a10], 1e-5))  - log.(val[1:a10])).^2 .* weight[1:a10]) * temp_weights[n] + 0.1 * sum((R[1:a10] - val[1:a10]).^2 .* weight[1:a10]) * temp_weights[n]
#            err += sum((log.(R[a20:end]) - log.(val[a20:end])).^2 .* weight[a20:end]) * temp_weights[n] + 0.1 * sum((R[a20:end] - val[a20:end]).^2 .* weight[a20:end]) * temp_weights[n]
#            err = 0.0
            
        end

        println("err $err")

        c += 1
        if c > 100
            c = 0
            println("x $x")
            if doplot
                sum_plot()
            end
        end
        
        return err

    end


#    x0 = [coupling, E1, E2, E3, gamma_phonon, gamma_phonon2, f1,f2,f3,A1,A2,A3,Amag,background]
    x0 = [coupling, E1, gamma_phonon, f1,f2,f3,A1, Amag, background, D1, Dmag, f1_lin, f1_quad, f1_cube, f2_lin, A1_lin, E1_lin ]

    err = f(x2)
    if doplot
        sum_plot()
    end
        
    
    if false
        res = optimize(f, x2 )
        x = Optim.minimizer(res)
    else
        x = x2
        err = f(x2)
    end    
    
    println("final x $x")
    coupling, E1, gamma_phonon, f1,f2,f3,A1, Amag, background, D1, Dmag, f1_lin, f1_quad, f1_cube, f2_lin, A1_lin, E1_lin  = x0
    E2, A2, D2, A2_lin, E3, A3, D3, A3_lin, gamma_phonon2, back2 = x

    
    println("coupling $coupling")
    println("E1 $E1")
    println("gamma_phonon $gamma_phonon")
    println("f1 $f1")
    println("f2 $f2")
    println("f3 $f3")
    println("A1 $A1")
    println("Amag $Amag")
    println("background $background")
    println("D1 $D1")
    println("Dmag $Dmag")
    println("f1_lin $f1_lin")
    println("f1_quad $f1_quad")
    println("f1_cube $f1_cube")
    println("f2_lin $f2_lin")
    println("A1_lin $A1_lin")
    println("E1_lin $E1_lin")
    
    
    if doplot
        sum_plot()
    end

    return x
    
end


