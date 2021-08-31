"""
    function main_fit( data ; mode=:plot_only, npts=400, start=missing, t1=70, t2=151 , doplot=true, poff = 1.0)

This is the main fitting code for the global fit. It should be called by
the helper function below main_fit_VV and main_fit_VH for VV and VH data.

"""
function main_fit( data ; mode=:plot_only, npts=400, start=missing, t1=70, t2=151 , doplot=true, poff = 1.0)

    println("TEMPs ", temps)
    
    if ismissing(start)
        #VH global 2
        x0 = [0.6835049066786283 * sqrt(200), 108.358106047425, 0.4543220022001645, 128.3650860901054, -0.04105257349414344, 0.02797865938260505, 6.488216937453528, 3.005566008569702, 0.1206692152863844, 95.40560187276805, 14.520817580740925, -0.027582651488746123, 0.005891826941019795, -0.00019900159172185888, 0.29434884664864885, -0.009780255038509683, 0.03281605040568809]
        x2 = [82.80035101271274, 4.065091207163782, 26.833584014264982, 0.015159946747169563, 146.96607541862306, 31.471772591116984, 24.935170319025428, 0.05872107317439383, 0.41007703434565074, -0.11729931180742176]
        
    else
        x0 = deepcopy(start[1])
        x2 = deepcopy(start[2])
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
    temp_weights[12:14] .= 2.0
    

    #helper function for plotting
    function sum_plot()

        L = plot(layout=(2,1), framestyle=:box, size=(500,1000))

        ylims!(0, 25)
        xlabel!("Freq (cm-1)")
        ylabel!("Raman Intensity")
        
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

                    
                    
                    plot!(L[1], freq, val .+ offset, color="blue", line=:solid, label="Experiment", lw=2)
                    display(plot!(L[1], freq, R .+ offset , color="orange", line=:solid, label="Model"))

                    display(plot!(L[2], freq, R2 .+ offset , color="green", line=:solid, label="phonon1"))
                    display(plot!(L[2], freq, R1 .+ offset , color="magenta", line=:solid, label="phonon2"))
                    display(plot!(L[2], freq, R3 .+ offset , color="yellow", line=:solid, label="phonon3"))
                    display(plot!(L[2], freq, Rm .+ offset , color="cyan", line=:solid, label="2mag"))
                else
                    display(plot!(L[2], freq, R3 .+ offset , color="yellow", line=:solid, label=""))
                    display(plot!(L[2], freq, R1 .+ offset , color="magenta", line=:solid, label=""))
                    display(plot!(L[2], freq, Rm .+ offset , color="cyan", line=:solid, label=""))
                    display(plot!(L[2], freq, R2 .+ offset , color="green", line=:solid, label=""))

                    plot!(L[1],  freq, val .+ offset, color="blue", line=:solid, label="", lw=2)
                    display(plot!(L[1], freq, R .+ offset , color="orange", line=:solid, label=""))

                end                
                
                offset += poff
            end
            
        end
        
        display(plot!())

        #uncomment to write data
#        write("model_VH_2/model_all", FREQ, RR1 + RR2 + RR3 + RRm)
#        write("model_VH_2/model_phonon1", FREQ, RR2)
#        write("model_VH_2/model_phonon2", FREQ, RR1)
#        write("model_VH_2/model_phonon3", FREQ, RR3)
#        write("model_VH_2/model_2mag", FREQ, RRm)

#        write("model_VV_2/model_all", FREQ, RR1 + RR2 + RR3 + RRm)
#        write("model_VV_2/model_phonon1", FREQ, RR2)
#        write("model_VV_2/model_phonon2", FREQ, RR1)
#        write("model_VV_2/model_phonon3", FREQ, RR3)
#        write("model_VV_2/model_2mag", FREQ, RRm)
        
        
    end


    #this is the actual function you minimize
    function f(x)

        coupling_, E1_, gamma_phonon_, f1_,f2_,f3_,A1_, Amag_, background_, D1_, Dmag_, f1_lin_, f1_quad_, f1_cube_, f2_lin_, A1_lin_, E1_lin_,E2_, A2_, D2_, A2_lin_, E3_, A3_, D3_, A3_lin_, gamma_phonon2_, back2_ = x


        
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

            #downweight higher frequencies.
            if t2 > 140.0
                wcut = findall(x->x > 140 , freq )[1] - 1
                weight[wcut:end] .= 0.1
            end

            
            temp = temps[n]


            #amplitudes
            A1t_ = (A1_ + A1_lin_ * temp)  * exp(- temp / D1_)
            Amagt_ = Amag_ * exp(- temp / Dmag_)

            #resonances
            f1t_ = f1_ + f1_lin_ * temp + f1_quad_ * temp^2 + f1_cube_ * temp^3
            f2t_ = f2_ + f2_lin_ * temp           #+ f1_lin * temp^2            


            E1t_ = E1_ +  E1_lin_ * temp
#            E1t_ = E1_ + 0.0 * E1_lin_ * temp

            #            A2t_ = (A2_ + A2_lin_ * temp)  * exp(- temp / D2_)
            A2t_ = (A2_ + 0.0 *  A2_lin_ * temp)  * exp(- temp / D2_)
            A3t_ = (A3_ + A3_lin_ * temp)  * exp(- temp / D3_)

            
            #setup and solve the hamiltonian


            R1,R2,R3,Rm = ham(freq, npts, coupling_, E1t_, gamma_phonon_, E2_, gamma_phonon_, E3_, gamma_phonon2_ , gamma_mag, [f1t_,f2t_,f3_])

            R = A1t_ * R1  + Amagt_ * Rm .+ background_ .+ R2 * A2t_ .+ R3 * A3t_ 

            Ra1 =A1t_ * R1
            Ra2 =A2t_ * R2
            Ra3 =A3t_ * R3
            Ram = Amagt_ * Rm

#            Ra1 =  R1
#            Ra2 =  R2
#            Ra3 = R3
 #           Ram = Rm / 5
            

            
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
                #this is the actual quantity to minimize.
                
                err += sum((log.(R) - log.(val)).^2 .* weight) * temp_weights[n] + 0.05 * sum((R - val).^2 .* weight) * temp_weights[n]
                
            catch
                err = 10000.0
            end
            
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


    err = f([x0;x2])
    if doplot
        sum_plot()
    end
        
    
    if mode != :plot_only
        #        res = optimize(f, [x0;x2] , ConjugateGradient())
        #        res = optimize(f, [x0;x2] , ParticleSwarm())
        res = optimize(f, [x0;x2] )
        x = Optim.minimizer(res)
    else
        x = [x0;x2]
        err = f([x0;x2])
    end    
    
    println("final x $x")
    println("----")
    println()
    coupling, E1, gamma_phonon, f1,f2,f3,A1, Amag, background, D1, Dmag, f1_lin, f1_quad, f1_cube, f2_lin, A1_lin, E1_lin,E2, A2, D2, A2_lin, E3, A3, D3, A3_lin, gamma_phonon2, back2  = x

    #    E2, A2, D2, A2_lin, E3, A3, D3, A3_lin, gamma_phonon2, back2 = x


    
    println("coupling (lambda) $coupling")
    println("E1 (freq of phonon2) $E1")
    println("E1_lin (freq of phonon2 linear in temperature) $E1_lin")
    println()
    println("E2 (freq of phonon1) $E2")
    println("E3 (freq of phonon3) $E3")
    println()
    println("gamma_phonon (spread phonon1 and phonon2) $gamma_phonon")
    println("gamma_phonon2 (spread phonon3 ) $gamma_phonon")
    println()
    println("f1 2magnon resonance zero temperature $f1")
    println("f1_lin 2mag resonance linear $f1_lin")
    println("f1_quad 2mag resonance quadratic $f1_quad")
    println("f1_cube 2mag resonance cubic $f1_cube")
    println()
    println("f2 2magnon gamma zero temp $f2")
    println("f2_lin 2magnon gamma linear in temp $f2_lin")
    println("f3 magnon q $f3")
    println()
    println("A1 phonon2 amplitude $A1")
    println("A1_lin phonon2 amplitude linear temperature $A1_lin")
    println("D1 phonon2 decay constant $D1")
    println()

    println("A2 phonon1 amplitude $A2")
    println("D2 phonon1 decay constant $D2")
    println()

    println("A3 phonon3 amplitude $A3")
    println("A3_lin phonon3 amplitude linear temperature $A3_lin")
    println("D3 phonon3 decay constant $D3")
    println()
    
    println("Amag 2mag amplitude $Amag")
    println("Dmag 2mag decay constant $Dmag")
    println()
    println("background $background")
    println()

    
    
#    if doplot
#        sum_plot()
#    end

    return x
    
end


"""
    function main_fit_VH(;mode=:plot_only, npts=400, t1=70, t2=151 , doplot=true, poff = 1.0)

This is the main fitting code for VH data. If run with default values it will plot the global fit.
Otherwise, you can start from that fit and try to do optimize yourself.

- If `mode` is not equal to `:plot_only`, then the code will actually do the fitting.
- `npts=400` controls the number of points used to approximate the continuum.
- `t1=70` and `t2=151` are the lower and upper bounds in frequency, in cm-1
- `doplot=false` to avoid plotting
- `poff = 1.0` The offset between plots.


"""
function main_fit_VH(;mode=:plot_only, npts=200, t1=70, t2=151 , doplot=true, poff = 1.0)

    #starting data for VH

    x0 = [0.6835049066786283 * sqrt(200), 108.358106047425, 0.4543220022001645, 128.3650860901054, -0.04105257349414344, 0.02797865938260505, 6.488216937453528, 3.005566008569702, 0.1206692152863844, 95.40560187276805, 14.520817580740925, -0.027582651488746123, 0.005891826941019795, -0.00019900159172185888, 0.29434884664864885, -0.009780255038509683, 0.03281605040568809]
    x2 = [82.80035101271274, 4.065091207163782, 26.833584014264982, 0.015159946747169563, 146.96607541862306, 31.471772591116984, 24.935170319025428, 0.05872107317439383, 0.41007703434565074, -0.11729931180742176]

    x0 = [42.14301427991219, 108.96458245968437,   0.3671497196891554, 131.86018169096064,  -0.04237922781509023,   6.370350562998323,   8.066942154147569,   7.4310704587847995,   0.1853291908297635,  70.32697165332858,  47.74174812637775,  -0.021699828192152402,  -0.005632677375607418,   1.0483113989854774e-5,   0.17014831448192136,  -0.008659521568708695,   0.015334192361995843]
    x2 = [82.92777517212615,   3.235297036432314,  48.2473670137221,  -0.02819672188709332, 146.3509429922471, 219.23958292338347/ 10.0,  33.910871992806264,  -0.012360658504345236,   0.48356872680308777,   0.014446038148484778]


    x0 = [39.944688885126894, 108.96458245968437,   0.2979372943816094, 131.20872720753079,  -0.04412596473809375,  23.5541549058618,  10.850581418304865,  13.123961854865605,   0.17973813235992037,  53.01493514857207,  21.918482755921552,  -0.026140080811667223,  -0.004882761369163288,  -4.233132922994685e-6,   0.11225451245561427,  -0.015237482887206534,   0.020594513658826807]
    x2 = [82.92777517212615,   3.1340936113663957,  28.60796533780442,  -0.031488828522276864, 146.3509429922471,  25.4797608432849,  26.114016384967915,  -0.021084891415018635,   0.3552090536534399,   0.0018599038339241392]


    x0 = [  31.590260395890922,  108.32899192012053,    0.43230892593438747,  130.68450051192426,   -0.04545134821652863,   35.5644777334443,    7.253416542272829,   12.31108252740858,    0.08166540365643941,   55.24142556727246,   36.60222026697,   -0.03229474292474604,   -0.006221216638113039,    3.667839715727892e-5,    0.1880928386526172,   -0.01856441126577166,    0.03604905626556161]
    x2 = [82.77587614171492,    3.3283891673673267,   37.05016141027897,   -0.03774178272505245,  146.88999963854502,   27.107834582349586,   37.739432597987616,   -0.01778241680426505,    0.3000185023896063,   -0.013769888918831522]

    x0[1] = 42.0
    x0[14] = -6e-5 

    x0=[36.1303567246408, 108.72913153006768,   0.4415732761921268, 130.34937026581417,  -0.04457839993256533,  42.72236771512939,   7.013903245162495,  12.3228341381047,   0.08491707281584261,  60.911598142876834,  35.07721601149555,  -0.033109767338510426,  -0.004625809345159903,   5.425431999454724e-6,   0.18484697779281822,  -0.016800526790898596,   0.0275283522272423]
    x2 = [82.86403764704566,   3.0989285279819083,  37.69724851853046,  -0.037220650615825725, 146.70654701375233,  30.865248044467876,  36.28155557186203,  -0.016990159143966947,   0.29138977249961384,  -0.008315888884396126]

    x0[14] = -1e-5 
    
    
    x = main_fit(VH,  start = (x0, x2), mode=mode, npts=npts, t1=70, t2=151 , doplot=doplot, poff = 1.0)

    x0 = [34.07623603880653, 108.56185726011745,   0.4347700072986739, 130.09877716612502,  -0.0444677435538983,  58.09582747257295,   6.6718032436071955,  12.679099139302899,   0.0872064391505721,  66.95020673046649,  33.71605266088026,  -0.03154409674031596,  -0.004107110591766565,  -3.3958354658766757e-6,   0.18527437391825585,  -0.012682686748259531,   0.03009348989262881]
    x2 = [82.8251279510044,   3.1539452309207765,  37.07264984237419,  -0.037190798905708124, 146.82578800737767,  27.689866098499152,  44.27211310623068,  -0.017900057097098677,   0.26536618966240666,  -0.00950598863837015]


    
    return x
    
end



"""
    function main_fit_VV(;mode=:plot_only, npts=400, t1=70, t2=151 , doplot=true, poff = 1.0)

Same, but for VV data

"""
function main_fit_VV(;mode=:plot_only, npts=200, t1=70, t2=151 , doplot=true, poff = 1.0)


    #starting data for VV
#    x0 = [0.6742588756733324 * sqrt(200) , 107.9573152290177, 0.539405630710322, 131.7904841400483, -0.04116806432910838, 0.13216936944670651, 5.3202665619743295, 100*2.863376011682241, 0.21242265395655052, 93.10845023705731, 17.26445320337706, -0.02813861848258876, -0.003541611719747818, -3.3065858988794736e-5, 0.22643428157217202, -0.017448307360460862, 0.02838862190628346]

    x0 = [5* 0.6742588756733324 * sqrt(200) , 107.9573152290177, 0.539405630710322, 131.7904841400483, -0.04116806432910838, 10.0, 5.3202665619743295, 3*2.863376011682241, 0.21242265395655052, 93.10845023705731, 2.0*17.26445320337706, -0.02813861848258876, -0.003541611719747818, -3.3065858988794736e-5, 0.22643428157217202, -0.017448307360460862, 0.02838862190628346]

#    x0[1] = 0.0
    
    x0[14] = -6e-5 

    x2 = [82.63885069507708, 2.6884950644927446, 46.08455682641799, -0.026069940793432625, 146.83426194256313, 177.7018681987727, 35.86738948611662, -0.008194054132508007, 0.61068408400651, 0.01040673628717732]

    x0 = [  45.92307209251176,  107.9573152290177,   0.38882742860279007, 131.5423696487797,  -0.0413588132095163,   3.86998116865837,   7.109832364779594,   5.68318649979629,   0.19387872289653235,  92.22074914414094,  45.00807078232238,  -0.010032181271378025,  -0.004240510990579296,  -1.3788536278742538e-5,   0.14236958036963127,  -0.009656814275695962,   0.027154196271006206]
    x2 = [82.63885069507708,   3.3728569479938417,  43.04109774704145,  -0.027160335617558068, 146.83426194256313, 158.13580351443466,  35.87334479498228,  -0.020359313215578618,   0.6854830478876333,   0.018621545690937213]

    x0 = [42.14301427991219, 108.96458245968437,   0.3671497196891554, 131.86018169096064,  -0.04237922781509023,   6.370350562998323,   8.066942154147569,   7.4310704587847995,   0.1853291908297635,  70.32697165332858,  47.74174812637775,  -0.021699828192152402,  -0.005632677375607418,   1.0483113989854774e-5,   0.17014831448192136,  -0.008659521568708695,   0.015334192361995843]
    x2 = [82.92777517212615,   3.235297036432314,  48.2473670137221,  -0.02819672188709332, 146.3509429922471, 219.23958292338347,  33.910871992806264,  -0.012360658504345236,   0.48356872680308777,   0.014446038148484778]

    x0[14] = -1e-5 
    

    x0 = [40.56858406766799,  108.5849919608786,   0.36979237859061115, 130.84981701659277,  -0.0426226760320891,   7.707842815646611,   8.841667006800092,   8.332181658088437,   0.1959801354917807,  74.95412012862728,  40.66839418410081,  -0.019509037776960586,  -0.0032427479956535273,  -2.727191489644125e-5,   0.16737019616221638,  -0.013274567933279172,   0.022046342100459776]
    x2 = [82.8818958602397,   3.710974748308286,  52.28361185263288,  -0.028226951127929736, 146.40269400693694, 245.7452949835227,  33.96893225416177,  -0.016520206984805462,   0.42816024660732616,   0.01624922314433173]
    
    x = main_fit(VV,  start = (x0, x2), mode=mode, npts=npts, t1=70, t2=151 , doplot=doplot, poff = 1.0)

    return x
    
end

