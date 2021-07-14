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
#        write("model_VH/model_all", FREQ, RR1 + RR2 + RR3 + RRm)
#        write("model_VH/model_phonon1", FREQ, RR2)
#        write("model_VH/model_phonon2", FREQ, RR1)
#        write("model_VH/model_phonon3", FREQ, RR3)
#        write("model_VH/model_2mag", FREQ, RRm)

        
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
            R1,R2,R3,Rm = ham(freq, npts, coupling_, E1t_, gamma_phonon_, E2_, gamma_phonon_, E3_ , gamma_phonon2_ , gamma_mag, [f1t_,f2t_,f3_])

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
function main_fit_VH(;mode=:plot_only, npts=400, t1=70, t2=151 , doplot=true, poff = 1.0)

    #starting data for VH

    x0 = [0.6835049066786283 * sqrt(200), 108.358106047425, 0.4543220022001645, 128.3650860901054, -0.04105257349414344, 0.02797865938260505, 6.488216937453528, 3.005566008569702, 0.1206692152863844, 95.40560187276805, 14.520817580740925, -0.027582651488746123, 0.005891826941019795, -0.00019900159172185888, 0.29434884664864885, -0.009780255038509683, 0.03281605040568809]
    x2 = [82.80035101271274, 4.065091207163782, 26.833584014264982, 0.015159946747169563, 146.96607541862306, 31.471772591116984, 24.935170319025428, 0.05872107317439383, 0.41007703434565074, -0.11729931180742176]

    x = main_fit(VH,  start = (x0, x2), mode=:plot_only, npts=npts, t1=70, t2=151 , doplot=true, poff = 1.0)

    return x
    
end



"""
    function main_fit_VV(;mode=:plot_only, npts=400, t1=70, t2=151 , doplot=true, poff = 1.0)

Same, but for VV data

"""
function main_fit_VV(;mode=:plot_only, npts=400, t1=70, t2=151 , doplot=true, poff = 1.0)


    #starting data for VV
    x0 = [0.6742588756733324 * sqrt(200) , 107.9573152290177, 0.539405630710322, 131.7904841400483, -0.04116806432910838, 0.13216936944670651, 5.3202665619743295, 2.863376011682241, 0.21242265395655052, 93.10845023705731, 17.26445320337706, -0.02813861848258876, -0.003541611719747818, -3.3065858988794736e-5, 0.22643428157217202, -0.017448307360460862, 0.02838862190628346]
    x0[14] = -6e-5 

    x2 = [82.63885069507708, 2.6884950644927446, 46.08455682641799, -0.026069940793432625, 146.83426194256313, 177.7018681987727, 35.86738948611662, -0.008194054132508007, 0.61068408400651, 0.01040673628717732]

    
    x = main_fit(VV,  start = (x0, x2), mode=:plot_only, npts=npts, t1=70, t2=151 , doplot=true, poff = 1.0)

    return x
    
end

