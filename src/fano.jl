
function fano(x, E_res, G, q)

    return 1.0 .- (q * G / 2 .+ (x .- E_res)	 ).^2 ./ ( (G/2)^2 .+ (x .- E_res).^2)

end

x = 0:.1:160

f1a = fano(x,  130, 15 * 1, 0.1) ; 
f30a = fano(x, 125, 15 * 2, 0.1) ; 
f50a = fano(x, 110, 15 * 4, 0.1) ;
f70a = fano(x, 80, 15 * 8, 0.1) ;


int1 = sum(f1a);
int30 = sum(f30a);
int50 = sum(f50a);
int70 = sum(f70a);

f1a  = f1a ;
f30a = f30a * int1 / int30 * 0.65;
f50a = f50a * int1 / int50 * 0.65^2;
f70a = f70a * int1 / int70 *0.65^3;

plot(x, f1a, label="T=1K")
plot(x, f30a, label="T=30K")
plot(x, f50a, label="T=50K")
plot(x, f70a, label="T=70K")

legend(fontsize=12)
xlabel("Shift cm-1", fontsize=14)
ylabel("Raman (arb units)", fontsize=14)
