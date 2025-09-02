from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt

m = 18260
cp = 460
# ho = 15         # without insulation
ho = 0.09     # effective h when insulation is consisdered

hf_values = [60,80,100,500,1000,3000]

# l = 7.8
# ri = 1.6
# rw = 1.616
# ro = 1.816
# Tw = 263
A = 110.584
Tf = 110.1
To = 313


C1 = m*cp
C2 = ho * A

T0 = 313
t_span = (0,10000)
t_eval = np.linspace(*t_span, 10000)

times_below_111 = {}

for hf in hf_values:
    C3 = hf * A
    def dTdt(t,y):
        T = y[0]
        # return (C2 * (T -To) - C3 * (T - Tf)) / C1
        return (C2 * (To -T) + C3 * (Tf - T)) / C1
    
    sol = solve_ivp(dTdt, t_span, [T0], t_eval=t_eval)
    # plt.plot(sol.t, sol.y[0],label=fr"$h_{{f}} = {hf}$ W/m$^2$")


    # find the time when temperature drops below 115K
    below_111 = np.where(sol.y[0] < 111)[0]
    if below_111.size > 0:
        t_below_111 = (sol.t[below_111[0]])/60
        times_below_111[hf] = t_below_111
        label = fr"$h_{{f}} = {hf}$ W/m$^2$, $t_{{111}}$ = {t_below_111:.0f} min"

    else:
        times_below_111[hf] = None
        label = fr"$h_{{f}} = {hf}$ W/m$^2$, $t_{{111}}$ = N/A"
    
    plt.plot((sol.t)/60, sol.y[0], label=label)




plt.xlabel("Time [min]")
plt.ylabel("Temperature [K]")
plt.title("Temperature Change Over Time for Single Node Model")
plt.ylim(90,330)
plt.yticks(np.arange(90,331,20))

plt.grid()
plt.legend()

plt.savefig("single_node_variable_hf.png",dpi=300)

plt.show()