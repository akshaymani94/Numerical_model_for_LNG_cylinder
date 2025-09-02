from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt

# 80% gas volume , 1% gas going out , 5% spray evaporating 

As = 16.066         
Ain = 72.3822       # interior
Aiex = 87.083       # insulation exterior

Tf = 110.1
To = 313


ho = 15
hg = 500
hf = 500
Riw = 0.128

Cin = 1211612.057
Cw = 8399600
Cg = 46955.43

B1 = ho * Aiex
B2 = 1/Riw
B3 = 1/Riw
B4 = hg * Ain
B5 = hg * Ain
B6 = hf * As

def dTdt(t,y):
    Tin,Tw,Tg = y
    dT1dt = (B1 * (To - Tin) - B2 * (Tin - Tw)) / Cin 
    dT2dt = (B3 * (Tin - Tw) - B4 * (Tw - Tg)) / Cw 
    dT3dt = (B5 * (Tw - Tg) - B6 * (Tg - Tf)  - (0.4672*1005*(Tg - To)) - 95.2) / Cg
    return [dT1dt, dT2dt, dT3dt]

y0 = [313, 313, 313]

t_span = (0, 10000)
t_eval = np.linspace(*t_span, 10000)

times_below_111 = {}

sol =  solve_ivp(dTdt, t_span, y0, t_eval=t_eval)
# find the time when temperature drops below 111K
below_111 = np.where(sol.y[1] < 111)[0]
if below_111.size > 0:
    t_below_111 = (sol.t[below_111[0]])/60
    times_below_111[hf] = t_below_111
    text1 = fr"$h_{{f}} = {hf}$ W/m$^2$, $t_{{111}}$ = {t_below_111:.0f} min"

else:
    times_below_111[hf] = None
    text1 = fr"$h_{{f}} = {hf}$ W/m$^2$, $t_{{111}}$ = N/A"

plt.plot((sol.t)/60, (sol.y[0]), label="Tin (insulation)")
plt.plot((sol.t)/60, sol.y[1], label=fr"Tw (wall)")
plt.plot((sol.t)/60, sol.y[2], label="Tg (gas)")

plt.text(0.65, 0.8, 
         text1,
         transform=plt.gca().transAxes,
         fontsize=10,
         verticalalignment='top',
         bbox=dict(boxstyle="round", facecolor="#90ee90", alpha=1))


plt.xlabel("Time [min]")
plt.ylabel("Temperature [K]")
plt.title("Temperature Change Over Time for three node Model")
plt.ylim(90,330)
plt.yticks(np.arange(90,331,20))
plt.legend()
plt.grid()
plt.savefig("three_node_model.png",dpi=300)

plt.show()