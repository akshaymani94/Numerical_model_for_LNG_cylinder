from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt

# 80% gas volume , 1% gas going out , 5% spray evaporating 

As = 16.066
Aex = 72.3822
Aiex = 87.083

Tf = 112.25
To = 314.15



ho = 0.0624
hf = 500
hg = 500

C1 = 8399600
C2 = 46955.43
B1 = ho * Aiex
B2 = hg * Aex
B3 = hg * Aex
B4 = hf * As


def dTdt(t,y):
    Tw,Tg = y
    dT1dt = (B1 * (To - Tw) - B2 * (Tw - Tg)) / C1 
    dT2dt = (B3 * (Tw - Tg) - B4 * (Tg - Tf)  - (0.46721*1005*(Tg - To)) - 95.2) / C2
    return [dT1dt, dT2dt]

y0 = [314.15, 314.15]

t_span = (0, 10000)
t_eval = np.linspace(*t_span, 10000)

times_below_115 = {}

sol =  solve_ivp(dTdt, t_span, y0, t_eval=t_eval)
# find the time when temperature drops below 111K
below_115 = np.where(sol.y[0] < 115)[0]
if below_115.size > 0:
    t_below_115 = (sol.t[below_115[0]])/60
    times_below_115[hf] = t_below_115
    text1 = fr"$h_{{f}} = {hf}$ W/m$^2$, $t_{{115}}$ = {t_below_115:.0f} min"

else:
    times_below_115[hf] = None
    text1 = fr"$h_{{f}} = {hf}$ W/m$^2$, $t_{{115}}$ = N/A"

plt.plot((sol.t)/60, sol.y[0], label=fr"Tw (wall)")
plt.plot((sol.t)/60, sol.y[1], label="Tg (gas)")

plt.text(0.65, 0.8, 
         text1,
         transform=plt.gca().transAxes,
         fontsize=10,
         verticalalignment='top',
         bbox=dict(boxstyle="round", facecolor="#90ee90", alpha=1))


plt.xlabel("Time [min]")
plt.ylabel("Temperature [K]")
plt.title("Temperature Change Over Time for Gas explicit Model 80percent vol higher hg ")
plt.ylim(90,330)
plt.yticks(np.arange(90,331,20))
plt.legend()
plt.grid()
plt.savefig("gas_explicit_80vol_higher_hg.png",dpi=300)

plt.show()