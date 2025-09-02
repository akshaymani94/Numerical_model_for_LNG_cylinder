from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt

# 100 mm gas layer , 1% gas going out , 5% spray evaporating 

As = 65.973
Aex = 72.3822
Aiex = 87.083

Tf = 110.1
To = 313


cp = 460

ho = 0.0624
hg = 1400
hf = 1400

C1 = 83996000
C2 = 8553.12
B1 = ho * Aiex
B2 = hg * Aex
B3 = hg * Aex
B4 = hf * As

#
def dTdt(t,y):
    Tw,Tg = y
    dT1dt = (B1 * (To - Tw) - B2 * (Tw - Tg)) / C1 
    dT2dt = (B3 * (Tw - Tg) - B4 * (Tg - Tf)  - (0.0851*1005*(Tg - To)) - 95.2) / C2
    return [dT1dt, dT2dt]

y0 = [313, 313]

t_span = (0, 10000)
t_eval = np.linspace(*t_span, 10000)

sol =  solve_ivp(dTdt, t_span, y0, t_eval=t_eval)

plt.plot((sol.t)/60, sol.y[0], label=fr"Tw (wall) $h_{{f}} = {hf}$ W/m$^2$")
plt.plot((sol.t)/60, sol.y[1], label="Tg (gas)")
plt.xlabel("Time [min]")
plt.ylabel("Temperature [K]")
plt.title("Temperature Change Over Time for Gas explicit Model (100mm) gas layer higher hg ")
plt.ylim(90,330)
plt.yticks(np.arange(90,331,20))
plt.legend()
plt.grid()
plt.savefig("gas_explicit_100mm_higher_hg.png",dpi=300)

plt.show()