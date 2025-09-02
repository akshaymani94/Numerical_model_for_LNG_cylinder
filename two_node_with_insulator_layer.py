from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt


       
Ain = 72.3822       # interior
Aiex = 87.083       # insulation exterior

Tf = 112.25
To = 314.15


ho = 15
hf_values = [60,70,100,500,1000,3000]
Riw = 0.128

Cin = 1211612.057
Cw = 8399600


B1 = ho * Aiex
B2 = 1/Riw
B3 = 1/Riw



cmap = plt.get_cmap('tab10')
colors = [cmap(i) for i in range(len(hf_values))]


for idx,hf in enumerate(hf_values):
    B4 = hf * Ain
    def dTdt(t,y):
        Tin,Tw = y
        dT1dt = (B1 * (To - Tin) - B2 * (Tin - Tw)) / Cin 
        dT2dt = (B3 * (Tin - Tw) - B4 * (Tw - Tf)) / Cw 
        return [dT1dt, dT2dt]

    y0 = [314.15, 314.15]

    t_span = (0, 10000)
    t_eval = np.linspace(*t_span, 10000)

    times_below_115 = {}

    sol =  solve_ivp(dTdt, t_span, y0, t_eval=t_eval)
    # find the time when temperature drops below 111K
    below_115 = np.where(sol.y[1] < 115)[0]
    if below_115.size > 0:
        t_below_115 = (sol.t[below_115[0]])/60
        times_below_115[hf] = t_below_115
        text1 = fr"$h_{{f}} = {hf}$ W/m$^2$, $t_{{115}}$ = {t_below_115:.0f} min"

    else:
        times_below_115[hf] = None
        text1 = fr"$h_{{f}} = {hf}$ W/m$^2$, $t_{{115}}$ = N/A"

    plt.plot((sol.t)/60, (sol.y[0]), label="Tin (insulation) " + text1, linestyle=':', linewidth=2.5, color=colors[idx] )
    plt.plot((sol.t)/60, sol.y[1], label=fr"Tw (wall)", linestyle='-', linewidth=1, color=colors[idx])


# plt.text(0.65, 0.8, 
#          text1,
#          transform=plt.gca().transAxes,
#          fontsize=10,
#          verticalalignment='top',
#          bbox=dict(boxstyle="round", facecolor="#90ee90", alpha=1))


plt.xlabel("Time [min]")
plt.ylabel("Temperature [K]")
plt.title("Temperature Change Over Time for Two-node Model with Insulator Layer")
plt.ylim(90,330)
plt.yticks(np.arange(90,331,20))
plt.legend()
plt.grid()
plt.savefig("two_node_insulator.png",dpi=300)

plt.show()