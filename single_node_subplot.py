from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt

# SIngle node without insulation

m = 18260
cp = 460
ho = 15         # without insulation
hoi = 0.09     # effective h when insulation is consisdered

hf_values = [60,80,100,500,1000,1050,3000]

# l = 7.8
# ri = 1.6
# rw = 1.616
# ro = 1.816
# Tw = 263
A = 82.48
Tf = 112.25     # -160.9 deg celcius
To = 314.15        # 40 deg celcius


C1 = m*cp
C2 = ho * A

T0 = 314.15
t_span = (0,10000)
t_eval = np.linspace(*t_span, 10000)

times_below_115 = {}

fig, axes = plt.subplots(1,2,figsize=(14,6))

for hf in hf_values:
    C2 = ho * A
    C3 = hf * A
    def dTdt(t,y):
        T = y[0]
        # return (C2 * (T -To) - C3 * (T - Tf)) / C1
        return (C2 * (To -T) + C3 * (Tf - T)) / C1
    
    sol = solve_ivp(dTdt, t_span, [T0], t_eval=t_eval)
    # plt.plot(sol.t, sol.y[0],label=fr"$h_{{f}} = {hf}$ W/m$^2$")


    # find the time when temperature drops below 115K
    below_115 = np.where(sol.y[0] < 115)[0]
    if below_115.size > 0:
        t_below_115 = (sol.t[below_115[0]])/60
        times_below_115[hf] = t_below_115
        label = fr"$h_{{f}} = {hf}$ W/m$^2$, $t_{{115}}$ = {t_below_115:.0f} min"

    else:
        times_below_115[hf] = None
        label = fr"$h_{{f}} = {hf}$ W/m$^2$, $t_{{115}}$ = N/A"
    
    axes[0].plot((sol.t)/60, sol.y[0], label=label)

axes[0].set_xlabel("Time [min]")
axes[0].set_ylabel("Temperature [K]")
axes[0].set_title("Temperature Change Over Time without insulation")
axes[0].set_ylim(90,330)
axes[0].set_yticks(np.arange(90,331,20))
axes[0].grid()
axes[0].legend()


for hf in hf_values:
    C2 = hoi * A
    C3 = hf * A
    def dTdt(t,y):
        T = y[0]
        # return (C2 * (T -To) - C3 * (T - Tf)) / C1
        return (C2 * (To -T) + C3 * (Tf - T)) / C1
    
    sol = solve_ivp(dTdt, t_span, [T0], t_eval=t_eval)
    # plt.plot(sol.t, sol.y[0],label=fr"$h_{{f}} = {hf}$ W/m$^2$")


    # find the time when temperature drops below 115K
    below_115 = np.where(sol.y[0] < 115)[0]
    if below_115.size > 0:
        t_below_115 = (sol.t[below_115[0]])/60
        times_below_115[hf] = t_below_115
        label = fr"$h_{{f}} = {hf}$ W/m$^2$, $t_{{115}}$ = {t_below_115:.0f} min"

    else:
        times_below_115[hf] = None
        label = fr"$h_{{f}} = {hf}$ W/m$^2$, $t_{{115}}$ = N/A"
    
    axes[1].plot((sol.t)/60, sol.y[0], label=label)

axes[1].set_xlabel("Time [min]")
axes[1].set_ylabel("Temperature [K]")
axes[1].set_title("Temperature Change Over Time with insulation")
axes[1].set_ylim(90,330)
axes[1].set_yticks(np.arange(90,331,20))
axes[1].grid()
axes[1].legend()

plt.suptitle("Temperature Change Over Time for Single Node Model",fontweight='bold')
plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.tight_layout()
plt.savefig("single_node_subplot.png",dpi=300)

plt.show()