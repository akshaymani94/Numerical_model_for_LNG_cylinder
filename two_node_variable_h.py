from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt

A1 = 104.55
A2 = 105.922
Aex = 123.5959

K = 45
L = 0.008  # 16 mm steel wall thickness

Tf = 110.1
To = 313

m1 = 9130
m2 = 9130
cp = 460

ho = 0.0624
hf_values = [60,80,100,500,1000,3000]

C1 = C2 = m1 * cp

B2 = (K * A1) / L
B3 = ho * Aex
B4 = (K * A2) / L

y0 = [313, 313]

t_span = (0, 10000)
t_eval = np.linspace(*t_span, 10000)

times_below_111 = {}

cmap = plt.get_cmap('tab10')
colors = [cmap(i) for i in range(len(hf_values))]


for idx,hf in enumerate(hf_values):
    B1 = hf * A1
    def dTdt(t,y):
    
        T1,T2 = y
        dT1dt = (B1 * (Tf - T1) + B2 * (T2 - T1)) / C1 
        dT2dt = (B3 * (To - T2) + B4 * (T1 - T2)) / C2
        return [dT1dt, dT2dt]

    sol =  solve_ivp(dTdt, t_span, y0, t_eval=t_eval)
    below_111 = np.where(sol.y[0] < 111)[0]
    if below_111.size > 0:
        t_below_111 = (sol.t[below_111[0]])/60
        times_below_111[hf] = t_below_111
        label_T1 = fr"T1 (inner wall) $h_{{f}} = {hf}$ W/m$^2$, $t_{{111}}$ = {t_below_111:.0f} min"

    else:
        times_below_111[hf] = None
        label = fr"$h_{{f}} = {hf}$ W/m$^2$, $t_{{111}}$ = N/A"
    label_T2 = fr"T2 (outer wall) $h_{{f}} = {hf}$ W/m$^2$"
    plt.plot((sol.t)/60, sol.y[0], label=label_T1, linestyle=':', linewidth=2.5, color=colors[idx])
    plt.plot((sol.t)/60, sol.y[1], label=label_T2, linestyle='-', linewidth=1, color=colors[idx])


# plt.plot(sol.t, sol.y[0], label="T1 (inner wall)")
# plt.plot(sol.t, sol.y[1], label="T2 (outer wall)")
plt.xlabel("Time [min]")
plt.ylabel("Temperature [K]")
plt.title("Temperature Change Over Time for Two-Node Model")
plt.ylim(90,330)
plt.yticks(np.arange(90,331,20))
plt.legend()
plt.grid()
plt.savefig("two_node_variable_hf.png",dpi=300)

plt.show()