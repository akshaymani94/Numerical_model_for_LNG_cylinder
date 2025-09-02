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
hf = 60

C1 = C2 = m1 * cp
B1 = hf * A1
B2 = (K * A1) / L
B3 = ho * Aex
B4 = (K * A2) / L

def dTdt(t,y):
    T1,T2 = y
    dT1dt = (B1 * (Tf - T1) + B2 * (T2 - T1)) / C1 
    dT2dt = (B3 * (To - T2) + B4 * (T1 - T2)) / C2
    return [dT1dt, dT2dt]

y0 = [313, 313]

t_span = (0, 10000)
t_eval = np.linspace(*t_span, 10000)

sol =  solve_ivp(dTdt, t_span, y0, t_eval=t_eval)

plt.plot((sol.t)/60, sol.y[0], label=fr"T1 (inner wall) $h_{{f}} = {hf}$ W/m$^2$")
plt.plot((sol.t)/60, sol.y[1], label="T2 (outer wall)")
plt.xlabel("Time [min]")
plt.ylabel("Temperature [K]")
plt.title("Temperature Change Over Time for Two-Node Model")
plt.ylim(90,330)
plt.yticks(np.arange(90,331,20))
plt.legend()
plt.grid()
plt.savefig("two_node.png",dpi=300)

plt.show()