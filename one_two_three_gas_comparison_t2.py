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

fig, axes = plt.subplots(2,2,figsize=(14,10))

# First plot single node with insulation

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
        label = fr"$h_{{f}} = {hf}$ W/m$^2$, $t_{{115}}$ = {t_below_115:.0f} min , $T_{{lowest}}$ = {(np.min(sol.y[0])-273):.0f} $^oC$"  

    else:
        times_below_115[hf] = None
        label = fr"$h_{{f}} = {hf}$ W/m$^2$, $t_{{115}}$ = N/A , $T_{{lowest}}$ = {(np.min(sol.y[0])-273):.0f} $^oC$"  
    
    axes[0,0].plot((sol.t)/60, sol.y[0], label=label)

axes[0,0].set_xlabel("Time [min]")
axes[0,0].set_ylabel("Temperature [K]")
axes[0,0].set_title("One-Node Model")
axes[0,0].set_ylim(90,330)
axes[0,0].set_yticks(np.arange(90,331,20))
axes[0,0].grid()
axes[0,0].legend()




# Second plot two node


A1 = 72.38
A2 = 73.43


K = 45
L = 0.008  # 16 mm steel wall thickness

Tf = 112.25
To = 314.15

m1 = 9130
m2 = 9130
cp = 460

ho = 0.09
hf_values = [60,70,80,100,500,1000,3000]

C1 = C2 = m1 * cp

B2 = (K * A1) / L
B3 = ho * A2
B4 = (K * A2) / L

y0 = [314.15, 314.15]

t_span = (0, 10000)
t_eval = np.linspace(*t_span, 10000)

times_below_115 = {}

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
    below_115 = np.where(sol.y[0] < 115)[0]
    if below_115.size > 0:
        t_below_115 = (sol.t[below_115[0]])/60
        times_below_115[hf] = t_below_115
        label_T1 = fr"T1 (inner wall) $h_{{f}} = {hf}$ W/m$^2$, $t_{{115}}$ = {t_below_115:.0f} min, $T_{{lowest}}$ = {(np.min(sol.y[0])-273):.0f} $^oC$"  

    else:
        times_below_115[hf] = None
        label = fr"$h_{{f}} = {hf}$ W/m$^2$, $t_{{115}}$ = N/A, $T_{{lowest}}$ = {(np.min(sol.y[0])-273):.0f} $^oC$"
    label_T2 = fr"T2 (outer wall) $h_{{f}} = {hf}$ W/m$^2$"
    axes[0,1].plot((sol.t)/60, sol.y[0], label=label_T1, linestyle=':', linewidth=2.5, color=colors[idx])
    axes[0,1].plot((sol.t)/60, sol.y[1], label=label_T2, linestyle='-', linewidth=1, color=colors[idx])


# plt.plot(sol.t, sol.y[0], label="T1 (inner wall)")
# plt.plot(sol.t, sol.y[1], label="T2 (outer wall)")
axes[0,1].set_xlabel("Time [min]")
axes[0,1].set_ylabel("Temperature [K]")
axes[0,1].set_title("Two-Node Model")
axes[0,1].set_ylim(90,330)
axes[0,1].set_yticks(np.arange(90,331,20))
axes[0,1].legend(fontsize=8)
axes[0,1].grid()



# three node model (gas explicit)

# 80% gas volume , mass balance of spray and gas 
# complete spray evaporated till Tg reaches 113K

As = 16.066
Aex = 72.3822
Aiex = 87.083

Tf = 112.25
To = 314.15



ho = 0.0624
hf = 60
hg = 60

C1 = 8399600
C2 = 46955.43
B1 = ho * Aiex
B2 = hg * Aex
B3 = hg * Aex
B4 = hf * As


def dTdt(t,y):
    Tw,Tg = y
    if (y[1] >= 113.15):
        
        dT1dt = (B1 * (To - Tw) - B2 * (Tw - Tg)) / C1 
        dT2dt = (B3 * (Tw - Tg) - B4 * (Tg - Tf)  - (1.5922*1005*(Tg - To)) - ((18.39+(204.128*hg))*510000)) / C2
    else:
    
        dT1dt = (B1 * (To - Tw) - B2 * (Tw - Tg)) / C1 
        dT2dt = (B3 * (Tw - Tg) - B4 * (Tg - Tf)  - (1.5922*1005*(Tg - To))) / C2
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
    text1 = fr"$h_{{f}} = {hf}$ W/m$^2$, $t_{{115}}$ = {t_below_115:.0f} min ,$T_{{w,lowest}}$ = {(np.min(sol.y[0])-273):.0f} $^oC$"  

else:
    times_below_115[hf] = None
    text1 = fr"$h_{{f}} = {hf}$ W/m$^2$, $t_{{115}}$ = N/A ,$T_{{w,lowest}}$ = {(np.min(sol.y[0])-273):.0f} $^oC$" 

axes[1,0].plot((sol.t)/60, sol.y[0], label=fr"Tw (wall)")
axes[1,0].plot((sol.t)/60, sol.y[1], label="Tg (gas)")

axes[1,0].text(0.3, 0.8, 
         text1,
         transform=axes[1,0].transAxes,
         fontsize=10,
         verticalalignment='top',
         bbox=dict(boxstyle="round", facecolor="#90ee90", alpha=1))


# Add notes on the right side (x=0.99 is near the right edge, y=0.5 is vertical center)
axes[1,0].text(0.3, 0.7,
    "Assumptions: \n"
    "1:  80% Volume occupied by gas, 20% by spray\n"
    "2:  Mass of spray entering = Mass of air leaving\n"
    "3:  $m_{spray}$ evaporated calculated through energy balance\n"
    "    $Q_{req} = Q_{avail}$ upto $T_{g} \\geq 160$ $^\\circ$C\n",
         transform=axes[1,0].transAxes,
         fontsize=10,
         verticalalignment='top',bbox=dict(boxstyle="round", facecolor="#f0f0f0", alpha=0.9)
)


axes[1,0].set_xlabel("Time [min]")
axes[1,0].set_ylabel("Temperature [K]")
axes[1,0].set_title("Gas explicit Model (80 percent vol)")
axes[1,0].set_ylim(90,330)
axes[1,0].set_yticks(np.arange(90,331,20))
axes[1,0].legend()
axes[1,0].grid()


# three node model (with insulation explicitly)
# mass balance of spray and gas
# complete spray evaporated till Tg reaches 113K


As = 16.066         
Ain = 72.3822       # interior
Aiex = 87.083       # insulation exterior

Tf = 112.25
To = 314.15


ho = 15
hg = 60
hf = 60
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
    if(y[2] >= 113.15):
        
        dT1dt = (B1 * (To - Tin) - B2 * (Tin - Tw)) / Cin 
        dT2dt = (B3 * (Tin - Tw) - B4 * (Tw - Tg)) / Cw 
        dT3dt = (B5 * (Tw - Tg) - B6 * (Tg - Tf)  - (1.5922*1005*(Tg - To)) - ((18.39+(204.128*hg))*510000)) / Cg
    else:
        
        dT1dt = (B1 * (To - Tin) - B2 * (Tin - Tw)) / Cin 
        dT2dt = (B3 * (Tin - Tw) - B4 * (Tw - Tg)) / Cw 
        dT3dt = (B5 * (Tw - Tg) - B6 * (Tg - Tf)  - (1.5922*1005*(Tg - To))) / Cg
    return [dT1dt, dT2dt, dT3dt]

y0 = [314.15, 314.15, 314.15]

t_span = (0, 10000)
t_eval = np.linspace(*t_span, 10000)

times_below_115 = {}

sol =  solve_ivp(dTdt, t_span, y0, t_eval=t_eval)
# find the time when temperature drops below 111K
below_115 = np.where(sol.y[1] < 115)[0]
if below_115.size > 0:
    t_below_115 = (sol.t[below_115[0]])/60
    times_below_115[hf] = t_below_115
    text1 = fr"$h_{{f}} = {hf}$ W/m$^2$, $t_{{115}}$ = {t_below_115:.0f} min, $T_{{w,lowest}}$ = {(np.min(sol.y[1])-273):.0f} $^oC$"

else:
    times_below_115[hf] = None
    text1 = fr"$h_{{f}} = {hf}$ W/m$^2$, $t_{{115}}$ = N/A, $T_{{w,lowest}}$ = {(np.min(sol.y[1])-273):.0f} $^oC$"

axes[1,1].plot((sol.t)/60, (sol.y[0]), label="Tin (insulation)")
axes[1,1].plot((sol.t)/60, sol.y[1], label=fr"Tw (wall)")
axes[1,1].plot((sol.t)/60, sol.y[2], label="Tg (gas)")

axes[1,1].text(0.3, 0.9, 
         text1,
         transform=axes[1,1].transAxes,
         fontsize=10,
         verticalalignment='top',
         bbox=dict(boxstyle="round", facecolor="#90ee90", alpha=1))

# Add notes on the right side (x=0.99 is near the right edge, y=0.5 is vertical center)
axes[1,1].text(0.3, 0.8,
    "Assumptions: \n"
    "1:  80% Volume occupied by gas, 20% by spray\n"
    "2:  Mass of spray entering = Mass of air leaving\n"
    "3:  $m_{spray}$ evaporated calculated through energy balance\n"
    "    $Q_{req} = Q_{avail}$ upto $T_{g} \\geq 160$ $^\\circ$C\n",
         transform=axes[1,1].transAxes,
         fontsize=10,
         verticalalignment='top',bbox=dict(boxstyle="round", facecolor="#f0f0f0", alpha=0.9)
)

axes[1,1].set_xlabel("Time [min]")
axes[1,1].set_ylabel("Temperature [K]")
axes[1,1].set_title("Three-node Model")
axes[1,1].set_ylim(90,330)
axes[1,1].set_yticks(np.arange(90,331,20))
axes[1,1].legend(
    loc='lower right',   # places legend inside at right bottom
    fontsize=8,
    ncol=1
)
axes[1,1].grid()



plt.suptitle("Temperature Change Over Time",fontweight='bold')
plt.tight_layout(rect=[0, 0.03, 1, 0.95])
# plt.tight_layout()
plt.savefig("Comparison_of_different_models_t2.png",dpi=300)

plt.show()