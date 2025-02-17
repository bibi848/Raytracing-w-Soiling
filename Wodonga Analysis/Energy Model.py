#%%
# Simulating energy demand
import numpy as np
import matplotlib.pyplot as plt

# f = t time space, T = width, beta = flatness, peak = , f0 = time of peak
def raised_cosine(f, T, beta, peak, f0):
    f = f - f0
    if np.abs(f) <= (1-beta)/(2*T):
        return peak
    
    elif np.abs(f) > (1-beta)/(2*T) and np.abs(f) <= (1+beta)/(2*T):
        return peak * 0.5 * (1 + np.cos((np.pi * T/ beta) * (np.abs(f) - (1-beta)/(2*T))))
    
    else:
        return 0

t = np.linspace(0, 24, 288)

baseline_power = 3

random_variability = np.random.normal(0, 0.6, len(t))

power1 = np.array([raised_cosine(ts, T=0.16, beta=0.3, peak=14, f0=8.5) for ts in t])
power2 = np.array([raised_cosine(ts, T=0.2, beta=0.3, peak=13.5, f0=15) for ts in t])

total_power = baseline_power + power1 + power2 + random_variability


# Plot
plt.figure(figsize=(10, 5))
plt.plot(t, total_power)
plt.xlabel("Time of Day (hours)")
plt.ylabel("Power (MW)")
plt.title("Simulated Energy Usage of Plant")
plt.xticks(np.arange(0, 25, 2))
plt.yticks(np.arange(1, 20, 2))
plt.grid()
plt.show()

#%%
# Sampling the period and making an energy demand spreadsheet

energy_demand = []

num_timesteps = 183308
i = 0
u = 0
while i < num_timesteps:

    energy_demand.append(total_power[u]*1000/3.6)

    u += 1
    i += 1

    if u == len(total_power):
        u = 0

import os
import pandas as pd
current_script_path = os.path.abspath(__file__)
wodonga_path = os.path.dirname(current_script_path)

file_path = wodonga_path + "\\Wodonga Data\\Wodonga Soiled Data.csv"
df = pd.read_csv(file_path)
df['Energy Demand [kWh]'] = energy_demand

df.to_csv(file_path, index=False)