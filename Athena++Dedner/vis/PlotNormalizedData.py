import h5py
from athena_read import athdf
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import ticker

#Physical Constants
#We will use normalized variables to simplify the computations. Here are the scale constants that will be used as outlined in Hans document:
r0 = 1 #AU
v0 = 52.483 #km/s
t0 = 2.851e6 #s = 33 days
rho0 = 8.361e-27 #kg cm-3 (5 * m_p)
P0 = 2.3e-11 #N/m^2
T0 = 111111 #K

#Read HDF5 data from an Athena++ run
data = athdf('C:/Users/Luke/Desktop/athena/bin/GOODRUNSMR/3DHydrodynamicHelioCartesianSMR.out1.00347.athdf')

# Extract raw data from HDF5
x = np.array(data['x1v']) #AU
y = np.array(data['x2v']) #AU
z = np.array(data['x3v']) #AU
press = np.array(data['press']) 
rho = np.array(data['rho']) 
vel1 = np.array(data['vel1'])
vel2 = np.array(data['vel2']) 

# Save data over stagnation axis (y=0, z=0, x varies)
# Find indices where x2 and x3 are closest to zero (stagnation axis is at these values)
i2_0 = np.argmin(np.abs(y))
i3_0 = np.argmin(np.abs(z))

print(y[i3_0], y[i3_0+1], y[i3_0-1])

# Extract 1D slice along x1 at y=0, z=0 (shape = *[x3,x2,x1])
rho_axis = rho[i3_0, i2_0, :]
press_axis = press[i3_0, i2_0, :]
rad = np.sqrt(x*x + y[i2_0] * y[i2_0] + z[i3_0] * z[i3_0])

# Plot density (rho) data along the stagnation axis. Note rho=1 is the normalized value at 1AU (rho_1AU/rho0 = 1)
plt.figure(figsize=(8,6))
plt.plot(x, rho_axis, label='Rho', marker='o', mfc='none', markersize=4, color='black')
plt.plot(x, 1 / (rad * rad), label='$\\frac{1}{r^2}$ Profile', color='red')
plt.yscale('log')
plt.ylim([0.00001,0.75])
plt.xlabel('x (AU)')
plt.ylabel('Rho (Normalized Units)')
plt.title('Rho Along Stagnation Axis')
plt.legend()
plt.grid(True)
plt.tight_layout()
#plt.savefig('RhoStagAxisFullRange')
plt.close()

plt.figure(figsize=(8,6))
plt.plot(x, rho_axis, label='Rho', marker='o', mfc='none', markersize=4, color='black')
plt.plot(x, 1 / (rad * rad), label='$\\frac{1}{r^2}$ Profile', color='red')
plt.yscale('log')
plt.ylim([0.00001,0.75])
plt.xlim([-300,200])
plt.xlabel('x (AU)')
plt.ylabel('Rho (Normalized Units)')
plt.title('Rho Along Stagnation Axis')
plt.legend()
plt.grid(True)
plt.tight_layout()
#plt.savefig('RhoStagAxisShockRange')
plt.close()

plt.figure(figsize=(8,6))
plt.plot(x, rho_axis, label='Rho', marker='o', mfc='none', markersize=4, color='black')
plt.plot(x, 1 / (rad * rad), label='$\\frac{1}{r^2}$ Profile', color='red')
plt.yscale('log', base=2)
plt.xscale('log', base=2)
#Get current axes, set power of 10 for custom ticks
ax = plt.gca()
y_ticks = np.logspace(-9, 0, num=10, base=10)
x_ticks = np.logspace(0, 3, num=4, base=10)
ax.set_xticks(x_ticks)
ax.set_yticks(y_ticks)
ax.xaxis.set_major_formatter(ticker.LogFormatterMathtext(base=10))
ax.yaxis.set_major_formatter(ticker.LogFormatterMathtext(base=10))
plt.ylim([10e-6,0.75])
plt.xlabel('x (AU, Upstream)')
plt.ylabel('Rho (Normalized Units)')
plt.title('Rho Along Stagnation Axis (log-log, Base 2)')
plt.legend()
plt.grid(True)
#plt.savefig('RhoStagAxisloglog')
plt.close()

plt.figure(figsize=(8,6))
plt.plot(x, rho_axis*rad*rad, label='$\\rho*r^2$', marker='o', mfc='none', markersize=4, color='black')
plt.yscale('log')
plt.xlabel('x (AU)')
plt.ylabel('Rho (Normalized Units)')
plt.title('Rho Along Stagnation Axis')
plt.legend()
plt.grid(True)
#plt.savefig('RhoStagAxisRhoRSquared')
plt.close()


# Plot Thermal Pressure (press) data along the stagnation axis. Note 0.6 is the normalized 1AU value (therm_press_1AU / P0 = 0.6)
plt.figure(figsize=(8,6))
plt.plot(x, press_axis, label='Therm Press', marker='o', mfc='none', markersize=4, color='black')
plt.plot(x, 0.6 * np.power(rad,-3.3333), label='$\\frac{1}{r^{10/3}}$ Profile', color='red')
plt.yscale('log')
plt.ylim([10e-10,1])
plt.xlabel('x (AU)')
plt.ylabel('press (Normalized Units)')
plt.title('Press Along Stagnation Axis')
plt.legend()
plt.grid(True)
plt.tight_layout()
#plt.savefig('PressStagAxisFullRange')
plt.close()

# Plot Thermal Pressure (press) data along the stagnation axis
plt.figure(figsize=(8,6))
plt.plot(x, press_axis, label='Therm Press', marker='o', mfc='none', markersize=4, color='black')
plt.plot(x, 0.6 * np.power(rad,-3.3333), label='$\\frac{1}{r^{10/3}}$ Profile', color='red')
plt.yscale('log')
plt.ylim([10e-10,1])
plt.xlim([-300,200])
plt.xlabel('x (AU)')
plt.ylabel('press (Normalized Units)')
plt.title('Press Along Stagnation Axis')
plt.legend()
plt.grid(True)
plt.tight_layout()
#plt.savefig('PressStagAxisShockRange')
plt.close()

# Plot Thermal Pressure (press) data along the stagnation axis
plt.figure(figsize=(8,6))
plt.plot(x, press_axis, label='Therm Press', marker='o', mfc='none', markersize=4, color='black')
plt.plot(x, 0.6 * np.power(rad,-3.3333), label='$\\frac{1}{r^{10/3}}$ Profile', color='red')
plt.yscale('log', base=10.0/3.0)
plt.xscale('log', base=10.0/3.0)
#Get current axes, set power of 10 for custom ticks
ax = plt.gca()
y_ticks = np.logspace(-9, 0, num=10, base=10)
x_ticks = np.logspace(0, 3, num=4,base=10)
ax.set_xticks(x_ticks)
ax.set_yticks(y_ticks)
ax.xaxis.set_major_formatter(ticker.LogFormatterMathtext(base=10))
ax.yaxis.set_major_formatter(ticker.LogFormatterMathtext(base=10))
plt.ylim([10e-9, 1])
plt.xlim([1,1000])
plt.xlabel('x (AU)')
plt.ylabel('press (Normalized Units)')
plt.title('Press Along Stagnation Axis (log-log, Base 10/3)')
plt.legend()
plt.grid(True)
plt.tight_layout()
#plt.savefig('PressStagAxisLogLog')
plt.close()

# Plot Thermal Pressure (press) data along the stagnation axis
plt.figure(figsize=(8,6))
plt.plot(x, press_axis * np.power(rad,10.0/3.0), label='$Press*r^{\\frac{10}{3}}$', marker='o', mfc='none', markersize=4, color='black')
plt.yscale('log')
#Get Current Axes
plt.ylabel('press (Normalized Units)')
plt.title('Press Along Stagnation Axis')
plt.legend()
plt.grid(True)
plt.tight_layout()
#plt.savefig('PressStagAxisPressR10thirds')
plt.close()