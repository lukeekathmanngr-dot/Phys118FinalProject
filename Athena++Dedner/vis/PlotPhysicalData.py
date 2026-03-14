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
data = athdf('C:/Users/Luke/Desktop/athena/bin/3DHydrodynamicHelioCartesianSMRSphericalBC.out1.00039.athdf')
# Extract raw data from HDF5
x = np.array(data['x1v']) #AU
y = np.array(data['x2v']) #AU
z = np.array(data['x3v']) #AU
press = np.array(data['press']) * P0 
dens = np.array(data['rho']) * 5
temp = (np.array(data['press']) / np.array(data['rho'])) * T0 * (1/(5.0/3.0 - 1)) #gamma-1
vel1 = np.array(data['vel1']) * v0
vel2 = np.array(data['vel2']) * v0

# Save data over stagnation axis (y=0, z=0, x varies)
# Find indices where x2 and x3 are closest to zero (stagnation axis is at these values)
i2_0 = np.argmin(np.abs(y))
i3_0 = np.argmin(np.abs(z))

# Extract 1D slice along x1 at y=0, z=0 (shape = *[x3,x2,x1])
dens_axis = dens[i3_0, i2_0, :]
press_axis = press[i3_0, i2_0, :]
temp_axis = temp[i3_0, i2_0, :]

#Calculate radial values along the line we are plotting over (y and z are not quite exactly zero.)
rad = np.sqrt(x*x + y[i2_0] * y[i2_0] + z[i3_0] * z[i3_0])

# Plot density data along the stagnation axis. Note dens=5 (cm-3) is the boundary value at 1AU
plt.figure(figsize=(8,6))
plt.plot(x, dens_axis, label='Dens', marker='o', mfc='none', markersize=4, color='black')
plt.plot(x, 5 / (rad * rad), label='$5*\\frac{1}{r^2}$ Profile', color='red')
plt.yscale('log')
#plt.ylim([1*pow(10,-5),5])
plt.xlabel('x (AU)')
plt.ylabel('Dens (cm-3)')
plt.title('Dens Along Stagnation Axis')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('DensStagAxisFullRange')
plt.close()

plt.figure(figsize=(8,6))
plt.plot(x, dens_axis, label='Dens', marker='o', mfc='none', markersize=4, color='black')
plt.plot(x, 5 / (rad * rad), label='$5*\\frac{1}{r^2}$ Profile', color='red')
plt.yscale('log')
#plt.ylim([1*pow(10,-5),5])
plt.xlim([-300,200])
plt.xlabel('x (AU)')
plt.ylabel('Dens (cm-3)')
plt.title('Dens Along Stagnation Axis')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('DensStagAxisShockRange')
plt.close()

plt.figure(figsize=(8,6))
plt.plot(x, dens_axis, label='Dens', marker='o', mfc='none', markersize=4, color='black')
plt.plot(x, 5 / (rad * rad), label='$5*\\frac{1}{r^2}$ Profile', color='red')
plt.yscale('log', base=2)
plt.xscale('log', base=2)
#Get current axes, set power of 10 for custom ticks
ax = plt.gca()
y_ticks = np.logspace(-4, 0, num=5, base=10)
x_ticks = np.logspace(0, 3, num=4, base=10)
ax.set_xticks(x_ticks)
ax.set_yticks(y_ticks)
ax.xaxis.set_major_formatter(ticker.LogFormatterMathtext(base=10))
ax.yaxis.set_major_formatter(ticker.LogFormatterMathtext(base=10))
ax.autoscale(enable=True, axis='both', tight=True) 
plt.xlabel('x (AU, Upstream)')
plt.ylabel('Dens (cm-3)')
plt.title('Dens Along Stagnation Axis (log-log, Base 2)')
plt.legend()
plt.grid(True)
plt.savefig('DensStagAxisloglog')
plt.close()

plt.figure(figsize=(8,6))
plt.plot(x, dens_axis*rad*rad, label='$Dens*r^2$', marker='o', mfc='none', markersize=4, color='black')
plt.yscale('log')
plt.xlabel('x (AU)')
plt.ylabel('Dens (cm-3)')
plt.title('Dens Along Stagnation Axis')
plt.legend()
plt.grid(True)
plt.savefig('DensStagAxisDensRSquared')
plt.close()


# Plot Temp data along the stagnation axis. Note T = 111111K (press/rho) * (1/(gamma-1)) for the normalized constants we have picked
plt.figure(figsize=(8,6))
plt.plot(x, temp_axis, label='T', marker='o', mfc='none', markersize=4, color='black')
plt.plot(x, 100000 * np.power(rad,-1.3333), label='$100000*\\frac{1}{r^{4/3}}$ Profile', color='red')
plt.yscale('log')
#plt.ylim([10,5*pow(10,6)])
plt.xlabel('x (AU)')
plt.ylabel('T (K)')
plt.title('T Along Stagnation Axis')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('TStagAxisFullRange')
plt.close()

plt.figure(figsize=(8,6))
plt.plot(x, temp_axis, label='T', marker='o', mfc='none', markersize=4, color='black')
plt.plot(x, 100000 * np.power(rad,-1.3333), label='$100000*\\frac{1}{r^{4/3}}$ Profile', color='red')
plt.yscale('log')
plt.xlim([-300,200])
#plt.ylim([10,5*pow(10,6)])
plt.xlabel('x (AU)')
plt.ylabel('T (K)')
plt.title('T Along Stagnation Axis')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('TStagAxisShockRange')
plt.close()

plt.figure(figsize=(8,6))
plt.plot(x, temp_axis, label='T', marker='o', mfc='none', markersize=4, color='black')
plt.plot(x, 100000 * np.power(rad,-1.3333), label='$100000*\\frac{1}{r^{4/3}}$ Profile', color='red')
plt.yscale('log', base=4.0/3.0)
plt.xscale('log', base=4.0/3.0)
#Get current axes, set power of 10 for custom ticks
ax = plt.gca()
y_ticks = np.logspace(2, 6, num=5, base=10)
x_ticks = np.logspace(0, 3, num=4, base=10)
ax.set_xticks(x_ticks)
ax.set_yticks(y_ticks)
ax.xaxis.set_major_formatter(ticker.LogFormatterMathtext(base=10))
ax.yaxis.set_major_formatter(ticker.LogFormatterMathtext(base=10))
#plt.xlim([1,1000])
#plt.ylim([10e1,5*pow(10,6)])
plt.xlabel('x (AU)')
plt.ylabel('T (K)')
plt.title('T Along Stagnation Axis (log-log, Base 4/3)')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('TStagAxisLogLog')
plt.close()

# Plot Thermal Pressure (press) data along the stagnation axis
plt.figure(figsize=(8,6))
plt.plot(x, temp_axis * np.power(rad,4.0/3.0), label='$T*r^{\\frac{4}{3}}$', marker='o', mfc='none', markersize=4, color='black')
plt.yscale('log')
plt.ylabel('T (K)')
plt.title('T Along Stagnation Axis')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('TStagAxisPressR10thirds')
plt.close()


#Make dual axis temp/dens plot to compare to Hans
#Note: fig controls figure-level properties: size, overall title, saving the whole figure, spacing between subplots.
#Note: ax controls axes-level plotting: labels, limits, plotting data, legends for this specific subplot.

fig, ax1 = plt.subplots(figsize=(10, 8))
ax2 = ax1.twinx()
#Plot Data
ax1.plot(x, temp_axis, label='T', marker='o', mfc='none', markersize=4, color='red')
ax2.plot(x, dens_axis, label='Dens', marker='o', mfc='none', markersize=4, color='blue')
fig.suptitle('Density, Temperature Along Stagnation Axis \n(Log-Log, Log Base T=4/3, Dens=2)')
ax1.set_xlabel('x (AU, Upwind)')
ax1.set_ylabel('T (K)')
ax2.set_ylabel('Dens (cm-3)')
ax1.set_yscale('log', base=4.0/3.0)
ax1.set_xscale('log', base=4.0/3.0)
ax2.set_yscale('log', base=2.0)
ax2.set_xscale('log', base=2.0)
#Set T axes to scale by base 10 log
y1_ticks = np.logspace(2, 7, num=7, base=10)
x1_ticks = np.logspace(0, 3, num=4, base=10)
ax1.set_xticks(x_ticks, minor=False)
ax1.set_yticks(y_ticks, minor=False)
ax1.xaxis.set_major_formatter(ticker.LogFormatterMathtext(base=10))
ax1.yaxis.set_major_formatter(ticker.LogFormatterMathtext(base=10))
ax1.set_ylim([temp_axis[x>0].min(), temp_axis[x>0].max()])
ax1.autoscale(enable=True, axis='x') 
#Set dens axes to scale by base 10 log
y2_ticks = np.logspace(-4, 1, num=6, base=10)
x2_ticks = np.logspace(0, 3, num=4, base=10)
ax2.set_xticks(x2_ticks, minor=False)
ax2.set_yticks(y2_ticks, minor=False)
ax2.xaxis.set_major_formatter(ticker.LogFormatterMathtext(base=10))
ax2.yaxis.set_major_formatter(ticker.LogFormatterMathtext(base=10))
ax2.set_ylim([dens_axis[x>0].min(), dens_axis[x>0].max()])
ax2.autoscale(enable=True, axis='x') 
fig.legend()
fig.tight_layout()
plt.savefig('TDensDualAxisLogLog')