import numpy as np
import matplotlib.pyplot as plt

plt.ion()

path = "/data/amarchal/ROHSA_paper/ROHSA/"
data1 = np.genfromtxt(path+"Tb_reso_0.8km.s-1_Tmin_0_Tmax_inf_ROHSA_noise_0.05_K_beam_0_2_2_gauss_run_35_timestep.dat")
data2 = np.genfromtxt(path+"Tb_reso_0.8km.s-1_Tmin_0_Tmax_inf_ROHSA_noise_0.05_K_beam_0_2_2_gauss_run_36_timestep.dat")
data4 = np.genfromtxt(path+"Tb_reso_0.8km.s-1_Tmin_0_Tmax_inf_ROHSA_noise_0.05_K_beam_0_2_2_gauss_run_38_timestep.dat")
data6 = np.genfromtxt(path+"Tb_reso_0.8km.s-1_Tmin_0_Tmax_inf_ROHSA_noise_0.05_K_beam_0_2_2_gauss_run_40_timestep.dat")
data8 = np.genfromtxt(path+"Tb_reso_0.8km.s-1_Tmin_0_Tmax_inf_ROHSA_noise_0.05_K_beam_0_2_2_gauss_run_42_timestep.dat")

fig=plt.figure(figsize=(10,10))
ax = fig.add_subplot(111)
ax.plot(data1[:,0],data1[:,1]/60., linewidth=2., label="N=1")
ax.plot(data2[:,0],data2[:,1]/60., linewidth=2., label="N=2")
ax.plot(data4[:,0],data4[:,1]/60., linewidth=2., label="N=4")
ax.plot(data6[:,0],data6[:,1]/60., linewidth=2., label="N=6")
ax.plot(data8[:,0],data8[:,1]/60., linewidth=2., label="N=8")
ax.set_xlim([0., 300.])
ax.set_xlabel(r'Grid size', fontsize = 16)
ax.set_ylabel(r'Time (min)', fontsize = 16)
plt.legend(loc = 1, numpoints = 1)
leg = plt.gca().get_legend()
ltext  = leg.get_texts()
plt.setp(ltext, fontsize = 'small')
plt.savefig('plot/timestep.png', format='png', bbox_inches='tight', pad_inches=0.02)
