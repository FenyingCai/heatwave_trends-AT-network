
import numpy as np
import xarray as xr


##--##--##--##--  networks (concurrent days)  --##--##--##--##
ds1 = xr.open_dataset("data/PDF_heathigh_daily_1_5std.nc")
concurrent_heat_high_1_5std = ds1.pdf.data

ds2 = xr.open_dataset("data/PDF_heatlow_daily_1_5std.nc")
concurrent_heat_low_1_5std = ds2.pdf.data

ds3 = xr.open_dataset("data/PDF_heathigh_interannualR_1_5std.nc")
IntensityR_heat_high_1_5std = ds3.cor_pdf.data

ds4 = xr.open_dataset("data/PDF_heatlow_interannualR_1_5std.nc")
IntensityR_heat_low_1_5std = ds4.cor_pdf.data



ds5 = xr.open_dataset("data/PDF_heathigh_daily_90th.nc")
concurrent_heat_high_90th = ds5.pdf.data

ds6 = xr.open_dataset("data/PDF_heatlow_daily_10th.nc")
concurrent_heat_low_10th = ds6.pdf.data

ds7 = xr.open_dataset("data/PDF_heathigh_interannualR_90th.nc")
IntensityR_heat_high_90th = ds7.cor_pdf.data

ds8 = xr.open_dataset("data/PDF_heatlow_interannualR_10th.nc")
IntensityR_heat_low_10th = ds8.cor_pdf.data








## Plotting Import
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.use('agg')

## Define a Map Function for plotting
def create_map(ax):
  # plt.axes( xscale="log" )
  plt.xlim(0,9)
  plt.ylim(0.0,0.32)
  plt.xticks(np.linspace(0.0,9,10))
  # plt.yticks(np.linspace(0,0.3,7)); ax.set_yticklabels([0,5,10,15,20,25,30])
  # ax.set_yticklabels(np.linspace(0,30,7))
  ax.tick_params(axis='both', length=9.0, width=1.4, labelsize=15)
  ax.spines['bottom'].set_linewidth(1.1)
  ax.spines['top'].set_linewidth(1.1)
  ax.spines['left'].set_linewidth(1.1)
  ax.spines['right'].set_linewidth(1.1)
  return ax


##--##--##--  create figure  --##--##--##
fig = plt.figure(dpi=400,figsize=(12,12))




ax1 = plt.axes([0.0,0.60, 0.4,0.22])
create_map(ax1)
ax1.set_xlabel(' ', size=18)
ax1.set_ylabel('PDF', size=18)
ax1.set_ylim(0.0,0.04)
ax1.set_yticks([0,0.005,0.01,0.015,0.02,0.025,0.03,0.035,0.040]); ax1.set_yticklabels([0,"",0.01,"",0.02,"",0.03,"",0.04], color='r')
ax1.spines['left'].set_color('r')
ax1.tick_params(axis='y', right=False, length=9.0, width=1.4, labelsize=15, color='r')
ax1.plot(np.linspace(0.6,10.2,52-3), concurrent_heat_high_1_5std[3:], label='heat-high (>40, p<0.01)', color='r', marker='o', linewidth=0.6,  markersize=6.0, linestyle='-')
ax1.legend(loc=(0.408, 0.835), fontsize=12.5)
plt.title("(a)  T2m  &  H500   (concurrent days)", size=18)

ax2 = ax1.twinx()
create_map(ax2)
ax2.set_ylim(0.0,0.18)
ax2.set_yticks([0,0.02,0.04,0.06,0.08,0.10,0.12,0.14,0.16,0.18]); ax2.set_yticklabels([0,"",0.04,"",0.08,"",0.12,"",0.16,""], color='b')
ax2.spines['left'].set_color('r')
ax2.spines['right'].set_color('b')
ax2.tick_params(axis='y', left=False, length=9.0, width=1.4, labelsize=15, color='b')
ax2.plot(np.linspace(0.6,10.2,52-3), concurrent_heat_low_1_5std[3:], label='heat-low (>40, p<0.01)', color='b', markerfacecolor='none', markeredgecolor='b', markeredgewidth=1.2, marker='o', linewidth=0.85,  markersize=7.0, linestyle='--')
ax2.plot([1.5,1.5], [0.0,0.085], label='', color='g', linewidth=1.8,  linestyle='--')
ax2.plot([4.0,4.0], [0.0,0.039], label='', color='g', linewidth=1.8,  linestyle='--')
ax2.legend(loc=(0.425, 0.675), fontsize=12.5)

ax1.text(1.70, 0.025, s='~2600 km', fontsize=12.5, color='b')
ax1.text(4.75, 0.0195, s='~5700 km', fontsize=12.5, color='r')
ax1.text(3.0, 0.048, s='Pressure   extremes   (H500 anomaly > 1.5 std,  < -1.5 std)', fontsize=20.0, color='black')



ax3 = plt.axes([0.52,0.60, 0.4,0.22])
create_map(ax3)
ax3.set_xlabel(' ', size=18)
ax3.set_ylabel('', size=18)
ax3.set_ylim(0.0,0.04)
ax3.set_yticks([0,0.005,0.01,0.015,0.02,0.025,0.03,0.035,0.04]); ax3.set_yticklabels([0,"",0.01,"",0.02,"",0.03,"",0.04], color='r')
ax3.spines['left'].set_color('r')
ax3.tick_params(axis='y', right=False, length=9.0, width=1.4, labelsize=15, color='r')
ax3.plot(np.linspace(0.6,10.2,52-3), IntensityR_heat_high_1_5std[3:], label=' heat-high (p<0.05)', color='r', marker='o', linewidth=0.6,  markersize=6.0, linestyle='-')
ax3.legend(loc=(0.486, 0.835), fontsize=12.5)
plt.title("(b)  T2m  &  H500   (yearly intensity, R)", size=18)

ax4 = ax3.twinx()
create_map(ax4)
ax4.set_ylim(0.0,0.048)
ax4.set_yticks([0,0.005,0.01,0.015,0.02,0.025,0.03,0.035,0.04,0.045]); ax4.set_yticklabels([0,"",0.01,"",0.02,"",0.03,"",0.04,""], color='b')
ax4.spines['left'].set_color('r')
ax4.spines['right'].set_color('b')
ax4.tick_params(axis='y', left=False, length=9.0, width=1.4, labelsize=15, color='b')
ax4.plot(np.linspace(0.6,10.2,52-3), IntensityR_heat_low_1_5std[3:], label=' heat-low (p<0.05)', color='b', markerfacecolor='none', markeredgecolor='b', markeredgewidth=1.2, marker='o', linewidth=0.85,  markersize=7.0, linestyle='--')
ax4.plot([1.5,1.5], [0.0,0.018], label='', color='g', linewidth=1.8,  linestyle='--')
ax4.plot([4.0,4.0], [0.0,0.014], label='', color='g', linewidth=1.8,  linestyle='--')
ax4.legend(loc=(0.495, 0.675), fontsize=12.5)

ax3.text(1.85, 0.0245, s='~2800 km', fontsize=12.5, color='b')
ax3.text(4.75, 0.0215, s='~5800 km', fontsize=12.5, color='r')







ax5 = plt.axes([0.0,0.25, 0.4,0.22])
create_map(ax5)
ax5.set_xlabel('Distance  (*1000 km)', size=18)
ax5.set_ylabel('PDF', size=18)
ax5.set_ylim(0.0,0.04)
ax5.set_yticks([0,0.005,0.01,0.015,0.02,0.025,0.03,0.035,0.040]); ax5.set_yticklabels([0,"",0.01,"",0.02,"",0.03,"",0.04], color='r')
ax5.spines['left'].set_color('r')
ax5.tick_params(axis='y', right=False, length=9.0, width=1.4, labelsize=15, color='r')
ax5.plot(np.linspace(0.6,10.2,52-3), concurrent_heat_high_90th[3:], label='heat-high (>40, p<0.01)', color='r', marker='o', linewidth=0.6,  markersize=6.0, linestyle='-')
ax5.legend(loc=(0.408, 0.835), fontsize=12.5)
plt.title("(c)  T2m  &  H500   (concurrent days)", size=18)

ax6 = ax5.twinx()
create_map(ax6)
ax6.set_ylim(0.0,0.18)
ax6.set_yticks([0,0.02,0.04,0.06,0.08,0.10,0.12,0.14,0.16,0.18]); ax6.set_yticklabels([0,"",0.04,"",0.08,"",0.12,"",0.16,""], color='b')
ax6.spines['left'].set_color('r')
ax6.spines['right'].set_color('b')
ax6.tick_params(axis='y', left=False, length=9.0, width=1.4, labelsize=15, color='b')
ax6.plot(np.linspace(0.6,10.2,52-3), concurrent_heat_low_10th[3:], label='heat-low (>40, p<0.01)', color='b', markerfacecolor='none', markeredgecolor='b', markeredgewidth=1.2, marker='o', linewidth=0.85,  markersize=7.0, linestyle='--')
ax6.plot([1.5,1.5], [0.0,0.085], label='', color='g', linewidth=1.8,  linestyle='--')
ax6.plot([4.0,4.0], [0.0,0.045], label='', color='g', linewidth=1.8,  linestyle='--')
ax6.legend(loc=(0.425, 0.675), fontsize=12.5)

ax5.text(1.70, 0.0265, s='~2600 km', fontsize=12.5, color='b')
ax5.text(4.75, 0.022, s='~5700 km', fontsize=12.5, color='r')
ax5.text(3.0, 0.048, s='Pressure   extremes   (H500 > 90th,  < 10th, >= 3days)', fontsize=20.0, color='black')




ax7 = plt.axes([0.52,0.25, 0.4,0.22])
create_map(ax7)
ax7.set_xlabel('Distance  (*1000 km)', size=18)
ax7.set_ylabel('', size=18)
ax7.set_ylim(0.0,0.04)
ax7.set_yticks([0,0.005,0.01,0.015,0.02,0.025,0.03,0.035,0.04]); ax7.set_yticklabels([0,"",0.01,"",0.02,"",0.03,"",0.04], color='r')
ax7.spines['left'].set_color('r')
ax7.tick_params(axis='y', right=False, length=9.0, width=1.4, labelsize=15, color='r')
ax7.plot(np.linspace(0.6,10.2,52-3), IntensityR_heat_high_90th[3:], label=' heat-high (p<0.05)', color='r', marker='o', linewidth=0.6,  markersize=6.0, linestyle='-')
ax7.legend(loc=(0.486, 0.835), fontsize=12.5)
plt.title("(d)  T2m  &  H500   (yearly intensity, R)", size=18)

ax8 = ax7.twinx()
create_map(ax8)
ax8.set_ylim(0.0,0.048)
ax8.set_yticks([0,0.005,0.01,0.015,0.02,0.025,0.03,0.035,0.04,0.045]); ax8.set_yticklabels([0,"",0.01,"",0.02,"",0.03,"",0.04,""], color='b')
ax8.spines['left'].set_color('r')
ax8.spines['right'].set_color('b')
ax8.tick_params(axis='y', left=False, length=9.0, width=1.4, labelsize=15, color='b')
ax8.plot(np.linspace(0.6,10.2,52-3), IntensityR_heat_low_10th[3:], label=' heat-low (p<0.05)', color='b', markerfacecolor='none', markeredgecolor='b', markeredgewidth=1.2, marker='o', linewidth=0.85,  markersize=7.0, linestyle='--')
ax8.plot([1.5,1.5], [0.0,0.0185], label='', color='g', linewidth=1.8,  linestyle='--')
ax8.plot([4.0,4.0], [0.0,0.015], label='', color='g', linewidth=1.8,  linestyle='--')
ax8.legend(loc=(0.495, 0.675), fontsize=12.5)

ax7.text(1.85, 0.0245, s='~2800 km', fontsize=12.5, color='b')
ax7.text(4.75, 0.0215, s='~5800 km', fontsize=12.5, color='r')











fig.savefig('FigureS19_pdf__1_5stdhPa.png', bbox_inches = 'tight')
plt.show()

# %% [markdown]
# 


