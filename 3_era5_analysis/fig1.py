# %%
import numpy as np
import xarray as xr


##--##--##--##--  distance pdf of 6 networks  --##--##--##--##
ds1 = xr.open_dataset("data/PDF_heatheat_daily.nc")
concurrent_heat_heat = ds1.pdf.data

ds2 = xr.open_dataset("data/PDF_heatcold_daily.nc")
concurrent_heat_cold = ds2.pdf.data

ds3 = xr.open_dataset("data/PDF_heathigh_daily.nc")
concurrent_heat_high = ds3.pdf.data

ds4 = xr.open_dataset("data/PDF_heatlow_daily.nc")
concurrent_heat_low = ds4.pdf.data

ds5 = xr.open_dataset("data/PDF_heathigh_interannualR.nc")
IntensityR_heat_high = ds5.cor_pdf.data

ds6 = xr.open_dataset("data/PDF_heatlow_interannualR.nc")
IntensityR_heat_low = ds6.cor_pdf.data




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
  ax.tick_params(axis='both', length=6.0, width=1.0, labelsize=9)
  ax.spines['bottom'].set_linewidth(1.0)
  ax.spines['top'].set_linewidth(1.0)
  ax.spines['left'].set_linewidth(1.0)
  ax.spines['right'].set_linewidth(1.0)
  return ax


##--##--##--  create figure  --##--##--##
width_in_inches = 210 / 25.4  ## A4 page
height_in_inches = 297 / 25.4

fig = plt.figure(dpi=400,figsize=(width_in_inches,height_in_inches))

ax1 = plt.axes([0.27,0.535, 0.45,0.18])
create_map(ax1)
ax1.set_xlabel('', size=11)
ax1.set_ylabel('PDF', size=11)
ax1.set_ylim(0.0,0.045)
ax1.set_yticks([0,0.005,0.01,0.015,0.02,0.025,0.03,0.035,0.040,0.045]); ax1.set_yticklabels([0,"",0.01,"",0.02,"",0.03,"",0.04,""], color='r')
ax1.spines['left'].set_color('r')
ax1.tick_params(axis='y', right=False, length=6.0, width=1.0, labelsize=8.5, color='r')
ax1.plot(np.linspace(0.6,10.2,52-3), concurrent_heat_heat[3:], label='heat-heat   (p<0.01)', color='r', marker='o', linewidth=0.7,  markersize=3.5, linestyle='-')
ax1.legend(loc=(0.60, 0.735), fontsize=7.0)
plt.title("Distance   distribution   of   links", size=12)


ax2 = ax1.twinx()
create_map(ax2)
ax2.set_ylim(0.0,0.083)
ax2.set_yticks([0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08]); ax2.set_yticklabels([0,"",0.02,"",0.04,"",0.06,"",0.08], color='b')
ax2.spines['left'].set_color('r')
ax2.spines['right'].set_color('b')
ax2.tick_params(axis='y', left=False, length=6.0, width=1.0, labelsize=8.5, color='b')
ax2.plot(np.linspace(0.6,10.2,52-3), concurrent_heat_cold[3:], label='heat-cold   (p<0.01)', color='b', markerfacecolor='none', markeredgecolor='b', markeredgewidth=0.8, marker='o', linewidth=0.7,  markersize=5.0, linestyle='--')
ax2.plot([1.5,1.5], [0.0,0.041], label='', color='g', linewidth=1.2,  linestyle='--')
ax2.plot([4.0,4.0], [0.0,0.026], label='', color='g', linewidth=1.2,  linestyle='--')
ax2.legend(loc=(0.602, 0.615), fontsize=7.0)

ax1.text(1.55, 0.0405, s='Heat  &  Cold   (concurrent days)', fontsize=10)
ax1.text(2.08, 0.0285, s='~2700 km', fontsize=9, color='b')
ax1.text(4.45, 0.0225, s='~5200 km', fontsize=9, color='r')






ax3 = plt.axes([0.09,0.32, 0.365,0.15])
create_map(ax3)
ax3.set_xlabel('Distance  (*1000 km)', size=11)
ax3.set_ylabel('PDF', size=11)
ax3.set_ylim(0.0,0.04)
ax3.set_yticks([0,0.005,0.01,0.015,0.02,0.025,0.03,0.035,0.040]); ax3.set_yticklabels([0,"",0.01,"",0.02,"",0.03,"",0.04], color='r')
ax3.spines['left'].set_color('r')
ax3.tick_params(axis='y', right=False, length=6.0, width=1.0, labelsize=8.5, color='r')
ax3.plot(np.linspace(0.6,10.2,52-3), concurrent_heat_high[3:], label='heat-high (>80, p<0.01)', color='r', marker='o', linewidth=0.7,  markersize=3.5, linestyle='-')
ax3.legend(loc=(0.428, 0.835), fontsize=7)
plt.title("T2m  &  H500   (concurrent days)", size=10)

ax4 = ax3.twinx()
create_map(ax4)
ax4.set_ylim(0.0,0.149)
ax4.set_yticks([0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10,0.11,0.12,0.13,0.14]); ax4.set_yticklabels([0,"","",0.03,"","",0.06,"","",0.09,"","",0.12,"",""], color='b')
ax4.spines['left'].set_color('r')
ax4.spines['right'].set_color('b')
ax4.tick_params(axis='y', left=False, length=6.0, width=1.0, labelsize=8.5, color='b')
ax4.plot(np.linspace(0.6,10.2,52-3), concurrent_heat_low[3:], label='heat-low (>80, p<0.01)', color='b', markerfacecolor='none', markeredgecolor='b', markeredgewidth=0.8, marker='o', linewidth=0.7,  markersize=5.0, linestyle='--')
ax4.plot([1.5,1.5], [0.0,0.064], label='', color='g', linewidth=1.2,  linestyle='--')
ax4.plot([4.0,4.0], [0.0,0.036], label='', color='g', linewidth=1.2,  linestyle='--')
ax4.legend(loc=(0.445, 0.675), fontsize=7)

ax3.text(1.75, 0.0250, s='~2700 km', fontsize=9, color='b')
ax3.text(4.70, 0.0201, s='~5700 km', fontsize=9, color='r')







ax7 = plt.axes([0.575,0.32, 0.365,0.15])
create_map(ax7)
ax7.set_xlabel('Distance  (*1000 km)', size=11)
ax7.set_ylabel('', size=11)
ax7.set_ylim(0.0,0.04)
ax7.set_yticks([0,0.005,0.01,0.015,0.02,0.025,0.03,0.035,0.04]); ax7.set_yticklabels([0,"",0.01,"",0.02,"",0.03,"",0.04], color='r')
ax7.spines['left'].set_color('r')
ax7.tick_params(axis='y', right=False, length=6.0, width=1.0, labelsize=8.5, color='r')
ax7.plot(np.linspace(0.6,10.2,52-3), IntensityR_heat_high[3:], label=' heat-high (p<0.05)', color='r', marker='o', linewidth=0.7,  markersize=3.5, linestyle='-')
ax7.legend(loc=(0.506, 0.835), fontsize=7)
plt.title("T2m  &  H500   (yearly  correlations)", size=10)

ax8 = ax7.twinx()
create_map(ax8)
ax8.set_ylim(0.0,0.053)
ax8.set_yticks([0,0.005,0.01,0.015,0.02,0.025,0.03,0.035,0.04,0.045,0.05]); ax8.set_yticklabels([0,"",0.01,"",0.02,"",0.03,"",0.04,"",0.05], color='b')
ax8.spines['left'].set_color('r')
ax8.spines['right'].set_color('b')
ax8.tick_params(axis='y', left=False, length=6.0, width=1.0, labelsize=8.5, color='b')
ax8.plot(np.linspace(0.6,10.2,52-3), IntensityR_heat_low[3:], label=' heat-low (p<0.05)', color='b', markerfacecolor='none', markeredgecolor='b', markeredgewidth=0.8, marker='o', linewidth=0.7,  markersize=5.0, linestyle='--')
ax8.plot([1.5,1.5], [0.0,0.022], label='', color='g', linewidth=1.2,  linestyle='--')
ax8.plot([4.0,4.0], [0.0,0.0195], label='', color='g', linewidth=1.2,  linestyle='--')
ax8.legend(loc=(0.515, 0.675), fontsize=7)

ax7.text(1.8, 0.0265, s='~2700 km', fontsize=9, color='b')
ax7.text(4.5, 0.0215, s='~5600 km', fontsize=9, color='r')


ax7.text(-7.7, 0.109, s='a', fontsize=13, color='black', fontweight='bold')
ax7.text(-12.2, 0.044, s='b', fontsize=13, color='black', fontweight='bold')
ax7.text(-0.2, 0.044, s='c', fontsize=13, color='black', fontweight='bold')


fig.savefig('/public/home/fcai/extreme1_AT/NC2_1980_2010/3_era5_analysis/Figure1.png', bbox_inches = 'tight')
fig.savefig('/public/home/fcai/extreme1_AT/NC2_1980_2010/3_era5_analysis/Figure1.pdf')
plt.show()





