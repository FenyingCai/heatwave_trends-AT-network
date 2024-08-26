# %%
import numpy as np
import xarray as xr


##--##--##--##--  read pdf  --##--##--##--##
ds3 = xr.open_dataset("data/Concurrence_heathigh_gt100.nc")
concurrent_heat_high1 = ds3.pdf.data

ds4 = xr.open_dataset("data/Concurrence_heatlow_gt100.nc")
concurrent_heat_low1 = ds4.pdf.data

ds5 = xr.open_dataset("data/Concurrence_heathigh_between100_80.nc")
concurrent_heat_high2 = ds5.pdf.data

ds6 = xr.open_dataset("data/Concurrence_heatlow_between100_80.nc")
concurrent_heat_low2 = ds6.pdf.data

ds7 = xr.open_dataset("data/Concurrence_heathigh_between80_60.nc")
concurrent_heat_high3 = ds7.pdf.data

ds8 = xr.open_dataset("data/Concurrence_heatlow_between80_60.nc")
concurrent_heat_low3 = ds8.pdf.data






ds3 = xr.open_dataset("data/Correlation_heathigh_gt99.nc")
correlate_heat_high1 = ds3.cor_pdf.data

ds4 = xr.open_dataset("data/Correlation_heatlow_gt99.nc")
correlate_heat_low1 = ds4.cor_pdf.data

ds5 = xr.open_dataset("data/Correlation_heathigh_between99_975.nc")
correlate_heat_high2 = ds5.cor_pdf.data

ds6 = xr.open_dataset("data/Correlation_heatlow_between99_975.nc")
correlate_heat_low2 = ds6.cor_pdf.data

ds7 = xr.open_dataset("data/Correlation_heathigh_between975_95.nc")
correlate_heat_high3 = ds7.cor_pdf.data

ds8 = xr.open_dataset("data/Correlation_heatlow_between975_95.nc")
correlate_heat_low3 = ds8.cor_pdf.data







# %%
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




ax1 = plt.axes([0.0,0.65, 0.4,0.22])
create_map(ax1)
ax1.set_xlabel(' ', size=18)
ax1.set_ylabel('PDF', size=18)
ax1.set_ylim(0.0,0.065)
ax1.set_yticks([0,0.01,0.02,0.03,0.04,0.05,0.06]); ax1.set_yticklabels([0,"",0.02,"",0.04,"",0.06], color='r')
ax1.spines['left'].set_color('r')
ax1.tick_params(axis='y', right=False, length=9.0, width=1.4, labelsize=15, color='r')
ax1.plot(np.linspace(0.6,10.2,52-3), concurrent_heat_high1[3:], label='heat-high (>100, p<0.01)', color='r', marker='o', linewidth=0.6,  markersize=6.0, linestyle='-')
ax1.legend(loc=(0.385, 0.835), fontsize=12.5)
plt.title("(a)  High  significant  concurrences", size=18)

ax2 = ax1.twinx()
create_map(ax2)
ax2.set_ylim(0.0,0.24)
ax2.set_yticks([0,0.02,0.04,0.06,0.08,0.10,0.12,0.14,0.16,0.18,0.20,0.22,0.24]); ax2.set_yticklabels([0,"",0.04,"",0.08,"",0.12,"",0.16,"",0.20,"",0.24], color='b')
ax2.spines['left'].set_color('r')
ax2.spines['right'].set_color('b')
ax2.tick_params(axis='y', left=False, length=9.0, width=1.4, labelsize=15, color='b')
ax2.plot(np.linspace(0.6,10.2,52-3), concurrent_heat_low1[3:], label='heat-low (>100, p<0.01)', color='b', markerfacecolor='none', markeredgecolor='b', markeredgewidth=1.2, marker='o', linewidth=0.85,  markersize=7.0, linestyle='--')
ax2.plot([1.5,1.5], [0.0,0.085], label='', color='g', linewidth=1.8,  linestyle='--')
ax2.plot([4.0,4.0], [0.0,0.03], label='', color='g', linewidth=1.8,  linestyle='--')
ax2.legend(loc=(0.40, 0.675), fontsize=12.5)

ax1.text(1.6, 0.0355, s='~2500 km', fontsize=12.5, color='b')
ax1.text(5.0, 0.0183, s='~6000 km', fontsize=12.5, color='r')




ax3 = plt.axes([0.52,0.65, 0.4,0.22])
create_map(ax3)
ax3.set_xlabel(' ', size=18)
ax3.set_ylabel('', size=18)
ax3.set_ylim(0.0,0.05)
ax3.set_yticks([0,0.005,0.01,0.015,0.02,0.025,0.03,0.035,0.04]); ax3.set_yticklabels([0,"",0.01,"",0.02,"",0.03,"",0.04], color='r')
ax3.spines['left'].set_color('r')
ax3.tick_params(axis='y', right=False, length=9.0, width=1.4, labelsize=15, color='r')
ax3.plot(np.linspace(0.6,10.2,52-3), correlate_heat_high1[3:], label=' heat-high (p<0.01)', color='r', marker='o', linewidth=0.6,  markersize=6.0, linestyle='-')
ax3.legend(loc=(0.486, 0.835), fontsize=12.5)
plt.title("(b)  High  significant  correlations", size=18)

ax4 = ax3.twinx()
create_map(ax4)
ax4.set_ylim(0.0,0.055)
ax4.set_yticks([0,0.005,0.01,0.015,0.02,0.025,0.03,0.035,0.04,0.045,0.05,0.055]); ax4.set_yticklabels([0,"",0.01,"",0.02,"",0.03,"",0.04,"",0.05,""], color='b')
ax4.spines['left'].set_color('r')
ax4.spines['right'].set_color('b')
ax4.tick_params(axis='y', left=False, length=9.0, width=1.4, labelsize=15, color='b')
ax4.plot(np.linspace(0.6,10.2,52-3), correlate_heat_low1[3:], label=' heat-low (p<0.01)', color='b', markerfacecolor='none', markeredgecolor='b', markeredgewidth=1.2, marker='o', linewidth=0.85,  markersize=7.0, linestyle='--')
ax4.plot([1.5,1.5], [0.0,0.0205], label='', color='g', linewidth=1.8,  linestyle='--')
ax4.plot([4.0,4.0], [0.0,0.019], label='', color='g', linewidth=1.8,  linestyle='--')
ax4.legend(loc=(0.495, 0.675), fontsize=12.5)

ax3.text(1.8, 0.0363, s='~2700 km', fontsize=12.5, color='b')
ax3.text(4.4, 0.0223, s='~5200 km', fontsize=12.5, color='r')









ax5 = plt.axes([0.0,0.35, 0.4,0.22])
create_map(ax5)
ax5.set_xlabel('', size=18)
ax5.set_ylabel('PDF', size=18)
ax5.set_ylim(0.0,0.045)
ax5.set_yticks([0,0.005,0.01,0.015,0.02,0.025,0.03,0.035,0.040,0.045]); ax5.set_yticklabels([0,"",0.01,"",0.02,"",0.03,"",0.04,""], color='r')
ax5.spines['left'].set_color('r')
ax5.tick_params(axis='y', right=False, length=9.0, width=1.4, labelsize=15, color='r')
ax5.plot(np.linspace(0.6,10.2,52-3), concurrent_heat_high2[3:], label='heat-high (80-100, p<0.01)', color='r', marker='o', linewidth=0.6,  markersize=6.0, linestyle='-')
ax5.legend(loc=(0.358, 0.835), fontsize=12.5)
plt.title("(c)  Medium  significant  concurrences", size=18)

ax6 = ax5.twinx()
create_map(ax6)
ax6.set_ylim(0.0,0.16)
ax6.set_yticks([0,0.02,0.04,0.06,0.08,0.10,0.12,0.14,0.16]); ax6.set_yticklabels([0,"",0.04,"",0.8,"",0.12,"",0.16], color='b')
ax6.spines['left'].set_color('r')
ax6.spines['right'].set_color('b')
ax6.tick_params(axis='y', left=False, length=9.0, width=1.4, labelsize=15, color='b')
ax6.plot(np.linspace(0.6,10.2,52-3), concurrent_heat_low2[3:], label='heat-low (80-100, p<0.01)', color='b', markerfacecolor='none', markeredgecolor='b', markeredgewidth=1.2, marker='o', linewidth=0.85,  markersize=7.0, linestyle='--')
ax6.plot([1.5,1.5], [0.0,0.050], label='', color='g', linewidth=1.8,  linestyle='--')
ax6.plot([4.0,4.0], [0.0,0.037], label='', color='g', linewidth=1.8,  linestyle='--')
ax6.legend(loc=(0.375, 0.675), fontsize=12.5)

ax5.text(1.75, 0.026, s='~2600 km', fontsize=12.5, color='b')
ax5.text(4.60, 0.023, s='~5800 km', fontsize=12.5, color='r')







ax7 = plt.axes([0.52,0.35, 0.4,0.22])
create_map(ax7)
ax7.set_xlabel('', size=18)
ax7.set_ylabel('', size=18)
ax7.set_ylim(0.0,0.040)
ax7.set_yticks([0,0.005,0.01,0.015,0.02,0.025,0.03,0.035,0.04]); ax7.set_yticklabels([0,"",0.01,"",0.02,"",0.03,"",0.4], color='r')
ax7.spines['left'].set_color('r')
ax7.tick_params(axis='y', right=False, length=9.0, width=1.4, labelsize=15, color='r')
ax7.plot(np.linspace(0.6,10.2,52-3), correlate_heat_high2[3:], label=' heat-high (0.01<p<0.025)', color='r', marker='o', linewidth=0.6,  markersize=6.0, linestyle='-')
ax7.legend(loc=(0.366, 0.835), fontsize=12.5)
plt.title("(d)  Medium  significant  correlations", size=18)

ax8 = ax7.twinx()
create_map(ax8)
ax8.set_ylim(0.0,0.055)
ax8.set_yticks([0,0.005,0.01,0.015,0.02,0.025,0.03,0.035,0.04,0.045,0.05]); ax8.set_yticklabels([0,"",0.01,"",0.02,"",0.03,"",0.04,"",0.05], color='b')
ax8.spines['left'].set_color('r')
ax8.spines['right'].set_color('b')
ax8.tick_params(axis='y', left=False, length=9.0, width=1.4, labelsize=15, color='b')
ax8.plot(np.linspace(0.6,10.2,52-3), correlate_heat_low2[3:], label=' heat-low (0.01<p<0.025)', color='b', markerfacecolor='none', markeredgecolor='b', markeredgewidth=1.2, marker='o', linewidth=0.85,  markersize=7.0, linestyle='--')
ax8.plot([1.5,1.5], [0.0,0.0225], label='', color='g', linewidth=1.8,  linestyle='--')
ax8.plot([4.0,4.0], [0.0,0.018], label='', color='g', linewidth=1.8,  linestyle='--')
ax8.legend(loc=(0.380, 0.675), fontsize=12.5)

ax7.text(1.9, 0.0245, s='~2800 km', fontsize=12.5, color='b')
ax7.text(4.8, 0.0225, s='~5600 km', fontsize=12.5, color='r')













ax11 = plt.axes([0.0,0.05, 0.4,0.22])
create_map(ax11)
ax11.set_xlabel('Distance  (*1000 km)', size=18)
ax11.set_ylabel('PDF', size=18)
ax11.set_ylim(0.0,0.045)
ax11.set_yticks([0,0.005,0.01,0.015,0.02,0.025,0.03,0.035,0.040,0.045]); ax11.set_yticklabels([0,"",0.01,"",0.02,"",0.03,"",0.04,""], color='r')
ax11.spines['left'].set_color('r')
ax11.tick_params(axis='y', right=False, length=9.0, width=1.4, labelsize=15, color='r')
ax11.plot(np.linspace(0.6,10.2,52-3), concurrent_heat_high3[3:], label='heat-high (60-80, p=0.01)', color='r', marker='o', linewidth=0.6,  markersize=6.0, linestyle='-')
ax11.legend(loc=(0.383, 0.835), fontsize=12.5)
plt.title("(e)  Low  significant  concurrences", size=18)

ax12 = ax11.twinx()
create_map(ax12)
ax12.set_ylim(0.0,0.10)
ax12.set_yticks([0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10]); ax12.set_yticklabels([0,"",0.02,"",0.4,"",0.06,"",0.08,"",0.10], color='b')
ax12.spines['left'].set_color('r')
ax12.spines['right'].set_color('b')
ax12.tick_params(axis='y', left=False, length=9.0, width=1.4, labelsize=15, color='b')
ax12.plot(np.linspace(0.6,10.2,52-3), concurrent_heat_low3[3:], label='heat-low (60-80, p=0.01)', color='b', markerfacecolor='none', markeredgecolor='b', markeredgewidth=1.2, marker='o', linewidth=0.85,  markersize=7.0, linestyle='--')
ax12.plot([1.5,1.5], [0.0,0.02], label='', color='g', linewidth=1.8,  linestyle='--')
ax12.plot([4.0,4.0], [0.0,0.025], label='', color='g', linewidth=1.8,  linestyle='--')
ax12.legend(loc=(0.400, 0.675), fontsize=12.5)

ax11.text(1.85, 0.027, s='~2800 km', fontsize=12.5, color='b')
ax11.text(4.60, 0.024, s='~5800 km', fontsize=12.5, color='r')







ax13 = plt.axes([0.52,0.05, 0.4,0.22])
create_map(ax13)
ax13.set_xlabel('Distance  (*1000 km)', size=18)
ax13.set_ylabel('', size=18)
ax13.set_ylim(0.0,0.040)
ax13.set_yticks([0,0.005,0.01,0.015,0.02,0.025,0.03,0.035,0.04]); ax13.set_yticklabels([0,"",0.01,"",0.02,"",0.03,"",0.4], color='r')
ax13.spines['left'].set_color('r')
ax13.tick_params(axis='y', right=False, length=9.0, width=1.4, labelsize=15, color='r')
ax13.plot(np.linspace(0.6,10.2,52-3), correlate_heat_high3[3:], label=' heat-high (0.025<p<0.05)', color='r', marker='o', linewidth=0.6,  markersize=6.0, linestyle='-')
ax13.legend(loc=(0.367, 0.835), fontsize=12.5)
plt.title("(f)  Low  significant  correlations", size=18)

ax14 = ax13.twinx()
create_map(ax14)
ax14.set_ylim(0.0,0.055)
ax14.set_yticks([0,0.005,0.01,0.015,0.02,0.025,0.03,0.035,0.04,0.045,0.05]); ax14.set_yticklabels([0,"",0.01,"",0.02,"",0.03,"",0.04,"",0.05], color='b')
ax14.spines['left'].set_color('r')
ax14.spines['right'].set_color('b')
ax14.tick_params(axis='y', left=False, length=9.0, width=1.4, labelsize=15, color='b')
ax14.plot(np.linspace(0.6,10.2,52-3), correlate_heat_low3[3:], label=' heat-low (0.025<p<0.05)', color='b', markerfacecolor='none', markeredgecolor='b', markeredgewidth=1.2, marker='o', linewidth=0.85,  markersize=7.0, linestyle='--')
ax14.plot([1.5,1.5], [0.0,0.020], label='', color='g', linewidth=1.8,  linestyle='--')
ax14.plot([4.0,4.0], [0.0,0.018], label='', color='g', linewidth=1.8,  linestyle='--')
ax14.legend(loc=(0.380, 0.675), fontsize=12.5)

ax13.text(1.9, 0.0230, s='~2800 km', fontsize=12.5, color='b')
ax13.text(4.8, 0.0225, s='~5800 km', fontsize=12.5, color='r')








fig.savefig('FigureS6_pdf_different_threshold.png', bbox_inches = 'tight')
plt.show()

# %% [markdown]
# 


