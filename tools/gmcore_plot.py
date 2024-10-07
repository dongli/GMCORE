import argparse
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as ticker
import matplotlib.dates as mdates
import cartopy.crs as ccrs
import cartopy.util as cutil
from metpy.interpolate import interpolate_1d
from metpy.units import units

def vinterp(zi, var, zo):
	return interpolate_1d(zo, zi, var).squeeze()

def plot_contour(ax, lon, lat, var, cmap=None, levels=None, left_string=None, right_string=None, with_grid=True):
	ax.set_title(var.long_name)
	var, lon = cutil.add_cyclic_point(var, coord=lon)
	im = ax.contourf(lon, lat, var, transform=ccrs.PlateCarree(), cmap=cmap, levels=levels)
	ax.contour(lon, lat, var, transform=ccrs.PlateCarree(), levels=levels, linewidths=0.2, colors='k')
	if with_grid:
		gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, color='gray', alpha=0.5, linestyle='--')
		gl.top_labels = False
		gl.right_labels = False
		gl.xformatter = ccrs.cartopy.mpl.gridliner.LONGITUDE_FORMATTER
		gl.yformatter = ccrs.cartopy.mpl.gridliner.LATITUDE_FORMATTER
	if left_string is not None: ax.set_title(left_string)
	cax = make_axes_locatable(ax).append_axes('right', size='5%', pad=0.05, axes_class=plt.Axes)
	cbar = plt.colorbar(im, cax=cax, orientation='vertical')
	formatter = ticker.ScalarFormatter(useMathText=True)
	formatter.set_scientific(True)
	formatter.set_powerlimits((0, 0))
	cbar.formatter = formatter
	cbar.update_ticks()

def plot_time_series(ax, time, var, color='black'):
	ax.plot(time, var, '-', color=color)
	ax.set_title(var.long_name)
	# locator = mdates.AutoDateLocator(maxticks=4, interval_multiples=True)
	locator = mdates.DayLocator(interval=60)
	formatter = mdates.AutoDateFormatter(locator)
	ax.xaxis.set_major_locator(locator)
	ax.xaxis.set_major_formatter(formatter)
