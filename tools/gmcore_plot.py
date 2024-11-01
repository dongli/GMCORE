import argparse
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as ticker
import matplotlib.dates as mdates
import matplotlib.path as mpath
import cartopy.crs as ccrs
import cartopy.util as cutil
from metpy.interpolate import interpolate_1d
from metpy.units import units

def vinterp(zi, var, zo):
	plev = np.array([zo]) * units.hPa
	res = interpolate_1d(plev, zi, var)
	res = xr.DataArray(res, coords=var.coords, dims=var.dims)
	return res

def plot_contour(ax, lon, lat, var, cmap=None, levels=None, left_string=None, right_string=None, with_grid=True, font_size=12, add_cyclic_point=True, use_scientific=False, show_lat_labels=True, with_contour=True, contour_linewidth=0.5):
	if left_string is not None:
		ax.set_title(left_string)
		# How to set title font size?
		ax.title.set_fontsize(font_size)
	else:
		ax.set_title('')
	if add_cyclic_point: var, lon = cutil.add_cyclic_point(var, coord=lon)
	im = ax.contourf(lon, lat, var, transform=ccrs.PlateCarree(), cmap=cmap, levels=levels, extend='both')
	if with_contour: ax.contour(lon, lat, var, transform=ccrs.PlateCarree(), levels=levels, linewidths=contour_linewidth, colors='k')
	if with_grid:
		gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, color='gray', alpha=0.5, linestyle='--')
		gl.top_labels = False
		gl.right_labels = False
		gl.xformatter = ccrs.cartopy.mpl.gridliner.LONGITUDE_FORMATTER
		gl.yformatter = ccrs.cartopy.mpl.gridliner.LATITUDE_FORMATTER
		gl.xlabel_style = {'size': font_size}
		gl.ylabel_style = {'size': font_size if show_lat_labels else 0}
	cax = make_axes_locatable(ax).append_axes('right', size='2%', pad=0.05, axes_class=plt.Axes)
	cbar = plt.colorbar(im, cax=cax, orientation='vertical')
	formatter = ticker.ScalarFormatter(useMathText=True)
	if use_scientific:
		formatter.set_scientific(True)
		formatter.set_powerlimits((0, 0))
	cbar.formatter = formatter
	cbar.ax.tick_params(labelsize=font_size)
	cbar.update_ticks()

	if isinstance(ax.projection, ccrs.NorthPolarStereo):
		theta = np.linspace(0, 2 * np.pi, 100)
		center, radius = [0.5, 0.5], 0.5
		verts = np.vstack([np.sin(theta), np.cos(theta)]).T
		circle = mpath.Path(verts * radius + center)
		ax.set_boundary(circle, transform=ax.transAxes)

def plot_time_series(ax, time, var, ylim=None, color='black'):
	ax.plot(time, var, '-', color=color)
	if ylim is not None: ax.set_ylim(ylim)
	ax.set_title(var.long_name)
	# locator = mdates.AutoDateLocator(maxticks=4, interval_multiples=True)
	locator = mdates.DayLocator(interval=60)
	formatter = mdates.AutoDateFormatter(locator)
	ax.xaxis.set_major_locator(locator)
	ax.xaxis.set_major_formatter(formatter)
