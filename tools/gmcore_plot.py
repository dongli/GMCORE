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
import metpy.calc as calc
from scipy.interpolate import CubicSpline
import pint
import os

g = 9.80616 * units.m / units.s ** 2
cp = 1004.5 * units.J / units.kg / units.K

def vinterp_z(zi, var, zo, method='cubic'):
	res = None
	if method == 'cubic':
		if len(np.shape(var)) == 1:
			cubic = CubicSpline(zi[:], var[:])
			res = cubic(zo)
		elif len(np.shape(var)) == 2:
			res = np.zeros((len(zo), np.shape(var)[1]))
			for i in range(np.shape(var)[1]):
				cubic = CubicSpline(zi[::-1,i], var[::-1,i], bc_type='not-a-knot')
				res[:,i] = cubic(zo)
	elif method == 'linear':
		res = interpolate_1d(zo, zi, var)
	else:
		print(f'[Error]: Method {method} is not supported!')
		exit(1)
	if isinstance(var, xr.DataArray):
		res = xr.DataArray(res, coords=var.coords, dims=var.dims)
		res.attrs = var.attrs
		if 'lev' in var.dims:
			res['lev'] = zo
			res['lev'].attrs['long_name'] = 'Height'
			res['lev'].attrs['units'] = zo.units
		elif 'ilev' in var.dims:
			res['ilev'] = zo
			res['ilev'].attrs['long_name'] = 'Height'
			res['ilev'].attrs['units'] = zo.units
	return res

def vinterp_p(zi, var, zo):
	plev = np.array([zo]) * units.hPa
	res = interpolate_1d(plev, zi, var)
	res = xr.DataArray(res, coords=var.coords, dims=var.dims)
	return res

def plot_contour_zonal_height(ax, var,
	cmap=None,
	norm=None,
	levels=None,
	ticks=None,
	left_string=None,
	right_string=None,
	with_grid=True,
	font_size=8,
	use_scientific=False,
	with_contour=False,
	linewidth=0.5,
	cbar_orient='vertical'):
	if left_string is not None:
		ax.set_title(left_string)
		# How to set title font size?
		ax.title.set_fontsize(font_size)
	elif 'long_name' in var.attrs:
		ax.set_title(var.long_name)
	else:
		ax.set_title('')
	if 'lat' in var.dims:
		lat = var.lat.copy()
	elif 'ilat' in var.dims:
		lat = var.ilat.copy()
	else:
		print(f'[Error]: Variable {var.name} has no lat or ilat dimension!')
		exit(1)
	if 'lev' in var.dims:
		lev = var.lev.copy()
	else:
		print(f'[Error]: Variable {var.name} has no lev dimension!')
		exit(1)
	if levels is not None and norm is None:
		im = ax.contourf(lat, lev, var, cmap=cmap, levels=levels, extend='both')
		if with_contour: ax.contour(lat, lev, var, levels=levels, linewidths=linewidth, colors='k')
	elif levels is None and norm is not None:
		im = ax.contourf(lat, lev, var, cmap=cmap, norm=norm, levels=norm.boundaries, extend='both')
		if with_contour: ax.contour(lat, lev, var, norm=norm, levels=norm.boundaries, linewidths=linewidth, colors='k')
	else:
		im = ax.contourf(lat, lev, var, cmap=cmap, extend='both')
		if with_contour: ax.contour(lat, lev, var, levels=levels, linewidths=linewidth, colors='k')
	if with_grid:
		ax.grid(True)
	if cbar_orient == 'vertical':
		cax = make_axes_locatable(ax).append_axes('right', size='2%', pad=0.05, axes_class=plt.Axes)
	else:
		cax = make_axes_locatable(ax).append_axes('bottom', size='5%', pad=0.3, axes_class=plt.Axes)
	if levels is not None and norm is None:
		if ticks is not None:
			cbar = plt.colorbar(im, cax=cax, orientation=cbar_orient, boundaries=levels, ticks=ticks)
		else:
			cbar = plt.colorbar(im, cax=cax, orientation=cbar_orient, boundaries=levels, ticks=levels)
	elif levels is None and norm is not None:
		if ticks is not None:
			cbar = plt.colorbar(im, cax=cax, orientation=cbar_orient, boundaries=norm.boundaries[1::2], ticks=ticks)
		else:
			cbar = plt.colorbar(im, cax=cax, orientation=cbar_orient, boundaries=norm.boundaries[1::2], ticks=norm.boundaries[1::2])
	else:
		cbar = plt.colorbar(im, cax=cax, orientation=cbar_orient)
	formatter = ticker.ScalarFormatter(useMathText=True)
	if use_scientific:
		formatter.set_scientific(True)
		formatter.set_powerlimits((0, 0))
	cbar.formatter = formatter
	cbar.ax.tick_params(labelsize=font_size)
	cbar.update_ticks()

def plot_contour_map(ax, var,
	cmap=None,
	levels=None,
	left_string=None,
	right_string=None,
	with_grid=True,
	font_size=8,
	add_cyclic_point=True,
	use_scientific=False,
	show_lat_labels=True,
	with_contour=True,
	linewidth=0.5):
	if left_string is not None:
		ax.set_title(left_string)
		# How to set title font size?
		ax.title.set_fontsize(font_size)
	elif 'long_name' in var.attrs:
		ax.set_title(var.long_name)
	else:
		ax.set_title('')
	if 'lon' in var.dims:
		lon = var.lon.copy()
	elif 'ilon' in var.dims:
		lon = var.ilon.copy()
	else:
		print(f'[Error]: Variable {var.name} has no lon or ilon dimension!')
		exit(1)
	if 'lat' in var.dims:
		lat = var.lat.copy()
	elif 'ilat' in var.dims:
		lat = var.ilat.copy()
	else:
		print(f'[Error]: Variable {var.name} has no lat or ilat dimension!')
		exit(1)
	if add_cyclic_point: var, lon = cutil.add_cyclic_point(var, coord=lon)
	# Use transform_first as in https://scitools.org.uk/cartopy/docs/latest/gallery/scalar_data/contour_transforms.html to avoid failure in contouring.
	lon2d, lat2d = np.meshgrid(lon, lat) # Needs 2D coordinates for transform_first=True.
	im = ax.contourf(lon2d, lat2d, var, transform=ccrs.PlateCarree(), cmap=cmap, levels=levels, extend='both', transform_first=True)
	if with_contour: ax.contour(lon2d, lat2d, var, transform=ccrs.PlateCarree(), levels=levels, linewidths=linewidth, colors='k', transform_first=True)
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

def plot_time_series(ax, time, var, ylim=None, color='black', font_size=8):
	ax.plot(time, var, '-', color=color)
	ax.tick_params(axis='both', labelsize=font_size)
	if ylim is not None: ax.set_ylim(ylim)
	ax.set_title(var.long_name)
	if time.dtype == '<M8[ns]':
		# locator = mdates.AutoDateLocator(maxticks=4, interval_multiples=True)
		locator = mdates.DayLocator(interval=60)
		formatter = mdates.AutoDateFormatter(locator)
		ax.xaxis.set_major_locator(locator)
		ax.xaxis.set_major_formatter(formatter)
		ax.set_xlabel('Time')
	else:
		ax.set_xlabel('Sol')
