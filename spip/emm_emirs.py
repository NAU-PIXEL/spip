#!/usr/bin/env python3
# -*- coding: utf-8 -*-

## emm_emirs.py
## Created by Aurélien STCHERBININE
## Last modified by Aurélien STCHERBININE : 11/05/2022

##----------------------------------------------------------------------------------------
"""Projection of an EMM/EMIRS pixel fov on Mars.
"""
##----------------------------------------------------------------------------------------
##----------------------------------------------------------------------------------------
## Packages
# Global
import astropy.units as u
from copy import deepcopy
import os
import yaml
from tqdm import tqdm
# Local
from . import circular_pixels

# Name of the current file
_py_file = 'emm_emirs.py'
# Data path
data_path = os.path.join(
    os.path.abspath(os.path.dirname(__file__)),
    'data')

##-----------------------------------------------------------------------------------
## Load data
with open(os.path.join(data_path, 'instruments.yml')) as f_inst:
    instruments = yaml.safe_load(f_inst)
with open(os.path.join(data_path, 'celestial_objects.yml')) as f_obj:
    objects = yaml.safe_load(f_obj)

# EMIRS IFOV [rad]
ifov_emirs = u.Quantity(instruments['instruments']['emm-emirs']['ifov']).to('rad').value
# Mars radius [km]
Rmars = u.Quantity(objects['planets']['mars']['radius']).to('km').value

##-----------------------------------------------------------------------------------
## EMIRS pixel footprint projection
def emirs_ifov_px_projection(lon_grid, lat_grid, lon_px, lat_px, emer, lon_subsc, lat_subsc, 
                             d_Mars_sc, Rmars=Rmars, ifov=ifov_emirs):
    """Identify the pixels of a lat/lon grid that are within an EMIRS pixel footprint.

    Parameters
    ==========
    lon_grid : 2D ndarray
        The longitude grid [deg].
    lat_grid : 2D ndarray
        The latitude grid [deg].
    lon_px : float
        Center longitude of the EMIRS pixel FOV [deg].
    lat_px : float
        Center latitude of the EMIRS pixel FOV [deg].
    emer : float
        The emergence angle of the pixel as seen by EMIRS [deg].
    lon_subsc : float
        The sub-spacecraft longitude [deg].
    lat_subsc : float
        The sub-spacecraft latitude [deg].
    d_Mars_sc : float
        The distance between the spacraft and Mars [km].
    Rmars : float, optional (default 3390)
        Martian radius [km].
    ifov : float, optional (default 2.7e-3)
        The half angular size of the EMIRS pixel FOV [rad].

    Returns
    =======
    mask_ellipse : 2D ndarray
        Identify the pixels of the lat/lon grid within the EMIRS footprint.
        | 1 -> Within
        | 0 -> Outside
    """
    # Compute footprint ellipse parameters
    a, b, θ = circular_pixels.params_ellipse_fov(lon_px, lat_px, emer, lon_subsc, lat_subsc, d_Mars_sc, ifov)
    # Search for grid pixels within the EMIRS footprint
    mask_ellipse = circular_pixels.in_ellipse_spherical(lon_grid, lat_grid, Rmars, a, b, lon_px, lat_px, θ)
    mask_ellipse2 = deepcopy(mask_ellipse) * 1.  # bool to float
    # Output
    return mask_ellipse2

def emirs_ifov_multi_px_projection(lon_grid, lat_grid, lon, lat, data, emer, 
                                   lon_subsc, lat_subsc, d_Mars_sc, Rmars=Rmars, ifov=ifov_emirs, 
                                   method='avg', zero_values=True, negative_values=False, 
                                   hide_progressbar=False, leave_progressbar=True):
    """Identify the pixels of a lat/lon grid that are within a list of EMIRS pixels footprints.

    Parameters
    ==========
    lon_grid : 2D ndarray
        The longitude grid [deg].
    lat_grid : 2D ndarray
        The latitude grid [deg].
    lon : 1D array of floats
        Center longitude of the EMIRS pixels FOV [deg].
    lat : 1D array of floats
        Center latitude of the EMIRS pixels FOV [deg].
    data : 1D array of floats
        The data value associated with each pixel.
    emer : 1D array of floats
        The emergence angles of the pixels as seen by EMIRS [deg].
    lon_subsc : 1D array of floats
        The sub-spacecraft longitudes [deg].
    lat_subsc : 1D array of floats
        The sub-spacecraft latitudes [deg].
    d_Mars_sc : 1D array of floats
        The distance between the spacraft and Mars for each pixel [km].
    Rmars : float, optional (default 3390)
        Martian radius [km].
    ifov : float, optional (default 2.7e-3)
        The half angular size of the EMIRS pixel FOV [rad].
    method : str, optional (default 'avg')
        The method for overlaping pixels.
        | 'avg' / 'average' / 'mean' -> Average the data
        | 'max' -> Take the maxmal value.
        | 'min' -> Take the minimal value.
    zero_values : bool, optional (default True)
        Set if the null values are considered as relevant data or not.
    negative_values : bool, optional (default False)
        Set if the negative values are considered as relevant data or not.
    hide_progressbar : bool, optional (default False)
        Disable the diplay of the progress bar.
    leave_progressbar : bool, optional (default True)
        Leave the progress bar displayed on terminal.

    Returns
    =======
    data_grid : 2D array of floats
        The data values, sampled on the provided lon/lat grid.
    mask : 2D array
        The array indicating where the new grid has been filled by the data.
    density : 2D array
        The number of obs data for each pixel.
        (Similar to mask, but without clipping the data to (0,1).)
    lon2 : 1D array of floats
        The longitude values from the longitude grid with data,
        after footprints projection.
        (lon_grid[mask] reshaped to 1D array)
    lat2 : 1D array of floats
        The latitude values from the latitude grid with data,
        after footprints projection.
        (lat_grid[mask] reshaped to 1D array)
    data2 : 1D array of floats
        The data values.
        (data_grid[mask] reshaped to 1D array)
    """
    # Test size compatibility
    if not (len(data) == len(lat) == len(lon)):
        raise ValueError("data, lat, and lon arrays must have the same size")
    # Check method
    if not method in ['avg', 'average', 'mean', 'max', 'min']:
        raise ValueError("Invalid argument for method. "
                "Must be one of 'avg', 'average', 'mean', 'max', 'min'.")
    # Initialisation
    Nlat, Nlon = lon_grid.shape
    data_grid = np.zeros((Nlat, Nlon))
    mask = np.zeros((Nlat, Nlon))
    # Sampling on the new grid
    for i in tqdm(range(len(data)), desc='Pixels footprints projection', disable=hide_progressbar,
                  leave=leave_progressbar):
        longi, lati, datai = lon[i], lat[i], data[i]
        lon_subsc_i, lat_subsc_i, = lon_subsc[i], lat_subsc[i]
        emer_i, d_Mars_sc_i = emer[i], d_Mars_sc[i]
        ignore = (
            np.isnan(longi) or np.isnan(lati)                       # Test coords
            or np.isnan(datai)                                      # Test values
            or (not negative_values and (datai < 0))                # Ignore negative values if specified
            or (not zero_values and (datai == 0))                   # Ignore null values if specified
                )
        if ignore:
            continue
        mask_pixel_fov = emirs_ifov_px_projection(lon_grid, lat_grid, longi, lati, emer_i,
                                                  lon_subsc_i, lat_subsc_i, d_Mars_sc_i, Rmars, ifov)
        mask += mask_pixel_fov
        data_i = mask_pixel_fov * datai
        if method in ['avg', 'average', 'mean']:
            data_grid += data_i
        elif method == 'max':
            mask_pixel_fov[mask_pixel_fov==0] = np.nan
            data_grid = np.nanmax([grid_data, data_i], axis=0)
        elif method == 'min':
            mask_pixel_fov[mask_pixel_fov==0] = np.nan
            data_grid = np.nanmin([grid_data, data_i], axis=0)
    data_grid[mask==0] = np.nan
    if method in ['avg', 'average', 'mean']:
        data_grid = data_grid / mask       # Normalisation
    density = deepcopy(mask)
    mask2 = np.clip(mask, 0, 1).astype(bool)
    # 2D -> 1D arrays
    lon2 = deepcopy(lon_grid)[mask2].reshape(-1)
    lat2 = deepcopy(lat_grid)[mask2].reshape(-1)
    data2 = deepcopy(data_grid)[mask2].reshape(-1)
    return data_grid, mask2, density, lon2, lat2, data2

##-----------------------------------------------------------------------------------
## End of code
##-----------------------------------------------------------------------------------
