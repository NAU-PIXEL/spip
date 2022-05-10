#!/usr/bin/env python3
# -*- coding: utf-8 -*-

## emm_emirs.py
## Created by Aurélien STCHERBININE
## Last modified by Aurélien STCHERBININE : 09/05/2022

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

# Name of the current file
_py_file = 'emm_emirs.py'
# Data path
data_path = os.path.join(
    os.path.abspath(os.path.dirname(__file__)),
    'data')

##-----------------------------------------------------------------------------------
## Load data
with open(os.path.join(data_path, 'instruments.yml') as f_inst:
    instruments = yaml.safe_load(f_inst)
with open(os.path.join(data_path, 'celestial_objects.yml') as f_obj:
    objects = yaml.safe_load(f_obj)

# EMIRS IFOV [rad]
ifov_emirs = u.Quantity(instruments['instruments']['emm-emirs']['ifov']).to('rad').value
# Mars radius [km]
Rmars = u.Quantity(instruments['planets']['mars']['radius']).to('km').value

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
    emergence : float
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
    a, b, θ = params_ellipse_fov(lon_px, lat_px, emer, lon_subsc, lat_subsc, d_Mars_sc, ifov)
    # Search for grid pixels within the EMIRS footprint
    mask_ellipse = in_ellipse_spherical(lon_grid, lat_grid, Rmars, a, b, lon_px, lat_px, θ)
    mask_ellipse2 = deepcopy(mask_ellipse) * 1.  # bool to float
    # Output
    return mask_ellipse2

##-----------------------------------------------------------------------------------
## End of code
##-----------------------------------------------------------------------------------
