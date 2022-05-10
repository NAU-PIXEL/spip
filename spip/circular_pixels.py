#!/usr/bin/env python3
# -*- coding: utf-8 -*-

## circular_pixels.py
## Created by Aurélien STCHERBININE
## Last modified by Aurélien STCHERBININE : 09/05/2022

##----------------------------------------------------------------------------------------
"""Projection of circular pixels fov.
"""
##----------------------------------------------------------------------------------------
##----------------------------------------------------------------------------------------
## Packages
# Global
import numpy as np

# Name of the current file
_py_file = 'circular_pixels.py'

##-----------------------------------------------------------------------------------
## Distrance between two points in a spherical referential
def spherical_dist_AB(coordsA, coordsB):
    """Compute the distance between two points A and B on a sphere of radius R.

    Convention:
      * R: radius
      * θ: colatitude [rad]
      * φ: longitude [rad]

    Parameters
    ==========
    coordsA : tuple of floats
        Spherical coordinates of point A (R, θ, φ).
    coordsB : tuple of floats
        Spherical coordinates of point B (R, θ', φ').

    Returns
    =======
    dAB : float
        The distance between A and B.
    """
    RA, θA, φA = coordsA
    RB, θB, φB = coordsB
    if RA != RB:
        raise ValueError('A and B must have the same radius coordinate')
    dAB = RA * np.arccos( np.cos(θA)*np.cos(θB) + np.sin(θA)*np.sin(θB)*np.cos(φB - φA) )
    return dAB

def spherical_dist_AB_latlon(latA, lonA, latB, lonB, R):
    """Compute the distance between two points A and B on a sphere of radius R.

    Parameters
    ==========
    latA : float
        Latitude of point A.
    lonA : float
        Longitude of point A.
    latB : float
        Latitude of point B.
    lonB : float
        Longitude of point B.
    R : float
        Radius of the sphere.

    Returns
    =======
    dAB : float
        The distance between A and B.
    """
    deg2rad = np.pi / 180
    θA, φA = (90 - latA) * deg2rad, lonA * deg2rad
    θB, φB = (90 - latB) * deg2rad, lonB * deg2rad
    dAB = R * np.arccos( np.cos(θA)*np.cos(θB) + np.sin(θA)*np.sin(θB)*np.cos(φB - φA) )
    return dAB

##-----------------------------------------------------------------------------------
## Retrieve the parameters of the projected ellipse of the pixel footprint
def params_ellipse_fov(lon_px, lat_px, emergence, lon_subsc, lat_subsc, d_Mars_sc,
                       ifov=2.7e-3):
    """Retrieve the parameters of the projected ellipse of a pixel FOV on Mars.

    Parameters
    ==========
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
    ifov : float, optional (default 2.7e-3)
        The half angular size of the EMIRS pixel FOV [rad].

    Returns
    =======
    a : float
        The FOV ellipse semi-major axis.
    b : float
        The FOV ellipse semi-minor axis.
    θ : float
        The rotation angle between the semi-major axis of the ellipse and the x-axis [rad].
    """
    # IFOV radius in km
    r_ifov = d_Mars_sc * np.tan(ifov)
    # Ellipse axes
    a = r_ifov / np.cos(emergence * np.pi / 180)
    b = r_ifov
    # Ellipse rotation angle
    #> v10 - Spherical computation of distances - Formule des cos
    dist_lon = spherical_dist_AB_latlon(lat_subsc, lon_px, lat_subsc, lon_subsc, 1)
    dist_lat = spherical_dist_AB_latlon(lat_px, lon_px, lat_subsc, lon_px, 1)
    dist_pts = spherical_dist_AB_latlon(lat_px, lon_px, lat_subsc, lon_subsc, 1)
    Δlon = lon_px - lon_subsc
    Δlat = lat_px - lat_subsc
    # > Cosinus formula
    # > cos(a) = cos(b)cos(c) + sin(b)sin(c)cos(α)
    # >> α = Arccos( ( cos(a) - cos(b)cos(c) ) / ( sin(b)sin(c) ) )
    α = np.arccos(
            ( np.cos(dist_lon) - np.cos(dist_lat) * np.cos(dist_pts) )
            / ( np.sin(dist_lat) * np.sin(dist_pts) )
                )
    if (Δlon % 360) <= 180:
        if Δlat >= 0:
            θ = (np.pi / 2) - α
        else:
            θ = (np.pi / 2) + α
    else:
        if Δlat >= 0:
            θ = (np.pi / 2) + α
        else:
            θ = (np.pi / 2) - α
    return a, b, θ

##-----------------------------------------------------------------------------------
## Test if a point of coordinates (x, y) / (lon, lat) is within an ellipse
def in_ellipse(x, y, a, b, x0=0, y0=0, θ=0):
    """Test if a point of coordinates (x, y) is within an ellipse.

    Parameters
    ==========
    x : float or ndarray
        The x-coordinate(s).
    y : float or ndarray
        The y-coordinate(s).
    a : float
        The semi-major axis af the ellipse.
    b : float
        The semi-minor axis of the ellipse.
    x0 : float, optional (default 0)
        The x-coordinate of the center of the ellipse.
    y0 : float, optional (default 0)
        The y-coordinate of the center of the ellipse.
    θ : float, optional (default 0)
        The angle between the semi-major axis of the ellipse and the x-axis [rad].

    Returns
    =======
    testin : bool or ndarray
        True if the point is in the ellipse, False if outside.
    """
    testin = (
        ( (x - x0) * np.cos(θ) + (y - y0) * np.sin(θ) )**2 / a**2 +
        ( (x - x0) * np.sin(θ) - (y - y0) * np.cos(θ) )**2 / b**2
        )
    return (testin <= 1)

def in_ellipse_spherical(lon, lat, R, a, b, lon0=0, lat0=0, θ=0):
    """Test if a point of coordinates (lat, lon) on a sphere of radius R is within an ellipse.

    Parameters
    ==========
    lon : float or ndarray
        Longitude of the point [deg].
    lat : float or ndarray
        Latitude of the point [deg].
    R : float
        Radius of the sphere [km].
    a : float
        The semi-major axis af the ellipse [km].
    b : float
        The semi-minor axis of the ellipse [km].
    lon0 : float, optional (default 0)
        The longitude of the center of the ellipse [deg].
    lat0 : float, optional (default 0)
        The latitude of the center of the ellipse [deg].
    θ : float, optional (default 0)
        The angle between the semi-major axis of the ellipse and the x-axis [rad].

    Returns
    =======
    testin : bool or ndarray
        True if the point is in the ellipse, False if outside.
    """
    deg2rad = np.pi / 180
    # Compute distances to the center of the ellipse in km
    Δx = 2 * (R * np.cos(lat*deg2rad)) * np.tan( (lon - lon0)*deg2rad / 2 )
    Δy = 2 * R * np.tan( (lat - lat0)*deg2rad / 2 )
    # Test if in ellipse
    testin = (
        ( Δx * np.cos(θ) + Δy * np.sin(θ) )**2 / a**2 +
        ( Δx * np.sin(θ) - Δy * np.cos(θ) )**2 / b**2
        )
    return (testin <= 1)

##-----------------------------------------------------------------------------------
## End of code
##-----------------------------------------------------------------------------------
