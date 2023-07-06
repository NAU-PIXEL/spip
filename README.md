![version](https://img.shields.io/badge/version-1.1.3-blue)
![pythonversion](https://img.shields.io/badge/Python-3.6+-blue)
[![DOI](https://zenodo.org/badge/490546598.svg)](https://zenodo.org/badge/latestdoi/490546598)

# SPiP : Spacecraft Pixel footprint Projection

Projection of the pixel footprint from an instrument onboard an orbital spacecraft on a planetary surface.

## Installation & Update
### Method 1: from PyPI (recommended)
**Installation:** `pip3 install planetary-spip`

**Update:** `pip3 install planetary-spip --upgrade` 


### Method 2: from the GitHub repository (development version)
**Installation:** Clone the repository and install with pip:

~~~bash
git clone https://github.com/NAU-PIXEL/spip.git
cd spip
pip3 install .
~~~

**Update:** Go to the previously cloned repository, pull the last updates, and install them with pip:
~~~bash
cd spip
git pull
pip3 install .
~~~

## Usage
~~~python
# Package importation
import spip
import numpy as np

# Genaration of a longitude/latitude grid to project the footprint
lon_array = np.arange(0.25, 360, 0.5)   # 0° -> 360°, 2 points per degree
lat_array = np.arange(-89.75, 90, 0.5)  # -90° -> 90°, 2 points per degree
lon_grid, lat_grid = np.meshgrid(lon_array, lat_array)

# Pixel parameters
lon_px, lat_px = 58.2, 79.3         # Longitude/Latitude of the pixel center [deg]
emer = 61.7                         # Emergence angle of the pixel on the surface [deg]
lon_subsc, lat_subsc = 28.0, 24.9   # Longitude/Latitude sub-spacecraft [deg]
d_Mars_sc = 27212                   # Distance aeroid-spacraft [km]
Rmars = spip.emm_emirs.Rmars        # 3390 [km]
ifov = spip.emm_emirs.ifov_emirs    # 2.7e-3 [rad]

# Generation of the mask map of the pixel footprint
#-- EMM/EMIRS spacecraft
mask_footprint = spip.emm_emirs.emirs_ifov_px_projection(
    lon_grid,
    lat_grid,
    lon_px,
    lat_px,
    emer,
    lon_subsc,
    lat_subsc,
    d_Mars_sc,
    Rmars,
    ifov
    )

# To easily handle multiple pixels, there is the spip.emm_emirs.emirs_ifov_multi_px_projection() function.

#-- General case, cirular pixel
a, b, θ = spip.circular_pixels.params_ellipse_fov(
    lon_px,
    lat_px,
    emer,
    lon_subsc,
    lat_subsc,
    d_Mars_sc,
    ifov
    )

mask_footprint = spip.circular_pixels.in_ellipse_spherical(
    lon_grid,
    lat_grid,
    Rmars,
    a,
    b,
    lon_px,
    lat_px,
    θ
    ) * 1.
    
# > mask_footprint
# For each point of the longitude/latitude grid:
#   | 1 -> Within the pixel footprint ellipse.
#   | 0 -> Outside.
~~~

## Documentation

Planetary Data Workshop 2023: [abstract](https://github.com/NAU-PIXEL/spip/blob/main/docs/Stcherbinine_PDW2023_7008_spip.pdf) & [poster](https://github.com/NAU-PIXEL/spip/blob/main/docs/poster_PDW_SPiP_stcherbinine.pdf)


## Credits

© Aurélien Stcherbinine (2023) NAU-PIXEL

Department of Astronomy and Planetary Science, Northern Arizona University, Flagstaff, AZ, USA


## License
This package is released under a MIT open source license. See [`LICENSE`](https://github.com/NAU-PIXEL/spip/blob/main/LICENSE) for more details.
