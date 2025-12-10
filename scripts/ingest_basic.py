"""
Basic ingestion script for the Galaxy Local Group Database.
Uses a single reference for now.

Columns:
| Column Name   | Destination table | Description                                                                                                                                                                       |
|---------------|-------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| GalaxyName    |   Sources           | Galaxy name                                                                                                                                                                       |
| RA            |   Sources           | Right Ascension (hours minutes seconds)                                                                                                                                           |
| Dec           |   Sources           | Declination (degrees arcminutes arcseconds)                                                                                                                                       |
| EB-V          |  skip            | Foreground reddening, E(B-V), measured directly from the Schlegel et al 1998 maps (they do not include the recalibration by Schlafly & Finbiener 2011)                            |
| dmod          |    Modeled Parameters          | Distance Modulus, (m-M)o err+ err-                                                                                                                                                |
| dmod+         |              |                                                                                                                                                                                   |
| dmod-         |              |                                                                                                                                                                                   |
| vh            |    Radial Velocity          | Heliocentric radial velocity, vh(km/s) err+ err-                                                                                                                                  |
| vh+           |              |                                                                                                                                                                                   |
| vh-           |              |                                                                                                                                                                                   |
| Vmag          |   Photometry           | Apparent V magnitude in Vega mags, Vmag err+ err-                                                                                                                                 |
| Vmag+         |              |                                                                                                                                                                                   |
| Vmag-         |              |                                                                                                                                                                                   |
| Classification            |   Morphology           |  galaxy type (Irr, Sph, etc)                                                                                                   |
| PA            |   Morphology           | Position Angle of major axis in degrees measured east from north, PA err+ err-                                                                                                    |
| PA+           |              |                                                                                                                                                                                   |
| PA-           |              |                                                                                                                                                                                   |
| e=1-b/a       |    Morphology          | Projected ellipticity, e=1-b/a err+ err-                                                                                                                                          |
| e+            |              |                                                                                                                                                                                   |
| e-            |              |                                                                                                                                                                                   |
| muVo          |  Modeled Parameters          | Central V surface brightness, muVo(mag/sq.arcsec) err+ err-                                                                                                                       |
| muVo+         |              |                                                                                                                                                                                   |
| muVo-         |              |                                                                                                                                                                                   |
| rh            |     Morphology         | Half-light radius measured on major axis, rh(arcmins) err+ err-                                                                                                                   |
| rh+           |              |                                                                                                                                                                                   |
| rh-           |              |                                                                                                                                                                                   |
| sigma_s       |  Modeled Parameters            | Stellar radial velocity dispersion, sigma_s(km/s) err+ err-                                                                                                                       |
| sigma_s+      |              |                                                                                                                                                                                   |
| sigma_s-      |              |                                                                                                                    |
| vrot_s        |   Modeled Parameters           |       Stellar peak/max rotation velocity, vrot_s(km/s) err+ err-                                                                                                                                                                                 |
| vrot_s+       |              |                                                                                                                                                                                   |
| vrot_s-       |              |                                                                                                                                                                                   |
| MHI           |  Modeled Parameters            | Mass of HI (calculated for the adopted distance modulus), M_HI (10^6M_sun)                                                                                                        |
| sigma_g       |  Modeled Parameters            | HI radial velocity dispersion, sigma_g(km/s) err+ err-                                                                                                                            |
| sigma_g+      |              |                                                                                                                                                                                   |
| sigma_g-      |              |                                                                                                                                                                                   |
| vrot_g        |   Modeled Parameters           | HI peak/max rotation velocity, vrot_g(km/s) err+ err-                                                                                                                             |
| vrot_g+       |              |                                                                                                                                                                                   |
| vrot_g-       |              |                                                                                                                                                                                   |
| [Fe/H]        |   Modeled Parameters            | Stellar mean metallicity in dex, [Fe/H] err+ err-                                                                                                                                 |
| feh+          |              |                                                                                                                                                                                   |
| feh-          |              |                                                                                                                                                                                   |
| F             |              | Flag, F, for technique used to measure metallicity: 1=High res. spectroscopy; 2=Spectral synthesis; 3=CaT calibration; 4=CMD modelling; 5=Isochrones; 6=color of red giant branch |
| pmra          |   Proper Motions           | mu_alpha*cos(delta) - proper motion in RA direction in mas/yr,  +err, -err                                                                                                        |
| epmra+        |              |                                                                                                                                                                                   |
| epmra-        |              |                                                                                                                                                                                   |
| pmdec         |    Proper Motions          | mu_delta - proper motion in declination direction in mas/yr, +err, -err                                                                                                           |
| epmdec+       |              |                                                                                                                                                                                   |
| epmdec-       |              |                                                                                                                                                                                   |
| References    |   Publications           | References, numbers refer to list given in supplementary file "References.dat"                                                                                                    |
"""

import sqlalchemy as sa
import numpy as np
from astrodb_utils import build_db_from_json
from astrodb_utils.loaders import DatabaseSettings
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table

def sanitize_value(value):
    # Handle custom use of 999 for missing values
    if value is None:
        return None
    
    if isinstance(value, (float, np.float32, np.float64)) and bool(np.isclose(value, 999)):
        return None

    # Handle numpy MaskedConstant (masked value from astropy.table.Table reads)
    if type(value).__module__ == "numpy.ma.core" and type(value).__name__ == "MaskedConstant":
        return None
    
    return value


# Load FITS file
t = Table.read("NearbyGalaxies_Jan2021_PUBLIC.fits")

# Instantiate database from settings file
db_settings = DatabaseSettings(settings_file="database.toml")
db = build_db_from_json(settings_file=db_settings.settings_file)

# Reference to use
reference = "McConnachie12"

# Loop over rows of FITS file
for row in t:
    # Convert RA and Dec from hours, minutes, seconds to degrees
    coord = SkyCoord(ra=row["RA"], dec=row["Dec"], unit=(u.hourangle, u.deg))
    ra_deg = coord.ra.deg
    dec_deg = coord.dec.deg
    source = row["GalaxyName"].strip().replace("*", "")

    # Prepare source data
    source_data = {
        "source": source,
        "ra_deg": ra_deg,
        "dec_deg": dec_deg,
        "equinox": "ICRS",
        "reference": reference,
    }
    print(source_data)

    # Insert source data into database
    try:
        with db.engine.connect() as conn:
            conn.execute(db.Sources.insert().values(source_data))
            conn.commit()
    except Exception:
        # Try an update if an insert fails
        with db.engine.connect() as conn:
            conn.execute(db.Sources.update().where(db.Sources.c.source == source_data["source"]).values(source_data))
            conn.commit()

    # Morphology data
    morphology_data = {
        "source": source,
        "position_angle_deg": sanitize_value(row["PA"]),
        "position_angle_error_upper": sanitize_value(row["PA+"]),
        "position_angle_error_lower": sanitize_value(row["PA-"]),
        "ellipticity": sanitize_value(row["e=1-b/a"]),
        "ellipticity_error_upper": sanitize_value(row["e+"]),
        "ellipticity_error_lower": sanitize_value(row["e-"]),
        "half_light_radius_arcmin": sanitize_value(row["rh"]),
        "half_light_radius_error_upper": sanitize_value(row["rh+"]),
        "half_light_radius_error_lower": sanitize_value(row["rh-"]),
        "adopted": True,
        "reference": reference,
    }
    print(morphology_data)

    # Insert morphology data into database
    try:
        with db.engine.connect() as conn:
            conn.execute(db.Morphology.insert().values(morphology_data))
            conn.commit()
    except sa.exc.IntegrityError:
        # Try an update if an insert fails (there is existing data for this source and reference)
        with db.engine.connect() as conn:
            conn.execute(db.Morphology.update().where(db.Morphology.c.source == source).values(morphology_data))
            conn.commit()

    # Radial velocity data
    radial_velocity_data = {
        "source": source,
        "rv_kms": sanitize_value(row["vh"]),
        "rv_error_upper": sanitize_value(row["vh+"]),
        "rv_error_lower": sanitize_value(row["vh-"]),
        "adopted": True,
        "reference": reference,
    }

    if radial_velocity_data["rv_kms"] is not None:
        print(radial_velocity_data)
        # Insert radial velocity data into database
        with db.engine.connect() as conn:
            conn.execute(db.RadialVelocities.insert().values(radial_velocity_data))
            conn.commit()

    # Photometry data
    photometry_data = {
        "source": source,
        "band": "Generic/Johnson.V",
        "magnitude": sanitize_value(row["Vmag"]),
        "magnitude_error_upper": sanitize_value(row["Vmag+"]),
        "magnitude_error_lower": sanitize_value(row["Vmag-"]),
        "reference": reference,
    }

    if photometry_data["magnitude"] is not None:
        print(photometry_data)
        # Insert photometry data into database
        with db.engine.connect() as conn:
            conn.execute(db.Photometry.insert().values(photometry_data))
            conn.commit()

    # Proper motions data
    proper_motions_data = {
        "source": source,
        "pm_ra": sanitize_value(row["pmra"]),
        "pm_ra_error_upper": sanitize_value(row["epmra+"]),
        "pm_ra_error_lower": sanitize_value(row["epmra-"]),
        "pm_dec": sanitize_value(row["pmdec"]),
        "pm_dec_error_upper": sanitize_value(row["epmdec+"]),
        "pm_dec_error_lower": sanitize_value(row["epmdec-"]),
        "adopted": True,
        "reference": reference,
    }

    if proper_motions_data["pm_ra"] is not None and proper_motions_data["pm_dec"] is not None:
        print(proper_motions_data)
        # Insert proper motions data into database
        with db.engine.connect() as conn:
            conn.execute(db.ProperMotions.insert().values(proper_motions_data))
            conn.commit()

    # Modeled parameters data
    parameter_to_name = {
        "dmod": "distance modulus",
        "muVo": "muVo",
        "sigma_s": "sigma_s",
        "vrot_s": "vrot_s",
        "MHI": "MHI",
        "sigma_g": "sigma_g",
        "vrot_g": "vrot_g",
        "[Fe/H]": "metallicity",
    }

    # Loop over the various columns and load them if they have values
    for column in ["dmod", "muVo", "sigma_s", "vrot_s", "MHI", "sigma_g", "vrot_g", "[Fe/H]"]:
        if row.get(column) is not None:
            # Custom handling of error columns for [Fe/H]
            if column == "[Fe/H]":
                error_column = "feh"
            else:
                error_column = column

            modeled_parameters_data = {
                "source": source,
                "model": None,
                "parameter": parameter_to_name.get(column, column),
                "value": sanitize_value(row.get(column)),
                "error_upper": sanitize_value(row.get(error_column + "+")),
                "error_lower": sanitize_value(row.get(error_column + "-")),
                "reference": reference,
            }

            if modeled_parameters_data["value"] is not None:
                print(modeled_parameters_data)
                # Insert modeled parameters data into database
                with db.engine.connect() as conn:
                    conn.execute(db.ModeledParameters.insert().values(modeled_parameters_data))
                    conn.commit()
