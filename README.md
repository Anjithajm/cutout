# 🌕🌒 FITS Cutout Generator

This repository provides a Python-based script for generating cutouts of objects in FITS format— either from single FITS tiles or via coadded mosaics from multiple overlapping tiles — at specified sky positions.
The code is optimized for both serial and MPI-parallel execution.

---

## 📑 Overview

This reads an input table containing target celestial coordinates (RA, Dec) along with associated FITS tile filenames, and produces cutout images centered at those positions. 
If the object is located in multiple tiles, the code coadds them using interpolation-based reprojection.

**Key features:**
-cutout generation for objects located in single-tile and multiple-tile.

---

## 🛠️ Dependencies

The following Python packages are required:

- `numpy==2.1.3`
- `scipy==1.15.2`
- `astropy==7.0.1`
- `reproject==0.14.1`

Install them via:
```bash
pip install numpy==2.1.3 scipy==1.15.2 astropy==7.0.1 reproject==0.14.1
```


## Usage Example

The following command demonstrates how to run the cutout generation script with a sample input catalog, specifying the relevant columns for Right Ascension, Declination, and tile filenames, as well as output parameters:

```bash
python cutout.py \
  --input-table data.fits \
  --ra-col RAJ2000 \
  --dec-col DECJ2000 \
  --tile-col tile \
  --band g \
  --cutout-size 8 \
  --output-dir cutouts
```
To run the script using MPI with 8 parallel processes, use the following command:
 ```bash
mpirun -n 8 python cutout.py \
  --input-table data.fits \
  --ra-col RAJ2000 \
  --dec-col DECJ2000 \
  --tile-col tile \
  --cutout-size 8 \
  --band r \
  --output-dir ./cutouts \
  --use-mpi
```
- `--input-table test_catalog_input.fits`  
  Specifies the input FITS table containing the target objects with sky coordinates and associated tile filenames.

- `--ra-col RAJ2000`  
  Indicates the column name in the input table that stores the Right Ascension (RA) in degrees, here `RAJ2000`.

- `--dec-col DECJ2000`  
  Specifies the column name containing Declination (DEC) in degrees, here `DECJ2000`.

- `--tile-col tile`  
  Defines the column containing the tile filenames. The column format is String. These filenames may be a comma-separated list if multiple tiles cover the target. For example, see the `data.fits` file

- `--band g`  
  Sets the photometric band label used to tag output files and directories (e.g., 'g' band).

- `--cutout-size 8`  
  Sets the size of the cutout on each side (e.g., 8 arcsec).

- `--output-dir test_cutouts`  
   directory where the generated FITS cutout files will be saved.
