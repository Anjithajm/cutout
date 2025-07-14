# 📦 FITS Cutout Generator with MPI Support

This repository provides a Python-based script for generating cutouts of objects in FITS format— either from single FITS tiles or via coadded mosaics from multiple overlapping tiles — at specified sky positions.
The code is optimized for both serial and MPI-parallel execution.

---

## 📑 Overview

This reads an input table containing target celestial coordinates (RA, Dec) along with associated FITS tile filenames, and produces small cutout images centered at those positions. 
If the object located in multiple tiles the code coadds them using interpolation-based reprojection.

**Key features:**
-cutout generation for objects located in single-tile and multi-tiles

---

## 🛠️ Dependencies

The following Python packages are required:

- `numpy==2.1.3`
- `scipy==1.15.2`
- `astropy==7.0.1`
- `reproject==0.14.1`

Install them via:
```bash
pip install numpy==2.1.3 scipy==1.15.2 astropy==7.0.1 reproject==0.14.1```


