#!/usr/bin/env python

import os
import time
import argparse
import numpy as np

from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D
from astropy.nddata.utils import NoOverlapError
from astropy.wcs import WCS

from reproject import reproject_interp
from reproject.mosaicking import reproject_and_coadd

try:
    from mpi4py import MPI
    mpi_available = True
except ImportError:
    mpi_available = False


def parse_arguments():
    parser = argparse.ArgumentParser(description="Produce FITS cutouts (single or coadded) from tile data with RA/DEC positions.")

    parser.add_argument('--input-table', type=str, required=True,
                        help='Path to the input FITS table with RA, DEC, and tile columns')
    parser.add_argument('--ra-col', type=str, required=True,
                        help='Column name for Right Ascension in degrees')
    parser.add_argument('--dec-col', type=str, required=True,
                        help='Column name for Declination in degrees')
    parser.add_argument('--tile-col', type=str, required=True,
                        help='Column name in the table containing tile filenames (comma-separated if multiple)')
    parser.add_argument('--cutout-size', type=float, default=8,
                        help='Cutout size in arcseconds (default: 8 arcsec)')
    parser.add_argument('--band', type=str, required=True,
                        help='Band name (e.g., u,g, r, i) to include in output filenames and directory')
    parser.add_argument('--output-dir', type=str, required=True,
                        help='Base directory to save cutout FITS files')
    parser.add_argument('--use-mpi', action='store_true',
                        help='Use MPI for parallel processing (requires mpi4py)')

    
    return parser.parse_args()


def ensure_output_dir(path, comm=None):
    if comm is not None:
        rank = comm.Get_rank()
        if rank == 0 and not os.path.exists(path):
            os.makedirs(path, exist_ok=True)
        comm.Barrier()
    else:
        os.makedirs(path, exist_ok=True)


def single_tile_cutout(ra, dec, fits_file, cutout_size_arcsec, output_dir, band, rank=0, index=None):
    """
    Generate a cutout from a single FITS tile centered on (ra, dec).
    """
    if not os.path.isfile(fits_file):
        print(f"[Rank {rank}] File {fits_file} does not exist or is not a file. Skipping.")
        return False

    try:
        hdul = fits.open(fits_file)
        header = hdul[0].header
        data = hdul[0].data
        wcs = WCS(header)

        coord = SkyCoord(ra=ra, dec=dec, unit='deg')
        pixel_scale = np.abs(header.get('CD2_2', 1e-5)) * 3600  # arcsec/pixel fallback
        cutout_pixels = int(np.ceil(cutout_size_arcsec / pixel_scale))
        size = (cutout_pixels, cutout_pixels)

        cutout = Cutout2D(data, coord, size=size, wcs=wcs)
        cutout_header = cutout.wcs.to_header()

        hdu_out = fits.PrimaryHDU(data=cutout.data.astype(np.float64), header=cutout_header)

        out_fname = f"cutout_{band}_RA{ra}_DEC{dec}.fits"
        out_path = os.path.join(output_dir, out_fname)
        hdu_out.writeto(out_path, overwrite=True)

        print(f"[Rank {rank}] Saved single-tile cutout for object {index} at RA={ra:.6f}, DEC={dec:.6f}")
        hdul.close()
        return True

    except NoOverlapError:
        print(f"[Rank {rank}] No overlap for single-tile object {index} at RA={ra:.6f}, DEC={dec:.6f}")
        return False
    except Exception as e:
        print(f"[Rank {rank}] Error processing single-tile object {index} at RA={ra:.6f}, DEC={dec:.6f}: {e}")
        return False


def multi_tile_coadd_cutout(ra, dec, fits_files, cutout_size_arcsec, output_dir, band, rank=0, index=None):
    """
    Generate a coadded cutout from multiple FITS tiles centered on (ra, dec).
    """
    hdus = []
    for fpath in fits_files:
        if os.path.isfile(fpath):
            try:
                hdul = fits.open(fpath)
                hdus.append(hdul[0])
            except Exception as e:
                print(f"[Rank {rank}] Error opening {fpath}: {e}. Skipping this file.")
        else:
            print(f"[Rank {rank}] File {fpath} does not exist or is not a file. Skipping.")

    if not hdus:
        print(f"[Rank {rank}] No valid FITS files found for coadd object {index} at RA={ra}, DEC={dec}")
        return False

    try:
        ref_header = hdus[0].header
        wcs_ref = WCS(ref_header)
        pixel_scale = np.abs(ref_header.get('CD2_2', 1e-5)) * 3600  # arcsec/pixel fallback
        cutout_pixels = int(np.ceil(cutout_size_arcsec / pixel_scale))
        size = (cutout_pixels, cutout_pixels)

        coord = SkyCoord(ra=ra, dec=dec, unit='deg')
        cutout_dummy = Cutout2D(hdus[0].data, coord, size=size, wcs=wcs_ref)
        cutout_header = cutout_dummy.wcs.to_header()

        array_coadd, footprint = reproject_and_coadd(hdus, cutout_header,
                                                     shape_out=size,
                                                     reproject_function=reproject_interp,
                                                     combine_function='mean')

        hdu_out = fits.PrimaryHDU(data=array_coadd.astype(np.float64), header=cutout_header)

        out_fname = f"cutout_{band}_RA{ra}_DEC{dec}.fits"
        out_path = os.path.join(output_dir, out_fname)
        hdu_out.writeto(out_path, overwrite=True)

        print(f"[Rank {rank}] Saved coadded cutout for object {index} at RA={ra:.6f}, DEC={dec:.6f}")

        for hdu in hdus:
            hdu.file.close()
        return True

    except NoOverlapError:
        print(f"[Rank {rank}] No overlap for coadd object {index} at RA={ra:.6f}, DEC={dec:.6f}")
        return False
    except Exception as e:
        print(f"[Rank {rank}] Error processing coadd object {index} at RA={ra:.6f}, DEC={dec:.6f}: {e}")
        return False


def process_rows(table_slice, cutout_size_arcsec, output_dir, tile_col, band, ra_col, dec_col, rank=0):
    no_overlap_rows = []

    for i, row in enumerate(table_slice):
        ra = row[ra_col]
        dec = row[dec_col]
        tile_files_str = row[tile_col]

        if not tile_files_str:
            print(f"[Rank {rank}] Warning: Empty tile list for RA={ra}, DEC={dec}. Skipping.")
            continue

        tile_files = [f.strip() for f in tile_files_str.split(',') if f.strip()]
        if len(tile_files) == 1:
            success = single_tile_cutout(ra, dec, tile_files[0], cutout_size_arcsec, output_dir, band, rank=rank, index=i)
            if not success:
                no_overlap_rows.append(row)
        else:
            success = multi_tile_coadd_cutout(ra, dec, tile_files, cutout_size_arcsec, output_dir, band, rank=rank, index=i)
            if not success:
                no_overlap_rows.append(row)

    return no_overlap_rows


def main():
    args = parse_arguments()

    output_dir_band = os.path.join(args.output_dir, args.band)
    if args.use_mpi and mpi_available:
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()
    else:
        comm = None
        rank = 0
        size = 1

    full_table = Table.read(args.input_table)
    total_rows = len(full_table)

    chunk_size = (total_rows + size - 1) // size
    start_idx = rank * chunk_size
    end_idx = min(start_idx + chunk_size, total_rows)
    sub_table = full_table[start_idx:end_idx]

    if rank == 0:
        print(f"Total rows: {total_rows}, MPI size: {size}")
        print(f"Output directory: {output_dir_band}")

    ensure_output_dir(output_dir_band, comm=comm)

    start_time = time.time()
    print(f"[Rank {rank}] Processing rows {start_idx} to {end_idx - 1}...")

    no_overlap_rows = process_rows(sub_table, args.cutout_size, output_dir_band,
                                  args.tile_col, args.band, args.ra_col, args.dec_col,
                                  rank=rank)

    if comm:
        all_no_overlap = comm.gather(no_overlap_rows, root=0)
        if rank == 0:
            all_no_overlap_flat = [item for sublist in all_no_overlap for item in sublist]
            if all_no_overlap_flat:
                nomatch_table = Table(rows=all_no_overlap_flat, names=full_table.colnames)
                no_overlap_path = os.path.join(output_dir_band, "no_overlap_objects.fits")
                nomatch_table.write(no_overlap_path, overwrite=True)
                print(f"[Rank 0] Wrote {len(all_no_overlap_flat)} no-overlap entries to {no_overlap_path}")
    else:
        if no_overlap_rows:
            nomatch_table = Table(rows=no_overlap_rows, names=full_table.colnames)
            no_overlap_path = os.path.join(output_dir_band, "no_overlap_objects.fits")
            nomatch_table.write(no_overlap_path, overwrite=True)
            print(f"Wrote {len(no_overlap_rows)} no-overlap entries to {no_overlap_path}")

    elapsed_hours = (time.time() - start_time) / 3600.0
    print(f"[Rank {rank}] Elapsed time: {elapsed_hours:.3f} hours")


if __name__ == '__main__':
    main()
