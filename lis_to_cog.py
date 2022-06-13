import os
import numpy as np
import rioxarray as rxr
import xarray as xr

from rasterio.crs import CRS
from rasterio.io import MemoryFile
from rio_cogeo.cogeo import cog_translate
from rio_cogeo.profiles import cog_profiles
from rasterio.warp import calculate_default_transform



def _fix_coords(ds):

    lons, lats = ds.lon.values[0, :], ds.lat.values[:, 0]

    ds = ds.drop_vars(["lon", "lat"]) \
            .rename_dims({
                "east_west": "lon",
                "north_south": "lat"
                }) \
            .assign_coords(
                lon=lons,
                lat=lats[::-1]
                )

    return ds


def _flatten_layers(ds, vars=["SoilMoist_tavg", "SoilTemp_tavg"]):

    for var in vars:

        if var in list(ds.keys()):

            var_ds = ds[var]

            layers = {
                f"{var}_{layer_num}": var_ds[layer_num] for layer_num in range(len(var_ds))
            }
            ds = ds.drop(var)
            ds = xr.merge([ds, layers])

    return ds


def _reproject(ds):
    
    wgs84_crs = CRS.from_epsg(4326)
    web_mercator_crs = CRS.from_epsg(3857)

    ds = ds.rio.write_crs(wgs84_crs) \
            .rio.set_spatial_dims(x_dim="lon", y_dim="lat") \
            .rio.reproject(web_mercator_crs)
    
    return ds


def _write_cog(ds, inpath, outdir):

    filename = inpath.split('/')[-1][:-2] + ".tif"

    outpath = os.path.join(outdir, filename) if outdir else filename

    ds.rio.to_raster(outpath, driver="COG")


def _select_vars(ds, vars):

    return ds[[vars]]


def to_cog(inpath, vars=None, outdir=None):

    ds = xr.open_dataset(inpath)

    if vars:
        ds = _select_vars(ds, vars)

    ds = _fix_coords(ds)

    ds = _flatten_layers(ds)

    ds = _reproject(ds)

    _write_cog(ds, inpath, outdir)


def to_cog_v2(inpath, vars=None, outdir=None):

    ds = xr.open_dataset(inpath)

    if vars:
        ds = _select_vars(ds, vars)

    ds = _fix_coords(ds)

    ds = _flatten_layers(ds)

    ds_dtype = np.float32
    nodata_value = ds.attrs["missing_value"]

    x_variable, y_variable = "lon", "lat"

    src_height, src_width = ds.dims[y_variable], ds.dims[x_variable]

    lats, lons = ds[y_variable].values, ds[x_variable].values
    xmin, xmax = np.min(lons), np.max(lons)
    ymin, ymax = np.min(lats), np.max(lats)

    src_crs = CRS.from_epsg(4326)
    dst_crs = CRS.from_epsg(3857)

    dst_transform, dst_width, dst_height = calculate_default_transform(
        src_crs,
        dst_crs,
        src_width,
        src_height,
        left=xmin,
        bottom=ymin,
        right=xmax,
        top=ymax,
    )

    output_profile = dict(
        driver="GTiff",
        dtype=ds_dtype,
        count=len(ds),
        transform=dst_transform,
        crs=dst_crs,
        height=src_height,
        width=src_width,
        nodata=nodata_value,
        tiled=True,
        compress="lzw",
        blockxsize=256,
        blockysize=256
    )

    filename = inpath.split('/')[-1][:-2] + "tif"
    outpath = os.path.join(outdir, filename) if outdir else filename

    with MemoryFile() as memfile:
        with memfile.open(**output_profile) as mem:
            data = ds.drop_dims(["time"]).to_array().values
            mem.write(data)
        cog_translate(
            memfile,
            outpath,
            output_profile,
            config=dict(GDAL_NUM_THREADS="ALL_CPUS", GDAL_TIFF_OVR_BLOCKSIZE="128")
        )
    return {"filename": outpath}

