import os
import rioxarray as rxr
import xarray as xr

from rasterio.crs import CRS


def _fix_coords(ds):

    lons, lats = ds.lon.values[0, :], ds.lat.values[:, 0]

    ds = ds.drop_vars(["lon", "lat"]) \
            .rename_dims({
                "east_west": "lon",
                "north_south": "lat"
                }) \
            .assign_coords(
                lon=lons,
                lat=lats
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


def lis_to_cog(inpath, vars=None, outdir=None):

    ds = xr.open_dataset(inpath)

    if vars:
        ds = _select_vars(ds, vars)

    ds = _fix_coords(ds)

    ds = _flatten_layers(ds)

    ds = _reproject(ds)

    _write_cog(ds, inpath, outdir)




