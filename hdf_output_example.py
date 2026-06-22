"""Example readers for the HDF5 outputs produced by this pipeline.

Two outputs are covered:

* ``tree_bricks{iout:05d}.h5`` — the HaloMaker/GalaxyMaker catalog
  (``read_bricks_hdf`` plus the ``read_members`` / ``read_photometry`` helpers).
* ``gas_bricks{iout:05d}.h5`` — the optional GasMaker post-processing output
  (``read_gas_bricks_hdf``), row-aligned with the catalog and joined by ``id``
  (``join_gas_to_catalog``).

See ``CATALOG_FORMAT.md`` for the authoritative field list, groups, and units.
These functions are intentionally small and dependency-light (only ``h5py`` and
``numpy``) so they can be copied into an analysis script as a starting point.
"""
import h5py
import numpy as np


# Catalog position/radius fields are stored in physical Mpc (see
# CATALOG_FORMAT.md). The box spans [-box_physical, 0] Mpc, so dividing by
# box_physical = aexp * boxsize2 maps positions to [-1, 0); adding 1 shifts them
# to RAMSES code units [0, 1). Radii are lengths, so they are only scaled.
_POS_FIELDS = ("px", "py", "pz", "px*", "py*", "pz*")
_RADIUS_FIELDS = ("r", "r*", "r50", "r90", "rvir")


def read_bricks_hdf(fname, return_params=False, return_pids=False,
                    to_code_unit=True):
    '''
    Read a halo/galaxy catalog from a HaloMaker HDF5 file.

    Parameters
    ----------
    fname : str
        Path to the ``tree_bricks*.h5`` file.
    return_params : bool, optional
        Also return the merged ``/header`` and ``/input`` attributes. Default False.
    return_pids : bool, optional
        Also return the flat member particle IDs and their per-halo offsets.
        Default False. (Unlike older versions, this now works independently of
        ``return_params``.)
    to_code_unit : bool, optional
        Convert all position and radius fields from the stored physical Mpc to
        RAMSES code units [0, 1). Default True (matches the historical example
        default). Pass False to keep the raw physical-Mpc values as stored.

    Returns
    -------
    halo : numpy structured array
        The ``/catalog/halo`` table (all fields, see CATALOG_FORMAT.md).
    params : dict, optional
        Merged ``/header`` + ``/input`` attributes. Returned if ``return_params``.
    pids : numpy array, optional
        Flat member particle IDs (``/member/pids``). Returned if ``return_pids``.
    index : numpy array, optional
        Offsets into ``pids`` of length ``nhalo + 1``; members of catalog row
        ``i`` are ``pids[index[i]:index[i+1]]``. Returned if ``return_pids``.

    Examples
    --------
    >>> halo = read_bricks_hdf("tree_bricks00601.h5")
    >>> halo, params = read_bricks_hdf("tree_bricks00601.h5", return_params=True)
    >>> halo, pids, index = read_bricks_hdf("tree_bricks00601.h5", return_pids=True)
    >>> members_of_row0 = pids[index[0]:index[1]]
    '''
    with h5py.File(fname, 'r') as f:
        halo = f['catalog']['halo'][:]

        if to_code_unit:
            boxsize2 = f['header'].attrs['boxsize2']
            aexp = f['header'].attrs['aexp']
            box_physical = aexp * boxsize2
            names = halo.dtype.names
            for field in _POS_FIELDS:
                if field in names:
                    halo[field] = halo[field] / box_physical + 1.0
            for field in _RADIUS_FIELDS:
                if field in names:
                    halo[field] = halo[field] / box_physical

        out = [halo]
        if return_params:
            params = dict(f['header'].attrs)
            params.update(f['input'].attrs)
            out.append(params)
        if return_pids:
            member = f['member']
            out.append(member['pids'][:])
            out.append(member['index'][:])

    if len(out) == 1:
        return out[0]
    return tuple(out)


def read_members(fname):
    '''
    Read the full ``/member`` group as a dict.

    Returns a dict with whatever member datasets are present. ``index`` (offsets,
    length ``nhalo + 1``) and ``pids`` are always written. ``pos`` / ``vel`` /
    ``mass`` are present only when the run enabled ``dump_DMs``. ``count`` (members
    per halo) is present on recent outputs. Members of catalog row ``i`` are the
    slice ``arr[index[i]:index[i+1]]`` for any flat member array ``arr``.

    >>> mem = read_members("tree_bricks00601.h5")
    >>> i = 0
    >>> pids_i = mem['pids'][mem['index'][i]:mem['index'][i + 1]]
    >>> pos_i = mem['pos'][mem['index'][i]:mem['index'][i + 1]]  # if dump_DMs
    '''
    out = {}
    with h5py.File(fname, 'r') as f:
        member = f['member']
        out['attrs'] = dict(member.attrs)  # e.g. index_base, layout, fields
        for key in ('index', 'count', 'pids', 'pos', 'vel', 'mass'):
            if key in member:
                out[key] = member[key][:]
    return out


def read_photometry(fname, model=None):
    '''
    Read intrinsic rest-frame stellar photometry tables (if present).

    Each model group (``CB07`` / ``BC03`` / ``FSPS``) holds a ``data`` table
    row-aligned with ``/catalog/halo`` (joined by ``id``) plus provenance
    attributes (IMF, source, mass_source, K_filter, ...).

    Parameters
    ----------
    model : str, optional
        Return only this model. If None, return all available models.

    Returns
    -------
    dict
        ``{model_name: {"data": structured_array, "attrs": dict}}``, or for a
        single requested ``model`` the inner ``{"data": ..., "attrs": ...}``.

    >>> phot = read_photometry("tree_bricks00601.h5")
    >>> bc03 = phot["BC03"]["data"]      # umag, gmag, ..., Kmag, age_r, ...
    >>> phot["FSPS"]["attrs"]["K_filter"]
    '''
    out = {}
    with h5py.File(fname, 'r') as f:
        if 'photometry' not in f:
            return out
        group = f['photometry']
        names = [model] if model is not None else list(group.keys())
        for name in names:
            node = group[name]
            out[name] = {
                "data": node['data'][:],
                "attrs": dict(node.attrs),
            }
    if model is not None:
        return out.get(model)
    return out


def read_gas_bricks_hdf(fname):
    '''
    Read a GasMaker ``gas_bricks*.h5`` output.

    Returns
    -------
    summary : numpy structured array
        The ``/gas/summary`` table (gas masses, metallicity, chemistry,
        kinematics, spherical-overdensity quantities, overlap diagnostics).
        Row-aligned with ``/catalog/halo`` of the source catalog; join by ``id``.
        Note GasMaker masses are in Msun (not 10^11 Msun) and ``r200``/``r500``
        are in code units [0, 1) — see CATALOG_FORMAT.md.
    processed : numpy bool array
        Per-row completion mask (only the requested roots' descendants are
        computed; the rest stay NaN / False).
    header : dict
        ``/header`` attributes (overlap settings, schema_version, source_catalog,
        completed_root_count, ...).

    >>> gas, processed, ghdr = read_gas_bricks_hdf("gas_bricks00601.h5")
    >>> gas["mgas_rvir"][processed][:5]
    '''
    with h5py.File(fname, 'r') as f:
        summary = f['gas']['summary'][:]
        processed = f['gas']['processed'][:] if 'processed' in f['gas'] else None
        header = dict(f['header'].attrs)
    return summary, processed, header


def join_gas_to_catalog(halo, gas_summary):
    '''
    Align a GasMaker ``/gas/summary`` table to a catalog ``halo`` table by ``id``.

    Both tables are written row-aligned to the same ``/catalog/halo``, so for a
    matching catalog and gas file this is an identity check. When the arrays may
    differ in length/order (e.g. a filtered catalog), this matches by ``id`` and
    returns the gas rows reordered to follow ``halo`` (missing ids get -1).

    Returns
    -------
    gas_index : numpy int array
        For each row of ``halo``, the index into ``gas_summary`` with the same
        ``id``, or -1 if absent. Use ``gas_summary[gas_index]`` after masking -1.
    '''
    gas_ids = gas_summary['id']
    order = np.argsort(gas_ids)
    pos = np.searchsorted(gas_ids[order], halo['id'])
    pos = np.clip(pos, 0, len(gas_ids) - 1)
    gas_index = np.where(gas_ids[order][pos] == halo['id'], order[pos], -1)
    return gas_index


if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        print("usage: python hdf_output_example.py "
              "<tree_bricks*.h5> [gas_bricks*.h5]")
        raise SystemExit(1)

    catalog_path = sys.argv[1]
    halo, params = read_bricks_hdf(catalog_path, return_params=True)
    print(f"catalog: {catalog_path}")
    print(f"  halos: {halo.shape[0]}  fields: {len(halo.dtype.names)}")
    print(f"  nb_of_halos={params.get('nb_of_halos')} "
          f"nb_of_subhalos={params.get('nb_of_subhalos')} "
          f"zoomin={params.get('zoomin')}")

    phot = read_photometry(catalog_path)
    if phot:
        print(f"  photometry models: {sorted(phot)}")

    if len(sys.argv) >= 3:
        gas_path = sys.argv[2]
        gas, processed, ghdr = read_gas_bricks_hdf(gas_path)
        ncomplete = int(processed.sum()) if processed is not None else gas.shape[0]
        print(f"gas: {gas_path}")
        print(f"  rows: {gas.shape[0]}  completed: {ncomplete}  "
              f"schema_version={ghdr.get('schema_version')}")
        gas_index = join_gas_to_catalog(halo, gas)
        matched = int((gas_index >= 0).sum())
        print(f"  rows matched to catalog by id: {matched}/{halo.shape[0]}")
