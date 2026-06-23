from halo_defs import mem, frange
import halo_defs as H
import numpy as np
import json
import atexit, os, signal, sys
from scipy.io import FortranFile
from tqdm import tqdm
# from multiprocessing import Pool
import multiprocessing as mp
ctx = mp.get_context('fork')
Pool = ctx.Pool

import faulthandler
faulthandler.enable()
faulthandler.register(signal.SIGUSR1, file=sys.stderr, all_threads=True)

UNITS_VERSION = "halomaker_units_v2"

CATALOG_FIELD_UNITS = {
    "id": "int",
    "timestep": "int",
    "nmem": "count",
    "ndm": "count",
    "nstar": "count",
    "nbsub": "count",
    "hosthalo": "int",
    "hostsub": "int",
    "level": "int",
    "nextsub": "int",
    "px": "code_unit",
    "py": "code_unit",
    "pz": "code_unit",
    "px*": "code_unit",
    "py*": "code_unit",
    "pz*": "code_unit",
    "vx": "km/s",
    "vy": "km/s",
    "vz": "km/s",
    "Lx": "Msun Mpc km/s",
    "Ly": "Msun Mpc km/s",
    "Lz": "Msun Mpc km/s",
    "Lx*": "Msun Mpc km/s",
    "Ly*": "Msun Mpc km/s",
    "Lz*": "Msun Mpc km/s",
    "sha": "code_unit",
    "shb": "code_unit",
    "shc": "code_unit",
    "m": "Msun",
    "mdm": "Msun",
    "m*": "Msun",
    "r": "code_unit",
    "r*": "code_unit",
    "r50": "code_unit",
    "r90": "code_unit",
    "age": "Gyr",
    "metal": "mass_fraction",
    "SFR": "Msun/yr",
    "SFR_r50": "Msun/yr",
    "SFR_r90": "Msun/yr",
    "SFR10": "Msun/yr",
    "SFR10_r50": "Msun/yr",
    "SFR10_r90": "Msun/yr",
    "spin": "dimensionless",
    "sigma": "km/s",
    "sigma_dm": "km/s",
    "sigma*": "km/s",
    "vrot": "km/s",
    "sig3d": "km/s",
    "sigcyl": "km/s",
    "vrot_r50": "km/s",
    "sig3d_r50": "km/s",
    "sigcyl_r50": "km/s",
    "vrot_r90": "km/s",
    "sig3d_r90": "km/s",
    "sigcyl_r90": "km/s",
    "ek": "Msun (km/s)^2",
    "ep": "Msun (km/s)^2",
    "et": "Msun (km/s)^2",
    "rvir": "code_unit",
    "mvir": "Msun",
    "tvir": "K",
    "cvel": "km/s",
    "rho_0": "Msun/kpc^3",
    "r_c": "code_unit",
    "cNFW": "dimensionless",
    "cNFWerr": "dimensionless",
    "vmaxcir": "km/s",
    "rmaxcir": "code_unit",
    "inslope": "dimensionless",
    "inslopeerr": "dimensionless",
    "mcontam": "Msun",
}

PHOTOMETRY_FIELD_UNITS = {
    "id": "int",
    "umag": "mag",
    "gmag": "mag",
    "rmag": "mag",
    "imag": "mag",
    "zmag": "mag",
    "Umag": "mag",
    "Bmag": "mag",
    "Vmag": "mag",
    "Kmag": "mag",
    "age_r": "Gyr",
    "metal_r": "mass_fraction",
    "r50_r": "code_unit",
    "r90_r": "code_unit",
}


def _json_attr(mapping):
    return json.dumps(mapping, sort_keys=True)


def _convert_halo_catalog_units(cat, box_physical_mpc):
    out = cat.copy()
    names = set(out.dtype.names)
    for field in ("m", "mdm", "m*", "mvir", "mcontam"):
        if field in names:
            out[field] *= 1.0e11
    for field in ("px", "py", "pz", "px*", "py*", "pz*"):
        if field in names:
            out[field] = np.mod(out[field] / box_physical_mpc, 1.0)
    for field in (
        "r", "r*", "r50", "r90", "rvir", "r_c", "rmaxcir",
        "sha", "shb", "shc",
    ):
        if field in names:
            out[field] /= box_physical_mpc
    for field in ("Lx", "Ly", "Lz", "Lx*", "Ly*", "Lz*", "ek", "ep", "et"):
        if field in names:
            out[field] *= 1.0e11
    if "rho_0" in names:
        out["rho_0"] *= 100.0
    return out


def _convert_photometry_units(photo, box_physical_mpc):
    out = photo.copy()
    for field in ("r50_r", "r90_r"):
        if field in out.dtype.names:
            out[field] /= box_physical_mpc
    return out

#///////////////////////////////////////////////////////////////////////
#***********************************************************************
def read_data_10():
    '''
     This routine read the output of N-body simulations (particles positions and speeds, 
     cosmological and technical parameters)
    @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
     WARNING: this routine just reads the data and converts positions       
              and velocities from CODE units to these units                 
              -- positions are between -0.5 and 0.5                         
              -- velocities are in km/s                                     
                 in units of Hubble velocity accross the box for SIMPLE (SN)
              -- total box mass is 1.0                                      
                 for simulation with hydro (only with -DRENORM) flag        
              -- initial (beg of simulation) expansion factor is ai=1.0      
    @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    '''

    print(f"\n> In read_data: timestep  ---> {H.numero_step}")

    if(H.numero_step == 1):
       print(f"> output_dir: `{H.output_dir}`")
       # contains the number of snapshots to analyze and their names, type and number (see below)
       f12 = open('inputfiles_HaloMaker.dat','r')

    # then read name of snapshot, its type (pm, p3m, SN, Nzo, Gd), num of procs used and number of snapshot
    iscomment = True
    while iscomment:
        name_of_file, H.simtype, H.nbPes, H.numstep = f12.readline().split()
        # skip comment lines starting with # or !
        if name_of_file[0] not in ['#', '!']:
            iscomment = False
    H.nbPes = int(H.nbPes); H.numstep = int(H.numstep)
    if(name_of_file[0]=="'")or(name_of_file[0]=='"'):
        name_of_file = name_of_file[1:-1]
    H.file_num = f"{int(H.numstep):05d}"
    if(name_of_file[0] != '/'):
        name_of_file = f'{H.output_dir}/{name_of_file}'
    print(f"name_of_file: `{name_of_file}`")

    # Note 1: old treecode SNAP format has to be converted [using SNAP_to_SIMPLE (on T3E)] 
    #     into new treecode SIMPLE (SN) format.
    # Note 2: of the five format (pm, p3m, SN, Nzo, Gd) listed above, only  SN, Nzo and Gd 
    #     are fully tested so the code stops for pm and p3m

    if(H.numero_step == H.nsteps): f12.close()

    if(H.simtype=='SN'): raise NotImplementedError("`SN` format is not implemented yet")

    elif(H.simtype=='Ra'):
        read_ramses_100(name_of_file)
        # Computation of omega_t = omega_matter(t)
        #
        #                            omega_f*(1+z)^3
        # omega(z)   = ------------------------------------------------------------------
        #              omega_f*(1+z)^3+(1-omega_f-omega_lambda_f)*(1+z)^2+omega_lambda_f
        #
        #                              omega_lambda_0
        # omega_L(z) = ----------------------------------------------------------------
        #              omega_f*(1+z)^3+(1-omega_f-omega_lambda_f)*(1+z)^2+omega_lambda_f
        H.omega_t  = H.omega_f*(H.af/H.aexp)**3
        H.omega_t  = H.omega_t/(H.omega_t+(1.-H.omega_f-H.omega_lambda_f)*(H.af/H.aexp)**2+H.omega_lambda_f)
    elif(H.simtype[:2]=='Ra'):
        read_ramses_new_101(name_of_file, rver=H.simtype)
        H.omega_t  = H.omega_f*(H.af/H.aexp)**3
        H.omega_t  = H.omega_t/(H.omega_t+(1.-H.omega_f-H.omega_lambda_f)*(H.af/H.aexp)**2+H.omega_lambda_f)
    elif(H.simtype=='Nzo'): raise NotImplementedError("`Nzo` format is not implemented yet")
    elif(H.simtype=='Gd'): raise NotImplementedError("`Gd` format is not implemented yet")
    else: raise NotImplementedError(f"> Don''t know the snapshot format: `{H.simtype}`")

    pos = mem['pos_10']
    print(f"> min max position (in box units)   : {np.min(pos)},{np.max(pos)}")
    vel = mem['vel_10']
    print(f"> min max velocities (in km/s)      : {np.min(vel)},{np.max(vel)}")
    print(f"> Reading done.")
    print(f"> aexp = {H.aexp}")

def skip_records(f, skip_num=1):
    """
    Skips a record from current position, faster than read_ints.
    """
    for _ in range(skip_num):
        first_size = f._read_size()

        f._fp.seek(first_size, 1)

        second_size = f._read_size()
        if first_size != second_size:
            raise IOError(f'Sizes do not agree in the header({first_size}) and footer({second_size}) for '
                        'this record - check header dtype')

#***********************************************************************
def read_ramses_100(repository):
    ''' This routine reads DM particles dumped in the RAMSES format.
    # NB: repository is directory containing output files
    # e.g. /horizon1/teyssier/ramses_simu/boxlen100_n256/output_00001/
    '''
    atexit.unregister(H.flush)
    signal.signal(signal.SIGINT, signal.SIG_DFL)#, H.flush)
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)#, H.flush)
    signal.signal(signal.SIGTERM, H.flush)
    # read cosmological params in header of amr file
    ipos    = repository.find("output_")
    nchar   = repository[ipos+7:ipos+12]
    nomfich = f"{repository}/amr_{nchar}.out00001"
    with FortranFile(nomfich, 'r') as f:
        H.ncpu, = f.read_ints()
        H.ndim, = f.read_ints()   
        nx,ny,nz = f.read_ints()
        H.nlevelmax, = f.read_ints()
        ngridmax, = f.read_ints()
        nstep_coarse, = f.read_ints()
        boxlen, = f.read_reals()
        # temps conforme tau, expansion factor, da/dtau
        tco,aexp_ram,hexp = f.read_reals()
        omega_m,omega_l,omega_k,omega_b = f.read_reals()
        # to get units cgs multiply by these scale factors
        scale_l,scale_d,scale_t = f.read_reals()
    # use approximate comv from cm to Mpc to match Romain's conversion... 
    H.Lboxp          = boxlen*scale_l/3.08e24/aexp_ram # converts cgs to Mpc comoving
    #write(errunit,*) 'af,hf,lboxp,ai,aexp',af,h_f,lboxp,ai,aexp_ram
    H.aexp           = aexp_ram*H.af  
    H.omega_f        = omega_m+omega_b
    H.omega_lambda_f = omega_l
    H.omega_c_f      = omega_k
    H.Lf             = boxlen*scale_l/3.08e24/aexp_ram
    H.mboxp          = np.float64(2.78782)*(H.Lf**3)*(H.H_f/100.)**2*H.omega_f
    print(f"> From AMR file: `{nomfich}`")
    print(f">     ncpu={str(H.ncpu):>6}, ndim={H.ndim:1d}, nstep_coarse={nstep_coarse:6d}")
    print(f">     nlevelmax={H.nlevelmax}, ngridmax={ngridmax}")
    print(f">     t={tco:.3e}, aexp={aexp_ram:.3e}, hexp={hexp:.3e}")
    print(f">     omega_m={omega_m:.3f}, omega_l={omega_l:.3f}, omega_k={omega_k:.3f}, omega_b={omega_b:.3f}")
    print(f">     boxlen={H.Lboxp:.3e} h-1 Mpc")
    print(boxlen*scale_l/3.08e24)
    print(boxlen*scale_l/3.08e24/aexp_ram)
    print(boxlen*scale_l/3.08e24*aexp_ram)
    stop()

    # now read the particle data files
    nomfich = f"{repository}/part_{nchar}.out00001"
    print(f"\n> From part file: `{nomfich}`")
    with FortranFile(nomfich, 'r') as f:
        H.ncpu, = f.read_ints()
        H.ndim, = f.read_ints()

    H.npart = 0
    for icpu1 in range(1,H.ncpu+1):
       nomfich = f"{repository}/part_{nchar}.out{icpu1:05d}"
       with FortranFile(nomfich, 'r') as f:
           ncpu2, = f.read_ints()
           ndim2, = f.read_ints()
           npart2, = f.read_ints()
       H.npart = H.npart+npart2
    
    H.nusedpart = H.npart
    print(f"> Found {H.npart} particles")
    print(f"> Reading positions and masses...")
    
    H.allocate('pos_10', (H.npart, H.ndim), dtype=np.float64)
    H.allocate('vel_10', (H.npart, H.ndim), dtype=np.float64)
    H.allocate('mass_10', (H.npart,), dtype=np.float64)
    H.massalloc = True
  
    iterobj = range(1,H.ncpu+1)
    if(H.TQDM)and(H.megaverbose):
        iterobj = tqdm(range(1,H.ncpu+1), desc="Reading particles", unit="cpu", ncols=100)
    for icpu1 in iterobj:
        nomfich = f"{repository}/part_{nchar}.out{icpu1:05d}"
        with FortranFile(nomfich, 'r') as f:
            ncpu2, = f.read_ints()
            ndim2, = f.read_ints()
            npart2, = f.read_ints()
            tmpp = np.empty((npart2, ndim2), dtype=np.float64)
            tmpv = np.empty((npart2, ndim2), dtype=np.float64)
            tmpm = np.empty(npart2, dtype=np.float64)
            idp = np.empty(npart2, dtype=np.int32)
            
            # read all particle positions
            for idim0 in range(H.ndim):
                tmpp[:,idim0] = f.read_reals()
            # read all particle velocities
            for idim0 in range(H.ndim):
                tmpv[:,idim0] = f.read_reals()
            # read all particle masses
            tmpm[:] = f.read_reals()
            # read all particle ids
            idp[:] = f.read_ints()

        # now sort DM particles in ascending id order
        for idim0 in range(H.ndim):
            # put all positions between -0.5 and 0.5
            mem['pos_10'][idp-1,idim0] = tmpp[:,idim0] - 0.5
            # convert code units to km/s 
            mem['vel_10'][idp-1,idim0] = tmpv[:,idim0]*scale_l/scale_t*1e-5
            mem['mass_10'][idp-1] = tmpm[:]
        del tmpp; del tmpv; del tmpm; del idp
    
        mtot = np.sum(mem['mass_10'])
        # that is for the dark matter so let's add baryons now if there are any 
        # and renormalization flag is on ##
        massres = np.min(mem['mass_10'])*H.mboxp*1e11
        H.massp   = np.min(mem['mass_10'])
        print(f"> particle mass (in M_sun)               = {massres}")
        if(H.RENORM):
            massres /= mtot
            H.massp /= mtot
            print(f"> particle mass (in M_sun) after renorm  = {massres}")
        if(H.BIG_RUN):
            H.deallocate('mass_10')


def _ra4_initial_mass_skip(repository):
    descriptor = os.path.join(repository, "part_file_descriptor.txt")
    if not os.path.exists(descriptor):
        return None

    fields = []
    with open(descriptor, "r") as stream:
        for line in stream:
            if line.lstrip().startswith("#"):
                continue
            columns = [column.strip() for column in line.split(",")]
            if len(columns) >= 2:
                fields.append(columns[1])

    if "initial_mass" not in fields:
        return None
    if "metallicity" not in fields:
        raise RuntimeError(
            f"Ra4 descriptor has initial_mass but no metallicity: {descriptor}"
        )
    metal_index = fields.index("metallicity")
    initial_mass_index = fields.index("initial_mass")
    if initial_mass_index <= metal_index:
        raise RuntimeError(
            f"Unsupported Ra4 initial_mass record order: {descriptor}"
        )
    return initial_mass_index - metal_index - 1


#***********************************************************************
def _read_ramses_new_1010(icpu, kwargs):
    repository = kwargs['repository']
    rver = kwargs['rver']
    nchar = kwargs['nchar']
    H.ndim = kwargs['ndim']
    scale_l = kwargs['scale_l']
    scale_t = kwargs['scale_t']
    dmcount = kwargs['dmcount']

    nomfich = f"{repository}/part_{nchar}.out{icpu:05d}"
    with FortranFile(nomfich, 'r') as f:
        ncpu2, = f.read_ints()
        ndim2, = f.read_ints()
        npart2, = f.read_ints()
        skip_records(f, 1)
        nstar, = f.read_ints()
        skip_records(f, 2)
        nsink, = f.read_ints()
        # assert nsize[icpu-1] == npart2

        if dmcount:
            skip_records(f,H.ndim)
            skip_records(f,H.ndim)
            skip_records(f,1)
        else:
            tmpp = np.empty((npart2,3), dtype=np.float64)
            tmpv = np.empty((npart2,3), dtype=np.float64)
            tmpm = np.empty(npart2, dtype=np.float64)
            
            # read all particle positions
            for idim0 in range(H.ndim):
                tmpp[:,idim0] = f.read_reals()
            # read all particle velocities
            for idim0 in range(H.ndim):
                tmpv[:,idim0] = f.read_reals()
            # read all particle masses
            tmpm = f.read_reals()
        # read all particle ids
        idp = f.read_ints()
        # read grid level of particles
        skip_records(f, 1)
        tmpt = None
        tmpmetal = None
        tmpm0 = None
        if(rver=='Ra4'):
            # read particle family
            fam = f.read_ints(dtype=np.int8)
            # read particle tag
            skip_records(f, 1)
        if((nstar>0)or(nsink>0)):
            if(rver != 'Ra4' or not dmcount):
                tmpt = f.read_reals()
            if not dmcount:
                tmpmetal = f.read_reals()
                if kwargs['m0_skip'] is not None:
                    skip_records(f, kwargs['m0_skip'])
                    tmpm0 = f.read_reals()
            elif rver != 'Ra4':
                skip_records(f, 1)
    # now sort DM particles in ascending id order and get rid of stars
    if(rver=='Ra4'):
        dmmask = fam==1
        if dmcount:
            mask = dmmask # DM particles only
        else:
            starmask = (fam==2)
            mask = dmmask | starmask
    else:
        if tmpt is None:
            tmpt = np.zeros(npart2, dtype=np.float64)
        dmmask = (idp>0)&(tmpt==0)
        if dmcount:
            mask = dmmask
        else:
            starmask = ((tmpt < 0) & (idp > 0)) | ((tmpt != 0) & (idp < 0))
            mask = dmmask | starmask
    if not dmcount:
        idp = np.abs(idp)
        star_slots = idp[starmask] - 1
        idp = np.where(starmask, idp+H.ndm, idp)
    npart_tmp = np.sum(mask)
    if not dmcount:
        ind = idp[mask]-1
        pos_tmp_101 = H.maccess('pos_tmp_101')
        for idim0 in range(H.ndim):
            # put all positions between -0.5 and 0.5
            pos_tmp_101[ind,idim0] = tmpp[mask,idim0]-0.5
            # convert code units to km/s 
            mem['vel_tmp_101'][ind,idim0] = tmpv[mask,idim0]*scale_l/scale_t*1e-5
        mem['mass_tmp_101'][ind] = tmpm[mask]
        if starmask.any():
            age_10 = H.maccess('age_10')
            metal_10 = H.maccess('metal_10')
            age_10[star_slots] = (
                kwargs['snapshot_age'] - H.epoch_to_age(tmpt[starmask]))
            metal_10[star_slots] = tmpmetal[starmask]
            if tmpm0 is not None:
                H.maccess('m0_10')[star_slots] = tmpm0[starmask]
    return npart_tmp
#***********************************************************************
import time
def read_ramses_new_101(repository, rver='Ra3'):
    ''' This routine reads DM particles dumped in the RAMSES format.
    # NB: repository is directory containing output files
    # e.g. /horizon1/teyssier/ramses_simu/boxlen100_n256/output_00001/
    '''
    atexit.unregister(H.flush)
    signal.signal(signal.SIGINT, signal.SIG_DFL)#, H.flush)
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)#, H.flush)
    signal.signal(signal.SIGTERM, H.flush)
    if(H.verbose): print()
    if(H.verbose): print(f"\t------------------------------------------------------------------")
    if(H.verbose): print(f"\t| Reading RAMSES version {rver}  ")
    if(H.verbose): print(f"\t------------------------------------------------------------------")
    m0_skip = _ra4_initial_mass_skip(repository) if rver == 'Ra4' else None

    # read cosmological params in header of amr file
    ipos    = repository.find("output_")
    nchar   = repository[ipos+7:ipos+12]
    nomfich = f"{repository}/amr_{nchar}.out00001"
    with FortranFile(nomfich, 'r') as f:
        H.ncpu, = f.read_ints()
        H.ndim, = f.read_ints()   
        nx,ny,nz = f.read_ints()
        H.nlevelmax, = f.read_ints()
        ngridmax, = f.read_ints()
        skip_records(f, 2) # nboundary, ngrid_current
        boxlen, = f.read_reals()
        nout,idum,idum = f.read_ints()
        skip_records(f, 2) # tout, aout
        tco, = f.read_reals()
        skip_records(f, 2) # dtold, dtnew
        idum, nstep_coarse = f.read_ints() # nstep, nstep_coarse
        skip_records(f, 1) # einit, mass_tot_0, rho_tot
        temp = f.read_reals()
        omega_m,omega_l,omega_k,omega_b,dummy = temp[:5]
        temp = f.read_reals()
        aexp_ram, hexp = temp[:2]
    if(H.verbose): print(f"\t|> From AMR file: `{nomfich}`")
    if(H.verbose): print(f"\t|>     ncpu={H.ncpu:6d}, ndim={H.ndim:1d}, nstep_coarse={nstep_coarse:6d}")
    if(H.verbose): print(f"\t|>     nlevelmax={H.nlevelmax:3d}, ngridmax={ngridmax:8d}")
    if(H.verbose): print(f"\t|>     t={tco:.3E}, aexp={aexp_ram:.3E}, hexp={hexp:.3E}")
    if(H.verbose): print(f"\t|>     omega_m={omega_m:.3f}, omega_l={omega_l:.3f}, omega_k={omega_k:.3f}, omega_b={omega_b:.3f}")

    nomfich = f"{repository}/info_{nchar}.txt"
    with open(nomfich, 'r') as f:
        for line in f:
            line = line.strip()  # Remove leading and trailing whitespace
            if not line or line.startswith('#'):
                continue  # Skip empty lines and lines starting with '#'

            # Split the line at the '=' character
            name, value = map(str.strip, line.split('='))

            # Check for a comment at the end of the value
            if '!' in value:
                value = value.split('!')[0].strip()

            # Process the name and value
            if name in ('unit_l', 'scale_l'):
                scale_l = np.float64(value)
            elif name in ('unit_d', 'scale_d'):
                scale_d = np.float64(value)
            elif name in ('unit_t', 'scale_t'):
                scale_t = np.float64(value)
                break

    H.Lboxp          = boxlen*scale_l/np.float64(3.08e24)/aexp_ram # converts cgs to Mpc comoving
    H.aexp           = aexp_ram*H.af  
    H.omega_f        = omega_m # Override
    H.omega_lambda_f = omega_l # Override
    H.omega_c_f      = omega_k # Override
    if(H.verbose): print(f"\t|>     boxlen={boxlen*scale_l/np.float64(3.08e24):.3e} h-1 Mpc")
    H.Lf = boxlen*scale_l/3.08e24/aexp_ram # Override
    H.mboxp     = np.float64(2.78782)*(H.Lf**3)*(H.H_f/100.)**2*H.omega_f # Override
    H.set_cosmology()
    snapshot_age = H.epoch_to_age(tco)

    # now read the particle data files
    nomfich = f"{repository}/part_{nchar}.out00001"
    if(H.verbose): print(f"\t|> From Part file: `{nomfich}`")
    with FortranFile(nomfich, 'r') as f:
        H.ncpu, = f.read_ints()
        H.ndim, = f.read_ints()

    H.npart = 0
    # nsize = np.zeros(H.ncpu, dtype=np.int32)
    for icpu1 in range(1,H.ncpu+1):
        nomfich = f"{repository}/part_{nchar}.out{icpu1:05d}"
        with FortranFile(nomfich, 'r') as f:
            ncpu2, = f.read_ints()
            ndim2, = f.read_ints()
            npart2, = f.read_ints()
            idum = f.read_ints()
            nstar, = f.read_ints()
            idum = f.read_ints()
            idum = f.read_ints()
            nsink, = f.read_ints()
        H.npart += npart2

    if(H.verbose): print(f"\t|> Found {H.npart} Total particles")
    H.nstar = nstar
    H.nusedpart = H.npart
    if(H.verbose): print(f"\t|        {H.npart-H.nstar} other particles")
    if(H.verbose): print(f"\t|        {nstar} star particles")
    if(H.verbose): print(f"\t|> Reading positions, velocities and masses...")
    H.allocate('pos_tmp_101', (H.nusedpart, H.ndim), dtype=np.float64)
    H.allocate('vel_tmp_101', (H.nusedpart, H.ndim), dtype=np.float64)
    H.allocate('mass_tmp_101', (H.nusedpart,), dtype=np.float64)
    H.massalloc = True
  
    # Count only DM particles
    kwargs = {
        'repository':repository, 'rver':rver, 'nchar':nchar,
        'ndim':H.ndim, 'scale_l':scale_l, 'scale_t':scale_t,
        'snapshot_age':snapshot_age, 'dmcount':True,
        'm0_skip':m0_skip}
    iterobj = range(1,H.ncpu+1)
    if(H.nbPes==1): # Sequential reading
        if(H.TQDM)and(H.megaverbose): pbar = tqdm(total=H.ncpu, desc=f"\t|  Count DMs (nbPes={H.nbPes})", unit="cpu", file=sys.stdout, disable=(not H.megaverbose), ncols=100)
        ndm = 0
        for icpu1 in iterobj:
            ndm += _read_ramses_new_1010(icpu1, kwargs)
            if(H.TQDM)and(H.megaverbose): pbar.update(1)
        if(H.TQDM)and(H.megaverbose): pbar.close()
    else: # Multiprocessing
        signal.signal(signal.SIGTERM, signal.SIG_DFL)
        with Pool(processes=H.nbPes) as pool:
            async_results = []
            for icpu1 in iterobj:
                r = pool.apply_async(_read_ramses_new_1010, (icpu1, kwargs))
                async_results.append((icpu1, r))
            ndm = 0
            for icpu1, r in tqdm(async_results, total=H.ncpu, desc=f"\t|  Count DMs (nbPes={H.nbPes})", unit="cpu", disable=(not H.megaverbose), ncols=100):
                try:
                    ndm += r.get(timeout=300)  # 300 sec
                except TimeoutError:
                    print(f"\n[HANG?] icpu={icpu1} still not finished after 300s")
                    raise
        signal.signal(signal.SIGTERM, H.flush)
    H.ndm =ndm
    if(H.verbose): print(f"\t|> Found {H.ndm} DM particles after masking")

    H.nusedpart = H.ndm + H.nstar
    if H.nstar > 0:
        H.allocate('age_10', (H.nstar,), dtype=np.float64)
        H.allocate('metal_10', (H.nstar,), dtype=np.float64)
        mem['age_10'][:] = np.nan
        mem['metal_10'][:] = np.nan
        if m0_skip is not None:
            H.allocate('m0_10', (H.nstar,), dtype=np.float64)
            mem['m0_10'][:] = np.nan
    # Read all parts
    kwargs['dmcount'] = False
    iterobj = range(1,H.ncpu+1)
    if(H.nbPes==1): # Sequential reading
        if(H.TQDM)and(H.megaverbose): pbar = tqdm(total=H.ncpu, desc=f"\t|  Reading parts(nbPes={H.nbPes})", unit="cpu", file=sys.stdout, ncols=100)
        npart = 0
        for icpu1 in iterobj:
            npart += _read_ramses_new_1010(icpu1, kwargs)
            if(H.TQDM)and(H.megaverbose): pbar.update(1)
        if(H.TQDM)and(H.megaverbose): pbar.close()
    else: # Multiprocessing
        signal.signal(signal.SIGTERM, signal.SIG_DFL)
        with Pool(processes=H.nbPes) as pool:
            async_results = []
            for icpu1 in iterobj:
                r = pool.apply_async(_read_ramses_new_1010, (icpu1, kwargs))
                async_results.append((icpu1, r))
            npart = 0
            for icpu1, r in tqdm(async_results, total=H.ncpu, desc=f"\t|  Reading parts(nbPes={H.nbPes})", unit="cpu", disable=(not H.megaverbose), ncols=100):
                try:
                    npart += r.get(timeout=300)  # 300 sec
                except TimeoutError:
                    print(f"\n[HANG?] icpu={icpu1} still not finished after 300s")
                    raise
        signal.signal(signal.SIGTERM, H.flush)
    H.npart = npart
    if(H.verbose): print(f"\t|> Reading parts done", flush=True)
    if(H.verbose): print(f"\t|> Found {H.npart} (DM+Star) particles after masking")
    if H.nstar > 0:
        n_age = np.sum(np.isfinite(mem['age_10']))
        n_metal = np.sum(np.isfinite(mem['metal_10']))
        n_m0 = (
            np.sum(np.isfinite(mem['m0_10']) & (mem['m0_10'] > 0))
            if H.allocated('m0_10') else H.nstar
        )
        if n_age != H.nstar or n_metal != H.nstar or n_m0 != H.nstar:
            raise RuntimeError(
                f"Incomplete stellar properties: age={n_age}, metal={n_metal}, "
                f"m0={n_m0}, "
                f"expected={H.nstar}")
    if(H.verbose and H.nstar > 0):
        print(f"\t|> Stellar age range (Gyr)          = "
              f"{np.min(mem['age_10']):.6g} .. {np.max(mem['age_10']):.6g}")
        print(f"\t|> Stellar metallicity range        = "
              f"{np.min(mem['metal_10']):.6g} .. {np.max(mem['metal_10']):.6g}")
        if H.allocated('m0_10'):
            print(f"\t|> Stellar initial-mass range        = "
                  f"{np.min(mem['m0_10']):.6g} .. {np.max(mem['m0_10']):.6g}")
    H.allocate('pos_10', (H.npart, H.ndim), dtype=np.float64)
    H.allocate('vel_10', (H.npart, H.ndim), dtype=np.float64)
    H.allocate('mass_10', (H.npart,), dtype=np.float64)
    # H.allocate('id_10', (H.npart,), dtype=np.int32)
    mem['pos_10'][:H.npart, :] = mem['pos_tmp_101'][:H.npart, :]
    mem['vel_10'][:H.npart, :] = mem['vel_tmp_101'][:H.npart, :]
    mem['mass_10'][:H.npart] = mem['mass_tmp_101'][:H.npart]
    # mem['id_10'][:] = np.arange(1,H.npart+1)
    H.deallocate('pos_tmp_101','vel_tmp_101','mass_tmp_101')

    mtot = np.sum(mem['mass_10'])
    # that is for the dark matter so let's add baryons now if there are any 
    # and renormalization flag is on ##
    dmmassres = np.min(mem['mass_10'][:H.ndm])*H.mboxp*1e11
    starmassres = np.min(mem['mass_10'][H.ndm:])*H.mboxp*1e11
    H.massp   = np.min(mem['mass_10'][:H.ndm])
    if(H.verbose): print(f"\t|> DM particle mass (in M_sun)               = {dmmassres}")
    if(H.verbose): print(f"\t|> Star particle mass (in M_sun)             = {starmassres}")
    if(H.RENORM):
        massres /= mtot
        H.massp /= mtot
        if(H.verbose): print(f"\t> particle mass (in M_sun) after renorm  = {massres}")
    if(H.BIG_RUN):
        H.deallocate('mass_10')

    if(H.verbose): print(f"\t------------------------------------------------------------------\n", flush=True)

#***********************************************************************
def write_tree_brick_1d():
    '''
    This subroutine writes the information relevant to building a halo 
    merging tree (using the build_tree program) i.e. for each halo:
      1/ the list of all the particles it contains (this enables us --- as  
         particle numbers are time independent --- to follow the halo history) 
      2/ its properties which are independent of merging history (mass ...)
    '''
    import os
    nchar   = f'{int(H.file_num):05d}'
    if(H.BIG_RUN):
        if(H.write_resim_masses):
            f44 = FortranFile(f'{H.output_dir}/resim_masses{H.prefix}.dat', 'w')
            f44.write_record(H.nusedpart)
            f44.write_record(mem['mass_10'])
            f44.close()
            full_path = os.path.abspath(f'{H.output_dir}/resim_masses{H.prefix}.dat')
            os.chmod(full_path, H.fchmod); os.chown(full_path, H.uid, H.gid)
            H.write_resim_masses = False

    whereIam_idxs = mem['whereIam_idxs']
    whereIam_counts = mem['whereIam_counts']
    pids0_groupsorted = mem['pids0_groupsorted']
    if H.dump_dms:
        mass_10 = mem['mass_10']
        pos_10 = mem['pos_10']
        vel_10 = mem['vel_10']

    if(not H.fsub):
        filename = f"{H.output_dir}/tree_brick_{nchar}{H.prefix}"
    else:
        filename = f"{H.output_dir}/tree_bricks{nchar}{H.prefix}"
    f44 = FortranFile(filename, 'w')
    print()
    print('> Output data to build halo merger tree to: ',filename)
    f44.write_record(H.nusedpart)
    f44.write_record(H.massp)
    f44.write_record(H.aexp)
    f44.write_record(H.omega_t)
    f44.write_record(H.age_univ)
    f44.write_record(H.nb_of_halos, H.nb_of_subhalos)
    # nb_of_parts_o0_1 = mem['nb_of_parts_o0_1']
    # first_part_oo_1 = mem['first_part_oo_1']
    timerecord = [0 for _ in range(10)]
    iterator = range(H.nb_of_halos + H.nb_of_subhalos)
    if H.verbose:
        iterator = tqdm(iterator, desc="Writing halos", unit="halo", ncols=100)
    for i0 in iterator:
        ith=0; ref = time.time()
        me = H.liste_halos_o0[i0+1]
        timerecord[ith]+=time.time()-ref; ith+=1; ref = time.time()
        # nmem = nb_of_parts_o0_1[i0+1]
        nmem = whereIam_counts[i0+1]
        timerecord[ith]+=time.time()-ref; ith+=1; ref = time.time()
        idx = whereIam_idxs[i0+1]
        timerecord[ith]+=time.time()-ref; ith+=1; ref = time.time()
        count = whereIam_counts[i0+1]
        timerecord[ith]+=time.time()-ref; ith+=1; ref = time.time()
        indexps = pids0_groupsorted[idx:idx+count]
        timerecord[ith]+=time.time()-ref; ith+=1; ref = time.time()
        # write list of particles in each halo
        if(H.dump_dms):
            mass_memb = mass_10[indexps]
            pos_memb = pos_10[indexps]
            vel_memb = vel_10[indexps]
        members = indexps + 1 # particle ids are 1-based in Fortran
        timerecord[ith]+=time.time()-ref; ith+=1; ref = time.time()
        f44.write_record(nmem)
        f44.write_record(members)
        timerecord[ith]+=time.time()-ref; ith+=1; ref = time.time()

        if(H.dump_dms):
            ncharg = f"{me['id']:07d}"
            nomfich = f"HAL_{nchar}/halo_dms_{ncharg}"
            f9 = FortranFile(nomfich, 'w')
            f9.write_record(me['id'])
            f9.write_record(me['level'])
            f9.write_record(me['m'])
            f9.write_record(me['px'],me['py'],me['pz'])
            f9.write_record(me['vx'],me['vy'],me['vz'])
            f9.write_record(me['Lx'],me['Ly'],me['Lz'])
            f9.write_record(nmem)
            for idim0 in range(H.ndim):
                # mdump[:]=pos_memb[:,idim0]
                f9.write_record( pos_memb[:,idim0] )
            for idim0 in range(H.ndim):
                # mdump[:]=vel_memb[:,idim0]
                f9.write_record( vel_memb[:,idim0] )
            f9.write_record( mass_memb )
            f9.write_record( members )
            del mass_memb; del pos_memb; del vel_memb; del mdump
            f9.close()
            full_path = os.path.abspath(nomfich)
            os.chmod(full_path, H.fchmod); os.chown(full_path, H.uid, H.gid)

        # del members
        # write each halo properties [MAIN BOTTLENECK]
        # write_halo_1d0(me,f44)
        f44.write_record( me['id'] )
        f44.write_record( me['timestep']  )
        f44.write_record( me['level'],me['hosthalo'],me['hostsub'],me['nbsub'],me['nextsub'] )
        f44.write_record( me['m'] )
        f44.write_record( me['px'],me['py'],me['pz'] )
        f44.write_record( me['vx'],me['vy'],me['vz'] )
        f44.write_record( me['Lx'],me['Ly'],me['Lz']  )
        f44.write_record( me['r'], me['sha'], me['shb'], me['shc'] )
        f44.write_record( me['ek'],me['ep'],me['et'] )
        f44.write_record( me['spin'] )
        f44.write_record( me['sigma'] )
        f44.write_record( me['rvir'],me['mvir'],me['tvir'],me['cvel'] )
        f44.write_record( me['rho_0'],me['r_c'],me['cNFW'] )
        f44.write_record( me['mcontam'] )
        timerecord[ith]+=time.time()-ref; ith+=1; ref = time.time()
    f44.close()
    full_path = os.path.abspath(filename)
    os.chmod(full_path, H.fchmod); os.chown(full_path, H.uid, H.gid)


#***********************************************************************
def write_tree_brick_hdf():
    '''
    This subroutine writes the information relevant to building a halo 
    merging tree (using the build_tree program) i.e. for each halo:
      1/ the list of all the particles it contains (this enables us --- as  
         particle numbers are time independent --- to follow the halo history) 
      2/ its properties which are independent of merging history (mass ...)
    '''
    import os, h5py
    nchar   = f'{int(H.file_num):05d}'
    if(H.BIG_RUN):
        if(H.write_resim_masses):
            f44 = FortranFile(f'{H.output_dir}/resim_masses{H.prefix}.dat', 'w')
            f44.write_record(H.nusedpart)
            f44.write_record(mem['mass_10'])
            f44.close()
            full_path = os.path.abspath(f'{H.output_dir}/resim_masses{H.prefix}.dat')
            os.chmod(full_path, H.fchmod); os.chown(full_path, H.uid, H.gid)
            H.write_resim_masses = False

    whereIam_idxs = mem['whereIam_idxs']
    whereIam_counts = mem['whereIam_counts']
    pids0_groupsorted = mem['pids0_groupsorted']
    if H.dump_dms:
        mass_10 = mem['mass_10']
        pos_10 = mem['pos_10']
        vel_10 = mem['vel_10']

    if(not H.fsub):
        filename = f"{H.output_dir}/tree_brick_{nchar}{H.prefix}.h5"
    else:
        filename = f"{H.output_dir}/tree_bricks{nchar}{H.prefix}.h5"
    # f44 = FortranFile(filename, 'w')
    print()
    print('> Output data to build halo merger tree to: ',filename)
    with h5py.File(filename, 'w') as f44:
        #---------------------------------
        # HEADER
        #---------------------------------
        f44.create_group('header')
        header = f44['header']
        # Snapshot data
        header.attrs['npart'] = H.npart
        header.attrs['nusedpart'] = H.nusedpart
        header.attrs['ndm'] = H.ndm
        header.attrs['nstar'] = H.nstar
        header.attrs['massp'] = H.massp
        header.attrs['aexp'] = H.aexp
        header.attrs['omega_t'] = H.omega_t
        header.attrs['age_univ'] = H.age_univ
        header.attrs['boxsize2']=H.boxsize2
        header.attrs['box_comoving_mpc'] = H.boxsize2
        header.attrs['box_physical_mpc'] = H.aexp * H.boxsize2
        header.attrs['units_version'] = UNITS_VERSION
        header.attrs['hubble']=H.hubble
        header.attrs['mboxp']=H.mboxp
        # HaloMaker data
        header.attrs['nb_of_halos'] = H.nb_of_halos
        header.attrs['nb_of_subhalos'] = H.nb_of_subhalos
        # User data
        header.attrs['createtime'] = H.mprefix[2]

        #---------------------------------
        # INPUT
        #---------------------------------
        # input_HaloMaker.dat
        f44.create_group('input')
        finput = f44['input']
        finput.attrs['omega_f']=H.omega_f
        finput.attrs['omega_lambda_f']=H.omega_lambda_f
        finput.attrs['af']=H.af
        finput.attrs['Lf']=H.Lf
        finput.attrs['H_f']=H.H_f
        finput.attrs['FlagPeriod']=H.FlagPeriod
        finput.attrs['nMembers']=H.nMembers
        finput.attrs['cdm']=H.cdm
        finput.attrs['method']=H.method
        finput.attrs['b_init']=H.b_init
        finput.attrs['nvoisins']=H.nvoisins
        finput.attrs['nhop']=H.nhop
        finput.attrs['rho_threshold']=H.rho_threshold
        finput.attrs['fudge']=H.fudge
        finput.attrs['fudgepsilon']=H.fudgepsilon
        finput.attrs['alphap']=H.alphap
        finput.attrs['verbose']=H.verbose
        finput.attrs['megaverbose']=H.megaverbose
        finput.attrs['DPMMC']=H.DPMMC
        finput.attrs['SC']=H.SC
        finput.attrs['dcell_min']=H.dcell_min
        finput.attrs['eps_SC']=H.eps_SC
        finput.attrs['nsteps']=H.nsteps
        finput.attrs['dump_dms']=H.dump_dms
        #---------------------------------
        # Catalog
        #---------------------------------
        box_physical_mpc = H.aexp * H.boxsize2
        cat = _convert_halo_catalog_units(H.liste_halos_o0[1:], box_physical_mpc)
        grp = f44.create_group('catalog')
        halo_ds = grp.create_dataset('halo', shape=cat.shape, dtype=cat.dtype, data=cat, compression='lzf')
        halo_ds.attrs['field_units'] = _json_attr(CATALOG_FIELD_UNITS)
        halo_ds.attrs['units_version'] = UNITS_VERSION

        from ssp_photometry import MODELS, model_metadata
        photometry = f44.create_group('photometry')
        for requested_model in MODELS:
            model, metadata = model_metadata(requested_model)
            model_group = photometry.create_group(model)
            for key, value in metadata.items():
                model_group.attrs[key] = value
            model_group.attrs['mass_source'] = (
                'initial_mass' if H.allocated('m0_10')
                else 'current_mass_fallback'
            )
            model_group.attrs['interpolation'] = 'bilinear in log-age/log-metallicity, magnitude space'
            model_group.attrs['frame'] = 'intrinsic rest-frame'
            model_group.attrs['dust'] = False
            model_group.attrs['nebular_emission'] = False
            model_group.attrs['SDSS_system'] = 'AB'
            model_group.attrs['Johnson_system'] = 'Vega'
            model_group.attrs['row_alignment'] = '/catalog/halo'
            model_group.attrs['fields'] = ','.join(H.photometry_ssp_dtype.names)
            photo = _convert_photometry_units(
                mem[f'photometry_{model.lower()}'][1:],
                box_physical_mpc,
            )
            photo_ds = model_group.create_dataset(
                'data', shape=photo.shape, dtype=photo.dtype,
                data=photo, compression='lzf',
            )
            photo_ds.attrs['field_units'] = _json_attr(PHOTOMETRY_FIELD_UNITS)
            photo_ds.attrs['units_version'] = UNITS_VERSION

        #---------------------------------
        # Member
        #---------------------------------
        # cumsum index, 1d member IDs, ...
        grp = f44.create_group('member')
        grp.attrs['index_base'] = 1
        grp.attrs['layout'] = 'flat'
        grp.create_dataset('index', data=whereIam_idxs, compression='lzf')
        grp.create_dataset('count', data=whereIam_counts, compression='lzf')
        pids = pids0_groupsorted+1 # 0-based to 1-based
        
        member_chunk = min(max(1, len(pids)), 65536)
        grp.create_dataset('pids', data=pids, compression='lzf', chunks=(member_chunk,))
        if H.dump_dms:
            grp.attrs['fields'] = 'pids,pos,vel,mass'
            pos_ds = grp.create_dataset(
                'pos', shape=(len(pids), 3), dtype=pos_10.dtype,
                compression='lzf', chunks=(member_chunk, 3))
            vel_ds = grp.create_dataset(
                'vel', shape=(len(pids), 3), dtype=vel_10.dtype,
                compression='lzf', chunks=(member_chunk, 3))
            mass_ds = grp.create_dataset(
                'mass', shape=(len(pids),), dtype=mass_10.dtype,
                compression='lzf', chunks=(member_chunk,))
            pos_ds.attrs['unit'] = 'code_unit'
            pos_ds.attrs['coordinate_frame'] = 'RAMSES box [0,1), periodic'
            vel_ds.attrs['unit'] = 'km/s peculiar'
            mass_ds.attrs['unit'] = 'Msun'
            mass_ds.attrs['units_version'] = UNITS_VERSION
            pos_ds.attrs['units_version'] = UNITS_VERSION
            vel_ds.attrs['units_version'] = UNITS_VERSION
            iterator = range(0, len(pids), member_chunk)
            if H.verbose:
                iterator = tqdm(iterator, desc="Writing flat halo members", unit="chunk", ncols=100)
            for start in iterator:
                end = min(start + member_chunk, len(pids))
                src = pids0_groupsorted[start:end]
                pos_ds[start:end, :] = np.mod(pos_10[src, :] / box_physical_mpc, 1.0)
                vel_ds[start:end, :] = vel_10[src, :]
                mass_ds[start:end] = mass_10[src] * 1.0e11

    full_path = os.path.abspath(filename)
    os.chmod(full_path, H.fchmod); os.chown(full_path, H.uid, H.gid)

#***********************************************************************
def write_halo_1d0(h:np.void,unitfile:FortranFile):
    # Masses (h['m'],h['mvir']) are in units of 10^11 Msol, and 
    # Lengths (h['x'],h['y'],h['z'],h['r'],h['rvir']) are in units of Mpc
    # Velocities (h['vx'],h['vy'],h['vz'],h['cvel']) are in km/s
    # Energies (h['ek'],h['ep'],h['et'])
    # Temperatures (h['tvir']) are in K
    # Angular Momentum (h['Lx'],h['Ly'],h['Lz']) are in
    # Other quantities are dimensionless (h['id'],h['timestep'],h['spin'])  

    unitfile.write_record( h['id'] )
    unitfile.write_record( h['timestep']  )
    unitfile.write_record( h['level'],h['hosthalo'],h['hostsub'],h['nbsub'],h['nextsub'] )
    unitfile.write_record( h['m'] )
    unitfile.write_record( h['px'],h['py'],h['pz'] )
    unitfile.write_record( h['vx'],h['vy'],h['vz'] )
    unitfile.write_record( h['Lx'],h['Ly'],h['Lz']  )
    unitfile.write_record( h['r'], h['sha'], h['shb'], h['shc'] )
    unitfile.write_record( h['ek'],h['ep'],h['et'] )
    unitfile.write_record( h['spin'] )
    unitfile.write_record( h['sigma'] )
    unitfile.write_record( h['rvir'],h['mvir'],h['tvir'],h['cvel'] )
    unitfile.write_record( h['rho_0'],h['r_c'],h['cNFW'] )
    unitfile.write_record( h['mcontam'] )
