from halo_defs import mem, frange
import halo_defs as H
import numpy as np
import atexit, signal, sys
from scipy.io import FortranFile
from tqdm import tqdm
# from multiprocessing import Pool
import multiprocessing as mp
ctx = mp.get_context('fork')
Pool = ctx.Pool

import faulthandler, signal, sys
faulthandler.enable()
faulthandler.register(signal.SIGUSR1, file=sys.stderr, all_threads=True)

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

    if(H.zoomin):
        refmask = mem['refmask_10']
        if np.sum(refmask) < H.nMembers:
            H.nusedpart = np.sum(refmask)
            print()
            print('ERROR in this snapshot:')
            print('> nusedpart=',H.nusedpart)
            print('> nMembers=',H.nMembers)
            print('> nusedpart < nMembers_threshold !')
            print('> Skip this step')
            if(H.allocated('pos_10')): H.deallocate('pos_10')
            if(H.allocated('vel_10')): H.deallocate('vel_10')
            if(H.allocated('mass_10')): H.deallocate('mass_10')
            if(H.allocated('whereIam_parts')): H.deallocate('whereIam_parts')
            if(H.allocated('refmask_10')): H.deallocate('refmask_10')
            if(H.allocated('id_10')): H.deallocate('id_10')
            if(len(H.liste_halos_o0)>0): H.liste_halos_o0 = np.empty(0, dtype=H.halo_dtype)
            return
    
    pos = mem['pos_10']
    print(f"> min max position (in box units)   : {np.min(pos)},{np.max(pos)}")
    if H.zoomin:
        # refmask = mem['refmask_10']
        rpos = pos[refmask]
        print(f">                  (zoom-in)        : {np.min(rpos)},{np.max(rpos)}")
    vel = mem['vel_10']
    print(f"> min max velocities (in km/s)      : {np.min(vel)},{np.max(vel)}")
    if H.zoomin:
        rvel = vel[refmask]
        print(f">                    (zoom-in)      : {np.min(rvel)},{np.max(rvel)}")
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
        if(rver=='Ra4'):
            # read particle family
            fam = f.read_ints(dtype=np.int8)
            # read particle tag
            skip_records(f, 1)
        else:
            # read all particle creation times if necessary
            if((nstar>0)or(nsink>0)):
                tmpt = f.read_reals()
                if(H.METALS):
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
        dmmask = (idp>0)&(tmpt==0)
        if dmcount:
            mask = dmmask
        else:
            starmask = ((tmpt < 0) & (idp > 0)) | ((tmpt != 0) & (idp < 0))
            mask = dmmask | starmask
    if not dmcount:
        idp = np.abs(idp)
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
    kwargs = {'repository':repository, 'rver':rver, 'nchar':nchar, 'ndim':H.ndim, 'scale_l':scale_l, 'scale_t':scale_t, 'dmcount':True}
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

    if (H.zoomin):
        # TODO: Reduce pos, vel, mass, id arrays to zoom-in region to save memory and speed up subsequent computations?
        if(H.verbose): print(f"\t|> Applying zoom-in mask to particles...")
        H.allocate('refmask_10', (H.npart,), dtype=np.bool_)
        goodmask = mem['mass_10'] < 1.0001*H.massp
        # if(H.verbose): print(f"\t|> Mass mask applied: {np.sum(goodmask)} particles kept out of {H.npart} ({100*np.sum(goodmask)/H.npart:.2f}%)", flush=True)
        if H.zoombox is None:
            goodpos = mem['pos_10'][goodmask]
            xmin, xmax = np.min(goodpos[:,0]), np.max(goodpos[:,0])
            ymin, ymax = np.min(goodpos[:,1]), np.max(goodpos[:,1])
            zmin, zmax = np.min(goodpos[:,2]), np.max(goodpos[:,2])
            H.zoombox = np.array([xmin,xmax,ymin,ymax,zmin,zmax])
        xmin,xmax,ymin,ymax,zmin,zmax = H.zoombox
        if(H.verbose): print(f"\t|> Zoom-in box (in box units): ")
        if(H.verbose): print(f"\t|    x=[{xmin:.3f}, {xmax:.3f}]")
        if(H.verbose): print(f"\t|    y=[{ymin:.3f}, {ymax:.3f}]")
        if(H.verbose): print(f"\t|    z=[{zmin:.3f}, {zmax:.3f}]")
        H.zoombox = np.array([xmin,xmax,ymin,ymax,zmin,zmax])
        # goodmask = goodmask & (mem['pos_10'][:,0] >= xmin) & (mem['pos_10'][:,0] <= xmax) & (mem['pos_10'][:,1] >= ymin) & (mem['pos_10'][:,1] <= ymax) & (mem['pos_10'][:,2] >= zmin) & (mem['pos_10'][:,2] <= zmax)
        mem['refmask_10'][:] = goodmask
        if(H.verbose): print(f"\t|> Zoom-in mask applied: {np.sum(goodmask)} particles kept out of {H.npart} ({100*np.sum(goodmask)/H.npart:.2f}%)", flush=True)
        

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
    if(H.dump_dms):
        #    call system('mkdir HAL_'//TRIM(nchar))
        os.mkdir(f'HAL_{nchar}')    
        full_path = os.path.abspath(f'HAL_{nchar}')
        os.chmod(full_path, H.dchmod); os.chown(full_path, H.uid, H.gid)

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
    if(H.dump_dms):
        #    call system('mkdir HAL_'//TRIM(nchar))
        os.mkdir(f'HAL_{nchar}')    
        full_path = os.path.abspath(f'HAL_{nchar}')
        os.chmod(full_path, H.dchmod); os.chown(full_path, H.uid, H.gid)

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
        finput.attrs['zoomin']=H.zoomin
        if H.zoomin:
            finput.attrs['zoombox']=H.zoombox

        #---------------------------------
        # Catalog
        #---------------------------------
        cat = H.liste_halos_o0[1:]
        grp = f44.create_group('catalog')
        grp.create_dataset('halo', shape=cat.shape, dtype=cat.dtype, data=cat, compression='lzf')

        #---------------------------------
        # Member
        #---------------------------------
        # cumsum index, 1d member IDs, ...
        grp = f44.create_group('member')
        grp.create_dataset('index', data=whereIam_idxs, compression='lzf')
        pids = pids0_groupsorted+1 # 0-based to 1-based
        
        grp.create_dataset('pids', data=pids, compression='lzf')
        if H.dump_dms:
            iterator = range(H.nb_of_halos + H.nb_of_subhalos)
            if H.verbose:
                iterator = tqdm(iterator, desc="Writing halo members", unit="halo", ncols=100)
            for i0 in iterator:
                idx = whereIam_idxs[i0+1]
                count = whereIam_counts[i0+1]
                indexps = pids[idx:idx+count] # 1-based
                # write list of particles in each halo
                if(H.dump_dms):
                    mass_memb = mass_10[indexps-1]
                    pos_memb = pos_10[indexps-1]
                    vel_memb = vel_10[indexps-1]

                mgrp=grp.create_group(f"{i0+1:07d}")
                mgrp.create_dataset('id', data=indexps, compression='gzip')
                mgrp.create_dataset('pos', data=pos_memb, compression='gzip')
                mgrp.create_dataset('vel', data=vel_memb, compression='gzip')
                mgrp.create_dataset('m', data=mass_memb, compression='gzip')

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
