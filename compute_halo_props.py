from input_output import *
import halo_defs as H
from halo_defs import mem,vector,frange
import time, os
import numpy as np
from scipy.io import FortranFile
from itertools import combinations
from ssp_photometry import (
    BANDS, MODELS, interpolation_coordinates, interpolate_magnitude,
    load_all_tables,
)

#//////////////////////////////////////////////////////////////////////////
#**************************************************************************
def init_0():
    # initialize gravitationalsoftening 
    initgsoft_00()

    # initialize cosmological and technical parameters of the simulation
    init_cosmo_01()

    print("Open inputfiles_HaloMaker.dat")
    nlines = 0
    with open('inputfiles_HaloMaker.dat', 'r') as f:
        for line in f:
            if line[0] == '#': continue
            if line[0] == '!': continue
            if (len(line.strip())==0): continue

            nlines += 1
    H.nsteps = nlines
    print(f'> Number of snapshots to analyze =   {H.nsteps}')

    return True

#*************************************************************************
def initgsoft_00():
    '''
    subroutine to initialize the required arrays for the gravitational 
    field smoothing interpolation in routine SOFTGRAV. This is only
    required for a cubic spline kernel; interpolation is performed in 
    distance.
    '''

    tiny = np.float64(1.e-19)
    one  = np.float64(1.0)
    two  = np.float64(2.0)

    if(H.gravsoft=='harmonic'): return None

    deldrg      = 2./H.ninterp
    H.phsmooth_oo[0] = 7.*np.sqrt(tiny)/5.

    for i0 in range(H.ninterp):
       xw  = (i0+1)*deldrg
       xw2 = xw*xw
       xw3 = xw2*xw
       xw4 = xw2*xw2
       if(xw <= one):
          H.phsmooth_oo[i0+1] = -2.*xw3*(one/3.-3.*xw2/20.+xw3/20.)+7.*xw/5.
       else:
          H.phsmooth_oo[i0+1] = -one/15.+8.*xw/5.-xw3*(4./3.-xw+3.*xw2/10.-xw3/30.)
       if (xw >= two):
          H.phsmooth_oo[i0+1] = one
    return None

#*************************************************************************
def init_cosmo_01():
    '''
    This routine reads in the `input_HaloMaker.dat` file which contains the cosmological 
    and technical parameters of the N-Body simulation to analyze.
    '''
    # cannot do otherwise than setting all the strings to a value larger than that of a line 
    # in the input file and trim them whenever it is needed

    # Initial (beginning of the simulation) expansion factor
    # H.ai             = np.float64(1)#/(1+50))
    cosmo_given={}
    H.ai             = np.float64(1/(1+50))
    # The file is based on a keyword list. Keywords are optional, so we set defaults here.  
    H.omega_f        = np.float64(0.3333)      #   mass density at final timestep
    cosmo_given['omega_f'] = False
    H.omega_lambda_f = np.float64(0.6667)      #   lambda at final timestep 
    cosmo_given['omega_lambda_f'] = False
    # H.af             = np.float64(36.587)      #   expansion factor of the final timestep  
    H.af             = np.float64(1)      #   expansion factor of the final timestep  
    H.Lf             = np.float64(150.0)       #   final length of box in physical Mpc 
    cosmo_given['Lf'] = False
    H.H_f            = np.float64(66.667)      #   Hubble constant in km/s/Mpc, at final timestep
    H.b_init         = np.float64(0.2)         #   linking length friend-of-friend parameter @ z=0
    H.nMembers       = 20          #   minimum number of particles for a halo
    H.nsteps         = 1           #   number of files to analyse (listed in 'files.dat')
    H.method         = "FOF"
    if(H.ANG_MOM_OF_R): H.agor_file      = "ang_mom_of_r"

    print('\n> Values of input parameters:  ')
    print('> ---------------------------- ')
    print("")
    print(f'> Looking for `input_HaloMaker.dat` in directory: `{H.output_dir}`')
    # Parse input_HaloMaker.dat: each key maps to an H attribute via
    # H.ALIAS_TO_ATTR / H.PARAMS (defined in halo_defs.py).
    f20 = open('input_HaloMaker.dat', 'r')

    for line in f20:
        i = line.find('=')
        if i == -1 or line.strip().startswith('#') or line.strip().startswith('!'):
            continue

        name  = line[:i].strip()
        value = line[i+1:].strip()

        # check for a comment at end of line
        i = value.find('!')
        if i != -1:
            value = value[:i].strip()

        print(f'>{name:>15} : {str(value):>10}')

        if name == 'agor_file':
            if H.ANG_MOM_OF_R:
                H.agor_file = f"{H.output_dir}/{value}"

        elif name in H.ALIAS_TO_ATTR:
            attr = H.ALIAS_TO_ATTR[name]
            dtype = H.PARAMS[attr][1]
            setattr(H, attr, dtype(value))
            try: cosmo_given[attr] = True
            except KeyError: pass

        elif name[0] in ["#", "!"]:
            pass

        elif name in H.IGNORED_INPUT_KEYS:
            pass
        else:
            print(f'dont recognise parameter: {name}')

    f20.close()


    # command-line override
    for attr, value in vars(H.args).items():
        if value is None:
            continue

        if attr == 'output_dir':
            continue

        if attr == 'agor_file':
            if H.ANG_MOM_OF_R:
                H.agor_file = f"{H.output_dir}/{value}"
                print(f'[override] H.agor_file = {H.agor_file}')
            continue
        elif attr == 'prefix':
            H.prefix = H.args.prefix
            if H.prefix[0] != "_":
                H.prefix = f"_{H.prefix}"
            print(f'[override] H.{attr} = {H.prefix}')
        elif attr in H.PARAMS:
            setattr(H, attr, value)
            print(f'[override] H.{attr} = {value}')


    header_keys = ('omega_f', 'omega_lambda_f', 'Lf')
    missing_header_values = [key for key in header_keys if not cosmo_given[key]]
    if missing_header_values:
        joined = ', '.join(missing_header_values)
        print(f'[postprocess] {joined} not provided; RAMSES header values will be used after read_data')
    cosmo_inputs_ready = cosmo_given['omega_f'] and cosmo_given['omega_lambda_f']
    box_input_ready = cosmo_given['Lf']
    # ==========================================================
    # ==========================================================
    # ==========================================================

    if cosmo_inputs_ready:
        if( ((H.omega_f+H.omega_lambda_f) != 1.0)and(H.omega_lambda_f != 0.0) ):
            raise NotImplementedError('> lambda + non flat Universe not implemented yet')

    # In the most general of cases:
    #     af/aexp         = 1+z
    #     omega(z)        = (H_f/H(z))**2 * (1+z)**3  * omega_f
    #     omega_lambda(z) = omega_lambda_f * (H_f/H(z))**2
    #     omega_c(z)      = omega_c_f * (1+z)**2 * (H_f/H(z))**2
    #     H(z)**2         = H_f**2*( omega_f*(1+z)**3 + omega_c_f*(1+z)**2 + omega_lambda_f)

        H.omega_c_f = 1. - H.omega_f - H.omega_lambda_f
        H.H_i       = H.H_f * np.sqrt( H.omega_f*(H.af/H.ai)**3 + H.omega_c_f*(H.af/H.ai)**2 + H.omega_lambda_f)

    if cosmo_inputs_ready and box_input_ready:
        # initial size of the box in physical Mpc (NB: f index stands for final quantities)
        H.Lboxp = H.Lf*(H.ai/H.af)
        # rho_crit = 2.78782 h^2  (critical density in units of 10**11 M_sol/Mpc^3)
        H.mboxp = np.float64(2.78782)*(H.Lf**3)*(H.H_f/100.)**2*H.omega_f
        print()
        print( f'> Initial/Final values of parameters:  ')
        print( f'> -----------------------------------  ')
        print( f' At (initial) redshift              :   ',H.af/H.ai-1.)
        print( f' > box size (Mpc)                   :   ',H.Lboxp)
        print( f' > Hubble parameter  (km/s/Mpc)     :   ',H.H_i)
        # print( f' > box mass (10^11 Msol)            :   ',H.mboxp)
        print( f' At (final) redshift                :   ',H.af/H.af-1.)
        print( f' > box size (Mpc)                   :   ',H.Lf)
        print( f' > Hubble parameter  (km/s/Mpc)     :   ',H.H_f)
        print( f' box mass (10^11 Msol)              :   ',H.mboxp)
        print()
    elif cosmo_inputs_ready:
        print()
        print( f'> Initial/Final values of parameters:  ')
        print( f'> -----------------------------------  ')
        print( f' At (initial) redshift              :   ',H.af/H.ai-1.)
        print( f' > Hubble parameter  (km/s/Mpc)     :   ',H.H_i)
        print( f' At (final) redshift                :   ',H.af/H.af-1.)
        print( f' > Hubble parameter  (km/s/Mpc)     :   ',H.H_f)
        print( ' > box size and box mass are finalized after read_data')
        print()
    else:
        print()
        print( f'> Initial/Final values of parameters:  ')
        print( f'> -----------------------------------  ')
        print( ' > RAMSES header will set           :   ', joined)
        print( ' > box size, box mass, and snapshot cosmology are finalized after read_data')
        print()

    H.write_resim_masses = True
    H.mprefix[1] = f"u{os.getuid()}"
    H.mprefix[2] = time.strftime("t%Y%m%d_%H%M%S", time.localtime())
    if H.prefix != "":
        H.mprefix[2] = f"{H.mprefix[2]}{H.prefix}"


def _compute_halo_props(ih1, member, fagor, printdatacheckhalo):
    # member is a tuple of (count, indexps, mypos, myvel, mymass, mydensity)
    h = mem['liste_halos_o0'][ih1]
    h['nmem'] = member[0]
    if(printdatacheckhalo):
        print('> halo:', h['id'],'nb_of_parts_o0_1',h['nmem'])
        t0 = time.time()
    det_mass_17(h, member=member)
    det_center_18(h, member=member)
    refine_centers_1b(h, member=member)
    compute_ang_mom_19(h, member=member)
    r_halos_1a(h, member=member)
    det_vir_1b(h, fagor=fagor, member=member)
    compute_spin_parameter_1c(h)
    compute_stellar_1d(h, member=member)
    compute_ssp_photometry_1f(ih1, h, member=member)
    compute_extended_profiles_1e(h, member=member)

    if(printdatacheckhalo): print('> mass:', h['m'], 'mhalo:', h['mdm'], 'mstar:', h['m*'])
    if(printdatacheckhalo): print('> center:', h['px'],h['py'],h['pz'])
    if(printdatacheckhalo): print('> angular momentum:', h['Lx'],h['Ly'],h['Lz'])
    if(printdatacheckhalo): print('> angular momentum (stars):', h['Lx*'],h['Ly*'],h['Lz*'])
    if(printdatacheckhalo): print('> radius:', h['r'])
    if(printdatacheckhalo): print('> mvir,rvir:', h['mvir'], h['rvir'])
    if(printdatacheckhalo): print('> spin:', h['spin'])
    if(printdatacheckhalo):
        t1 = time.time()
        print('> halo computation took:',int(t1- t0) ,'s')
        print()

#*************************************************************************
def new_step_1():
    '''
    This is the main subroutine: it builds halos from the particle simulation and 
    computes their properties ...
    '''
    print(f'\n\n> Timestep  --->    {H.numero_step}')
    print('> -------------------')

    read_time_ini = time.time()
    timereport = []

    # read N-body info for this new step
    if(H.verbose): print(f"\n\n\n$$ Read data...", flush=True)
    _ref = time.time()
    read_data_10()
    if(H.nusedpart<H.nMembers):
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
        if(H.allocated('id_10')): H.deallocate('id_10')
        if(H.allocated('age_10')): H.deallocate('age_10')
        if(H.allocated('metal_10')): H.deallocate('metal_10')
        if(H.allocated('m0_10')): H.deallocate('m0_10')
        if(len(H.liste_halos_o0)>0): H.liste_halos_o0 = np.empty(0, dtype=H.halo_dtype)
        return
    if(H.verbose): print(f"\n$$ Read data done ({time.time()-_ref:.2f} sec)", flush=True)
    timereport.append(('- read_data', time.time()-_ref))

    # determine the age of the universe and current values of cosmological parameters
    if(H.verbose): print(f"\n$$ Determine the Age...", flush=True)
    _ref = time.time()
    det_age_11()
    if(H.verbose): print(f"\n$$ Determine the Age done ({time.time()-_ref:.2f} sec)", flush=True)
    timereport.append(('- det_age', time.time()-_ref))
    
    # first compute the virial overdensity (with respect to the average density of the universe) 
    # predicted by the top-hat collapse model at this redshift 
    if(H.verbose): print(f"\n$$ Compute the virial overdensity...", flush=True)
    _ref = time.time()
    virial_12()
    if(H.verbose): print(f"\n$$ Compute the virial overdensity done ({time.time()-_ref:.2f} sec)", flush=True)
    timereport.append(('- virial', time.time()-_ref))

    if(H.verbose): print(f"\n$$ Make halos...", flush=True)
    _ref = time.time()
    if(H.method!="FOF"):
       H.allocate('whereIam_parts',H.nusedpart,dtype=np.int32)
    timerecords = make_halos_13()  
    if(H.verbose): print(f"\n$$ Make halos done ({time.time()-_ref:.2f} sec)", flush=True)
    timereport.append(('- make_halos', time.time()-_ref))
    if timerecords is not None:
        for itr in timerecords:
            timereport.append(itr)

    if(H.Test_FOF):
        filelisteparts = f"whereIam_parts_{H.numstep}"
        with FortranFile(filelisteparts, 'w') as f:
            f.write_record(H.nusedpart,H.nb_of_halos)
            f.write_record(mem['whereIam_parts'][:H.nusedpart])
        full_path = os.path.abspath(filelisteparts)
        os.chmod(full_path, H.fchmod); os.chown(full_path, H.uid, H.gid)
        # reset nb of halos not to construct halos
        H.nb_of_halos = 0
    # if there are no halos go to the next timestep
    if (H.nb_of_halos == 0):
        print('no halos deallocating')
        H.deallocate('pos_10','vel_10')
        if(H.allocated('density_1312')): H.deallocate('density_1312')
        if(H.allocated('mass_10')): H.deallocate('mass_10')
        if(H.allocated('age_10')): H.deallocate('age_10')
        if(H.allocated('metal_10')): H.deallocate('metal_10')
        if(H.allocated('m0_10')): H.deallocate('m0_10')
        H.deallocate('whereIam_parts')
        return

    ### Move this to inside `make_linked_list` function
    # H.allocate('first_part_oo_1', H.nb_of_halos+H.nb_of_subhalos+1, dtype=np.int32)
    # H.allocate('nb_of_parts_o0_1', H.nb_of_halos+H.nb_of_subhalos+1, dtype=np.int32)
    # # make a linked list of the particles so that each member of a halo points to the next 
    # # until there are no more members (last particles points to -1)
    # H.allocate('linked_list_oo_1', 1+H.nusedpart+1, dtype=np.int32)

    if(H.verbose): print(f"\n$$ Make linked list...", flush=True)
    _ref = time.time()
    make_linked_list_14()

    # H.deallocate whereIam_parts bc it has been replaced by linked_list_oo
    H.deallocate('whereIam_parts')

    # in case of resimulation (or individual particle masses) count how many halos are contaminated 
    # (i.e. contain "low res" particles) 
    # NB: for numerical reasons, it can happen that mass is sligthly larger 
    #     than massp. Thus, for detecting LR particles, the criterion 
    #     if (mass(indexp) > massp) ...  may not work (it can get ALL parts, also HR parts!).
    #     Thus, to be conservative, use: if (mass(indexp) > massp * (1+1e-5))) ...
    whereIam_idxs = mem['whereIam_idxs']
    whereIam_counts = mem['whereIam_counts']
    pids0_groupsorted = mem['pids0_groupsorted']
    pos_10 = mem['pos_10']
    vel_10 = mem['vel_10']
    mass_10 = mem['mass_10']
    density_1312 = mem['density_1312']
    level_1319 = mem['level_1319']

    n_halo_contam = 0 
    n_subs_contam = 0  
    for i0 in range(H.nb_of_halos+H.nb_of_subhalos):
        idx = whereIam_idxs[i0+1]
        count = whereIam_counts[i0+1]
        indexps = pids0_groupsorted[idx:idx+count]
        mymass = mass_10[indexps]

        iscontam = (mymass > H.massp*1.00001).any()
        if(iscontam):
            if(H.fsub):
                if(level_1319[i0] == 1): 
                    n_halo_contam += 1
                else:
                    n_subs_contam += 1
            else:
                n_halo_contam += 1
    print('> # of halos, # of CONTAMINATED halos :',H.nb_of_halos,n_halo_contam,H.nb_of_subhalos,n_subs_contam)
    if(H.megaverbose):
        f222 = open(f'ncontam_halos{H.prefix}.dat', 'a+')
        f222.write(f'{H.numero_step:6d} {H.nb_of_halos:6d} {n_halo_contam:6d} {H.nb_of_subhalos:6d} {n_subs_contam:6d}\n')
        f222.close()
        full_path = os.path.abspath(f'ncontam_halos{H.prefix}.dat')
        os.chmod(full_path, H.fchmod); os.chown(full_path, H.uid, H.gid)

    if(H.verbose): print(f"\n$$ Make linked list done ({time.time()-_ref:.2f} sec)", flush=True)
    timereport.append(('- make_linked_list', time.time()-_ref))

    # until now we were using code units for positions and velocities
    # this routine changes that to physical (non comoving) coordinates for positions 
    # and peculiar (no Hubble flow) velocities in km/s
    # The masses are changed from code units into 10^11 M_sun as well.
    change_units_15()

    # allocation and initialization of the halo list
    H.allocate('liste_halos_o0', H.nb_of_halos+H.nb_of_subhalos+1, dtype=H.halo_dtype)
    H.liste_halos_o0 = mem['liste_halos_o0']
    for model in MODELS:
        photometry = H.allocate(
            f'photometry_{model.lower()}',
            H.nb_of_halos+H.nb_of_subhalos+1,
            dtype=H.photometry_ssp_dtype,
        )
        photometry['id'] = 0
        for field in H.photometry_ssp_fields:
            photometry[field] = np.nan
    load_all_tables()
    # (liste_halos_o0 is list of `halo` class, and defined in `halo_defs`)
    # (It may should be changed to shared memory if you want to implement the multiprocessing)

    init_halos_16()


    if(H.verbose): print(f"\n$$ Calculate halo properties...", flush=True)
    _ref = time.time()
    fagor=None
    if(H.ANG_MOM_OF_R):
        filename = f"{H.agor_file}.{H.numstep:03d}"
        fagor = FortranFile(filename, 'w')
        fagor.write_record(H.nb_of_halos,H.nb_of_subhalos)
        fagor.write_record(H.nshells)

    printdatacheckhalo = False
    # Let's make below loop parallel with multiprocessing, but for now, it is serial.
    iterator = frange(1,H.nb_of_halos + H.nb_of_subhalos)

    DEBUG=False
    if H.nbPes==1 or DEBUG:
        if H.megaverbose:
            iterator = tqdm(iterator, 
                total = H.nb_of_halos + H.nb_of_subhalos,
                desc = "[Ncpu=1] Calc halo props", ncols=100)
        for i1 in iterator:
            me = H.liste_halos_o0[i1]
            # member particles
            my_number = me['id']
            idx = whereIam_idxs[my_number]
            count = whereIam_counts[my_number]
            indexps = pids0_groupsorted[idx:idx+count]
            mypos = pos_10[indexps]
            myvel = vel_10[indexps]
            mymass = mass_10[indexps]
            mydensity = density_1312[indexps]
            member = (count, indexps, mypos, myvel, mymass, mydensity)
            
            _compute_halo_props(i1, member=member, fagor=fagor, printdatacheckhalo=printdatacheckhalo)
    else:
        callback=None
        if H.megaverbose:
            pbar = tqdm(total = H.nb_of_halos + H.nb_of_subhalos, desc = f"[Ncpu={H.nbPes}] Calc halo props", ncols=100)
            callback = lambda *args: pbar.update()
        signal.signal(signal.SIGTERM, signal.SIG_DFL)
        with Pool(processes=H.nbPes) as pool:
            async_results = []
            for i1 in iterator:
                me = H.liste_halos_o0[i1]
                my_number = me['id']
                idx = whereIam_idxs[my_number]
                count = whereIam_counts[my_number]
                indexps = pids0_groupsorted[idx:idx+count]
                mypos = pos_10[indexps]
                myvel = vel_10[indexps]
                mymass = mass_10[indexps]
                mydensity = density_1312[indexps]
                member = (count, indexps, mypos, myvel, mymass, mydensity)
                args = (i1, member, fagor, printdatacheckhalo)
                r = pool.apply_async(_compute_halo_props, args=args, callback=callback)
                async_results.append(r)
            for r in async_results:
                r.get()
        signal.signal(signal.SIGTERM, H.flush)
        if H.megaverbose: pbar.close()



    if(H.ANG_MOM_OF_R):
        fagor.close()
        full_path = os.path.abspath(filename)
        os.chmod(full_path, H.fchmod); os.chown(full_path, H.uid, H.gid)
    if(H.verbose): print(f"\n$$ Calculate halo properties done ({time.time()-_ref:.2f} sec).", flush=True)
    timereport.append(('- calc_halo_props', time.time()-_ref))

    if(H.verbose): print(f"\n$$ Write tree_bricks...", flush=True)
    _ref = time.time()
    # write_tree_brick_1d()
    write_tree_brick_hdf()
    if(H.verbose): print(f"\n$$ Write tree_bricks done ({time.time()-_ref:.2f} sec)", flush=True)
    timereport.append(('- write_tree_bricks', time.time()-_ref))

    # H.liste_halos_o0 = []
    H.liste_halos_o0 = np.empty(0, dtype=H.halo_dtype)
    H.deallocate('liste_halos_o0')
    for model in MODELS:
        H.deallocate(f'photometry_{model.lower()}')
    # H.deallocate('nb_of_parts_o0_1','first_part_oo_1','linked_list_oo_1')
    H.deallocate('pos_10','vel_10')
    H.deallocate('whereIam_idxs','whereIam_counts','pids0_groupsorted')
    if(H.allocated('mass_10')): H.deallocate('mass_10')
    if(H.allocated('age_10')): H.deallocate('age_10')
    if(H.allocated('metal_10')): H.deallocate('metal_10')
    if(H.allocated('m0_10')): H.deallocate('m0_10')
    if(not H.cdm): H.deallocate('density_1312')

    read_time_end = time.time()

    print('> time_step computations took : ',round(read_time_end - read_time_ini),' seconds', flush=True)
    if(H.megaverbose):
        ttotal = sum([t for _, t in timereport if t>=0])
        for name, t in timereport:
            if t >= 0:
                print(f"    {name:25s} : {t:.2f} sec ({t/ttotal*100:.2f}%)")
            else:
                print(f"    {name:25s} : {-t:.2f} sec")
        print()
    H.mlist()
    return True

#***********************************************************************
def make_halos_13():
    '''
    subroutine which builds the halos from particle data using fof or adaptahop
    '''
    import compute_neiKDtree_mod as neiKDtree

    print('> In routine make_halos ')
    print('> ----------------------')
    print(f"> npart={H.npart}")
    
    print( )
    print( '_______________________________________________________________________'  )
    print( )
    print( '          Compute neiKDtree'  )
    print( '          -----------------'  )
    print( )
    print( '_______________________________________________________________________'  )

    H.fPeriod[:]    = H.FlagPeriod
    if(H.FlagPeriod==1):
        print('> WARNING: Assuming PERIODIC boundary conditions --> make sure this is correct', flush=True)
        periodic = True
    else:
        print('> WARNING: Assuming NON PERIODIC boundary conditions --> make sure this is correct', flush=True)
        periodic = False

    if(H.numero_step==1):
        if(H.cdm):
            print('> Center of haloes and subhaloes are defined as the particle the closest to the cdm')
        else:
            print('> Center of haloes and subhaloes are defined as the particle with the highest density' )

        if(H.method == "FOF"):
            raise NotImplementedError('FOF is not implemented yet')
            print('> HaloMaker is using Friend Of Friend algorithm', flush=True)
            # fof_init()
            H.fsub = False
        elif(H.method == "HOP"):
            print('> HaloMaker is using Adaptahop in order to' )
            print('> Detect halos, and subhaloes will not be selected'    , flush=True)
            neiKDtree.init_adaptahop_130()
            H.fsub = False
        elif(H.method == "DPM"):
            print('> HaloMaker is using Adaptahop in order to' )
            print('> Detect halos, and subhaloes with the Density Profile Method'    , flush=True)
            neiKDtree.init_adaptahop_130()
            H.fsub = True
        elif(H.method == "MSM"):
            print('> HaloMaker is using Adaptahop in order to' )
            print('> Detect halos, and subhaloes with the Most massive Subhalo Method', flush=True)
            neiKDtree.init_adaptahop_130()
            H.fsub = True
        elif(H.method == "BHM"):
            print('> HaloMaker is using Adaptahop in order to' )
            print('> Detect halos, and subhaloes with the Branch History Method'    , flush=True)
            neiKDtree.init_adaptahop_130()
            H.fsub = True
        else:
            print('> Selection method: ',H.method,' is not included')
            print('> Please check input file:, input_HaloMaker.dat', flush=True)
    
    read_time_ini = time.time()
    timerecords = None
    if(H.method =='FOF'):
        # compute density when most dense particle is choosen as center
        raise NotImplementedError('FOF is not implemented yet')
        fof_main()
        if( (not H.cdm)and(H.nb_of_halos>0) ): compute_density_for_fof()
        H.nb_of_subhalos = 0       
    else:
        timerecords = neiKDtree.compute_adaptahop_131()
        # no nead to keep information about density when choosing particle closest to the cdm as center
        if(H.cdm): H.deallocate('density_1312')

    read_time_end = time.time()

    print('> Number of halos with more than     ', H.nMembers,' particles:',H.nb_of_halos)
    print('> Number of sub-halos with more than ', H.nMembers,' particles:',H.nb_of_subhalos)
    if(H.megaverbose): print('> time_step computations took : ',round(read_time_end - read_time_ini),' seconds', flush=True)
  
    return timerecords

#***********************************************************************
def make_linked_list_14():
    '''
    Subroutine builds a linked list of parts for each halo which 
    contains all its particles.
    '''
    # H.allocate('first_part_oo_1', H.nb_of_halos+H.nb_of_subhalos+1, dtype=np.int32)
    # H.allocate('nb_of_parts_o0_1', H.nb_of_halos+H.nb_of_subhalos+1, dtype=np.int32)
    # make a linked list of the particles so that each member of a halo points to the next 
    # until there are no more members (last particles points to -1)
    # H.allocate('linked_list_oo_1', 1+H.nusedpart+1, dtype=np.int32)

    # current_ptr_o0 = np.zeros(1+H.nb_of_halos+H.nb_of_subhalos, dtype=np.int32)-1
    # # initialization of linked list
    # mem['first_part_oo_1'][:]  = -1
    # mem['nb_of_parts_o0_1'][:] =  0
    # mem['linked_list_oo_1'][:] = -1
 
    # # make linked list: a few (necessary) explanations ....
    # #   1/ hindex1 (whereIam_parts(i)) is the number of the halo to which belongs particle i (i is in [1..nusedpart]) 
    # #   2/ first_part_oo(j) is the number of the first particle of halo j  (j is in [1..nhalo])
    # #   3/ current_ptr_o0(hindex1) is the number of the latest particle found in halo number hindex1 
    # # to sum up, first_part_oo(i) contains the number of the first particle of halo i,
    # # linked_list_oo(first_part_oo(i)) the number of the second particle of halo i, etc ... until 
    # # the last particle which points to number -1

    # for pindex0 in range(H.nusedpart):
    #     hindex1 = mem['whereIam_parts'][pindex0]
    #     if(hindex1>(H.nb_of_halos+H.nb_of_subhalos)): raise IndexError('error in whereIam_parts')
    #     if (mem['first_part_oo_1'][hindex1] == -1):
    #         mem['first_part_oo_1'][hindex1]  = pindex0
    #         mem['nb_of_parts_o0_1'][hindex1] = 1
    #         current_ptr_o0[hindex1] = pindex0
    #     else:
    #         ind              = current_ptr_o0[hindex1]
    #         mem['linked_list_oo_1'][ind] = pindex0
    #         current_ptr_o0[hindex1] = pindex0
    #         mem['nb_of_parts_o0_1'][hindex1] += 1

    # # close linked list
    # for i0 in range(0,H.nb_of_halos+H.nb_of_subhalos+1):
    #     if(current_ptr_o0[i0] == -1): continue
    #     index2              = current_ptr_o0[i0]
    #     mem['linked_list_oo_1'][index2] = -1
    return


#*************************************************************************
def init_halos_16():
    level_1319 = mem['level_1319']
    first_daughter_1319 = mem['first_daughter_1319']
    mother_1319 = mem['mother_1319']
    first_sister_1319 = mem['first_sister_1319']
    for ih1 in frange(1, H.nb_of_halos + H.nb_of_subhalos):
        me = H.liste_halos_o0[ih1]
        H.clear_halo(me)
        me['id'] = ih1
        me['timestep'] = H.numstep
        if(H.fsub):
            if(level_1319[ih1-1]==1):
                if(first_daughter_1319[ih1-1]>0):
                    me['nextsub'] = first_daughter_1319[ih1-1]
                me['hosthalo']   = ih1
            else:
                me['level']    = level_1319[ih1-1]
                me['hostsub'] =mother_1319[ih1-1]
                imother1 =me['hostsub']
                mother = H.liste_halos_o0[imother1]
                mother['nbsub'] += 1
                if(first_daughter_1319[ih1-1]>0):
                    me['nextsub'] = first_daughter_1319[ih1-1]
                elif(first_sister_1319[ih1-1]>0):
                    me['nextsub'] = first_sister_1319[ih1-1]
                else:
                    ihtmp1 = ih1
                    while((first_sister_1319[ihtmp1-1]<=0)and(level_1319[ihtmp1-1]>1)):
                        ihtmp1 =mother_1319[ihtmp1-1]
                    if(level_1319[ihtmp1-1]>1):
                        ihtmp1 = first_sister_1319[ihtmp1-1]
                        me['nextsub'] = ihtmp1
                imother1 = me['hostsub']
                imother1 = mother['hosthalo']
                mother = H.liste_halos_o0[imother1]
                if(mother['level']!=1): raise ValueError('wrong id for halo host')
                me['hosthalo'] = imother1
    
    if(H.fsub): H.deallocate('mother_1319','first_daughter_1319','first_sister_1319','level_1319')

#*************************************************************************
def det_mass_17(h:np.void, member:tuple):
    '''
    adds up masses of particles to get total mass of the halo. 
    '''

    # member = (count, indexps, mypos, myvel, mymass, mydensity)
    npch = member[0]
    imass = member[4]
    indexps = member[1]

    mtot = np.sum(imass)
    dmmask = indexps < H.ndm
    idm = imass[dmmask]
    masshalo = np.sum(idm)
    mcontam = np.sum(idm[idm > H.massp* 1.00001])
    istar = imass[~dmmask]
    massgal = np.sum(istar)

    h['ndm'] = len(idm)
    h['nstar'] = len(istar)
    h['m'] = mtot  # in 10^11 M_sun
    h['mdm'] = masshalo  # in 10^11 M_sun
    h['m*'] = massgal  # in 10^11 M_sun
    h['mcontam'] = mcontam

#***********************************************************************
def compute_ang_mom_19(h:np.void, member:tuple):
    '''
    compute angular momentum of all halos
    we compute r * m * v, where r & v are pos and vel of halo particles relative to center of halo
    (particle closest to center of mass or most dense particle)
    '''
    # (count, indexps, mypos, myvel, mymass, mydensity)
    indexps, ipos, ivel, imass = member[1], member[2], member[3], member[4]
    dmmask = indexps < H.ndm

    drxs = correct_for_periodicity_1d(ipos[:,0] - h['px'])
    drys = correct_for_periodicity_1d(ipos[:,1] - h['py'])
    drzs = correct_for_periodicity_1d(ipos[:,2] - h['pz'])
    pvxs = imass * (ivel[:,0] - h['vx'])
    pvys = imass * (ivel[:,1] - h['vy'])
    pvzs = imass * (ivel[:,2] - h['vz'])
    lxs = drys*pvzs - drzs*pvys
    lys = drzs*pvxs - drxs*pvzs
    lzs = drxs*pvys - drys*pvxs
    h['Lx'] = np.sum(lxs)
    h['Ly'] = np.sum(lys)
    h['Lz'] = np.sum(lzs)
    if((~dmmask).any() and np.isfinite(h['px*'])):
        spos = ipos[~dmmask]
        svel = ivel[~dmmask]
        smass = imass[~dmmask]
        svx = np.sum(smass * svel[:,0]) / h['m*']
        svy = np.sum(smass * svel[:,1]) / h['m*']
        svz = np.sum(smass * svel[:,2]) / h['m*']
        sdrx = correct_for_periodicity_1d(spos[:,0] - h['px*'])
        sdry = correct_for_periodicity_1d(spos[:,1] - h['py*'])
        sdrz = correct_for_periodicity_1d(spos[:,2] - h['pz*'])
        spvx = smass * (svel[:,0] - svx)
        spvy = smass * (svel[:,1] - svy)
        spvz = smass * (svel[:,2] - svz)
        h['Lx*'] = np.sum(sdry*spvz - sdrz*spvy)
        h['Ly*'] = np.sum(sdrz*spvx - sdrx*spvz)
        h['Lz*'] = np.sum(sdrx*spvy - sdry*spvx)
    else:
        h['Lx*'] = np.float64(np.nan)
        h['Ly*'] = np.float64(np.nan)
        h['Lz*'] = np.float64(np.nan)


#***********************************************************************
def r_halos_1a(h:np.void, member:tuple):
    '''
    compute distance of the most remote particle (with respect to center of halo, which
    is either center of mass or most bound particle)
    '''
    # (count, indexps, mypos, myvel, mymass, mydensity)
    indexps = member[1]
    ipos = member[2]
    drxs = correct_for_periodicity_1d(ipos[:,0] - h['px'])
    drys = correct_for_periodicity_1d(ipos[:,1] - h['py'])
    drzs = correct_for_periodicity_1d(ipos[:,2] - h['pz'])
    dr2s = drxs**2 + drys**2 + drzs**2
    dr2max = np.max(dr2s)

    h['r'] = np.sqrt(dr2max)
    starmask = indexps >= H.ndm
    if(starmask.any() and np.isfinite(h['px*'])):
        spos = ipos[starmask]
        drxs = correct_for_periodicity_1d(spos[:,0] - h['px*'])
        drys = correct_for_periodicity_1d(spos[:,1] - h['py*'])
        drzs = correct_for_periodicity_1d(spos[:,2] - h['pz*'])
        h['r*'] = np.sqrt(np.max(drxs**2 + drys**2 + drzs**2))
    else:
        h['r*'] = np.float64(np.nan)

#***********************************************************************
def correct_for_periodicity(dr:vector, copy=False):
    '''
    subroutine corrects for the fact that if you have periodic boundary conditions,
    then groups of particles that sit on the edge of the box can have part of their
    particles on one side and part on the other. So we have to take out a box-length
    when measuring the distances between group members if needed.
    '''

    if (H.FlagPeriod == 0): return  #--> NO PERIODIC BCs 
    
    if(copy):
        temp = vector()
        temp.x = dr.x*1
        temp.y = dr.y*1
        temp.z = dr.z*1
        if (temp.x > + H.Lbox_pt2): temp.x = temp.x - H.Lbox_pt
        if (temp.x < - H.Lbox_pt2): temp.x = temp.x + H.Lbox_pt 

        if (temp.y > + H.Lbox_pt2): temp.y = temp.y - H.Lbox_pt
        if (temp.y < - H.Lbox_pt2): temp.y = temp.y + H.Lbox_pt 

        if (temp.z > + H.Lbox_pt2): temp.z = temp.z - H.Lbox_pt
        if (temp.z < - H.Lbox_pt2): temp.z = temp.z + H.Lbox_pt 
        return temp
    else:
        if (dr.x > + H.Lbox_pt2): dr.x = dr.x - H.Lbox_pt
        if (dr.x < - H.Lbox_pt2): dr.x = dr.x + H.Lbox_pt 

        if (dr.y > + H.Lbox_pt2): dr.y = dr.y - H.Lbox_pt
        if (dr.y < - H.Lbox_pt2): dr.y = dr.y + H.Lbox_pt 

        if (dr.z > + H.Lbox_pt2): dr.z = dr.z - H.Lbox_pt
        if (dr.z < - H.Lbox_pt2): dr.z = dr.z + H.Lbox_pt 

def correct_for_periodicity_1d(dr:np.ndarray, copy=False):
    if (H.FlagPeriod == 0): return dr #--> NO PERIODIC BCs

    # In-place masked correction (avoids the full-size temporaries that two
    # np.where calls allocate). copy=False mutates the passed array; all callers
    # pass throwaway deltas such as (ipos[:,0] - h['px']), so this is safe.
    out = dr.copy() if copy else dr
    out[out > + H.Lbox_pt2] -= H.Lbox_pt
    out[out < - H.Lbox_pt2] += H.Lbox_pt
    return out


#***********************************************************************                
def det_center_18(h:np.void, member:tuple):
    '''
    compute position of center of mass of halo, and its velocity.
    '''
    # (count, indexps, mypos, myvel, mymass, mydensity)

    icenter = -1
    pos_10 = mem['pos_10']
    indexps = member[1]
    dmmask = indexps < H.ndm
    ipos = member[2] # shape (N, 3)
    ivel = member[3] # shape (N, 3)
    imass = member[4]
    
    hmy_number = h['id']
    hm = h['m']
    if(H.cdm):
        # compute cdm
        pc = vector()
        pos_10_ifirst = ipos[0]
        drxs = correct_for_periodicity_1d(ipos[:,0] - pos_10_ifirst[0])
        drys = correct_for_periodicity_1d(ipos[:,1] - pos_10_ifirst[1])
        drzs = correct_for_periodicity_1d(ipos[:,2] - pos_10_ifirst[2])
        # It's simply the center of mass.
        pc.x = np.sum(imass * drxs) / hm + pos_10_ifirst[0]
        pc.y = np.sum(imass * drys) / hm + pos_10_ifirst[1]
        pc.z = np.sum(imass * drzs) / hm + pos_10_ifirst[2]
        correct_for_periodicity(pc)
        # search particule closest to the cdm
        distmin = H.Lbox_pt
        drxs = correct_for_periodicity_1d(ipos[:,0] - pc.x)
        drys = correct_for_periodicity_1d(ipos[:,1] - pc.y)
        drzs = correct_for_periodicity_1d(ipos[:,2] - pc.z)
        drs2 = drxs**2 + drys**2 + drzs**2
        min_index = np.argmin(drs2)
        if drs2[min_index] < distmin*distmin:
            cpos = ipos[min_index]
    else:
        dens = member[5]
        max_index = np.argmax(dens)
        cpos = ipos[max_index]

    h['px']  = cpos[0] # in Mpc
    h['py']  = cpos[1] # in Mpc
    h['pz']  = cpos[2] # in Mpc

    # velocity of center is set equal velocity of center of mass:
    vcx = np.sum(imass * ivel[:,0])
    vcy = np.sum(imass * ivel[:,1])
    vcz = np.sum(imass * ivel[:,2])
    
    h['vx'] = vcx/hm # in km/s
    h['vy'] = vcy/hm # in km/s
    h['vz'] = vcz/hm # in km/s

    # vmean=np.sqrt(vcx**2+vcy**2+vcz**2)/hm
    # vnorms = np.sqrt(ivel[:,0]**2 + ivel[:,1]**2 + ivel[:,2]**2)
    # sigma2 = np.sum(imass * (vnorms - vmean)**2)

    # h['sigma'] = np.sqrt(sigma2/hm)
    dv2 = imass * ((ivel[:,0] - h['vx'])**2 + (ivel[:,1] - h['vy'])**2 + (ivel[:,2] - h['vz'])**2)
    sigma2 = np.sum(dv2)
    h['sigma'] = np.sqrt(sigma2 / hm)
    h['sigma_dm'] = np.sqrt(np.sum(dv2[dmmask]) / h['mdm']) if(h['mdm']>0) else 0
    h['sigma*'] = np.sqrt(np.sum(dv2[~dmmask]) / h['m*']) if(h['m*']>0) else 0

#***********************************************************************
def tab_props_inside_1b40(h:np.void,nr:int,v, member=None):
    '''
    returns the cumulative mass contained in concentric ellipsoids centered on the center of the
    halo (cdm or mbp)
    '''
    tabm2_o0 = np.zeros(nr); tabk2_o0 = np.zeros(nr)
    epsilon = 1e-2
    from num_rec import rf
    count, _, ipos, ivel, imass, _ = member

    # rescale to get ellipsoid  concentric to principal ellipsoid
    # which contains all the particles of the halo
    drxs = correct_for_periodicity_1d(ipos[:,0] - h['px']) # [Mpc]
    drys = correct_for_periodicity_1d(ipos[:,1] - h['py']) # [Mpc]
    drzs = correct_for_periodicity_1d(ipos[:,2] - h['pz']) # [Mpc]
    # project position vector along the principal ellipsoid axis
    dras = drxs*v[0,0] + drys*v[1,0] + drzs*v[2,0] # [Mpc]
    drbs = drxs*v[0,1] + drys*v[1,1] + drzs*v[2,1] # [Mpc]
    drcs = drxs*v[0,2] + drys*v[1,2] + drzs*v[2,2] # [Mpc]
    r_ells = np.sqrt((dras / h['sha'])**2 + (drbs / h['shb'])**2 + (drcs / h['shc'])**2) # [dimensionless]
    rmax = np.max(r_ells) # [dimensionless]

    amax = rmax * h['sha'] * (1.0 + epsilon) # [Mpc]
    bmax = rmax * h['shb'] * (1.0 + epsilon) # [Mpc]
    cmax = rmax * h['shc'] * (1.0 + epsilon) # [Mpc]

    # initialize loop quantities
    tabm2_o0[:] = 0
    tabk2_o0[:] = 0
    louped_parts = 0

    # compute velocities in the halo frame adding in the Hubble flow
    vtxs = ivel[:,0] - h['vx'] + drxs*H.Hub_pt # [km/s]
    vtys = ivel[:,1] - h['vy'] + drys*H.Hub_pt # [km/s]
    vtzs = ivel[:,2] - h['vz'] + drzs*H.Hub_pt # [km/s]
    v2 = vtxs**2 + vtys**2 + vtzs**2
    # biggest ellipsoid is divided in nr concentric ellipsoid shells: we
    # calculate below the ellipsoid bin in which each particle falls and fill up the 
    # mass and energy tables accordingly
    # NB: if by chance the most distant particle from the center is ON the shortest
    #     axis, then r_ell is equal to 1-epsilon (one of the terms below is one and the others 0)
    #     otherwise r_ell is between 0 and 1-epsilon and so we just multiply it by nr to 
    #     find the ellipsoid shell containing the particle.
    r_ells = np.sqrt((dras / amax)**2 + (drbs / bmax)**2 + (drcs / cmax)**2) # [dimensionless]
    i_ells = (r_ells*nr).astype(int) # [dimensionless]
    valid = i_ells < nr
    np.add.at(tabm2_o0, i_ells[valid], imass[valid])
    np.add.at(tabk2_o0, i_ells[valid], 0.5*imass[valid]*v2[valid])
    louped_parts = np.sum(~valid)

    if (louped_parts  >  0):
       raise ValueError('> Problem in tab_props_inside : missed ',louped_parts,' particles\n')

    np.cumsum(tabm2_o0, out=tabm2_o0)
    np.cumsum(tabk2_o0, out=tabk2_o0)
    # approximation based on appendix B of paper GALICS 1:
    # better than 10-15 . accuracy on average
    tabp2_o0 = -3/5 * H.gravconst * tabm2_o0**2 * rf(h['sha']**2,h['shb']**2,h['shc']**2)

    # correct potential energy estimate for small halos which are calculated by direct summation 
    if (h['ep'] != tabp2_o0[nr-1]): tabp2_o0 = tabp2_o0/tabp2_o0[nr-1]*h['ep']
    return tabm2_o0,tabk2_o0,tabp2_o0,amax,bmax,cmax

#***********************************************************************

def det_vir_props_1b4(h:np.void,v,amax=0,bmax=0,cmax=0,ttab=1000, member=None):
    '''
    computes the virial properties (radius, mass) of a halo
    '''
    # compute properties inside ttab concentric principal ellipsoids centered on center of halo
    tab_mass_o0,tab_ekin_o0,tab_epot_o0,amax,bmax,cmax = tab_props_inside_1b40(h,ttab,v, member=member)


    # find the outermost ellipsoid bin where virial theorem is either satisfied better than 20 .
    # or satisfied best ... 
    mvir      = tab_mass_o0[ttab-1]
    # initialize rvir to be the geometric average of the axis radii of the outermost ellipsoid shell
    # which contains at least one particle
    for i1 in frange(ttab-1,1,-1):
        if (tab_mass_o0[i1-1] < tab_mass_o0[i1]): break

    avir      = i1/(ttab-1)*amax
    bvir      = i1/(ttab-1)*bmax
    cvir      = i1/(ttab-1)*cmax
    rvir      = (avir*bvir*cvir)**(1./3.)
    # assume initial departure from virialization is 100 .
    virth_old = 1.0
    virth = 100.0
    virths = np.abs((2.0*tab_ekin_o0+tab_epot_o0)/(tab_ekin_o0+tab_epot_o0))
    good_virth = virths <= 0.2
    if good_virth.any():
        i0 = np.where(good_virth)[0][-1]
        virth = virths[i0]
        if virth < virth_old:
            mvir      = tab_mass_o0[i0] 
            # take the min here bc initialization throws away all the empty outer shells
            avir      = min(avir,i0/(ttab-1)*amax)
            bvir      = min(bvir,i0/(ttab-1)*bmax)
            cvir      = min(cvir,i0/(ttab-1)*cmax)
            rvir      = (avir*bvir*cvir)**(1./3.)
            # kvir      = tab_ekin_o0[i0] 
            # pvir      = tab_epot_o0[i0]
    for i0 in frange(ttab-1,0,-1):
        # if region is unbound, it cannot be virialized in the same time 
        if( (tab_ekin_o0[i0]+tab_epot_o0[i0]) >= 0.0): continue
        # the region is bound so compute relative virialization |2*K+P|/|K+P|
        virth = abs((2.0*tab_ekin_o0[i0]+tab_epot_o0[i0])/(tab_ekin_o0[i0]+tab_epot_o0[i0]))
        # if region is better virialized then update virial quantities
        if (virth < virth_old):
            mvir      = tab_mass_o0[i0] 
            # take the min here bc initialization throws away all the empty outer shells
            avir      = min(avir,i0/(ttab-1)*amax)
            bvir      = min(bvir,i0/(ttab-1)*bmax)
            cvir      = min(cvir,i0/(ttab-1)*cmax)
            rvir      = (avir*bvir*cvir)**(1./3.)
            # kvir      = tab_ekin_o0[i0] 
            # pvir      = tab_epot_o0[i0]
            virth_old = virth
            # if virial theorem holds with better than 20 . accuracy, exit do loop
            if (virth <= 0.20): break

 
    # for small halos it may happen that virial theorem is not enforced to within 15%.
    # bc the halo is not fully virialized yet or that it is valid by fluke (right combination
    # of potential and kinetic energy) ... so .... 
    # 1/ in the latter case, we further check that the halo density is high enough 
    # (order of what predicted by the spherical top-hat model, vir_overdens) 
    # for us to believe in the measurement of the virial theorem. 
    # 2/ in the former case, as an inner region measurement of the virialization would be too noisy for 
    # lack of mass resolution (not enough particles) we only use the overdensity criterion.
    # NB: this criterion is similar to NFW but with vir_overdens * average density and NOT 200. * critical density
    #     which does not makes sense when Omega_matter(z) != 1 (not critical) and we take ellipsoids
    #     NOT spheres bc they yield more accurate volumes and average halo densities
    # average density of the universe at current timestep (in 10^11 M_sun / Mpc^3)
 
    # volume of the smallest concentric ellipsoid
    volmin   = 4./3.*H.pi*(amax/(ttab-1))*(bmax/(ttab-1))*(cmax/(ttab-1))
    if( (virth > 0.20)or(mvir < H.vir_overdens*H.rho_mean*4./3.*H.pi*avir*bvir*cvir)):
        for ii1 in frange(ttab-1,1,-1):
            # assume that the virial mass and radii are obtained when the density inside the ellipsoid 
            # is greater than vir_overdens * rho_mean AND there is at least one particle inside the outermost
            # ellipsoid shell
            mvir = H.vir_overdens * H.rho_mean * volmin * ii1**3
            #mvir = 200d0 * 3d0*Hub_pt**2/8d0/acos(-1d0)/gravconst * volmin * ii1**3
            if( (tab_mass_o0[ii1] >= mvir)and(tab_mass_o0[ii1-1] < tab_mass_o0[ttab-1]) ): break
        mvir   = tab_mass_o0[ii1]
        # kvir   = tab_ekin_o0[ii1]
        # pvir   = tab_epot_o0[ii1]
        avir   = ii1/(ttab-1)*amax
        bvir   = ii1/(ttab-1)*bmax
        cvir   = ii1/(ttab-1)*cmax
        rvir   = (avir*bvir*cvir)**(1./3.)

    # check if virialization conditions were met --> if not set relevant quantities to zero: this is a 
    # non-virialized halo ....
    if( (mvir > 0.0)and(rvir > 0.0) ):
       # it may happen (although it is very rare bc vector linking center of halo to most distant 
       # particle has to be roughly parallel to minor axis of the halo in such a case) that the virial 
       # radius estimated from the geometric mean of the 3 principal axis is slightly larger than the 
       # distance of the most distant particle to the center of the halo (what we call r)
       # when this happens we set the virial radius to be r
       h['rvir']         = min(rvir,h['r'])    # in Mpc
       h['mvir']         = mvir             # in 10^11 M_sun
       # circular velocity at r_vir in km/s
       h['cvel']         = np.sqrt(H.gravconst*h['mvir']/h['rvir']) 
       # temperature at r_vir in K
       h['tvir']         = 35.9*H.gravconst*h['mvir']/h['rvir']
       # compute halo density profile within the virialized region
       compute_halo_profile_1b41(h)
    else:
       print('halo bugged (ID, Mvir, Rvir)',h['id'],mvir,rvir)
       raise ValueError('at ',h['px'],h['py'],h['pz'])
    return amax,bmax,cmax

#***********************************************************************
def compute_halo_profile_1b41(h:np.void):
    if (H.profile == 'TSIS'):
       # for the singular isothermal sphere the profile is defined @ rvir for it is singular at r=0.
       h['rho_0'] = h['mvir'] / (4.0 * H.pi * h['rvir']**3)
       h['r_c']   = h['rvir']
    else:
       raise NotImplementedError('Other profiles than TSIS not yet fully implemented')

#***********************************************************************
def det_vir_1b(h:np.void, fagor:FortranFile=None, member=None):
    # (count, indexps, mypos, myvel, mymass, mydensity)
    '''
    determine virial properties of the halos, energies and profiles
    '''
    # compute principal axis of halo
    d,v = det_main_axis_1b2(h, member=member) # d [Mpc], v [dimensionless]
 
    # compute halo energies if necessary i.e. in the case where the center of the halo is the 
    # center of mass bc if it is the most bound particle then energies are already available  
    det_halo_energies_1b3(h, member=member)

    # compute virial properties based on conditions necessary for the virial theorem to apply
    amax,bmax,cmax = det_vir_props_1b4(h,v, member=member)
    if(H.ANG_MOM_OF_R):
        det_ang_momentum_per_shell_1b5(h,amax,bmax,cmax,v,fagor=fagor, member=member)
        

#***********************************************************************
def det_age_11():
    '''
    subroutine which determines the current age of the Universe (Gyr) and
    values of current cosmological parameters 
    '''
    from num_rec import qromo, age_temps
    from functools import partial
    import scipy.integrate as integrate

    omm  = H.omega_f
    oml  = H.omega_lambda_f
    # age0 is the age of the universe @ the beginning of the simulation (in Gyr)
    # somme0 = qromo(age_temps,(H.af/H.ai),10001,omm=omm,oml=oml)
    func = partial(age_temps, oml=oml, omm=omm)
    somme0 = integrate.quad(func, (H.af/H.ai), 10001)[0]
    age0 = 977.78*somme0/H.H_f            

    # age1 is the time between the beginning of the simulation and the current timestep (in Gyr)
    if(H.aexp != H.ai):
    #    somme1 = qromo(age_temps,(H.af/H.aexp),(H.af/H.ai),omm=omm,oml=oml)
       somme1 = integrate.quad(func, (H.af/H.aexp), (H.af/H.ai))[0]
    else:
       somme1 = 0.0

    age1 = 977.78*somme1/H.H_f 

    # finally the age of the universe is the sum 
    H.age_univ = age0+age1  

    print()
    print( '> Current values of parameters: ')
    print( '> ----------------------------- ')
    print( '> aexp                        : ',H.aexp)
    print( '> Redshift                    : ',H.af/H.aexp-1.)
    print( '> Age of the Universe (Gyr)   : ',H.age_univ)

    H.Hub_pt   = H.H_f * np.sqrt(H.omega_f*(H.af/H.aexp)**3 +  H.omega_c_f*(H.af/H.aexp)**2 + H.omega_lambda_f)
    H.Lbox_pt  = H.Lboxp*(H.aexp/H.af)
    H.Lbox_pt2 = H.Lbox_pt / 2.0

    print( '> Hubble Parameter  (km/s/Mpc): ',H.Hub_pt)
    print( '> Box Length (Mpc)            : ',H.Lbox_pt)
    print()

    return

#***********************************************************************
def virial_12():
    '''
    compute the overdensity factor for virialization in a tophat collapse model at a given redshift 
    for a given cosmological model
    '''
    from num_rec import cubic
    
    # convert age of universe from Gyr back into inverse Hubble parameter units
    age       = H.age_univ/977.78*H.H_f
    # compute the overdensity needed to reach maximum expansion by half the age of the universe
    omega_maxexp = H.omega_f
    omega_maxexp = collapse_120(age/2.0,omega_maxexp,1.e-6)
    # calculate how far an object collapses to virial equilibrium
    eta = 2.0*H.omega_lambda_f/omega_maxexp*(H.af/H.aexp)**3
    if(eta == 0.0):
       reduce = 0.5
    else:
       a      = 2.0*eta
       b      = -(2.0+eta)
       reduce = cubic(a,0.0,b,1.0)
    
    H.vir_overdens = omega_maxexp/H.omega_f/reduce**3*(H.aexp/H.af)**3
    H.rho_mean     = H.mboxp/H.Lbox_pt**3
    print(f"> Virial overdensity             : {H.vir_overdens}")
    print(f"> Mean density (1e11 M_sun/Mpc3) : {H.rho_mean}")


#***********************************************************************
def collapse_120(age0,omm,acc):
    '''
    this subroutine performs a brute-force search using bracketing to find the value of the cosmic curvature 
    that gives turnaround at the specified expansion parameter. The expansion parameter is defined to
    be 1 at turnaround.

    note that for a constant cosmological constant, the age at a given
    expansion factor decreases monotonically with omegam (before turnaround;
    after, there is more than one age at a given expansion factor).

    third argument is the desired fractional accurracy.
    '''    
    # IMPORTANT NOTE: omax0 corresponds to a perturbation turning around at z=140 in a LCDM standard cosmology
    #                 and needs to be increased if you analyze outputs before this redshift ...
    from num_rec import age_temps_turn_around, qromo
    from functools import partial
    from scipy import integrate

    omax0=1.e7 ;omin0=1.e0; age_f=0
    age  = -1.0  # impossible value
    omax = omax0
    omin = omin0
    oml  = H.omega_lambda_f
    while( (abs(age-age0) > acc*age0)and((omax-omin) > acc*omm) ):
        omm = 0.5*(omax+omin)
        
        ### Replica of fortran subroutine `qromo`
        ### (More consistent with previous implementation)
        # age_f = qromo(age_temps_turn_around,101,10001,omm=omm,oml=oml)
        # age = qromo(age_temps_turn_around,1,101,omm=omm,oml=oml)
        
        ### Using scipy.integrate.quad
        ### (More standard way of integration in Python)
        func = partial(age_temps_turn_around, oml=oml, omm=omm)
        age_f = integrate.quad(func, 101, 10001)[0]
        age = integrate.quad(func, 1, 101)[0]
        
        if(age+age_f > age0):
            omin = omm
        else:
            omax = omm

    if( (omax==omax0)or(omin==omin0) ):
        print('WARNING: presumed bounds for omega are inadequate in collapse.')
        print('WARNING: omax,omax0,omin,omin0=',omax,omax0,omin,omin0)
    return omm

#***********************************************************************
def change_units_15():
    '''
    subroutine which goes from positions in code units to physical (non comoving)
    Mpc and from velocities in code units to peculiar (no Hubble flow) velocities 
    in km/s. masses are also changed from code units to 10^11 M_sun
    '''
    mem['pos_10']   *= H.Lbox_pt
    if(type == 'SN'): mem['vel_10'] *= H.Hub_pt*H.Lbox_pt
    H.massp = H.massp * H.mboxp
    if (H.allocated('mass_10')):
        mem['mass_10']  *= H.mboxp
    if H.allocated('m0_10'):
        mem['m0_10'] *= H.mboxp

#***********************************************************************
def det_halo_energies_1b3(h:np.void, full_PE=1000, member=None):
    from num_rec import rf
    count, _, ipos, ivel, imass, _ = member

    inp          = count#mem['nb_of_parts_o0_1'][h['id']]
    count_pairs = (inp < full_PE)
    # get potential energy 
    if (not count_pairs):
        # formula B1 of appendix of paper GalICS 1 :
        # EP = -3/5 G M^2 Rf, 
        # with Rf = 0.5 * Int_0^infty dt/sqrt((t+x)(t+y)(t+z)) = 0.5 * rf_numrec
        # NB : rf_numrec returns RF in inverse Mpc (because x is in Mpc^2)
        h['ep'] = (-3/5) * H.gravconst * h['m']**2 * rf(h['sha']**2,h['shb']**2,h['shc']**2)
        #     (km/s)^2 Mpc/(10^11 Msol) * (10^11 Msol)^2  / Mpc
        #   = (km/s)^2 * 10^11 Msol * Mpc / Mpc
        #   = (km/s)^2 * 10^11 Msol
        #   = 10^11 Msol * (km/s)^2 ---> in the same units as kinetic energy (see below)
    else:
        # Direct summation of the potential energy over all pairs of particles in the halo
        comb = list(combinations(np.arange(inp), 2)) # only count pairs once
        drs = np.squeeze(np.diff(ipos[comb], axis=1))
        dists = np.sqrt(np.sum(drs**2, axis=1))/H.Lbox_pt
        ms = np.prod(imass[comb], axis=1)
        ped = np.sum(-ms / dists)
        h['ep'] = ped * H.gravconst / H.Lbox_pt

    # get kinetic energy (in center-of-halo frame)
    drxs = correct_for_periodicity_1d(ipos[:,0] - h['px'])
    drys = correct_for_periodicity_1d(ipos[:,1] - h['py'])
    drzs = correct_for_periodicity_1d(ipos[:,2] - h['pz'])
    # add Hubble flow 
    vtsx = ivel[:,0] - h['vx'] + drxs * H.Hub_pt
    vtsy = ivel[:,1] - h['vy'] + drys * H.Hub_pt
    vtsz = ivel[:,2] - h['vz'] + drzs * H.Hub_pt
    v2s = vtsx**2 + vtsy**2 + vtsz**2
    ked = np.sum(imass * v2s)

    h['ek'] = 0.5*ked
    # get total energy 
    h['et'] = h['ek'] + h['ep']
    
    return

#***********************************************************************
def det_halo_center_multiscale_1b0(h:np.void,x0,y0,z0,r0,nx, member=None):
    '''
    The initial guesss of the C.O.M. was calculated by det_center().
    '''
    pos_10 = mem['pos_10']    
    count, indexps, ipos, _, imass, idensity = member

    mass_grid = np.zeros((nx,nx,nx), dtype=np.float64)

    xmin=x0-r0
    ymin=y0-r0
    zmin=z0-r0
    deltax=2*r0/nx

    # Assign mass to a uniform mesh with NGP
    iis = ((ipos[:,0] - xmin)/deltax).astype(int)
    jjs = ((ipos[:,1] - ymin)/deltax).astype(int)
    kks = ((ipos[:,2] - zmin)/deltax).astype(int)
    valid = (iis >= 0) & (iis < nx) & (jjs >= 0) & (jjs < nx) & (kks >= 0) & (kks < nx)
    iis = iis[valid]
    jjs = jjs[valid]
    kks = kks[valid]
    iids = indexps[valid]
    imass = imass[valid]
    mass_grid = np.zeros((nx,nx,nx), dtype=np.float64)
    np.add.at(mass_grid, (iis, jjs, kks), imass)

    # Search for the cell containing the maximum mass
    mass_max = np.max(mass_grid)
    imax,jmax,kmax = np.where(mass_grid == mass_max)
    imax = imax[0]; jmax = jmax[0]; kmax = kmax[0]

    nxnew=3
    if(deltax > nxnew*H.dcell_min):
        xc=xmin+(imax+0.5)*deltax
        yc=ymin+(jmax+0.5)*deltax
        zc=zmin+(kmax+0.5)*deltax
        det_halo_center_multiscale_1b0(h,xc,yc,zc,deltax,nxnew, member=member)
    else:
        # Find the particle with the maximum density within 
        # the cell with the maximum mass
        mass_max=0.0
        argmax = np.argmax(idensity)
        itarget = iids[argmax]

        # Assign the new halo center
        h['px']=pos_10[itarget,0]
        h['py']=pos_10[itarget,1]
        h['pz']=pos_10[itarget,2]

#***********************************************************************
def det_halo_center_sphere_1b1(h:np.void,x0,y0,z0,r0, member=None):
    count, indexps, ipos, _, imass, _ = member

    r02=r0*r0

    # compute cdm
    mtot=0
    drxs = correct_for_periodicity_1d(ipos[:,0] - x0)
    drys = correct_for_periodicity_1d(ipos[:,1] - y0)
    drzs = correct_for_periodicity_1d(ipos[:,2] - z0)
    dr2 = drxs**2 + drys**2 + drzs**2
    rmask = dr2 <= r02
    drxs = drxs[rmask]; drys = drys[rmask]; drzs = drzs[rmask]; imass = imass[rmask]
    pcx = np.sum(imass*drxs)
    pcy = np.sum(imass*drys)
    pcz = np.sum(imass*drzs)
    mtot = np.sum(imass)

    if(mtot > 0):
        xc  = pcx / mtot + x0
        yc  = pcy / mtot + y0
        zc  = pcz / mtot + z0
    else:
        xc = x0
        yc = y0
        zc = z0   

    if(r0 > H.dcell_min and mtot>0):
        det_halo_center_sphere_1b1(h,xc,yc,zc,(1-H.eps_SC)*r0, member=member)
    else:
        # search particule closest to the cdm
        drxs = correct_for_periodicity_1d(ipos[:,0] - xc)
        drys = correct_for_periodicity_1d(ipos[:,1] - yc)
        drzs = correct_for_periodicity_1d(ipos[:,2] - zc)
        dr2 = drxs**2 + drys**2 + drzs**2
        argmin = np.argmin(dr2)
        # itarget = indexps[argmin]
        distmin = dr2[argmin]
        if(distmin > 1.0): raise ValueError(f'distmin(={distmin}) > 1')

        # Assign the new halo center
        h['px']=ipos[argmin,0]
        h['py']=ipos[argmin,1]
        h['pz']=ipos[argmin,2]

def _radius_from_center(ipos, x0, y0, z0):
    if(len(ipos) == 0):
        return np.float64(np.nan)
    drxs = correct_for_periodicity_1d(ipos[:,0] - x0)
    drys = correct_for_periodicity_1d(ipos[:,1] - y0)
    drzs = correct_for_periodicity_1d(ipos[:,2] - z0)
    return np.sqrt(np.max(drxs**2 + drys**2 + drzs**2))

def _member_subset(member, mask):
    count, indexps, ipos, ivel, imass, idensity = member
    return (np.int32(np.sum(mask)), indexps[mask], ipos[mask], ivel[mask], imass[mask], idensity[mask])

def refine_centers_1b(h:np.void, member=None):
    count, indexps, ipos, ivel, imass, idensity = member
    if(H.method != "FOF"):
        r0 = _radius_from_center(ipos, h['px'], h['py'], h['pz'])
        if(np.isfinite(r0) and r0 > 0):
            if(H.DPMMC):
                nx=9
                det_halo_center_multiscale_1b0(h,h['px'],h['py'],h['pz'],r0,nx, member=member)
                r0 = _radius_from_center(ipos, h['px'], h['py'], h['pz'])
            if(H.SC):
                det_halo_center_sphere_1b1(h,h['px'],h['py'],h['pz'],r0, member=member)

    starmask = indexps >= H.ndm
    h['px*'] = np.float64(np.nan)
    h['py*'] = np.float64(np.nan)
    h['pz*'] = np.float64(np.nan)
    if(not starmask.any()):
        return

    star_member = _member_subset(member, starmask)
    _, _, spos, _, _, sdensity = star_member
    imax = np.argmax(sdensity)
    total_center = (np.float64(h['px']), np.float64(h['py']), np.float64(h['pz']))
    h['px'] = spos[imax,0]
    h['py'] = spos[imax,1]
    h['pz'] = spos[imax,2]
    if(H.method != "FOF"):
        r0 = _radius_from_center(spos, h['px'], h['py'], h['pz'])
        if(np.isfinite(r0) and r0 > 0):
            if(H.DPMMC):
                nx=9
                det_halo_center_multiscale_1b0(h,h['px'],h['py'],h['pz'],r0,nx, member=star_member)
                r0 = _radius_from_center(spos, h['px'], h['py'], h['pz'])
            if(H.SC):
                det_halo_center_sphere_1b1(h,h['px'],h['py'],h['pz'],r0, member=star_member)
    h['px*'] = h['px']
    h['py*'] = h['py']
    h['pz*'] = h['pz']
    h['px'], h['py'], h['pz'] = total_center

#***********************************************************************
def compute_spin_parameter_1c(h:np.void):
    hl                  = np.sqrt(h['Lx']**2 + h['Ly']**2 + h['Lz']**2)
    spin                = hl * np.sqrt(np.abs(h['et'])) / h['m']**2.5 / H.gravconst
    h['spin']              = spin

def _stellar_kinematics_1d(mass, dvx, dvy, dvz, vphi, vrad, vzyl, select):
    select = select & np.isfinite(mass) & (mass > 0)
    if np.count_nonzero(select) < 2:
        return np.nan, np.nan, np.nan

    weight = mass[select]
    weight_sum = np.sum(weight)
    if not np.isfinite(weight_sum) or weight_sum <= 0:
        return np.nan, np.nan, np.nan

    sig3d2 = np.average(
        dvx[select]**2 + dvy[select]**2 + dvz[select]**2,
        weights=weight,
    )

    cyl_select = select & np.isfinite(vphi) & np.isfinite(vrad) & np.isfinite(vzyl)
    if np.count_nonzero(cyl_select) < 2:
        return np.nan, np.sqrt(max(sig3d2, 0.0)), np.nan

    cyl_weight = mass[cyl_select]
    vrot = np.average(vphi[cyl_select], weights=cyl_weight)
    mean_vrad = np.average(vrad[cyl_select], weights=cyl_weight)
    mean_vzyl = np.average(vzyl[cyl_select], weights=cyl_weight)
    sigcyl2 = np.average(
        (vrad[cyl_select] - mean_vrad)**2
        + (vphi[cyl_select] - vrot)**2
        + (vzyl[cyl_select] - mean_vzyl)**2,
        weights=cyl_weight,
    )
    return vrot, np.sqrt(max(sig3d2, 0.0)), np.sqrt(max(sigcyl2, 0.0))

def compute_stellar_1d(h:np.void, member=None):
    count, indexps, mypos, myvel, mymass, mydensity = member
    h['r50'] = np.float64(np.nan)
    h['r90'] = np.float64(np.nan)
    for field in (
        'vrot', 'sig3d', 'sigcyl',
        'vrot_r50', 'sig3d_r50', 'sigcyl_r50',
        'vrot_r90', 'sig3d_r90', 'sigcyl_r90',
    ):
        h[field] = np.float64(np.nan)
    starmask = indexps >= H.ndm
    if starmask.any() and np.isfinite(h['px*']):
        star_indices = indexps[starmask] - H.ndm
        star_masses = mymass[starmask]
        star_vel = myvel[starmask]
        star_ages = mem['age_10'][star_indices]
        star_metals = mem['metal_10'][star_indices]
        drxs = correct_for_periodicity_1d(mypos[starmask,0] - h['px*'])
        drys = correct_for_periodicity_1d(mypos[starmask,1] - h['py*'])
        drzs = correct_for_periodicity_1d(mypos[starmask,2] - h['pz*'])
        dr2 = drxs**2 + drys**2 + drzs**2
        sorted_indices = np.argsort(dr2)
        sorted_masses = star_masses[sorted_indices]
        cumulative_mass = np.cumsum(sorted_masses)
        total_stellar_mass = cumulative_mass[-1]
        r50_index = np.searchsorted(cumulative_mass, 0.5 * total_stellar_mass)
        r90_index = np.searchsorted(cumulative_mass, 0.9 * total_stellar_mass)
        h['r50'] = np.sqrt(dr2[sorted_indices[r50_index]])
        h['r90'] = np.sqrt(dr2[sorted_indices[r90_index]])

        h['age'] = np.average(star_ages, weights=star_masses)
        h['metal'] = np.average(star_metals, weights=star_masses)
        inside_r50 = dr2 < h['r50']**2
        inside_r90 = dr2 < h['r90']**2

        young100 = star_ages < 0.1
        h['SFR'] = np.sum(star_masses[young100]) * 1.0e3
        h['SFR_r50'] = np.sum(star_masses[young100 & inside_r50]) * 1.0e3
        h['SFR_r90'] = np.sum(star_masses[young100 & inside_r90]) * 1.0e3

        young10 = star_ages < 0.01
        h['SFR10'] = np.sum(star_masses[young10]) * 1.0e4
        h['SFR10_r50'] = np.sum(star_masses[young10 & inside_r50]) * 1.0e4
        h['SFR10_r90'] = np.sum(star_masses[young10 & inside_r90]) * 1.0e4

        stellar_mass = np.sum(star_masses)
        angular_momentum = np.array(
            [h['Lx*'], h['Ly*'], h['Lz*']], dtype=np.float64
        )
        angular_momentum_norm = np.linalg.norm(angular_momentum)
        if (
            len(star_masses) >= 2
            and np.isfinite(stellar_mass)
            and stellar_mass > 0
        ):
            bulk_velocity = np.average(star_vel, axis=0, weights=star_masses)
            dvx = star_vel[:,0] - bulk_velocity[0]
            dvy = star_vel[:,1] - bulk_velocity[1]
            dvz = star_vel[:,2] - bulk_velocity[2]

            vphi = np.full(len(star_masses), np.nan, dtype=np.float64)
            vrad = np.full(len(star_masses), np.nan, dtype=np.float64)
            vzyl = np.full(len(star_masses), np.nan, dtype=np.float64)
            if np.isfinite(angular_momentum_norm) and angular_momentum_norm > 0:
                axis = angular_momentum / angular_momentum_norm
                zpos = drxs*axis[0] + drys*axis[1] + drzs*axis[2]
                rperp_x = drxs - zpos*axis[0]
                rperp_y = drys - zpos*axis[1]
                rperp_z = drzs - zpos*axis[2]
                rperp = np.sqrt(rperp_x**2 + rperp_y**2 + rperp_z**2)
                valid_cyl = rperp > np.finfo(np.float64).eps

                inv_rperp = np.zeros_like(rperp)
                inv_rperp[valid_cyl] = 1.0 / rperp[valid_cyl]
                erad_x = rperp_x * inv_rperp
                erad_y = rperp_y * inv_rperp
                erad_z = rperp_z * inv_rperp
                ephi_x = axis[1]*erad_z - axis[2]*erad_y
                ephi_y = axis[2]*erad_x - axis[0]*erad_z
                ephi_z = axis[0]*erad_y - axis[1]*erad_x

                vrad[valid_cyl] = (
                    dvx[valid_cyl]*erad_x[valid_cyl]
                    + dvy[valid_cyl]*erad_y[valid_cyl]
                    + dvz[valid_cyl]*erad_z[valid_cyl]
                )
                vphi[valid_cyl] = (
                    dvx[valid_cyl]*ephi_x[valid_cyl]
                    + dvy[valid_cyl]*ephi_y[valid_cyl]
                    + dvz[valid_cyl]*ephi_z[valid_cyl]
                )
                vzyl[valid_cyl] = (
                    dvx[valid_cyl]*axis[0]
                    + dvy[valid_cyl]*axis[1]
                    + dvz[valid_cyl]*axis[2]
                )

            all_stars = np.ones(len(star_masses), dtype=bool)
            for suffix, select in (
                ('', all_stars),
                ('_r50', inside_r50),
                ('_r90', inside_r90),
            ):
                vrot, sig3d, sigcyl = _stellar_kinematics_1d(
                    star_masses, dvx, dvy, dvz, vphi, vrad, vzyl, select
                )
                h[f'vrot{suffix}'] = vrot
                h[f'sig3d{suffix}'] = sig3d
                h[f'sigcyl{suffix}'] = sigcyl

def compute_ssp_photometry_1f(ih1, h:np.void, member=None):
    count, indexps, mypos, myvel, mymass, mydensity = member
    for model in MODELS:
        mem[f'photometry_{model.lower()}'][ih1]['id'] = h['id']

    starmask = indexps >= H.ndm
    if np.count_nonzero(starmask) == 0 or not np.isfinite(h['px*']):
        return

    star_indices = indexps[starmask] - H.ndm
    age = mem['age_10'][star_indices]
    metal = mem['metal_10'][star_indices]
    if H.allocated('m0_10'):
        mass_msol = mem['m0_10'][star_indices] * 1.0e11
    else:
        mass_msol = mymass[starmask] * 1.0e11
    valid = (
        np.isfinite(age) & (age > 0)
        & np.isfinite(metal) & (metal > 0)
        & np.isfinite(mass_msol) & (mass_msol > 0)
    )
    if np.count_nonzero(valid) == 0:
        return

    age = age[valid]
    metal = metal[valid]
    mass_msol = mass_msol[valid]
    star_pos = mypos[starmask][valid]
    drx = correct_for_periodicity_1d(star_pos[:,0] - h['px*'])
    dry = correct_for_periodicity_1d(star_pos[:,1] - h['py*'])
    drz = correct_for_periodicity_1d(star_pos[:,2] - h['pz*'])
    radius = np.sqrt(drx**2 + dry**2 + drz**2)
    order = np.argsort(radius)

    for model in MODELS:
        result = mem[f'photometry_{model.lower()}'][ih1]
        coordinates = interpolation_coordinates(model, age, metal)
        r_luminosity = None
        for output_name, _, _ in BANDS:
            magnitude_per_msol = interpolate_magnitude(
                model, output_name, coordinates
            )
            luminosity = mass_msol * 10.0**(-0.4 * magnitude_per_msol)
            total_luminosity = np.sum(luminosity)
            if np.isfinite(total_luminosity) and total_luminosity > 0:
                result[output_name] = -2.5 * np.log10(total_luminosity)
            if output_name == 'rmag':
                r_luminosity = luminosity

        total_r_luminosity = np.sum(r_luminosity)
        if not np.isfinite(total_r_luminosity) or total_r_luminosity <= 0:
            continue
        result['age_r'] = np.average(age, weights=r_luminosity)
        result['metal_r'] = np.average(metal, weights=r_luminosity)
        cumulative_luminosity = np.cumsum(r_luminosity[order])
        i50 = np.searchsorted(
            cumulative_luminosity, 0.5 * total_r_luminosity
        )
        i90 = np.searchsorted(
            cumulative_luminosity, 0.9 * total_r_luminosity
        )
        result['r50_r'] = radius[order[min(i50, len(order)-1)]]
        result['r90_r'] = radius[order[min(i90, len(order)-1)]]

def compute_extended_profiles_1e(h:np.void, member=None):
    count, indexps, ipos, ivel, imass, idensity = member
    h['vmaxcir'] = np.float64(np.nan)
    h['rmaxcir'] = np.float64(np.nan)
    h['cNFWerr'] = np.float64(np.nan)
    h['inslope'] = np.float64(np.nan)
    h['inslopeerr'] = np.float64(np.nan)

    drxs = correct_for_periodicity_1d(ipos[:,0] - h['px'])
    drys = correct_for_periodicity_1d(ipos[:,1] - h['py'])
    drzs = correct_for_periodicity_1d(ipos[:,2] - h['pz'])
    dist = np.sqrt(drxs**2 + drys**2 + drzs**2)

    nonzero = dist > 0
    if(np.sum(nonzero) > 0):
        order = np.argsort(dist[nonzero])
        sdist = dist[nonzero][order]
        smass = imass[nonzero][order]
        vcir = np.sqrt(H.gravconst * np.cumsum(smass) / sdist)
        imax = np.argmax(vcir)
        h['vmaxcir'] = vcir[imax]
        h['rmaxcir'] = sdist[imax]

    dmmask = (indexps < H.ndm) & nonzero
    ndm = np.sum(dmmask)
    if(ndm < 8 or h['rvir'] <= 0):
        h['cNFW'] = np.float64(np.nan)
        return

    dmdist = dist[dmmask]
    dmmass = imass[dmmask]
    order = np.argsort(dmdist)
    dmdist = dmdist[order]
    dmmass = dmmass[order]
    rvir = h['rvir']
    inside = dmdist <= rvir
    if(np.sum(inside) < 8):
        h['cNFW'] = np.float64(np.nan)
        return
    dmdist = dmdist[inside]
    dmmass = dmmass[inside]

    try:
        from scipy.optimize import curve_fit

        def log_nfw(r, rs, logrho0):
            x = r / rs
            return logrho0 - np.log10(x) - 2.0*np.log10(1.0 + x)

        nbins = max(8, min(80, int(np.sqrt(len(dmdist)))))
        bins = np.logspace(np.log10(dmdist[0]), np.log10(dmdist[-1]), nbins + 1)
        shell_mass, _ = np.histogram(dmdist, bins=bins, weights=dmmass)
        shell_num, _ = np.histogram(dmdist, bins=bins)
        shell_vol = 4.0/3.0*H.pi*(bins[1:]**3 - bins[:-1]**3)
        radius = 0.5*(bins[1:] + bins[:-1])
        valid = (shell_num > 0) & (shell_mass > 0) & (shell_vol > 0)
        if(np.sum(valid) >= 4):
            rho = shell_mass[valid] / shell_vol[valid]
            radius = radius[valid]
            sigma = 1.0 / np.sqrt(shell_num[valid])
            popt, pcov = curve_fit(
                log_nfw, radius, np.log10(rho),
                p0=[max(rvir/10.0, radius[0]), np.log10(rho[0])],
                sigma=sigma, absolute_sigma=False,
                bounds=([0.0, -np.inf], [rvir, np.inf]),
            )
            rs = popt[0]
            rs_err = np.sqrt(pcov[0,0]) if pcov.size else np.nan
            h['cNFW'] = rvir / rs
            h['cNFWerr'] = h['cNFW'] * np.sqrt((rs_err/rs)**2) if rs > 0 else np.nan
        else:
            h['cNFW'] = np.float64(np.nan)
    except Exception:
        h['cNFW'] = np.float64(np.nan)
        h['cNFWerr'] = np.float64(np.nan)

    try:
        rs = rvir / h['cNFW'] if np.isfinite(h['cNFW']) and h['cNFW'] > 0 else rvir
        rcut = min(rs, rvir)
        mask = dmdist < rcut
        if(np.sum(mask) >= 4):
            r_in = dmdist[mask]
            m_in = dmmass[mask]
            xx = np.log10(r_in / rvir)
            shell_vol = 4.0/3.0*H.pi*r_in**3
            shell_vol = np.diff(np.insert(shell_vol, 0, 0.0))
            # Drop roundoff-degenerate shells: when adjacent sorted radii differ
            # by ~machine epsilon, the per-particle shell volume collapses to a
            # near-zero value and mass/shell_vol spikes, corrupting the log-log
            # fit and making `inslope` platform-sensitive. Filtering these shells
            # removes the spike; halos with no degenerate shells are unaffected
            # (this is a no-op there), so the fix is surgical.
            rel_width = np.empty_like(r_in)
            rel_width[0] = 1.0  # first shell (0 -> r[0]) is never degenerate
            rel_width[1:] = (r_in[1:] - r_in[:-1]) / np.maximum(r_in[1:], np.finfo(np.float64).tiny)
            valid = (shell_vol > 0) & (m_in > 0) & (rel_width > 1.0e-9)
            xx = xx[valid]
            yy = np.log10(m_in[valid] / shell_vol[valid])
            if(len(xx) >= 4):
                from scipy.optimize import curve_fit
                def linefit(x, a, b):
                    return a*x + b
                popt, pcov = curve_fit(linefit, xx, yy)
                h['inslope'] = popt[0]
                h['inslopeerr'] = np.sqrt(pcov[0,0]) if pcov.size else np.nan
    except Exception:
        h['inslope'] = np.float64(np.nan)
        h['inslopeerr'] = np.float64(np.nan)




#***********************************************************************
def det_inertial_tensor_1b20(h:np.void, member=None):
    '''
    Compute inertial tensor with respect to center of halo (either cdm or mbp)
    '''
    _, _, ipos, _, imass, _ = member

    drxs = correct_for_periodicity_1d(ipos[:,0] - h['px']) # [Mpc]
    drys = correct_for_periodicity_1d(ipos[:,1] - h['py'])
    drzs = correct_for_periodicity_1d(ipos[:,2] - h['pz'])
    mat = np.zeros((3,3), dtype=np.float64)
    mat[0,0] = np.sum(imass*drxs*drxs) # (10^11 M_sun) * (Mpc)^2
    mat[0,1] = np.sum(imass*drxs*drys)
    mat[0,2] = np.sum(imass*drxs*drzs)
    mat[1,0] = np.sum(imass*drxs*drys)
    mat[1,1] = np.sum(imass*drys*drys)
    mat[1,2] = np.sum(imass*drys*drzs)
    mat[2,0] = np.sum(imass*drxs*drzs)
    mat[2,1] = np.sum(imass*drys*drzs)
    mat[2,2] = np.sum(imass*drzs*drzs)
    return mat

#***********************************************************************
def det_main_axis_1b2(h:np.void, member=None):
    '''
    determine the principal axis of the halo (h['sha'],b,c)
    '''
    from scipy.linalg import eigh

    mat = det_inertial_tensor_1b20(h, member=member) # (10^11 M_sun) * (Mpc)^2

    d, v = eigh(mat) # d in (10^11 M_sun) * (Mpc)^2, v is unitless
    d      = np.sqrt(d/h['m']) # in Mpc
    h['sha'] = d[0] # in Mpc
    h['shb'] = d[1]
    h['shc'] = d[2]
    return d,v

#***********************************************************************
def det_ang_momentum_per_shell_1b5(h,amax,bmax,cmax,v,fagor:FortranFile=None, member=None):    
    # NOT YET OPTIMIZED!!
    # ALSO, THERE COULD BE A RACE CONDITION FOR WRITING
    dr=vector();dv=vector()
    L = [vector() for _ in range(H.nshells)]
    m = np.zeros(H.nshells, dtype=np.float64)

    for i0 in range(H.nshells):
       L[i0].x = 0.0
       L[i0].y = 0.0
       L[i0].z = 0.0
    m[:]   = 0.0
    indexp = mem['first_part_oo_1'][h['id']]
    while(indexp != -1):
        # particle positions relative to center of halo
        dr.x   = mem['pos_10'][indexp,0] - h['px']
        dr.y   = mem['pos_10'][indexp,1] - h['py']
        dr.z   = mem['pos_10'][indexp,2] - h['pz']
        correct_for_periodicity(dr)
        # convert dr into ellipsoid coords
        dra    = dr.x*v[0,0]+dr.y*v[1,0]+dr.z*v[2,0]
        drb    = dr.x*v[0,1]+dr.y*v[1,1]+dr.z*v[2,1]
        drc    = dr.x*v[0,2]+dr.y*v[1,2]+dr.z*v[2,2]
        # index of shell containing particle
        r_ell  = np.sqrt((dra / amax)**2 + (drb / bmax)**2 + (drc / cmax)**2)
        i_ell  = int(r_ell*H.nshells)
        if (i_ell > H.nshells):
            print('> Problem in get_ang_momentum_per_shell : i_ell > nshells ')
            raise ValueError(f"{i_ell},{H.nshells}")
        # velocity relative to halo velocity (cdm)
        dv.x   = mem['vel_10'][indexp,0]-h['vx']
        dv.y   = mem['vel_10'][indexp,1]-h['vy']
        dv.z   = mem['vel_10'][indexp,2]-h['vz']
        # update mass and angular momentum of shell
        if(H.BIG_RUN):
            m[i_ell]   += H.massp
            L[i_ell].x += H.massp*(dr.y*dv.z - dr.z*dv.y) # in 10**11 Msun * km/s * Mpc
            L[i_ell].y += H.massp*(dr.z*dv.x - dr.x*dv.z)
            L[i_ell].z += H.massp*(dr.x*dv.y - dr.y*dv.x)
        else:
            m[i_ell]   += mem['mass_10'][indexp]
            L[i_ell].x += mem['mass_10'][indexp]*(dr.y*dv.z - dr.z*dv.y)
            L[i_ell].y += mem['mass_10'][indexp]*(dr.z*dv.x - dr.x*dv.z)
            L[i_ell].z += mem['mass_10'][indexp]*(dr.x*dv.y - dr.y*dv.x)
        indexp = mem['linked_list_oo_1'][indexp]   
    
    fagor.write_record(
        h['id'], h['px'], h['py'], h['pz'], h['vx'], h['vy'], h['vz'],
        h['rvir'], h['mvir'], h['m'], h['r'], h['spin'],
        amax, bmax, cmax, v[0,:3], v[1,:3], v[2,:3],
        m[0:H.nshells], L[0:H.nshells].x, L[0:H.nshells].y, L[0:H.nshells].z
        )
    
    return
