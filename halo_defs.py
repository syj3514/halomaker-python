import numpy as np
from multiprocessing import Pool, shared_memory
from collections import defaultdict
from scipy.io import FortranFile
from tqdm import tqdm
import time
import atexit, signal
import sys, os


#======================================================================
# Make Flags
#======================================================================
RENORM = False
BIG_RUN = False
METALS = True
TQDM = True
DEV = False
ANG_MOM_OF_R = False
Test_FOF = False
SCIPY = True
FORTRAN = True

#======================================================================
# Memory allocation
#======================================================================
mprefix = ["HaloMaker", "uNone", "tNone"] # Halo/Galaxy, uid, runtime
def collapse_mprefix():
    global mprefix
    return "_".join([i for i in mprefix])
mem = {}
mem_address = {}
globvars = {}
massalloc=True

pos_ref_0 = np.array([
    -1., -1., -1.,
     1., -1., -1.,
    -1.,  1., -1.,
     1.,  1., -1.,
    -1., -1.,  1.,
     1., -1.,  1.,
    -1.,  1.,  1.,
     1.,  1.,  1.], dtype=np.float64)
pos_ref_0 = pos_ref_0.reshape((3, 8), order='F')

def maccess(name):
    global mem_address
    temp, shape, dtype = mem_address[name]
    return np.ndarray(shape, dtype=dtype, buffer=temp.buf)


def mlist():
    global mem_address
    rows = []
    for name in mem_address.keys():
        shape = mem_address[name][1]
        dtype = mem_address[name][2]
        dtype_str = getattr(dtype, "name", str(dtype))
        filesize_MB = mem_address[name][0].size / (1024*1024)
        rows.append((name, filesize_MB, str(shape).replace(' ', ''), dtype_str))
    if len(rows)==0: return

    name_w  = max(len(r[0]) for r in rows)
    size_w  = max(len(f"{r[1]:.2f}") for r in rows)+3
    shape_w = max(len(r[2]) for r in rows)
    dtype_w = max(len(r[3]) for r in rows)

    header = (
        f"{'name':<{name_w}}  "
        f"{'MB':>{size_w}}  "
        f"{'shape':<{shape_w}}  "
        f"{'dtype':<{dtype_w}}"
    )
    print("-" * len(header))
    print("Memory Allocations:")
    print(header)
    print("-" * len(header))
    size_total = 0
    for name, size_mb, shape_str, dtype_str in rows:
        print(
            f"{name:<{name_w}}  "
            f"{size_mb:>{size_w}.2f}  "
            f"{shape_str:<{shape_w}}  "
            f"{dtype_str:<{dtype_w}}"
        )
        size_total += size_mb
    ustr = 'MB'
    if size_total > 1000:
        size_total /= 1024
        print(f" (Total size in GB)")
        ustr = 'GB'
    print()
    print(
        f"{'TOTAL':<{name_w}}  "
        f"{size_total:>{size_w}.2f} {ustr}  "
    )
    print("-" * len(header))


def allocate(name, shape, dtype:type|np.dtype=np.float64):
    global mem, mem_address
    if(DEV): print(f"\t@Allocating memory `{name}` | {shape}")
    arr = np.empty(shape, dtype=dtype)
    _name = collapse_mprefix() + f"_{name}"
    temp = shared_memory.SharedMemory(create=True, size=arr.nbytes, name=_name)
    if(name in mem_address.keys()):
        if(mem_address[name] is not None):
            raise Exception(f"\t@Memory address `{name}` already exists")
    mem_address[name] = [temp, arr.shape, arr.dtype]
    mem[name] = np.ndarray(arr.shape, dtype=arr.dtype, buffer=mem_address[name][0].buf)
    arr = None
    return mem[name]

def allocated(name):
    global mem
    return name in mem.keys()

def deallocate(*names):
    global mem, mem_address
    for name in names:
        if(name in mem.keys()):
            if(DEV): print(f"\t@Deallocating memory `{name}`")
            if mem_address[name] is not None:
                mem[name] = None
                mem_address[name][0].close()
                mem_address[name][0].unlink()
                mem_address[name][0] = None
            del mem[name]
            del mem_address[name]
        else:
            print(f"\t@Memory address `{name}` does not exist")
    if(list(mem.keys())==0):
        signal.signal(signal.SIGINT, signal.SIG_DFL)
        signal.signal(signal.SIGPIPE, signal.SIG_DFL)
        signal.signal(signal.SIGTERM, signal.SIG_DFL)

def flush(msg=True, parent=''):
    signal.signal(signal.SIGINT, signal.SIG_DFL)
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    signal.signal(signal.SIGTERM, signal.SIG_DFL)
    if(len(mem_address) > 0):
        if(msg or DEV): print(f"\t@{parent} Clearing memory")
        if(msg or DEV): print(f"\t   {[i[0].name for i in mem_address.values()]}")
        deallocate(*mem_address.keys())
    if(msg or DEV): print("\t@Memory Clear Done")
    sys.exit(0)
atexit.unregister(flush)

def adjust_size(name, new_shape, dtype, tmp=None, init=None):
    if isinstance(new_shape, int):
        new_shape = (new_shape,)
    old_shape = mem[name].shape
    if tmp is None:
        tmp = np.empty(old_shape, dtype=dtype)
    if len(old_shape)==1:
        tmp[:] = mem[name][:]
        deallocate(name)
        allocate(name, new_shape, dtype=dtype)
        leng = min(old_shape[0], new_shape[0])
        mem[name][:leng] = tmp[:leng]
        if init is not None:
            mem[name][leng:] = init
    elif len(old_shape)==2:
        nrow = old_shape[0]
        tmp[:,:] = mem[name][:,:]
        deallocate(name)
        allocate(name, new_shape, dtype=dtype)
        leng = min(old_shape[1], new_shape[1])
        mem[name][:nrow,:leng] = tmp[:nrow,:leng]
        if init is not None:
            mem[name][nrow:,leng:] = init



def datdump(data, path, msg=False):
    # assert isinstance(data[0], np.ndarray), "Data should be numpy.ndarray"
    # assert isinstance(data[1], str)
    leng = len(data)
    with open(path, "wb") as f:
        f.write(leng.to_bytes(4, byteorder='little'))
        f.write(data.tobytes())
    full_path = os.path.abspath(path)
    os.chmod(full_path, fchmod); os.chown(full_path, uid, gid)
    if(msg): print(f" `{path}` saved")

def datload(path, dtype='f8', msg=False):
    a = dtype[0]; b=int(dtype[1])
    with open(path, "rb") as f:
        leng = int.from_bytes(f.read(4), byteorder='little')
        data = np.frombuffer(f.read(b*leng), dtype=dtype)
    if(msg): print(f" `{path}` loaded")
    return data

#======================================================================
# useful types :
#======================================================================
class vector:
    __slots__ = ["x", "y", "z"]
    def __init__(self):
        self.x:float = 0
        self.y:float = 0
        self.z:float = 0

class shape:
    __slots__ = ["a", "b", "c"]
    def __init__(self):
        self.a:float = 0
        self.b:float = 0
        self.c:float = 0

class baryon:
    __slots__ = ["rvir", "mvir", "tvir", "cvel"]
    def __init__(self):
        self.rvir:float = 0
        self.mvir:float = 0
        self.tvir:float = 0
        self.cvel:float = 0

class extend:
    __slots__ = ['mcontam']
    def __init__(self):
        self.mcontam:float = 0

class hprofile:
    __slots__ = ["rho_0", "r_c", "cNFW"]
    def __init__(self):
        self.rho_0:float = 0
        self.r_c:float = 0
        self.cNFW:float = 0



halo_dtype = np.dtype([
        ('id', 'i4'),
        ('timestep', 'i4'),
        ('nmem', 'i4'),('ndm', 'i4'),('nstar', 'i4'),
        ('nbsub', 'i4'),
        ('hosthalo', 'i4'),
        ('hostsub', 'i4'),
        ('level', 'i4'),
        ('nextsub', 'i4'),
        ('px', 'f8'),('py', 'f8'),('pz', 'f8'),
        ('px*', 'f8'),('py*', 'f8'),('pz*', 'f8'),
        ('vx', 'f8'),('vy', 'f8'),('vz', 'f8'),
        ('Lx', 'f8'),('Ly', 'f8'),('Lz', 'f8'),
        ('Lx*', 'f8'),('Ly*', 'f8'),('Lz*', 'f8'),
        ('sha', 'f8'),('shb', 'f8'),('shc', 'f8'),
        ('m', 'f8'),('mdm', 'f8'),('m*', 'f8'),
        ('r', 'f8'),('r*', 'f8'),('r50', 'f8'),('r90', 'f8'),
        ('age', 'f8'),('metal', 'f8'),
        ('SFR', 'f8'),('SFR_r50', 'f8'),('SFR_r90', 'f8'),
        ('SFR10', 'f8'),('SFR10_r50', 'f8'),('SFR10_r90', 'f8'),
        ('spin', 'f8'),
        ('sigma', 'f8'),('sigma_dm', 'f8'),('sigma*', 'f8'),
        ('vrot', 'f8'),('sig3d', 'f8'),('sigcyl', 'f8'),
        ('vrot_r50', 'f8'),('sig3d_r50', 'f8'),('sigcyl_r50', 'f8'),
        ('vrot_r90', 'f8'),('sig3d_r90', 'f8'),('sigcyl_r90', 'f8'),
        ('ek', 'f8'),('ep', 'f8'),('et', 'f8'),
        ('rvir','f8'),('mvir','f8'),('tvir','f8'),('cvel','f8'),
        ('rho_0','f8'),('r_c','f8'),('cNFW','f8'),('cNFWerr','f8'),
        ('vmaxcir','f8'),('rmaxcir','f8'),
        ('inslope','f8'),('inslopeerr','f8'),
        ('mcontam','f8')
])

def clear_halo(h):
    h['id']=0
    h['timestep'] = 0
    h['nmem'] = 0; h['ndm'] = 0; h['nstar'] = 0
    h['nbsub'] = 0
    h['hosthalo'] = 0
    h['hostsub'] = 0
    h['level'] = 1
    h['nextsub'] = -1
    h['px'] = np.float64(0.0); h['py'] = np.float64(0.0); h['pz'] = np.float64(0.0)
    h['px*'] = np.float64(np.nan); h['py*'] = np.float64(np.nan); h['pz*'] = np.float64(np.nan)
    h['vx'] = np.float64(0.0); h['vy'] = np.float64(0.0); h['vz'] = np.float64(0.0)
    h['Lx'] = np.float64(0.0); h['Ly'] = np.float64(0.0); h['Lz'] = np.float64(0.0)
    h['Lx*'] = np.float64(0.0); h['Ly*'] = np.float64(0.0); h['Lz*'] = np.float64(0.0)
    h['sha'] = np.float64(0.0);h['shb'] = np.float64(0.0);h['shc'] = np.float64(0.0)
    h['m'] = np.float64(0.0); h['mdm'] = np.float64(0.0); h['m*'] = np.float64(0.0)
    h['r'] = np.float64(0.0); h['r*'] = np.float64(np.nan); h['r50'] = np.float64(np.nan); h['r90'] = np.float64(np.nan)
    h['age'] = np.float64(np.nan); h['metal'] = np.float64(np.nan)
    h['SFR'] = np.float64(0.0); h['SFR_r50'] = np.float64(0.0); h['SFR_r90'] = np.float64(0.0)
    h['SFR10'] = np.float64(0.0); h['SFR10_r50'] = np.float64(0.0); h['SFR10_r90'] = np.float64(0.0)
    h['spin'] = np.float64(0.0)
    h['sigma'] = np.float64(0.0); h['sigma_dm'] = np.float64(0.0); h['sigma*'] = np.float64(0.0)
    h['vrot'] = np.float64(np.nan); h['sig3d'] = np.float64(np.nan); h['sigcyl'] = np.float64(np.nan)
    h['vrot_r50'] = np.float64(np.nan); h['sig3d_r50'] = np.float64(np.nan); h['sigcyl_r50'] = np.float64(np.nan)
    h['vrot_r90'] = np.float64(np.nan); h['sig3d_r90'] = np.float64(np.nan); h['sigcyl_r90'] = np.float64(np.nan)
    h['ek'] = np.float64(0.0); h['ep'] = np.float64(0.0); h['et'] = np.float64(0.0)
    h['rvir'] = np.float64(0.0); h['mvir'] = np.float64(0.0); h['tvir'] = np.float64(0.0); h['cvel'] = np.float64(0.0)
    h['rho_0'] = np.float64(0.0); h['r_c'] = np.float64(0.0); h['cNFW'] = np.float64(np.nan); h['cNFWerr'] = np.float64(np.nan)
    h['vmaxcir'] = np.float64(np.nan); h['rmaxcir'] = np.float64(np.nan)
    h['inslope'] = np.float64(np.nan); h['inslopeerr'] = np.float64(np.nan)
    h['mcontam'] = np.float64(0.0)
    

liste_halos_o0 = np.empty(0, dtype=halo_dtype)

photometry_ssp_dtype = np.dtype([
    ('id', 'i4'),
    ('umag', 'f8'), ('gmag', 'f8'), ('rmag', 'f8'),
    ('imag', 'f8'), ('zmag', 'f8'),
    ('Umag', 'f8'), ('Bmag', 'f8'), ('Vmag', 'f8'), ('Kmag', 'f8'),
    ('age_r', 'f8'), ('metal_r', 'f8'),
    ('r50_r', 'f8'), ('r90_r', 'f8'),
])
photometry_ssp_fields = photometry_ssp_dtype.names[1:]
#======================================================================

cosmo_table = None
cosmo_table_key = None

def H_over_H0(aexp):
    return np.sqrt(omega_f * aexp ** -3 + omega_lambda_f)

def set_cosmology(n=5000):
    global cosmo_table, cosmo_table_key
    key = (float(omega_f), float(omega_lambda_f), float(H_f), int(n))
    if cosmo_table is not None and cosmo_table_key == key:
        return cosmo_table
    try:
        from scipy.integrate import cumulative_trapezoid as cumtrapz
    except ImportError:
        from scipy.integrate import cumtrapz

    km = 1e5
    Gyr = 3.154e16
    Mpc = 3.08e24
    aarr = np.linspace(0, 1, n)[1:] ** 2
    aarr_st = (aarr[:-1] + aarr[1:]) / 2
    duda = 1. / (aarr_st ** 3 * H_over_H0(aarr_st))
    dtda = 1. / (H_f * km * Gyr / Mpc * aarr_st * H_over_H0(aarr_st))
    aarr = aarr[1:]

    uarr = cumtrapz(duda[::-1], aarr[::-1], initial=0)[::-1]
    tarr = cumtrapz(dtda, aarr, initial=0)
    cosmo_table = np.core.records.fromarrays(
        [aarr, tarr, uarr],
        dtype=[('aexp', 'f8'), ('t', 'f8'), ('u', 'f8')])
    cosmo_table_key = key
    return cosmo_table

def interpolate_cosmo_table(value, name1, name2):
    table = set_cosmology()
    return np.interp(value, table[name1], table[name2])

def epoch_to_age(epoch):
    return interpolate_cosmo_table(epoch, 'u', 't')


#======================================================================
# parameters relative to the simulation analysis
#======================================================================
gravsoft = 'cubsplin' # type of gravitational softening
nbPes = 1 # obsolete vars for reading treecode format
nsteps = 1 # number of timesteps
alpha = np.float64(1.0); tnow = np.float64(1.0)
nusedpart = 1 # number of particles in the simulation
nMembers = 1 # minimal number of particles in a fof halo
b_init = np.float64(1.0) # linking length parameter of the fof at z=0
profile = 'TSIS' # type of halo profile (only isothermal_sphere yet)
ninterp = 1 # nb of bins for interp. of smoothed grav. field
FlagPeriod = 1 # flag for periodicity of boundary conditions
fPeriod = np.array([1.0, 1.0, 1.0], dtype=np.float64)
zoomin = False # flag for zoom-in simulation
# Manual Zoombox
# [x1,x2,y1,y2,z1,z2] code unit centered on 0. 
    # ex: np.array([-0.5,0.5,-0.5,0.5,-0.5,0.5])
# Or, you can override in the command
#   `--zoombox [-0.310,-0.290,-0.150,-0.110,-0.010,0.010]`
zoombox = None # 
#----------- For gadget format: ----------------------------------------
# nhr = 1 # to read only the selected HR particles
# minlrmrat = 1.0 # to recognize contaminated halos
#======================================================================


#======================================================================
# Definitions specific to input/output
#======================================================================
output_dir = "."
file_num = "00001"
numstep = 1
# errunit = 0
write_resim_masses = False # for writing resim_masses.dat file
dchmod = 0o755
fchmod = 0o644
uid = -1
gid = -1
prefix=""

#======================================================================


#======================================================================
# Constants
#======================================================================
# gravconst = np.float64(4.302e-6) # G in units of (km/s)^2 kpc/(10^11 Msol)
gravconst = np.float64(430.1) # G in units of (km/s)^2 kpc/(10^11 Msol)
pi = np.pi
#======================================================================


#======================================================================
# Global variables 
#======================================================================
ndim = 3
#   real(kind=8),allocatable         :: pos(:,:),vel(:,:)
massp = np.float64(0.0)
#   real(kind=8),allocatable         :: epsvect(:),mass(:)
omega_t,omega_lambda_t,omega_f,omega_lambda_f,omega_c_f = np.float64(0.0), np.float64(0.0), np.float64(0.0), np.float64(0.0), np.float64(0.0)
rho_crit,aexp,Lboxp,mboxp,af,ai,Lf,H_f,H_i = np.float64(0.0), np.float64(0.0), np.float64(0.0), np.float64(0.0), np.float64(0.0), np.float64(0.0), np.float64(0.0), np.float64(0.0), np.float64(0.0)
age_univ,Lbox_pt,Lbox_pt2,Hub_pt,omega_0,hubble,omega_lambda_0 = np.float64(0.0), np.float64(0.0), np.float64(0.0), np.float64(0.0), np.float64(0.0), np.float64(0.0),np.float64(0.0)
vir_overdens,rho_mean = np.float64(0.0), np.float64(0.0)
#   integer(kind=4),allocatable      :: linked_list_oo(:), whereIam_parts(:)
#   integer(kind=4),allocatable      :: first_part_oo(:),nb_of_parts_o0(:)
#   type (halo),allocatable          :: liste_halos_o0(:)
phsmooth_oo: np.ndarray[np.float64] = np.zeros(2+ninterp)
nb_of_halos, nb_of_subhalos = 0, 0
numero_step = 1
simtype = 'Ra3'
#======================================================================


#======================================================================
# defs for Adaptahop
#======================================================================
nparmax=512*512*512
nparbuffer=128**3
ncellbuffermin=128**3
nlevelmax=30; npartpercell=100
lin=10; lsin=11; lin2=12
loudis=12; lounei=13; lounode=14; loupartnode=15; lounodedyn=16
bignum=1.e30
# Physical constants (units : m s kg) ->
gravitational_constant=np.float64(6.6726e-11)
critical_density= np.float64(1.8788e-26)
mega_parsec=np.float64(3.0857e22)
solar_mass=np.float64(1.989e30)
convert_in_mps2=np.float64(1.e6)
sizeroot = np.float64(0.0)
xlong, ylong, zlong, boxsize, boxsize2 = np.float64(0.0), np.float64(0.0), np.float64(0.0), np.float64(0.0), np.float64(0.0)
boxarr = np.array([xlong, ylong, zlong], dtype=np.float64)
xlongs2, ylongs2, zlongs2 = np.float64(0.0), np.float64(0.0), np.float64(0.0)
omega0,omegaL,mass_in_kg,GMphys = np.float64(0.0), np.float64(0.0), np.float64(0.0), np.float64(0.0)
aexp_hubble_const,aexp_mega_parsec_sqrt3 = np.float64(0.0), np.float64(0.0)
aexp_mega_parsec,aexp_max = np.float64(0.0), np.float64(0.0)
npart = 0
nstar = 0
ndm = 0
nvoisins,nhop,ntype = 0, 0, 0
ncellmx = 0
ngroups,nmembthresh,nnodes,nnodesmax = 0, 0, 0, 0
ncpu,nmpi,niterations = 0, 0, 0
ncellbuffer = 0
rho_threshold = np.float64(0.0)
verbose,megaverbose,periodic = False, False, False
fgas = np.float64(0.0)
fudge,alphap,epsilon,fudgepsilon = np.float64(0.0), np.float64(0.0), np.float64(1e-2), np.float64(0.0)
pos_shift = np.zeros(3, dtype=np.float64)

class grp:
    __slots__ = ["nhnei", "isad_gr", "rho_saddle_gr", "tmp"]
    def __init__(self):
        self.nhnei:np.int32 = 0
        self.isad_gr:np.ndarray[np.int32] = None # integer(kind=4), dimension(:),pointer
        self.rho_saddle_gr:np.ndarray[np.float64] = None # real(kind=8), dimension(:),pointer
        self.tmp:tuple = None

    def __str__(self):
        return f"grp(nhnei={self.nhnei}, isad_gr={self.isad_gr}, rho_saddle_gr={self.rho_saddle_gr}) (tmp={self.tmp})"
# class supernode:
#     __slots__ = ["level","mother","firstchild","nsisters","sister","rho_saddle","density","densmax","radius","mass","truemass","position"]
#     def __init__(self):
#         # Int32
#         self.mass = 0
#         self.level = 0
#         self.mother = 0
#         self.firstchild = 0
#         self.nsisters = 0
#         self.sister = 0
#         # Float64
#         self.rho_saddle = 0.0
#         self.density = 0.0
#         self.densmax = 0.0
#         self.radius = 0.0
#         self.truemass = 0.0
#         self.position = np.empty(3, dtype=np.float64)

# node_0:list['supernode'] = []
node_dtype = np.dtype([
                    ('rho_saddle','f8'), ('density','f8'), ('densmax','f8'), ('radius','f8'),
                    ('truemass','f8'), ('px','f8'), ('py','f8'), ('pz','f8'),
                    ('mass','i4'), ('level','i4'), ('mother','i4'), ('firstchild','i4'),
                    ('nsisters','i4'), ('sister','i4')
                    ])
node_0:np.ndarray = np.empty(0, dtype=node_dtype)
group:list['grp'] = []
#    type (grp), allocatable       :: group(:)
#    type (supernode), allocatable :: node_0(:)

#======================================================================
# Flags for halo finder selection
#======================================================================
method = 'MSM' # flag to notify which and how the halofinder is to be used
fsub = False # flag to notify whether subhaloes are included
cdm = False # flag to select particle closest to the cdm instead of the one with the highest density
DPMMC = False # flag to select the densest particle in the most massive cell of the halo (not with FOF)
SC = False # flag to select the com within concentric spheres (not with FOF)
dcell_min = np.float64(0.0)
eps_SC = np.float64(0.0)
dump_dms = False
dump_stars = False
nchem=0

#======================================================================
# array to build the structure tree
#======================================================================
nstruct = 0

# used for the merger history method
ex_nb_of_structs = 0


if(ANG_MOM_OF_R):
    agor_unit = 19
    agor_file = ""
    nshells = 100
  

def frange(first, last, step=1):
    if(step>0):
        return range(first, last+1, step)
    else:
        return range(first, last-1, step)





###################################################
# Variables translators from Fortran to Python
###################################################
# I (Seyoung Jeon) changed some variable names to be more descriptive, 
# but I will keep the original names as comments for reference. 
# Please refer to the Fortran code for the original variable names and their usage.
# (Hopefully, useful for search)

# Fortran            | Python
# liste_parts        | whereIam_parts


def fbool(value):
    v = str(value).strip().lower()
    if v in ('.true.', 'true', 't', '1', 'yes', 'y'):
        return True
    if v in ('.false.', 'false', 'f', '0', 'no', 'n'):
        return False
    raise ValueError(f"Cannot parse bool: {value}")


def chmod(value):
    return int(f"0o{int(value)}", 8)

def getzoombox(values):
    str_to_parse = str(values).strip()
    if str_to_parse.startswith('[') and str_to_parse.endswith(']'):
        str_to_parse = str_to_parse[1:-1]
    values = [float(v.strip()) for v in str_to_parse.split(',')]
    if len(values) != 6:
        raise ValueError(f"Expected 6 values for zoombox, got {len(values)}: {values}")
    return np.array([*values])

PARAMS = {
    # name in H        aliases in input file / CLI                 type
    'omega_f':        (['omega_0', 'Omega_0', 'omega_f'],          np.float64),
    'omega_lambda_f': (['omega_l', 'lambda_0', 'lambda_f'],        np.float64),
    'af':             (['af', 'afinal', 'a_f'],                    np.float64),
    'Lf':             (['Lf', 'lf', 'lbox'],                       np.float64),
    'H_f':            (['H_f', 'H_0', 'H'],                        np.float64),
    'FlagPeriod':     (['FlagPeriod'],                             np.int32),

    'nMembers':       (['n', 'N', 'npart'],                        np.int32),
    'cdm':            (['cdm'],                                    fbool),
    'method':         (['method'],                                 str),
    'b_init':         (['b'],                                      np.float64),
    'nvoisins':       (['nvoisins'],                               np.int32),
    'nhop':           (['nhop'],                                   np.int32),
    'rho_threshold':  (['rhot'],                                   np.float64),
    'fudge':          (['fudge'],                                  np.float64),
    'fudgepsilon':    (['fudgepsilon'],                            np.float64),
    'alphap':         (['alphap'],                                 np.float64),

    'verbose':        (['verbose'],                                fbool),
    'megaverbose':    (['megaverbose'],                            fbool),
    'DPMMC':          (['DPMMC'],                                  fbool),
    'SC':             (['SC'],                                     fbool),
    'dcell_min':      (['dcell_min'],                              np.float64),
    'eps_SC':         (['eps_SC'],                                 np.float64),

    'nsteps':         (['nsteps', 'nsteps_do'],                    np.int32),
    'dump_dms':       (['dump_DMs'],                               fbool),
    'dump_stars':     (['dump_stars'],                             fbool),
    'nchem':          (['nchem'],                                  np.int32),

    'dchmod':         (['dchmod'],                                 chmod),
    'fchmod':         (['fchmod'],                                 chmod),
    'uid':            (['uid'],                                    int),
    'gid':            (['gid'],                                    int),

    'zoomin':         (['zoomin'],                                 fbool),
    'prefix':         (['prefix'],                                 str),
    'zoombox':         (['zoombox'],                               getzoombox),
}


ALIAS_TO_ATTR = {
    alias: attr
    for attr, (aliases, _) in PARAMS.items()
    for alias in aliases
}
