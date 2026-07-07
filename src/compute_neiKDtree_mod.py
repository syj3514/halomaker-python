import halo_defs as H
from halo_defs import mem,frange,maccess, datdump, datload
import numpy as np
from tqdm import tqdm
from num_rec import spline, splines,splines_numba, icellid, icellids, counting_argsort_8, assign_struct_ids
import multiprocessing as mp
ctx = mp.get_context('fork')
Pool = ctx.Pool
import time, os
import faulthandler, signal, sys
from scipy.spatial import cKDTree
from functools import partial

def stop(): raise ValueError("Stop")

#   use halo_defs
#   integer::ncall         #YDdebug
# ncall_create_nodes = 0
# ncall = 0
print_prefix = "    "
ncall = 0  # debug counter for do_colorize_131430; plain int avoids import-time Manager() server (sandbox-safe)
ntest = 0
inccell = 0
ind_deepcount = 0
switch = 0
#   integer::icolor_select #YDdebug
icolor_select = 0
from collections import defaultdict
timers = defaultdict(float)
#   contains

if H.FORTRAN:
    # import create_nodes
    # from create_nodes import neikdtree
    # from create_nodes import fhalo_defs as fH
    from compute_adaptahop import neikdtree

#=======================================================================
def compute_adaptahop_131():
#=======================================================================
    timerecords = []
    ref = time.time()
    change_pos_1310()
    timerecords.append(('    change_pos', -time.time()+ref)); ref = time.time()
    if H.FORTRAN:
        sync_fortran()
        if(H.megaverbose): print(f'{print_prefix}Go into Fortran for full-box...')
        neikdtree.compute_adaptahop(
            np.asfortranarray(mem['pos_10'].T),
            np.asfortranarray(mem['mass_10']))
        timerecords.append(('    compute_adaptahop', -time.time()+ref)); ref = time.time()
        
        H.allocate('density_1312',H.npart, dtype=np.float64)
        mem['whereIam_parts'][:] = neikdtree.liste_parts[:]
        mem['density_1312'][:] = neikdtree.density[:]
        H.nnodes = int(neikdtree.nnodes)
        H.node_0 = np.zeros(H.nnodes+1, dtype=H.node_dtype)
        arr = [neikdtree.real_table.T, neikdtree.integer_table.T]
        node_0 = fromndarrays(arr, dtype=H.node_dtype)
        H.node_0[1:] = node_0
        
        # Close real_table, integer_table, liste_parts, density
        del arr
        neikdtree.close()
        if(H.megaverbose): print(f'{print_prefix}Fortran closed')
        timerecords.append(('    f2py', -time.time()+ref)); ref = time.time()
    else:
        # action neighbors
        if H.SCIPY:
            tree = create_tree_structure_1311_scipy()
            timerecords.append(('    create_tree', -time.time()+ref)); ref = time.time()
            compute_mean_density_and_np_1312_scipy(tree)
            timerecords.append(('    compute_density', -time.time()+ref)); ref = time.time()
        else:
            create_tree_structure_1311()
            timerecords.append(('    create_tree', -time.time()+ref)); ref = time.time()
            compute_mean_density_and_np_1312()
            timerecords.append(('    compute_density', -time.time()+ref)); ref = time.time()
        # action adaptahop
        find_local_maxima_1313()
        timerecords.append(('    find_local_maxima', -time.time()+ref)); ref = time.time()
        create_group_tree_1314()      # MAIN BOTTLENECT
    timerecords.append(('    create_group_tree', -time.time()+ref)); ref = time.time()
    change_pos_back_1315()
    timerecords.append(('    change_pos_back', -time.time()+ref)); ref = time.time()
    # check that we have halos
    count_halos_1316()
    timerecords.append(('    count_halos', -time.time()+ref)); ref = time.time()
    if(H.nb_of_halos>0):
        if(H.verbose): print()
        if(H.verbose): print(f"{print_prefix}Select halos and subhalos...")
        _ref = time.time()
        # reinit halo and subhalo count
        H.nb_of_halos    = 0
        H.nb_of_subhalos = 0
        # H.node_0 structure tree -> halo structure tree: 3 methods avalable
        if(H.method=="HOP"):
            select_halos_1317()
        elif(H.method=="DPM"):
            select_with_DP_method_1318()
        elif(H.method=="MSM"):
            select_with_MS_method_1319()
        elif(H.method=="BHM"):
            select_with_BH_method_131a()
        else:
            print(f'{print_prefix}> could not recognise selection method:',H.method)

        # group whereIam_parts by halo
        whereIam_parts = mem['whereIam_parts']
        pids0 = np.arange(H.npart)
        argsort = np.argsort(whereIam_parts)
        whereIam_parts_sorted = whereIam_parts[argsort]
        H.allocate('pids0_groupsorted', len(whereIam_parts_sorted), dtype=np.int32)
        mem['pids0_groupsorted'] = pids0[argsort]
        unis, idxs, counts = np.unique(whereIam_parts_sorted, return_index=True, return_counts=True)
        H.allocate('whereIam_idxs', len(idxs), dtype=np.int32)
        H.allocate('whereIam_counts', len(counts), dtype=np.int32)
        mem['whereIam_idxs'] = idxs
        mem['whereIam_counts'] = counts
        if(H.megaverbose): print(f"{print_prefix}--> {time.time()-_ref:.2f} seconds to select halos")
        timerecords.append(('    select_halos', -time.time()+_ref))
    else:
        H.node_0 = np.empty(0, dtype=H.node_dtype)
    return timerecords

#=======================================================================
def list_parameters_1300():
#=======================================================================


    print(f'{print_prefix}============================================================')
    print(f'{print_prefix}Compute_neiKDtree was called with the following attributes :')
    print(f'{print_prefix}verbose         :',H.verbose)
    print(f'{print_prefix}megaverbose     :',H.megaverbose)
    print(f'{print_prefix}boxsize2        :',H.boxsize2)
    print(f'{print_prefix}hubble          :',H.hubble)
    print(f'{print_prefix}omega0          :',H.omega0)
    print(f'{print_prefix}omegaL          :',H.omegaL)
    print(f'{print_prefix}aexp_max        :',H.aexp_max)
    print(f'{print_prefix}action          : neighbors')
    print(f'{print_prefix}nvoisins        :',H.nvoisins)
    print(f'{print_prefix}nhop            :',H.nhop)
    print(f'{print_prefix}action          : adaptahop')
    print(f'{print_prefix}rho_threshold   :',H.rho_threshold)
    print(f'{print_prefix}nmembthresh     :',H.nmembthresh)
    print(f'{print_prefix}alphap          :',H.alphap)
    print(f'{print_prefix}fudge           :',H.fudge)
    print(f'{print_prefix}fudgepsilon     :',H.fudgepsilon)
    print(f'{print_prefix}epsilon [kpc]   :',1000*H.fudgepsilon*H.xlong/H.npart**(1/3))
    print(f'{print_prefix}Selection method:',' ',H.method)
    print(f'{print_prefix}============================================================', flush=True)


#=======================================================================
def init_adaptahop_130():
#=======================================================================
    H.omegaL   = H.omega_lambda_f
    H.omega0   = H.omega_f
    H.aexp_max = H.af
    H.hubble   = H.H_f*1e-2
    H.boxsize2 = H.Lf
    H.xlong    = H.boxsize2
    H.ylong    = H.xlong
    H.zlong    = H.xlong
    H.boxarr   = np.array([H.xlong, H.ylong, H.zlong])
    H.xlongs2  = H.xlong*0.5
    H.ylongs2  = H.ylong*0.5
    H.zlongs2  = H.zlong*0.5

    pos_renorm  =H.xlong
    H.pos_shift[0]=0.0
    H.pos_shift[1]=0.0
    H.pos_shift[2]=0.0
    H.mass_in_kg=(H.xlong*H.ylong*H.zlong/H.npart)*H.mega_parsec**3 \
        *H.omega0*H.critical_density*H.hubble**2 
    H.Hub_pt = 100.*H.hubble * np.sqrt(H.omega0*(H.aexp_max/H.aexp)**3 + (1-H.omegaL-H.omega0)*(H.aexp_max/H.aexp)**2 + H.omegaL) 
    
    H.nmembthresh = H.nMembers  
    
    list_parameters_1300()

    if H.FORTRAN:
        neikdtree.sync_from_init_adaptahop(
            # i4
            H.npart,H.nmembthresh,H.nMembers,
            # f8
            H.omegaL, H.omega_lambda_f, 
            H.omega0, H.omega_f, 
            H.aexp_max, H.af, 
            H.hubble, H.H_f, 
            H.boxsize2, H.Lf, 
            H.xlong, H.ylong, H.zlong, 
            H.xlongs2, H.ylongs2, H.zlongs2,
            H.Hub_pt,H.aexp,
            # arr
            H.pos_shift,
        )

#=======================================================================
def change_pos_1310():
#=======================================================================
    H.npart    = H.nusedpart
    H.epsilon  = H.fudgepsilon*H.xlong/H.npart**(1/3)
    if(H.SCIPY)and(not H.FORTRAN): mem['pos_10'] += 0.5
    mem['pos_10']      *= H.boxsize2
    if H.FORTRAN:
        neikdtree.sync_from_change_pos(
            H.npart, H.nusedpart, # i4
            H.epsilon, H.fudgepsilon, H.xlong, H.boxsize2 # f8
            )


#=======================================================================
def change_pos_back_1315():
#=======================================================================
    mem['pos_10'] /= H.boxsize2
    if H.SCIPY: mem['pos_10'] -= 0.5

#=======================================================================
def count_halos_1316():
#=======================================================================
    # H.nb_of_halos = 0
    # for inode1 in frange(1,H.nnodes):
    #     if(H.node_0[inode1]['level']==1): H.nb_of_halos += 1
    H.nb_of_halos = np.count_nonzero(H.node_0['level'][1:]==1)

#=======================================================================
def select_halos_1317():
#=======================================================================
    node_to_halo = np.full(H.nnodes, -1, dtype=np.int32)

    # counting number of halos
    H.nb_of_halos = 0
    for inode1 in frange(1,H.nnodes):
        if(H.node_0[inode1]['mother']<=0): H.nb_of_halos += 1
        ih1 = inode1
        while(H.node_0[ih1]['mother'] !=0):
            ih1 = H.node_0[ih1]['mother'] if(node_to_halo[ih1-1]<=0) else node_to_halo[ih1-1]
        node_to_halo[inode1-1] = ih1 
    print(f'{print_prefix}number of nodes :',H.nnodes)
    print(f'{print_prefix}number of haloes:', H.nb_of_halos)

    for ipar0 in range(H.npart):
        inode1 = mem['whereIam_parts'][ipar0]
        if(inode1>0):
            ih1 = node_to_halo[inode1-1]
            if((ih1<=H.nb_of_halos)and(ih1>0)):
                mem['whereIam_parts'][ipar0] = ih1
            else:
                raise ValueError('ihalo > nb_of_halos or nil')

    H.node_0 = np.empty(0, dtype=H.node_dtype)

#=======================================================================
def select_with_DP_method_1318():
#=======================================================================
    raise NotImplementedError('`select_with_DP_method` not implemented')

#=======================================================================
def select_with_MS_method_1319(mass_acc = 1.0e-2, strict=False):
#=======================================================================
    mostmasssub = np.zeros(H.nnodes, dtype=np.int32)
    node_to_struct = np.zeros(H.nnodes, dtype=np.int32)
    if(H.verbose): print(f'{print_prefix}{print_prefix}Using MS method')

    # First recompute H.node_0 mass accordingly to their particle list
    where = np.where(H.node_0['firstchild']>0)[0]
    if where[0]==0: where=where[1:]
    for inode1 in where:
        me = H.node_0[inode1]
        isub1 = me['firstchild']
        while(isub1>0):
            # Initially, each `mass` includes all particles in the node and its subnodes. 
            # We need to subtract the mass of the subnodes to get the mass of the node itself.
            me['mass'] -= H.node_0[isub1]['mass']
            me['truemass'] -= H.node_0[isub1]['truemass']
            isub1 = H.node_0[isub1]['sister']
        if(me['mass']<=0 or me['truemass']<=0.0):
            print(f'{print_prefix}Error in computing H.node_0',inode1,'mass_10')
            raise ValueError(' mass, truemass:', me['mass'], me['truemass'])

    if(H.verbose)and(strict):
        # check that the new masses are correct
        npartcheck_0 = np.zeros(H.nnodes+1, dtype=np.int32)
        for ip1 in frange(1,H.nusedpart):
            if(mem['whereIam_parts'][ip1-1]<0): raise ValueError('whereIam_parts is smaller than 0')
            npartcheck_0[mem['whereIam_parts'][ip1-1]] += 1
        if(sum(npartcheck_0)!=H.nusedpart): raise ValueError('Error in particles count')
        for inode1 in frange(1,H.nnodes):
            if(H.node_0[inode1]['mass']!=npartcheck_0[inode1]):
                print(f'{print_prefix}Error in H.node_0 particle count, for H.node_0',inode1)
                print(f'{print_prefix}it first subnode is:',H.node_0[inode1]['firstchild'])
                if(H.node_0[inode1]['firstchild']>0): print(f'{print_prefix}it has:',H.node_0[H.node_0[inode1]['firstchild']]['nsisters'],'subnodes' )
                raise ValueError("")
        del npartcheck_0
    
    mostmasssub[:]    = -1
    firstchild = H.node_0['firstchild']
    for inode1 in frange(H.nnodes, 1, -1):
        # Search inode1 subnodes for the most massive one
        isub1 = firstchild[inode1]
        if(isub1>0):
            # init all on the first subnode
            maxmasssub         = H.node_0[isub1]['truemass']
            mostmasssub[inode1-1] = isub1
            isub1               = H.node_0[isub1]['sister']
        while(isub1>0):
            if(H.node_0[isub1]['truemass']>maxmasssub):
                maxmasssub         = H.node_0[isub1]['truemass']
                mostmasssub[inode1-1] = isub1
            isub1               = H.node_0[isub1]['sister']
        # add mostmasssub mass to inode1 mass and inode1's densmax = mostmasssub's densmax
        if (mostmasssub[inode1-1] > 0):
            H.node_0[inode1]['mass']     += H.node_0[mostmasssub[inode1-1]]['mass']
            H.node_0[inode1]['truemass'] += H.node_0[mostmasssub[inode1-1]]['truemass']
            H.node_0[inode1]['densmax']  = H.node_0[mostmasssub[inode1-1]]['densmax']

    # Second run to write node_to_struct array and count the number of substructures
    mothers = H.node_0['mother'][1:H.nnodes+1]
    indices = np.arange(1, H.nnodes + 1)
    is_halo = (mothers <= 0)
    valid_mother_idx = np.where(~is_halo, mothers - 1, 0) 
    is_most_mass = (~is_halo) & (indices == mostmasssub[valid_mother_idx])
    is_new_sub = (~is_halo) & (~is_most_mass)
    nb_halos = np.sum(is_halo)
    nb_sub = np.sum(is_new_sub)

    # Assign a struct ID to each node
    is_starter = is_halo | is_new_sub
    starter_ids = np.zeros(H.nnodes, dtype=int)
    starter_ids[is_starter] = np.cumsum(is_starter)[is_starter]

    # Define the parent-child relationship: the most massive child points to its mother
    parent = indices-1 # 0-based
    parent[is_most_mass] = mothers[is_most_mass] - 1

    # Pointer jumping
    # : All nodes will point to the 'starter node' of their branch
    for _ in range(int(np.ceil(np.log2(H.nnodes)))):
        new_parent = parent[parent]
        if np.array_equal(parent, new_parent):
            break
        parent = new_parent
    node_to_struct[:] = starter_ids[parent]
    
    # check halos and subhalos number 
    H.nstruct = nb_halos + nb_sub
    if strict:
        istruct = np.max(starter_ids)
        if(H.nstruct!=istruct):
            raise ValueError('Error in subroutines count,',H.nstruct, istruct)

    if(H.verbose):
        # Check that the structure tree is ordered as it should be 
        assert (H.node_0['level'][1:nb_halos+1]==1).all(), 'error in H.node_0 order'
        assert (node_to_struct[:nb_halos] == indices[:nb_halos]).all(), 'error in H.node_0 order'
        assert (H.node_0['level'][nb_halos+1:]!=1).all(), 'error shouldn''t be any halo here'
        # for inode1 in range(1,H.nnodes):
        #     if(inode1<=nb_halos):
        #         if( not( (H.node_0[inode1]['level']==1)and(node_to_struct[inode1-1]==inode1) )):
        #             raise ValueError('error in H.node_0 order' )
        #     else:
        #         if(H.node_0[inode1]['level']==1):
        #             raise ValueError('error shouldn''t be any halo here' )
        print(f"{print_prefix}{print_prefix}")
        print(f'{print_prefix}{print_prefix}> number of nodes               :', H.nnodes)
        print(f'{print_prefix}{print_prefix}> number of halos               :', nb_halos)
        print(f'{print_prefix}{print_prefix}> number of substructures       :', nb_sub)
        print(f'{print_prefix}{print_prefix}> number of node removed        :', H.nnodes - H.nstruct)
        print(f"{print_prefix}{print_prefix}")
        print(f'{print_prefix}{print_prefix}> Cleaning whereIam_parts')
        if(strict):
            npartcheck_0 = np.zeros(H.nstruct, dtype=np.int32)#; npartcheck_0[:] = 0

    # Cleaning whereIam_parts
    ip0s = np.where(mem['whereIam_parts']>0)[0]
    inode1s = mem['whereIam_parts'][ip0s]
    mem['whereIam_parts'][ip0s] = node_to_struct[inode1s-1]

    if(H.verbose)and(strict):
        # Check H.node_0[inode]['mass'] it should now correspond to npartcheck_0 count
        for inode0 in range(H.nnodes):
            if(node_to_struct[inode0]<=0): raise ValueError('node_to_struct is nil')
            imother1 = H.node_0[inode0+1]['mother']
            if((imother1<=0) or (imother1>0 and (node_to_struct[imother1-1]!=node_to_struct[inode0]))):
                if(H.node_0[inode0+1]['mass']!=npartcheck_0[node_to_struct[inode0]-1]):
                    print(f'{print_prefix}Wrong nb of part in struct: ', node_to_struct[inode0])
                    print(f"{print_prefix}inode0,H.node_0[inode0+1]['mass'],istruct,npartcheck_0(istruct)",inode0, \
                        H.node_0[inode0+1]['mass'],node_to_struct[inode0],npartcheck_0[node_to_struct[inode0]-1])
                    raise ValueError("")
        del npartcheck_0

    if(H.verbose): print(f'{print_prefix}{print_prefix}> Creating new structure tree')
    # creating new structure tree
    H.allocate('mother_1319', H.nstruct, dtype=np.int32)
    H.allocate('first_sister_1319', H.nstruct, dtype=np.int32)
    H.allocate('first_daughter_1319', H.nstruct, dtype=np.int32)
    H.allocate('level_1319', H.nstruct, dtype=np.int32)

    first_sister_1319 = mem['first_sister_1319']
    mother_1319 = mem['mother_1319']
    first_daughter_1319 = mem['first_daughter_1319']
    level_1319 = mem['level_1319']
    mother_1319[:]         = -1
    first_sister_1319[:]   = -1
    first_daughter_1319[:] = -1
    level_1319[:]          =  0
    ihalo1          = -1 

    for inode0 in range(H.nnodes):
        istruct1 = node_to_struct[inode0]
        if(istruct1<=0): raise ValueError('index nil for istrut')
        # If not touched
        if(mother_1319[istruct1-1]<0):
            # I'm a main halo
            if(H.node_0[inode0+1]['mother']<=0):
                mother_1319[istruct1-1] = 0
                if(ihalo1>0): first_sister_1319[ihalo1-1] = istruct1
                ihalo1 = istruct1
                level_1319[istruct1-1] = 1
            # I'm a subhalo
            else:
                imother1          = node_to_struct[H.node_0[inode0+1]['mother']-1]
                mother_1319[istruct1-1]  = imother1
                level_1319[istruct1-1]   = level_1319[imother1-1] + 1
                # I'm first daughter of my mother
                if(first_daughter_1319[imother1-1]<=0):
                    first_daughter_1319[imother1-1] = istruct1
                # I have already a sister, I need to be the last sister
                else:
                    isub1 = first_daughter_1319[imother1-1]
                    while(first_sister_1319[isub1-1]>0):
                        isub1 = first_sister_1319[isub1-1]
                    first_sister_1319[isub1-1] = istruct1

    if(H.megaverbose):
        # For test we shall output the structure tree in file struct_tree.dat'
        f120 = open("struct_tree.dat", mode='w+')
        for istruct0 in range(H.nstruct):
            if(mem['mother_1319'][istruct0]<=0):
                f120.write('---------\n')
                f120.write(
                    f"halo:{istruct0+1:6d} first_child:{mem['first_daughter_1319'][istruct0]:6d} sister:{mem['first_sister_1319'][istruct0]:6d}\n")
                isub1 = mem['first_daughter_1319'][istruct0]
                while(isub1>0):
                    f120.write(
                        f"sub:{isub1:6d} mother:{mem['mother_1319'][isub1-1]:6d} first_child:{mem['first_daughter_1319'][isub1-1]:6d} sister:{mem['first_sister_1319'][isub1-1]:6d}\n")
                    isub1 = mem['first_sister_1319'][isub1-1]
            else:
                isub1 = mem['first_daughter_1319'][istruct0]
                if(isub1!=0):
                    f120.write(
                        f"sub:{istruct0+1:6d} mother:{mem['mother_1319'][istruct0]:6d} first_child:{mem['first_daughter_1319'][istruct0]:6d} sister:{mem['first_sister_1319'][istruct0]:6d}\n")
                    while(isub1>0):
                        f120.write(
                            f"sub:{isub1:6d} mother:{mem['mother_1319'][isub1-1]:6d} first_child:{mem['first_daughter_1319'][isub1-1]:6d} sister:{mem['first_sister_1319'][isub1-1]:6d}\n")
                        isub1 = mem['first_sister_1319'][isub1-1]
        f120.close()
        full_path = os.path.abspath('struct_tree.dat')
        os.chmod(full_path, H.fchmod); os.chown(full_path, H.uid, H.gid)

    H.nb_of_halos    = nb_halos
    H.nb_of_subhalos = nb_sub
    H.node_0 = np.empty(0, dtype=H.node_dtype)

#=======================================================================
def select_with_BH_method_131a():
#=======================================================================
    raise NotImplementedError('`select_with_BH_method` not implemented')

#=======================================================================
def main_prog(inode,mainprog):
#=======================================================================
    raise NotImplementedError("`main_prog` is used in BHM which is not implemented")

#=======================================================================
def main_prog_node_sub(inode,isub,mainprog):
#=======================================================================
    raise NotImplementedError("`main_prog_node_sub` is used in BHM which is not implemented")

#=======================================================================
def make_linked_node_list():
#=======================================================================
    raise NotImplementedError("`make_linked_node_list` is used in BHM which is not implemented")

def _compute_mean_density_worker(ipar1, poshere):
    dist2_0, iparnei = walk_tree_norec_131200(ipar1, poshere) # <---- Main bottleneck
    if(np.sum(iparnei)==0):
        print(ipar1)
        print(dist2_0)
        print(iparnei)
        raise ValueError('No neighbors found for particle', ipar1)
    compute_density_13121(ipar1,dist2_0,iparnei)
    mem['iparneigh_1312'][:H.nhop,ipar1-1]=iparnei[:H.nhop]

def _compute_mean_density_chunk(ipar1start, ipar1end):
    poss = mem['pos_10'][ipar1start-1:ipar1end]
    for ipar1, ipos in zip(range(ipar1start, ipar1end+1), poss):
        _compute_mean_density_worker(ipar1, ipos)
    
#=======================================================================
def compute_mean_density_and_np_1312_scipy(tree):
#=======================================================================
    if (H.verbose):
        print(f'{print_prefix}Compute mean density for each particle...')
        ref = time.time()

    H.allocate('iparneigh_1312',(H.nhop,H.npart), dtype=np.int32)

    # Neighbors
    _ref = time.time()
    dist1_0s, idx0 = tree.query(mem['pos_10'], k=H.nhop+1, workers=H.nbPes)
    idx0 = idx0[:,1:]
    mem['iparneigh_1312'][:,:] = idx0.T+1
    print(f"iparneigh ({time.time()-_ref:.6f} sec)")

    # Density
    _ref = time.time()
    H.allocate('density_1312',H.npart, dtype=np.float64)
    rs = dist1_0s[:,H.nvoisins]*0.5
    imass = mem['mass_10']

    tmp = dist1_0s[:,1:H.nvoisins+1] / rs[:,None]
    splined = splines_numba(tmp)
    contrib = np.sum(imass[idx0] * splined, axis=1) + imass
    mem['density_1312'][:] = (H.xlong*H.ylong*H.zlong)*contrib /(H.pi*rs**3)
    print(f"density ({time.time()-_ref:.6f} sec)")

    
    # Check for average density
    if (H.verbose):
        _ref = time.time()
        print(f"{print_prefix}    Calc average density...")
        densav = mem['density_1312'].mean()
        print(f'{print_prefix}--> Average density :',densav)
        print(f"avgdensity ({time.time()-_ref:.6f} sec)")
        print(f"{print_prefix}--> {time.time()-ref:.2f} seconds to compute mean density")
    
#=======================================================================
def compute_mean_density_and_np_1312():
#=======================================================================
    if (H.verbose): print(f'{print_prefix}Compute mean density for each particle...')
    _ref = time.time()

    H.allocate('iparneigh_1312',(H.nhop,H.npart), dtype=np.int32)
    H.allocate('density_1312',H.npart, dtype=np.float64)

    DEBUG=False
    if (DEBUG)and(os.path.exists('density_1312.tmp')):
        print(f'{print_prefix}Loading density and iparneigh from tmp files...')
        tmp = datload('density_1312.tmp', dtype='f8')
        mem['density_1312'][:] = tmp[:]
        tmp = datload('iparneigh_1312.tmp', dtype='i4').reshape((H.nvoisins, H.npart))
        mem['iparneigh_1312'][:,:] = tmp[:,:]
        
        # Check for average density
        if (H.verbose):
            print(f"{print_prefix}    Calc average density...")
            densav = mem['density_1312'].mean()
            print(f'{print_prefix}--> Average density :',densav)
        
        H.deallocate('mass_cell_1311')
        H.deallocate('size_cell_1311')
        H.deallocate('pos_cell_1311')
        H.deallocate('sister_1311')
        H.deallocate('firstchild_1311')
    else:
        CHUNK = True
        pos = mem['pos_10']
        if (H.nbPes == 1)or(DEBUG): # Single process
            iterator = tqdm(zip(frange(1,H.npart), pos), total=H.npart, desc=f'{print_prefix}    [Ncpu=1] Compute density', ncols=100)
            for ipar1, ipos in iterator:
                _compute_mean_density_worker(ipar1, ipos)
            iterator.close()
        else:
            if CHUNK:
                pbar = tqdm(total=H.npart, desc=f'{print_prefix}    [Ncpu={H.nbPes}] Compute density', mininterval=0.5, ncols=100)
                signal.signal(signal.SIGTERM, signal.SIG_DFL)
                # Split the work into chunks to reduce overhead
                tmp = np.linspace(1, H.npart+1, H.nbPes*100+1, dtype=np.int32)
                ipar1starts = tmp[:-1]
                ipar1ends   = tmp[1:]
                with Pool(processes=H.nbPes) as pool:
                    async_results = []
                    for ipar1start, ipar1end in zip(ipar1starts, ipar1ends):
                        r = pool.apply_async(_compute_mean_density_chunk, args=(ipar1start, ipar1end), callback=lambda _: pbar.update(ipar1end - ipar1start+1))
                        async_results.append(r)
                    for r in async_results:
                        r.get()
                signal.signal(signal.SIGTERM, H.flush)
                pbar.close()
            else:
                pbar = tqdm(total=H.npart, desc=f'{print_prefix}    [Ncpu={H.nbPes}] Compute density', mininterval=0.5, ncols=100)
                signal.signal(signal.SIGTERM, signal.SIG_DFL)
                with Pool(processes=H.nbPes) as pool:
                    async_results = []
                    for ipar1, ipos in zip(range(1, H.npart+1), pos):
                        r = pool.apply_async(_compute_mean_density_worker, args=(ipar1, ipos), callback=lambda _: pbar.update())
                        async_results.append(r)
                    for r in async_results:
                        r.get()
                signal.signal(signal.SIGTERM, H.flush)
                pbar.close()

        # Check for average density
        if (H.verbose):
            print(f"{print_prefix}    Calc average density...")
            densav = mem['density_1312'].mean()
            print(f'{print_prefix}Average density :',densav)
        
        H.deallocate('mass_cell_1311')
        H.deallocate('size_cell_1311')
        H.deallocate('pos_cell_1311')
        H.deallocate('sister_1311')
        H.deallocate('firstchild_1311')

        # datdump(mem['iparneigh_1312'].flatten(), 'iparneigh_1312.tmp')
        # datdump(mem['density_1312'], 'density_1312.tmp')
    print(f"{print_prefix}--> {time.time()-_ref:.2f} seconds to compute mean density")
    return


#=======================================================================
def find_local_maxima_1313():
#=======================================================================
    if (H.verbose): print(f'{print_prefix}Now Find local maxima...')

    mem['idpart_1311'][:]=0
    H.ngroups=0
    density = mem['density_1312']
    iparneigh = mem['iparneigh_1312']
    arange0 = np.arange(H.npart)
    arange1 = arange0+1
    denmask = density > H.rho_threshold
    darr = np.empty((H.nhop+1, H.npart), dtype=density.dtype)
    darr[0,:] = density; darr[1:,:] = density[iparneigh - 1]
    narr = np.empty((H.nhop+1, H.npart), dtype=iparneigh.dtype)
    narr[0,:] = arange1; narr[1:,:] = iparneigh
    argmax = np.argmax(darr, axis=0)
    densest = darr[argmax, arange0]
    iparsel = narr[argmax, arange0]
    iparsel[~denmask] = 0

    whoistie = darr == densest
    nwhoistie = np.sum(whoistie, axis=0) > 1
    if nwhoistie.any():
        if (H.verbose):
            print(f'{print_prefix}WARNING : equal densities in find_local_maxima.')
        where = np.where(nwhoistie)[0]
        for iw in where:
            candidates = narr[whoistie[:, iw], iw]
            iparsel[iw] = np.min(candidates)
    
    iam_max = (iparsel==arange1)&(denmask)
    H.ngroups = np.sum(iam_max)
    mem['idpart_1311'][iam_max] = -np.arange(1, H.ngroups+1)
    mem['idpart_1311'][~iam_max] = iparsel[~iam_max]
    if (H.verbose): print(f'{print_prefix}--> Number of local maxima found :',H.ngroups)



    # Now Link the particles associated to the same maximum
    if (H.verbose): print(f'{print_prefix}Create link list...')

    H.allocate('densityg_1313',H.ngroups, dtype=np.float64)
    nmemb = np.zeros(H.ngroups, dtype=np.int32)
    H.allocate('firstpart_1313',H.ngroups, dtype=np.int32)
    H.allocate('igrouppart_1313', H.npart, dtype=np.int32)

    mem['igrouppart_1313'][:] = 0
    wherepeak = np.where(mem['idpart_1311']<0)[0]
    mem['densityg_1313'][-mem['idpart_1311'][wherepeak]-1] = density[wherepeak]
    iparid1s = mem['idpart_1311'].copy()
    active = iparid1s > 0
    while active.any():
        nxt = mem['idpart_1311'][iparid1s[active] - 1]
        iparid1s[active] = nxt
        active = iparid1s > 0
    mem['igrouppart_1313'][:] = -iparid1s
    mem['igrouppart_1313'][~denmask] = 0

    nmemb[:] = 0; mem['firstpart_1313'][:] = 0
    argsort = np.argsort(mem['igrouppart_1313'][denmask], kind='stable')
    tmp = mem['igrouppart_1313'][denmask][argsort]
    _ipar0 = arange0[denmask][argsort]
    ind = np.flatnonzero(np.r_[True, tmp[1:] != tmp[:-1]])
    counts = np.diff(np.r_[ind, tmp.size])
    uni = tmp[ind]
    mem['idpart_1311'][_ipar0] = 0
    mem['idpart_1311'][_ipar0[1:]] = _ipar0[:-1] + 1
    mem['idpart_1311'][_ipar0[ind]] = 0
    nmemb[uni-1] = counts
    mem['firstpart_1313'][uni - 1] = _ipar0[ind + counts - 1] + 1

    nmembmax=np.max(nmemb)
    nmembtot=np.sum(nmemb)
    if (H.verbose):
        print(f'{print_prefix}--> Number of particles of the largest group :',nmembmax)
        print(f'{print_prefix}--> Total number of particles in groups ',nmembtot)
    del nmemb


getgroup_material = None
def get_group_members(igroup1, verbose=False):
    global getgroup_material
    if getgroup_material is None:
        igrouppart = mem['igrouppart_1313']
        density = mem['density_1312']
        denmask = (density > H.rho_threshold)&(igrouppart>0)
        tarr = igrouppart[denmask]
        argsort = np.argsort(tarr, kind='stable')
        tmp = tarr[argsort]
        _ipar1 = (np.arange(H.npart)+1)[denmask][argsort]
        ind = np.flatnonzero(np.r_[True, tmp[1:] != tmp[:-1]])
        counts = np.diff(np.r_[ind, tmp.size])
        getgroup_material = (_ipar1, ind, counts)
    _ipar1, ind, counts = getgroup_material
    i0 = ind[igroup1-1]; i1 = i0+counts[igroup1-1]
    if verbose:
        igrouppart = mem['igrouppart_1313']
        alltrue = (igroup1==igrouppart[_ipar1[i0:i1]-1]).all()
        if alltrue: print(f"{print_prefix}# Good")
        else:
            print()
            print(f"{print_prefix}# {igroup1}")
            print(f"{print_prefix}# {igrouppart[_ipar1[i0:i1]-1]}")
    return _ipar1[i0:i1]

def set_group_members(ipar1s, igroup1s):
    global getgroup_material
    mem['igrouppart_1313'][ipar1s-1] = igroup1s

def init_group_members():
    global getgroup_material
    getgroup_material = None


def sync_fortran():
    assert H.FORTRAN
    neikdtree.sync_others(
        H.verbose, H.megaverbose, H.npart, H.nbPes,
        H.rho_threshold, H.massp, H.boxsize,
        H.nhop, H.nvoisins,
        H.fudge, H.alphap,
        H.method, H.nlevelmax,
    )

def get_bytes_data(array):
    # works as non-copy view of structured array if array is contiguous
    barr = array.view('b').reshape((array.size, array.itemsize))
    return barr

def fromndarrays(ndarrays, dtype):
    """
    convert list of ndarray to structured array with given dtype
    faster than np.rec.fromarrays
    only works for 2d arrays for now
    """
    descr = np.dtype(dtype)

    itemsize = 0
    nitem = None
    for nda in ndarrays:
        if (nitem is None):
            nitem = nda.shape[0]
        elif (nitem != nda.shape[0]):
            raise ValueError("Array shape does not match")
        itemsize += nda.shape[1] * nda.dtype.itemsize
    if (descr.itemsize != itemsize):
        raise ValueError(f"Sum of itemsize ({itemsize}) does not match with desired dtype ({descr.itemsize})")

    array = np.empty(nitem, descr)
    barr = get_bytes_data(array)
    col = 0
    for nda in ndarrays:
        bnda = nda.view('b')
        barr[:, col:col + bnda.shape[1]] = bnda
        col += bnda.shape[1]
    return array


#=======================================================================
def create_group_tree_1314():
#=======================================================================
    if (H.verbose): print(f'{print_prefix}Create the tree of structures of structures')

    # End of the branches of the tree
    ref = time.time()
    compute_saddle_list_13140()
    # Output: H.group
    if (H.verbose): print(f"{print_prefix}--> {time.time()-ref:.2f} seconds to saddle_list")

    if (H.verbose): print(f'{print_prefix}Build the hierarchical tree')

    H.nnodesmax=2*H.ngroups

    # Allocations
    H.node_0 = np.zeros(H.nnodesmax+1, dtype=H.node_dtype)
    H.allocate('idgroup_1314',H.ngroups, dtype=np.int32)
    H.allocate('color_1314',H.ngroups, dtype=np.int32)
    H.allocate('igroupid_1314',H.ngroups, dtype=np.int32)
    H.allocate('idgroup_tmp_1314',H.ngroups, dtype=np.int32)

    # Initializations
    mem['whereIam_parts'][:] = 0

    # Iterative loop to build the rest of the tree
    inode               = 0
    H.nnodes            = 0
    rhot                = H.rho_threshold
    H.node_0[inode]['mother']  = 0
    mass_loc            = 0
    truemass            = 0
    igroupref           = 0
    masstmp             = 0
    rsquare             = 0
    densmoy             = 0
    truemasstmp         = 0

    igrouppart = mem['igrouppart_1313']
    mass = mem['mass_10']
    groupmask = igrouppart > 0
    mass_loc = np.count_nonzero(groupmask)
    truemass = np.sum(mass[groupmask])

    H.node_0[inode]['mass']          = mass_loc
    H.node_0[inode]['truemass']      = truemass
    H.node_0[inode]['radius']        = 0
    H.node_0[inode]['density']       = 0
    H.node_0[inode]['px']  = 0
    H.node_0[inode]['py']  = 0
    H.node_0[inode]['pz']  = 0
    H.node_0[inode]['densmax']       = np.max(mem['densityg_1313'])
    H.node_0[inode]['rho_saddle']    = 0.
    H.node_0[inode]['level']         = 0
    H.node_0[inode]['nsisters']      = 0
    H.node_0[inode]['sister']        = 0
    igrA = 1
    igrB = H.ngroups
    mem['idgroup_1314'][:] = np.arange(1, H.ngroups+1)
    mem['igroupid_1314'][:] = np.arange(1, H.ngroups+1)

    # H.allocate('timer_1314', 10, dtype=np.float64)
    # mem['timer_1314'][:]=0
    signal.signal(signal.SIGTERM, signal.SIG_DFL)
    print(f"{print_prefix}Create nodes ...")
    cref = time.time()
    ### MAIN BOTTLENECK!!!
    with Pool(processes=H.nbPes) as global_pool:
        create_nodes_13143(rhot, inode, igrA, igrB, pool=global_pool)
    print(f"{print_prefix}--> {time.time()-cref:.6f} seconds to create nodes")
    signal.signal(signal.SIGTERM, H.flush)
    # for i, itime in enumerate(mem['timer_1314']):
    #     print(f"{print_prefix}    Timer {i}: {itime:.6f} seconds")

    # H.deallocate('timer_1314')
    H.deallocate('igrouppart_1313')
    H.deallocate('idgroup_1314')
    H.deallocate('color_1314')
    H.deallocate('igroupid_1314')
    H.deallocate('idgroup_tmp_1314')
    H.deallocate('idpart_1311')
    H.deallocate('densityg_1313')
    H.deallocate('firstpart_1313')

def _create_nodes_worker(args):
    icolor0, igrA, igrB, rhot = args
    posg=np.zeros(3)
    posref = None
    mass_loc=0
    truemass=0
    rsquareg=0
    densmoyg=0
    densmaxgroup=-1.
    mass_comp=0
    densmoy_comp_max=-1.
    igroupref=0

    idgroup = mem['idgroup_1314']
    densityg = mem['densityg_1313']
    for igr1 in frange(igrA, igrB):
        igroup1=idgroup[igr1-1]
        densmaxgroup=max(densmaxgroup,densityg[igroup1-1])
        posgtmp,masstmp,igroupref,posref,rsquaretmp,densmoytmp,truemasstmp = treat_particles_new_13141(
            igroup1,rhot,
            igroupref,posref)
        posg[:] += posgtmp[:]
        rsquareg += rsquaretmp
        mass_loc += masstmp
        truemass += truemasstmp
        densmoyg += densmoytmp
        densmoytmp /= masstmp
        mass_comp=max(mass_comp,masstmp)
        if (masstmp > 0): densmoy_comp_max=max(densmoy_comp_max, densmoytmp/(1+H.fudge/np.sqrt(masstmp)))

    rsquare_val=np.sqrt(np.abs((truemass*rsquareg-(posg[0]**2+posg[1]**2+posg[2]**2) )/ truemass**2 ))
    densmoy_val=densmoyg/mass_loc

    ifok_val=(
        # 1) Minimum number of particles
        (mass_loc>=H.nmembthresh) and 
        # 2) Density at least equal to the threshold (with fudge factor)
        (densmoy_val > rhot*(1+H.fudge/np.sqrt(mass_loc)) or densmoy_comp_max>rhot) and
        # 3) Density at least equal to the alpha factor times the density of the mother node
        (densmaxgroup >= H.alphap*densmoy_val) and
        # 4) Radius at least equal to the softening length (to avoid numerical artifacts)
        (rsquare_val >= H.epsilon)
    )
    return (icolor0, mass_loc, truemass, posg, densmaxgroup, densmoy_comp_max, mass_comp, rsquare_val, densmoy_val, ifok_val)

#=======================================================================
def create_nodes_13143(rhot,inode,igrA,igrB, pool=None):
#=======================================================================
    global ncall, icolor_select#, ncall_create_nodes
    # timer=mem['timer_1314']; ref = time.time(); ith = 0
    mem['color_1314'][igrA-1:igrB] = 0
    densityg = mem['densityg_1313']
    # Percolate the groups
    icolor_select = 0
    queue = []
    idgroup = mem['idgroup_1314']
    igroupid = mem['igroupid_1314']
    color = mem['color_1314']
    for igr1 in range(igrA, igrB+1):
        igroup1 = idgroup[igr1-1]
        if color[igr1-1] == 0:
            icolor_select += 1
            queue.append(igr1); color[igr1-1]=-1
        while queue:
            igr1 = queue.pop(0)
            igroup1 = idgroup[igr1-1]
            me = H.group[igroup1-1]
            my_isad_gr = me.isad_gr
            my_rho_saddle_gr = me.rho_saddle_gr
            color[igr1-1] = icolor_select
            for inA0 in range(me.nhnei):
                if my_rho_saddle_gr[inA0] > rhot:
                    igroup2 = my_isad_gr[inA0]
                    igr2 = igroupid[igroup2-1]
                    color2 = color[igr2-1]
                    if color2 == 0:
                        queue.append(igr2); color[igr2-1] = -1
                    elif (color2>0) and (color2 != icolor_select):
                        raise ValueError('The connections are not symmetric.')
                    else:
                        pass
                else:
                    # We do not need this saddle anymore 
                    # (and use the fact that saddles are ranked in decreasing order)
                    me.nhnei -= 1
    # timer[ith] += time.time()-ref; ith+=1; ref=time.time()
    # # 12.54 sec

    # We select only groups where we are sure of having at least one
    # particle above the threshold density rhot
    # Then sort them to gather them on the list
    igrpos_0 = np.empty(icolor_select+1, dtype=np.int32)
    igrinc = np.empty(icolor_select, dtype=np.int32)
    igrpos_0[0]=igrA-1
    igrpos_0[1:icolor_select+1]=0
    fslice = slice(igrA-1, igrB)
    dmask = densityg[idgroup[fslice]-1] > rhot
    icolors = color[fslice][dmask]
    np.add.at(igrpos_0, icolors, 1)    
    igrpos_0 = np.cumsum(igrpos_0)    
    if (igrpos_0[icolor_select]-igrA+1 == 0):
        print(f'{print_prefix}ERROR in create_nodes :')
        raise ValueError('All subgroups are below the threshold.')    
    # timer[ith] += time.time()-ref; ith+=1; ref=time.time()
    # # 1.39 sec

    ## Original Loop
    igrinc[:icolor_select]=0
    for igr1 in frange(igrA, igrB):
        icolor1=color[igr1-1]
        igroup1=idgroup[igr1-1]
        if (densityg[igroup1-1]>rhot):
            igrinc[icolor1-1] += 1
            igr_eff1=igrinc[icolor1-1]+igrpos_0[icolor1-1]
            mem['idgroup_tmp_1314'][igr_eff1-1]=igroup1
            mem['igroupid_1314'][igroup1-1]=igr_eff1
    igrB=igrpos_0[icolor_select]
    mem['idgroup_1314'][igrA-1:igrB]=mem['idgroup_tmp_1314'][igrA-1:igrB]
    idgroup = mem['idgroup_1314']
    inc_color_tot = np.sum(igrinc[:icolor_select]>0)

    igrposnew_0 = np.empty(1+inc_color_tot, dtype=np.int32)
    igrposnew_0[0] = igrpos_0[0]
    igrmask = igrinc>0
    nmask = np.sum(igrmask)
    igrposnew_0[1:1+nmask] = igrpos_0[1:][igrmask]
    # timer[ith] += time.time()-ref; ith+=1; ref=time.time()
    # # 3.50 sec

    del igrpos_0
    del igrinc

    isisters=0
    posgg = np.empty((3, inc_color_tot), dtype=np.float64)
    massg = np.empty(inc_color_tot, dtype=np.int32)
    truemassg = np.empty(inc_color_tot, dtype=np.float64)
    densmaxg = np.empty(inc_color_tot, dtype=np.float64)
    densmoy_comp_maxg = np.empty(inc_color_tot, dtype=np.float64)
    mass_compg = np.empty(inc_color_tot, dtype=np.int32)
    rsquare = np.empty(inc_color_tot, dtype=np.float64)
    densmoy = np.empty(inc_color_tot, dtype=np.float64)
    ifok = np.empty(inc_color_tot, dtype=np.bool_)
    ifok[:inc_color_tot]=False
    # timer[ith] += time.time()-ref; ith+=1; ref=time.time()
    # # 0.17 sec

    ### MAIN BOTTLENECK!!! (~96.66 sec)
    DEBUG = False
    if H.nbPes==1 or DEBUG or (pool is None):
        iterator = range(inc_color_tot)
        # if H.verbose: iterator = tqdm(range(inc_color_tot), desc=f'{print_prefix}[nbPes=1] Gathering same color')
        for icolor0 in iterator:
            igrA=igrposnew_0[icolor0]+1
            igrB=igrposnew_0[icolor0+1]
            args = (icolor0, igrA, igrB, rhot)
            res = _create_nodes_worker(args)
            idx = res[0]
            massg[idx], truemassg[idx], posgg[:3, idx], densmaxg[idx], densmoy_comp_maxg[idx], mass_compg[idx], rsquare[idx], densmoy[idx], ifok[idx] = res[1:]
        # if H.verbose: iterator.close()
    else:
        callback=None
        # if H.verbose:
        #     pbar = tqdm(total=inc_color_tot, desc=f'{print_prefix}[nbPes={H.nbPes}] Gathering same color', mininterval=0.5, leave=True)
        #     callback = lambda _: pbar.update()
        async_results = []
        for icolor0 in range(inc_color_tot):
            igrA=igrposnew_0[icolor0]+1
            igrB=igrposnew_0[icolor0+1]
            args = (icolor0, igrA, igrB, rhot)
            async_result = pool.apply_async(_create_nodes_worker, args=(args,), callback=callback)
            async_results.append(async_result)
        for r in async_results:
            res = r.get()
            idx = res[0]
            massg[idx], truemassg[idx], posgg[:3, idx], densmaxg[idx], densmoy_comp_maxg[idx], mass_compg[idx], rsquare[idx], densmoy[idx], ifok[idx] = res[1:]
        # if H.verbose: pbar.close()
    nsisters = np.sum(ifok)
    icolor_ref = np.where(ifok)[0][-1] if nsisters > 0 else 0
    # timer[ith] += time.time()-ref; ith+=1; ref=time.time()
    # # 96.66 sec

    if (nsisters>1):
        isisters=0
        inodetmp=H.nnodes+1
        for icolor0 in range(inc_color_tot):
            if (ifok[icolor0]):
                isisters += 1
                H.nnodes += 1
                if (H.nnodes>H.nnodesmax):
                    print(f'{print_prefix}ERROR in create_nodes :')
                    raise ValueError(f'H.nnodes({H.nnodes}) > H.nnodes max({H.nnodesmax})')
                if((H.nnodes%max(H.nnodesmax/10000,1))==0 and H.megaverbose): print(f'{print_prefix}H.nnodes=',H.nnodes)
                H.node_0[H.nnodes]['mother']=inode
                H.node_0[H.nnodes]['densmax']=densmaxg[icolor0]
                if (isisters>1): H.node_0[H.nnodes]['sister']=H.nnodes-1
                else: H.node_0[H.nnodes]['sister']=0

                H.node_0[H.nnodes]['nsisters']=nsisters
                H.node_0[H.nnodes]['mass']=massg[icolor0]
                H.node_0[H.nnodes]['truemass']=truemassg[icolor0]
                if (massg[icolor0]==0):
                    print(f'{print_prefix}ERROR in create_nodes :')
                    raise ValueError(f'NULL mass for H.nnodes={H.nnodes}')
                posfin = posgg[:3,icolor0]/truemassg[icolor0]
                H.node_0[H.nnodes]['radius']=rsquare[icolor0]
                H.node_0[H.nnodes]['density']=densmoy[icolor0]
                posfin -= np.round(posfin / H.boxarr) * H.boxarr
                H.node_0[H.nnodes]['px']=posfin[0]
                H.node_0[H.nnodes]['py']=posfin[1]
                H.node_0[H.nnodes]['pz']=posfin[2]
                H.node_0[H.nnodes]['rho_saddle']=rhot
                H.node_0[H.nnodes]['level'] = H.node_0[inode]['level']+1
                # if (H.megaverbose and (H.node_0[H.nnodes]['mass']>=H.nmembthresh)):
                #     print(f'{print_prefix}*****************************************')
                #     print(f'{print_prefix}new H.node_0 :',H.nnodes)
                #     print(f'{print_prefix}level    :',H.node_0[H.nnodes]['level'])
                #     print(f'{print_prefix}nsisters :',H.node_0[H.nnodes]['nsisters'])
                #     print(f'{print_prefix}mass     :',H.node_0[H.nnodes]['mass'])
                #     print(f'{print_prefix}true mass:',H.node_0[H.nnodes]['truemass'])
                #     print(f'{print_prefix}radius   :',H.node_0[H.nnodes]['radius'])
                #     print(f'{print_prefix}position : ({H.node_0[H.nnodes]['px']}, {H.node_0[H.nnodes]['py']}, {H.node_0[H.nnodes]['pz']})')
                #     print(f'{print_prefix}rho_saddl:',H.node_0[H.nnodes]['rho_saddle'])
                #     print(f'{print_prefix}rhomax   :',H.node_0[H.nnodes]['densmax'])
                #     print(f'{print_prefix}*****************************************')
        H.node_0[inode]['firstchild']=H.nnodes
        inodeout=inodetmp
        # timer[ith] += time.time()-ref; ith+=1; ref=time.time()
        # # 0.21 sec
        for icolor0 in range(inc_color_tot):
            if (ifok[icolor0]):
                igrAout=igrposnew_0[icolor0]+1
                igrBout=igrposnew_0[icolor0+1]
                for igr1 in frange(igrAout, igrBout):
                    paint_particles_new_131431(idgroup[igr1-1],inodeout,rhot)
                rhotout=rhot*(1+H.fudge/np.sqrt(mass_compg[icolor0]))
                # timer[-1] += time.time()-ref
                # # 5.69 sec
                if (igrBout!=igrAout):
                    create_nodes_13143(rhotout,inodeout,igrAout,igrBout, pool=pool)
                else:
                    H.node_0[inodeout]['firstchild']=0
                inodeout += 1
    elif (nsisters==1):
        inodeout=inode
        rhotout=rhot*(1+H.fudge/np.sqrt(mass_compg[icolor_ref]))
        igrAout=igrposnew_0[0]+1
        igrBout=igrposnew_0[inc_color_tot]
        # timer[ith] += time.time()-ref; ith+=1; ref=time.time()
        # # 0.21 sec
        if (igrBout!=igrAout): create_nodes_13143(rhotout,inodeout,igrAout,igrBout, pool=pool)
        else: H.node_0[inode]['firstchild']=0
    else:
        H.node_0[inode]['firstchild']=0
        # timer[ith] += time.time()-ref; ith+=1; ref=time.time()
        # # 0.21 sec

#=======================================================================
def paint_particles_131431(igroup1,inode,rhot):
#=======================================================================
    ipar1=mem['firstpart_1313'][igroup1-1]
    while(ipar1>0):
        if(mem['density_1312'][ipar1-1]>rhot):
            mem['whereIam_parts'][ipar1-1]=inode
        ipar1=mem['idpart_1311'][ipar1-1]

#=======================================================================
def paint_particles_new_131431(igroup1,inode,rhot):
#=======================================================================
    ipar1s = get_group_members(igroup1)
    denmask = mem['density_1312'][ipar1s-1]>rhot
    mem['whereIam_parts'][ipar1s[denmask]-1] = inode

#=======================================================================
def treat_particles_13141(igroup1,rhot,posg,imass,igroupref,posref,rsquare,densmoy,truemass, strict=False):
#=======================================================================
    boxarr = H.boxarr
    posg[:] = 0; imass=0; 
    rsquare=0; densmoy=0; truemass=0
    
    ipar1=mem['firstpart_1313'][igroup1-1]
    first_good=False
    density = mem['density_1312']
    idpart = mem['idpart_1311']
    firstpart = mem['firstpart_1313']
    pos = mem['pos_10']
    mass = mem['mass_10']
    while (ipar1>0):
        if (density[ipar1-1] > rhot):
            if ( not first_good):
                if (igroupref==0):
                    posref[:3]=pos[ipar1-1,:3]
                    igroupref=igroup1
                first_good=True
                firstpart[igroup1-1]=ipar1
                if strict:
                    densmin=density[ipar1-1]
                    densmax=densmin
            else:
                idpart[iparold1-1]=ipar1
                
            iparold1=ipar1
            imass += 1
            xmasspart = mass[ipar1-1]
            truemass += xmasspart
            dpos = pos[ipar1-1,:3] - posref[:3]
            dpos = np.where(dpos > boxarr/2, boxarr - dpos, dpos) + posref
            posg[:3] += dpos[:3]*xmasspart
            rsquare += xmasspart*np.sum(dpos[:3]**2)
            densmoy += density[ipar1-1]
            if strict:
                densmax = max(densmax,density[ipar1-1])
                densmin = min(densmin,density[ipar1-1])
        ipar1=idpart[ipar1-1]
    if ( not first_good): firstpart[igroup1-1]=0

    if strict:
        if ( (densmin<=rhot)or(densmax!=mem['densityg_1313'][igroup1-1]) ):
            print(f'{print_prefix}ERROR in treat_particles')
            print(f'{print_prefix}igroup1, densmax, rhot=',igroup1,mem['densityg_1313'][igroup1-1],rhot)
            raise ValueError('denslow, denshigh    =',densmin,densmax)
    return posg,imass,igroupref,posref,rsquare,densmoy,truemass

#=======================================================================
# def treat_particles_new_13141(igroup1,rhot,posg,imass,igroupref,posref,rsquare,densmoy,truemass):
def treat_particles_new_13141(igroup1,rhot,igroupref,posref, strict=False):
#=======================================================================
    density = mem['density_1312']
    firstpart = mem['firstpart_1313']
    
    
    
    _ipar1s = get_group_members(igroup1)
    dens = density[_ipar1s-1]
    mask = dens > rhot
    if (mask.any())and(firstpart[igroup1-1]>0):
        idpart = mem['idpart_1311']
        pos = mem['pos_10']
        mass = mem['mass_10']
        dens = dens[mask]
        ipar1s = _ipar1s[mask]
        firstpart[igroup1-1] = ipar1s[-1]
        if igroupref==0:
            igroupref = igroup1
            posref=pos[ipar1s[-1]-1,:3]
        if strict:
            densmin, densmax = np.min(dens), np.max(dens)

        imass = len(ipar1s)
        truemass = np.sum(mass[ipar1s-1])
        dpos = pos[ipar1s-1,:3] - posref[:3]
        boxarr = H.boxarr
        dpos = np.where(dpos > boxarr/2, boxarr - dpos, dpos) + posref
        posg = np.sum(dpos[:,:3]*mass[ipar1s-1][:,None], axis=0)
        rsquare = np.sum(mass[ipar1s-1][:,None]*np.sum(dpos[:,:3]**2, axis=1)[:,None])
        densmoy = np.sum(dens)

        idpart[ipar1s[1:]-1] = ipar1s[:-1] # <--BOOKMARK: not sure if this is correct
    else:
        posg = np.zeros(3)
        imass=0; rsquare=0; densmoy=0; truemass=0
        firstpart[igroup1-1] = 0
    
    if not mask.all():
        set_group_members(_ipar1s[~mask], 0)

    if strict:
        if ( (densmin<=rhot)or(densmax!=mem['densityg_1313'][igroup1-1]) ):
            print(f'{print_prefix}ERROR in treat_particles')
            print(f'{print_prefix}igroup1, densmax, rhot=',igroup1,mem['densityg_1313'][igroup1-1],rhot)
            raise ValueError('denslow, denshigh    =',densmin,densmax)
    return posg,imass,igroupref,posref,rsquare,densmoy,truemass

#=======================================================================
def do_colorize_131430(igroup1,igr1,rhot): #YDdebug
#=======================================================================
    global ncall, icolor_select

    ncall += 1
    print(f"{print_prefix}do_colorize : igroup1={igroup1}, igr1={igr1}, color={icolor_select}, rhot={rhot}, ncall={ncall}")
    mem['color_1314'][igr1-1]=icolor_select
    neig=H.group[igroup1-1].nhnei
    for ineig0 in range(H.group[igroup1-1].nhnei):
        if (H.group[igroup1-1].rho_saddle_gr[ineig0] > rhot):
            # We connect this group to its neighbourg
            igroup2=H.group[igroup1-1].isad_gr[ineig0]
            igrB=mem['igroupid_1314'][igroup2-1]
            if (mem['color_1314'][igrB-1]==0):
                do_colorize_131430(igroup2,igrB,rhot) #YDdebug
            elif (mem['color_1314'][igrB-1]!=icolor_select):
                print(f"{print_prefix}ERROR in do_colorize : color(igrB)({mem['color_1314'][igrB-1]}) <> icolor_select({icolor_select})")
                raise ValueError('The connections are not symmetric.')
            else:
                pass
        else:
            # We do not need this saddle anymore (and use the fact that saddles 
            # are ranked in decreasing order)
            neig -= 1
    H.group[igroup1-1].nhnei=neig


def _compute_saddle_list_worker_legacy(_ipar0, mygids, iinds, icounts):
    _nhneis = np.zeros_like(mygids, dtype=np.int32)
    count = 0
    iparneigh = mem['iparneigh_1312'].view()
    igrouppart = mem['igrouppart_1313'].view()
    for mygroup, iind, icount in zip(mygids, iinds, icounts):
        allneigh = np.unique(iparneigh[:, _ipar0[iind:iind+icount]].flatten())
        theirgroups = np.unique(igrouppart[allneigh - 1])
        ineig = np.sum((theirgroups != mygroup)&(theirgroups > 0))
        _nhneis[count] = ineig
        count += 1
    return _nhneis

def _compute_saddle_list_worker(_ipar0, mygids, iinds, icounts):
    iparneigh = mem['iparneigh_1312'].view()
    igrouppart = mem['igrouppart_1313'].view()
    density = mem['density_1312'].view()
    _nhneis = np.zeros_like(mygids, dtype=np.int32)
    _isad_grs = [None for _ in range(len(mygids))]
    _rho_saddle_grs = [None for _ in range(len(mygids))]
    count = 0
    for mygroup, iind, icount in zip(mygids, iinds, icounts):
        allneigh = iparneigh[:, _ipar0[iind:iind+icount]]
        theirgroups =igrouppart[allneigh - 1]
        theirdens = density[allneigh - 1]
        ourdens = density[_ipar0[iind:iind+icount]]
        theirdens = np.where((theirdens <= ourdens), theirdens, ourdens)
        validmask = (theirgroups != mygroup)&(theirgroups > 0)
        theirdens = theirdens[validmask]
        theirgroups = theirgroups[validmask]
        sort_idx = np.argsort(theirgroups)
        g_sorted = theirgroups[sort_idx]
        d_sorted = theirdens[sort_idx]
        isad_gr, split_indices = np.unique(g_sorted, return_index=True)
        _nhneis[count] = len(isad_gr)
        rho_saddle_gr = np.maximum.reduceat(d_sorted, split_indices)
        if len(isad_gr)>0:
            argsort = np.argsort(-rho_saddle_gr)
            _isad_grs[count] = isad_gr[argsort]
            _rho_saddle_grs[count] = rho_saddle_gr[argsort]
        else:
            _isad_grs[count] = np.array([0], dtype=np.int32)
            _rho_saddle_grs[count] = np.array([0.0], dtype=np.float64)
        count += 1
    return _nhneis, _isad_grs, _rho_saddle_grs

#=======================================================================
def compute_saddle_list_13140():
#=======================================================================
# Compute the lowest density threshold below which each group is 
# connected to an other one
#=======================================================================
    if (H.verbose): print(f'{print_prefix}    Fill the end of the branches of the group tree')

    # Allocate the array of nodes
    H.group = [H.grp() for _ in range(H.ngroups)]

    if (H.verbose): print(f'{print_prefix}    First count the number of neighbourgs of each elementary group...')

    touch = np.zeros(H.ngroups, dtype=np.bool_)#; touch[:H.ngroups] = False
    listg = np.zeros(H.ngroups, dtype=np.int32)#; listg[:H.ngroups] = 0

    firstpart = mem['firstpart_1313']
    igrouppart = mem['igrouppart_1313']
    iparneigh = mem['iparneigh_1312']
    density = mem['density_1312']
    groupmask = igrouppart > 0

    DEBUG=False
    if DEBUG: ref = time.time()
    _ipar0 = np.arange(H.npart)
    argsort = np.argsort(igrouppart[groupmask], kind='stable')
    _igrouppart = igrouppart[groupmask][argsort]
    _ipar0 = _ipar0[groupmask][argsort]
    ind = np.flatnonzero(np.r_[True, _igrouppart[1:] != _igrouppart[:-1]])
    counts = np.diff(np.r_[ind, _igrouppart.size])
    uni = _igrouppart[ind]
    # Single Process
    if H.nbPes==1 or DEBUG:
        for mygroup, iind, icount in zip(uni, ind, counts):
            igroupA0 = mygroup - 1
            ourdens = density[_ipar0[iind:iind+icount]]
            allneigh = iparneigh[:, _ipar0[iind:iind+icount]]
            theirdens = density[allneigh - 1]; theirgroups =igrouppart[allneigh - 1]
            theirdens = np.where((theirdens <= ourdens), theirdens, ourdens)
            validmask = (theirgroups != mygroup)&(theirgroups > 0)
            theirdens = theirdens[validmask]; theirgroups = theirgroups[validmask]
            sort_idx = np.argsort(theirgroups)
            g_sorted = theirgroups[sort_idx]
            d_sorted = theirdens[sort_idx]
            isad_gr, split_indices = np.unique(g_sorted, return_index=True)
            H.group[igroupA0].nhnei = len(isad_gr)
            rho_saddle_gr = np.maximum.reduceat(d_sorted, split_indices)
            if len(rho_saddle_gr)==0:
                H.group[igroupA0].isad_gr = np.array([0], dtype=np.int32)
                H.group[igroupA0].rho_saddle_gr = np.array([0.0], dtype=np.float64)
            else:
                argsort = np.argsort(-rho_saddle_gr)
                H.group[igroupA0].isad_gr = isad_gr[argsort]
                H.group[igroupA0].rho_saddle_gr = rho_saddle_gr[argsort]
    # Multiprocessing
    else:
        nhneis = [0 for _ in range(H.nbPes)]
        isad_grs = [None for _ in range(H.nbPes)]
        rho_saddle_grs = [None for _ in range(H.nbPes)]
        
        def _callback(result, icpu):
            nhneis[icpu] = result[0]
            isad_grs[icpu] = result[1]
            rho_saddle_grs[icpu] = result[2]
        signal.signal(signal.SIGTERM, signal.SIG_DFL)
        with Pool(processes=H.nbPes) as pool:
            results = []
            for idx in range(H.nbPes):
                istart = idx * (len(uni) // H.nbPes)
                iend = (idx + 1) * (len(uni) // H.nbPes) if idx != (H.nbPes-1) else len(uni)
                mygids = uni[istart:iend]
                iinds = ind[istart:iend]
                icounts = counts[istart:iend]
                my_ipar0s = _ipar0[iinds[0] : iinds[-1] + icounts[-1]]
                res = pool.apply_async(
                    _compute_saddle_list_worker,
                    args=(my_ipar0s, mygids, iinds-iinds[0], icounts),
                    callback=lambda r, icpu=idx: _callback(r, icpu))
                results.append(res)
            for res in results:
                res.get()
        signal.signal(signal.SIGTERM, H.flush)

        icpu = 0; i1 = 0; cursor = 0
        _nhneis = nhneis[icpu]
        _isad_grs = isad_grs[icpu]
        _rho_saddle_grs = rho_saddle_grs[icpu]
        for igroupA0 in range(H.ngroups):
            ineig = _nhneis[cursor]
            H.group[igroupA0].nhnei = ineig
            i2 = i1 + ineig
            H.group[igroupA0].rho_saddle_gr = np.array(_rho_saddle_grs[cursor], dtype=np.float64)
            H.group[igroupA0].isad_gr = np.array(_isad_grs[cursor], dtype=np.int32)
            
            i1 = i2
            cursor += 1
            if cursor >= len(_nhneis):
                icpu += 1; i1 = 0; cursor = 0
                if icpu < len(nhneis):
                    _nhneis = nhneis[icpu]
                    _isad_grs = isad_grs[icpu]
                    _rho_saddle_grs = rho_saddle_grs[icpu]

    del touch
    del listg
    

    if (H.verbose): print(f'{print_prefix}    Establish symmetry in connections...')
    whereIam_dict = [{} for _ in range(H.ngroups)]
    for igroupA0 in range(H.ngroups):
        me = H.group[igroupA0]
        if me.nhnei > 0:
            myneighs = me.isad_gr
            for inA0 in range(me.nhnei):
                igroupB1 = myneighs[inA0]
                whereIam_dict[igroupB1-1][igroupA0+1] = inA0
    # Total number of connections count
    icon_count=0

    # Destroy the connections between 2 groups which are not symmetric
    # This might be rather slow and might be discarded later
    idestroy=0

    for igroupA0 in range(H.ngroups):
        me = H.group[igroupA0]
        mynhnei = me.nhnei
        if (mynhnei>0):
            igroupB1s = me.isad_gr
            whereIam = np.fromiter((whereIam_dict[igroupA0].get(igroupB1, -1) for igroupB1 in igroupB1s), dtype=np.int32, count=len(igroupB1s))
            mutualmask = whereIam>=0
            new_nhnei = np.sum(mutualmask)
            new_isad_gr = igroupB1s[mutualmask]
            if mutualmask.any():
                whereIam = whereIam[mutualmask]
                mutuals = igroupB1s[mutualmask]
                my_rho_sad12 = me.rho_saddle_gr[mutualmask]
                their_rho_sad12 = np.array([H.group[igroupB1-1].rho_saddle_gr[iw] for iw,igroupB1 in zip(whereIam,mutuals)], dtype=np.float64)
                TooHigh = my_rho_sad12 >= their_rho_sad12
                # Choose the smaller one
                new_rho_saddle_gr = np.where(TooHigh, their_rho_sad12, my_rho_sad12)
                # Change theirs
                for iw, igroupB1, rsad, win in zip(whereIam, mutuals, new_rho_saddle_gr, ~TooHigh):
                    if win:
                        H.group[igroupB1-1].rho_saddle_gr[iw] = rsad
            else:
                new_isad_gr = np.array([0], dtype=np.int32)
                new_rho_saddle_gr = np.array([0.0], dtype=np.float64)    
            me.tmp = (new_nhnei, new_isad_gr, new_rho_saddle_gr)
            idestroy += mynhnei - new_nhnei
            icon_count += new_nhnei

    for igroupA0 in range(H.ngroups):
        me = H.group[igroupA0]
        if me.nhnei > 0:
            new_nhnei, new_isad_gr, new_rho_saddle_gr = me.tmp
            argsort = np.argsort(-new_rho_saddle_gr)
            me.isad_gr = new_isad_gr[argsort]
            me.rho_saddle_gr = new_rho_saddle_gr[argsort]
            me.nhnei = new_nhnei
            me.tmp = None

    if (H.verbose): print(f'{print_prefix}--> Number of connections removed :',idestroy)
    if (H.verbose): print(f'{print_prefix}--> Total number of connections remaining :',icon_count)
    if (H.verbose): print(f'{print_prefix}Rebuild groups with undesired connections removed...')
    # ---------------------------------------------------------

    H.deallocate('iparneigh_1312')
    # H.deallocate('igrouppart_1313')

#=======================================================================
def compute_density_13121(ipar1,dist2_0,iparnei):
#=======================================================================
    # Vectorized version
    r = np.sqrt(dist2_0[H.nvoisins]) * 0.5
    mneighs = mem['mass_10'][iparnei - 1]
    mmy = mem['mass_10'][ipar1 - 1]
    contrib = np.sum(mneighs * splines(np.sqrt(dist2_0[1:H.nvoisins+1]) / r)) + mmy
    mem['density_1312'][ipar1-1]=(H.xlong*H.ylong*H.zlong)*contrib /(H.pi*r**3)

#=======================================================================
def find_nearest_parts_13120(ipar1, poshere):
#=======================================================================
    dist2_0, iparnei = walk_tree_norec_131200(ipar1,poshere)
    return dist2_0, iparnei

def walk_tree_norec_131200(iparid1,poshere):
    '''
    Non-recursive version of walk_tree_131200.
    '''
    dist2_0 = np.full(H.nvoisins, H.bignum, dtype=np.float64)
    iparnei = np.zeros(H.nvoisins, dtype=np.int32)
    icellidin1 = 1
    # Current Information
    maxdist = dist2_0[H.nvoisins-1]
    boxarr = H.boxarr
    boxarrT = boxarr[:, np.newaxis]
    poshereT = poshere[:, np.newaxis]

    firstchild = mem['firstchild_1311']
    pos_cell = mem['pos_cell_1311']
    size_cell = mem['size_cell_1311']
    sister = mem['sister_1311']
    mass_cell = mem['mass_cell_1311']
    pos = mem['pos_10']
    idpart = mem['idpart_1311']

    stack = [(icellidin1, 0.0)]
    while stack:
        icellid, idiscell = stack.pop()
        if idiscell >= maxdist: continue
        icell_identity1 = icellid
        ifirst = firstchild[icell_identity1-1]
        # Am I leaf or not??
        if ifirst < 0:
            # If still close
            if (idiscell < maxdist):
                
                first_pos_this_node=-ifirst-1 # First particle index
                mynpart = mass_cell[icell_identity1-1]
                myparts = idpart[first_pos_this_node : first_pos_this_node+mynpart]
                if iparid1 in myparts:
                    myparts = myparts[myparts != iparid1]
                    mynpart -= 1
                if mynpart==0: continue
                dpos = np.abs(pos[myparts-1,:] - poshere)
                dpos = np.where(dpos > boxarr/2, boxarr - dpos, dpos)
                distance2p = np.sum(dpos*dpos, axis=1)
                dmask = distance2p < maxdist
                if np.any(dmask):
                    # Merge and sort `dist2_0`
                    tmp = np.concatenate( (dist2_0, distance2p[dmask]) )
                    partitioned_indices = np.argpartition(tmp, H.nvoisins)[:H.nvoisins]
                    ordk = np.argsort(tmp[partitioned_indices])
                    partitioned_indices = partitioned_indices[ordk]
                    dist2_0[:] = tmp[partitioned_indices]; maxdist = dist2_0[H.nvoisins-1]
                    # Merge and sort `iparnei`
                    tmp = np.concatenate( (iparnei, myparts[dmask]) )
                    iparnei = tmp[partitioned_indices]
        else:
            icell_identity1 = ifirst
            sc = size_cell[icell_identity1-1]
            # STEP1) Split Child cells
            icell_identity1s = np.empty(8, dtype=np.int32)
            for inc in range(8):
                icell_identity1s[inc] = icell_identity1
                icell_identity1=sister[icell_identity1-1]
                if (icell_identity1 == 0):
                    break
            inc += 2
            if inc<9: icell_identity1s = icell_identity1s[:inc-1]
            
            dpos = np.abs(pos_cell[:, icell_identity1s-1] - poshereT)
            dpos = np.where(dpos > boxarrT/2, boxarrT - dpos, dpos)
            dpos -= sc
            dpos[dpos < 0] = 0.0  # <--- Inside the cell
            distance2s = np.sum(dpos*dpos, axis=0)
            
            sel = np.flatnonzero(maxdist > distance2s)
            if sel.size==0:
                continue
            inc = sel.size+1
            distance2s = distance2s[sel]
            icell_identity1s = icell_identity1s[sel]
            argsort = np.argsort(distance2s)
            discell2_0 = distance2s[argsort]
            icid = icell_identity1s[argsort]  # Cell IDs

            for i in range(inc-2, -1, -1):
                idiscell = discell2_0[i]
                if idiscell < maxdist:
                    icellid_out1 = icid[i]
                    stack.append((icellid_out1, idiscell))
    return np.insert(dist2_0,0,0), iparnei


#=======================================================================
def walk_tree_131200(iparid1,icellidin1,poshere,dist2_0,iparnei, inccellpart):
#=======================================================================
    '''
    The subroutine 
    calculates the distances between the particle and cells (or the particles within cells),
    stores the closest cells/particles in `discell2_0` and `icid` arrays, and
    updates `dist2_0` and `iparnei` arrays with the closest particles.
    The subroutine recursively calls itself with the children cells of the current cell until it reaches the leaf cells that contain particles.
    The algorithm is an implementation of a tree structure, where cells are nodes and the particles are leaves.
    
    - by ChatGPT Jan 9
    '''

    maxdist = dist2_0[H.nvoisins]
    icell_identity1 = mem['firstchild_1311'][icellidin1-1]
    sc = mem['size_cell_1311'][icell_identity1-1]  # Distance to cell face
    # lvl = np.log2(H.xlong/mem['size_cell_1311'][icellidin1-1]/2)
    # Until icell_identity1==0: (Final leaf of tree)
    # Calc distance (poshere <-> cells)
    sister = mem['sister_1311']
    boxarr = H.boxarr
    boxarrT = boxarr[:, np.newaxis]

    icell_identity1s = np.zeros(8, dtype=np.int32)
    for inc in range(8):
        icell_identity1s[inc] = icell_identity1
        icell_identity1=sister[icell_identity1-1]
        if (icell_identity1 == 0):
            break
    inc += 2
    if inc<9: icell_identity1s = icell_identity1s[:inc-1]

    poshereT = poshere[:, np.newaxis]
    dpos = np.abs(mem['pos_cell_1311'][:, icell_identity1s-1] - poshereT)
    dpos = np.where(dpos > boxarrT/2, boxarrT - dpos, dpos)
    # sc = mem['size_cell_1311'][icell_identity1s-1]  # Distance to cell face
    dpos -= sc
    dpos[dpos < 0] = 0.0  # <--- Inside the cell
    distance2s = np.sum(dpos**2, axis=0)
    
    mask = maxdist > distance2s
    if not np.any(mask):
        return dist2_0, iparnei
    inc = np.sum(mask)+1
    inccellpart += inc - 1
    distance2s = distance2s[mask]
    icell_identity1s = icell_identity1s[mask]
    argsort = np.argsort(distance2s)
    discell2_0 = np.empty(distance2s.size + 1, dtype=np.float64)
    discell2_0[0] = 0.0; discell2_0[1:] = distance2s[argsort]
    icid = icell_identity1s[argsort]  # Cell IDs

    # Loop for counted cells,
    # Update the closest particle ID(in iparnei) and the distance to that part(in dist2_0)
    pos = mem['pos_10']
    mass_cell = mem['mass_cell_1311']
    idpart = mem['idpart_1311']
    for ic0 in range(inc-1):
        icellid_out1=icid[ic0]
        # lvl = np.log2(H.xlong/sc/2)
        # If this cell is a leaf, calc dists from inside particles
        if (mem['firstchild_1311'][icellid_out1-1] < 0):
            # If still close
            if (discell2_0[ic0+1]<maxdist):
                # print(icellid_out1, discell2_0[ic0+1], maxdist)
                first_pos_this_node=-mem['firstchild_1311'][icellid_out1-1]-1 # First particle index
                mynpart = mass_cell[icellid_out1-1]
                # ----------------------------------------------------
                # [Development - vectorized version]
                myparts = idpart[first_pos_this_node : first_pos_this_node+mynpart]
                if iparid1 in myparts:
                    myparts = myparts[myparts != iparid1]
                    mynpart -= 1
                if mynpart==0: return dist2_0, iparnei
                dpos = np.abs(pos[myparts-1,:] - poshere)
                dpos = np.where(dpos > boxarr/2, boxarr - dpos, dpos)
                distance2p = np.sum(dpos**2, axis=1)
                if distance2p.min() < maxdist:
                    # Merge and sort `dist2_0`
                    tmp = np.concatenate( (dist2_0[1:], distance2p) )
                    partitioned_indices = np.argpartition(tmp, H.nvoisins)[:H.nvoisins]
                    ordk = np.argsort(tmp[partitioned_indices])
                    partitioned_indices = partitioned_indices[ordk]
                    tmp = tmp[partitioned_indices]
                    dist2_0 = np.empty(H.nvoisins+1, dtype=np.float64)
                    dist2_0[0] = 0.0; dist2_0[1:] = tmp; maxdist = dist2_0[H.nvoisins]
                    # Merge and sort `iparnei`
                    tmp = np.concatenate( (iparnei, myparts) )
                    iparnei = tmp[partitioned_indices]
                # ----------------------------------------------------
        # If (not a leaf) and (still close)
        elif (discell2_0[ic0+1] < maxdist):
            # walk_tree(icellid_out1,poshere,dist2_0,iparid1,inccellpart,iparnei)
            dist2_0, iparnei = walk_tree_131200(iparid1,icellid_out1,poshere,dist2_0,iparnei, inccellpart)
            maxdist = dist2_0[H.nvoisins]
        else: pass
    return dist2_0, iparnei

#=======================================================================
def create_tree_structure_1311_scipy():
#=======================================================================
    from scipy.spatial import KDTree

    H.allocate('idpart_1311',H.npart, dtype=np.int32)
    if (H.verbose): print(f'{print_prefix}Create KDtree structure...')
    if (H.verbose): ref = time.time()
    tree = KDTree(mem['pos_10'], 
                  leafsize=H.npartpercell, compact_nodes=True,
                  copy_data=False, balanced_tree=True,
                  boxsize=H.xlong)
    if (H.verbose): print(f'{print_prefix}--> total number of cells =',tree.size)
    if (H.verbose): print(f"{print_prefix}--> {time.time()-ref:.2f} seconds to create KDtree structure")
    return tree



#=======================================================================
def create_tree_structure_1311():
#=======================================================================
    global inccell, timers
    pos_this_node = np.empty(3, dtype=np.float64)

    if (H.verbose): print(f'{print_prefix}Create tree structure...')

    # we modified to put 2*H.npart-1 instead of 2*H.npart so that AdaptaHOP can work on a 1024^3, 2*(1024^3)-1 is still an integer(kind=4), 2*(1024^3) is not 
    H.ncellmx=int( (2*H.npart -1)/H.nvoisins*4 )
    H.ncellbuffer=max(round(0.1*H.npart),H.ncellbuffermin)
    H.allocate('idpart_1311',H.npart, dtype=np.int32)
    H.allocate('mass_cell_1311',H.ncellmx, dtype=np.int32)
    H.allocate('size_cell_1311',H.ncellmx, dtype=np.float64)
    H.allocate('pos_cell_1311',(3,H.ncellmx), dtype=np.float64)
    H.allocate('sister_1311',H.ncellmx, dtype=np.int32)
    H.allocate('firstchild_1311',H.ncellmx, dtype=np.int32)

    # -----------------------------------------------------
    # [ITERATIVE MODE]
    mem['idpart_1311'][:] = np.arange(H.npart, dtype=np.int32)+1

    nlevel=0
    inccell=0
    idmother=0
    pos_this_node[:]=np.float64(0.)
    npart_this_node=H.npart
    first_pos_this_node=0
    mem['pos_cell_1311'][:]=0
    mem['size_cell_1311'][:]=0
    mem['mass_cell_1311'][:]=0
    mem['sister_1311'][:]=0
    mem['firstchild_1311'][:]=0
    H.sizeroot = np.float64(max(H.xlong, H.ylong, H.zlong))

    ref = time.time()

    # -------------------------------------------------------
    # Calculate root cells (<= baselevel) first
    # -------------------------------------------------------
    baselevel = min(int(np.ceil(np.log2(H.nbPes)/3)), 3) + 1
    tmp = 8**np.arange(baselevel+1)
    nchunk = 8**baselevel
    nrootcells = np.sum(tmp)
    if (H.verbose): print(f"{print_prefix}    Base level = {baselevel} ({nchunk} cells) for {H.nbPes} PEs")
    H.allocate('rmass_cell_1311',nrootcells, dtype=np.int32)
    H.allocate('rsize_cell_1311',nrootcells, dtype=np.float64)
    H.allocate('rpos_cell_1311',(3,nrootcells), dtype=np.float64)
    H.allocate('rsister_1311',nrootcells, dtype=np.int32)
    H.allocate('rfirstchild_1311',nrootcells, dtype=np.int32)
    mem['rpos_cell_1311'][:]=0
    mem['rsize_cell_1311'][:]=0
    mem['rmass_cell_1311'][:]=0
    mem['rsister_1311'][:]=0
    mem['rfirstchild_1311'][:]=0
    stack_chunks, npartchunks = create_KDtree_root(baselevel=baselevel)
    if (H.verbose): print(f"{print_prefix}    Root builds {np.sum(tmp)} cells ({[int(t) for t in tmp]})")
    if (H.verbose): print(f"{print_prefix}    Nparts per chunk: {npartchunks.min()}-{npartchunks.max()} (total {npartchunks.sum()})")

    nroot_counts = np.zeros(nchunk, dtype=np.int32)
    for ilvl in range(baselevel+1):
        tmp = np.arange(0, nchunk, 8**(baselevel-ilvl))
        nroot_counts[tmp] += 1
    nroot_cumsum = np.cumsum(nroot_counts)

    # -------------------------------------------------------
    # Calculate sons of root cells in parallel
    # -------------------------------------------------------
    H.allocate('ncell_chunks_1311',nchunk, dtype=np.int32)
    if H.nbPes == 1: # Single process
        if (H.verbose): iterator = tqdm(range(nchunk), desc=f"{print_prefix}    [Ncpu=1] Creating KD-tree chunks", ncols=100)
        else: iterator = range(nchunk)
        for ichunk in iterator:
            _create_KDtree_worker(ichunk, stack_chunks[ichunk], inccell)
        if (H.verbose): iterator.close()
    else: # Multi processes
        if (H.verbose):
            pbar = tqdm(total=nchunk, desc=f"{print_prefix}    [Ncpu={H.nbPes}] Creating KD-tree chunks", ncols=100)
            callback = lambda _: pbar.update()
        else:
            pbar = None
            callback = None
        signal.signal(signal.SIGTERM, signal.SIG_DFL)
        with Pool(processes=H.nbPes) as pool:
            async_results = []
            for icpu1 in range(nchunk):
                r = pool.apply_async(_create_KDtree_worker, args=(icpu1, stack_chunks[icpu1], inccell), callback=callback)
                async_results.append(r)
            for r in async_results:
                r.get()
        signal.signal(signal.SIGTERM, H.flush)
        if (H.verbose): pbar.close()
        
    # -------------------------------------------------------
    # Rearrange
    # -------------------------------------------------------
    ncell_chunks = mem['ncell_chunks_1311']
    inccell += np.sum(ncell_chunks)
    ncell_cumul = np.cumsum(ncell_chunks)
    _ncell_cumul = np.empty(nchunk+1, dtype=np.int32)
    _ncell_cumul[0] = 0; _ncell_cumul[1:] = ncell_cumul
    insert_starts = _ncell_cumul[:-1] + nroot_cumsum
    _nroot_cumsum = np.empty(nchunk+1, dtype=np.int32)
    _nroot_cumsum[0] = 0; _nroot_cumsum[1:] = nroot_cumsum
    insert_ends   = _nroot_cumsum[1:]  + ncell_cumul
    nroot_remains = nrootcells - nroot_cumsum
    root_index = np.arange(nrootcells) + np.repeat(_ncell_cumul[:-1], nroot_counts)
    

    # -------------------------------------------------------
    # Gathering
    # -------------------------------------------------------
    cellkeys = [('mass_cell', 'i4'), ('size_cell', 'f8'),
                ('pos_cell', 'f8'), ('sister', 'i4'),
                ('firstchild', 'i4')]
    if H.nbPes == 1: # Single process
        if (H.verbose): 
            iterator = tqdm(range(nchunk), desc=f"{print_prefix}    [Ncpu=1] Gathering KD-tree chunks", ncols=100)
        else:
            iterator = range(nchunk)
        for ichunk in iterator:
            willbe_shifted = nroot_cumsum[ichunk]
            nroot_remain = nroot_remains[ichunk]
            insert_start = insert_starts[ichunk]
            insert_end   = insert_ends[ichunk]
            ncell_ichunk = ncell_chunks[ichunk]
            ncell_allchunk = ncell_cumul[ichunk]
            iroot1 = _nroot_cumsum[ichunk]
            iroot2 = _nroot_cumsum[ichunk+1]
            _gather_KDtree_worker(
            ichunk,
            willbe_shifted, nroot_remain, insert_start, insert_end,
            ncell_ichunk, ncell_allchunk, iroot1, iroot2,
            root_index, cellkeys, nroot_counts
            )
        if (H.verbose): iterator.close()
    else: # Multi processes
        if (H.verbose): 
            pbar = tqdm(total=nchunk, desc=f"{print_prefix}    [Ncpu={H.nbPes}] Gathering KD-tree chunks", ncols=100)
            callback = lambda _: pbar.update()
        else:
            pbar = None
            callback = None
        signal.signal(signal.SIGTERM, signal.SIG_DFL)
        with Pool(processes=H.nbPes) as pool:
            async_results = []
            for ichunk in range(nchunk):
                willbe_shifted = nroot_cumsum[ichunk]
                nroot_remain = nroot_remains[ichunk]
                insert_start = insert_starts[ichunk]
                insert_end   = insert_ends[ichunk]
                ncell_ichunk = ncell_chunks[ichunk]
                ncell_allchunk = ncell_cumul[ichunk]
                iroot1 = _nroot_cumsum[ichunk]
                iroot2 = _nroot_cumsum[ichunk+1]
                r = pool.apply_async(
                    _gather_KDtree_worker,
                    args=(
                        ichunk,
                        willbe_shifted, nroot_remain, insert_start, insert_end,
                        ncell_ichunk, ncell_allchunk, iroot1, iroot2,
                        root_index, cellkeys, nroot_counts
                    ),
                    callback=callback
                )
                async_results.append(r)
            for r in async_results:
                r.get()
        signal.signal(signal.SIGTERM, H.flush)
        if (H.verbose): pbar.close()

    # Change root (0 <= lvl < baselvl) firstchild
    firstchild = mem["firstchild_1311"]
    sizezero = H.sizeroot/2
    lvls = np.int32(np.log2(sizezero / mem['size_cell_1311'][root_index]))
    mask = lvls < baselevel
    need_to_change = root_index[mask]
    for idx in need_to_change:
        old_fc = firstchild[idx]
        new_fc = root_index[old_fc -1]+1  # -1 for Fortran->C index
        firstchild[idx] = new_fc

    # Change sister sign
    sister_1311 = mem["sister_1311"]
    sister_1311 = np.where(sister_1311 < 0, -sister_1311, sister_1311)
    mem["sister_1311"] = sister_1311

    ncell = inccell
    if (H.verbose): print(f'{print_prefix}--> total number of cells =',ncell)
    if (H.verbose): print(f"{print_prefix}--> {time.time()-ref:.2f} seconds to create the tree structure")

    H.deallocate(
        'rmass_cell_1311', 'rsize_cell_1311', 'rpos_cell_1311',
        'rsister_1311', 'rfirstchild_1311', 'ncell_chunks_1311',
        )
    
def _create_KDtree_worker(icpu, stack, inccell_root):
    pos_ref_0 = H.pos_ref_0

    idpart_1311      = maccess('idpart_1311')#mem['idpart_1311']
    pos_10           = maccess('pos_10')
    
    ncellmx=int( (2*H.npart -1)/H.nvoisins )
    pos_cell   = np.empty((3,ncellmx), dtype=np.float64)
    mass_cell  = np.empty(ncellmx, dtype=np.int32)
    size_cell  = np.empty(ncellmx, dtype=np.float64)
    sister     = np.zeros(ncellmx, dtype=np.int32)
    firstchild = np.zeros(ncellmx, dtype=np.int32)

    pos_cell[:, :inccell_root] = mem['rpos_cell_1311']
    mass_cell[:inccell_root]   = mem['rmass_cell_1311']
    size_cell[:inccell_root]   = mem['rsize_cell_1311']
    sister[:inccell_root]      = mem['rsister_1311']
    firstchild[:inccell_root]  = mem['rfirstchild_1311']

    sizeroot     = H.sizeroot
    npartpercell = H.npartpercell
    nlevelmax    = H.nlevelmax

    _inccell = inccell_root
    while stack:
        nlevel, pos_this_node, npart_this_node, first_pos_this_node, idmother = stack.pop()

        # Stop refinement due to no particle in leaf cell
        if npart_this_node <= 0:
            continue

        # ---------------- Enter node / allocate cell id ----------------
        _inccell += 1
        # ---------------- Increase Array -------------------
        if _inccell > ncellmx:
            ncellmx_old = ncellmx
            ncellmx += H.ncellbuffer
            
            new_shape = (ncellmx,)
            tmp = np.zeros(ncellmx_old, dtype=np.int32); old_shape = tmp.shape
            tmp[:] = mass_cell[:]
            mass_cell = np.empty(new_shape, dtype=np.int32)
            mass_cell[:old_shape[0]] = tmp[:]
            tmp[:] = sister[:]
            sister = np.zeros(new_shape, dtype=np.int32)
            sister[:old_shape[0]] = tmp[:]
            tmp[:] = firstchild[:]
            firstchild = np.zeros(new_shape, dtype=np.int32)
            firstchild[:old_shape[0]] = tmp[:]
            del tmp

            tmp = np.zeros(ncellmx_old, dtype=np.float64)
            size_cell = np.empty(new_shape, dtype=np.float64)
            size_cell[:old_shape[0]] = tmp[:]
            del tmp

            tmp = np.zeros((3, ncellmx_old), dtype=np.float64); old_shape = tmp.shape
            pos_cell = np.empty((3, ncellmx), dtype=np.float64)
            pos_cell[:, :old_shape[1]] = tmp[:, :]
            del tmp
        # ---------------------------------------------------

        # meta write
        pos_cell[:3, _inccell-1] = pos_this_node
        mass_cell[_inccell-1] = npart_this_node
        size_cell[_inccell-1] = (2.0 ** (-nlevel)) * sizeroot * 0.5

        # link to mother
        if idmother > 0:
            sister[_inccell-1] = firstchild[idmother-1]
            firstchild[idmother-1] = _inccell # <--------could be race condition
        else:
            pass

        # leaf condition
        if (npart_this_node <= npartpercell) or (nlevel == nlevelmax):
            firstchild[_inccell-1] = -(first_pos_this_node + 1)
            continue

        # ---------------- Split / sort within this node ----------------
        # Count the number of particles in each subcell of this node
        idpart1s = idpart_1311[first_pos_this_node:first_pos_this_node+npart_this_node].view()

        # gather positions then compute icids
        tmp = np.take(pos_10, idpart1s - 1, axis=0)
        icids = icellids(tmp, pos_this_node=pos_this_node, mode=1)
        incsubcell_0 = np.bincount(icids, minlength=8).astype(np.int32)
        nsubcell_0 = np.empty(8, dtype=np.int32)
        nsubcell_0[0] = 0
        nsubcell_0[1:] = np.cumsum(incsubcell_0[:-1])
        argsort = np.argsort(icids, kind='mergesort')

        idpart_1311[first_pos_this_node:first_pos_this_node+npart_this_node] = idpart1s[argsort]

        # ---------------- Push children (reverse order for same traversal) ----------------
        idmother_out = _inccell
        scale = sizeroot * (2.0 ** (-nlevel - 2))

        # push j0=7..0 so that pop() processes 0..7 like recursive loop
        for j0 in range(7, -1, -1):
            npart_this_node_out = int(incsubcell_0[j0])
            if npart_this_node_out <= 0:
                continue
            first_pos_this_node_out = int(first_pos_this_node + nsubcell_0[j0])
            pos_this_node_out = pos_this_node[:] + (pos_ref_0[:, j0] * scale)
            stack.append((
                nlevel + 1,
                pos_this_node_out,
                npart_this_node_out,
                first_pos_this_node_out,
                idmother_out,
            ))

    # Dump
    datdump(pos_cell[:, inccell_root:_inccell].flatten(), f'pos_cell_i{icpu:04d}.tmp')
    datdump(mass_cell[inccell_root:_inccell], f'mass_cell_i{icpu:04d}.tmp')
    datdump(size_cell[inccell_root:_inccell], f'size_cell_i{icpu:04d}.tmp')
    datdump(sister[inccell_root:_inccell], f'sister_i{icpu:04d}.tmp')
    datdump(firstchild[inccell_root:_inccell], f'firstchild_i{icpu:04d}.tmp')
    datdump(firstchild[:inccell_root], f'firstchild_r{icpu:04d}.tmp')

    ncell_ichunk = _inccell - inccell_root
    ncell_chunks_1311 = maccess('ncell_chunks_1311')
    ncell_chunks_1311[icpu] = ncell_ichunk

def _gather_KDtree_worker(
        # variables
        ichunk,
        willbe_shifted, nroot_remain, insert_start, insert_end,
        ncell_ichunk, ncell_allchunk, iroot1, iroot2,
        # static
        root_index, cellkeys, nroot_counts
        ):
    for key, dtype in cellkeys:
        # Shift root cells      
        _root = mem[f"r{key}_1311"][...,iroot1:iroot2]
        if ('sister' in key):
            _rootold = _root.copy()
            _root = np.where(_rootold > 0, root_index[_rootold-1]+1, _rootold)
        mem[f"{key}_1311"][..., insert_start-nroot_counts[ichunk]:insert_start] = _root

        # Insert chunk cells
        fname = f"{key}_i{ichunk:04d}.tmp"
        arr = datload(fname, dtype=dtype)
        if 'pos' in key:
            arr = arr.reshape(3,-1)
        if key == 'firstchild':
            # Correct firstchild indices of chunks
            arr = arr.copy()
            arr[arr > 0] += ncell_allchunk - ncell_ichunk - nroot_remain
        if key=='sister':
            # Correct sister indices of chunks
            arr = arr.copy()
            arr[arr > nroot_remain] += ncell_allchunk - ncell_ichunk - nroot_remain
            arr *= -1
        mem[f"{key}_1311"][..., insert_start:insert_end] = arr
        os.remove(fname)

    # Modify firstchild of {ichunk}th Root (lvl=baselevel)
    rname = f"firstchild_r{ichunk:04d}.tmp"
    rarr = datload(f"firstchild_r{ichunk:04d}.tmp", dtype=dtype)
    our_mothder = rarr[willbe_shifted-1]
    mem[f"firstchild_1311"][root_index[willbe_shifted-1]] = our_mothder + ncell_allchunk - ncell_ichunk - nroot_remain
    os.remove(rname)

def create_KDtree_root(baselevel=0):
    global inccell
    stacks = []
    chunks = [[] for _ in range(8**baselevel)]
    npartchunks = [np.int32(0) for _ in range(8**baselevel)]
    ichunk = 0

    pos_ref_0 = H.pos_ref_0

    # Level 0 Root cell
    nlevel=0
    idmother=0
    pos_this_node=np.zeros(3, dtype=np.float64)
    npart_this_node=H.npart
    first_pos_this_node=0

    # local refs to avoid repeated dict lookups
    idpart_1311      = mem['idpart_1311']
    pos_10           = mem['pos_10']
    pos_cell_1311    = mem['rpos_cell_1311']
    mass_cell_1311   = mem['rmass_cell_1311']
    size_cell_1311   = mem['rsize_cell_1311']
    sister_1311      = mem['rsister_1311']
    firstchild_1311  = mem['rfirstchild_1311']
    
    sizeroot     = H.sizeroot
    npartpercell = H.npartpercell
    nlevelmax    = H.nlevelmax

    stacks.append( (nlevel, pos_this_node.copy(), npart_this_node, first_pos_this_node, idmother) )

    tmp = 8**np.arange(baselevel+1)
    nrootcells = np.sum(tmp)
    if (H.verbose): 
        pbar = tqdm(total=nrootcells, desc=f"{print_prefix}    Creating KD-tree root", ncols=100)

    # for ilvl in range(baselevel+1):
    while stacks:
        nlevel, pos_this_node, npart_this_node, first_pos_this_node, idmother = stacks.pop()
        ilvl = nlevel
        inccell += 1
        pos_cell_1311[:3, inccell-1] = pos_this_node
        mass_cell_1311[inccell-1] = npart_this_node
        size_cell_1311[inccell-1] = (2.0 ** (-nlevel)) * sizeroot * 0.5
    
        # link to mother
        if idmother > 0:
            sister_1311[inccell-1] = firstchild_1311[idmother-1]
            firstchild_1311[idmother-1] = inccell
        else:
            pass

        # leaf condition
        if (npart_this_node <= npartpercell):
            firstchild_1311[inccell-1] = -(first_pos_this_node + 1)
            continue

        # ---------------- Split / sort within this node ----------------
        # Count the number of particles in each subcell of this node
        idpart1s = idpart_1311[first_pos_this_node:first_pos_this_node+npart_this_node].view()

        # gather positions then compute icids
        icids = icellids(np.take(pos_10, idpart1s - 1, axis=0), pos_this_node=pos_this_node, mode=1)
        argsort = np.argsort(icids, kind='mergesort')
        incsubcell_0 = np.bincount(icids, minlength=8).astype(np.int32)
        nsubcell_0 = np.empty(8, dtype=np.int32); nsubcell_0[0] = 0
        nsubcell_0[1:] = np.cumsum(incsubcell_0[:-1])
        
        idpart_1311[first_pos_this_node:first_pos_this_node+npart_this_node] = idpart1s[argsort]

        # ---------------- Push children (reverse order for same traversal) ----------------
        idmother_out = inccell
        scale = sizeroot * (2.0 ** (-nlevel - 2))

        # push j0=7..0 so that pop() processes 0..7 like recursive loop
        for j0 in range(7, -1, -1):
            npart_this_node_out = int(incsubcell_0[j0])
            if npart_this_node_out <= 0:
                continue
            first_pos_this_node_out = int(first_pos_this_node + nsubcell_0[j0])
            pos_this_node_out = pos_this_node[:3] + (pos_ref_0[:3, j0] * scale)
            if ilvl < baselevel:
                stacks.append((
                    nlevel + 1,
                    pos_this_node_out,
                    npart_this_node_out,
                    first_pos_this_node_out,
                    idmother_out,
                ))
            else:
                chunks[ichunk].append((
                    nlevel + 1,
                    pos_this_node_out,
                    npart_this_node_out,
                    first_pos_this_node_out,
                    idmother_out,
                ))
        if ilvl == baselevel:
            npartchunks[ichunk] = np.sum(incsubcell_0)
            ichunk += 1
        if (H.verbose): pbar.update(1)
                
    if (H.verbose): pbar.close()
    return chunks, np.asarray(npartchunks, dtype=np.int32)


#=======================================================================
def create_KDtree_13110(
        nlevel:np.int32,
        pos_this_node:np.ndarray[np.float64],
        npart_this_node, 
        first_pos_this_node,
        idmother,
        pos_ref_0=None):
#=======================================================================
#  nlevel : level of the H.node_0 in the octree. Level zero corresponds to 
#           the full box
#  pos_this_node : position of the center of this H.node_0
#  H.npart  : total number of particles 
#  idpart : array of dimension H.npart containing the id of each
#           particle. It is sorted such that neighboring particles in
#           this array belong to the same cell H.node_0.
#  npart_this_node : number of particles in the considered H.node_0
#  first_pos_this_node : first position in idpart of the particles 
#           belonging to this H.node_0
#  pos :    array of dimension 3.H.npart giving the positions of 
#           each particle belonging to the halo
#  inccell : cell id number for the newly created structured grid site
#  pos_cell : array of dimension 3.H.npart (at most) corresponding
#           to the positions of each H.node_0 of the structured grid
#  mass_cell : array of dimension H.npart (at most) giving the 
#           number of particles in each H.node_0 of the structured grid
#  size_cell : array of dimension H.npart (at most) giving half the
#           size of the cube forming each H.node_0 of the structured grid
#  H.sizeroot : size of the root cell (nlevel=0)
#  idmother : id of the mother cell
#  sister   : sister of a cell (at the same level, with the same mother)
#  firstchild : first child of a cell (then the next ones are found 
#             with the array sister). If it is a cell containing only
#           one particle, it gives the id of the particle.
#  ncellmx : maximum number of cells
#  H.megaverbose : detailed H.verbose mode
#=======================================================================
    ref = time.time()
    global inccell, timers, switch
    
    if(pos_ref_0 is None):
        pos_ref_0 = np.array([-1., -1., -1.,
                             1., -1., -1.,
                            -1.,  1., -1.,
                             1.,  1., -1.,
                            -1., -1.,  1.,
                             1., -1.,  1.,
                            -1.,  1.,  1.,
                             1.,  1.,  1.], dtype=np.float64)
        pos_ref_0 = pos_ref_0.reshape((3, 8), order='F')  # 'F' for Fortran-style ordering
        timers[-1] += time.time() - ref; ref = time.time()

    header=0; stage=0
    if (npart_this_node>0):
        inccell += 1
        if ( ((inccell%1000000)==0)and(H.megaverbose) ):
            print(f'{print_prefix}    inccell=',inccell)
        timers[-1] += time.time() - ref; ref = time.time()
        
        # ---------------- Increase Array -------------------
        if (inccell>H.ncellmx):
            # If we have reached the maximum number of cells, we increase
            # the size of the arrays and reallocate them
            ncellmx_old=H.ncellmx
            H.ncellmx += H.ncellbuffer
            if(H.megaverbose):
                print(f'{print_prefix}ncellmx({ncellmx_old}) is too small. Increase(+{H.ncellbuffer}) it and reallocate arrays accordingly')
            tmp = np.zeros(ncellmx_old, dtype=np.int32)
            H.adjust_size('mass_cell_1311', H.ncellmx, np.int32, tmp=tmp)
            H.adjust_size('sister_1311', H.ncellmx, np.int32, tmp=tmp)
            H.adjust_size('firstchild_1311', H.ncellmx, np.int32, tmp=tmp, init=0)
            del tmp
            tmp = np.zeros(ncellmx_old, dtype=np.float64)
            H.adjust_size('size_cell_1311', H.ncellmx, np.float64, tmp=tmp)
            del tmp
            pos_cell_tmp = np.zeros((3,ncellmx_old), dtype=np.float64)
            H.adjust_size('pos_cell_1311', (3,H.ncellmx), np.float64, tmp=pos_cell_tmp)
            del pos_cell_tmp
            timers[-1] += time.time() - ref; ref = time.time()
        # ---------------------------------------------------
        
        ref = time.time()
        mem['pos_cell_1311'][:3,inccell-1]=pos_this_node
        timers[header+stage] += time.time() - ref; ref = time.time(); header+=0; stage+=0.1
        mem['mass_cell_1311'][inccell-1]=npart_this_node
        timers[header+stage] += time.time() - ref; ref = time.time(); header+=0; stage+=0.1
        mem['size_cell_1311'][inccell-1]=2.**(-nlevel)*H.sizeroot*0.5
        timers[header+stage] += time.time() - ref; ref = time.time(); header+=1; stage=0 # 0->1
        if (idmother>0):
            # If this is not the root cell, we link it to its mother
            mem['sister_1311'][inccell-1]=mem['firstchild_1311'][idmother-1]
            timers[header+stage] += time.time() - ref; ref = time.time(); header+=0; stage+=0.1
            mem['firstchild_1311'][idmother-1]=inccell
            timers[header+stage] += time.time() - ref; ref = time.time(); header+=1; stage=0 # 1->2
        else:
            header += 1; stage=0
        if ((npart_this_node <= H.npartpercell) or (nlevel==H.nlevelmax)):
            # If there is only `H.npartpercell` particles in the `H.node_0` or we have reach
            # maximum level of refinement, we are done
            mem['firstchild_1311'][inccell-1]=-(first_pos_this_node+1)
            timers[-1] += time.time() - ref; ref = time.time() # 1->2
            return
    else:
        # Stop refinement due to no particle in leaf cell
        return
   
    #  Count the number of particles in each subcell of this H.node_0
    idpart1s = mem['idpart_1311'][first_pos_this_node : first_pos_this_node+npart_this_node].view()
    timers[header+stage] += time.time() - ref; ref = time.time(); header+=1; stage=0 # 2->3

    # tmp = mem['pos_10'][idpart1s-1]
    tmp = np.take(mem['pos_10'], idpart1s-1, axis=0)
    timers[header+stage] += time.time() - ref; ref = time.time(); header+=0; stage+=0.1
    icids = icellids(tmp, pos_this_node=pos_this_node, mode=1)
    timers[header+stage] += time.time() - ref; ref = time.time(); header+=1; stage=0 # 3->4

    incsubcell_0 = np.bincount(icids, minlength=8).astype(np.int32)
    timers[header+stage] += time.time() - ref; ref = time.time(); header+=1; stage=0 # 4->5

    #  Create the array of positions of the first particle of the lists
    #  of particles belonging to each subnode
    nsubcell_0 = np.empty(8, dtype=np.int32)
    nsubcell_0[0] = 0
    nsubcell_0[1:] = np.cumsum(incsubcell_0[:-1])
    timers[header+stage] += time.time() - ref; ref = time.time(); header+=1; stage=0 # 5->6

    #  Sort the array of ids (idpart) to gather the particles belonging
    #  to the same subnode. Put the result in `idpart_tmp`.
    argsort = np.argsort(icids, kind='mergesort')
    timers[header+stage] += time.time() - ref; ref = time.time(); header+=1; stage=0 # 6->7

    #  Put back the sorted ids in idpart
    mem['idpart_1311'][first_pos_this_node : first_pos_this_node+npart_this_node] = idpart1s[argsort]
    timers[header+stage] += time.time() - ref; ref = time.time(); header+=1; stage=0 # 7->8

    #  Call again the routine for the 8 subnodes:
    #  Compute positions of subnodes, new level of refinement, 
    #  positions in the array idpart corresponding to the subnodes,
    #  and call for the treatment recursively.
    # nlevel_out=nlevel+1
    idmother_out=inccell
    scale = H.sizeroot * 2**(-nlevel-2)
    timers[-1] += time.time() - ref; ref = time.time()
    for j0 in range(7+1):
        create_KDtree_13110(
            nlevel+1, # nlevel_out,
            # pos_this_node[:3] + H.sizeroot*pos_ref_0[:3,j0]*2**(-nlevel-2), #pos_this_node_out,
            pos_this_node[:3] + scale*pos_ref_0[:3,j0], #pos_this_node_out,
            incsubcell_0[j0], #npart_this_node_out,
            first_pos_this_node+nsubcell_0[j0], #first_pos_this_node_out,
            idmother_out,
            pos_ref_0=pos_ref_0)

#=======================================================================
def remove_degenerate_particles():
#=======================================================================
    raise NotImplementedError('remove_degenerate_particles is not implemented yet')
