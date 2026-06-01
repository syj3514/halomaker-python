import sys
import argparse
# from halo_defs import *
sys.stdout.reconfigure(line_buffering=True)
sys.stderr.reconfigure(line_buffering=True)

def main():
    # this program reads in the particle data from N-Body simulation snapshots and 
    # puts them into halos (using a friend-of-friend algorithm). It then goes on to
    # compute various properties for these halos (mass, spin, energy ...)


    print()
    print( '_______________________________________________________________________')
    print()
    print( '              HalaxyMaker                                                ')
    print( '              ---------                                                ' )
    print()
    print()
    print( ' first version    :        S. Ninin                  (1999)            ')
    print( ' modif.           :        J. Devriendt              (1999-2002)       ')
    print( ' modif.           :        B. Lanzoni & S. Hatton    (2000-2001)       ')
    print( ' galics v1.0      :        J. Blaizot & J. Devriendt (2002)            ')
    print( ' horizon v1.0     :        J. Devriendt & D. Tweed   (2006)            ')
    print( ' horizon v2.0     :        J. Devriendt & D. Tweed   (2007)            ')
    print( ' horizon v2.0.2   :        J. Devriendt, D. Tweed & J.Blaizot (2008)   ')
    print( ' horizon python   :        S. Jeon                   (2026)            ')
    print( )
    print( '_______________________________________________________________________'  )
    print( )



    import halo_defs as H
    import compute_halo_props
    # get directory where to input/output the data
    print("Usage: python HaloMaker.py [output_dir]")
    parser = argparse.ArgumentParser(description='Run HaloMaker')
    parser.add_argument('output_dir', nargs='?', default='.', help='Directory for output data (default: current directory)')
    # parser.add_argument('--prefix', default='', help='Prefix for output files (default: "")')
    for attr, (aliases, dtype) in H.PARAMS.items():
        flags = [f'--{attr}']

        for alias in aliases:
            flag = f'--{alias}'
            if flag not in flags:
                flags.append(flag)
        
        parser.add_argument(
            *flags,
            dest=attr,
            type=dtype,
            default=None,
        )
    H.args = parser.parse_args()
    H.output_dir = H.args.output_dir
    # if(len(sys.argv)<2): H.output_dir = '.'
    # else: H.output_dir = sys.argv[1]
    # initialize cosmological and technical parameters of the N_Body simulation 
    print( )
    print( '_______________________________________________________________________'  )
    print( )
    print("          Initialization")
    print("          --------------")
    print( )
    print( '_______________________________________________________________________'  )
    compute_halo_props.init_0()

    # loop over snapshots
    print( )
    print( '_______________________________________________________________________'  )
    print( )
    print("          Loop over snapshots")
    print("          -------------------")
    print( )
    print( '_______________________________________________________________________'  )
    for istep1 in range(1,H.nsteps+1,1):
        H.numero_step = istep1
        compute_halo_props.new_step_1()
        print(H.numero_step,' of ',H.nsteps)
        print()
        print()

    print()
    print( '> Bricks now available to build halo merger tree with TreeMaker')
    print()
    print( '_______________________________________________________________________'  )
    print()
    print( '      End of HalaxyMaker                                        ')
    print()
    print( '_______________________________________________________________________'  )

if __name__ == "__main__":
    main()