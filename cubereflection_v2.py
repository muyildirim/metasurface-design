import meep as mp
import math
import cmath
import argparse



def main(args):

    # default unit length is 1 um
    um_scale = 1.0



    resolution = 40     # pixels/um

    a = args.aa         # lattice periodicity
    d = args.dd         # cube side length
    h = 1.7             # metal rod height

    tcell= 1.7     	# unitcell thickness
    tsub = 8.85          # substrate thickness
    tpml = 5.0          # PML thickness
    tair = 8.0         # air thickness

    sz = 2*tpml+tair+tcell+tsub
    cell_size = mp.Vector3(a,a,sz)

    pml_layers = [mp.PML(thickness=tpml,direction=mp.Z)]


    lmin = 5.0         # source min wavelength
    lmax = 10.0         # source max wavelength
    fmin = 1/lmax       # source min frequency
    fmax = 1/lmin       # source max frequency
    fcen = 0.5*(fmin+fmax)
    df = fmax-fmin


        # Generate
    nTe = 5.3
    Te = mp.Medium(index=nTe)

    nbaf2=1.4
    BaF2=mp.Medium(index=nbaf2)

    if args.empty:
        geometry = []
    else:
        geometry = [mp.Block(material=Te, size=mp.Vector3(d,d,d),
                     center=mp.Vector3(0,0,0)),
                     mp.Block(material=BaF2, size=mp.Vector3(mp.inf,mp.inf,tsub+tpml),
                              center=mp.Vector3(0,0,-tsub/2.0-d/2.0)) ]

        # CCW rotation angle (degrees) about Y-axis of PW current source; 0 degrees along -z axis
    theta = math.radians(args.theta)

        # k with correct length (plane of incidence: XZ)

    k = mp.Vector3(math.sin(theta),0,math.cos(theta)).scale(fcen)


    def pw_amp(k, x0):
        def _pw_amp(x):

            return cmath.exp(1j * 2 * math.pi * k.dot(x + x0))

        return _pw_amp

    src_pos = d/2.0+tair*0.75
    sources = [ mp.Source(mp.GaussianSource(fcen, fwidth=df), component=mp.Ey, center=mp.Vector3(0,0,src_pos),
                          size=mp.Vector3(a,a,0),
                          amp_func=pw_amp(k, mp.Vector3(0,0,src_pos))) ]

    sim = mp.Simulation(cell_size=cell_size,
                        geometry=geometry,
                        sources=sources,
                        boundary_layers=pml_layers,
                        k_point = k,
                        resolution=resolution)

    nfreq = 50
    refl = sim.add_flux(fcen, df, nfreq, mp.FluxRegion(center=mp.Vector3(0,0,+d/2.0+0.5*tair),size=mp.Vector3(a,a,0)))

    if not args.empty:
        sim.load_minus_flux('refl-flux', refl)

    sim.run(mp.at_beginning(mp.output_epsilon),
                  until_after_sources=mp.stop_when_fields_decayed(25, mp.Ey, mp.Vector3(0,0,d/2.0+0.5*tair), 1e-3))

    if args.empty:
        sim.save_flux('refl-flux', refl)

    sim.display_fluxes(refl)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-empty', action='store_true', default=False, help="empty? (default: False)")
    parser.add_argument('-aa', type=float, default=3.4, help='lattice periodicity (default: 4.5 um)')
    parser.add_argument('-dd', type=float, default=1.7, help='cube side length (default: 1.7 um)')
    parser.add_argument('-theta', type=float, default=0, help='angle of planewave current source (default: 0 degrees)')
    args = parser.parse_args()
    main(args)
