import numpy as np
import astropy.units as u
from scipy.spatial import ConvexHull
from nexoclom2.data_simulation.ModelResult import ModelResult


class MESSENGERModel(ModelResult):
    """ Simulated MESSENGER UVVS result from a model Output. Returns the
    simulated radiance, column density, and number of packets sampled.
    
    Parameters
    ----------
    scdata: MESSENGERdata object
        Subset of the MESSENGER UVVS data
        
    output: nexoclom2 Output object
        Result from a model run.
        
    params: dict
        Possible fields
        
        * wavelength: astropy distance Quantity. There is some logic to
        determine default radiance for species
        
        * dphi: astropy angle Quantity. Angular size of aperture opening.
        Default = 1º
        
    Attributes
    ----------
    species: nexoclom2 Atom object
    query: str
    radiance
    column
    
    Notes
    -----
    * The field of view includes all packets within the angle dphi of the
    UVVS boresight
    
    """
    def __init__(self, scdata, output, params=None):
        if params is None:
            params = {'quantity': 'radiance'}
        else:
            params['quantity'] = 'radiance'
        super().__init__(output, params)
        
        assert output.frame == 'MERCURYSOLAR', 'output.frame must be MERCURYSOLAR'
        assert output.species == scdata.species, 'output.species =/= scdata.species'
        
        self.species = scdata.species
        self.query = scdata.comparisons
        self.type = 'LineOfSight'
        self.dphi = params.get('dphi', 1*u.deg)
        self.radiance = np.zeros(len(scdata))*u.kR
        self.column = np.zeros(len(scdata))/u.cm**2
        self.packets = np.zeros(len(scdata), dtype=int)
        self.compute_radiance(scdata, output)
        
    def compute_radiance(self, data, output, chunksize=5000000):
        st = 0
        while st < output.n_final_packets:
            fin = st + np.min([st+chunksize, output.n_final_packets])
            st += fin
            packets = output.final_state(which=range(st, fin))
            
            mercury = output.objects['Mercury']
            dist_from_plan = np.sqrt(data.x**2 + data.y**2 + data.z**2)
            
            # Angle between look direction and planet
            ang = np.arccos((-data.x*data.xbore - data.y*data.ybore -
                             data.z*data.zbore)/dist_from_plan)
            
            asize_plan = np.arcsin(1*mercury.unit/dist_from_plan)
            
            # If LOS doesn't hit planet, integrate to infinity
            dist_from_plan[ang > asize_plan] = 1e30*mercury.unit
            
            sc_xyz = np.column_stack([data.x, data.y, data.z])
            sc_bore = np.column_stack([data.xbore, data.ybore, data.zbore])

            for i in range(len(data)):
                # angle between points and s/c boresight
                x_rel_sc = packets.X() - sc_xyz[i,:]
                r_rel_sc = np.sqrt(np.sum(x_rel_sc**2, axis=1))
                cos_theta = np.sum(x_rel_sc*sc_bore[i,:], axis=1)/r_rel_sc
                
                # Determine the projection of the packet onto the LOS in view of s/c
                in_view = ((cos_theta > np.cos(self.dphi)) &
                           (r_rel_sc*cos_theta < dist_from_plan[i]))
                
                if np.any(in_view):
                    # distance from s/c along los
                    los_r = r_rel_sc[in_view]*cos_theta[in_view]
                    
                    # locations of packets projected onto los rel Mercury
                    los_xyz = sc_xyz[i,:] + los_r[:,np.newaxis]*sc_bore[i,:]
                    sub = packets[in_view]
                    sub.x = los_xyz[:,0]
                    sub.y = los_xyz[:,1]
                    sub.z = los_xyz[:,2]
                    sub.X = sub.X()
                    sub.V = sub.V()
                    rad_per_atom = self.radiance_per_atom(sub, output)
                    
                    # compute column density
                    area_pix = np.pi * (r_rel_sc[in_view]*np.sin(self.dphi))**2
                    col = sub.frac * output.atoms_per_packet/area_pix
                    self.column[i] += col.sum()
                    self.radiance[i] += np.sum(col*rad_per_atom)
                    self.packets[i] += in_view.sum()
                else:
                    pass
                
            del packets
