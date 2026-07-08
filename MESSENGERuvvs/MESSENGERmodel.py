import numpy as np
from sklearn.neighbors import KDTree
import astropy.units as u
from astropy.modeling import models, fitting
from astropy.visualization import PercentileInterval
from astropy.time import Time
import h5py
from nexoclom2.data_simulation.ModelResult import ModelResult


class MMPacket:
    def __init__(self):
        # self.x = None
        # self.y = None
        # self.z = None
        self.time = None
        self.X = None
        self.vx = None
        self.vy = None
        self.vz = None
        self.frac = None
        
    def __getitem__(self, q):
        if self.time is not None:
            new = MMPacket()
            new.time = self.time[q]
            new.X = self.X[q]
            new.vx = self.vx[q]
            new.vy = self.vy[q]
            new.vz = self.vz[q]
            new.frac = self.frac[q]
            return new
        else:
            return None


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
        
        assert output.species == scdata.species, 'output.species =/= scdata.species'
        
        self.inputs = output.inputs
        self.species = scdata.species
        self.query = scdata.comparisons
        self.type = 'LineOfSight'
        self.dphi = params.get('dphi', 1*u.deg)
        self.radiance = np.zeros(len(scdata))*u.kR
        self.column = np.zeros(len(scdata))/u.cm**2
        self.packets = np.zeros(len(scdata), dtype=int)
        self.compute_radiance(scdata, output)
        fitfactor, goodness = self.determine_source_rate(scdata)
        self.radiance *= fitfactor
        self.sourcerate = fitfactor * output.sourcerate
        self.goodness_of_fit = goodness
        
    def compute_radiance(self, data, output, chunksize=1_000_000):
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
        
        # Distance along line of sight to the edge of the model
        # (boresight * t + x_sc)**2 = (outer_edge)**2
        a = np.sum((sc_bore*mercury.unit)**2, axis=1)  ## == 1
        b = 2*np.sum(sc_xyz*sc_bore*mercury.unit, axis=1)
        c = np.sum(sc_xyz**2, axis=1) - output.inputs.options.outer_edge**2
        t_edge = (-b + np.sqrt(b**2-4*a*c))/(2*a)
        
        print('Determining radiance')
        with h5py.File(output.savefile) as store:
            final_state = store['final_state']
            
            t0 = Time.now()
            packets = MMPacket()
            packets.time = final_state['time'][:]*u.s
            packets.X = np.column_stack([final_state['x'][:],
                                         final_state['y'][:],
                                         final_state['z'][:]])*output.unit
            packets.vx = final_state['vx'][:]*output.unit/u.s
            packets.vy = final_state['vy'][:]*output.unit/u.s
            packets.vz = final_state['vz'][:]*output.unit/u.s
            packets.frac = final_state['frac'][:]
            # frac = final_state['frac'][:]
            # vx = final_state['vx'][:]*output.unit/u.s
            # vy = final_state['vy'][:]*output.unit/u.s
            # vz = final_state['vz'][:]*output.unit/u.s
            tree = KDTree(packets.X)
            t1 = Time.now()
            print((t1-t0).to(u.s))
            
            for i in range(len(data)):
                # Find points along line of sight
                xyz, bore = sc_xyz[i,:], sc_bore[i,:]
                
                t = np.linspace(0, t_edge[i], 100)*mercury.unit
                pts_los = (xyz[np.newaxis,:] + np.outer(t, bore))
                width = 2*t*np.sin(self.dphi)
                inds = np.unique(np.concatenate(tree.query_radius(pts_los, width)))
                packets_inds = packets[inds]
                # packets = allpackets[inds]
                # t = final_state['time'][inds]*u.s
                # x = final_state['x'][inds]*output.unit
                # y = final_state['y'][inds]*output.unit
                # z = final_state['z'][inds]*output.unit
                # frac = final_state['frac'][inds]
                # vx = final_state['vx'][inds]*output.unit/u.s
                # vy = final_state['vy'][inds]*output.unit/u.s
                # vz = final_state['vz'][inds]*output.unit/u.s
                
                # angle between points and s/c boresight
                x_rel_sc = packets.X[inds,:] - sc_xyz[i,:]
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
                    sub = packets_inds[in_view]
                    # sub = MMPacket()
                    # sub.time = t[in_view]
                    # sub.x = los_xyz[:,0]
                    # sub.y = los_xyz[:,1]
                    # sub.z = los_xyz[:,2]
                    # sub.X = np.column_stack([sub.x, sub.y, sub.z])
                    # sub.vx = vx[in_view]
                    # sub.vy = vx[in_view]
                    # sub.vz = vx[in_view]
                    # sub.frac = frac[in_view]
                    
                    # compute column density
                    rad_per_atom = self.radiance_per_atom(sub, output)
                    area_pix = np.pi * (r_rel_sc[in_view]*np.sin(self.dphi))**2
                    col = sub.frac * output.atoms_per_packet/area_pix.to(u.cm**2)
                    self.column[i] += col.sum()
                    self.radiance[i] += np.sum(col*rad_per_atom)
                    self.packets[i] += in_view.sum()
                else:
                    pass
                
            
    def make_mask(self, scdata):
        mask = np.array([True for _ in scdata.radiance])
        return mask

    def determine_source_rate(self, scdata, use_weights=True, mask=None):
        if mask is None:
            mask = self.make_mask(scdata)
        else:
            pass
        if use_weights:
            weights = 1./scdata.data.sigma.values[mask]**2
        else:
            weights = np.ones_like(scdata.sigma[mask])
            
        linmodel = models.Multiply()
        fitter = fitting.LinearLSQFitter()
        if not np.all(self.radiance == 0*u.kR):
            best_fit = fitter(linmodel, self.radiance[mask],
                              scdata.radiance[mask],
                              weights=weights)
            
            return best_fit.factor.value, None
        else:
            return 0, None
