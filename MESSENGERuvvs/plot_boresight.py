

oedge = min([output.inputs.options.outer_edge*2, 100*mercury.unit])
far_pt = sc_xyz + oedge*sc_bore
radius = oedge*np.sin(self.dphi)

vec_in_circle = np.column_stack([np.ones_like(data.xbore),
                                 np.ones_like(data.ybore),
                                 -(data.xbore+data.ybore)/data.zbore])
vec_in_circle /= np.sqrt(np.sum(vec_in_circle**2, axis=1))[:,np.newaxis]
bore_cross_vec = np.cross(sc_bore, vec_in_circle)
cone = np.append(circle, [sc_xyz[0,:]], axis=0)

circle = (radius*np.cos(phi)[:,np.newaxis]*vec_in_circle[i,:] +
          radius*np.sin(phi)[:,np.newaxis]*bore_cross_vec[i,:] +
          far_pt[i,:])
