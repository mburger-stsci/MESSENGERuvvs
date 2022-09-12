import numpy as np
import pickle
import PIL
import dash
import plotly.io as pio
import plotly.graph_objects as go


class MESSENGERview:
    def __init__(self, data, layer='solar'):
        self.mdata = data
        self.mdata.set_frame('MSO')
        self.layer = layer

        # Create a blank figure
        self.create_figure()
        self.app = self.create_app()
        
    def create_mercury_globe(self, layer):
        # Draw a sphere
        longitude_ = np.linspace(0, 2 * np.pi, 73)
        latitude_ = np.linspace(-np.pi / 2, np.pi / 2, 37)
        longitude, latitude = np.meshgrid(longitude_, latitude_)
        loctime = (longitude * 12 / np.pi + 12) % 24
        ptsx = np.cos(longitude) * np.cos(latitude)
        ptsy = np.sin(longitude) * np.cos(latitude)
        ptsz = np.sin(latitude)
        if layer == 'solar':
            # Solar shading of surface
            ptsc = np.copy(ptsx)
            ptsc[ptsc < 0] = 0
        elif layer == 'abundance':
            # Source map on surface
            sourcemap = pickle.load(open('sourcemap.pkl', 'rb'))
            longitude, latitude = np.meshgrid(sourcemap['longitude'].value,
                                              sourcemap['latitude'].value)
            ptsc = np.zeros_like(longitude)
            ptsc[:-1, :-1] = sourcemap['abundance'].T
            ptsc[:, -1] = ptsc[:, 0]
            ptsc[-1, :] = ptsc[0, :]
        elif layer == 'map':
            # Load in Mercury's surface map
            assert 0
            mercmap = np.array(PIL.Image.open('MercuryMap.tiff'))
            indx = np.interp(latitude.flatten(),
                             np.linspace(-np.pi / 2, np.pi / 2, num=mercmap.shape[0], endpoint=False),
                             np.arange(mercmap.shape[0]))
            indy = np.interp(longitude.flatten(),
                             np.linspace(0, 2 * np.pi, num=mercmap.shape[1], endpoint=False),
                             np.arange(mercmap.shape[1]))
            ptsc = mercmap[indx.astype(int), indy.astype(int), :]
        else:
            assert 0, 'Not set up'

        # Mercury sphere
        customdata = np.array([loctime,
                               longitude*180/np.pi,
                               latitude*180/np.pi,
                               ptsc]).T
        hovertemp = (#'x, y, z: %{x:.1f}, %{y:.1f}, %{z:.1f} <br>'
            'Local Time: %{customdata[0]:.1f} hr <br>'
            'Longitude: %{customdata[1]:.1f}°<br>'
            'Latitude: %{customdata[2]:.1f}°')
        mercury = go.Surface(name='Mercury',
                             x=ptsx, y=ptsy, z=ptsz,
                             surfacecolor=ptsc*256/ptsc.max(),
                             colorscale='Greys',
                             customdata=customdata,
                             hovertemplate=hovertemp,
                             showscale=False, reversescale=True,
                             lighting={'diffuse': 1,
                                       'ambient': 1,
                                       'specular': 0})
        
        self.mercury_figure.add_trace(mercury)
        
        
    def orientation_arrows(self):
        sundir = go.Scatter3d(name='SunDir', x=[0, 4], y=[0, 0], z=[0, 0],
                              mode='lines', hovertext='To Sun',
                              hoverinfo='text',
                              line={'width':5, 'color':'red'}, )
        duskdir = go.Scatter3d(name='DuskDir', x=[0, 0], y=[0, 4], z=[0, 0],
                               mode='lines', hovertext='To Dusk',
                               hoverinfo='text',
                               line={'width':5, 'color':'green'})
        northdir = go.Scatter3d(name='NorthDir', x=[0, 0], y=[0, 0], z=[0, 4],
                                mode='lines', hovertext='To North',
                                hoverinfo='text',
                                line={'width':5, 'color':'blue'})
        
        self.mercury_figure.add_trace(sundir)
        self.mercury_figure.add_trace(duskdir)
        self.mercury_figure.add_trace(northdir)

    def messenger_orbit(self):
        customdata = np.array([self.mdata.data.radiance,
                               self.mdata.data.alttan,
                               self.mdata.data.loctimetan,
                               self.mdata.data.lattan * 180 / np.pi]).T
        hovertemp = ('S/C x, y, z: %{x:.1f}, %{y:.1f}, %{z:.1f} <br>'
                     'Radiance: %{customdata[0]:.1f} kR<br>'
                     'Tangent Altitude: %{customdata[1]:.1f} km<br>'
                     'Tangent Local Time: %{customdata[2]:.1f} hr<br>'
                     'Tangent Latitude: %{customdata[3]:.1f}°')
        
        # Plot the orbit in black
        orbit = go.Scatter3d(name='Orbit', x=self.mdata.data.x, y=self.mdata.data.y,
                                 z=self.mdata.data.z, mode='markers',
                                 marker={'size':3, 'color':'black'},
                                 customdata=customdata, hovertemplate=hovertemp)
        self.mercury_figure.add_trace(orbit)
        self.mercury_figure.update_layout(showlegend=False)
     
    def add_lines_of_sight(self):
        # Add lines to tangent point
        lines_of_sight = []
        hovertemp = ('S/C x, y, z: %{customdata[0]:.1f}, '
                     '%{customdata[1]:.1f}, %{customdata[2]:.1f} <br>'
                     'Radiance: %{customdata[3]:.1f} kR<br>'
                     'Tangent Altitude: %{customdata[4]:.1f} km<br>'
                     'Tangent Local Time: %{customdata[5]:.1f} hr<br>'
                     'Tangent Latitude: %{customdata[6]:.1f}°')
        for i, row in self.mdata.data.iterrows():
            customdata = np.array([[row.x, row.x],
                                   [row.y, row.y],
                                   [row.z, row.z],
                                   [row.radiance, row.radiance],
                                   [row.alttan, row.alttan],
                                   [row.loctimetan, row.loctimetan],
                                   [row.lattan*180/np.pi, row.lattan*180/np.pi]]).T
            line_of_sight = go.Scatter3d(name=f'Line of Sight',
                                         x=[row.x, row.xtan],
                                         y=[row.y, row.ytan],
                                         z=[row.z, row.ztan],
                                         line={'color':[row.radiance, row.radiance],
                                               'cmin':0,
                                               'cmax':self.mdata.data.radiance.max(),
                                               'width':3},
                                         mode='lines',
                                         customdata=customdata,
                                         hovertemplate=hovertemp)
            lines_of_sight.append(line_of_sight)
            self.mercury_figure.add_trace(line_of_sight)
            
        lines_of_sight[0].update(line={'colorbar':{'title':'Radiance (kR)',
                                       'len':0.75, 'lenmode':'fraction'}})
        self.mercury_figure.update_layout(showlegend=False)
    
    def create_figure(self):
        pio.templates.default = 'plotly_white'
    
        self.mercury_figure = go.Figure(layout={'height': 1000,
                                        'width': 1000})
        
        # Add Mercury sphere
        self.create_mercury_globe(self.layer)
        self.orientation_arrows()
        self.messenger_orbit()
        self.add_lines_of_sight()
        
        self.mercury_figure.update_layout(
            {'scene':{'xaxis':{'range':[-6, 6],
                               'title':'Sunward (R_M)'},
                      'yaxis':{'range':[-6, 6],
                               'title':'Duskward (R_M)'},
                      'zaxis':{'range':[-6, 6],
                               'title':'Northward (R_M)'},
                                'aspectratio':{'x':1,
                                               'y':1,
                                               'z':1}
                      }
             })
        
        self.mercury_figure.update_layout(showlegend=False)
    
    def create_app(self):
        title_info = f'''
# MESSENGER UVVS
### {self.mdata.species}, {self.mdata.query}
'''
        app = dash.Dash(f'MESSENGER UVVS Orbit Viewer',
                        external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css'])
        
        app.layout = dash.dash.html.Div([
            dash.dcc.Markdown(title_info),
            dash.dcc.Graph(id='OrbitFigure', figure=self.mercury_figure)])

        return app
