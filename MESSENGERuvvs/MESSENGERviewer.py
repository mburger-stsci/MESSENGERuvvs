import dash
from MESSENGERuvvs.MESSENGERview import MESSENGERview
from MESSENGERuvvs.MESSENGERdata import MESSENGERdata


class MESSENGERviewer(MESSENGERview):
    def __init__(self, species='Na', query=None, layer='solar', show=True):
        self.species = species
        self.query = query
        self.load_data()
        super().__init__(self.mdata)
        self.layer = layer

        # Create a blank figure
        self.create_figure()

        self.app = self.create_app()
        if show:
            self.show()
        
    def load_data(self):
        try:
            data = MESSENGERdata(self.species, self.query)
            usedefault = data.data is None
        except:
            usedefault = True
            
        if usedefault:
            data = MESSENGERdata(self.species, 'orbit=22')
            print(f'No data for species={self.species} and query="{self.query}".')
            print('Using default query.')
            self.query = 'orbit=22'
            data.set_frame('MSO')
        else:
            pass
        
        self.mdata = data
        
    def remove_trace(self, name):
        data = list(self.mercury_figure.data)
        data_ = []
    
        for i, part in enumerate(data):
            if part['name'] != name:
                data_.append(part)
            else:
                pass
        self.mercury_figure.data = data_

    def create_app(self):
        title_info = f'''
# MESSENGER UVVS
### {self.mdata.species}, {self.mdata.query}
'''
        app = dash.Dash(f'MESSENGER UVVS Orbit Viewer',
                        external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css'])
        app.layout = dash.dash.html.Div([
            dash.dcc.Markdown(title_info, id='TitleInfo'),
            dash.html.Table([
                dash.html.Tr([  # Top row: Species, show MESSENGER orbit, Show Model
                    dash.html.Td([
                        dash.html.Label('Species'),
                        dash.dcc.Dropdown(id='SpeciesDropdown',
                                          options=[{'label': 'Calcium', 'value': 'Ca'},
                                                   {'label':'Sodium', 'value':'Na'},
                                                   {'label':'Magnesium', 'value':'Mg'}],
                                          searchable=False, clearable=False,
                                          value=self.species, style={'width': 140}),
                    ]),
                    dash.html.Td([
                        dash.html.Label('MESSENGER Query'),
                        dash.dcc.Input(id='MesQuery',
                                       value=self.query,
                                       style={'width': 300},
                                       debounce=True,
                                       type='text'),
                        ], id='MesSquare'),
                ]),
            ]),
            dash.dcc.Graph(id='OrbitFigure', figure=self.mercury_figure)
        ])

        # Reveal Model panel
        @app.callback(
             dash.Output('TitleInfo', 'text'),
             dash.Output('MesQuery', 'value'),
             dash.Output('OrbitFigure', 'figure'),
             
             dash.Input('SpeciesDropdown', 'value'),
             dash.Input('MesQuery', 'value'))
        def update_figure(species, query):
            if (species != self.species) or (query != self.query):
                self.species = species
                self.query = query
                self.load_data()
                self.remove_trace('Orbit')
                self.messenger_orbit()

                self.remove_trace('Line of Sight')
                self.add_lines_of_sight()
            else:
                pass
            
            title_info = f'''
# MESSENGER UVVS
### {self.mdata.species}, {self.mdata.query}
'''
            self.mercury_figure.update_layout()
            print('figure updated')
            return title_info, self.query, self.mercury_figure

        return app
    
    def show(self):
        self.app.run_server()
