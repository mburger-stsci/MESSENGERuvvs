"""MESSENGER UVVS data class"""
import numpy as np
from scipy.optimize import minimize_scalar
import pandas as pd
import bokeh.plotting as bkp
from bokeh.models import HoverTool, Whisker
from bokeh.palettes import Set1
from bokeh.io import export_png
from astropy import units as u
from astropy.visualization import PercentileInterval
from .database_setup import database_connect
from nexoclom import Input


class InputError(Exception):
    """Raised when a required parameter is not included in the inputfile."""
    def __init__(self, expression, message):
        self.expression = expression
        self.message = message


class MESSENGERdata:
    """Retrieve MESSENGER data from database.
    Given a species and set of comparisons, retrieve MESSSENGER UVVS
    data from the database. The list of searchable fields is given at
    :doc:`database_fields`.
    
    Returns a MESSENGERdata Object.
    
    **Parameters**
    
    species
        Species to search. This is required because the data from each
        species is stored in a different database table.
        
    query
        A SQL-style list of comparisons.
        
    The data in the object created is extracted from the database tables using
    the query:
    
    ::
    
        SELECT *
        FROM <species>uvvsdata, <species>pointing
        WHERE <query>
    
    See examples below.
    
    **Class Atributes**
    
    species
        The object can only contain a single species.
        
    frame
        Coordinate frame for the data, either MSO or Model.
        
    query
        SQL query used to search the database and create the object.
    
    data
        Pandas dataframe containing result of SQL query. Columns in the
        dataframe are the same as in the database except *frame* and
        *species* have been dropped as they are redundant. If models have been
        run, there are also columns in the form modelN for the Nth model run.
        
    taa
        Median true anomaly for the data in radians.
        
    model_label
        If *N* models have been run, this is a dictionary in the form
        `{'model0':label0, ..., 'modelN':labelN}` containing descriptions for
        the models.
        
    model_strength
        If *N* models have been run, this is a dictionary in the form
        `{'model0':strength0, ..., 'modelN':strengthN}` containing modeled
        source rates in units of :math:`10^{26}` atoms/s.
        
    **Examples**
    
    1. Loading data
    
    ::
    
        >>> from MESSENGERuvvs import MESSENGERdata
        
        >>> CaData = MESSENGERdata('Ca', 'orbit = 36')
        
        >>> print(CaData)
        Species: Ca
        Query: orbit = 36
        Frame: MSO
        Object contains 581 spectra.
        
        >>> NaData = MESSENGERdata('Na', 'orbit > 100 and orbit < 110')
        
        >>> print(NaData)
        Species: Na
        Query: orbit > 100 and orbit < 110
        Frame: MSO
        Object contains 3051 spectra.
        
        >>> MgData = MESSENGERdata('Mg',
                'loctimetan > 5.5 and loctimetan < 6.5 and alttan < 1000')
        
        >>> print(len(MgData))
        45766
        
    2. Accessing data.
    
    * The observations are stored within the MESSENGERdata object in a
      `pandas <https://pandas.pydata.org>`_ dataframe attribute called *data*.
      Please see the `pandas documentation <https://pandas.pydata.org>`_ for
      more information on how to work with dataframes.
    
    ::
    
        >>> print(CaData.data.head(5))
                                 utc  orbit  merc_year  ...  loctimetan         slit               utcstr
        unum                                            ...
        3329 2011-04-04 21:24:11.820     36          0  ...   14.661961  Atmospheric  2011-04-04T21:24:11
        3330 2011-04-04 21:25:08.820     36          0  ...   12.952645  Atmospheric  2011-04-04T21:25:08
        3331 2011-04-04 21:26:05.820     36          0  ...   12.015670  Atmospheric  2011-04-04T21:26:05
        3332 2011-04-04 21:27:02.820     36          0  ...   12.007919  Atmospheric  2011-04-04T21:27:02
        3333 2011-04-04 21:27:59.820     36          0  ...   12.008750  Atmospheric  2011-04-04T21:27:59
        
        [5 rows x 29 columns]

    * Individual observations can be extracted using standard Python
      slicing techniques:
     
    ::
        
        >>> print(CaData[3:8])
        Species: Ca
        Query: orbit = 36
        Frame: MSO
        Object contains 5 spectra.

        >>> print(CaData[3:8].data['taa'])
        unum
        3332    1.808107
        3333    1.808152
        3334    1.808198
        3335    1.808243
        3336    1.808290
        Name: taa, dtype: float64

    3. Modeling data
    
    ::
    
        >>> inputs = Input('Ca.spot.Maxwellian.input')
        >>> CaData.model(inputs, 1e5, label='Model 1')
        >>> inputs..speeddist.temperature /= 2.  # Run model with different temperature
        >>> CaData.model(inputs, 1e5, label='Model 2')
        
    4. Plotting data
    
    ::
    
        >>> CaData.plot('Ca.orbit36.models.html')
    
    5. Exporting data to a file
    
    ::
    
        >>> CaData.export('modelresults.csv')
        >>> CaData.export('modelresults.html', columns=['taa'])

    
    """
    def __init__(self, species=None, comparisons=None):
        allspecies = ['Na', 'Ca', 'Mg']
        self.species = None
        self.frame = None
        self.query = None
        self.data = None
        self.taa = None
        self.inputs = None
        self.model_info = None
        
        if species is None:
            pass
        elif species not in allspecies:
            # Return list of valid species
            print(f"Valid species are {', '.join(allspecies)}")
        elif comparisons is None:
            # Return list of queryable fields
            with database_connect() as con:
                columns = pd.read_sql(
                    f'''SELECT * from {species}uvvsdata, {species}pointing
                        WHERE 1=2''', con)
            print('Available fields are:')
            for col in columns.columns:
                print(f'\t{col}')
        else:
            # Run the query and try to make the object
            query = f'''SELECT * from {species}uvvsdata, {species}pointing
                        WHERE unum=pnum and {comparisons}
                        ORDER BY unum'''
            try:
                with database_connect() as con:
                    data = pd.read_sql(query, con)
            except Exception:
                raise InputError('MESSENGERdata.__init__',
                                 'Problem with comparisons given.')

            if len(data) > 0:
                self.species = species
                self.frame = data.frame[0]
                self.query = comparisons
                data.drop(['species', 'frame'], inplace=True, axis=1)
                self.data = data
                self.data.set_index('unum', inplace=True)
                self.taa = np.median(data.taa)
            else:
                print(query)
                print('No data found')
                
    def __str__(self):
        result = (f'Species: {self.species}\n'
                  f'Query: {self.query}\n'
                  f'Frame: {self.frame}\n'
                  f'Object contains {len(self)} spectra.')
        return result

    def __repr__(self):
        result = ('MESSENGER UVVS Data Object\n'
                  f'Species: {self.species}\n'
                  f'Query: {self.query}\n'
                  f'Frame: {self.frame}\n'
                  f'Object contains {len(self)} spectra.')
        return result

    def __len__(self):
        try:
            return len(self.data)
        except Exception:
            return 0

    def __getitem__(self, q_):
        if isinstance(q_, int):
            q = slice(q_, q_+1)
        elif isinstance(q_, slice):
            q = q_
        elif isinstance(q_, pd.Series):
            q = np.where(q_)[0]
        else:
            raise TypeError

        new = MESSENGERdata()
        new.species = self.species
        new.frame = self.frame
        new.query = self.query
        new.taa = self.taa
        new.data = self.data.iloc[q].copy()
        new.model_info = self.model_info
        new.inputs = self.inputs

        return new

    def __iter__(self):
        for i in range(len(self.data)):
            yield self[i]

    def keys(self):
        """Return all keys in the object, including dataframe columns"""
        keys = list(self.__dict__.keys())
        keys.extend([f'data.{col}' for col in self.data.columns])
        return keys

    def set_frame(self, frame=None):
        """Convert between MSO and Model frames.

        More frames could be added if necessary.
        If Frame is not specified, flips between MSO and Model."""
        if (frame is None) and (self.frame == 'MSO'):
            frame = 'Model'
        elif (frame is None) and (self.frame == 'Model'):
            frame = 'MSO'
        else:
            pass

        allframes = ['Model', 'MSO']
        if frame not in allframes:
            print('{} is not a valid frame.'.format(frame))
            return None
        elif frame == self.frame:
            pass
        elif (self.frame == 'MSO') and (frame == 'Model'):
            # Convert from MSO to Model
            self.data.x, self.data.y = self.data.y.copy(), -self.data.x.copy()
            self.data.xbore, self.data.ybore = (self.data.ybore.copy(),
                                                -self.data.xbore.copy())
            self.data.xtan, self.data.ytan = (self.data.ytan.copy(),
                                              -self.data.xtan.copy())
            self.frame = 'Model'
        elif (self.frame == 'Model') and (frame == 'MSO'):
            self.data.x, self.data.y = -self.data.y.copy(), self.data.x.copy()
            self.data.xbore, self.data.ybore = (-self.data.ybore.copy(),
                                                self.data.xbore.copy())
            self.data.xtan, self.data.ytan = (-self.data.ytan.copy(),
                                              self.data.xtan.copy())
            self.frame = 'MSO'
        else:
            assert 0, 'You somehow picked a bad combination.'

    def model(self, inputs_, npackets, quantity='radiance',
              fit_method='middle50', dphi=3*u.deg, overwrite=False,
              filenames=None, label=None):

        if isinstance(inputs_, str):
            inputs = Input(inputs_)
        elif hasattr(inputs_, 'line_of_sight'):
            inputs = inputs_
        else:
            raise InputError('MESSENGERdata.model', 'Problem with the inputs.')

        # TAA needs to match the data
        oldtaa = inputs.geometry.taa
        if len(self.data.orbit.unique()) == 1:
            inputs.geometry.taa = np.median(self.data.taa)*u.rad
        elif self.data.taa.max()-self.data.taa.min() < 3*np.pi/180.:
            inputs.geometry.taa = self.data.taa.median()*u.rad
        else:
            assert 0, 'Too wide a range of taa'
            
        # If using a planet-fixed source map, need to set subsolarlon
        if ((inputs.spatialdist.type == 'surface map') and
            (inputs.spatialdist.coordinate_system == 'planet-fixed')):
            inputs.spatialdist.subsolarlon = self.data.subslong.median()*u.rad
        else:
            pass

        # Run the model
        inputs.run(npackets, overwrite=overwrite)

        # Simulate the data
        if self.inputs is None:
            self.inputs = [inputs]
        else:
            self.inputs.append(inputs)

        self.set_frame('Model')
        model_result = inputs.line_of_sight(self.data, quantity,
                                            dphi=dphi, filenames=filenames,
                                            overwrite=overwrite)
        
        # modkey is the number for this model
        modkey = f'model{len(self.inputs)-1:00d}'
        packkey = f'packets{len(self.inputs)-1:00d}'
        self.data[modkey] = model_result.radiance/1e3 # Convert to kR
        self.data[packkey] = model_result.packets
        
        strength, goodness_of_fit = self.fit_model(modkey, fit_method)
        self.data[modkey] = self.data[modkey]*strength.value

        if label is None:
            label = modkey.capitalize()
        else:
            pass

        model_info = {'fit_method': fit_method,
                      'goodness-of-fit': goodness_of_fit,
                      'strength': strength,
                      'label': label}
        if self.model_info is None:
            self.model_info = {modkey: model_info}
        else:
            self.model_info[modkey] = model_info
        
        print(f'Model strength for {label} = {strength}')

        # Put the old TAA back in.
        inputs.geometry.taa = oldtaa
    
    def fit_model(self, modkey, fit_method):
        def chisq(x):
            return np.sum((self.data[mask].radiance -
                           x * self.data[mask][modkey])**2 /
                          self.data[mask].sigma[mask]**2)
        
        if fit_method.lower().startswith('middle'):
            perinterval = float(fit_method[6:])
            # Estimate model strength (source rate) by fitting middle %
            interval = PercentileInterval(perinterval)
            lim = interval.get_limits(self.data.radiance)
            mask = ((self.data.radiance >= lim[0]) &
                    (self.data.radiance <= lim[1]))
    
            strunit = u.def_unit('10**26 atoms/s', 1e26/u.s)
            model_strength = minimize_scalar(chisq)
            return model_strength.x * strunit, model_strength.fun
        else:
            raise InputError('MESSENGERdata.fit_model',
                             f'fitmethod = {fit_method} not defined.')

    def plot(self, filename=None, show=True, **kwargs):
        if filename is not None:
            if not filename.endswith('.html'):
                filename += '.html'
            else:
                pass
            bkp.output_file(filename)
        else:
            pass

        # Format the date correction
        self.data['utcstr'] = self.data['utc'].apply(
            lambda x: x.isoformat()[0:19])

        # Put the dataframe in a useable form
        self.data['lower'] = self.data.radiance - self.data.sigma
        self.data['upper'] = self.data.radiance + self.data.sigma
        source = bkp.ColumnDataSource(self.data)

        # Make the figure
        fig = bkp.figure(plot_width=1200, plot_height=800,
                         x_axis_type='datetime',
                         title=f'{self.species}, {self.query}',
                         x_axis_label='UTC',
                         y_axis_label='Radiance (kR)',
                         tools=['pan', 'box_zoom', 'reset', 'save'])

        # plot the data
        dplot = fig.circle(x='utc', y='radiance', size=7, color='black',
                           legend_label='Data', hover_color='white',
                           source=source)
        
        # Add error bars
        fig.add_layout(Whisker(source=source, base='utc', upper='upper',
                               lower='lower'))

        # tool tips
        tips = [('index', '$index'),
                ('UTC', '@utcstr'),
                ('Radiance', '@radiance{0.2f} kR')]
        if self.model_info is not None:
            for modkey, info in self.model_info.items():
                tips.append((info['label'], f'@{info["strength"]}{{0.2f}} kR'))
        datahover = HoverTool(tooltips=tips,
                              renderers=[dplot])
        fig.add_tools(datahover)

        # Plot the model
        col = (c for c in Set1[9])
        if self.model_info is not None:
            for modkey, info in self.model_info.items():
                try:
                    c = next(col)
                except StopIteration:
                    col = (c for c in Set1[9])
                    c = next(col)

                f0 = fig.line(x='utc', y=modkey, source=source,
                               legend_label=info['label'], color=c)
                f1 = fig.circle(x='utc', y=modkey, size=7, source=source,
                               legend_label=info['label'], color=c)
                datahover.renderers.append(f0)
                datahover.renderers.append(f1)

        # Labels, etc.
        fig.title.align = 'center'
        fig.title.text_font_size = '16pt'
        fig.axis.axis_label_text_font_size = '16pt'
        fig.axis.major_label_text_font_size = '16pt'
        fig.legend.label_text_font_size = '16pt'
        fig.legend.click_policy = 'hide'

        if filename is not None:
            bkp.output_file(filename)
            export_png(fig, filename=filename.replace('.html', '.png'))
            bkp.save(fig)
            
        else:
            pass

        if show:
            bkp.show(fig)

    def export(self, filename, columns=('utc', 'radiance')):
        """Export data and models to a file.
        **Parameters**
        
        filename
            Filename to export model results to. The file extension determines
            the format. Formats available: csv, pkl, html, tex
            
        columns
            Columns from the data dataframe to export. Available columns can
            be found by calling the `keys()` method on the data object.
            Default = ['utc', 'radiance'] and all model result columns. Note:
            The default columns are always included in the output
            regardless of whether they are specified.
        
        **Returns**
        
        No outputs.
        
        """
        columns_ = list(columns)
        if self.model_info is not None:
            columns_.extend(self.model_info.keys())
        else:
            pass
        
        # Make sure radiance is in there
        if 'radiance' not in columns_:
            columns_.append('radiance')
        else:
            pass

        # Make sure UTC is in there
        if 'utc' not in columns_:
            columns_.append('utc')
        else:
            pass

        if len(columns_) != len(set(columns_)):
            columns_ = list(set(columns_))
        else:
            pass

        for col in columns_:
            if col not in self.data.columns:
                columns_.remove(col)
            else:
                pass

        subset = self.data[columns_]
        if filename.endswith('.csv'):
            subset.to_csv(filename)
        elif filename.endswith('.pkl'):
            subset.to_pickle(filename)
        elif filename.endswith('.html'):
            subset.to_html(filename)
        elif filename.endswith('.tex'):
            subset.to_latex(filename)
        else:
            print('Valid output formats = csv, pkl, html, tex')


