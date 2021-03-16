# -*- coding: utf-8 -*-
#
# (c) Copyright 2003, 2004, 2005
#     Author: Ola Skavhaug
#     Simula Research Laboratory AS
#     
#     This file is part of PyPDE.
#     As PyPDE is unfinshed and not distributed, this is distributed as part
#     of PyFDM / PySE.
#
#     PyPDE is free software; you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation; either version 2 of the License, or
#     (at your option) any later version.
#
#     PyPDE is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with PyPDE; if not, write to the Free Software
#     Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA



"""A general and simple plotting facility for pyPDE and friends.
VizCurve is equipped with a Matlab-like interface, and plots using Gnuplot, BTL
or pyx.

Some hints:
plot(x, y, hardcopy=True) turns on hardcopy with a default filename
plot(x, y, hardcopy='myfilename.eps) turns on hardcopy to file 'myfilename.eps'

"""

import Gnuplot
from plot_utils import *
import plot_utils
try:
    from pyx import graph as pyxgraph
    from pyx import data as pyxdata
    from pyx import canvas as pyxcanvas
    from pyx import text as pyxtext
except:
    pass
from pyPDE import seq, concatenate as _concatenate, NewAxis as _NewAxis
from os import path as _path

class FieldData(object):
    """Curve data. Needed to save the state of plots."""
    def __init__(self, plotter):
        data_dict = plotter._field_attrs
        self.legend = data_dict["legend"]
        self.color = data_dict["color"]
        self.line_type = data_dict["line_type"]
        self.x = data_dict["x"]
        self.y = data_dict["y"]

def manual():
    print "Write me!"

class VizCurveBase(object): 
    """
Base class for curve-plots. Dictionaries are used to store plotting
parameters:
    self._figure_attrs: Label info, scaling, movie, hardcopy, hold
    self._axes_attrs: min/max of axes.
    self._field_attrs: curve information; x, y, legend, and line-type"""
    def __init__(self, *args, **kwargs):
        self._set_defaults()

        self.configure_field_attrs(*args, **kwargs)
        self.configure_axes_attrs(**kwargs)
        self.configure_figure_attrs(**kwargs)
        self.fieldDataList = []

    def _set_defaults(self):
        self._figure_attrs = {
            'xlabel':   "x",
            'ylabel':   "y",
            'title':    "y(x)",
            'scale':    "linear",
            'filename': None,
            'movie':     "off",
            'frame_count': 0,
            'hardcopy': False,
            'default_filename': "tmp.eps",
            'screen_plot': True,
            'hold': 'off'
        }

        self._axes_attrs = {
# Axes used for plotting:
            'axmin':        None,
            'axmax':        None,
            'aymin':        None,
            'aymax':        None,
# User specified axes ranges (if given):
            'xmin':         None,
            'xmax':         None,
            'ymin':         None,
            'ymax':         None,
            'fixed_view':   False
        }

        self._field_attrs = {
            'legend':       "",
            'color':        "blue",
            'line_type':    "solid",
            'x':            None,
            'y':            None
        }

    def set_field_data(self, *args):
        if len(args) == 2:
            self._field_attrs['x'] = args[0]
            self._field_attrs['y'] = args[1]
        elif len(args) == 1:
            if not self._field_attrs['x']:
                self._field_attrs['x'] = seq(len(args[0])-1) # From 0 to len(args[1])-1
            else:
                if self._figure_attrs['hold'] == 'on':
                    print "Warning: x values not given, this may be dangerous"
            self._field_attrs['y'] = args[0]

    def configure_field_attrs(self, *args, **kwargs):
        """Configure parameters for the field (i.e. the "graph")"""
        self.set_field_data(*args)
        for k in kwargs:
            if self._field_attrs.has_key(k):
                if isinstance(self._field_attrs[k], (type(kwargs[k]), type(None))):
                    self._field_attrs[k] = kwargs[k]

    def configure_axes_attrs(self, **kwargs):

#        # Do nothing when axes are fixed
#        if self._axes_attrs['fixed_view'] == 'on': 
#            return 
#        # Hold is 'on', check several cases
#        elif self._figure_attrs['hold'] == 'on':
#            x = self._field_attrs['x']
#            y = self._field_attrs['y']
#            self._axes_attrs['xmin'] = min(min(x), self._axes_attrs['xmin'])
#            self._axes_attrs['xmax'] = max(max(x), self._axes_attrs['xmax'])
#            self._axes_attrs['ymin'] = min(min(y), self._axes_attrs['ymin'])
#            self._axes_attrs['ymax'] = max(max(y), self._axes_attrs['ymax'])
#
#        # Hold is 'off', we are saved
#        else:
#            x = self._field_attrs['x']
#            y = self._field_attrs['y']
#
#            if not (x and y):
#                return
#            
#            self._axes_attrs['xmin'] = min(x)
#            self._axes_attrs['xmax'] = max(x)
#            self._axes_attrs['ymin'] = min(y)
#            self._axes_attrs['ymax'] = max(y)

        for k in kwargs:
            if self._axes_attrs.has_key(k):
                if isinstance(self._axes_attrs[k], (type(kwargs[k]), type(None))):
                    self._axes_attrs[k] = kwargs[k]

    def set_figure_axes(self):
        # First plot case
        _fa = self._field_attrs # Short name
        _aa = self._axes_attrs  # Short name
        if _aa['axmin'] == None:
            _aa['axmin'] = min(_fa['x'])
            _aa['axmax'] = max(_fa['x'])
            _aa['aymin'] = min(_fa['y'])
            _aa['aymax'] = max(_fa['y'])
        else: # Old plot exists, find new axes range
            _aa['axmin'] = min(_aa['axmin'], min(_fa['x']))
            _aa['axmax'] = max(_aa['axmax'], max(_fa['x']))
            _aa['aymin'] = min(_aa['aymin'], min(_fa['y']))
            _aa['aymax'] = max(_aa['aymax'], max(_fa['y']))

        # Always use the user specified ranges to override
        if not _aa['xmin'] == None: _aa['axmin'] = _aa['xmin']
        if not _aa['xmax'] == None: _aa['axmax'] = _aa['xmax']
        if not _aa['ymin'] == None: _aa['aymin'] = _aa['ymin']
        if not _aa['ymax'] == None: _aa['aymax'] = _aa['ymax']

    def configure_figure_attrs(self, **kwargs):
        """Configure parameters for the figure"""
        for k in kwargs:
            if self._figure_attrs.has_key(k):
                if isinstance(self._figure_attrs[k], (type(kwargs[k]), type(None))):
                    self._figure_attrs[k] = kwargs[k]
                elif k == 'hardcopy': # Special case; can be str, None, bool, int
                    self._figure_attrs[k] = kwargs[k] # No type checking

    def set(self, **kwargs):
        self.configure_field_attrs( **kwargs)
        self.configure_axes_attrs(**kwargs)
        self.configure_figure_attrs(**kwargs)
        if self._figure_attrs['movie'] == 'on': 
            self._axes_attrs['fixed_view'] = 'on'
        else: 
            self._axes_attrs['fixed_view'] = 'off'

    
    def _clear_plot(self):
        self._set_defaults()

    def _next_frame_name(self, filename):
        (base, ext) = _path.splitext(filename)
        if not ext: ext='.eps'
        newfilename = "%s%.4d%s" % (base, self._figure_attrs['frame_count'], ext) 
        self._figure_attrs['frame_count'] += 1
        return newfilename

    def _get_file_type(self, filename):
        (base, ext) = _path.splitext(filename)
        if not ext: ext = '.eps'
        return ext


    def view(self, xmin, xmax, ymin, ymax):
        _aa = self._axes_attrs
        _aa['xmin'] = xmin
        _aa['xmax'] = xmax
        _aa['ymin'] = ymin
        _aa['ymax'] = ymax

    def saveField(self):
        self.fieldDataList.append(FieldData(self))

       

class VizCurveGnuplot(VizCurveBase):

    def __init__(self, *args, **kwargs):
        self.g = Gnuplot.Gnuplot(persist=1)

        self.d = []
        self.colors = {'black':    -1, 'red':   1, # lt
                       'green':     2, 'blue':  3,
                       'purple':    4, 'aqua':  5,
                       'brown':     6, 'orange':7,
                       'light brown':8}
# http://sparky.rice.edu/~hartigan/gnuplot.html
        self.point_types = { 'diamond': 1, '+': 2, # pt
                            'square':   3, 'x': 4,
                            'triangle': 5, '*': 6,
                            'solid': 'solid'}

        VizCurveBase.__init__(self, *args, **kwargs)
        self.init()


    def init(self):
        self.d = []
        if self._figure_attrs['movie'] == 'off': 
            self.gnuplot_color = None
            self.gnuplot_line_type = None


    def _clear_plot(self):
        if self._figure_attrs['movie'] == 'off': 
            VizCurveBase._clear_plot(self)
        self.init()

    def set(self, **kwargs):
        """Configure gnuplot specific settings"""
        if kwargs.has_key('hold'): self._figure_attrs['hold'] = kwargs['hold']
        if self._figure_attrs['hold'] == 'off':
            self._clear_plot()

        VizCurveBase.set(self, **kwargs)
        self.gnuplot_color = self.colors[self._field_attrs['color']]
        self.gnuplot_line_type = self.point_types[self._field_attrs['line_type']]
        self.g.xlabel(self._figure_attrs['xlabel'])
        self.g.ylabel(self._figure_attrs['ylabel'])
        self.g.title(self._figure_attrs['title'])
        scale = self._figure_attrs['scale']
        if scale == 'loglog':
            self.g('set logscale x')
            self.g('set logscale y')
            self.g('set autoscale')
        elif scale == 'logx':
            self.g('set logscale x')
            self.g('set nologscale y')
            self.g('set autoscale')
        elif scale == 'logy':
            self.g('set logscale y')
            self.g('set nologscale x')
            self.g('set autoscale')
        elif scale == 'linear':
            self.g('set nologscale y')
            self.g('set nologscale x')
            self.g('set nologscale x')
            self.g('set noautoscale')

    def update(self):
        self.g.plot(*self.d)

    def plot(self,*args,**kwargs):
        """ Plot the figure with Gnuplot."""
        self.set_field_data(*args)
        # Save the plot for later use
        self.saveField()
        self.set_figure_axes()

        self.g('set xrange[%g:%g]' % (self._axes_attrs['axmin'], self._axes_attrs['axmax']))
        self.g('set yrange[%g:%g]' % (self._axes_attrs['aymin'], self._axes_attrs['aymax']))

        if self.gnuplot_line_type == 'solid':
            withstring = "lines lt %d" % (self.gnuplot_color)
        else:
            withstring="linespoints lt %d pt %d" % (self.gnuplot_color, self.gnuplot_line_type)
        self.d.append(Gnuplot.Data(self._field_attrs['x'], 
                      self._field_attrs['y'], with=withstring, 
                      title=self._field_attrs['legend']))
        if self._figure_attrs['screen_plot']:
            self.g.plot(*self.d)
        if self._figure_attrs['hardcopy']:
            if isinstance(self._figure_attrs['hardcopy'], str): 
                filename = self._figure_attrs['hardcopy']
            else:
                filename = self._figure_attrs['default_filename']
            if self._figure_attrs['movie'] == 'on':
                filename = self._next_frame_name(filename)
            file_type = self._get_file_type(filename)
            self.g('set output "%s"' % (filename))
            type_str = ""
            if file_type == '.eps': type_str = 'postscript eps color'
            elif file_type == '.png': type_str = 'png'
            elif file_type == '.tex': type_str = 'latex' 
            elif file_type == '.fig': type_str = 'fig' 
            else: type_str = 'postscript eps color'
            self.g('set terminal %s' % (type_str))
            self.g.plot(*self.d)
            self.g('set terminal x11')


class VizCurvePyx(VizCurveBase):

    def __init__(self, *args, **kwargs):
        VizCurveBase.__init__(self, *args, **kwargs)
        self.init()

    def init(self):
        self.d = []
        self.line = []
        self.pyx_line_types = {
#            'solid': pyxgraph.line(), 
#            'solid':    pyxgraph.symbol(None, lineattrs=[pyxcanvas.linestyle.solid]), 
            'diamond':  pyxgraph.symbol.diamond, 
            '+':        pyxgraph.symbol.plus, 
            'square':   pyxgraph.symbol.square, 
            'x':        pyxgraph.symbol.cross, 
            'triangle': pyxgraph.symbol.triangle, 
            '*':        pyxgraph.symbol.circle}

        self.pyx_line_colors = {
            'black':    pyxgraph.color.cmyk.Black,
            'red':      pyxgraph.color.cmyk.Red,
            'green':    pyxgraph.color.cmyk.Green,
            'blue':     pyxgraph.color.cmyk.Blue,
            'purple':   pyxgraph.color.cmyk.Purple,
            'aqua':     pyxgraph.color.cmyk.Aquamarine,
            'brown':    pyxgraph.color.cmyk.Brown,
            'orange':   pyxgraph.color.cmyk.Orange,
            'light brown': pyxgraph.color.cmyk.BrickRed,
            }
    def _clear_plot(self):
        if self._figure_attrs['movie'] == 'off': 
            VizCurveBase._clear_plot(self)
        self.init()

    def set(self, **kwargs):
        if kwargs.has_key('hold'): self._figure_attrs['hold'] = kwargs['hold']
        if self._figure_attrs['hold'] == 'off': self._clear_plot()
        VizCurveBase.set(self, **kwargs)
        la = [pyxcanvas.linestyle.solid, \
            self.pyx_line_colors[self._field_attrs['color']]]
        if not self._field_attrs['line_type'] == 'solid':
            lt = self.pyx_line_types[self._field_attrs['line_type']]
            sa = [pyxcanvas.stroked(self.pyx_line_colors['black']),
                pyxcanvas.filled(self.pyx_line_colors[self._field_attrs['color']]),
                pyxcanvas.linewidth.thin]
            self.line.append(pyxgraph.symbol(lt, symbolattrs=sa, lineattrs=la))
        else:
            self.line.append(pyxgraph.line(lineattrs=la))

    def plot(self, *args, **kwargs):
        self.set_field_data(*args)
        # Save the plot for later use
        self.saveField()
        self.set_figure_axes()
        p = pyxgraph.axispainter(titledist="0.1 cm", labeldist="0.2 cm", titledirection=None)
        self.pyx_x = pyxgraph.linaxis(min=self._axes_attrs['axmin'],
                                      max=self._axes_attrs['axmax'],
                                      title=self._figure_attrs['xlabel'],
                                      painter=p)
        self.pyx_y = pyxgraph.linaxis(min=self._axes_attrs['aymin'],
                                      max=self._axes_attrs['aymax'],
                                      title=self._figure_attrs['ylabel'], 
                                      painter=p)

        self.d.append(_concatenate((self._field_attrs['x'][:, _NewAxis], self._field_attrs['y'][:, _NewAxis]), 1).tolist())
        h = 5; w = 10
        if self._field_attrs['legend']:
            self.g = pyxgraph.graphxy( width=w, x=self.pyx_x,
                y=self.pyx_y, key=pyxgraph.key(pos="tr"))
        else:
            self.g = pyxgraph.graphxy( width=w, x=self.pyx_x, y=self.pyx_y)

        for i in xrange(len(self.d)):
            if self._field_attrs['legend']:
                self.g.plot(pyxgraph.data(pyxdata.data(self.d[i]), x=0, y=1, 
                            title=self._field_attrs['legend']), style=self.line[i])
            else:
                self.g.plot(pyxgraph.data(pyxdata.data(self.d[i]), x=0, y=1), style=self.line[i])
        # Find out what filename to save to
            if self._figure_attrs['hardcopy'] and isinstance(self._figure_attrs['hardcopy'], str):
                filename = self._figure_attrs['hardcopy']
            else:
                filename = self._figure_attrs['default_filename']
#        [x, y] = self.g.pos(0.5*(self._axes_attrs['axmin']+self._axes_attrs['axmax']), 1.05*self._axes_attrs['aymax'])
        [x, y] = self.g.pos(self._axes_attrs['axmin'], 1.05*self._axes_attrs['aymax'])
        self.g.text(x, y, self._figure_attrs['title'] )
        self.g.dodata(); self.g.writetofile(filename)



#    def clear_plot(self):
#        """Reset Plot"""
#        self.d=[]
#        VizCurveBase.clear_plot(self)
#        self.fixedView = False

def plot(*args,**kwargs):
    """Matlabish plotting command. Factory method for choosing the underlying
    plot engine."""
    if kwargs.has_key('program'): program = kwargs['program']
    else: program = 'Gnuplot'
    if program == 'Gnuplot':
        set_plot_engine(VizCurveGnuplot)
    elif program == 'Pyx':
        set_plot_engine(VizCurvePyx)
    elif program == 'BLT':
        print "Sorry, not implemented, using Gnuplot instead"
        set_plot_engine(VizCurveGnuplot)
    else:
        print "No such plotting engine, using Gnuplot instead"
        set_plot_engine(VizCurveGnuplot)

    if kwargs.has_key('figure'): figure(kwargs['figure'])
    if not plot_utils.plot_exists: 
        figure()
    plot_utils.plotter.set(**kwargs)
    plot_utils.plotter.plot(*args, **kwargs)
    return plot_utils.plotter

def view(xmin, xmax, ymin, ymax):
    plot_utils.plotter.view(xmin, xmax, ymin, ymax)
    plot_utils.plotter.update()
 
def semilogx(*args, **kwargs):
    kwargs['scale'] = 'logx'
    plot(*args, **kwargs)
 
def semilogy(*args, **kwargs):
    kwargs['scale'] = 'logy'
    plot(*args, **kwargs)

def loglog(*args, **kwargs):
    kwargs['scale'] = 'loglog'
    plot(*args, **kwargs)

def linear(*args, **kwargs):
    kwargs['scale'] = 'linear'
    plot(*args, **kwargs)

def update():
    plot_utils.plotter.update()

def clear_plot():
    """ Broken"""
    global plot_exists
#    plot_exists = 0

# set Gnuplot as default plot-engine
set_plot_engine(VizCurveGnuplot)
