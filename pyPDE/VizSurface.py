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


""" Surface plots in 3D, commonly referred to as 3D plots, but we call it 2D
plots... """

import Gnuplot
import sys,os
if os.environ.get('NUMPYARRAY','') == 'Numeric' or sys.modules.has_key('Numeric'):
    from Numeric import *
else:
    from numarray import *
from os import path
import threading
from time import sleep
from os import system as s
from pyPDE import amin,amax


_plot_exists=0

def manual():
    print "Write me!"

class VizSurfaceBase(object):
    def __init__(self, *args, **kwargs):
        self._set_defaults()

        self.configure_field_attrs(*args, **kwargs)
        self.configure_axes_attrs(**kwargs)
        self.configure_figure_attrs(**kwargs)

    def _set_defaults(self):
        self._figure_attrs = {
            'xlabel':   "x",
            'ylabel':   "y",
            'zlabel':   "z",
            'title':    "z(x,y)",
#            'scale':    "linear",
            'filename': None,
            'movie':     "off",
            'frame_count': 0,
            'hardcopy': False,
            'default_filename': "tmp.eps",
            'screen_plot': True #,
#            'hold': 'off'
        }

        self._axes_attrs = {
# Axes used for plotting:
            'axmin':     None,
            'axmax':     None,
            'aymin':     None,
            'aymax':     None,
            'azmin':     None,
            'azmax':     None,
# User specified axes ranges (if given):
            'xmin':     None,
            'xmax':     None,
            'ymin':     None,
            'ymax':     None,
            'zmin':     None,
            'zmax':     None,
            'fixed_view': False
        }

        self._field_attrs = {
            'legend':   "",
            'color':    "blue",
            'line_type': "solid",
            'line_width': 0.3,
            'color_map': "def1",
            'x':        None,
            'y':        None,
            'z':        None
        }

    def set_field_data(self, *args):
        if len(args) == 3:
            self._field_attrs['x'] = args[0]
            self._field_attrs['y'] = args[1]
            self._field_attrs['z'] = args[2]
        elif len(args) == 1:
            (xlen, ylen) = args[0].shape
            self._field_attrs['x'] = arange(xlen)
            self._field_attrs['y'] = arange(ylen)
            self._field_attrs['z'] = args[0]
#        else:
#            raise ValueError, "Error: 1 or 3 args must be supplied"

    def configure_field_attrs(self, *args, **kwargs):
        """Configure parameters for the field (i.e. the "graph")"""
        self.set_field_data(*args)
        for k in kwargs:
            if self._field_attrs.has_key(k):
                self._field_attrs[k] = kwargs[k]

    def configure_axes_attrs(self, **kwargs):
        for k in kwargs:
            if self._axes_attrs.has_key(k):
                self._axes_attrs[k] = kwargs[k]

    def set_figure_axes(self):
        # First plot case
        _fa = self._field_attrs # Short name
        _aa = self._axes_attrs # Short name
        if _aa['axmin'] == None:
            _aa['axmin'] = min(_fa['x'])
            _aa['axmax'] = max(_fa['x'])
            _aa['aymin'] = min(_fa['y'])
            _aa['aymax'] = max(_fa['y'])
            _aa['azmin'] = amin(_fa['z'])
            _aa['azmax'] = amax(_fa['z'])
        else: # Old plot exists, find new axes range
            _aa['axmin'] = min(_aa['axmin'],min(_fa['x']))
            _aa['axmax'] = max(_aa['axmax'],max(_fa['x']))
            _aa['aymin'] = min(_aa['aymin'],min(_fa['y']))
            _aa['aymax'] = max(_aa['aymax'],max(_fa['y']))
            _aa['azmin'] = min(_aa['azmin'],amin(_fa['z']))
            _aa['azmax'] = max(_aa['azmax'],amax(_fa['z']))

        # Always use the user specified ranges to override
        if not _aa['xmin'] == None: _aa['axmin'] = _aa['xmin']
        if not _aa['xmax'] == None: _aa['axmax'] = _aa['xmax']
        if not _aa['ymin'] == None: _aa['aymin'] = _aa['ymin']
        if not _aa['ymax'] == None: _aa['aymax'] = _aa['ymax']
        if not _aa['zmin'] == None: _aa['azmin'] = _aa['zmin']
        if not _aa['zmax'] == None: _aa['azmax'] = _aa['zmax']

    def configure_figure_attrs(self, **kwargs):
        """Configure parameters for the figure"""
        for k in kwargs:
            if self._figure_attrs.has_key(k):
                self._figure_attrs[k] = kwargs[k]

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
        (base, ext) = path.splitext(filename)
        if not ext: ext='.eps'
        newfilename = "%s%.4d%s" % (base, self._figure_attrs['frame_count'], ext) 
        self._figure_attrs['frame_count'] +=1
        return newfilename

    def _get_file_type(self, filename):
        (base, ext) = path.splitext(filename)
        if not ext: ext='.eps'
        return ext


    def view(self, xmin, xmax, ymin, ymax):
        _aa = self._axes_attrs
        _aa['xmin'] = xmin
        _aa['xmax'] = xmax
        _aa['ymin'] = ymin
        _aa['ymax'] = ymax
        _aa['zmin'] = zmin
        _aa['zmax'] = zmax
       

class VizSurfaceGnuplot(VizSurfaceBase):

    def __init__(self,*args,**kwargs):
        self.g = Gnuplot.Gnuplot(persist=1)

        self.d = []
        self.colors = {'black': -1, 'red': 1, # lt
                       'green': 2, 'blue':3,
                       'purple': 4, 'aqua': 5,
                       'brown':6, 'orange':7,
                       'light brown':8}
# http://sparky.rice.edu/~hartigan/gnuplot.html
        self.point_types = { 'diamond': 1, '+': 2, # pt
                            'square': 3, 'x': 4,
                            'triangle': 5, '*': 6,
                            'solid': 'solid'}

        self.gnuplot_color_map = {'def1': '22,-13,-27',
                                  'def2': '4,4,4'}
        VizSurfaceBase.__init__(self, *args, **kwargs)
        self.init()
        self.gnuplot_color = self.colors[self._field_attrs['color']]
        self.gnuplot_line_type = self.point_types[self._field_attrs['line_type']]
#        self.gnuplot_line_line_width = self.point_types[self._field_attrs['line_width']]


    def init(self):
        self.d = []
        if self._figure_attrs['movie'] == 'off': 
            self.gnuplotColor = None
            self.gnuplotLineType = None


    def _clear_plot(self):
        if self._figure_attrs['movie'] == 'off': 
            VizSurfaceBase._clear_plot(self)
        self.init()

    def set(self, **kwargs):
        """Configure gnuplot specific settings"""
#        if kwargs.has_key('hold'): self._figure_attrs['hold'] = kwargs['hold']
#        if self._figure_attrs['hold'] == 'off':
        self._clear_plot()
        VizSurfaceBase.set(self, **kwargs)
        self.gnuplot_color = self.colors[self._field_attrs['color']]
        self.gnuplot_line_type = self.point_types[self._field_attrs['line_type']]
        self.g.xlabel(self._figure_attrs['xlabel'])
        self.g.ylabel(self._figure_attrs['ylabel'])
#        self.g.zlabel(self._figure_attrs['zlabel'])
        self.g.title(self._figure_attrs['title'])
#        self.g('set palette rgbformulae %s' % (self.gnuplot_color_map[self._field_attrs['color_map']]))

#        scale = self._figure_attrs['scale']
#        if scale == 'loglog':
#            self.g('set logscale x')
#            self.g('set logscale y')
#            self.g('set autoscale')
#        elif scale == 'logx':
#            self.g('set logscale x')
#            self.g('set nologscale y')
#            self.g('set autoscale')
#        elif scale == 'logy':
#            self.g('set logscale y')
#            self.g('set nologscale x')
#            self.g('set autoscale')
#        elif scale == 'linear':
#            self.g('set nologscale y')
#            self.g('set nologscale x')
#            self.g('set nologscale x')
#            self.g('set noautoscale')
#


    def update(self):
        self.g.plot(*self.d)


    def plot(self,*args,**kwargs):
        """ Plot the figure with Gnuplot."""

        self.set_field_data(*args)
        self.set_figure_axes()

        self.g('set xrange[%g:%g]' % (self._axes_attrs['axmin'], self._axes_attrs['axmax']))
        self.g('set yrange[%g:%g]' % (self._axes_attrs['aymin'], self._axes_attrs['aymax']))
        self.g('set zrange[%g:%g]' % (self._axes_attrs['azmin'], self._axes_attrs['azmax']))

        if self.gnuplot_line_type == 'solid': withstring = "lines"
        else: withstring="linespoints pt %d" % (self.gnuplot_line_type)

        withstring += " lt %d lw %g" % ( self.gnuplot_color, self._field_attrs['line_width'] )
        self.d.append(Gnuplot.GridData(self._field_attrs['z'], self._field_attrs['x'], self._field_attrs['y'], with=withstring, title=self._field_attrs['legend']))

        if self._figure_attrs['screen_plot']:
#            self.g('set pm3d')
            self.g('set contour')
            self.g('set nokey')
            self.g.splot(*self.d)

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
            else: type_str = 'postscript eps'
            self.g('set terminal %s' % (type_str))
            self.g.splot(*self.d)
#            sleep(1)
#            s("cp %s %s" % (filename, filename+'~'))
#            s("cp %s %s" % (filename+'~', filename))
#            s("rm %s" % (filename+'~'))
            self.g('set terminal x11')
#

#    def clear_plot(self):
#        """Reset Plot"""
#        self.d=[]
#        VizCurveBase.clear_plot(self)
        
#        self.fixedView = False

def plot(*args,**kwargs):
    """Matlabish plotting command."""
    global _g
    global _plot_exists
    if not _plot_exists:
        _plot_exists = 1
        _g = VizSurfaceGnuplot(*args, **kwargs)
    _g.set(**kwargs)
    _g.plot(*args, **kwargs)
    return _g

def view(xmin, xmax, ymin, ymax):
    _g.view(xmin, xmax, ymin, ymax)
    _g.update()
 
def update():
    _g.update()

def clear_plot():
    global _plot_exists
#    _plot_exists = 0

