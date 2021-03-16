"""Common functionality for all types of plotting.

"""

plot_exists = 0
fig_dict = {}
fig_count = -1
_plot_engine = lambda : None
plotter = None

from time import sleep

def set_plot_engine(_pe):
    global _plot_engine
    #print "_pe is: %s" % (_pe)
    _plot_engine = _pe

def figure(*args):
    """Change the active plot figure. When called without arguments, a new figure
    is made with a default name (incremented integer name). If called with a
    non-exisiting figure name, make a new figure with this name. Otherwise,
    activate the named figure."""
    global plotter, fig_dict, fig_count, plot_exists, _plot_engine
    plot_exists = 1

    print "plot_engine is %s" % (_plot_engine)

    if len(args) == 0:
        plotter = _plot_engine()
        fig_count += 1
        fig_dict[fig_count] = plotter
        sleep(0.1) # Must sleep or else Gnuplot fails.
        return plotter
    try:
        plotter = fig_dict[args[0]]
    except: 
        plotter = _plot_engine()
        if not plotter == None:
            fig_dict[args[0]] = plotter
            sleep(0.1) # Must sleep or else Gnuplot fails.
        else:
            print "You must provide a plot engine!"
    return plotter


def get_figure_names():
    """Return a list of defined figure names."""
    return fig_dict.keys()
