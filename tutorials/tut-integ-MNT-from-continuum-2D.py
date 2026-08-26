#-------------------------------------------------------------------------------
# Name:        tut-integ-MNT-from-continuum-2D
# Purpose:     compute and plot MNT diagrams for a beam discretized using Q4 elements (continuum for structures)
#
# Author:      Andrzej Truty
#
# Created:     2020
# Copyright:   (c) ZACE Services Ltd 2020
# Licence:     <your licence>
#-------------------------------------------------------------------------------
import matplotlib.pyplot  as pyplot
import matplotlib.pylab   as pylab
from math import *
import numpy as np
import xlwt
import sys,os

#miktex is needed to use the following package
from matplotlib import rc
# those who have installed miktex environment can deactivate the next statement to have nice
# Latex style axes descriptions, equations etc...
#rc('text',usetex=True)


#path to SDK
# Only needed if you want to use a local clone of zsoil_output instead of the pip-installed package:
# import sys
# sys.path.insert(0, '..\\src')   # the folder that CONTAINS zsoil_output\

from zsoil_output import Mesh, HistoryOfExecution, RCF_info, Continuum_EleResults, Diagram, Element
from zsoil_output.C_ContinuumResults import (postprocess_MNQ_time_history_for_fict_uniform_beam,
                                              smooth_jumps_in_stress_resultant_time_history,
                                              get_MNQ_envelopes_for_fict_uniform_beam)


# ==========================================================
def beam_MNT_from_continuum_2D (zsoil_data, out_png):
# ==========================================================
    # goal:
    # this function will plot bending moment integrated from stresses in continuum elements 2D
    #
    # zsoil_data     - name of the project (without file extension !)
    # real_time      - time at which the
    # out_png        - name of the resulting png file
    # ------------------------------------------------------------------

    # ------------------------------------------------------------------
    #next 3 calls must always be done to be able to analyze ZSOIL results
    # (1) launch history of execution for this project
    his = HistoryOfExecution (data_file)

    # (2) read mesh from *.dat file
    mesh = Mesh (data_file)

    # (3) read configuration file on how the element and nodal results are stored in external binary files
    rcf = RCF_info (data_file)
    # ------------------------------------------------------------------

    #take all time instances
    solution_indices = his.give_converged_time_solutions()

    # create instance of diagram class which will contain plot contact stresses
    dg = Diagram()
    dg.set_figure_size((12, 12))
    markers = dg.get_markers(len(solution_indices))
    colors = dg.get_colors(len(solution_indices))

    zoom = [[], [], []]
    zoom_filter = [mesh.ZOOM_FILTER_ENOUGH_ONE_NODE_IS_IN, zoom]

    cont_rsl = Continuum_EleResults(mesh, his, rcf)
    mat_filter = [1]
    sel_elements = mesh.get_list_of_elements(Element.GROUP_CONTINUUM, \
                                             0.0, False, mat_filter=mat_filter, \
                                             ele_class_filter=[], zoom_filter=zoom_filter)

    x_dir = np.array([1.0, 0.0, 0.0]) # direction of beam X-axis
    pt    = np.array([0.0, 0.5, 0.0]) # starting central point of the beam cross section
    y_dir = np.array([0.0, 1.0, 0.0]) # local beam Y-axis

    #here we will get {NX, QY, QZ, MX, MY, MZ} at center of each elements layer along given drection and for each time instance --> array [n_ele_layrs,n_time_instances,6]
    statusOK, x_tab, x_size, Forces = cont_rsl.get_MNQ_time_history_for_fict_uniform_beam( \
        sel_elements, solution_indices, \
        pt, x_dir, y_dir)

    # postprocess the above result to have it along full beam length including starting and end cross sections
    x_tab_ex, Forces_ex = postprocess_MNQ_time_history_for_fict_uniform_beam(Forces, x_tab, x_size)

    # smooth bending moment diagram Mz
    smooth_jumps_in_stress_resultant_time_history(Forces_ex, 5)

    label = 'Mz'
    i = 0
    dg.add_plot(x_tab_ex, Forces_ex [:,0,5], marker_=markers[i], markevery_=1, markersize_=0, \
            markerfacecolor_=colors[i], markeredgecolor_='black', \
            linestyle_='-', label_=label, color_=colors[i])

    markers = dg.get_markers(3)
    colors = dg.get_colors  (3)

    # make envelope of bending moment MZ as a leading force
    env_min, env_max = get_MNQ_envelopes_for_fict_uniform_beam(Forces_ex, 5)  # MZ envelope

    i = 1
    label = 'Mz-env-min'
    dg.add_plot(x_tab_ex, env_min[:,5], marker_=markers[i], markevery_=1, markersize_=3, \
            markerfacecolor_=colors[i], markeredgecolor_='black', \
            linestyle_='--', label_=label, color_=colors[i])

    i = 2
    label = 'Mz-env-max'
    dg.add_plot(x_tab_ex, env_max[:, 5], marker_=markers[i], markevery_=1, markersize_=5, \
            markerfacecolor_=colors[i], markeredgecolor_='black', \
            linestyle_='', label_=label, color_=colors[i])


    # dg.set_x_axes_tick_params (labelsize=16)
    # dg.set_x_min   ( 30.0)
    # dg.set_x_max   (120.0)
    # dg.set_y_max   (  0.0)
    # dg.set_y_min   (-20.0)
    # #dg.set_x_log_scale ()

    dg.set_x_label (r'$x$ [m]')
    dg.set_y_label (r'$Mz$ [kNm]')
    dg.add_grid         ()
    dg.add_legend       ()
    dg.save_fig(out_png)


if __name__ == '__main__':
    data_file = "..\\tests\continuum-MNT\\beam-2D"
    plots_dir = "..\\tutorials\\plots-tut-integ-MNT-from-continuum-2D\\"

    beam_MNT_from_continuum_2D (data_file, plots_dir+"beam_MNT_from_continuum_2D.png")

    pyplot.show ()
