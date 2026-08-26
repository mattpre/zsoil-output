#-------------------------------------------------------------------------------
# Name:        tut-simple-beams-flex-2D-MNT
# Purpose:     this tutorial describes on how to extract results from flexibility based beam elements
#              bending moments, shear forces and axial forces are plot for selected time instances
#              note that these elements have 5 gauss points (the 1-st and 5-th coincide with element nodes)
#
# Author:      Andrzej Truty
#
# Created:     2020
# Copyright:   (c) user 2020
# Licence:     <your licence>
#-------------------------------------------------------------------------------

import os
import numpy
import sys
# those who have installed miktex environment can deactivate the next statement to have nice
# Latex style axes descriptions, equations etc...
#rc('text',usetex=True)

# Only needed if you want to use a local clone of zsoil_output instead of the pip-installed package:
# import sys
# sys.path.insert(0, '..\\src')   # the folder that CONTAINS zsoil_output\

from zsoil_output import Mesh, HistoryOfExecution, RCF_info, Beam_S_EleResults, Element



import pylab
from matplotlib.widgets import Slider, Button, RadioButtons

#===================================
def update_column (value):
#===================================
    print(value)
    return


#===================================
def make_MNT_for_column (indices_of_sorted_sel_ele,list_of_time_instances,his,mesh,rcf,out_png):
#===================================
    #create an instance of object ==> results
    beams_rsl = Beam_S_EleResults (mesh,his,rcf)
    #quit if there is nothing to do
    if len(list_of_time_instances) <= 0:
        print("null list of time steps -> nothing to do")
        return

    tmp = len(indices_of_sorted_sel_ele)
    #create array of vertical coordinates of beam nodal points
    y_tab = []
    for i in indices_of_sorted_sel_ele:
        e = mesh.elements [i-1]
        ngaus = e.get_ngaus()
        for igaus in range (ngaus):
            xyz = e.get_xyz_gp (igaus)
            y_tab.append (xyz [1])
    y_tab = np.array (y_tab)
    #y_tab.sort()


    #let's make a diagram  of bending moments
    # Two subplots, unpack the axes array immediately
    fig, (axes_Nx,axes_Qy,axes_Mz) = pylab.subplots (1, 3, sharey=True,figsize=(16,8))
    fig.subplots_adjust(hspace=0.0,wspace=0.1,top=0.95, right=0.95)

    #axcolor = 'lightgoldenrodyellow'
    #axamp   = pylab.axes([0.15, 0.96, 0.75, 0.02], axisbg='red')
    #samp    = Slider (axamp, 'Column',1,15,valinit=1)
    #samp.on_changed (update_column)

    colors     = ['b','r','c','m','k']
    linestyles_dashed = ['--','--','--','--','--']
    linestyles_solid  = ['-','-','-','-','-']
    linestyles = []
    if (len(list_of_time_instances) < len (linestyles_solid)):
        linestyles = linestyles + linestyles_solid
        linestyles = linestyles + linestyles_dashed
    else:
        linestyles = linestyles + linestyles_dashed
        linestyles = linestyles + linestyles_solid

    while (len(list_of_time_instances) > len (linestyles)):
        linestyles = linestyles + linestyles

    while (len(list_of_time_instances) > len (colors)):
        colors = colors + colors

    #make loop over time instances
    for i,solution_index in enumerate(list_of_time_instances):

        row = his.data [solution_index]
        time = row [his.DATA_TIME]

        N = beams_rsl.get_force_for_sel_elements(indices_of_sorted_sel_ele, \
                                                 solution_index,'X')
        N = N.reshape (N.size)
        Q = beams_rsl.get_force_for_sel_elements(indices_of_sorted_sel_ele, \
                                                 solution_index,'Y')
        Q = Q.reshape (Q.size)
        M = beams_rsl.get_moment_for_sel_elements(indices_of_sorted_sel_ele, \
                                               solution_index,'Z')
        M = M.reshape (M.size)

        #check the Y limit for the vertical axis
        ylim = (min(y_tab),max(y_tab))

        #make grid on all diagrams and put labels
        axes_Nx.grid (True,which="both",ls="-",lw=0.4)
        axes_Nx.plot (N,y_tab,color=colors [i],lw=1.75,linestyle=linestyles [i],label=str(time))
        #check the Y limit for the vertical axis
        #ylim = axes_Nx.get_ylim()
        axes_Nx.plot ([0.0,0.0],ylim,color='black',lw=1.5)
        axes_Nx.set_xlabel ('N [kN]',fontsize=18)
        axes_Nx.legend(loc='best')

        axes_Qy.grid (True,which="both",ls="-",lw=0.4)
        axes_Qy.plot (Q,y_tab,color=colors [i],lw=1.75,linestyle=linestyles [i],label=str(time))
        axes_Qy.plot ([0.0,0.0],ylim,color='black',lw=1.5)
        axes_Qy.set_xlabel ('Q [kN]',fontsize=18)
        axes_Qy.legend(loc='best')

        axes_Mz.grid (True,which="both",ls="-",lw=0.4)
        axes_Mz.plot (M*(-1.0),y_tab,color=colors [i],lw=1.75,linestyle=linestyles [i],label=str(time))
        axes_Mz.plot ([0.0,0.0],ylim,color='black',lw=1.5)
        axes_Mz.set_xlabel ('M [kNm]',fontsize=18)
        axes_Mz.legend(loc='best')

    #put x coordinate of a selected column on the plot
    #x_coord_str = '% 6.2f' % x_coord
    #pylab.title ("Column at x = "+x_coord_str)
    fig.subplots_adjust(hspace=0)

    #save the figure as external file in png format
    fig.savefig (out_png)


#===================================
def main():
#===================================
    project =     project = '..\\tests\\beams\\BEAM-STD-2D-1a'
    plot_dir = '..\\tutorials\\plots-tut-simple-beams-flex-2D-MNT\\'

    #create an instance of  history of execution object
    his     = HistoryOfExecution (project)
    #read mesh from *.dat file
    mesh    = Mesh     (project)
    #read configuration file on how the element and nodal results are stored in external binary files
    rcf     = RCF_info (project)

    #colect all time instances for which converged (!) results are stored (use time filter if needed)
    list_of_time_instances = his.give_converged_time_solutions (time_filter = [-1.0,1.0e38])
    #
    mat_filter = [] #put user material indices if needed
    #one can set also a zoom in x,y,z directions to extract only some set of beams
    zoom_filter = [mesh.ZOOM_FILTER_ENOUGH_ONE_NODE_IS_IN,\
                  [[], [], []]]
    # get beam elements with given filters
    sel_elements = mesh.get_list_of_elements(Element.GROUP_BEAM, \
                                              0.0, False, mat_filter=mat_filter, \
                                              ele_class_filter=[], zoom_filter=zoom_filter)
    # sort elements in given direction for plotting purposes (use direction that coincides with the beams local x axis direction)
    point = [0.0, 0.0, 0.0]
    dir   = [0.0, 1.0, 0.0]
    #get sorted element's list as list of indices rather than list of elements (objects)
    sorted_sel_ele = mesh.sort_sel_elements_by_dist_along_dir         (sel_elements,point,dir)
    indices_of_sorted_sel_ele = mesh.get_ele_indices_for_sel_elements (sorted_sel_ele)

    out_png =plot_dir+ "MNT-column.png"
    #plot results N, Q, M
    make_MNT_for_column (indices_of_sorted_sel_ele, list_of_time_instances, his, mesh, rcf, \
                         out_png)
    #show the result on the screen
    pylab.show()
    pylab.close()


#here we call main routine that does the job
if __name__ == '__main__':
    main()

