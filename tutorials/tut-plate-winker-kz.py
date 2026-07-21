#-------------------------------------------------------------------------------
# Name:        tutorial-2-plate
# Purpose:     this tutorial describes on how to extract results for a plate
#              then we will plot bending moments, shear forces for selected time instances
#              then will make an envelope of bending moments
#
# Author:      Andrzej Truty
#
# Created:     20-08-2015
# Copyright:   (c) user 2015
# Licence:     <your licence>
#-------------------------------------------------------------------------------

import os
import sys
# those who have installed miktex environment can deactivate the next statement to have nice
# Latex style axes descriptions, equations etc...
#rc('text',usetex=True)

sys.path.append ('..\\')

from C_Mesh import *
from C_HistoryOfExecution import *
from C_Rcf_info import *
from C_ShellResults import *

import numpy

import pylab
from scipy.interpolate import griddata

#===================================
def update_column (value):
#===================================
    print (value)
    return

#===================================
def make_color_contour_kz (mesh,his,rcf,solution_index_ref,solution_index_curr, \
                           sel_shells,adj_continuum,title,out_png):
#===================================
    nodal_results     = NodalResults (mesh, his, rcf)
    continuum_results = Continuum_EleResults(mesh, his, rcf, nodal_results)

    sorted_adj_continuum = adj_continuum [:]
    indices_of_sorted_adj_continuum = mesh.get_ele_indices_for_sel_elements (sorted_adj_continuum)

    if len(sel_shells) <= 0 or len(adj_continuum) <= 0:
        print ("shells list size",len(sel_shells))
        print ("continuum list size",len(adj_continuum))
        print ("error: one of the lists is empty")
        return

    x = np.zeros (len(sel_shells))
    y = np.zeros (len(sel_shells))

    if solution_index_ref != None:
        sigYY_tot_ref = continuum_results.get_rsl_for_sel_elements_ex (indices_of_sorted_adj_continuum, \
                                          solution_index_ref,Continuum_EleResults.TSTRESS_KEY_EXTRA,'YY')
        uY_ref = np.zeros (len(sel_shells))
        for i, shell in enumerate(sel_shells):
            uY_at_nodes = nodal_results.get_nodes_displacements (shell.nodes,solution_index_ref,'Y')
            uY_value = 0.0
            for k in range (uY_at_nodes.size):
                uY_value = uY_value + uY_at_nodes [k]
            uY_value = uY_value / uY_at_nodes.size
            uY_ref [i] = uY_value

    if solution_index_curr != None:
        sigYY_tot_curr = continuum_results.get_rsl_for_sel_elements_ex (indices_of_sorted_adj_continuum, \
                                          solution_index_curr,Continuum_EleResults.TSTRESS_KEY_EXTRA,'YY')
        uY_curr = np.zeros (len(sel_shells))
        for i, shell in enumerate(sel_shells):
            uY_at_nodes = nodal_results.get_nodes_displacements (shell.nodes,solution_index_curr,'Y')
            uY_value = 0.0
            for k in range (uY_at_nodes.size):
                uY_value = uY_value + uY_at_nodes [k]
            uY_value = uY_value / uY_at_nodes.size
            uY_curr [i] = uY_value
            xyz = shell.get_center()
            x [i] = xyz [0]
            y [i] = xyz [2]

    kz = np.zeros (x.size)
    for i in range (x.size):
        dUy = uY_curr [i]
        if solution_index_ref != None:
            dUy = dUy - uY_ref [i]
        dSigY = sigYY_tot_curr [i]
        if solution_index_ref != None:
            dSigY = dSigY - sigYY_tot_ref [i]
        if abs(dUy) < 1.0e-14:
            dUy = copysign (1.0e-14,dUy)
        kz [i] = dSigY / dUy

    xmin = min (x)
    xmax = max (x)
    ymin = min (y)
    ymax = max (y)

    fig    = pylab.figure (figsize=(8,8),dpi=80,facecolor='white') #RO
    axes   = fig.add_subplot(1,1,1)
    #pyplot.subplots_adjust(left=0.17, bottom=0.1, right=0.95, top=0.85, hspace = 0.05) #RO
    #
    xi = numpy.linspace(xmin,xmax,50)
    yi = numpy.linspace(ymin,ymax,50)
    # # grid the data.
    zi_max = griddata((x, y), kz, (xi[None,:], yi[:,None]), method='linear')
    # # contour the gridded data, plotting dots at the randomly spaced data points.
    CS = pylab.contour(xi,yi,zi_max,20,colors='k',)
    CS = pylab.contourf(xi,yi,zi_max,20,cmap=pylab.cm.jet)
    pylab.colorbar() # draw colorbar
    # # plot data points.
    # #pylab.scatter(x,y,marker='o',c='b',s=1)
    pylab.xlim (xmax,xmin)
    pylab.ylim (ymin,ymax)
    pylab.title(title)
    fig.savefig (out_png)


#===================================
def main():
#===================================
    plot_dir = '../tutorials/plots-tut-plate-winker-kz/'
    rsl_directory = '../tests/kz-plate/'
    project = rsl_directory + 'winkler'

    #create an instance of  history of execution object
    his     = HistoryOfExecution (project)
    #read mesh from *.dat file
    mesh    = Mesh     (project)
    #read configuration file on how the element and nodal results are stored in external binary files
    rcf     = RCF_info (project)

    #create a graph in which we collect elements of class continuum/shells attached to a given node
    mesh.set_adj_elements_to_nodes ([Element.GROUP_CONTINUUM,Element.GROUP_SHELL])

    #empty zoom
    zoom = []
    # user defined vector to set up local shell x-axis
    user_vec = [1.0,0.0,0.0]
    #colect all time instances for which we have stored results and  steps were converged with given time filtering
    list_of_time_instances = his.give_converged_time_solutions (time_filter =[-1.0,1.0e99])

    #take shell material with index 2 only
    mat_shells  = [2]
    zoom_filter = [mesh.ZOOM_FILTER_ENOUGH_ONE_NODE_IS_IN, zoom]
    sel_shells  = mesh.get_list_of_elements(Element.GROUP_SHELL, \
                                         0.0, False, mat_filter=mat_shells, \
                                         ele_class_filter=[], zoom_filter=zoom_filter)
    # set current local bases in all shell elements for the given vector
    # by doing this we can get T_VEC or T_TNS transformation matrices for the vector
    # and for the plane stress tensor (moments and membrane forces)
    for element in sel_shells:
        element.set_transf_matrices(user_vec)

    # collect continuum elements adjacent to selected shells
    adj_continuum = mesh.get_adj_elements_to (sel_shells,Element.GROUP_CONTINUUM,None)

    # kz - secant for step 3 (indexing from zero !)
    make_color_contour_kz (mesh, his, rcf, None, list_of_time_instances [3], \
                           sel_shells,adj_continuum,"kz-0-3",plot_dir+"kz-0-3.png")

    # kz - incremental between steps 2 and 3
    make_color_contour_kz(mesh, his, rcf, list_of_time_instances[1], list_of_time_instances[2], \
                          sel_shells, adj_continuum, "kz-1-2",plot_dir+"kz-1-2.png")

    # kz - incremental between steps 3 and 4
    make_color_contour_kz(mesh, his, rcf, list_of_time_instances[2], list_of_time_instances[3], \
                        sel_shells, adj_continuum, "kz-2-3",plot_dir+"kz-2-3.png")

    #
    # #show the result on the screen
    pylab.show()
    pylab.close()

#here we call main routine that does the job
if __name__ == '__main__':
    main()

