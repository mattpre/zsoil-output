#-------------------------------------------------------------------------------
# Name:        tut-plate
# Purpose:     this tutorial describes on how to extract results for a plate
#              then we will plot bending moments, shear forces for selected time instances
#              then will make an envelope of bending moments
#
# Author:      Andrzej Truty
#
# Created:     2020
# Copyright:   (c) ZACE Services Ltd 2020
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
    print(value)
    return

#===================================
def make_color_contour_M (mesh,his,rcf,solution_index,sel_elements,comp,out_png,title,
                          ArmerWood_flag = False,fiber_flag=None):
#===================================
    shell_results = Shell_EleResults(mesh, his, rcf)
    sorted_sel_ele = sel_elements [:]
    indices_of_sorted_sel_ele = mesh.get_ele_indices_for_sel_elements(sorted_sel_ele)

    if ArmerWood_flag:
       M_comp = shell_results.get_moment_AW_for_sel_elements (indices_of_sorted_sel_ele,solution_index,comp,fiber_flag)
    else:
        M_comp = shell_results.get_moment_for_sel_elements (indices_of_sorted_sel_ele,solution_index,comp)

    x = []
    y = []
    M = []

    for i, element in enumerate(sorted_sel_ele):
        ngaus = element.get_ngaus()
        #ele_coord = e.get_ele_coord ()
        for igaus in range(ngaus):
            xyz = element.get_xyz_gp (igaus)
            x.append (xyz [0])
            y.append (xyz [2])
            M.append (M_comp [i,igaus])

    xmin = min (x)
    xmax = max (x)
    ymin = min (y)
    ymax = max (y)

    fig    = pylab.figure (figsize=(8,8),dpi=80,facecolor='white') #RO
    axes   = fig.add_subplot(1,1,1)
    #pyplot.subplots_adjust(left=0.17, bottom=0.1, right=0.95, top=0.85, hspace = 0.05) #RO

    xi = numpy.linspace(xmin,xmax,100)
    yi = numpy.linspace(ymin,ymax,100)
    # grid the data.
    zi_max = griddata((x, y), M, (xi[None,:], yi[:,None]), method='cubic')
    # contour the gridded data, plotting dots at the randomly spaced data points.
    CS = pylab.contour(xi,yi,zi_max,20,colors='k')
    CS = pylab.contourf(xi,yi,zi_max,20,cmap=pylab.cm.jet)
    pylab.colorbar() # draw colorbar
    # plot data points.
    #pylab.scatter(x,y,marker='o',c='b',s=1)
    pylab.xlim (xmax,xmin)
    pylab.ylim (ymin,ymax)
    pylab.title(title)
    fig.savefig (out_png)

#===================================
def make_color_contour_envelope_M (mesh,his,rcf,solution_indices,sel_elements,comp,out_png,title,
                                 make_min=True,ArmerWood_flag = False,fiber_flag=None):
#===================================
    shell_results = Shell_EleResults(mesh, his, rcf)
    sorted_sel_ele = sel_elements [:]
    indices_of_sorted_sel_ele = mesh.get_ele_indices_for_sel_elements(sorted_sel_ele)

    if ArmerWood_flag:
        M_comp_min,M_comp_max = shell_results.get_moment_AW_envelope_for_sel_elements (indices_of_sorted_sel_ele,solution_indices,comp,fiber_flag)
        comp_in_all = Shell_EleResults.comps_all_AW.index ('M'+comp)
    else:
        M_comp_min,M_comp_max = shell_results.get_moment_envelope_for_sel_elements (indices_of_sorted_sel_ele,solution_indices,comp)
        comp_in_all = Shell_EleResults.comps_all.index ('M'+comp)

    x = []
    y = []
    M = []

    xmin =  1.0e99
    xmax = -1.0e99
    ymin =  1.0e99
    ymax = -1.0e99

    for i, element in enumerate(sorted_sel_ele):
        ngaus = element.get_ngaus()
        ele_coord = element.get_ele_coord ()
        for xyz in ele_coord:
            xmin = min(xmin,xyz [0])
            xmax = max(xmax,xyz [0])
            ymin = min(ymin,xyz [2])
            ymax = max(ymax,xyz [2])

        for igaus in range(ngaus):
            xyz = element.get_xyz_gp (igaus)
            x.append (xyz [0])
            y.append (xyz [2])
            if make_min:
                M.append (M_comp_min [i,igaus,comp_in_all])
            else:
                M.append (M_comp_max [i,igaus,comp_in_all])

    # xmin = min (x)
    # xmax = max (x)
    # ymin = min (y)
    # ymax = max (y)

    fig    = pylab.figure (figsize=(8,8),dpi=80,facecolor='white') #RO
    axes   = fig.add_subplot(1,1,1)
    #pyplot.subplots_adjust(left=0.17, bottom=0.1, right=0.95, top=0.85, hspace = 0.05) #RO

    xi = numpy.linspace(xmin,xmax,100)
    yi = numpy.linspace(ymin,ymax,100)
    # grid the data.
    zi_max = griddata((x, y), M, (xi[None,:], yi[:,None]), method='linear')
    # contour the gridded data, plotting dots at the randomly spaced data points.
    CS = pylab.contourf(xi,yi,zi_max,20,cmap=pylab.cm.jet,)
    pylab.colorbar() # draw colorbar
    # plot data points.
    #pylab.scatter(x,y,marker='o',c='b',s=1)
    pylab.xlim (xmax,xmin)
    pylab.ylim (ymin,ymax)
    pylab.title(title)
    fig.savefig (out_png)

#===================================
def main():
#===================================
    plot_dir = '../tutorials/plots-tut-plate/'
    rsl_directory = '../tests/plate/'
    project = rsl_directory + 'plate-2-dense'
    #create an instance of  history of execution object
    his     = HistoryOfExecution (project)
    #read mesh from *.dat file
    mesh    = Mesh     (project)
    #read configuration file on how the element and nodal results are stored in external binary files
    rcf     = RCF_info (project)

    #empty zoom filter
    zoom = []
    #user vector to set local shell x-axis
    user_vec = [1.0,0.0,0.0]
    #colect all time instances for which we have converged results (use time filter if needed)
    list_of_time_instances = his.give_converged_time_solutions (time_filter =[-1.0,1.0e99])

    zoom_filter = [mesh.ZOOM_FILTER_ENOUGH_ONE_NODE_IS_IN, zoom]
    sel_elements = mesh.get_list_of_elements(Element.GROUP_SHELL, \
                                         0.0, False, mat_filter=[], \
                                         ele_class_filter=[], zoom_filter=zoom_filter)
    # set current local bases in all shell elements for the given vector
    # by doing this we can get T_VEC or T_TNS transformation matrices for the vector
    # and for the plane stress tensor (moments and membrane forces)
    for element in sel_elements:
        element.set_transf_matrices(user_vec)

    #plot M-XX bending moments
    make_color_contour_M (mesh,his,rcf,list_of_time_instances[-1],sel_elements,"XX",\
                          plot_dir+'plate-MXX.png',\
                          "Moment-XX")

    #plot M-XY moments
    make_color_contour_M(mesh,his,rcf,list_of_time_instances[-1],sel_elements,"XY", \
                         plot_dir+'plate-MXY.png', \
                         "Moment-XY")

    # plot M-YY bending moments for top (last parameter is + which is used for the top fibers)
    make_color_contour_M(mesh, his, rcf, list_of_time_instances[-1], sel_elements, "XX", \
                         plot_dir + 'plate-MXX-AW-top.png', \
                        "Moment-XX(A-W)(top)",True,'+')

    # plot M-XX Wood-Armer bending moments bottom (-)
    make_color_contour_M(mesh, his, rcf, list_of_time_instances[-1], sel_elements, "XX", \
                         plot_dir + 'plate-MXX-AW-bot.png', \
                         "Moment-XX(A-W)(bot)", True, '-')

    # plot M-YY Wood-Armer bending moments top(+)
    make_color_contour_M(mesh, his, rcf, list_of_time_instances[-1], sel_elements, "YY", \
                         plot_dir + 'plate-MYY-AW-top.png', \
                         "Moment-YY(A-W)(top)", True, '+')

    # plot M-YY Wood-Armer bending moments top(-)
    make_color_contour_M(mesh, his, rcf, list_of_time_instances[-1], sel_elements, "YY", \
                         plot_dir + 'plate-MYY-AW-bot.png', \
                         "Moment-YY(A-W)(bot)", True, '-')

    # plot M-XX envelope
    make_color_contour_envelope_M (mesh, his, rcf, list_of_time_instances, sel_elements, "XX", \
                                   plot_dir + 'plate-MXX-envelope.png', \
                                   "Moment-XX-envelope")

    # plot M-XX Wood Armer envelope bottom
    make_color_contour_envelope_M(mesh, his, rcf, list_of_time_instances, sel_elements, "XX", \
                                  plot_dir + 'plate-MXX-envelope-bot.png', \
                                  "Moment-XX-envelope-bot",True,True,'-')

    # plot M-XX Wood Armer envelope top
    make_color_contour_envelope_M(mesh, his, rcf, list_of_time_instances, sel_elements, "XX", \
                                  plot_dir + 'plate-MXX-envelope-top.png', \
                                  "Moment-XX-envelope-top", True, True, '+')

    #show the result on the screen
    pylab.show()
    pylab.close()

#here we call main routine that does the job
if __name__ == '__main__':
    main()

