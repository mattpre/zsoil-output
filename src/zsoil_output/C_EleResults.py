import os
import numpy
import numpy as np
from .C_Element import *
from .C_RslFileCache import RslFileCache


#=====================================================
class EleResults ():
#=====================================================

    BUFFER_F_HANDLE    = 0
    BUFFER_SUM_ALL_GP  = 1
    BUFFER_GP_SIZE_INT = 2
    BUFFER_MAX         = 3

    #====================================================
    def __init__ (self,mesh,his,rcf):
    #=====================================================
        self.mesh = mesh
        self.his  = his
        self.rcf  = rcf

        self.rsl_files = RslFileCache ()
        self.group_buffers = []
        for i in range (Element.GROUP_MAX):
            buffer = []
            for j in range (self.BUFFER_MAX):
                buffer.append (None)
            self.group_buffers.append (buffer)


    #======================================
    def get_element_results_ex (self,sel_elements,solution_indices):
    #======================================
	# all elements must belong to the same group

        mesh = self.mesh
        his  = self.his
        rcf  = self.rcf

        out = []

        if len(sel_elements) <= 0:
            return out

        ext    = sel_elements [0].rsl_ext
        group  = sel_elements [0].group
        group_index = sel_elements[0].group_dict_inv [group]

        if self.group_buffers [group_index][self.BUFFER_F_HANDLE] == None:
            sum_of_all_gp_in_rsl_file = mesh.get_sum_of_gauss_points (group)
            gp_size_int = rcf.give_one_gp_size(group)
            self.group_buffers [group_index][self.BUFFER_F_HANDLE   ] = mesh.my_project+ext
            self.group_buffers [group_index][self.BUFFER_SUM_ALL_GP ] = sum_of_all_gp_in_rsl_file
            self.group_buffers [group_index][self.BUFFER_GP_SIZE_INT] = gp_size_int
        path                      = self.group_buffers [group_index][self.BUFFER_F_HANDLE   ]
        sum_of_all_gp_in_rsl_file = self.group_buffers [group_index][self.BUFFER_SUM_ALL_GP ]
        gp_size_int               = self.group_buffers [group_index][self.BUFFER_GP_SIZE_INT]

        gp_size = gp_size_int * 4 #in bytes

        # record index of every (element, gauss point) at solution 0;
        # element seek positions are multiples of one gauss point record
        gp_rows = []
        gp_ids  = []
        for e in sel_elements:
            first = e.seek_in_rsl // gp_size
            for igaus in range (len(e.xsiGP)):
                gp_rows.append (first + igaus)
                gp_ids.append ((e.index,igaus+1))
        sol_offs = np.asarray (solution_indices,dtype=np.int64) * sum_of_all_gp_in_rsl_file
        rec_indices = np.asarray (gp_rows,dtype=np.int64)[:,None] + sol_offs[None,:]

        data = self.rsl_files.read_records (path,gp_size_int,rec_indices.ravel())
        data = data.reshape (len(gp_rows),len(solution_indices)*gp_size_int)

        for k in range (len(gp_rows)):
            gauss_rsl = [gp_ids[k][0],gp_ids[k][1]]
            gauss_rsl.extend (data [k])
            out.append (gauss_rsl)
        return out

    #======================================
    def close (self):
    #======================================
        # releases the result files kept open for fast repeated reads
        self.rsl_files.close ()

    #======================================
    def get_element_results (self,group,solution_indices,mat_filter = [],zoom_filter = [0,[[],[],[]]]):
    #======================================
        mesh = self.mesh
        his  = self.his
        rcf  = self.rcf

        out = []

        sel_elements = mesh.get_list_of_elements (group,0.0,False,mat_filter=mat_filter,\
                                                  ele_class_filter=[],zoom_filter=zoom_filter)
        print(len(sel_elements))
        if len(sel_elements) <= 0:
            return out
        out = self.get_element_results_ex (sel_elements,solution_indices)
        return out
