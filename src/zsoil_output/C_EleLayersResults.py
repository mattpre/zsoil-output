import os
import numpy
from .C_RslFileCache import RslFileCache


#=====================================================
class EleLayersResults ():
#=====================================================

    #====================================================
    def __init__ (self,mesh,his,rcf,rcf_lay):
    #=====================================================
        self.mesh = mesh
        self.his  = his
        self.rcf  = rcf
        self.rcf_lay = rcf_lay
        self.rsl_files = RslFileCache ()

    #======================================
    def get_layers_results (self,group,solution_indices,mat_filter = [],zoom_filter = [0,[[],[],[]]]):
    #======================================
        mesh = self.mesh
        his  = self.his
        rcf  = self.rcf
        rcf_lay = self.rcf_lay

        out = []

        sel_elements = mesh.get_list_of_elements (group,0.0,False,mat_filter,[],zoom_filter)
        if len(sel_elements) <= 0:
            return out

        ext    = sel_elements [0].rsl_lay_ext
        if ext == None:
            return out

        sum_of_layers_in_rsl_file = mesh.get_sum_of_layers (group)
        lay_size_int = rcf_lay.give_one_lay_size (group)
        lay_size= lay_size_int * 4 #in bytes

        # record index of every (element, gauss point, layer) at solution 0;
        # element seek positions are multiples of one layer record
        lay_rows = []
        lay_ids  = []
        for e in sel_elements:
            ngaus = len(e.xsiGP)
            nlayers_all, nlayers_reinf = e.get_nlayers()
            first = e.seek_in_lay_rsl // lay_size
            for igaus in range (ngaus):
                for layer in range (nlayers_all):
                    lay_rows.append (first + igaus * nlayers_all + layer)
                    lay_ids.append ((e.index,igaus+1,layer+1))
        if len(lay_rows) <= 0:
            return out

        sol_offs = numpy.asarray (solution_indices,dtype=numpy.int64) * sum_of_layers_in_rsl_file
        rec_indices = numpy.asarray (lay_rows,dtype=numpy.int64)[:,None] + sol_offs[None,:]

        data = self.rsl_files.read_records (mesh.my_project+ext,lay_size_int,rec_indices.ravel())
        data = data.reshape (len(lay_rows),len(solution_indices)*lay_size_int)

        for k in range (len(lay_rows)):
            lay_rsl = list (lay_ids[k])
            lay_rsl.extend (data [k])
            out.append (lay_rsl)
        return out

    #======================================
    def close (self):
    #======================================
        # releases the result files kept open for fast repeated reads
        self.rsl_files.close ()
