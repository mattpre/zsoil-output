import os
import numpy
import numpy as np
from   math import *
from .C_RslFileCache import RslFileCache


#=====================================================
class NodalResults ():
#=====================================================

    #====================================================
    def __init__ (self,mesh,his,rcf):
    #=====================================================
        self.mesh = mesh
        self.his  = his
        self.rcf  = rcf

        self.dict_exts          = {"NODAL_RSL":".s00","NODAL_RES":".r00","NODAL_VEL":".v00","NODAL_ACC":".a00"}
        self.dict_node_size_int = {"NODAL_RSL":rcf.give_one_node_size(),"NODAL_RES":rcf.give_one_node_size(),"NODAL_VEL":rcf.give_one_node_size_v(),"NODAL_ACC":rcf.give_one_node_size_a()}
        self.rsl_files = RslFileCache ()

    #======================================
    def get_nodal_results (self,selected_nodes,solution_indices,dict_key):
    #======================================
        mesh = self.mesh
        his  = self.his
        rcf  = self.rcf

        out = []

        if len(selected_nodes) <= 0:
            return out

        node_size_int= self.dict_node_size_int[dict_key]
        node_indices_1 = [node.index for node in selected_nodes]
        data = self._read_nodal_records (node_indices_1,solution_indices,dict_key)
        data = data.reshape (len(selected_nodes),len(solution_indices)*node_size_int)

        for k in range (len(selected_nodes)):
            node_rsl = [node_indices_1[k]]
            node_rsl.extend (data [k])
            out.append (node_rsl)
        return out

    #======================================
    def _read_nodal_records (self,node_indices_1,solution_indices,dict_key):
    #======================================
        # returns float32 array of shape (len(node_indices_1),len(solution_indices),node_size_int)
        rsl_file_ext = self.dict_exts [dict_key]
        node_size_int= self.dict_node_size_int[dict_key]
        n_nodes = len(self.mesh.nodes)
        rows = np.asarray (node_indices_1,dtype=np.int64) - 1
        sol_offs = np.asarray (solution_indices,dtype=np.int64) * n_nodes
        rec_indices = rows[:,None] + sol_offs[None,:]
        data = self.rsl_files.read_records (self.mesh.my_project+rsl_file_ext,node_size_int,rec_indices.ravel())
        return data.reshape (len(node_indices_1),len(solution_indices),node_size_int)

    #======================================
    def close (self):
    #======================================
        # releases the result files kept open for fast repeated reads
        self.rsl_files.close ()



    #======================================
    def _get_node_vectorial_rsl_ths (self,node_index_1,solution_indices,dict_key,rsl_key,comp=None):
    #======================================
        node = self.mesh.get_node(node_index_1)
        selected_nodes = [node]
        res_all = self.get_nodal_results(selected_nodes, solution_indices, dict_key)
        vec = res_all[0]  # only one node is considered here
        size = len(solution_indices)

        offs = 1 # as the first position is occupied by the node number in the vector
        node_size_int = self.dict_node_size_int [dict_key]

        if comp == None:
            col_1, n_comps = self.rcf.give_col_for_nod_rsl(rsl_key, 'X')
            col_2, n_comps = self.rcf.give_col_for_nod_rsl(rsl_key, 'Y')
            col_3, n_comps = self.rcf.give_col_for_nod_rsl(rsl_key, 'Z')
            out = np.zeros(size * 3).reshape (size,3)
        else:
            col_1, ncomps = self.rcf.give_col_for_nod_rsl(rsl_key, comp)
            col_2 = -1
            col_3 = -1
            out = np.zeros(size)

        for i in range(len(solution_indices)):
            if comp == None:
                if col_1 != -1:
                    out[i][0] = vec[col_1+offs]
                if col_2 != -1:
                    out[i][1] = vec[col_2+offs]
                if col_3 != -1:
                    out[i][2] = vec[col_3+offs]
            else:
                if col_1 != -1:
                    out[i] = vec[col_1+offs]

            offs = offs + node_size_int

        return out


    # ======================================
    def _get_node_scalar_rsl_ths(self, node_index_1, solution_indices, dict_key, rsl_key):
    # ======================================
        node = self.mesh.get_node(node_index_1)
        selected_nodes = [node]
        res_all = self.get_nodal_results(selected_nodes, solution_indices, dict_key )
        vec = res_all[0]  # only one node is considered here
        size = len(solution_indices)


        offs = 1  # as the first position is occupied by the node number in the vector
        node_size_int = self.dict_node_size_int [dict_key]

        col, n_comps = self.rcf.give_col_for_nod_rsl(rsl_key, '')
        out = np.zeros(size)

        for i in range(len(solution_indices)):
            if col != -1:
                out[i] = vec[col+offs]
            offs = offs + node_size_int

        return out


    # ======================================
    def _get_nodes_vectorial_rsl (self, node_indices_1, solution_index, dict_key, rsl_key, comp=None):
    # ======================================
        # get displacements for list of nodes at one time instance
        # node_index_1 starts from 1
        # solution indices are verified in class HistoryOfExecution, from there import time instances or other quantities
        # comp is one of 'X', 'Y', 'Z'; if None take full vector
        # this functions returns numpy array of size (nr_time_instances,3) when comp==None
        # this functions returns numpy array of size (nr_time_instances) when comp != None
        if comp == None:
            ret = np.zeros(len(node_indices_1) * 3).reshape (len(node_indices_1),3)
            cols = [self.rcf.give_col_for_nod_rsl(rsl_key, c)[0] for c in ['X', 'Y', 'Z']]
        else:
            ret = np.zeros(len(node_indices_1))
            cols = [self.rcf.give_col_for_nod_rsl(rsl_key, comp)[0]]
        if len(node_indices_1) <= 0:
            return ret

        node_indices = [self.mesh.get_node(i).index for i in node_indices_1]
        data = self._read_nodal_records (node_indices,[solution_index],dict_key)[:,0,:]

        for k in range(len(cols)):
            if cols[k] != -1:
                if comp == None:
                    ret[:,k] = data[:,cols[k]]
                else:
                    ret[:] = data[:,cols[k]]

        return ret


    # ======================================
    def _get_nodes_scalar_rsl (self, node_indices_1, solution_index,dict_key, rsl_key):
    # ======================================
        # get pore pressures for list of nodes at one time instance
        # node_index_1 starts from 1
        # solution indices are verified in class HistoryOfExecution, from there import time instances or other quantities
        # this functions returns numpy array of size (nr_time_instances)
        ret = np.zeros(len(node_indices_1))
        col, n_comps = self.rcf.give_col_for_nod_rsl(rsl_key, '')
        if len(node_indices_1) <= 0 or col == -1:
            return ret

        node_indices = [self.mesh.get_node(i).index for i in node_indices_1]
        data = self._read_nodal_records (node_indices,[solution_index],dict_key)[:,0,:]
        ret[:] = data[:,col]

        return ret

    # ---------------------------------
    # public  functions
    #functions to be used by the user
    # ---------------------------------

    #======================================
    def get_node_displacements_time_history (self,node_index_1,solution_indices,comp=None):
    #======================================
        #get displacement for one node but several time instances
        #node_index_1 is a node index and it starts from 1
        # solution indices is a Python list; these can be verified in class HistoryOfExecution
        #comp is one of 'X', 'Y', 'Z'; if None ==> get full vector [3]
        #this functions returns numpy array of size (nr_time_instances,3) when comp==None
        # this functions returns numpy array of size (nr_time_instances) when comp != None
        return self._get_node_vectorial_rsl_ths (node_index_1,solution_indices,'NODAL_RSL','DISP_TRA',comp)

    #======================================
    def get_node_reaction_forces_time_history (self,node_index_1,solution_indices,comp=None):
    #======================================
        #get displacement for one node but several time instances
        #node_index_1 is a node index and it starts from 1
        # solution indices is a Python list; these can be verified in class HistoryOfExecution
        #comp is one of 'X', 'Y', 'Z'; if None ==> get full vector [3]
        #this functions returns numpy array of size (nr_time_instances,3) when comp==None
        # this functions returns numpy array of size (nr_time_instances) when comp != None
        return self._get_node_vectorial_rsl_ths (node_index_1,solution_indices,'NODAL_RES','DISP_TRA',comp)

    #======================================
    def get_node_veloc_time_history (self,node_index_1,solution_indices,comp=None):
    #======================================
        #get velocity for one node but several time instances
        #node_index_1 is a node index and it starts from 1
        # solution indices is a Python list; these can be verified in class HistoryOfExecution
        #comp is one of 'X', 'Y', 'Z'; if None ==> get full vector [3]
        #this functions returns numpy array of size (nr_time_instances,3) when comp==None
        # this functions returns numpy array of size (nr_time_instances) when comp != None
        return self._get_node_vectorial_rsl_ths (node_index_1,solution_indices,'NODAL_VEL','DISP_TRA',comp)

    #======================================
    def get_node_accel_time_history (self,node_index_1,solution_indices,comp=None):
    #======================================
        #get accelerations for one node but several time instances
        #node_index_1 is a node index and it starts from 1
        # solution indices is a Python list; these can be verified in class HistoryOfExecution
        #comp is one of 'X', 'Y', 'Z'; if None ==> get full vector [3]
        #this functions returns numpy array of size (nr_time_instances,3) when comp==None
        # this functions returns numpy array of size (nr_time_instances) when comp != None
        return self._get_node_vectorial_rsl_ths (node_index_1,solution_indices,'NODAL_ACC','DISP_TRA',comp)

    #======================================
    def get_nodes_displacements (self,node_indices_1,solution_index,comp=None):
    #======================================
        # get displacements for list of nodes at one time instance
        #node_index_1 starts from 1
        #solution indices are verified in class HistoryOfExecution, from there import time instances or other quantities
        #comp is one of 'X', 'Y', 'Z'; if None take full vector
        #this functions returns numpy array of size (nr_time_instances,3) when comp==None
        # this functions returns numpy array of size (nr_time_instances) when comp != None
        return self._get_nodes_vectorial_rsl (node_indices_1,solution_index,'NODAL_RSL','DISP_TRA',comp)

    #======================================
    def get_nodes_reaction_forces (self,node_indices_1,solution_index,comp=None):
    #======================================
        # get displacements for list of nodes at one time instance
        #node_index_1 starts from 1
        #solution indices are verified in class HistoryOfExecution, from there import time instances or other quantities
        #comp is one of 'X', 'Y', 'Z'; if None take full vector
        #this functions returns numpy array of size (nr_time_instances,3) when comp==None
        # this functions returns numpy array of size (nr_time_instances) when comp != None
        return self._get_nodes_vectorial_rsl (node_indices_1,solution_index,'NODAL_RES','DISP_TRA',comp)

    #======================================
    def get_nodes_veloc (self,node_indices_1,solution_index,comp=None):
    #======================================
        # get velocities for list of nodes at one time instance
        #node_index_1 starts from 1
        #solution indices are verified in class HistoryOfExecution, from there import time instances or other quantities
        #comp is one of 'X', 'Y', 'Z'; if None take full vector
        #this functions returns numpy array of size (nr_time_instances,3) when comp==None
        # this functions returns numpy array of size (nr_time_instances) when comp != None
        return self._get_nodes_vectorial_rsl (node_indices_1,solution_index,'NODAL_VEL','DISP_TRA',comp)


    #======================================
    def get_nodes_accel (self,node_indices_1,solution_index,comp=None):
    #======================================
        # get accelerations for list of nodes at one time instance
        #node_index_1 starts from 1
        #solution indices are verified in class HistoryOfExecution, from there import time instances or other quantities
        #comp is one of 'X', 'Y', 'Z'; if None take full vector
        #this functions returns numpy array of size (nr_time_instances,3) when comp==None
        # this functions returns numpy array of size (nr_time_instances) when comp != None
        return self._get_nodes_vectorial_rsl (node_indices_1,solution_index,'NODAL_ACC','DISP_TRA',comp)



    # ======================================
    def get_node_rotations_time_history (self, node_index_1, solution_indices, comp=None):
    # ======================================
        # get displacement for one node but several time instances
        # node_index_1 is a node index and it starts from 1
        # solution indices is a Python list; these can be verified in class HistoryOfExecution
        # comp is one of 'X', 'Y', 'Z'; if None ==> get full vector [3]
        # this functions returns numpy array of size (nr_time_instances,3) when comp==None
        # this functions returns numpy array of size (nr_time_instances) when comp != None
        return self._get_node_vectorial_rsl_ths(node_index_1, solution_indices,'NODAL_RSL','DISP_ROT', comp)

    # ======================================
    def get_node_reaction_moments_time_history (self, node_index_1, solution_indices, comp=None):
    # ======================================
        # get displacement for one node but several time instances
        # node_index_1 is a node index and it starts from 1
        # solution indices is a Python list; these can be verified in class HistoryOfExecution
        # comp is one of 'X', 'Y', 'Z'; if None ==> get full vector [3]
        # this functions returns numpy array of size (nr_time_instances,3) when comp==None
        # this functions returns numpy array of size (nr_time_instances) when comp != None
        return self._get_node_vectorial_rsl_ths(node_index_1, solution_indices,'NODAL_RES','DISP_ROT', comp)


    # ======================================
    def get_nodes_rotations (self, node_indices_1, solution_index, comp=None):
    # ======================================
        # get rotations for list of nodes at one time instance
        # node_index_1 starts from 1
        # solution indices are verified in class HistoryOfExecution, from there import time instances or other quantities
        # comp is one of 'X', 'Y', 'Z'; if None take full vector
        # this functions returns numpy array of size (nr_time_instances,3) when comp==None
        # this functions returns numpy array of size (nr_time_instances) when comp != None
        return self._get_nodes_vectorial_rsl(node_indices_1, solution_index,'NODAL_RSL','DISP_ROT', comp)

    # ======================================
    def get_nodes_reaction_moments (self, node_indices_1, solution_index, comp=None):
    # ======================================
        # get rotations for list of nodes at one time instance
        # node_index_1 starts from 1
        # solution indices are verified in class HistoryOfExecution, from there import time instances or other quantities
        # comp is one of 'X', 'Y', 'Z'; if None take full vector
        # this functions returns numpy array of size (nr_time_instances,3) when comp==None
        # this functions returns numpy array of size (nr_time_instances) when comp != None
        return self._get_nodes_vectorial_rsl(node_indices_1, solution_index,'NODAL_RES','DISP_ROT', comp)



    # ======================================
    def get_node_pore_pressures_time_history (self, node_index_1, solution_indices):
    # ======================================
        # get displacement for one node but several time instances
        # node_index_1 is a node index and it starts from 1
        # solution indices is a Python list; these can be verified in class HistoryOfExecution
        # this functions returns numpy array of size (nr_time_instances)
        return self._get_node_scalar_rsl_ths(node_index_1, solution_indices,'NODAL_RSL','PPRESS')


    # ======================================
    def get_nodes_pore_pressures(self, node_indices_1, solution_index):
    # ======================================
        # get pore pressures for list of nodes at one time instance
        # node_index_1 starts from 1
        # solution indices are verified in class HistoryOfExecution, from there import time instances or other quantities
        # this functions returns numpy array of size (nr_time_instances)
        return  self._get_nodes_scalar_rsl (node_indices_1, solution_index,'NODAL_RSL','PPRESS')

    # ======================================
    def get_node_pressure_heads_time_history (self, node_index_1, solution_indices):
    # ======================================
        # get displacement for one node but several time instances
        # node_index_1 is a node index and it starts from 1
        # solution indices is a Python list; these can be verified in class HistoryOfExecution
        # this functions returns numpy array of size (nr_time_instances)
        return self._get_node_scalar_rsl_ths(node_index_1, solution_indices,'NODAL_RSL','PRES_HEAD')


    # ======================================
    def get_nodes_pressure_heads(self, node_indices_1, solution_index):
    # ======================================
        # get pressure heads for list of nodes at one time instance
        # node_index_1 starts from 1
        # solution indices are verified in class HistoryOfExecution, from there import time instances or other quantities
        # this functions returns numpy array of size (nr_time_instances)
        return self._get_nodes_scalar_rsl(node_indices_1, solution_index,'NODAL_RSL','PRES_HEAD')

    # ======================================
    def get_node_temperatures_time_history (self, node_index_1, solution_indices):
    # ======================================
        # get displacement for one node but several time instances
        # node_index_1 is a node index and it starts from 1
        # solution indices is a Python list; these can be verified in class HistoryOfExecution
        # this functions returns numpy array of size (nr_time_instances)
        return self._get_node_scalar_rsl_ths(node_index_1, solution_indices,'NODAL_RSL','TEMP')


    # ======================================
    def get_nodes_temperatures(self, node_indices_1, solution_index):
    # ======================================
        # get temperatures for list of nodes at one time instance
        # node_index_1 starts from 1
        # solution indices are verified in class HistoryOfExecution, from there import time instances or other quantities
        # this functions returns numpy array of size (nr_time_instances)
        return self._get_nodes_scalar_rsl(node_indices_1, solution_index,'NODAL_RSL','TEMP')

    # ======================================
    def get_node_humidities_time_history (self, node_index_1, solution_indices):
    # ======================================
        # get displacement for one node but several time instances
        # node_index_1 is a node index and it starts from 1
        # solution indices is a Python list; these can be verified in class HistoryOfExecution
        # this functions returns numpy array of size (nr_time_instances)
        return self._get_node_scalar_rsl_ths(node_index_1, solution_indices,'NODAL_RSL','HUMID')


    # ======================================
    def get_nodes_humidities(self, node_indices_1, solution_index):
    # ======================================
        # get humidities for list of nodes at one time instance
        # node_index_1 starts from 1
        # solution indices are verified in class HistoryOfExecution, from there import time instances or other quantities
        # this functions returns numpy array of size (nr_time_instances)
        return self._get_nodes_scalar_rsl (node_indices_1, solution_index,'NODAL_RSL','HUMID')


    # =====================================================
    def _get_envelope_vec_for_sel_nodes (self, node_indices_1, solution_indices, dict_key, rsl_key, comp_index_4):
    # =====================================================
        # returns 2 matrices in which each row consists of [NX,QY,QZ,MX,MY,MZ] generalized force components
        out_max = np.zeros(len(node_indices_1) * 4).reshape(len(node_indices_1),4)
        out_min = np.zeros(len(node_indices_1) * 4).reshape(len(node_indices_1),4)
        aux     = np.zeros(len(node_indices_1) * 4).reshape(len(node_indices_1),4)

        for solution_index in solution_indices:
            vecx = self._get_nodes_vectorial_rsl(node_indices_1, solution_index, dict_key, rsl_key)
            aux[:,0:3] = vecx
            aux[:,3] = np.sqrt(vecx[:,0] * vecx[:,0] + vecx[:,1] * vecx[:,1] + vecx[:,2] * vecx[:,2])

            # keep the whole row (corresponding components) where the envelope component is exceeded
            mask = aux[:,comp_index_4] > out_max[:,comp_index_4]
            out_max[mask] = aux[mask]
            mask = aux[:,comp_index_4] < out_min[:,comp_index_4]
            out_min[mask] = aux[mask]

        return out_min, out_max

    # =====================================================
    def _get_envelope_scalar_for_sel_nodes (self, node_indices_1, solution_indices, dict_key, rsl_key):
    # =====================================================
        # returns 2 matrices in which each row consists of [NX,QY,QZ,MX,MY,MZ] generalized force components
        out_max = np.zeros(len(node_indices_1) )
        out_min = np.zeros(len(node_indices_1) )

        for solution_index in solution_indices:
            aux = self._get_nodes_scalar_rsl(node_indices_1, solution_index, dict_key, rsl_key)
            out_max = np.where(aux > out_max, aux, out_max)
            out_min = np.where(aux < out_min, aux, out_min)
        return out_min, out_max


    # =====================================================
    def  get_displ_envelope_for_sel_nodes (self,node_indices_1,solution_indices,comp):
    # =====================================================
        # returns 2 matrices of dimension (N_nodes,4) in which each row consists of [UX,UY,UZ,U-ABS]
        # it keeps envelope for the given displ. component and corresponding displacements
        # first matrix is for min envelope and second for the max one
        vecComps = ['X', 'Y', 'Z','ABS']
        comp_index = vecComps.index (comp)
        return self._get_envelope_vec_for_sel_nodes ( node_indices_1, solution_indices,\
                                                      'NODAL_RSL',"DISP_TRA",comp_index)

    # =====================================================
    def  get_reaction_forces_envelope_for_sel_nodes (self,node_indices_1,solution_indices,comp):
    # =====================================================
        # returns 2 matrices of dimension (N_nodes,4) in which each row consists of [UX,UY,UZ,U-ABS]
        # it keeps envelope for the given displ. component and corresponding displacements
        # first matrix is for min envelope and second for the max one
        vecComps = ['X', 'Y', 'Z','ABS']
        comp_index = vecComps.index (comp)
        return self._get_envelope_vec_for_sel_nodes ( node_indices_1, solution_indices,\
                                                      'NODAL_RES',"DISP_TRA",comp_index)

    # =====================================================
    def  get_rot_envelope_for_sel_nodes (self,node_indices_1,solution_indices,comp):
    # =====================================================
        # returns 2 matrices of dimension (N_nodes,4) in which each row consists of [RX,RY,RZ,R-ABS]
        # it keeps envelope for the given displ. component and corresponding displacements
        # first matrix is for min envelope and second for the max one
        vecComps = ['X', 'Y', 'Z','ABS']
        comp_index = vecComps.index (comp)
        return self._get_envelope_vec_for_sel_nodes ( node_indices_1, solution_indices,\
                                                      'NODAL_RSL',"DISP_ROT",comp_index)

    # =====================================================
    def  get_reaction_moments_envelope_for_sel_nodes (self,node_indices_1,solution_indices,comp):
    # =====================================================
        # returns 2 matrices of dimension (N_nodes,4) in which each row consists of [RX,RY,RZ,R-ABS]
        # it keeps envelope for the given displ. component and corresponding displacements
        # first matrix is for min envelope and second for the max one
        vecComps = ['X', 'Y', 'Z','ABS']
        comp_index = vecComps.index (comp)
        return self._get_envelope_vec_for_sel_nodes ( node_indices_1, solution_indices,\
                                                      'NODAL_RES',"DISP_ROT",comp_index)

    # =====================================================
    def  get_pressure_envelope_for_sel_nodes (self,node_indices_1,solution_indices):
    # =====================================================
        # returns 2 matrices of dimension (N_nodes)
        # first matrix is for min envelope and second for the max one
        return self._get_envelope_scalar_for_sel_nodes ( node_indices_1, solution_indices,\
                                                         'NODAL_RSL','PPRESS')

    # =====================================================
    def  get_pressure_head_envelope_for_sel_nodes (self,node_indices_1,solution_indices):
    # =====================================================
        # returns 2 matrices of dimension (N_nodes)
        # first matrix is for min envelope and second for the max one
        return self._get_envelope_scalar_for_sel_nodes ( node_indices_1, solution_indices,\
                                                         'NODAL_RSL','PRES_HEAD')

    # =====================================================
    def  get_temperature_envelope_for_sel_nodes (self,node_indices_1,solution_indices):
    # =====================================================
        # returns 2 matrices of dimension (N_nodes)
        # first matrix is for min envelope and second for the max one
        return self._get_envelope_scalar_for_sel_nodes ( node_indices_1, solution_indices,\
                                                         'NODAL_RSL','TEMP')

    # =====================================================
    def  get_humidity_envelope_for_sel_nodes (self,node_indices_1,solution_indices):
    # =====================================================
        # returns 2 matrices of dimension (N_nodes)
        # first matrix is for min envelope and second for the max one
        return self._get_envelope_scalar_for_sel_nodes ( node_indices_1, solution_indices,\
                                                         'NODAL_RSL','HUMID')



    # =====================================================
    def  get_veloc_envelope_for_sel_nodes (self,node_indices_1,solution_indices,comp):
    # =====================================================
        # returns 2 matrices of dimension (N_nodes,4) in which each row consists of [UX,UY,UZ,U-ABS]
        # it keeps envelope for the given displ. component and corresponding displacements
        # first matrix is for min envelope and second for the max one
        vecComps = ['X', 'Y', 'Z','ABS']
        comp_index = vecComps.index (comp)
        return self._get_envelope_vec_for_sel_nodes ( node_indices_1, solution_indices,\
                                                      'NODAL_VEL',"DISP_TRA",comp_index)

    # =====================================================
    def  get_accel_envelope_for_sel_nodes (self,node_indices_1,solution_indices,comp):
    # =====================================================
        # returns 2 matrices of dimension (N_nodes,4) in which each row consists of [UX,UY,UZ,U-ABS]
        # it keeps envelope for the given displ. component and corresponding displacements
        # first matrix is for min envelope and second for the max one
        vecComps = ['X', 'Y', 'Z','ABS']
        comp_index = vecComps.index (comp)
        return self._get_envelope_vec_for_sel_nodes ( node_indices_1, solution_indices,\
                                                      'NODAL_ACC',"DISP_TRA",comp_index)

