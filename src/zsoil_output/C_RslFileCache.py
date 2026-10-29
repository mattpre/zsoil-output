import os
import numpy


#=====================================================
class RslFileCache ():
#=====================================================
    # Memory-maps binary result files (float32 records of fixed size) once and
    # keeps them mapped, so many small reads do not each pay for open/seek/read.
    # A file is re-mapped automatically when its size or modification time changes.

    #====================================================
    def __init__ (self):
    #=====================================================
        self.maps = {}  # path -> (stat_key, memmap or None)

    #====================================================
    def _get_map (self,path):
    #=====================================================
        st = os.stat (path)
        stat_key = (st.st_size,st.st_mtime_ns)
        entry = self.maps.get (path)
        if entry is not None and entry[0] == stat_key:
            return entry[1]
        if st.st_size < 4:
            mm = None
        else:
            mm = numpy.memmap (path,dtype=numpy.float32,mode="r")
        self.maps [path] = (stat_key,mm)
        return mm

    #====================================================
    def read_records (self,path,rec_size_int,rec_indices):
    #=====================================================
        # returns float32 array of shape (len(rec_indices),rec_size_int)
        # rec_indices are 0-based record numbers (record = rec_size_int float32 values)
        rec_indices = numpy.asarray (rec_indices,dtype=numpy.int64)
        mm = self._get_map (path)
        n_rec = 0 if mm is None or rec_size_int <= 0 else len(mm) // rec_size_int
        if rec_indices.size > 0 and (rec_indices.min() < 0 or rec_indices.max() >= n_rec):
            raise IndexError ("record index out of range in %s (%d complete records of %d floats)"
                              % (path,n_rec,rec_size_int))
        if rec_indices.size == 0:
            return numpy.zeros ((0,rec_size_int),dtype=numpy.float32)
        view = mm [:n_rec*rec_size_int].reshape (n_rec,rec_size_int)
        return numpy.array (view [rec_indices])

    #====================================================
    def close (self):
    #=====================================================
        # releases all mapped files (needed on Windows before the files can be overwritten);
        # read_records only returns copies, so dropping the references unmaps the files
        self.maps = {}
