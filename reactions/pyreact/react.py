import os
import ctypes as ct
import numpy as np

def array_ctype(ndim, dtype=np.float64, flags="C"):
    return [ct.POINTER(ct.c_int)]*ndim + [np.ctypeslib.ndpointer(ndim=ndim, dtype=dtype, flags=flags)]

def array_arg(a):
    arr = a
    return (*(ct.c_int(s) for s in arr.shape), arr)

class ReactionModule:
    libname = "../lib/libreact.so"
    module_name = "reaction_module"

    def __init__(self):
        self.load_lib()

    def load_lib(self, path=None):
        if path is None:
            path = os.path.dirname(__file__)
        libpath = os.path.abspath(os.path.join(path, self.libname))
        self.lib = ct.CDLL(libpath)

    def get_function(self, name, c_bind=True):
        if c_bind:
            return getattr(self.lib, name)
        else:
            return getattr(self.lib, f"__{self.module_name}_MOD_{name}")

    def test_func(self, a):
        f = self.get_function("test_func")
        f.restype = np.int
        f.argtypes = [*array_ctype(ndim=2, dtype=np.float64)]

        r = f(*array_arg(a))
        return r

    def compute_reaction(self, h, n_s, omega_m, omega_b, sigma_8, mg1, mass_loop,
                               z, k, Pk, is_transfer=False,
                               verbose=True):
        f = self.get_function("compute_reaction")
        f.restype = np.int
        f.argtypes = [*array_ctype(ndim=1, dtype=np.float64), # P(k, z=0)
                      *array_ctype(ndim=1, dtype=np.float64), # k
                      *array_ctype(ndim=1, dtype=np.float64), # z
                      ct.POINTER(ct.c_bool),       # is_transfer
                      ct.POINTER(ct.c_double),     # h
                      ct.POINTER(ct.c_double),     # n_s
                      ct.POINTER(ct.c_double),     # omega_m
                      ct.POINTER(ct.c_double),     # omega_b
                      ct.POINTER(ct.c_double),     # sigma_8
                      ct.POINTER(ct.c_double),     # modified gravity param
                      ct.POINTER(ct.c_int),        # mass_loop
                      *array_ctype(ndim=2, dtype=np.float64), # reaction (output)
                      *array_ctype(ndim=2, dtype=np.float64), # linear MG power spectrum (output)
                     ]
        reaction = np.zeros((len(k), len(z)), dtype=Pk.dtype, order="C")
        p_lin = np.zeros((len(k), len(z)), dtype=Pk.dtype, order="C")

        r =   f(*array_arg(np.ascontiguousarray(Pk, dtype=np.float64)),
                *array_arg(np.ascontiguousarray(k, dtype=np.float64)),
                *array_arg(np.ascontiguousarray(z[::-1], dtype=np.float64)),
                ct.c_bool(is_transfer),
                ct.c_double(h), ct.c_double(n_s), ct.c_double(omega_m), ct.c_double(omega_b), ct.c_double(sigma_8),
                ct.c_double(mg1),
                ct.c_int(mass_loop),
                *array_arg(reaction),
                *array_arg(p_lin),
                )
        if r != 0:
            raise RuntimeError("Reaction code terminated with an error.")
        return reaction[:,::-1].T, p_lin[:,::-1].T
