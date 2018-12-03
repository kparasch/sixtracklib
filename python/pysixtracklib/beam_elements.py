#!/usr/bin/env python
# -*- coding: utf-8 -*-

from cobjects import CBuffer, CObject, CField
import numpy as np


class Drift(CObject):
    _typeid = 2
    length = CField(0, 'real', default=0.0, alignment=8)


class DriftExact(CObject):
    _typeid = 3
    length = CField(0, 'real', default=0.0, alignment=8)


class Multipole(CObject):
    _typeid = 4
    order = CField(0, 'int64',   default=0,    alignment=8)
    length = CField(1, 'real',    default=0.0,  alignment=8)
    hxl = CField(2, 'real',    default=0.0,  alignment=8)
    hyl = CField(3, 'real',    default=0.0,  alignment=8)
    bal = CField(4, 'real',    default=0.0,
                 length='2 * order + 2', pointer=True, alignment=8)

    def _factorial(self, x):
        if not isinstance(x, int):
            return 0
        return (x > 0) and (x * self._factorial(x - 1)) or 1

    def __init__(self, order=None, knl=None, ksl=None, bal=None, **kwargs):

        if bal is None and (not(knl is None) or not(ksl is None)):
            if knl is None:
                knl = []
            if ksl is None:
                ksl = []
            if order is None:
                order = 0

            n = max((order + 1), max(len(knl), len(ksl)))
            _knl = np.array(knl)
            nknl = np.zeros(n, dtype=_knl.dtype)
            nknl[:len(knl)] = knl
            knl = nknl
            del(_knl)

            _ksl = np.array(ksl)
            nksl = np.zeros(n, dtype=_ksl.dtype)
            nksl[:len(ksl)] = ksl
            ksl = nksl
            del(_ksl)

            assert(n > 0)
            order = n - 1

            bal = np.zeros(2 * order + 2)
            assert(len(knl) == len(ksl))

            for ii in range(0, len(knl)):
                inv_factorial = 1.0 / float(self._factorial(ii))
                jj = 2 * ii
                bal[jj] = knl[ii] * inv_factorial
                bal[jj + 1] = ksl[ii] * inv_factorial

        elif not(bal is None) and bal and \
                len(bal) > 2 and ((len(bal) % 2) == 0):

            order = (len(bal) - 2) / 2
            assert(order > 0)

        elif bal is None and knl is None and ksl is None and \
                not(order is None) and order > 0:
            bal = np.zeros(2 * order + 2)

        if not(bal is None or order is None):
            CObject.__init__(self, bal=bal, order=order, **kwargs)
        else:
            CObject.__init__(self, bal=[], order=0, **kwargs)


class Cavity(CObject):
    _typeid = 5
    voltage = CField(0, 'real', default=0.0,  alignment=8)
    frequency = CField(1, 'real', default=0.0,  alignment=8)
    lag = CField(2, 'real', default=0.0,  alignment=8)


class XYShift(CObject):
    _typeid = 6
    dx = CField(0, 'real',   default=0.0,  alignment=8)
    dy = CField(1, 'real',   default=0.0,  alignment=8)


class SRotation(CObject):
    _typeid = 7
    cos_z = CField(0, 'real',   default=1.0,  alignment=8)
    sin_z = CField(1, 'real',   default=0.0,  alignment=8)

    def __init__(self, angle=0, **nargs):
        anglerad = angle/180*np.pi
        cos_z = np.cos(anglerad)
        sin_z = np.sin(anglerad)
        CObject.__init__(self,
                         cos_z=cos_z, sin_z=sin_z, **nargs)


class BeamBeam4D(CObject):
    _typeid = 8
    size = CField(0, 'uint64', const=True, default=0)
    data = CField(1, 'float64',   default=0.0,
                  length='size', pointer=True)

    def __init__(self, data=None, **kwargs):
        if data is None:
            slots = ('q_part', 'N_part', 'sigma_x', 'sigma_y', 'beta_s',
                     'min_sigma_diff', 'Delta_x', 'Delta_y', 'Dpx_sub', 'Dpy_sub', 'enabled')
            data = [kwargs[ss] for ss in slots]
            CObject.__init__(self, size=len(data), data=data, **kwargs)
        else:
            CObject.__init__(self, **kwargs)


##### BB6D interanl objects ######
class ParBoost(CObject):
    _typeid = 903
    
    sphi = CField(0, 'real', default=0.0, alignment=8)
    cphi = CField(1, 'real', default=0.0, alignment=8)
    tphi = CField(2, 'real', default=0.0, alignment=8)
    salpha = CField(3, 'real', default=0.0, alignment=8)
    calpha  = CField(4, 'real',  default=0.0, alignment=8)

class Sigmas(CObject):
    _typeid = 902

    Sig_11_0 = CField(0, 'real', default=0.0, alignment=8)
    Sig_12_0 = CField(1, 'real', default=0.0, alignment=8)
    Sig_13_0 = CField(2, 'real', default=0.0, alignment=8)
    Sig_14_0 = CField(3, 'real', default=0.0, alignment=8)
    Sig_22_0 = CField(4, 'real', default=0.0, alignment=8)
    Sig_23_0 = CField(5, 'real', default=0.0, alignment=8)
    Sig_24_0 = CField(6, 'real', default=0.0, alignment=8)
    Sig_33_0 = CField(7, 'real', default=0.0, alignment=8)
    Sig_34_0 = CField(8, 'real', default=0.0, alignment=8)
    Sig_44_0 = CField(9, 'real', default=0.0, alignment=8)

##### End BB6D internal obj ######

class BeamBeam6D(CObject):
    _typeid = 9

    q_part = CField(0, 'real', default=0.0, alignment=8)
    parboost = CField(1, ParBoost)
    Sigmas_0_star = CField(2, Sigmas)
    min_sigma_diff = CField(3, 'real', default=0.0, alignment=8)
    threshold_singular = CField(4, 'real', default=0.0, alignment=8)
    N_slices = CField(5, 'int64', const=True, alignment=8)

    delta_x = CField(6, 'real', default=0.0, alignment=8)
    delta_y = CField(7, 'real', default=0.0, alignment=8)
    x_CO = CField(8, 'real', default=0.0, alignment=8)
    px_CO = CField(9, 'real', default=0.0, alignment=8)
    y_CO = CField(10, 'real', default=0.0, alignment=8)
    py_CO = CField(11, 'real', default=0.0, alignment=8)
    sigma_CO = CField(12, 'real', default=0.0, alignment=8)
    delta_CO = CField(13, 'real', default=0.0, alignment=8)
    Dx_sub = CField(14, 'real', default=0.0, alignment=8)
    Dpx_sub = CField(15, 'real', default=0.0, alignment=8)
    Dy_sub = CField(16, 'real', default=0.0, alignment=8)
    Dpy_sub = CField(17, 'real', default=0.0, alignment=8)
    Dsigma_sub = CField(18, 'real', default=0.0, alignment=8)
    Ddelta_sub = CField(19, 'real', default=0.0, alignment=8)
    enabled = CField(20, 'int64', default=0, alignment=8)

    N_part_per_slice = CField(21, 'real', default=0.0,
                 length='N_slices', pointer=True, alignment=8)
    x_slices_star = CField(22, 'real', default=0.0,
                 length='N_slices', pointer=True, alignment=8)
    y_slices_star = CField(23, 'real', default=0.0,
                 length='N_slices', pointer=True, alignment=8)
    sigma_slices_star = CField(24, 'real', default=0.0,
                 length='N_slices', pointer=True, alignment=8)


    def __init__(self, data=None, **kwargs):

        import pysixtrack
        data = pysixtrack.BB6Ddata.BB6D_init(
            **{kk: kwargs[kk] for kk in kwargs.keys() if kk != 'cbuffer'})

        CObject.__init__(self, N_slices=data.N_slices)
        
        self.q_part = data.q_part

        for attr in "sphi cphi tphi salpha calpha".split():
            setattr(self.parboost, attr, getattr(data.parboost, attr))

        for attr in "Sig_11_0 Sig_12_0 Sig_13_0 Sig_14_0 Sig_22_0 Sig_23_0 Sig_24_0 Sig_33_0 Sig_34_0 Sig_44_0".split(): 
            setattr(self.Sigmas_0_star, attr, getattr(data.Sigmas_0_star, attr))

        self.min_sigma_diff = data.min_sigma_diff
        self.threshold_singular = data.threshold_singular
        self.N_slices = data.N_slices

        self.delta_x = data.delta_x
        self.delta_y = data.delta_y
        self.x_CO = data.x_CO
        self.px_CO = data.px_CO
        self.y_CO = data.y_CO
        self.py_CO = data.py_CO
        self.sigma_CO = data.sigma_CO
        self.delta_CO = data.delta_CO
        self.Dx_sub = data.Dx_sub
        self.Dpx_sub = data.Dpx_sub
        self.Dy_sub = data.Dy_sub
        self.Dpy_sub = data.Dpy_sub
        self.Dsigma_sub = data.Dsigma_sub
        self.Ddelta_sub = data.Ddelta_sub
        self.enabled = data.enabled

        self.N_part_per_slice[:self.N_slices] = data.N_part_per_slice[:self.N_slices]
        self.x_slices_star[:self.N_slices] = data.x_slices_star[:self.N_slices]
        self.y_slices_star[:self.N_slices] = data.y_slices_star[:self.N_slices]
        self.sigma_slices_star[:self.N_slices] = data.sigma_slices_star[:self.N_slices]




class Elements(object):
    element_types = {'Cavity': Cavity,
                     'Drift': Drift,
                     'DriftExact': DriftExact,
                     'Multipole': Multipole,
                     #                     'RFMultipole': RFMultipole,
                     'SRotation': SRotation,
                     'XYShift': XYShift,
                     'BeamBeam6D': BeamBeam6D,
                     'BeamBeam4D': BeamBeam4D,
                     #                     'Line': Line,
                     #                     'Monitor': Monitor,
                     }

    def _mk_fun(self, buff, cls):
        def fun(*args, **nargs):
            # print(cls.__name__,nargs)
            return cls(cbuffer=buff, **nargs)
        return fun

    @classmethod
    def fromfile(cls, filename):
        cbuffer = CBuffer.fromfile(filename)
        return cls(cbuffer=cbuffer)

    @classmethod
    def fromline(cls, line):
        self = cls()
        for label, element_name, element in line:
            getattr(self, element_name)(**element._asdict())
        return self

    def tofile(self, filename):
        self.cbuffer.tofile(filename)

    def __init__(self, cbuffer=None):
        if cbuffer is None:
            self.cbuffer = CBuffer()
        else:
            self.cbuffer = cbuffer
        for name, cls in self.element_types.items():
            setattr(self, name, self._mk_fun(self.cbuffer, cls))
            self.cbuffer.typeids[cls._typeid] = cls

    def gen_builder(self):
        out = {}
        for name, cls in self.element_types.items():
            out[name] = getattr(self, name)
        return out

    def get_elements(self):
        n = self.cbuffer.n_objects
        return [self.cbuffer.get_object(i) for i in range(n)]

    def get(self, objid):
        return self.cbuffer.get_object(objid)
