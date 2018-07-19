#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 10 11:10:39 2018

@author: tam10

    This file is part of oniom_ee_resp.

    oniom_ee_resp is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    oniom_ee_resp is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with oniom_ee_resp.  If not, see <http://www.gnu.org/licenses/>.
"""

class Parameters(object):
    """
    The Parameters object is a class to organise the different types of 
    parameters used in an Amber calculation.
    """
    def __init__(self, VdWs = None, stretches = None, bends = None, torsions = None, improper_torsions = None):
        self.VdWs = VdWs
        self.stretches = stretches
        self.bends = bends
        self.torsions = torsions
        self.improper_torsions = improper_torsions
        self.nonbon_string = ""
        
    def update_parameters(self, other_params, replace = True):
        """
        Update this Parameter object with a new one.
        With replace=True, parameters with the same type will be replaced
        With replace=False, parameters with the same type will be ignored
        """
        for other_VdW in other_params.VdWs:
            self.update_VdW(other_VdW, replace)
            
        for other_stretch in other_params.stretches:
            self.update_stretch(other_stretch, replace)
            
        for other_bend in other_params.bends:
            self.update_bend(other_bend, replace)
            
        for other_torsion in other_params.torsions:
            self.update_torsion(other_torsion, replace)
            
        for other_improper_torsion in other_params.improper_torsions:
            self.update_improper_torsion(other_improper_torsion, replace)
        
    def update_VdW(self, other_VdW, replace = True):
        """
        Update or add a VdW parameter.
        With replace=True, parameters with the same type will be replaced
        With replace=False, parameters with the same type will be ignored
        """
        for i, VdW in enumerate(self.VdWs):
            if VdW.type_equivalent(other_VdW):
                if replace:
                    self.VdWs[i] = other_VdW
                return
        self.VdWs.append(other_VdW)
        
    def update_stretch(self, other_stretch, replace = True):
        """
        Update or add a Stretch parameter.
        With replace=True, parameters with the same type will be replaced
        With replace=False, parameters with the same type will be ignored
        """
        for i, stretch in enumerate(self.stretches):
            if stretch.type_equivalent(other_stretch)[0]:
                if replace:
                    self.stretches[i] = other_stretch
                return
        self.stretches.append(other_stretch)
        
    def update_bend(self, other_bend, replace = True):
        """
        Update or add a Bend parameter.
        With replace=True, parameters with the same type will be replaced
        With replace=False, parameters with the same type will be ignored
        """
        for i, bend in enumerate(self.bends):
            if bend.type_equivalent(other_bend)[0]:
                if replace:
                    self.bends[i] = other_bend
                return
        self.bends.append(other_bend)
        
    def update_torsion(self, other_torsion, replace = True):
        """
        Update or add a Torsion parameter.
        With replace=True, parameters with the same type will be replaced
        With replace=False, parameters with the same type will be ignored
        """
        for i, torsion in enumerate(self.torsions):
            if torsion.type_equivalent(other_torsion)[0]:
                if replace:
                    self.torsions[i] = other_torsion
                return
        self.torsions.append(other_torsion)
        
    def update_improper_torsion(self, other_improper_torsion, replace = True):
        """
        Update or add an Improper Torsion parameter.
        With replace=True, parameters with the same type will be replaced
        With replace=False, parameters with the same type will be ignored
        """
        for i, improper_torsion in enumerate(self.improper_torsions):
            if improper_torsion.type_equivalent(other_improper_torsion)[0]:
                if replace:
                    self.improper_torsions[i] = other_improper_torsion
                return
        self.improper_torsions.append(other_improper_torsion)
        
    @property
    def VdWs(self):
        return self._VdWs
    @VdWs.setter
    def VdWs(self, value):
        if value is None:
            self._VdWs = []
        elif hasattr(value, "__len__") and all([isinstance(VdWs, VdWParam) for VdWs in value]):
            self._VdWs = list(value)
        else:
            raise ValueError("'VdWs' must be a list of VdWParam objects")
        
    @property
    def stretches(self):
        return self._stretches
    @stretches.setter
    def stretches(self, value):
        if value is None:
            self._stretches = []
        elif hasattr(value, "__len__") and all([isinstance(stretch, StretchParam) for stretch in value]):
            self._stretches = list(value)
        else:
            raise ValueError("'stretches' must be a list of StretchParam objects")
        
    @property
    def bends(self):
        return self._bends
    @bends.setter
    def bends(self, value):
        if value is None:
            self._bends = []
        elif hasattr(value, "__len__") and all([isinstance(bend, BendParam) for bend in value]):
            self._bends = list(value)
        else:
            raise ValueError("'bends' must be a list of BendParam objects")
        
    @property
    def torsions(self):
        return self._torsions
    @torsions.setter
    def torsions(self, value):
        if value is None:
            self._torsions = []
        elif hasattr(value, "__len__") and all([isinstance(torsion, TorsionParam) for torsion in value]):
            self._torsions = list(value)
        else:
            raise ValueError("'torsions' must be a list of TorsionParam objects")
        
    @property
    def improper_torsions(self):
        return self._improper_torsions
    @improper_torsions.setter
    def improper_torsions(self, value):
        if value is None:
            self._improper_torsions = []
        elif hasattr(value, "__len__") and all([isinstance(improper_torsion, ImproperTorsionParam) for improper_torsion in value]):
            self._improper_torsions = list(value)
        else:
            raise ValueError("'improper_torsions' must be a list of ImproperTorsionParam objects")

    @property
    def nonbon_string(self):
        return self._nonbon_string
    @nonbon_string.setter
    def nonbon_string(self, value):
        if value is None:
            self._nonbon_string = ""
        elif isinstance(value, str):
            self._nonbon_string = value
        else:
            raise ValueError("'nonbon_string' must be a string, not {}".format(type(value)))
            
    def get_string(self):
        """
        Write all parameters in a Gaussian-readable string.
        """
        string = self.nonbon_string + "\n"
        for param in self.VdWs + self.stretches + self.bends + self.torsions + self.improper_torsions:
            string += param.get_string() + "\n"
            
        return string
            
    def __repr__(self):
        return "Parameters({:d} VdWs, {:d} Stretches, {:d} Bends, {:d} Torsions, {:d} ImproperTorsions".format(
            len(self.VdWs),
            len(self.stretches),
            len(self.bends),
            len(self.torsions),
            len(self.improper_torsions),
        )
    
class VdWParam(object):
    """
    Van der Waal Parameter for Amber calcutions
    type0: Amber type of Atom
    radius: Atomic radius
    potential: Potential well depth
    """
    def __init__(self, type0, radius, potential):
        self.type0 = type0
        self.radius = radius
        self.potential = potential
        self._header = "HrmStr1"
        self._type0_f = "{:<3s}"
        self._radius_f = "{:7.4f}"
        self._v_f = "{:7.4f}"
    
    def get_string(self):
        """
        Write this parameter in a Gaussian-readable string.
        """
        return " ".join([
            self._header,
            self._type0_f.format(self.type0),
            self._radius_f.format(self.radius),
            self._potential_f.format(self.potential),
        ])
        
    def type_equivalent(self, other):
        """
        Check whether this parameter has the same types as another
        returns:
            (equivalent<bool>)
        """
        if not isinstance(other, VdWParam):
            raise ValueError("'other' must be a VdWParam object, not {}".format(type(other)))
        
        return self.type0 == other.type0
        
    def value_equivalent(self, other):
        """
        Check whether this parameter is equivalent to another
        returns:
            (equivalent<bool>)
        """
        type_equivalent = self.type_equivalent(other)
        
        if not type_equivalent:
            return False
            
        return (self.radius == other.radius) and (self.v == other.v)
        
    @property
    def type0(self):
        return self._type0
    @type0.setter
    def type0(self, value):
        self._type0 = str(value)
        
    @property
    def radius(self):
        return self._radius
    @radius.setter
    def radius(self, value):
        self._radius = float(value)
        
    @property
    def potential(self):
        return self._potential
    @potential.setter
    def potential(self, value):
        self._potential = float(value)
        
    def __repr__(self):
        return "VdWParam({})".format(
            ", ".join(
                [s + "=" + getattr(self, "_{}_f".format(s)).format(getattr(self, s)) for s in ["type0", "radius", "potential"]]
            )
        )
    
class StretchParam(object):
    """
    Stretching Parameter for Amber calcutions
    type0: Amber type of Atom 0
    type1: Amber type of Atom 1
    k: Force constant
    r: Equilibrium distance
    """
    def __init__(self, type0, type1, k, r):
        self.type0 = type0
        self.type1 = type1
        self.k = k
        self.r = r
        self._header = "HrmStr1"
        self._type0_f = "{:<3s}"
        self._type1_f = "{:<3s}"
        self._k_f = "{:7.4f}"
        self._r_f = "{:7.4f}"
    
    def get_string(self):
        """
        Write this parameter in a Gaussian-readable string.
        """
        return " ".join([
            self._header,
            self._type0_f.format(self.type0),
            self._type1_f.format(self.type1),
            self._k_f.format(self.k),
            self._r_f.format(self.r),
        ])
        
    def type_equivalent(self, other):
        """
        Check whether this parameter has the same types as another
        returns:
            (equivalent<bool>, isreversed<bool>)
        """
        if not isinstance(other, StretchParam):
            raise ValueError("'other' must be a StretchParam object, not {}".format(type(other)))
            
        at0 = self.type0
        at1 = self.type1
        bt0 = other.type0
        bt1 = other.type1
        
        if (at0 == bt0) and (at1 == bt1):
            return (True, False)
        
        if (at0 == bt1) and (at1 == bt0):
            return (True, True)
            
        return (False, False)
        
    def value_equivalent(self, other):
        """
        Check whether this parameter is equivalent to another
        returns:
            (equivalent<bool>, isreversed<bool>)
        """
        type_equivalent, isreversed = self.type_equivalent(other)
        
        if not type_equivalent:
            return (False, isreversed)
            
        ak = self.k
        ar = self.r
        bk = other.k
        br = other.r
            
        if (ak == bk) and (ar == br):
            return (True, isreversed)
        
        return (False, isreversed)
        
    @property
    def type0(self):
        return self._type0
    @type0.setter
    def type0(self, value):
        self._type0 = str(value)
        
    @property
    def type1(self):
        return self._type1
    @type1.setter
    def type1(self, value):
        self._type1 = str(value)
        
    @property
    def k(self):
        return self._k
    @k.setter
    def k(self, value):
        self._k = float(value)
        
    @property
    def r(self):
        return self._r
    @r.setter
    def r(self, value):
        self._r = float(value)
        
    def __repr__(self):
        return "StretchParam({})".format(
            ", ".join(
                [s + "=" + getattr(self, "_{}_f".format(s)).format(getattr(self, s)) for s in ["type0", "type1", "k", "r"]]
            )
        )
    
class BendParam(object):
    """
    Bending Parameter for Amber calcutions
    type0: Amber type of Atom 0
    type1: Amber type of Atom 1
    type2: Amber type of Atom 2
    k: Force constant
    angle: Equilibrium angle
    """
    def __init__(self, type0, type1, type2, k, angle):
        self.type0 = type0
        self.type1 = type1
        self.type2 = type2
        self.k = k
        self.angle = angle
        self._header = "HrmBnd1"
        self._type0_f = "{:<3s}"
        self._type1_f = "{:<3s}"
        self._type2_f = "{:<3s}"
        self._k_f = "{:7.4f}"
        self._angle_f = "{:7.4f}"
    
    def get_string(self):
        """
        Write this parameter in a Gaussian-readable string.
        """
        return " ".join([
            self._header,
            self._type0_f.format(self.type0),
            self._type1_f.format(self.type1),
            self._type2_f.format(self.type2),
            self._k_f.format(self.k),
            self._angle_f.format(self.angle),
        ])
        
    def type_equivalent(self, other):
        """
        Check whether this parameter has the same types as another
        returns:
            (equivalent<bool>, isreversed<bool>)
        """
        if not isinstance(other, BendParam):
            raise ValueError("'other' must be a BendParam object, not {}".format(type(other)))
            
        at0 = self.type0
        at1 = self.type1
        at2 = self.type2
        bt0 = other.type0
        bt1 = other.type1
        bt2 = other.type2
        
        if (at0 == bt0) and (at1 == bt1) and (at2 == bt2):
            return (True, False)
        
        if (at0 == bt2) and (at1 == bt1) and (at2 == bt0):
            return (True, True)
            
        return (False, False)
        
    def value_equivalent(self, other):
        """
        Check whether this parameter is equivalent to another
        returns:
            (equivalent<bool>, isreversed<bool>)
        """
        type_equivalent, isreversed = self.type_equivalent(other)
        
        if not type_equivalent:
            return (False, isreversed)
            
        ak = self.k
        ar = self.angle
        bk = other.k
        br = other.angle
            
        if (ak == bk) and (ar == br):
            return (True, isreversed)
        
        return (False, isreversed)
        
    @property
    def type0(self):
        return self._type0
    @type0.setter
    def type0(self, value):
        self._type0 = str(value)
        
    @property
    def type1(self):
        return self._type1
    @type1.setter
    def type1(self, value):
        self._type1 = str(value)
        
    @property
    def type2(self):
        return self._type2
    @type2.setter
    def type2(self, value):
        self._type2 = str(value)
        
    @property
    def k(self):
        return self._k
    @k.setter
    def k(self, value):
        self._k = float(value)
        
    @property
    def angle(self):
        return self._angle
    @angle.setter
    def angle(self, value):
        self._angle = float(value)
        
    def __repr__(self):
        return "BendParam({})".format(
            ", ".join(
                [s + "=" + getattr(self, "_{}_f".format(s)).format(getattr(self, s)) for s in ["type0", "type1", "type2", "k", "angle"]]
            )
        )
    
class TorsionParam(object):
    """
    Torsion Parameter for Amber calcutions
    type[0-4]: Amber type of Atom 0-4
    v[0-4]: barrier for cosine multiplicity 1-5
    gamma[0-4]: phase offset for cosine periodicity 1-5
    n_paths: Number of paths
    """
    def __init__(
            self, 
            type0, 
            type1, 
            type2, 
            type3, 
            v0 = 0., 
            v1 = 0., 
            v2 = 0., 
            v3 = 0., 
            gamma0 = 0., 
            gamma1 = 0., 
            gamma2 = 0., 
            gamma3 = 0., 
            n_paths = 0
        ):
        
        self.type0 = type0
        self.type1 = type1
        self.type2 = type2
        self.type3 = type3
        self.v0 = v0
        self.v1 = v1
        self.v2 = v2
        self.v3 = v3
        self.gamma0 = gamma0
        self.gamma1 = gamma1
        self.gamma2 = gamma2
        self.gamma3 = gamma3
        self.n_paths = n_paths
        
        self._header = "AmbTrs"
        self._type0_f = "{:<3s}"
        self._type1_f = "{:<3s}"
        self._type2_f = "{:<3s}"
        self._type3_f = "{:<3s}"
        self._v0_f = "{:7.4f}"
        self._v1_f = "{:7.4f}"
        self._v2_f = "{:7.4f}"
        self._v3_f = "{:7.4f}"
        self._gamma0_f = "{:7.4f}"
        self._gamma1_f = "{:7.4f}"
        self._gamma2_f = "{:7.4f}"
        self._gamma3_f = "{:7.4f}"
        self._n_paths_f = "{:3.1f}"
    
    def get_string(self):
        """
        Write this parameter in a Gaussian-readable string.
        """
        return " ".join([
            self._header,
            self._type0_f.format(self.type0 if self.type0 != 'X' else '*'),
            self._type1_f.format(self.type1 if self.type0 != 'X' else '*'),
            self._type2_f.format(self.type2 if self.type0 != 'X' else '*'),
            self._type3_f.format(self.type3 if self.type0 != 'X' else '*'),
            self._v0_f.format(self.v0),
            self._v1_f.format(self.v1),
            self._v2_f.format(self.v2),
            self._v3_f.format(self.v3),
            self._gamma0_f.format(self.gamma0),
            self._gamma1_f.format(self.gamma1),
            self._gamma2_f.format(self.gamma2),
            self._gamma3_f.format(self.gamma3),
            self._n_paths_f.format(self.n_paths)
        ])
        
    def type_equivalent(self, other):
        """
        Check whether this parameter has the same types as another
        returns:
            (equivalent<bool>, isreversed<bool>)
        """
        if not isinstance(other, TorsionParam):
            raise ValueError("'other' must be a TorsionParam object, not {}".format(type(other)))
            
        at0 = self.type0
        at1 = self.type1
        at2 = self.type2
        at3 = self.type3
        bt0 = other.type0
        bt1 = other.type1
        bt2 = other.type2
        bt3 = other.type3
        
        if (at0 == bt0) and (at1 == bt1) and (at2 == bt2) and (at3 == bt3):
            return (True, False)
        
        if (at0 == bt3) and (at1 == bt2) and (at2 == bt1) and (at3 == bt0):
            return (True, True)
            
        return (False, False)
        
    def value_equivalent(self, other):
        """
        Check whether this parameter is equivalent to another
        returns:
            (equivalent<bool>, isreversed<bool>)
        """
        type_equivalent, isreversed = self.type_equivalent(other)
        
        if not type_equivalent:
            return (False, isreversed)
        
        a_vs = [getattr(self, "v{}".format(str(i))) for i in range(4)]
        b_vs = [getattr(other, "v{}".format(str(i))) for i in range(4)]
        
        a_gammas = [getattr(self, "gamma{}".format(str(i))) for i in range(4)]
        b_gammas = [getattr(other, "gamma{}".format(str(i))) for i in range(4)]
        
        
        if isreversed:
            if a_vs == list(reversed(b_vs)) and a_gammas == list(reversed(b_gammas)):
                return (True, True)
        else:
            if a_vs == b_vs and a_gammas == b_gammas:
                return (True, False)
        
        return (False, isreversed)
            
    @property
    def type0(self):
        return self._type0
    @type0.setter
    def type0(self, value):
        self._type0 = str(value)
        
    @property
    def type1(self):
        return self._type1
    @type1.setter
    def type1(self, value):
        self._type1 = str(value)
        
    @property
    def type2(self):
        return self._type2
    @type2.setter
    def type2(self, value):
        self._type2 = str(value)
        
    @property
    def type3(self):
        return self._type3
    @type3.setter
    def type3(self, value):
        self._type3 = str(value)
        
    @property
    def v0(self):
        return self._v0
    @v0.setter
    def v0(self, value):
        self._v0 = float(value)
        
    @property
    def v1(self):
        return self._v1
    @v1.setter
    def v1(self, value):
        self._v1 = float(value)
        
    @property
    def v2(self):
        return self._v2
    @v2.setter
    def v2(self, value):
        self._v2 = float(value)
        
    @property
    def v3(self):
        return self._v3
    @v3.setter
    def v3(self, value):
        self._v3 = float(value)
        
    @property
    def gamma0(self):
        return self._gamma0
    @gamma0.setter
    def gamma0(self, value):
        self._gamma0 = float(value)
        
    @property
    def gamma1(self):
        return self._gamma1
    @gamma1.setter
    def gamma1(self, value):
        self._gamma1 = float(value)
        
    @property
    def gamma2(self):
        return self._gamma2
    @gamma2.setter
    def gamma2(self, value):
        self._gamma2 = float(value)
        
    @property
    def gamma3(self):
        return self._gamma3
    @gamma3.setter
    def gamma3(self, value):
        self._gamma3 = float(value)
        
    @property
    def n_paths(self):
        return self._n_paths
    @n_paths.setter
    def n_paths(self, value):
        self._n_paths = float(value)
        
    def __repr__(self):
        return "TorsionParam({})".format(
            ", ".join(
                [
                    s + "=" + getattr(self, "_{}_f".format(s)).format(getattr(self, s)) for s in [
                        "type0", 
                        "type1", 
                        "type2", 
                        "type3", 
                        "v0",
                        "v1",
                        "v2",
                        "v3",
                        "gamma0",
                        "gamma1",
                        "gamma2",
                        "gamma3",
                        "n_paths"
                    ]
                ]
            )
        )
    
class ImproperTorsionParam(object):
    """
    Improper Torsion Parameter for Amber calcutions
    type[0-4]: Amber type of Atom 0-4
    v: barrier 
    gamma: phase offset
    nt: Cosine periodicity
    """
    def __init__(
            self, 
            type0, 
            type1, 
            type2, 
            type3, 
            v = 0.,  
            gamma = 0., 
            nt = 0
        ):
        
        self.type0 = type0
        self.type1 = type1
        self.type2 = type2
        self.type3 = type3
        self.v = v
        self.gamma = gamma
        self.nt = nt
        
        self._header = "ImpTrs"
        self._type0_f = "{:<3s}"
        self._type1_f = "{:<3s}"
        self._type2_f = "{:<3s}"
        self._type3_f = "{:<3s}"
        self._v_f = "{:7.4f}"
        self._gamma_f = "{:7.4f}"
        self._nt_f = "{:3.1f}"
    
    def get_string(self):
        """
        Write this parameter in a Gaussian-readable string.
        """
        return " ".join([
            self._header,
            self._type0_f.format(self.type0 if self.type0 != 'X' else '*'),
            self._type1_f.format(self.type1 if self.type0 != 'X' else '*'),
            self._type2_f.format(self.type2 if self.type0 != 'X' else '*'),
            self._type3_f.format(self.type3 if self.type0 != 'X' else '*'),
            self._v_f.format(self.v),
            self._gamma_f.format(self.gamma),
            self._nt_f.format(self.nt)
        ])
        
    def type_equivalent(self, other):
        """
        Check whether this parameter has the same types as another
        returns:
            (equivalent<bool>, isreversed<bool>)
        """
        if not isinstance(other, ImproperTorsionParam):
            raise ValueError("'other' must be a ImproperTorsionParam object, not {}".format(type(other)))
            
        at0 = self.type0
        at1 = self.type1
        at2 = self.type2
        at3 = self.type3
        bt0 = other.type0
        bt1 = other.type1
        bt2 = other.type2
        bt3 = other.type3
        
        if (at0 == bt0) and (at1 == bt1) and (at2 == bt2) and (at3 == bt3):
            return (True, False)
        
        if (at0 == bt3) and (at1 == bt2) and (at2 == bt1) and (at3 == bt0):
            return (True, True)
            
        return (False, False)
        
    def value_equivalent(self, other):
        """
        Check whether this parameter is equivalent to another
        returns:
            (equivalent<bool>, isreversed<bool>)
        """
        type_equivalent, isreversed = self.type_equivalent(other)
        
        if not type_equivalent:
            return (False, isreversed)
            
        if (self.v == other.v) and (self.gamma == other.gamma) and (self.nt == other.nt):
            return (True, isreversed)
        else:
            return (False, isreversed)
        
        return (False, isreversed)
        
    @property
    def type0(self):
        return self._type0
    @type0.setter
    def type0(self, value):
        self._type0 = str(value)
        
    @property
    def type1(self):
        return self._type1
    @type1.setter
    def type1(self, value):
        self._type1 = str(value)
        
    @property
    def type2(self):
        return self._type2
    @type2.setter
    def type2(self, value):
        self._type2 = str(value)
        
    @property
    def type3(self):
        return self._type3
    @type3.setter
    def type3(self, value):
        self._type3 = str(value)
        
    @property
    def v(self):
        return self._v
    @v.setter
    def v(self, value):
        self._v = float(value)
        
    @property
    def gamma(self):
        return self._gamma
    @gamma.setter
    def gamma(self, value):
        self._gamma = float(value)
        
    @property
    def nt(self):
        return self._nt
    @nt.setter
    def nt(self, value):
        self._nt = float(value)
        
    def __repr__(self):
        return "ImproperTorsionParam({})".format(
            ", ".join(
                [
                    s + "=" + getattr(self, "_{}_f".format(s)).format(getattr(self, s)) for s in [
                        "type0", 
                        "type1", 
                        "type2", 
                        "type3", 
                        "v",
                        "gamma",
                        "nt"
                    ]
                ]
            )
        )
    