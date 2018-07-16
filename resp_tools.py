#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 10 11:04:57 2018

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
    along with LepsPy.  If not, see <http://www.gnu.org/licenses/>.
"""

import numpy as np
from numba import jit
import sys

class RESP_Optimiser(object):
    def __init__(self):
        
        self._atomic_positions = None
        self._atom_potentials = None
        self._esp_positions = None
        self._qm_esps = None
        self._mm_esps = None
        self._esp_derived_charges = None
        self._model_esp_derived_charges = None
        self._resp_derived_charges = None
        self._model_resp_derived_charges = None
        self._model_indices = None
        
        self._total_charge = 0
        
        self._from_bohr = 0.529177249
        self._len_i = 0 #Number of ESP Positions
        self._len_j = 0 #Number of Atomic Positions
        
        self._convergence_threshold = 0.001
        self._restraint_strength = 0.01
        self._tightness = 0.1
        self._max_iterations = 20
        
    def log(self, message):
        print(message)
        sys.stdout.flush()
        
    @property
    def atomic_positions(self):
        return self._atomic_positions
    @atomic_positions.setter
    def atomic_positions(self, atomic_positions):
        self._atomic_positions = np.array(atomic_positions, dtype=float)
        self.model_indices = np.array(range(len(self._atomic_positions)))
        self.model_positions = self._atomic_positions

    @property
    def atom_potentials(self):
        return self._atom_potentials
    @atom_potentials.setter
    def atom_potentials(self, atom_potentials):
        self._atom_potentials = np.array(atom_potentials, dtype=float)

    @property
    def esp_positions(self):
        return self._esp_positions
    @esp_positions.setter
    def esp_positions(self, esp_positions):
        self._esp_positions = np.array(esp_positions, dtype=float)
        self._len_i = len(self._esp_positions)

    @property
    def qm_esps(self):
        return self._qm_esps
    @qm_esps.setter
    def qm_esps(self, qm_esps):
        self._qm_esps = np.array(qm_esps, dtype=float)

    @property
    def mm_esps(self):
        return self._mm_esps
    @mm_esps.setter
    def mm_esps(self, mm_esps):
        self._mm_esps = np.array(mm_esps, dtype=float)

    @property
    def esp_derived_charges(self):
        return self._esp_derived_charges
    @esp_derived_charges.setter
    def esp_derived_charges(self, esp_derived_charges):
        self._esp_derived_charges = np.array(esp_derived_charges, dtype=float)
        self.total_charge = sum(self._esp_derived_charges)

    @property
    def model_esp_derived_charges(self):
        return self._model_esp_derived_charges
    @model_esp_derived_charges.setter
    def model_esp_derived_charges(self, model_esp_derived_charges):
        self._model_esp_derived_charges = np.array(model_esp_derived_charges, dtype=float)
        self.model_resp_derived_charges = self._model_esp_derived_charges.copy()

    @property
    def resp_derived_charges(self):
        return self._resp_derived_charges
    @resp_derived_charges.setter
    def resp_derived_charges(self, resp_derived_charges):
        self._resp_derived_charges = np.array(resp_derived_charges, dtype=float)

    @property
    def model_resp_derived_charges(self):
        return self._model_resp_derived_charges
    @model_resp_derived_charges.setter
    def model_resp_derived_charges(self, model_resp_derived_charges):
        self._model_resp_derived_charges = np.array(model_resp_derived_charges, dtype=float)
        self.total_charge = sum(self._model_resp_derived_charges)

    @property
    def model_indices(self):
        return self._model_indices
    @model_indices.setter
    def model_indices(self, model_indices):
        self._model_indices = np.array(model_indices, dtype=int)
        self._len_j = len(self._model_indices)
        self.model_positions_from_indices()
        self.model_esp_charges_from_indices()

    @property
    def model_positions(self):
        return self._model_positions
    @model_positions.setter
    def model_positions(self, model_positions):
        self._model_positions = np.array(model_positions, dtype=float)
              
    @property
    def convergence_threshold(self):
        return self._convergence_threshold
    @convergence_threshold.setter
    def convergence_threshold(self, convergence_threshold):
        self._convergence_threshold = convergence_threshold
        
    @property
    def restraint_strength(self):
        return self._restraint_strength
    @restraint_strength.setter
    def restraint_strength(self, restraint_strength):
        self._restraint_strength = restraint_strength
        
    @property
    def tightness(self):
        return self._tightness
    @tightness.setter
    def tightness(self, tightness):
        self._tightness = tightness
        
    @property
    def max_iterations(self):
        return self._max_iterations
    @max_iterations.setter
    def max_iterations(self, max_iterations):
        self._max_iterations = int(max_iterations)
        
    @property
    def total_charge(self):
        return self._total_charge
    @total_charge.setter
    def total_charge(self, total_charge):
        self._total_charge = float(total_charge)
            
    def from_log(self, log_path):
        with open(log_path, "r") as log_obj:
            lines = log_obj.readlines()
            
        atomic_positions = []
        atom_potentials = []
        esp_positions = []
        qm_esps = []
        esp_derived_charges = []

        esp_section = False
        read_esps = False
        fsize = len(lines)
        i = 0
        while i < fsize:
            l = lines[i]
            if "Electrostatic Properties" in l:
                esp_section = True
            if esp_section:
                if "Atomic Center" in l:
                    x = float(l[32:42].strip())
                    y = float(l[42:52].strip())
                    z = float(l[52:62].strip())
                    atomic_positions.append(np.array([x, y, z]))
                elif "ESP Fit Center" in l:
                    x = float(l[32:42].strip())
                    y = float(l[42:52].strip())
                    z = float(l[52:62].strip())
                    esp_positions.append(np.array([x, y, z]))
                elif "Atom   " in l:
                    atom_potentials.append(float(l.split()[-1]))
                elif "Fit    " in l:
                    qm_esps.append(float(l.split()[-1]))
                elif read_esps:
                    if "Sum of ESP" in l:
                        read_esps = False
                    else:
                        esp_derived_charges.append(float(l.split()[-1]))
                elif "ESP charges:" in l:
                    read_esps = True
                    i += 1

            i += 1

        if len(esp_derived_charges) == 0:
            raise RuntimeError("ESP charges not found in {:s}".format(log_path))

        if len(atomic_positions) == 0:
            raise RuntimeError("Atomic Centers not found in {:s}".format(log_path))

        if len(esp_positions) == 0:
            raise RuntimeError("ESP Fit Centers not found in {:s}".format(log_path))

        if len(atom_potentials) == 0:
            raise RuntimeError("Atom Potentials not found in {:s}".format(log_path))

        if len(qm_esps) == 0:
            raise RuntimeError("ESP Potentials not found in {:s}".format(log_path))

        self.esp_derived_charges = np.array(esp_derived_charges)
        self.atomic_positions    = np.array(atomic_positions   ) / self._from_bohr
        self.esp_positions       = np.array(esp_positions      ) / self._from_bohr
        self.atom_potentials     = np.array(atom_potentials    )
        self.qm_esps             = np.array(qm_esps            )
        
        self.model_indices = []
        
    def to_resp_file(self, filename, atoms):
        """
        Write a .resp file that contains:
            atom positions, esp- and resp-derived charges.
            ESP positions and values
            connectivity
        """
        
        def get_model_neighbours(neighbours, index, model_dict):
            model_neighbours = []
            for neighbour in neighbours:
                if neighbour in model_dict:
                    model_neighbours.append(model_dict[neighbour])
            return model_neighbours

        if not filename.endswith(".resp"):
            raise Exception("RESP file must end with .resp")
        
        atom_f = "ATOM {:3s} {:12.6f} {:12.6f} {:12.6f} {:10.6f} {:10.6f} {:10.6f}\n"
        esp_f  = "ESP      {:12.6f} {:12.6f} {:12.6f} {:10.6f} {:10.6f}\n"
        conn_f = "CONN {:d} {:s}\n"

        mis = self.model_indices
        misd = {mis[j]: j for j in range(len(mis))}

        with open(filename, "w") as esp_file:
            for j in range(self._len_j):
                m = mis[j] 
                esp_file.write(
                    atom_f.format(
                        atoms[m].element,
                        self.model_positions[j, 0],
                        self.model_positions[j, 1],
                        self.model_positions[j, 2],
                        self.atom_potentials[m],
                        self.model_esp_derived_charges[j],
                        self.model_resp_derived_charges[j]
                    )
                )
            for i in range(self._len_i):
                esp_file.write(
                    esp_f.format(
                        self.esp_positions[i, 0],
                        self.esp_positions[i, 1],
                        self.esp_positions[i, 2],
                        self.qm_esps[i],
                        self.mm_esps[i]
                    )
                )

            for j in range(self._len_j):
                m = mis[j]
                ns = atoms.connectivity[m]
                esp_file.write(
                    conn_f.format(
                        j,
                        " ".join([str(n) for n in get_model_neighbours(ns, m, misd)])
                    )
                )

    def model_positions_from_indices(self):
        self.model_positions = np.zeros((self._len_j, 3))
        _model_positions_from_indices(
            self.model_positions,
            self.atomic_positions, 
            self.model_indices,
            self._len_j
        )
        
    def model_esp_charges_from_indices(self):
        self.model_esp_derived_charges = np.zeros((self._len_j))
        _model_esp_charges_from_indices(
            self.model_esp_derived_charges,
            self.esp_derived_charges, 
            self.model_indices,
            self._len_j
        )
        
    def calculate_mm_esps(self):
        
        self.mm_esps = np.zeros((self._len_i))
        _calculate_mm_esps(
            self.mm_esps,
            self.model_resp_derived_charges, 
            self.esp_positions, 
            self.model_positions,
            self._len_i,
            self._len_j
        )
        
    def check_charge_convergence(self, new_resp_charges):
        assert len(self.model_resp_derived_charges) == len(new_resp_charges)

        rms = np.zeros((self._len_j))

        for q in range(self._len_j):
            rms[q] = ((self.model_resp_derived_charges[q] - new_resp_charges[q]) ** 2) ** 0.5

        self.log("Max charge deviation: {} Average deviation: {}" .format(max(rms), np.average(rms)))
        if max(rms) < self.convergence_threshold:
            return True
        return False


    def get_new_charges(self):
        
        self.resp_derived_charges = np.zeros(len(self.atomic_positions))

        B = np.zeros((self._len_j + 1))
        _form_B(
            B,
            self.qm_esps, 
            self.esp_positions, 
            self.model_positions, 
            self.total_charge,
            self._len_i,
            self._len_j
        )

        for step in range(self.max_iterations):

            self.log("Step {:d}".format(step + 1))

            self.calculate_mm_esps()
                
            A = np.zeros((self._len_j + 1, self._len_j + 1))
            _form_A(
                A,
                self.model_resp_derived_charges, 
                self.esp_positions, 
                self.model_positions, 
                self.restraint_strength, 
                self.tightness,
                self._len_i,
                self._len_j
            )

            A_inv = np.linalg.inv(A)

            q = np.dot(B, A_inv)

            new_charges = q[:-1]
            lagrange = q[-1]
            
            converged = self.check_charge_convergence(new_charges)

            self.log("Least Square Fit: {}. Lambda: {}".format(_get_least_square_fit(self.qm_esps, self.mm_esps, self._len_i), lagrange))
            sys.stdout.flush()

            if converged:
                self.log("Converged")
                break

            self.model_resp_derived_charges = new_charges.copy()

        for j, model_index in enumerate(self.model_indices):
            self.resp_derived_charges[model_index] = new_charges[j]

        if not converged:
            self.log("Not converged")
            
@jit(signature_or_function="(float64[:,:], float64[:,:], int64[:], int64)", nopython=True)
def _model_positions_from_indices(model_positions, atomic_positions, model_indices, len_j):
    for j in range(len_j):
        mi = model_indices[j]
        model_positions[j,0] = atomic_positions[mi,0]
        model_positions[j,1] = atomic_positions[mi,1]
        model_positions[j,2] = atomic_positions[mi,2]

@jit(signature_or_function="(float64[:], float64[:], int64[:], int64)", nopython=True)
def _model_esp_charges_from_indices(esp_model_charges, esp_atomic_charges, model_indices, len_j):
    for j in range(len_j):
        mi = model_indices[j]
        esp_model_charges[j] = esp_atomic_charges[mi]
        

@jit(signature_or_function="float64(float64[:], float64[:])", nopython=True)
def _distance(p_i, p_j):
    r_ij = (
        (p_i[0] - p_j[0]) ** 2 +
        (p_i[1] - p_j[1]) ** 2 +
        (p_i[2] - p_j[2]) ** 2
    ) ** 0.5
    return r_ij

@jit(signature_or_function="(float64[:], float64[:], float64[:,:], float64[:,:], int64, int64)", nopython=True)
def _calculate_mm_esps(mm_esps, esp_atomic_potentials, esp_positions, atomic_positions, len_i, len_j):
    q_js = esp_atomic_potentials
    p_is = esp_positions
    p_js = atomic_positions
    for i in range(len_i):
        for j in range(len_j):
            q = q_js[j]
            r_ij = (
                (p_is[i,0] - p_js[j,0]) ** 2 +
                (p_is[i,1] - p_js[j,1]) ** 2 +
                (p_is[i,2] - p_js[j,2]) ** 2
            ) ** 0.5

            mm_esps[i] += q / r_ij

@jit(signature_or_function="float64(float64[:], float64[:], int64)", nopython=True)
def _get_least_square_fit(arr1, arr2, len_arr):
    lsf = 0.
    for i in range(len_arr):
        lsf += (arr1[i] - arr2[i]) ** 2
        
    return lsf

@jit(signature_or_function="float64(float64[:,:], float64[:,:], int64, int64, int64)", nopython=True)
def _form_a_jk(p_is, p_js, len_i, j, k):
    sum_i = 0.
    for i in range(len_i):
        r_ij = (
            (p_is[i,0] - p_js[j,0]) ** 2 +
            (p_is[i,1] - p_js[j,1]) ** 2 +
            (p_is[i,2] - p_js[j,2]) ** 2
        ) ** 0.5

        r_ik = (
            (p_is[i,0] - p_js[k,0]) ** 2 +
            (p_is[i,1] - p_js[k,1]) ** 2 +
            (p_is[i,2] - p_js[k,2]) ** 2
        ) ** 0.5
        
        sum_i += (1. / (r_ij * r_ik))
        
    return sum_i

@jit(signature_or_function="float64(float64[:,:], float64[:,:], int64, int64)", nopython=True)
def _form_a_jj(p_is, p_js, len_i, j):
    sum_i = 0.
    for i in range(len_i):
        
        r_ij = (
            (p_is[i,0] - p_js[j,0]) ** 2 +
            (p_is[i,1] - p_js[j,1]) ** 2 +
            (p_is[i,2] - p_js[j,2]) ** 2
        ) ** 0.5
        sum_i += 1. / (r_ij * r_ij)
        
    return sum_i

@jit(signature_or_function="(float64[:,:], float64[:], float64[:,:],  float64[:,:],    float64,      float64,  int64, int64)", nopython=True)
def _form_A(                 A,           old_charges, esp_positions, model_positions, scale_factor, tightness, len_i, len_j):
    
    # 
    # Christopher I. Bayly, Piotr Cieplak, Wendy D. Cornell, and Peter A. Kollman
    # J. Phys. Chem. 1993, 97, 10269-10280
    #
    # Brent H. Besler, Kenneth M. Merz, Jr. and Peter A. Kollman
    # J. Comp. Chem, 1990, 11, 431-439
    #
    # Aq = B
    #
    # | A11  A12 . . A1n  1 | |   q1   |   |   B1   | 
    # | A21  A22 . . A2n  1 | |   q2   |   |   B2   |
    # |  .    .  .    .     | |   .    |   |   .    |
    # |  .    .   .   .   1 | |   .    | = |   .    |
    # |  .    .    .  .     | |   .    |   |   .    |
    # | An1  An2 . . Ann  1 | |   qn   |   |   Bn   |
    # |  1    1   1   1   0 | | lambda |   |  q_tot |
    #
    #  Where Ajk (j!=k) = sum(i){1 / (r_ij * r_ik)}
    #        Ajk (j==k) = sum(i){1 / (r_ij ** 2)} + delta(chi_rstr**2)/delta(q_j)
    #
    #        delta(chi_rstr**2)/delta(q_j) = a * q_j / ((q_j ** 2 + b ** 2) ** 0.5)
    #        a = Scale factor / restraint strength
    #        b = Tightness of restraint hyperbola
    #
    #        Bj = sum(i){qm_esps[i] / r_ij}
    # 
    
    q_js = old_charges
    p_is = esp_positions
    p_js = model_positions
    a = scale_factor
    b = tightness
    
    for j in range(len_j + 1):
        for k in range(len_j + 1):
                
            if j > k:
                if j == len_j:
                    A[j,k] = A[k,j] = 1.
                else:
                    A[j,k] = A[k,j] = _form_a_jk(p_is, p_js, len_i, j, k)
                
            elif j == k:
                if j == len_j:
                    A[j,j] = 0.
                else:
                    q_j = q_js[j]
                    #A[j,j] = _form_a_jj(p_is, p_js[j], len_i) + a * q_j / ((q_j ** 2 + b ** 2) ** 0.5)
                    A[j,j] = _form_a_jj(p_is, p_js, len_i, j) + a / ((q_j ** 2 + b ** 2) ** 0.5)

@jit(signature_or_function="float64(float64[:,:], float64[:,:], int64, int64, float64[:])", nopython=True)
def _form_b_j(p_is, p_js, len_i, j, V_is):
    
    b_j = 0.
    for i in range(len_i):
        V_i = V_is[i]
        r_ij = (
            (p_is[i,0] - p_js[j,0]) ** 2 +
            (p_is[i,1] - p_js[j,1]) ** 2 +
            (p_is[i,2] - p_js[j,2]) ** 2
        ) ** 0.5
        b_j += V_i / (r_ij)

    return b_j

@jit(signature_or_function="(float64[:], float64[:], float64[:,:], float64[:,:], float64, int64, int64)", nopython=True)
def _form_B(B, qm_esps, esp_positions, model_positions, total_charge, len_i, len_j):
    
    V_is = qm_esps
    p_is = esp_positions
    p_js = model_positions
    q_tot = total_charge
    
    for j in range(len_j):
        b_ij = _form_b_j(p_is, p_js, len_i, j, V_is)
        B[j] =  b_ij
        
    B[len_j] = q_tot

