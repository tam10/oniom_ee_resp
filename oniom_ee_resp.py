#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 10 11:12:32 2018

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

from amber import Parameters, StretchParam, BendParam, ImproperTorsionParam
from atoms import Atoms
from resp_tools import RESP_Optimiser

import numpy as np
from subprocess import Popen, PIPE
import time
import shutil
import traceback as tb
import os

class ONIOM_Optimiser():
    def __init__(self, atoms): 
        
        self.settings = _ONIOM_Settings()
        
        self.resp_opt = RESP_Optimiser()
        self.resp_opt.log = self.log

        self._initialise_log()

        #Check atoms
        if type(atoms) != Atoms:
            raise RuntimeError("atoms must be a Protein Extensions Atoms instance")
        self.initial_atoms = self.atoms = atoms
        self.size = len(atoms)
        
        self.params = Parameters()

    def set_model_region(self, model_indices):

        if type(model_indices) != list:
            raise RuntimeError("model_indices must be a list")
        if len(model_indices) == 0:
            raise RuntimeError("model_indices is empty")
        self.initial_model_indices = model_indices
        self.get_model_region()

        #Have to update this if model region changed
        self._get_atom_info_strings()
     
    def set_params(self, params):
        
        if not self.links:
            raise RuntimeError("Model region not set")

        if params is not None:
            if not isinstance(params, Parameters):
                raise RuntimeError("'params' must be a Parameters object")
        self.params = self.initial_params = params
        
    def _initialise_log(self):
        with open(self.settings.global_settings.home + "/" + self.settings.global_settings.log_fn, "w") as log_obj:
            log_obj.write("{}: {}\n".format(time.ctime(), "Initialised"))

    def log(self, message):
        with open(self.settings.global_settings.home + "/" + self.settings.global_settings.log_fn, "a") as log_obj:
            log_obj.write("{}: {}\n".format(time.ctime(), message))

    def copy_from_node(self, fn, destination):
        try:
            shutil.copy(fn, destination)
            self.log("Copied {} to {}/{}".format(fn, destination, fn))
        except:
            self.log("Failed to copy {} to {}/{}\n{}".format(fn, destination, fn, tb.format_exc()))
            
    def copy_to_node(self, fn):
        try:
            shutil.copy(fn, "./")
            self.log("Copied {} to node".format(fn))
        except:
            self.log("Failed to copy {} to node\n{}".format(fn, tb.format_exc()))
        
    def _run_gauss(self, base_fn, copy_com = False, copy_log = False, copy_chk = False):
       
        global_settings = self.settings.global_settings
        
        if copy_com:
            self.copy_from_node(base_fn + ".com", global_settings.home)

        process = Popen(["{} < {}.com > {}.log".format(global_settings.gaussian_run_command, base_fn, base_fn)], stdout=PIPE, stderr=PIPE, shell=True)
        out, err = process.communicate()
        
        if err:
            raise Exception("Failed to run: " + err)

        if copy_log:
            self.copy_from_node(base_fn + ".log", global_settings.work)
        
        if copy_chk:
            self.copy_from_node(base_fn + ".chk", global_settings.work)
            
    def update_resps_from_log(self, log_fn, write_resp_file=False):
        
        global_settings = self.settings.global_settings
        resp_settings = self.settings.resp_settings
        
        self.log("Calculating RESP charges from: {}.log".format(log_fn))
        self.resp_opt.from_log(log_fn + ".log")

        models_and_links = self.model_indices + self.links.keys()
        self.resp_opt.model_indices = models_and_links
        
        log_str = "RESP on {} atoms using {} ESP Points\n".format(len(self.resp_opt.model_indices), self.resp_opt._len_i)

        self.resp_opt.total_charge = global_settings.charges[-1]
        self.resp_opt.convergence_threshold = resp_settings.resp_convergence_threshold
        self.resp_opt.restraint_strength = resp_settings.resp_restraint_strength
        self.resp_opt.tightness = resp_settings.resp_restraint_tightness
        self.resp_opt.max_iterations = resp_settings.resp_max_iterations
        self.resp_opt.get_new_charges()

        esp_charges = self.resp_opt.esp_derived_charges
        resp_charges = self.resp_opt.resp_derived_charges

        #Apply RESP charges
        for j in models_and_links:
            self.atoms[j].partial_charge = resp_charges[j]
            self.atoms[j].resp_charge = resp_charges[j]
            self.atoms[j].esp_charge = esp_charges[j]

        #Apply link atom restraints
        link_rcd   = resp_settings.resp_charge_distribution["layer1"]
        layer1_rcd = resp_settings.resp_charge_distribution["layer1"]
        layer2_rcd = resp_settings.resp_charge_distribution["layer2"]
        layer3_rcd = resp_settings.resp_charge_distribution["layer3"]

        for link_i, link in self.links.items():
            link_atom = self.atoms[link_i]
            link_charge = link_atom.partial_charge

            layer1 = link["layer1"]
            layer2 = link["layer2"]
            layer3 = link["layer3"]

            layer1_size = len(layer1)
            layer2_size = len(layer2)
            layer3_size = len(layer3)

            layer1_charge_per_atom = link_charge * layer1_rcd / float(layer1_size)
            layer2_charge_per_atom = link_charge * layer2_rcd / float(layer2_size)
            layer3_charge_per_atom = link_charge * layer3_rcd / float(layer3_size)

            log_str += "\nRedistributing charge from link atom (index: {}, charge {:7.4f}):\n".format(link_i, link_charge)

            log_str += "Layer 1: Weight: {:7.4f}, Distribution: {:7.4f}\n".format(layer1_rcd, layer1_charge_per_atom)
            for layer1_i in layer1:
                log_str += "Atom {:5d}. Old charge: {:7.4f}, New charge {:7.4f}\n".format(
                    layer1_i, 
                    self.atoms[layer1_i].partial_charge,
                    self.atoms[layer1_i].partial_charge + layer1_charge_per_atom
                )
                
                self.atoms[layer1_i].partial_charge += layer1_charge_per_atom

            log_str += "Layer 2: Weight: {:7.4f}, Distribution: {:7.4f}\n".format(layer2_rcd, layer2_charge_per_atom)
            for layer2_i in layer2:
                log_str += "Atom {:5d}. Old charge: {:7.4f}, New charge {:7.4f}\n".format(
                    layer2_i, 
                    self.atoms[layer2_i].partial_charge,
                    self.atoms[layer2_i].partial_charge + layer2_charge_per_atom
                )
                
                self.atoms[layer2_i].partial_charge += layer2_charge_per_atom

            log_str += "Layer 3: Weight: {:7.4f}, Distribution: {:7.4f}\n".format(layer3_rcd, layer3_charge_per_atom)
            for layer3_i in layer3:
                log_str += "Atom {:5d}. Old charge: {:7.4f}, New charge {:7.4f}\n".format(
                    layer3_i, 
                    self.atoms[layer3_i].partial_charge,
                    self.atoms[layer3_i].partial_charge + layer3_charge_per_atom
                )
                
                self.atoms[layer3_i].partial_charge += layer3_charge_per_atom

            link_atom.partial_charge *= link_rcd

        log_str += "RESP Calculation complete"

        self.log(log_str)

        if write_resp_file:
            self.resp_opt.to_resp_file(log_fn + ".resp", self.atoms)
            
    def _update_geometry_from_log(self, log_fn):
        with open(log_fn + ".log", "r") as log_obj:
            log_str = log_obj.read()
            
        positions_lines = log_str.split('Center     Atomic      Atomic             Coordinates (Angstroms)')[-1].split('--------------\n')[1].split('\n')[:-1]
        try:
            positions = [[float(p) for p in pl.split()[-3:]] for pl in positions_lines]
        except ValueError:
            raise Exception("No geometry in {}".format(log_fn))

        self.atoms.positions = np.array(positions)
        
    def write_com(self, base_fn, keywords_string, title, oldchk_fn = "", write_connectivity = True, write_params = True):
        
        global_settings = self.settings.global_settings
        
        with open(base_fn + ".com", "w") as com_obj:
            com_obj.write("%mem={:d}MB\n".format(global_settings.mem_mb))
            com_obj.write("%nprocshared={:d}\n".format(global_settings.ncpus))
            com_obj.write("%chk={:s}.chk\n".format(base_fn))
            if global_settings._last_calc:
                com_obj.write("%oldchk={:s}.chk\n".format(global_settings._last_calc))
            
            com_obj.write("{}\n\n".format(keywords_string))
            com_obj.write("{}\n\n".format(title))
            com_obj.write(" ".join(["{:d} {:d}".format(global_settings.charges[i], global_settings.multiplicities[i]) for i in range(len(global_settings.charges))]) + "\n")

            com_obj.write(self._get_atoms_string() + "\n")
            
            if write_connectivity:
                com_obj.write("\n".join(self._connectivity_strings) + "\n\n\n")
                
            if write_params:
                com_obj.write(self._get_params_string() + "\n")
        
        
    def get_resp(self, copy_com = False, copy_log = False, copy_chk = False, copy_resp = False):
        
        global_settings = self.settings.global_settings
        resp_settings = self.settings.resp_settings
        
        base_fn = "{:s}{:s}{:d}".format(global_settings.base_name_prefix, resp_settings.base_name, resp_settings.step)
        title = "ESP Calculation {:d}".format(resp_settings.step)
        
        keywords = resp_settings.keywords.copy()
        keywords.update(global_settings.additional_keywords)
        keywords.update(resp_settings.additional_keywords)
        if global_settings._last_calc:
            keywords.update({"guess": "read"})
            
        write_connectivity = "connectivity" in keywords["geom"]
        
        iops = resp_settings.iops.copy()
        
        keywords_string = self._get_keywords_string(keywords, iops, embed = resp_settings.use_embed)
        oldchk_fn = "%oldchk={:s}.chk\n".format(global_settings._last_calc) if global_settings._last_calc else ""
        
        self.write_com(base_fn, keywords_string, title, oldchk_fn, write_connectivity)
        
        self._run_gauss(base_fn, copy_com, copy_log, copy_chk)
        self.update_resps_from_log(base_fn, copy_resp)
        if copy_resp:
            self.copy_from_node(base_fn + ".resp", self._work)
        
        resp_settings.step += 1
        global_settings._last_calc = base_fn
    
    def optimise(self, copy_com = False, copy_log = False, copy_chk = False):
        
        global_settings = self.settings.global_settings
        opt_settings = self.settings.opt_settings
        
        base_fn = "{:s}_{:d}".format(opt_settings.base_name, opt_settings.step)
        title = "Optimisation {:d}".format(opt_settings.step)
        
        keywords = opt_settings.keywords.copy()
        if not "opt" in keywords:
            keywords["opt"] = dict()
            
        keywords.update(global_settings.additional_keywords)
        keywords.update(opt_settings.additional_keywords)
        keywords["opt"]["maxcycles"] = opt_settings.max_steps

        recompute_hessian = False
        
        if opt_settings.steps_per_hessian != 0:
            if opt_settings.step % opt_settings.steps_per_hessian == 0:
                recompute_hessian = True

        if opt_settings.step == 0:
            recompute_hessian = True

        if opt_settings.restart:
            recompute_hessian = True
            opt_settings.restart = False

        if recompute_hessian:
            keywords["opt"]["calcfc"] = None
        else:
            keywords["opt"]["readfc"] = None
            keywords.update({"guess": "read"})
            
        write_connectivity = "connectivity" in keywords["geom"]
        
        iops = opt_settings.iops.copy()
        iops.update({
            1: {
                139: opt_settings.max_mm_steps,
                8: opt_settings.step_size
            }
        })
        
        keywords_string = self._get_keywords_string(keywords, iops)
        oldchk_fn = "" if recompute_hessian else "%oldchk={:s}_{:d}.chk\n".format(opt_settings.base_name, opt_settings.step - 1)
        
        self.write_com(base_fn, keywords_string, title, oldchk_fn, write_connectivity)
            
        self._run_gauss(base_fn, copy_com, copy_log, copy_chk)
        self._update_geometry_from_log(base_fn)

        with open("{:s}.log".format(base_fn), "r") as log_obj:
            last_line = log_obj.readlines()[-1]
            if "Normal termination" in last_line:
                opt_settings.converged = True
   
        opt_settings.step += 1

    def get_missing_parameters(self, copy_com = False, copy_log = False, copy_chk = False):
        
        global_settings = self.settings.global_settings
        params_settings = self.settings.params_settings
        
        base_fn = "{:s}_{:d}".format(params_settings.base_name, params_settings.step)
        title = "Parameter determination step {:d}".format(params_settings.step)
        
        keywords = params_settings.keywords.copy()
        keywords.update(global_settings.additional_keywords)
        keywords.update(params_settings.additional_keywords)
        
        iops = params_settings.iops.copy()
        
        write_connectivity = "connectivity" in keywords["geom"]
        
        keywords_string = self._get_keywords_string(keywords, iops, oniom=False, low_method=True)
        
        self.write_com(base_fn, keywords_string, title, "", write_connectivity, False)
            
        self._run_gauss(base_fn, copy_com, copy_log, copy_chk)
        self._update_params_from_log(base_fn)
        
        params_settings.step += 1
        
    def irc_optimise(self, copy_com = False, copy_log = False, copy_chk = False):
        
        global_settings = self.settings.global_settings
        irc_settings = self.settings.params_settings

        base_fn = "{:s}_{:d}".format(irc_settings.base_name, irc_settings.step)
        title = "IRC Optimisation {:d}".format(irc_settings.step)

        keywords = irc_settings.keywords.copy()
        keywords.update(global_settings.additional_keywords)
        keywords.update(irc_settings.additional_keywords)

        keywords["irc"]["maxpoints"] = irc_settings.max_steps
        if irc_settings.step == 0:
            keywords["irc"]["calcfc"] = None
        else:
            if "calcfc" in keywords["irc"]:
                del keywords["irc"]["calcfc"]
            keywords["irc"]["rcfc"] = None
            keywords.update({"guess": "read"})
        
        iops = irc_settings.iops.copy()
        iops.update({
            1: {
                139: irc_settings.max_mm_steps,
                8: irc_settings.step_size
                }
        })
        
        write_connectivity = "connectivity" in keywords["geom"]
        
        keywords_string = self._get_keywords_string(keywords, iops)
        oldchk_fn = "%oldchk={:s}_{:d}.chk\n".format(irc_settings.base_name, irc_settings.step - 1) if irc_settings.step > 0 else ""
        
        self.write_com(base_fn, keywords_string, title, oldchk_fn, write_connectivity)
            
        self._run_gauss(base_fn, copy_com, copy_log, copy_chk)
        self._update_geometry_from_log(base_fn)

        with open("{:s}_{:d}.log".format(irc_settings.base_name, irc_settings.step), "r") as log_obj:
            lines = log_obj.readlines()
            for line in lines:
                if "PES minimum detected" in line:
                    irc_settings.converged = True
   
        irc_settings.step += 1

    def initialise(self):
    
        sum_charge_distribution = sum([v for k,v in self.settings.resp_settings.resp_charge_distribution.items()])

        if not np.abs(sum_charge_distribution - 1) < 0.0001:
            self.log("WARNING: Sum of charge distribution is not 1. Partial charge will not be conserved")

        if not self.charges or not self.multiplicities:
            raise Exception("Charges and/or multiplicities not set")

        if not self.settings.global_settings.high_method:
            raise Exception("High method not set")

        if not self.settings.global_settings.low_method:
            raise Exception("Low method not set")

        self._get_atom_info_strings()
    
    def run(self, copy_com = False, copy_log = False, copy_chk = False, copy_resp = False):
        
        self.initialise()

        if self.settings.global_settings.irc_mode == True:
            for step_num in range(self.settings.global_settings.max_overall_steps):
                self.get_resp(copy_com, copy_log, copy_chk, copy_resp)
                self.irc_optimise(copy_com, copy_log, copy_chk)
        else:
            for step_num in range(self.settings.global_settings.max_overall_steps):
                self.get_resp(copy_com, copy_log, copy_chk, copy_resp)
                self.optimise(copy_com, copy_log, copy_chk)
        
    def get_model_region(self):
        #Get the model region indices depending on whether:
        #    the initial indices were tags (in the case selection from non-protonated protein)
        #    the initial indices are the real indices
        #
        #Expand selection to include any neighbouring hydrogens and link atoms
        #Determine types of link atoms
        
        if not self.initial_model_indices:
            raise Exception("Model indices are not set")
            
        self.links = dict()
        self.model_indices = model_indices = self.initial_model_indices

        neighbours = self.atoms.expand_selection_by_bonds(model_indices, include_current_selection = False)
        h_neighbours  = [n for n in neighbours if self.atoms[n].element == "H"]
        links = [n for n in neighbours if n not in h_neighbours]

        for link_index in links:
        
            #Get atoms in the model region that neighbour the link atoms for RESP
            layer1 = self.atoms.expand_selection_by_bonds([link_index], include_current_selection = False)
            layer1 = [n for n in layer1 if n in model_indices]
            if not len(layer1) == 1:
                raise Exception("Link atom {} has wrong number of model region neighbours: {} ({})".format(link_index, len(layer1), layer1))
            
            link_type = self.settings.params_settings.link_type_dict[self.atoms[layer1[0]].element]
            if type(link_type) == dict:
                try:
                    link_type = link_type[self.atoms[layer1[0]].amber_type]
                except:
                    link_type = link_type["*"]

            layer2 = self.atoms.expand_selection_by_bonds([link_index], include_current_selection = False)
            layer2 = [n for n in layer2 if n in model_indices and n not in layer1]
            
            layer3 = self.atoms.expand_selection_by_bonds([link_index], include_current_selection = False)
            layer3 = [n for n in layer3 if n in model_indices and n not in layer1 and n not in layer2]
            

            self.links[link_index] = {"type": link_type, "layer1": layer1, "layer2": layer2, "layer3": layer3}

    def set_high_method(self, method, basis=''):
        self.settings.global_settings.set_high_method(method, basis)

    def set_low_method(self, method, basis=''):
        self.settings.global_settings.set_low_method(method, basis)
        
    def _get_low_method_keyword(self):
        low_method_keyword = "{}{}".format(
            self.settings.global_settings.low_method,
            "/{}".format(self.settings.global_settings.low_basis) if self.settings.global_settings.low_basis else ""
        )
        return low_method_keyword
        
    def _get_high_method_keyword(self):
        high_method_keyword = "{}{}".format(
            self.settings.global_settings.high_method,
            "/{}".format(self.settings.global_settings.high_basis) if self.settings.global_settings.high_basis else ""
        )
        return high_method_keyword
            
    def _get_oniom_keyword(self, embed=False):
        oniom_str = "oniom({}:{}){}".format(
            self._get_high_method_keyword(),
            self._get_low_method_keyword(),
            "=embed" if embed else ""
        )
        return oniom_str
        
    def _get_iop_string(self, iops):
        strings = []
        for overlay, options in iops.items():
            for option, value in options.items():
                strings.append("{:d}/{:d}={:d}".format(overlay, option, value))
        return "iop({})".format(",".join(strings))
        
    def _get_keywords_string(self, keywords_dict, iops, embed=False, oniom=True, low_method=True):
    
        #Keywords is supplied as a dictionary to allow easy replacement and updating
        #keywords_dict = {
        #    "keyword1":  
        #        {
        #            "option1": None,
        #            "option2": "value1",
        #        },
        #    "keyword2": "value2"
        #    "keyword3": None
        #}
        #
        # is Gaussian equivalent of:
        # keyword1=(option1,option2=value1) keyword2=value2 keyword3
        #
        #
        # Setting 'oniom' to True returns the full ONIOM keyword, and 'low_method' is ignored
        # If 'oniom' is False, setting 'low_method' to True returns only the low method keyword
        # If 'oniom' is False, setting 'low_method' to False returns only the high method keyword

        strings = []

        if self.settings.global_settings.additional_print:
            strings.append("#p")
        else:
            strings.append("#")

        for keyword, options in keywords_dict.items():
            keyword_string = "{}".format(keyword)

            if options != None:
                keyword_string += "="

                if type(options) == dict:
                    keyword_string += "("

                    options_strings = []

                    for option, value in options.items():
                        if value:
                            options_strings.append("{}={}".format(option, value))
                        elif value == None:
                            options_strings.append("{}".format(option))
                        else:
                            raise Exception("keywords not properly formatted")

                    keyword_string += ",".join(options_strings) + ")"
                else:
                    keyword_string += "{}".format(options)

            strings.append(keyword_string)

        if oniom:
            keywords_string = " ".join(strings + [self._get_oniom_keyword(embed)]) +  "\n" + self._get_iop_string(iops)
        elif low_method:
            keywords_string = " ".join(strings + [self._get_low_method_keyword()]) +  "\n" + self._get_iop_string(iops)
        else:
            keywords_string = " ".join(strings + [self._get_high_method_keyword()]) +  "\n" + self._get_iop_string(iops)
            
        return keywords_string
    
    def _get_atom_info_strings(self):
        #ONIOM(QM:AMBER) atom line
        #Build a string with formatting gaps
        #Charge gaps (f1) to be filled by RESP calc
        #Position gaps (f3) to be filled by Opt calc
        
        #Atoms formatting
        f0 = " {:s}-{:s}-"
        f1 = "{:.6f}"
        f2 = "(PDBName={:s},ResName={:s},ResNum={:d})"
        f3 = " 0 {:13.8f} {:13.8f} {:13.8f} "
        f4 = "{:s} {:s}-{:s} {:d}"
        
        link_keys = self.links.keys()
        
        self._atom_info_strings = []
        self._connectivity_strings = []
        
        for i in range(self.size):
            a = self.atoms[i]
            
            s0 = f0.format(a.element, a.amber_type)
            s2 = f2.format(a.pdb_type, a.residue_name, a.residue_number)
            
            if i in self.model_indices:
                s4 = "H"
            elif i in link_keys:
                l = self.links[i]
                s4 = f4.format("L", "H", l["type"], l["layer1"][0] + 1)
            else:
                s4 = "L"
                
            self._atom_info_strings.append([s0 + f1 + s2, f3 + s4])
            
            s = " {:d}".format(i + 1)
            ns = self.atoms.connectivity[i]
            
            for n in ns:
                s += " {:d} 1.0".format(n + 1)
                
            self._connectivity_strings.append(s)
            
    
    def _get_atoms_string(self):
        #Build the entire atoms section of com file by filling in gaps in atoms_info_strings
        
        atoms_string = ""
        
        f2 = "{:52s}"

        for i in range(self.size):
            
            a = self.atoms[i]
            f0, f1 = self._atom_info_strings[i]
            p = a.position
            
            atoms_string += f2.format(f0.format(a.amber_charge)) + f1.format(p[0],p[1],p[2]) + "\n"
            
        return atoms_string


    def _update_params_from_log(self, log_fn):

        def get_torsion_energy(torsionV, theta, phase_offset, nodes):
            return torsionV * (1 + np.cos(np.deg2rad(((nodes + 1)*theta-phase_offset))))

        def append_length_types_from_ambers(a0, a1, lengthAtoms):

            key1 = "{}-{}".format(a0, a1)
            key2 = "{}-{}".format(a1, a0)

            key = key1 if key1 > key2 else key2

            if lengthAtoms.get(key):
                lengthAtoms[key].append([i0, i1])
            else:
                lengthAtoms[key] = [[i0, i1]]

        def append_length_types_from_indices(atoms, i0, i1, lengthAtoms, links):

            a0 = atoms[i0].amber_type
            a1 = atoms[i1].amber_type

            append_length_types_from_ambers(a0, a1, lengthAtoms)
            if i0 in links:
                link_type = links[i0]["type"]
                append_length_types_from_ambers(link_type, a1, lengthAtoms)
            if i1 in links:
                #link_type = links[i0]["type"]
                link_type = links[i1]["type"]
                append_length_types_from_ambers(a0, link_type, lengthAtoms)

        def append_angle_types_from_ambers(a0, a1, a2, angleAtoms):

            key1 = "{}-{}-{}".format(a0, a1, a2)
            key2 = "{}-{}-{}".format(a2, a1, a0)

            key = key1 if key1 > key2 else key2

            if angleAtoms.get(key):
                angleAtoms[key].append([i0, i1, i2])
            else:
                angleAtoms[key] = [[i0, i1, i2]]

        def append_angle_types_from_indices(atoms, i0, i1, i2, angleAtoms, links):

            a0 = atoms[i0].amber_type
            a1 = atoms[i1].amber_type
            a2 = atoms[i2].amber_type

            append_angle_types_from_ambers(a0, a1, a2, angleAtoms)
            if i0 in links:
                link_type = links[i0]["type"]
                append_angle_types_from_ambers(link_type, a1, a2, angleAtoms)
            if i2 in links:
                link_type = links[i2]["type"]
                append_angle_types_from_ambers(a0, a1, link_type, angleAtoms)

        def append_torsion_types_from_ambers(a0, a1, a2, a3, torsionAtoms):

            key1 = "{}-{}-{}-{}".format(a0, a1, a2, a3)
            key2 = "{}-{}-{}-{}".format(a3, a2, a1, a0)

            key = key1 if key1 > key2 else key2

            if torsionAtoms.get(key):
                if [i0, i1, i2, i3] not in torsionAtoms[key]:
                    torsionAtoms[key].append([i0, i1, i2, i3])
            else:
                torsionAtoms[key] = [[i0, i1, i2, i3]]

        def append_torsion_types_from_indices(atoms, i0, i1, i2, i3, torsionAtoms, links):

            a0 = atoms[i0].amber_type
            a1 = atoms[i1].amber_type
            a2 = atoms[i2].amber_type
            a3 = atoms[i3].amber_type

            append_torsion_types_from_ambers(a0, a1, a2, a3, torsionAtoms)
            if i0 in links:
                link_type = links[i0]["type"]
                append_torsion_types_from_ambers(link_type, a1, a2, a3, torsionAtoms)
            if i3 in links:
                link_type = links[i3]["type"]
                append_torsion_types_from_ambers(a0, a1, a2, link_type, torsionAtoms)

        atoms = self.atoms

        with open("{:s}.log".format(log_fn), "r") as mp_obj:
            mp_lines = mp_obj.readlines()

        lengthKeq = str(self._lengthKeq)
        angleKeq  = str(self._angleKeq)
        torsionV  = str(self._torsionV)

        lengthAtoms  = {}
        angleAtoms  = {}
        torsionAtoms = {}

        lengthIndices = []
        angleIndices = []
        torsionIndices = []

        #Read output to see what's missing in bonds and angles
        #This is only for real, so we have to fill in the links later
        for mp_line in mp_lines:
            try:
                if "Bondstretch undefined" in mp_line:
                    i0 = int(mp_line.split()[4]) - 1
                    i1 = int(mp_line.split()[5]) - 1

                    lengthIndices.append([i0, i1])


                elif "Angle bend  undefined" in mp_line:
                    i0 = int(mp_line.split()[5]) - 1
                    i1 = int(mp_line.split()[6]) - 1
                    i2 = int(mp_line.split()[7]) - 1

                    angleIndices.append([i0, i1, i2])

            except(ValueError):
                print("Could not resolve: " + mp_line)
                raise

        for link_i in [i for i,a in enumerate(atoms.links) if atoms.layers[a] == "H"]:
            i1s = atoms.expand_selection_by_bonds([link_i], include_current_selection=False)
            for i1 in i1s:
                lengthIndices.append([link_i, i1])
                i2s = atoms.expand_selection_by_bonds([i1], include_current_selection=False)
                for i2 in i2s:
                    if i2 != link_i:
                        angleIndices.append([link_i, i1, i2])


        #Get the amber atom types for all lengths
        #Sort the amber length specifications as e.g. CH-HC == HC-CH
        for i0, i1 in lengthIndices:
            append_length_types_from_indices(atoms, i0, i1, lengthAtoms, self.links)
            #Get angles involved with links
            i2s = atoms.expand_selection_by_bonds([i0], include_current_selection=False)
            for i2 in i2s:
                if i2 != i1 and i2 in self.links:
                    angleIndices.append([i2, i0, i1])

            i2s = atoms.expand_selection_by_bonds([i2], include_current_selection=False)
            for i2 in i2s:
                if i2 != i0 and i2 in self.links:
                    angleIndices.append([i0, i1, i2])

        #Get amber types for all angles
        #Torsions aren't included in the output so add them using angles
        for i0, i1, i2 in angleIndices:
            #Get all torsions from angles
            i3s = atoms.expand_selection_by_bonds([i0], include_current_selection=False)
            for i3 in i3s:
                if i3 != i1:
                    torsionIndices.append([i3, i0, i1, i2])

            i3s = atoms.expand_selection_by_bonds([i2], include_current_selection=False)
            for i3 in i3s:
                if i3 != i1:
                    torsionIndices.append([i0, i1, i2, i3])
            append_angle_types_from_indices(atoms, i0, i1, i2, angleAtoms, self.links)

        for i0, i1, i2, i3 in torsionIndices:
            append_torsion_types_from_indices(atoms, i0, i1, i2, i3, torsionAtoms, self.links)


        for k, lengthIndicesList in lengthAtoms.items():
            lengths = [atoms.get_distance(*atomIndices) for atomIndices in lengthIndicesList]

            a0, a1 = k.split('-')
            stretch = StretchParam(a0, a1, lengthKeq, np.average(lengths))
            self.params.update_stretch(stretch)

        for k, angleIndicesList in angleAtoms.items():
            angles = [atoms.get_angle(*atomIndices) for atomIndices in angleIndicesList]

            a0, a1, a2 = k.split('-')
            bend = BendParam(a0, a1, a2, angleKeq, np.average(angles))
            self.params.update_bend(bend)

        for k, dihedralIndicesList in torsionAtoms.items():
            dihedrals = [np.rad2deg(atoms.get_dihedral(*atomIndices)) for atomIndices in dihedralIndicesList]

            nodes_phase_offsets = [
                [1, 0],
                [1, 180],
                [2, 0],
                [2, 180],
                [3, 0],
                [3, 180]
            ]

            scores = []
            for nodes, phase_offset in nodes_phase_offsets:
                score = 0.
                for dihedral in dihedrals:
                    score += get_torsion_energy(float(torsionV), dihedral, phase_offset, nodes)

                scores.append(score)

            nodes, phase_offset = nodes_phase_offsets[scores.index(min(scores))]


            a0, a1, a2, a3 = k.split('-')
            torsion = ImproperTorsionParam(a0, a1, a2, a3, gamma = phase_offset, v = torsionV, nt = nodes)
            self.params.update_torsion(torsion)
            
            
            
class _ONIOM_Settings(object):
    def __init__(self):
        self.global_settings = _Global_Settings()
        self.opt_settings    = _Opt_Settings()
        self.irc_settings    = _IRC_Settings()
        self.resp_settings   = _RESP_Settings()
        self.params_settings = _Params_Settings()
    
class _Global_Settings(object):
    def __init__(self):
        
        self.max_steps = 20
        self.irc_mode = False
        self.additional_keywords = dict()
        
        self.base_name_prefix  = ""
        self._last_calc = ""
        
        #_atom_info_strings needs to be filled out with charges and positions every update
        #form it with static info (symbols, ambers, oniom etc) to speed up
        self._atom_info_strings = []
        self._connectivity_strings = []

        self.charges = []
        self.multiplicities = []

        self.mem_mb = 4000
        self.ncpus = 4
        self.overhead_mem_mb_per_cpu = 100
        
        self._set_mem_mb_ncpus_from_env()
        
        self.home = "."
        self.work = "."
        self.set_directories(None, None)
        
        self.log_fn = "optimiser.log"
        self.gaussian_run_command = "module load gaussian/g16-a03; g16"
        
        self.high_method = ""
        self.high_basis = ""
        
        self.low_method = ""
        self.low_basis = ""

    def set_high_method(self, method, basis=''):
        self.high_method = method
        self.high_basis = basis
        
    def set_low_method(self, method, basis=''):
        self.low_method = method
        self.low_basis = basis

    def set_directories(self, home, work):

        if not home:
            home = os.environ.get("PBS_O_WORKDIR")
        if not home:
            home = os.environ["PWD"]

        if not work:
            work = home.replace("/home/", "/work/")

        self.home = home
        self.work = work
        
    def _set_mem_mb_ncpus_from_env(self):
        command = 'qstat -f $PBS_JOBID'
        p = Popen([command], stdout=PIPE, stderr=PIPE, shell=True)
        out, err = p.communicate()

        if out:
            ncpus_str = [s for s in out.split("\n") if "Resource_List.ncpus" in s][0].strip().split()[-1]
            ncpus = int(ncpus_str)

            self.ncpus = ncpus

            mem_str = [s for s in out.split("\n") if "Resource_List.mem" in s][0].strip().split()[-1]
            mem_int = int("".join([c for c in mem_str if c.isdigit()]))
            mem_ord = "".join([c for c in mem_str if c.isalpha()]).upper()

            if mem_ord == "GB":
                mem_int *= 1024
            elif mem_ord == "MB":
                    pass
            else:
                mem_int /= 1024

            self.mem_mb = mem_int - self.overhead_mem_mb_per_cpu * self.ncpus
            if self.mem_mb < 0:
                raise Exception("Negative memory allocated")
    
class _Opt_Settings(object):
    def __init__(self):
        
        self.restart = False
        self.max_steps = 2
        self.max_mm_steps = 100
        self.step_size = 10
        self.steps_per_hessian = 0
        
        self.keywords  = {
            "geom": "connectivity",
            "opt": None
        }
        self.additional_keywords  = dict()
        self.iops = dict()
        self.additional_print = False
        
        self.base_name  = "opt_"
        self.step  = 0
        
        self.converged = False
    
class _IRC_Settings(object):
    def __init__(self):
        
        self.max_steps = 1
        self.max_mm_steps = 100
        self.step_size = 10
        
        self.keywords = {
            "geom": "connectivity", 
            "irc": {"downhill": None}
        }
        self.additional_keywords = dict()
        self.iops = dict()
        self.additional_print = False
        
        self.base_name  = "irc_"
        self.step  = 0
        
        self.converged = False
    
class _RESP_Settings(object):
    def __init__(self):
        
        self.keywords = {
            "geom": "connectivity",
            "scf": "tight", 
            "pop": {
                "saveesp": None,
                "mk": None
            },
            "density": "current"
        }
        self.additional_keywords = dict()
        self.iops = {
            6: {
                33: 2,
                41: 10,
                42: 17,
            }
        }
        self.additional_print = False
        
        self.base_name  = "resp_"
        self.step  = 0
        
        self.converged = False
        
        self.use_embed = True
        self.resp_convergence_threshold = 0.001
        self.resp_restraint_strength = 0.01
        self.resp_restraint_tightness = 0.1
        self.resp_max_iterations = 20
        
        self.resp_charge_distribution = {
            "link"  : 0.0,
            "layer1": 0.0,
            "layer2": 0.5,
            "layer3": 0.5
        }
        
    def set_charge_distributions(self, layer_1_weight, layer_2_weight, layer_3_weight):

        layer_1_weight = float(layer_1_weight)
        self.resp_charge_distribution["layer1"] = layer_1_weight

        layer_2_weight = float(layer_2_weight)
        self.resp_charge_distribution["layer2"] = layer_2_weight

        layer_3_weight = float(layer_3_weight)
        self.resp_charge_distribution["layer3"] = layer_3_weight
    
class _Params_Settings(object):
    def __init__(self):
                    
        self.keywords = {
            "geom": "connectivity"
        }
        self.additional_keywords = dict()
        self.iops = dict()
        self.additional_print = False
        
        self.base_name  = "params_"
        self.step  = 0
        
        self.link_type_dict = {
            "C": {
                "C" : "HC",
                "CA": "HA",
                "CM": "HA",
                "CJ": "HA",
                "CT": "HC",
                "CR": "H5",
                "CW": "H4",
                "*" : "HC"
            },
            "N": "H",
            "O": "HO",
            "S": "HS"
        }
    
        self.nonbon = "NonBon 3 1 0 0 0.000 0.000 0.500 0.000 0.000 -1.000"
        self.lengthKeq = 400.00
        self.angleKeq  = 100.00
        self.torsionV  = 4.000
        
    