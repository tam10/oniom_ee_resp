#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 10 11:11:40 2018

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

class Atoms(object):
    def __init__(
                self, 
                atoms = None,
                elements = None,
                positions = None,
                partial_charges = None,
                esp_charges = None,
                resp_charges = None,
                tags = None,
                pdb_types = None,
                amber_types = None,
                residue_names = None,
                residue_numbers = None,
                layers = None,
                connectivity = None,
                links = None
            ):
        """
    Assemble and Atoms object using a list of Atom objects or lists/arrays.
    Using lists or arrays, positions and elements are essential.
    
    partial_charges: Charges used to write Gaussian input files
    esp_charges:     ESP derived charges read from a Gaussian output file.
    resp_charges:    RESP derived charges calculated from a RESP_Optimiser
    tags:            Integers used for tagging atoms 
                        (useful for differentiating atoms)
    pdb_types:       PDB atom identifier
    amber_types:     Amber atom identifier
    residue_names:   PDB residue name
    residue_numbers: PDB residue number
    layers:          List of strings representing ONIOM Layer.
                        Must be one of ["H", "M", "L"]
    connectivity:    Reciprocal list of neighbours by index
    links:           A subset of connectivity 
                        Only link connections are included
        """
        
        self._positions = np.array([])
        self.elements = np.array([])
        self.partial_charges = np.array([])
        self.esp_charges = np.array([])
        self.resp_charges = np.array([])
        self.tags = np.array([])
        self.pdb_types = np.array([])
        self.amber_types = np.array([])
        self.residue_names = np.array([])
        self.residue_numbers = np.array([])
        self.layers = np.array([])
        self.links = np.array([])
        
        if atoms:
            if not hasattr(atoms, "__len__"):
                raise ValueError("atoms argument must have a length")
            
            positions = []
            elements = []
            partial_charges = []
            esp_charges = []
            resp_charges = []
            tags = []
            pdb_types = []
            amber_types = []
            residue_names = []
            residue_numbers = []
            layers = []
            links = []
                
            for atom in atoms:
                if not type(atom) == Atom:
                    raise ValueError("Elements in atoms argument must be Atom types")
                        
                positions.append(atom.position)
                elements.append(atom.element)
                partial_charges.append(atom.partial_charge)
                esp_charges.append(atom.esp_charge)
                resp_charges.append(atom.resp_charge)
                tags.append(atom.tag)
                pdb_types.append(atom.pdb_type)
                amber_types.append(atom.amber_type)
                residue_names.append(atom.residue_name)
                residue_numbers.append(atom.residue_number)
                layers.append(atom.layer)
                links.append(atom.link)
                
            self._positions = np.array(positions)
            self.elements = elements
            self.partial_charges = partial_charges
            self.esp_charges = esp_charges
            self.resp_charges = resp_charges
            self.tags = tags
            self.pdb_types = pdb_types
            self.amber_types = amber_types
            self.residue_names = residue_names
            self.residue_numbers = residue_numbers
            self.layers = layers
            self.links = links
            self.connectivity = connectivity

        else:
            positions = np.array(positions)
            if not positions.shape[1] == 3:
                raise ValueError("'positions' must be an array of shape (,3), not {}".format(positions.shape))
            
            self._positions = positions
            self.elements = elements
            self.partial_charges = partial_charges
            self.esp_charges = esp_charges
            self.resp_charges = resp_charges
            self.tags = tags
            self.pdb_types = pdb_types
            self.amber_types = amber_types
            self.residue_names = residue_names
            self.residue_numbers = residue_numbers
            self.layers = layers
            self.links = links
            self.connectivity = connectivity
            
    @property
    def elements(self):
        return self._elements
    @elements.setter
    def elements(self, value):
        if not len(value) == self.size:
            raise ValueError("Wrong size of 'elements' ({}) for atoms ({})".format(len(value), self.size))
        self._elements = np.array(value, dtype="str")
        
    @property
    def positions(self):
        return self._positions
    @positions.setter
    def positions(self, value):
        value = np.array(value)
        if not value.shape == (self.size, 3):
            raise ValueError("'positions' must be an array of shape ({},3), not {}".format(self.size, value.shape))
        self._positions = value
        
    @property
    def partial_charges(self):
        return self._partial_charges
    @partial_charges.setter
    def partial_charges(self, value):
        if value is None:
            self._partial_charges = np.zeros(self.size, dtype="float")
        elif not len(value) == self.size:
            raise ValueError("Wrong size of 'partial_charges' ({}) for atoms ({})".format(len(value), self.size))
        else:
            self._partial_charges = np.array(value, dtype = "float")
        
    @property
    def esp_charges(self):
        return self._esp_charges
    @esp_charges.setter
    def esp_charges(self, value):
        if value is None:
            self._esp_charges = np.zeros(self.size, dtype="float")
        elif not len(value) == self.size:
            raise ValueError("Wrong size of 'esp_charges' ({}) for atoms ({})".format(len(value), self.size))
        else:
            self._esp_charges = np.array(value, dtype = "float")
        
    @property
    def resp_charges(self):
        return self._resp_charges
    @resp_charges.setter
    def resp_charges(self, value):
        if value is None:
            self._resp_charges = np.zeros(self.size, dtype="float")
        elif not len(value) == self.size:
            raise ValueError("Wrong size of 'resp_charges' ({}) for atoms ({})".format(len(value), self.size))
        else:
            self._resp_charges = np.array(value, dtype = "float")
        
    @property
    def tags(self):
        return self._tags
    @tags.setter
    def tags(self, value):
        if value is None:
            self._tags = np.zeros(self.size, dtype="int")
        elif not len(value) == self.size:
            raise ValueError("Wrong size of 'tags' ({}) for atoms ({})".format(len(value), self.size))
        else:
            self._tags = np.array(value, dtype = "int")
        
    @property
    def pdb_types(self):
        return self._pdb_types
    @pdb_types.setter
    def pdb_types(self, value):
        if value is None:
            self._pdb_types = np.array([""] * self.size, dtype = "str")
        elif not len(value) == self.size:
            raise ValueError("Wrong size of 'pdb_types' ({}) for atoms ({})".format(len(value), self.size))
        else:
            self._pdb_types = np.array(value, dtype = "str")
        
    @property
    def amber_types(self):
        return self._amber_types
    @amber_types.setter
    def amber_types(self, value):
        if value is None:
            self._amber_types = np.array([""] * self.size, dtype = "str")
        elif not len(value) == self.size:
            raise ValueError("Wrong size of 'amber_types' ({}) for atoms ({})".format(len(value), self.size))
        else:
            self._amber_types = np.array(value, dtype = "str")
        
    @property
    def residue_names(self):
        return self._residue_names
    @residue_names.setter
    def residue_names(self, value):
        if value is None:
            self._residue_names = np.array([""] * self.size, dtype = "str")
        elif not len(value) == self.size:
            raise ValueError("Wrong size of 'residue_names' ({}) for atoms ({})".format(len(value), self.size))
        else:
            self._residue_names = np.array(value, dtype = "str")
        
    @property
    def residue_numbers(self):
        return self._residue_numbers
    @residue_numbers.setter
    def residue_numbers(self, value):
        if value is None:
            self._residue_numbers = np.zeros(self.size,dtype="int")
        elif not len(value) == self.size:
            raise ValueError("Wrong size of 'residue_numbers' ({}) for atoms ({})".format(len(value), self.size))
        else:
            self._residue_numbers = np.array(value, dtype = "int")
        
    @property
    def layers(self):
        return self._layers
    @layers.setter
    def layers(self, value):
        if value is None:
            self._layers = np.array(["H"] * self.size, dtype = "str")
        elif not len(value) == self.size:
            raise ValueError("Wrong size of 'residue_names' ({}) for atoms ({})".format(len(value), self.size))
        else:
            allowed_layer_names = Atom._allowed_layer_names
            for v in value:
                if v not in allowed_layer_names:
                    raise ValueError("Each element of 'layers' must be one of {}".format(allowed_layer_names))
            self._layers = np.array(value, dtype = "str")
        
    @property
    def links(self):
        return self._links
    @links.setter
    def links(self, value):
        if value is None:
            self._links = [[] * self.size]
        elif not len(value) == self.size:
            raise ValueError("Wrong size of 'links' ({}) for atoms ({})".format(len(value), self.size))
        else:
            size = self.size
            for index, neighbours in enumerate(value):
                for neighbour in neighbours:
                    if not neighbour >=0 and neighbour < size:
                        raise IndexError("Neighbour atom {} of atom {} is out of range of atoms ({})".format(neighbour, index, size))
            self._links = value
        
    @property
    def connectivity(self):
        return self._connectivity
    @connectivity.setter
    def connectivity(self, value):
        if value is None:
            self._connectivity = [[] * self.size]
        elif not len(value) == self.size:
            raise ValueError("Wrong size of 'connectivity' ({}) for atoms ({})".format(len(value), self.size))
        else:
            size = self.size
            for index, neighbours in enumerate(value):
                for neighbour in neighbours:
                    if not neighbour >=0 and neighbour < size:
                        raise IndexError("Neighbour atom {} of atom {} is out of range of atoms ({})".format(neighbour, index, size))
            self._connectivity = value
        
    @property
    def formula(self):
        """
        Hill system of ordering elements.
        """
        def e_count(e_list, e):
            count = e_list.count(e)
            return e + (str(count) if count > 1 else '') if count > 0 else ''
        
        es = self.elements.tolist()
        e_set = set(es)
        
        if 'C' in e_set:
            e_set.remove('C')
        if 'H' in e_set:
            e_set.remove('H')
        
        sorted_es = sorted(list(e_set))
        
        formula = e_count(es, 'C')
        formula += e_count(es, 'H')
        
        for e in sorted_es:
            formula += e_count(es, e)
            
        return formula
        
    @property
    def atomic_centre(self):
        return np.average(self.positions, axis=0)
        
    @property
    def size(self):
        return len(self.positions)
        
    def get_distance(self, index0, index1):
        p0 = self.positions[index0]
        p1 = self.positions[index1]
        return np.linalg.norm(p0 - p1)
        
    def get_angle(self, index0, index1, index2):
        p0 = self.positions[index0]
        p1 = self.positions[index1]
        p2 = self.positions[index2]

        v10 = p0 - p1
        v12 = p2 - p1
        
        c = np.dot(v10, v12) / (np.linalg.norm(v10) * np.linalg.norm(v12))
        a = np.arccos(c)
        return np.degrees(a)
        
    def get_dihedral(self, index0, index1, index2, index3):
        p0 = self.positions[index0]
        p1 = self.positions[index1]
        p2 = self.positions[index2]
        p3 = self.positions[index3]
        
        v10 = p0 - p1
        v12 = p2 - p2
        v23 = p3 - p2
        
        v10 /= np.linalg.norm(v10)
        
        w10 = v10 - v12 * np.dot(v10, v12)
        w23 = v23 - v12 * np.dot(v23, v12)
        
        x = np.dot(w10, w23)
        y = np.dot(np.cross(v12, w10), w23)
        
        a = np.arctan2(y, x)
        return np.degrees(a)
        
    def expand_selection_by_bonds(self, current_selection, expand_by = 1, include_current_selection = True, exclude = None):
        """
        Starting with a list of atomic indices, expand the list using the connectivity.
        
        expand_by:                  Number of bonds to grow selection by.
        include_current_selection:  Include the current_selection list in the output.
        exclude:                    Remove these indices at every step.
        """
        exclude = [] if exclude is None else exclude

        connectivity = self.connectivity
        expanded_selection = np.copy(current_selection)
        
        for level in range(expand_by):
            expanded_selection = [
                index for index, neighbours in enumerate(connectivity)
                if (
                    index in expanded_selection 
                    or any([(neighbour in expanded_selection) for neighbour in neighbours])
                )
                and index not in exclude                      
            ]
            
        if include_current_selection:
            return expanded_selection
        else:
            return [i for i in expanded_selection if i not in current_selection]
        
    def __len__(self):
        return self.size
        
    def __getitem__(self, i):
        """
        Returns a new Atom object for an integer input.
        Returns a new Atoms object for a list or a slice.
        """
        if isinstance(i, int):
            if i < -self.size or i > self.size:
                raise IndexError("Index out of range")
                
            atom = Atom(
                element  = self.elements[i],
                position = self.positions[i],
                partial_charge = self.partial_charges[i],
                esp_charge = self.esp_charges[i],
                resp_charge = self.resp_charges[i],
                index = i,
                tag = self.tags[i],
                pdb_type = self.pdb_types[i],
                amber_type = self.amber_types[i],
                residue_name = self.residue_names[i],
                residue_number = self.residue_numbers[i],
                layer = self.layers[i]
            )
            return atom
        elif hasattr(i, "__len__"):
            atoms = Atoms([self[index] for index in i])
            connectivity = [sorted([i.index(n) for n in self.connectivity[old_index] if n in i]) for old_index in i]
            atoms.connectivity = connectivity
            return atoms
        elif isinstance(i, slice):
            return self.__getitem__(range(self.size)[i])
            
    def __repr__(self):
        repr_str = "Atoms({}, atomic_centre={})".format(self.formula, self.atomic_centre)
        return repr_str
        
class Atom(object):
    """
    The Atom object contains the following properties:
        
    partial_charge:  Charge used to write Gaussian input files
    esp_charge:      ESP derived charge read from a Gaussian output file.
    resp_charge:     RESP derived charge calculated from a RESP_Optimiser
    tag:             Integer used for tagging
                         (useful for differentiating atoms)
    pdb_type:        PDB atom identifier
    amber_type:      Amber atom identifier
    residue_name:    PDB residue name
    residue_number:  PDB residue number
    layer:           String representing ONIOM Layer.
                         Must be one of ["H", "M", "L"]
    """
    
    #Prevent the wrong strings being used for layers
    _allowed_layer_names = ["H", "M", "L"]

    def __init__(
                 self, 
                 element,
                 position,
                 partial_charge = 0.,
                 esp_charge = 0.,
                 resp_charge = 0.,
                 index = 0,
                 tag = 0,
                 pdb_type = "",
                 amber_type = "",
                 residue_name = "",
                 residue_number = 0,
                 layer = "H"
            ):
        
        self.element = element
        self.position = position
        self.partial_charge = partial_charge
        self.esp_charge = esp_charge
        self.resp_charge = resp_charge
        self.index = index
        self.tag = tag
        self.pdb_type = pdb_type
        self.amber_type = amber_type
        self.residue_name = residue_name
        self.residue_number = residue_number
        self.layer = layer
        
    @property
    def element(self):
        return self._element
    @element.setter
    def element(self, value):
        if not type(value) in [str, np.str_]:
            raise ValueError("'element' must be a string, not {}".format(type(value)))
        self._element = str(value)
        
    @property
    def position(self):
        return self._position
    @position.setter
    def position(self, value):
        if not hasattr(value, "__len__"):
            raise ValueError("'position' must have a length of 3")
        if len(value) != 3:
            raise ValueError("'position' must have a length of 3")
        self._position = np.array(value, dtype = "float")
        
    @property
    def partial_charge(self):
        return self._partial_charge
    @partial_charge.setter
    def partial_charge(self, value):
        self._partial_charge = float(value)
        
    @property
    def esp_charge(self):
        return self._esp_charge
    @esp_charge.setter
    def esp_charge(self, value):
        self._esp_charge = float(value)
        
    @property
    def resp_charge(self):
        return self._resp_charge
    @resp_charge.setter
    def resp_charge(self, value):
        self._resp_charge = float(value)
        
    @property
    def index(self):
        return self._index
    @index.setter
    def index(self, value):
        if not type(value) in [int, np.int64]:
            raise ValueError("'index' must be a integer, not {}".format(type(value)))
        self._index = int(value)
        
    @property
    def tag(self):
        return self._tag
    @tag.setter
    def tag(self, value):
        if not type(value) in [int, np.int64]:
            raise ValueError("'tag' must be a integer, not {}".format(type(value)))
        self._tag = int(value)
        
    @property
    def pdb_type(self):
        return self._pdb_type
    @pdb_type.setter
    def pdb_type(self, value):
        if not type(value) in [str, np.str_]:
            raise ValueError("'pdb_type' must be a string, not {}".format(type(value)))
        self._pdb_type = str(value)
        
    @property
    def amber_type(self):
        return self._amber_type
    @amber_type.setter
    def amber_type(self, value):
        if not type(value) in [str, np.str_]:
            raise ValueError("'amber_type' must be a string, not {}".format(type(value)))
        self._amber_type = str(value)
        
    @property
    def residue_name(self):
        return self._residue_name
    @residue_name.setter
    def residue_name(self, value):
        if not type(value) in [str, np.str_]:
            raise ValueError("'residue_name' must be a string, not {}".format(type(value)))
        self._residue_name = str(value)
        
    @property
    def residue_number(self):
        return self._residue_number
    @residue_number.setter
    def residue_number(self, value):
        if not type(value) in [int, np.int64]:
            raise ValueError("'residue_number' must be a integer, not {}".format(type(value)))
        self._residue_number = int(value)
        
    @property
    def layer(self):
        return self._layer
    @layer.setter
    def layer(self, value):
        if not value in self._allowed_layer_names:
            raise ValueError("'layer' must be one of: {}, not {}".format(self._allowed_layer_names, value))
        self._layer = str(value)
        
    def __repr__(self):
        strings = ["Atom("]
        names = [
            "element",
            "position",
            "partial_charge",
            "esp_charge",
            "resp_charge",
            "index",
            "tag",
            "pdb_type",
            "amber_type",
            "residue_name",
            "residue_number",
            "layer"
        ]
        
        for name in names:
            strings.append("{}={},".format(name, getattr(self, name)))
        
        strings.append(")")
            
        return " ".join(strings)
        
class ComReader(object):
    def __init__(self, filename):
        
        self.filename = filename
        self.chk = ""
        self.oldchk = ""
        self.memstr = ""
        self.procstr = ""
        self.keywords = ""
        self.title = ""
        self.charges = []
        self.multiplicities = []
        self.additional = []

        self.atoms = None

        self.parse()

    def parse_atom_string(self, atom_string):
        symbol = ""
        amber = ""
        pdb = ""
        charge = 0.
        charge_str = ""
        residue = ""
        resnum = 0
        resnum_str = ""
        
        phase = 0

        
        for s in atom_string:
            if phase == 0:
                if s == "-":
                    phase = 1
                else:
                    symbol += s
            elif phase == 1:
                if s == "-":
                    phase = 2
                else:
                    amber += s
            elif phase == 2:
                if s == "(":
                    phase = 3
                    charge = float(charge_str)
                else:
                    charge_str += s
            elif phase == 3:
                if s == "=":
                    phase = 4
            elif phase == 4:
                if s == ",":
                    phase = 5
                else:
                    pdb += s
            elif phase == 5:
                if s == "=":
                    phase = 6
            elif phase == 6:
                if s == ",":
                    phase = 7
                else:
                    residue += s
            elif phase == 7:
                if s == "=":
                    phase = 8
            elif phase == 8:
                if s == ")":
                    phase = 9
                    resnum = int(resnum_str)
                else:
                    resnum_str += s

        return symbol, amber, charge, pdb, residue, resnum

    def parse(self):

        with open(self.filename, "r") as com_obj:
            com_lines = com_obj.readlines()


        elements = []
        ambers = []
        pdbs = []
        amber_charges = []
        residues = []
        resnums = []
        positions = []
        layers = []
        link_dict = dict()

        phase = 0
        atom_num = 0
        for l in com_lines:
            l = l.strip()
            if phase == 0:
                if l.startswith("%mem"):
                    self.memstr = l.split("=")[1]
                elif l.startswith(("%nproc", "%ncpus")):
                    self.procstr = l.split("=")[1]
                elif l.startswith("%chk"):
                    self.chk = l.split("=")[1]
                elif l.startswith("%oldchk"):
                    self.oldchk = l.split("=")[1]
                elif l.startswith("#"):
                    self.keywords = l
                    phase += 1
            elif phase == 1:
                if l:
                    self.keywords += " " + l
                else:
                    phase += 1
            elif phase == 2:
                if l:
                    self.title += " " + l
                else:
                    if not self.title:
                        raise RuntimeError("No title section in com")
                    phase += 1
            elif phase == 3:
                sl = l.split()
                num_layers = int(len(sl) / 2)
                self.charges = [int(sl[i * 2    ]) for i in range(num_layers)]
                self.mults   = [int(sl[i * 2 + 1]) for i in range(num_layers)]
                phase += 1
            elif phase == 4:
                if l:

                    sl = l.strip().split()
                    atom_string = sl[0]

                    element, amber, charge, pdb, residue, resnum = self.parse_atom_string(atom_string)

                    try:
                        int(sl[1])
                        o = 1
                    except:
                        o = 0

                    position = [float(p) for p in sl[1+o:4+o]]

                    elements.append(element)
                    ambers.append(amber)
                    pdbs.append(pdb)
                    amber_charges.append(charge)
                    residues.append(residue)
                    resnums.append(resnum)
                    positions.append(position)
                    
                    if atom_num not in link_dict:
                        link_dict[atom_num] = []

                    try:
                        layer = sl[4+o]
                        layers.append(layer)
                    except:
                        layers.append("H")
                        
                    try:
                        link = int(sl[6+o]) - 1
                        link_dict[atom_num].append(link)
                        
                        if not link in link_dict:
                            link_dict[link] = []
                        link_dict[link].append(atom_num)
                        
                    except:
                        pass

                    atom_num += 1
                else:
                    if not elements:
                        raise RuntimeError("No atoms section in com")
                    phase += 1
            elif phase == 5:
                self.additional.append(l)

        links = [link_dict[i] if i in link_dict else [] for i in range(len(positions))]
        
        self.atoms = Atoms(
            elements = elements,
            amber_types = ambers,
            pdb_types = pdbs,
            partial_charges = amber_charges,
            esp_charges = amber_charges,
            residue_names = residues,
            residue_numbers = resnums,
            positions = positions,
            layers = layers,
            links = links
        )

        self.get_neighbours()

    def get_neighbours(self):

        size = len(self.atoms)
        neighbours = [[] for i in range(size)]

        blanks = 0
        for i in range(len(self.additional)):
            a = self.additional[i].strip()
            sa = a.split()
    
            if blanks == 2:
                break
            elif a == "":
                blanks += 1
            elif sa[0].isdigit():
                index = int(sa[0]) - 1
                connections = sa[1:]

                ns = [int(connections[2*cn]) - 1 for cn in range(int(len(connections) / 2))]

                for n in ns:
                    if not n in neighbours[index]:
                        neighbours[index] += [n]
                    if not index in neighbours[n]:
                        neighbours[n] += [index]

        self.atoms.connectivity = neighbours

def read_com(fn):
    com_reader = ComReader(fn)
    return com_reader.atoms