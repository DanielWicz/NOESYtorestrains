import csv
import re
import os
from itertools import product, tee
from cardinality import count, at_least

class TblFileMaker:
    """A class that's pruporse is to create tbl restrain files from assigned
       NOESY signals. The assigment is read from a csv file with a format shown
       below:
           HydrogenPair,AssignedNOESYIntensity
           NH1 - NH2,m
           NH2 - NH5,s
           NH7 - NH1,w
       
       To use the class you've to create its object instance with 
       two arguments `noecsv_fname`, `sequence_fname` and `restwrite_fname` as in 
       the example below:
       
       ```python
       from NOESYtorestrains import TblFileMaker

       tblfm = TblFileMaker(noecsv_fname='noe.csv',
                            restwrite_fname='rest.tbl',
                            sequence_fname='sequence.seq')
       tblfm.read_noecsv()
       tblfm.save_tbl()
    ```
    The only files needed to exist are sequence file and csv file. The tbl file
    is provided by the program.
    Arguments:
        noecsv_fname: It is a filename for the csv file with assigned
                      NOESY signals. The CSV is arranged in the fallowing way:
                      `hydrogen, intensity
                       1HA - 2HN,s
                       2HA - 3HN,s
                       2HB - 3HN,w
                       3HA - 4HN,s
                       3HB - 4HN,m
                       3HG1 - 4HN,m
                       4HG12 - 5HN,w`
                       Where s, m, w are acronyms for strong, medium, weak intensities. 

        restwrite_fname: It is a filename for restrains file written by the class.
                         The file is going to be written in the fallowing format:
                         `assign (resid 2 and name HN)(resid 1 and name HA) 2.5 0.7 0.4 
                          assign (resid 3 and name HN)(resid 2 and name HA) 3.0 1.2 0.5 
                          assign (resid 4 and name HN)(resid 3 and name HA) 4.0 2.2 1.1`
                          for every restrain interaction based on the CSV and the sequence file.
        sequence_fname: It is a filename for file with three letter sequence (e.g. 
                        ALA ARG ASN ASP ...).
        one_atom_max_limit: How many combinations to calculate per atom. If you have ambigious restrains (restrains with respect
                            to hydrogen atoms, that can be assigned in various ways), then the number of different combinations
                            goes as factorial of n (n!).

    """

    def __init__(self, noecsv_fname=None, restwrite_fname=None, sequence_fname=None, one_atom_max_limit=1):
        # The values are d, d_minus, d_plus range: (d - d_minus, d + d_plus)
        self.basedict = {'w': [4.0, 2.2, 1.1], 'm': [3.0, 1.2, 0.5], 's': [2.5, 0.7, 0.4]}
        # Selecting force field types by hydrogen force field name has a small ambiguity
        # e.g. hydrogens HE1 and HE2 in hystidyne are for carbon and nitrogen aswell
        # It should be implemented with a dictionary of carbon types (CG: [hydrogen type ...] or other way.
        self.forcefield_types = {'ALA': ['HN', 'HA', 'HB1rot', 'HB2rot', 'HB3rot'],
                                 'ARG': ['HN', 'HA', 'HB1', 'HB2', 'HG1', 'HG2', 'HH11', 'HH12', 'HH21', 'HH22'],
                                 'ASN': ['HN', 'HA', 'HB1', 'HB2', 'HD21', 'HD22'],
                                 'ASP': ['HN', 'HA', 'HB1', 'HB2'],
                                 'CYS': ['HN', 'HA', 'HB1', 'HB2', 'HG'],
                                 'GLU': ['HN', 'HA', 'HB1', 'HB2', 'HG1', 'HG2'],
                                 'GLN': ['HN', 'HA', 'HB1', 'HB2', 'HG1', 'HG2', 'HH21', 'HH22'],
                                 'GLY': ['HN', 'HA1', 'HA2'],
                                 'HIS': ['HN', 'HA', 'HB1', 'HB2', 'HD1', 'HD2', 'HE1', 'HE2'],
                                 'ILE': ['HN', 'HA', 'HB', 'HG11', 'HG12', 'HG21rot', 'HG22rot', 'HG23rot', 'HD11rot', 'HD12rot', 'HD13rot'],
                                 'LEU': ['HN', 'HA', 'HB', 'HG11', 'HG12', 'HG21rot', 'HG22rot', 'HG23rot', 'HD11rot', 'HD12rot', 'HD13rot'],
                                 'LYS': ['HN', 'HA', 'HB1', 'HB2', 'HG1', 'HG2', 'HD1', 'HD2', 'HE1', 'HE2', 'HZ11rot', 'HZ12rot', 'HZ13rot'],
                                 'MET': ['HN', 'HA', 'HB1', 'HB2', 'HG1', 'HG2', 'HE1rot', 'HE2rot', 'HE3rot'],
                                 'PHE': ['HN', 'HA', 'HB1', 'HB2', 'HD1', 'HD2', 'HE1', 'HE2', 'HZ'],
                                 'PRO': ['HA', 'HB1', 'HB2', 'HG1', 'HG2', 'HD1', 'HD2'],
                                 'SER': ['HN', 'HA', 'HB1', 'HB2', 'HG'],
                                 'THR': ['HN', 'HA', 'HB', 'HG1', 'HG21rot', 'HG22rot', 'HG23rot'],
                                 'TRP': ['HN', 'HA', 'HB1', 'HB2', 'HD1', 'HE1', 'HE3', 'HZ2', 'HZ3', 'HH2'],
                                 'TYR': ['HN', 'HA', 'HB1', 'HB2', 'HD1', 'HD2', 'HE1', 'HE2', 'HH'],
                                 'VAL': ['HN', 'HA', 'HB', 'HG11rot', 'HG12rot', 'HG13rot', 'HG21rot', 'HG22rot', 'HG23rot'],
                                 'CPC': ['HN', 'HA', 'HB', 'HG11', 'HG12', 'HG21', 'HG22', 'HD1', 'HD2'],
                                 'TPC': ['HN', 'HA', 'HB', 'HG11', 'HG12', 'HG21', 'HG22', 'HD1', 'HD2'],
                                 'ACE': ['HA1rot', 'HA2rot', 'HA3rot']}

        self.basestring = 'assign (resid {0} and name {1})(resid {2} and name {3}) {4} {5} {6}'
        self.csv_regexp = '([0-9]+)([A-Z]+)([0-9]*)\ *[\–\-]\ *([0-9]+)([A-Z]+)([0-9]*)'
        self.regexp_comp = re.compile(self.csv_regexp)
        self.read_file_lines = []
        self.parameters_to_save = []
        self.aminoacid_sequence = []
        self.noecsv_fname = noecsv_fname
        self.restwrite_fname = restwrite_fname
        self.sequence_fname = sequence_fname
        self.one_atom_max_limit = one_atom_max_limit

    def read_sequence(self):
        with open(self.sequence_fname, 'r') as f:
            for line in f:
                for aa in line.strip().split(' '):
                    if any((aa == '', aa == ' ', aa is None)):
                        raise ValueError('The file does have a'+\
                                         'bad structure, check it')
                    self.aminoacid_sequence.append(aa)

    def read_noecsv(self):
        with open(self.noecsv_fname, 'r') as csvfile:
            reader = csv.DictReader(csvfile)
            for i, row in enumerate(reader):
                keys = list(row.keys())
                strength = row[keys[1]]
                d, d_min, d_plus = self.basedict[strength]
                matched = self.regexp_comp.match(row[keys[0]])
                # It is going as fallows:
                #    group(1) is name of the residue that belongs
                #    to firsts paired atom
                #    group(2) is first atom name that is paired with a
                #    second atom
                #    group(3) is name of the residue that belongs
                #    to second paired atom
                #    group(5) is second atom name that is paired
                #    with the first atom
                #    d, d_min, d_plus were described previously.
                self.parameters_to_save.append((int(matched.group(1)),
                                               matched.group(2),
                                               self.to_int_else_none(matched.group(3)),
                                               int(matched.group(4)),
                                               matched.group(5),
                                               self.to_int_else_none(matched.group(6)),
                                               float(d), float(d_min), float(d_plus)))

    def to_int_else_none(self, str_to_int):
        try:
            number = int(str_to_int)
        except:
            number = None

        return number

    def atom_to_forcefield(self, atom_name=None, aminoacid=None,
                           atom_num=None):
        """For given atom and residue name, returns all possible
           force field atom names from the topology file and
           considers rotation of the hydrogens around a bond,
           if all of them were marked as rot in the force field
           dictionary.
           Arguments:
               atom_name: atom name passed to the function from
                          the csv file/distance restrain entry list.
               aminoacid: One aminoacid that corresponds to 
                          the atom.
           Returns: List of atom names for corresponding aminoacid together
                    with rotatory mark, indicating if given hydrogens
                    are structurally equivalent..
        """
        ff_types_for_aa = self.forcefield_types[aminoacid]
        regexp = re.compile('([A-Z]+)([0-9]*)(rot)?')
        atom_names_list = []
        for ff_type_for_atom in ff_types_for_aa:
            matched = regexp.match(ff_type_for_atom)
            ff_atom_name = matched.group(1)
            ff_atom_num = self.to_int_else_none(matched.group(2))
            ff_atom_rot = matched.group(3)
            if ff_atom_rot == 'rot':
                rotational = True
            else:
                rotational = False
            if atom_num is not None and ff_atom_num is not None:
                if atom_num != ff_atom_num and len(str(atom_num)) == len(str(ff_atom_num)):
                    continue
                elif int(str(atom_num)[0]) != int(str(ff_atom_num)[0]):
                    continue
            if ff_atom_name == atom_name:
                if ff_atom_num is not None:
                    atom_names_list.append((ff_atom_name+str(ff_atom_num), rotational))
                else:
                    atom_names_list.append((ff_atom_name, rotational))

        return atom_names_list

    def distance_restrain_pairs(self, entry_pair=None):
        """Returns list of distance restrains for given distance restrain entry.
           
           Arguments:
               entry_pair: One distance restrain entry with its residue numbers
                          , distances, atom names. So one
                          entry of the self.parameters_to_save
                          field.
               aminoacid: One aminoacid that corresponds to the entry.
           Returns: List of distance restrains for the given distance entry.
        """
        res1, atom1, atomnum1, res2, atom2, atomnum2, d, d_min, d_max = entry_pair
        atom1_f = self.atom_to_forcefield(atom_name=atom1, atom_num=atomnum1,
                                          aminoacid=self.aminoacid_sequence[res1-1])
        atom2_f = self.atom_to_forcefield(atom_name=atom2, atom_num=atomnum2,
                                          aminoacid=self.aminoacid_sequence[res2-1])
        rotational_i_1, rotational_i_2 = 0, 0
        for atom1 in atom1_f:
            if atom1[1]:
                rotational_i_1 += 1
        for atom2 in atom2_f:
            if atom2[1]:
                rotational_i_2 += 1
        if rotational_i_1 == len(atom1_f):
            atom1_f = [atom1_f[0]]
        if rotational_i_2 == len(atom2_f):
            atom2_f = [atom2_f[0]]
        atom1_f_names_only = list(zip(*atom1_f))[0]
        atom2_f_names_only = list(zip(*atom2_f))[0]
        atom_product, atom_product_tee = tee(product(atom1_f_names_only, atom2_f_names_only))
        if at_least(self.one_atom_max_limit + 1, atom_product_tee): 
            print('possible pairs for aminoacids {0} and {1} exceeding {2}'.format(res1, res2, self.one_atom_max_limit))
            return None
        all_atom_pairs = list(atom_product)
        if all_atom_pairs == []:
            return []
        distance_restrains = []
        for pair in all_atom_pairs:
            distance_restrains.append((res1, pair[0], res2, pair[1], d, d_min, d_max))

        return distance_restrains

    def restrain_postprocess(self):
        """Post processes generated distance restrains, by
           replacing atom names with force field names and
           if there is several possibilities (due to equivalent
           hydrogens) for restrains, then generates few distance 
           restrains as a list.

           Arguments:
               None
           Returns: List of distance restrains.
        """
        
        repeated_entries = []
        for entry in self.parameters_to_save:
            distance_restrains_pairs = self.distance_restrain_pairs(entry_pair=entry)
            if distance_restrains_pairs == [] or distance_restrains_pairs is None:
                continue
            repeated_entries.append(distance_restrains_pairs)
        # Generate tuples of protein restrains
        protein_restrains, protein_restrains_tee = tee(product(*repeated_entries))
        print('Number of files to save is: {0}'.format(count(protein_restrains_tee)))
        
        return protein_restrains

    def save_tbl(self, tbl_dir='tblfiles'):
        try:
            os.mkdir(tbl_dir)
        except Exception as e:
            print('Directory with name {} already exists'.format(tbl_dir))
        else:
            print('Directory with name {} was just created'.format(tbl_dir))
        print('Generating different combinations for files to save, if it is too long, kill the app')
        protein_restrains = self.restrain_postprocess()
        for i, protein_restrain in enumerate(protein_restrains):
            to_save_string_list = []
            for x in protein_restrain:
                to_save_string_list.append(self.basestring.format(*x))
            with open('./'+tbl_dir+'/'+str(i)+'_'+self.restwrite_fname, 'w') as f:
                whole_tbl = '\n'.join(to_save_string_list)
                f.write(whole_tbl)

if __name__ == '__main__':
    tblmaker = TblFileMaker(noecsv_fname='noe.csv', restwrite_fname='protein.tbl',
                            sequence_fname='sequence.seq')
    tblmaker.read_noecsv()
    tblmaker.read_sequence()
    tblmaker.save_tbl()
