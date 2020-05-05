import csv
import re
from itertools import product

class TblFileMaker:
    """A class that's pruporse is to create tbl restrain files from assigned
       NOESY signals. The assigment is saved into a csv file with a format shown
       below:
           HydrogenPair,AssignedNOESYIntensity
           NH1 - NH2,m
           NH2 - NH5,s
           NH7 - NH1,w
       To use the class you've to create its object instance with 
       two arguments `noecsv_fname` and `restwrite_fname` as in 
       the example below:
       
       ```python
       from NOESYtorestrains import TblFileMaker

       tblfm = TblFileMaker(noecsv_fname='noe.csv',
                            restwrite_fname='rest.tbl')
       tblfm.read_noecsv()
       tblfm.save_tbl()
    ```
    Arguments:
        noecsv_fname: It is a filename for the csv file with assigned
                      NOESY signals.
        restwrite_fname: It is a filename for restrains file written by the class.
        sequence_fname: It is a filename for file with three letter sequence (e.g. 
                        ALA ARG ASN ASP ...).

    """

    def __init__(self, noecsv_fname=None, restwrite_fname=None, sequence_fname=None):
        # The values are d, d_minus, d_plus range: (d - d_minus, d + d_plus)
        self.basedict = {'w': [4.0, 2.2, 1.1], 'm': [3.0, 1.2, 0.5], 's': [2.5, 0.7, 0.4]}
        self.forcefield_types = {'ALA': ['HN', 'HA', 'HB1rot', 'HB2rot', 'HB3rot'],
                                 'ARG': [],
                                 'ASN': [],
                                 'ASP': [],
                                 'CYS': [],
                                 'GLU': ['HN', 'HA', 'HB1', 'HB2', 'HG1', 'HG2'],
                                 'GLN': [],
                                 'GLY': ['HN', 'HA1', 'HA2'],
                                 'HIS': [],
                                 'ILE': ['HN', 'HA', 'HB', 'HG11', 'HG12', 'HG21rot', 'HG22rot', 'HG23rot', 'HD11rot', 'HD12rot', 'HD13rot'],
                                 'LEU': [],
                                 'LYS': ['HN', 'HA', 'HB1', 'HB2', 'HG1', 'HG2', 'HD1', 'HD2', 'HE1', 'HE2', 'HZ1rot', 'HZ2rot', 'HZ3rot'],
                                 'MET': [],
                                 'PHE': [],
                                 'PRO': [],
                                 'SER': [],
                                 'THR': [],
                                 'TRP': ['HN', 'HA', 'HB1', 'HB2', 'HD1', 'HE1', 'HE3', 'HZ2', 'HZ3', 'HH2'],
                                 'TYR': ['HN', 'HA', 'HB1', 'HB2', 'HD1', 'HD2', 'HE1', 'HE2', 'HH'],
                                 'VAL': [],
                                 'CPC': ['HN', 'HA', 'HB', 'HG11', 'HG12', 'HG21', 'HG22', 'HD1', 'HD2'],
                                 'ACE': ['HA1rot', 'HA2rot', 'HA3rot']}

        self.basestring = 'assign (resid {0} and name {1})(resid {2} and name {3}) {4} {5} {6}'
        self.csv_regexp = '([A-Z]+)([0-9]+)\ *[\–\-]\ *([A-Z]+)([0-9]+)'
        self.regexp_comp = re.compile(self.csv_regexp)
        self.read_file_lines = []
        self.parameters_to_save = []
        self.aminoacid_sequence = []
        self.noecsv_fname = noecsv_fname
        self.restwrite_fname = restwrite_fname
        self.sequence_fname = sequence_fname

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
                #    group(2) is name of the residue that belongs
                #    to firsts paired atom
                #    group(1) is first atom name that is paired with a
                #    second atom
                #    group(4) is name of the residue that belongs
                #    to second paired atom
                #    group(3) is second atom name that is paired
                #    with the first atom
                #    d, d_min, d_plus were described previously.
                self.parameters_to_save.append((int(matched.group(2)),
                                               matched.group(1),
                                               int(matched.group(4)),
                                               matched.group(3),
                                               float(d), float(d_min), float(d_plus)))

    def atom_to_forcefield(self, atom_name=None, aminoacid=None):
        """For given atom and residue name, returns all possible
           force field atom names from the topology file
           Arguments:
               atom_name: atom name passed to the function from
                          the csv file/distance restrain entry list.
               aminoacid: One aminoacid that corresponds to 
                          the atom.
           Returns: List of atom names for corresponding aminoacid.
        """
        pass 

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
        pass


           
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
        pass

    def save_tbl(self):
        with open(self.restwrite_fname, 'w') as f:
            whole_tbl = '\n'.join(self.read_file_lines)
            f.write(whole_tbl)

if __name__ == '__main__':
    tblmaker = TblFileMaker(noecsv_fname='noe.csv', restwrite_fname='protein.tbl',
                            sequence_fname='sequence.seq')
    tblmaker.read_noecsv()
    tblmaker.read_sequence()
