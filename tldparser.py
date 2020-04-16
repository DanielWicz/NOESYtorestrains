import csv
import re

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

    """
    def __init__(self, noecsv_fname=None, restwrite_fname=None):
        # The values are d, d_minus, d_plus range: (d - d_minus, d + d_plus)
        self.basedict = {'w': [4.0, 2.2, 1.1], 'm': [3.0, 1.2, 0.5], 's': [2.5, 0.7, 0.4]}
        self.basestring = 'assign (resid {0} and name {1})(resid {2} and name {3}) {4} {5} {6}'
        self.csv_regexp = '([A-Z]+)([0-9]+)\ *[\–\-]\ *([A-Z]+)([0-9]+)'
        self.regexp_comp = re.compile(self.csv_regexp)
        self.read_file_lines = [] 
        self.noecsv_fname = noecsv_fname
        self.restwrite_fname = restwrite_fname

    def read_noecsv(self):
        with open(self.noecsv_fname, 'r') as csvfile:
            reader = csv.DictReader(csvfile)
            for i, row in enumerate(reader):
                keys = list(row.keys())
                strength = row[keys[1]]
                d, d_min, d_plus = self.basedict[strength]
                matched = self.regexp_comp.match(row[keys[0]])
                line_to_save = self.basestring.format(matched.group(2),
                                                      matched.group(1),
                                                      matched.group(4),
                                                      matched.group(3),
                                                      d, d_min, d_plus)
                self.read_file_lines.append(line_to_save)


    def save_tbl(self):
        with open(self.restwrite_fname, 'w') as f:
            whole_tbl = '\n'.join(self.read_file_lines)
            f.write(whole_tbl)


if __name__ == '__main__':
    tblmaker = TblFileMaker(noecsv_fname='noe.csv', restwrite_fname='protein.tbl')
    tblmaker.read_noecsv()
    tblmaker.save_tbl()
