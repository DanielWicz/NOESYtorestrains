# NOESYtorestrains


### Introduction
   A script that's pruporse is to create tbl restrain files from assigned
   NOESY signals. The assigment is read from a csv file with a format shown
   below:
   
   ```
       HydrogenPair,AssignedNOESYIntensity
       NH1 - NH2,m
       NH2 - NH5,s
       NH7 - NH1,w
   ```
   
   To use the class you've to create its object instance with
   two arguments `noecsv_fname`, `sequence_fname` and `restwrite_fname` as in
   the example below.
   
   The script also allows for generation of `.tbl` restrain files that have ambigiously defined
   atoms (e.g. from signals that can be asssigned either between hydrogens `H1, H2` or `H1, H3`).
   To use this option, in the CSV file don't specify its number e.g. if you take the example at the beginning,
   then instead of writting `NH1 - NH2,m`, you write `NH - NH2,m` or `NH1 - NH,m`.
   To limit number of such combinations, specify `one_atom_max_limit=N`, where `N` is number of combinations
   for one atom - atom pair. E.g. if there would be 9 combinations for a pair and `N=3`, then the program
   would crash with appropiate message. Atoms with such a large combinations should be avoided, especially
   with other atoms with such a large number of combinations, because the total number of comb. grows as `N!` (or `N` factorial).  

### Example

```
   from noesytorestrains import TblFileMaker

   tblfm = TblFileMaker(noecsv_fname='noe.csv',
                        restwrite_fname='rest.tbl',
                        sequence_fname='sequence.seq',
                        one_atom_max_limit=1)
   tblfm.read_noecsv()
   tblfm.save_tbl()
```

### Usage 
The script is used as it is, imported directly from the directory, where the repository was cloned into.

### Requipments
- Python 3.5 or above
- Cardinality package (`https://pypi.org/project/cardinality/`)
