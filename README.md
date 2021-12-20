# NOESYtorestrains


### Introduction
   A script that's pruporse is to create tbl restrain files from assigned
   NOESY signals. The assigment is read from a csv file with a format shown
   below:
       HydrogenPair,AssignedNOESYIntensity
       NH1 - NH2,m
       NH2 - NH5,s
       NH7 - NH1,w

   To use the class you've to create its object instance with
   two arguments `noecsv_fname`, `sequence_fname` and `restwrite_fname` as in
   the example below:

### Example

```
   from NOESYtorestrains import TblFileMaker

   tblfm = TblFileMaker(noecsv_fname='noe.csv',
                        restwrite_fname='rest.tbl',
                        sequence_fname='sequence.seq')
   tblfm.read_noecsv()
   tblfm.save_tbl()
```

### Usage 
The script is used as it is, imported directly from the directory, where the repository was cloned into.
