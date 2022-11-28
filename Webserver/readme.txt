The GReg algorithm expects 4 input:
One file path that contains one pure DNA sequence inside
Maximum loop length allowed for the G-quadruplex (suggested from 1 to 10)
Maximum bulge length allowed for the G-quadruplex (suggested from 0 to 4)
Minimum folding temperatures as a threshold (suggested from 0 to 89)

The output will show each of the G4 containing regions like the example below: sequence, multiplicity and more details

GGGCGCGCGGGGCGCGGGTGGGG
[5, 5, 5, 0, 0, 0, 0, 0, 2, 4, 5, 3, 0, 1, 0, 5, 5, 5, 0, 3, 5, 5, 2]
> position: 282-304 | number of nucleotides: 23 | percentage G content: 74 | total number of quadruplex(es) could form 5 | number of tandem repeat(s) 1 | maximum estimated melting temperature: 65.8 | median estimated melting temperature: 60.0 | minimum estimated melting temperature:53.7