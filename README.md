# Free Homology Generators

This is a python module which computes generaing cycles for simplicial homology.

### Compatibility

The motivation for this module is to not only compute the (free) generators for homology, but to adapt to different simplicial complex data types from different topological data analysis software. To do this, the program has functions to turn a simplicial complex into a "companion array" which is a 2d list whose ith row is the list of all i-simplices. To add compatibility with new data types, one should add the code to convert their desired simplicial complex into this array according to its type.

Currently, this module is only compatible with simplex trees in Ghudi.


### How to Use

There is a jupyter notebook file which walks through the usage of this module (using Ghudi).


### Future Considerations

- Add compatibility with other data types.
- Add functionality to detect generators for torsion.
