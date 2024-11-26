import numpy as np
import gudhi as gd
import sympy as sp
from IPython.display import display, Math, Latex


# process a simplex tree into its entire rips filtration
def gudhi_rips_filtration(simplex_tree):
    """
    Returns a dictionary whose keys are a filtration value and the corresponding dictionary value is the simplicial complex in the rips filtration corresponding to this filtration value
    -----------------
    Args: 
    - simplex_tree: the simplex tree with filtration values to be interpreted
    
    Returns: 
    - rips_filtration: the dictionary described above
    """
    # get the generator of faces in the tree ordered by increasing filtration value
    rf = simplex_tree.get_filtration()

    # set a current_filtration value to compare future ones to
    current_filtration_value = 0.0
    # make an empty dictionary to be the rips filtration
    rips_filtration = {}

    # loop through the faces returned
    for term in rf:
        # unpack filtration filtration value
        filtration_value = term[1]
        # if the filtration value of the current simplex has changed
        if filtration_value != current_filtration_value:
            # make a copy of the original simplex tree
            # and prune the tree above the previous filtration value
            new_simplex = gd.SimplexTree(simplex_tree)
            new_simplex.prune_above_filtration(current_filtration_value)
            # add the new simplex to the list
            rips_filtration[current_filtration_value] = new_simplex

            # update the current filtration and add it to the list
            current_filtration_value = filtration_value

    # if you want to access parts of the filtration by index, you must cast the values (or whole dictionary) of the returned object to a lsit
    return rips_filtration


# process a simplex tree into part of its rips filtration
def gudhi_rips_filtration_index(simplex_tree, index, max_filtration_value=Ellipsis):
    """
    Returns a specific term in the rips filtration by its index (starting at 0).
    If max_filtration is specified, then this returns the simplicial complex which has this as its highest filtration value.
    -----------------
    Args: 
    - simplex_tree: the simplex tree with filtration values to be interpreted
    - index: the index in the rips filtration which is to be computed and returned
    - max_filtration_value: -1 by default to be ignored. if specified to be non negative, we will ignore index and instead return the term in the rips filtration which has all simplices with this filtration value or lower

    Returns: 
    - new_simplex: the simplex tree object based on the index or max_filtration_value as described above
    """
    # get the generator of faces in the tree ordered by increasing filtration value
    rf = simplex_tree.get_filtration()

    # check if the max filtration value is specified
    if max_filtration_value != Ellipsis:
        # create a copy of the given simplex and remove all simplices at higher filtration values than max_filtration_value
        new_simplex = gd.SimplexTree(simplex_tree)
        new_simplex.prune_above_filtration(max_filtration_value)
        return new_simplex

    # we must find the filtration value which falls at the desired index
    # initialize a list of filtraiton values
    current_filtration_value = 0.0
    filtration_values = [current_filtration_value]
    
    # add all distinct filtration values to the list
    for term in rf:
        if term[1] != current_filtration_value:
            filtration_values.append(term[1])

    # create a copy of the given simlex and remove all simplices at higher filtration values than the one at the desired index
    new_simplex = gd.SimplexTree(simplex_tree)
    new_simplex.prune_above_filtration(filtration_values[index])

    return new_simplex



def companion_array(simplicial_complex):
    """
    Converts a simplicial complex into a 2d python list whose ith row is a list of all i-simplices in the given simplicial complex
    -----------------
    Args: 
    - simplicial_complex: the simplicial complex which we want to convert to an array

    Returns: 
    - simplex_array: the 2d list whose ith row is a list of all i-simplices in the given simplicial_complex
    """
    
    # This function can accept different simplicial complex objects if you add the code in an if type(simplicial_complex) ==  statement below
    # print(type(simplicial_complex))
    # GHUDI SIMPLEX TREE
    if type(simplicial_complex) == gd.simplex_tree.SimplexTree:
        # get the list of simplices as a generator object
        gen = simplicial_complex.get_simplices()
        # initialize an empty 2d list with one row for each dimension of simplices in the complex
        simplex_array = [[] for i in range(1+simplicial_complex.dimension())]
        
        # loop through the terms in the generator object
        for term in gen:
            # unpack the term into a simplex and a filtration value
            simplex,f=term

            # add the simplex to the row of the simplex array corresponding to the dimension of the simplex
            simplex_array[len(simplex)-1].append(simplex)
    else:
        # change this to raise an error which tells the user that this version of the function does not support the given type of simplicial_complex but they can add this functionality themselves.
        raise ValueError(f"This function does not support the given simplicial complex data type. The current supported data types are: Gudhi simplex tree. You can add the functionality for the type {type(simplicial_complex)} by modifying this function.")
        
    return simplex_array


def boundary_matrix(complex_c, deg):
    """
    Given a simplicial complex companion array and a degree, computes the matrix associated to the boundary map in this degree of the simplicial chain complex using as a basis the simplices of the corresponding dimension in the array
    -----------------
    Args: 
    - complex_c: a 2d list whose ith row is a list of all i-simplices in the given simplex tree
    - deg: the degree (or dimension) which we want to compute the boundary map for (from degree to degree-1)

    Returns: 
    - matrix: a sympy matrix object which corresponds to the boundary map on the basis obtained from the given simplicial complex companion array
    """
    # set the basis to the set of simplices in the desired degree
    if deg >= len(complex_c):
        return 0
    else:
        basis = complex_c[deg]
    
    # if the degree is 0, the boundary is the 0 map
    # so return the matrix which sends each 0 simplex to 0
    if deg == 0:
        return sp.Matrix([[0]*len(basis)])
    
    # set the basis for the target space (deg-1) to the set of simplices of one dimension lower
    basis_next = complex_c[deg-1]
    # initialize a matrix with as many rows as basis elements in the target space and 0 columns
    matrix = sp.zeros(len(basis_next), 0)

    # loop through each simplex in the basis
    for simplex in basis:
        # initialize a new column with 0 entries
        column = [0]*len(basis_next)
        # for each facet of the simplex
        for i in range(deg+1):
            face = simplex.copy()
            face.pop(i)
            # find the index that this face appears in the basis for the target space
            index = basis_next.index(face)
            # at this index, add the appropriate coefficient of the face (a power of -1 to give us the alternating sum)
            column[index] = (-1)**i
        
        # join this column to the matrix
        matrix = matrix.row_join(sp.Matrix(column))
    
    return matrix




def homology(simplicial_complex, start, end = Ellipsis):
    """
    Computes the free rank of the homology with Z coefficients of a simplicial complex at every degree between two given indices (the ending index is set to the given starting index by default). Additionally, finda a list of free generators for the homology.
    --------------------
    Args:
    - simplicial_complex: the simplicial complex which we want to compute the homology of
    - start: the lowest degree to compute the homology in
    - end: the highest degree to compute the homology in
    Returns:
    - generators: a list of lists of generators
    - start: the degree at which the list of generators starts
    - complex_c: the companion array of the simplicial complex
    """
    
    # get the companion array for the simplex tree
    complex_c = companion_array(simplicial_complex)
    # dimension of the simplicial complex
    dimension = len(complex_c) - 1
    
    
    # if end is the placeholder object, Ellipsis, then set it to be the given starting index
    if end == Ellipsis:
        end = start
    # cap the last homology to compute at the top degree
    if end > dimension:
        end = dimension
    
    
    # 2d list whose kth row is a list of "vectors" (lists) in the (i+k) basis (the (i+k)th row of the companion array) which are generators for the free part of the (i+k)th homology group
    generators = []
    
    # first outgoing boundary
    boundary_out = boundary_matrix(complex_c, start)
    
    # we want to cut off the computation at dimension to handle this separately
    # if our ending index is equal to the dimension, we stop one loop early
    if end < dimension:
        cutoff = end+1
    else:
        cutoff = end
    
    # for each index from start to cutoff-1 (inclusive)
    for i in np.arange(start, cutoff):
        boundary_in = boundary_matrix(complex_c, i+1)
        # matrix which we will use to test if cycles are linearly independent, initialized to the incoming boundary
        linear_independence_test = sp.Matrix([list(b) for b in boundary_in.columnspace()]).T
        # create list of cycles
        cycles = boundary_out.nullspace()
        # initialize empty list for generators
        generating_cycles = []
        # check for each cycle...
        for cycle in cycles:
            # temporary copy of linear_independence_test which we add a cycle to
            test = linear_independence_test.col_insert(0, cycle)
            reduced,pivots = test.rref()
            # if the number of pivots columns is equal to the number of columns in the matrix
            if len(pivots)==sp.shape(test)[1]:
                # convert cycle to list
                cycle_as_list = list(cycle)
                
                # replace the linear_independence_test with the new set of linearly independent boundaries and cycles
                linear_independence_test = test

                generating_cycles.append(cycle_as_list)
        generators.append(generating_cycles)
        # update outgoing boundary for the next loop
        boundary_out = boundary_in
    
    if end == dimension:
        generators.append([list(c) for c in boundary_matrix(complex_c, dimension).nullspace()])
    
    
    # return the list of generators, the starting dimension for the homology, and the companion array for the given complex
    return generators,start,complex_c




def display_generators(generator_list, start, companion_array):
    """
    Interprets a list of generators starting at a certain index and a companion array and prints the generators in a readable format.
    ------------------
    Args:
    - generator_list: a list given by the homology function. Each element of the list is a list of vectors which correspond to generating cycles for the homology in a certain degree
    - start: the degree at which generator_list starts (this is returned by the homology function)
    - companion_array: the 2d python list whose ith row is a list of all i-simplices corresponding to the simplicial complex in question
    """
    # for each set of generators
    for i in range(len(generator_list)):
        # set the degree to the index of this set of generators plus the starting degree
        degree = i + start
        # if there are generators, print a header
        if len(generator_list[i]):
            # the size is the size of the basis of simplices for this degree
            size = len(generator_list[i][0])
            print(f"The free generators in degree {degree} are:")
        else:
            size = 0
            print(f"There are no free generators in degree {degree}.")
            print()
            # display(Math("0"))
        
        # for each generator in the set of generators
        for generator in generator_list[i]:
            # initialize an empty string which will become the formatted generator
            the_generator = ""
            # for each basis simplex
            for basis_index in range(size):
                # if the_generator is not empty (i.e. this is not the first thing being added to the string)
                if the_generator:
                    # check the coefficient
                    match generator[basis_index]:
                        # if the coefficient is 0, do nothing
                        case 0:
                            pass
                        # if the coefficient is 1, use a +
                        case 1:
                            the_generator += f"+{companion_array[degree][basis_index]}"
                        # if the coefficient is -1, use a -
                        case -1:
                            the_generator += f"-{companion_array[degree][basis_index]}"
                        # if the coefficient is negative, don't put a +
                        case _ if generator[basis_index] < 0:
                            the_generator += f"{generator[basis_index]}{companion_array[degree][basis_index]}"
                        # otherwise, just use + 
                        case _:
                            the_generator += f"+{generator[basis_index]}{companion_array[degree][basis_index]}"
                # if the generator is empty (i.e. this is the first thing being added to the list)
                else:
                    # check the coefficient
                    match generator[basis_index]:
                        # if the coefficient is 0, do nothing
                        case 0:
                            pass
                        # if the coefficient is 1, omit it
                        case 1:
                            the_generator += f"{companion_array[degree][basis_index]}"
                        # if the coefficient is -1, use a -
                        case -1:
                            the_generator += f"-{companion_array[degree][basis_index]}"
                        # otherwise, just put the coefficient with no sign
                        case _:
                            the_generator += f"{generator[basis_index]}{companion_array[degree][basis_index]}"
            
            # we have now populated the generator string with the linear combination of simplices
            # display the generator in latex
            display(Math(the_generator))
    return None