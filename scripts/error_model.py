#!/usr/bin/python

from optparse import OptionParser

import numpy as np
import sys


usage = """
A script for creating an error model from biased sequencing data.

Information about the biases should be formatted as a 4x4 table
indicating the probability of an experimentally-derived change
from the rows to the columns (i.e. the cell in row i and column j
should indicate the probability that base i was changed to base j.

A more concrete example is provided on the bwa-pssm web page:

    http://bwa-pssm.binf.ku.dk

Example:

    1   0       0   0
    0   0.8     0   0
    0   0       1   0
    0   0.2     0   1

The output is a table where the first column indicates the quality
to be replaced, the second the identity of the base and the last
four are numbers corresponding to the columns of the PSSM to be built.

Example:

    ...
    39 A 2.00 -12.77 -12.54 -12.34
    39 C -12.54 1.77 -12.54 -0.74
    39 G -12.54 -12.77 2.00 -12.34
    39 T -12.54 -12.77 -12.54 2.00
    40 A 2.00 -13.10 -12.87 -12.67
    40 C -12.87 1.77 -12.87 -0.74
    40 G -12.87 -13.10 2.00 -12.67
    40 T -12.87 -13.10 -12.87 2.00
"""

def phred_prob(q):
    '''
    Return the error probability corresponding to a quality value of 'q'.
    
    @param q: The PHRED quality value.
    @return: The error probability.
    '''
    return 10 ** (-q / 10.)

def main():
    parser = OptionParser(usage=usage)

    (options, args) = parser.parse_args()
    if len(args) < 1:
        parser.print_help()
        sys.exit(1)

    in_filename = args[0]

    if in_filename == '-':
        in_file = sys.stdin
    else:
        in_file = open(in_filename, 'r')

    bias_matrix = np.genfromtxt(in_file)
    
    if bias_matrix.shape != (4,4):
        print >>sys.stderr, "Invalid dimensions for the input table:", t.shape
        sys.exit(1)

    colsums = np.sum(bias_matrix, axis=0)
    if not np.allclose(colsums, np.array([1., 1., 1., 1.])):
        print >>sys.stderr, """
The probabilities along the columns should sum to \
1. The probabilities in the supplied table have \
the following sums: %s
""" % (colsums)

    bases = ['A', 'C', 'G', 'T']
    qual_matrix = np.zeros((4,4))

    for q in xrange(0, 42):
        error_prob = phred_prob(q) / 3.   #there's three different bases it could error to
        match_prob = 1 - 3. * error_prob

        qual_matrix.fill(error_prob)
        np.fill_diagonal(qual_matrix, match_prob)

        # pseudocounts
        qual_matrix += 10 ** -8

        bg = 0.25
        new_probs = bias_matrix.dot(qual_matrix)
        pssm_scores = np.log(new_probs / bg) / np.log(2.)

        for i, b in enumerate(bases):
            print "%d %s %s" % (q, b, " ".join(map("{:.4f}".format, pssm_scores[:,i])))

    in_file.close()

if __name__ == "__main__":
    main()
