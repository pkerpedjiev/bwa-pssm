#!/usr/bin/env python
''' Example use : 
cat file.fq | ./fastq2wm33.py > ./bwa pssm ... - 
Significant speedup can be achieved with using the pypy interpreter: http://pypy.org/
After installation, change interpreter from python to pypy in first line or run
cat file.fq | pypy fastq2wm33.py > ./bwa pssm ... - '''

import array
import sys
import math

from optparse import OptionParser

# ------------------------------------------------------------
# Global parameters
# ------------------------------------------------------------

# IMPORTANT: Set max read length $maxlength
# IMPORTANT: Set base for qual scores $qbase (33 or 64)

maxlength=76

# Specify quality scores used
qbase     = 33         # Base for transforming to ascii
qbase_chr = chr(qbase) # Convert above to character repr.
qmax      = 40         # Max quality score

# Genome base composition
Q = array.array( 'f', [0.25, 0.25, 0.25, 0.25] )

# Damage
# Set the probability of a C->T and G->A at any position in the
# 5'end of the read up to position $Ndamage.
# Cytosine to thymine misincorporation rates are highest (30.7%) at
# the first position of the sequences and decrease by approximately
# twofold per position as the read progresses (Fig. 4, bottom).
# This rate was reduced to 3.2% at the fifth nucleotide.
Ndamage = 6

CT5 = array.array( 'f', [0.307,0.16,0.067,0.043,0.032,0.024] )
GA5 = array.array( 'f', [0.,0.,0.,0.,0.,0.] )

# Similar for the 3' end of the read (starting with last base of read)
CT3 = array.array( 'f', [0.,0.,0.,0.,0.,0.] )
GA3 = array.array( 'f', [0.307,0.16,0.067,0.043,0.032,0.024] )

# Mutations
# If we have an overall mutation frequency A and assume that
# transversions and transitions have same prob,
# P(a|g)=A/3 for a different from g

ptransition = 0.005/3
ptransversion = ptransition

# Alphabet
alph = array.array( 'c', ["A", "C", "G", "T"] )

logConstant = 1.0/math.log(2.0);


# ------------------------------------------------------------
# define matrix functions
# ------------------------------------------------------------


# ------------------------------------------------------------
# Define mutation matrix
# This is the prob P(a|g) for sample base a and genome g.
def mutation_matrix() : 

    #Create a matrix (list of lists) 
    mut = [ [ptransversion for i in xrange(4)] for j in range(4) ]

    # Fill transitions
    mut[0][2] = mut[2][0] = mut[1][3] = mut[3][1] = ptransition

    # Fill diagonals
    mut[0][0] = mut[1][1] = mut[2][2] = mut[3][3] = 1.0 - 2.0*ptransversion - ptransition

    return mut

# ------------------------------------------------------------
# Define the damage matrix
# This is the product matrix P(b|g) = \sum_a P(b|a)P(a|g)
# for damaged b, sample a and genome g.
# For each position in the 5' and 3' end it differs
#
# For a given set of parameters, it returns a reference to a matrix
def damage_matrix( ct, ga, mut ) : 

    # Construct as identity matrix (list of lists)
    dam = [ [1.0 if i == j else 0 for i in xrange(4)] for j in range(4) ]

    # C -> T
    dam[3][1] = ct        # P(T|C)
    dam[1][1] = 1.0 - ct  # P(C|C)

    # G -> A 
    dam[0][2] = ga        # P(A|G)
    dam[2][2] = 1.0 - ga  # P(G|G)


    # Product \sum_a P(b|a)P(a|g) where P(a|g) is $mut[$a][$g]
    r = [ [0.0 for i in xrange(4)] for j in range(4) ]

    for b in xrange( 4 ) : 
        for g in xrange( 4 ) : 
            for a in xrange( 4 ) :
                r[b][g] += dam[b][a] * mut[a][g]


    return r


# ------------------------------------------------------------
# Define correction matrix
# Correction to qual score probability
#    A(b,g) = P(b,g)/\sum_g' P(b,g')
# where P(b,g)=P(b|g)P(g)
def correction_matrix( r, q ) : # arguments ( P(b|g), P(g) )
    
    A = [ [0.0 for i in xrange(4)] for j in range(4) ]

    for b in xrange( 4 ) : 

        # Sum over r(b|g) * q(g) for g = 0,1,2,3
        s = sum( iter(r[b][g] * q[g] for g in xrange(4) ) )

        for g in xrange(4) : 
            A[b][g] = r[b][g] * q[g] / s
    
    return A

# ------------------------------------------------------------
# Define geno probabilities
# The actual genomic probabilities for x - the called base and qual
# Returns the 4 base probs P(g|x)
def genome_probs( A, p ) : 

    r = [0.0 for i in xrange(4)]

    for g in xrange(4) : 
        for b in xrange(4) : 
            r[g] += p[b] * A[b][g]

    return r

# ------------------------------------------------------------
# Calculate adjusted scores
# For a given qual score (number) and called base (number),
# return the scores of the four bases using the correction matrix $A
def adjusted_scores( A, base, qual, q ) : 

    e = math.pow(10, -0.1 * qual )
    if e > 0.75 : e = 0.75

    p = [e/3, e/3, e/3, e/3] 
    p[base] = 1.0 - e

    r = genome_probs(A, p)

    for g in xrange(4) : 
        r[g] = logConstant * math.log( r[g] / q[g] )


    # Convert output to str and return
    return [ '%.2f'%ir for ir in r ]



def main():

    usage = """
    Example use : 
    cat file.fq | ./fastq2wm33.py > ./bwa pssm ... - 
    Significant speedup can be achieved with using the pypy interpreter: http://pypy.org/
    After installation, change interpreter from python to pypy in first line or run

    cat file.fq | pypy fastq2wm33.py > ./bwa pssm ... -
    """
    num_args= 0
    parser = OptionParser(usage=usage)

    #parser.add_option('-o', '--options', dest='some_option', default='yo', help="Place holder for a real option", type='str')
    #parser.add_option('-u', '--useless', dest='uselesss', default=False, action='store_true', help='Another useless option')

    (options, args) = parser.parse_args()

    if len(args) < num_args:
        parser.print_help()
        sys.exit(1)

    # ScoreHash will be a dictionary of string lists : [key][base] = str(value)
    scoreHash = {}

    
    mut = mutation_matrix()
    A = correction_matrix( mut, Q ) 



    # scoreHash away from ends: No damage; Rate = mutation matrix
    for i in xrange( qmax ) : 
        for x in xrange( 4 ) : 
            key = alph[x] + chr( qbase+i )
            scoreHash[key] = adjusted_scores( A,x,i,Q )




    # Do the same for 5' positions
    for k in xrange( Ndamage ) : 
        # Matrices here
        R = damage_matrix( CT5[k], GA5[k], mut )
        A = correction_matrix( R, Q )

        for i in xrange( qmax ) : 
            for x in xrange( 4 ) : 

                key = alph[x] + chr( qbase+i ) + str(k)
                scoreHash[key] = adjusted_scores( A, x, i, Q )



    # Do the same for 3' positions
    for k in xrange( Ndamage ) : 
        # Matrices here
        R = damage_matrix( CT3[k], GA3[k], mut )
        A = correction_matrix( R,Q )

        for i in xrange( qmax ) : 
            for x in xrange( 4 ) : 

                key = alph[x] + chr( qbase+i ) + '-' + str(k)
                scoreHash[key] = adjusted_scores(A, x, i, Q)

    # ------------------------------------------------------------
    # Read from STDIN, write to STDOUT
    # ------------------------------------------------------------

    count = 0
    seq  = array.array('c')
    qual = array.array('c')

    for line in sys.stdin : 

        # New read
        if line[0] == '@' : 
            print_str = ''
            count = 0

        if count < 2 : 
            print line,

        if count == 1 : 
            seq = array.array('c', list( line.strip() ) )

        elif count == 3 : 
            print '&'

            qual = array.array('c', list( line.strip() ) )
            length = len( qual )
 
            if length < 2*Ndamage + 2 : continue

            score_arr = [[] for x in xrange(length)]

            for i, (iseq, iqual) in enumerate( zip( seq, qual ) ) : 

                # Not recognised base
                if iseq not in 'ACGT' : 
                    iseq = 'A'
                    iqual = qbase_chr

                key = iseq + iqual

                # If in beginning of read
                if i < Ndamage : key += str(i)

                # If in end of read
                elif length < maxlength : 
                    j = length-i-1
                    if j < Ndamage : 
                        key += '-'+str(j)


                score_arr[i] = scoreHash[key]

            # Transpose list-list
            score_arr = map(list, zip(*score_arr))

            print '\n'.join( [' '.join( irow ) for irow in score_arr ] )

        # Increment where in read line is
        count+=1


if __name__ == '__main__':
    sys.exit( main() )

