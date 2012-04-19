#!/usr/bin/perl -w


# IMPORTANT: Set max read length $maxlength
# IMPORTANT: Set base for qual scores $qbase (33 or 64)

$maxlength=76;

# Specify quality scores used
$qbase = 33;   # Base for transforming to ascii
$qmax = 40;    # Max quality score

# Genome base composition
@Q = ( 0.25, 0.25, 0.25, 0.25 );

# Damage
# Set the probability of a C->T and G->A at any position in the
# 5'end of the read up to position $Ndamage.
# Cytosine to thymine misincorporation rates are highest (30.7%) at
# the first position of the sequences and decrease by approximately
# twofold per position as the read progresses (Fig. 4, bottom).
# This rate was reduced to 3.2% at the fifth nucleotide.
$Ndamage = 6;
@CT5 = (0.307,0.16,0.067,0.043,0.032,0.024);
@GA5 = (0.,0.,0.,0.,0.,0.);

# Similar for the 3' end of the read (starting with last base of read)
@CT3 = (0.,0.,0.,0.,0.,0.);
@GA3 = (0.307,0.16,0.067,0.043,0.032,0.024);


# Mutations
# If we have an overall mutation frequency A and assume that
# transversions and transitions have same prob,
# P(a|g)=A/3 for a different from g
$ptransition = 0.005/3;
$ptransversion = $ptransition;

# Alphabet
@alph = ( "A", "C", "G", "T");
# for ($i=0;$i<4;++$i) { $letternum{$alph[$i]} = $i; }



$logConstant = 1.0/log(2.0);


# Print 4x4 matrix
sub print_mat {
    my $mat=$_[0];
    for ($i=0; $i<4; ++$i) {
	print "$alph[$i]";
	for ($j=0; $j<4; ++$j) {
	    print " $mat->[$i][$j]";
	}
	print "\n";
    }
}

# This is the prob P(a|g) for sample base a and genome g.
sub mutation_matrix {
    # Initialize as transversions
    for ($i=0; $i<4; ++$i) {
	for ($j=0; $j<4; ++$j) { $mut[$i][$j] = $ptransversion; }
    }
    # Transitions
    $mut[0][2]=$mut[2][0]=$mut[1][3]=$mut[3][1]=$ptransition;
    # Diagonal
    $mut[0][0]=$mut[1][1]=$mut[2][2]=$mut[3][3]
	=1.-2.*$ptransversion-$ptransition;
    return @mut;
}


# This is the product matrix P(b|g) = \sum_a P(b|a)P(a|g)
# for damaged b, sample a and genome g.
# For each position in the 5' and 3' end it differs
#
# For a given set of parameters, it returns a reference to a matrix
sub damage_matrix {
    my $ct=$_[0];
    my $ga=$_[1];
    my @mut = @{$_[2]};
    my @dam = ();

    # First calculate the damage P(b|a) = $dam[$b][$a]

    # Initialize as diagonal
    for ($b=0; $b<4; ++$b) {
	for ($a=0; $a<4; ++$a) { $dam[$b][$a] = 0; }
	$dam[$b][$b] = 1.;
    }

    # C->T
    $dam[3][1]=$ct;     # P(T|C)
    $dam[1][1]=1.-$ct;  # P(C|C)
    # G->A
    $dam[0][2]=$ga;     # P(A|G)
    $dam[2][2]=1.-$ga;  # P(G|G)

#    print "Damage matrix\n";
#    print_mat(\@dam);

    # Product \sum_a P(b|a)P(a|g) where P(a|g) is $mut[$a][$g]
    my $r = [([(0,0,0,0)],[(0,0,0,0)],[(0,0,0,0)],[(0,0,0,0)])];
    for ($b=0; $b<4; ++$b) {
	for ($g=0; $g<4; ++$g) {
	    for ($a=0; $a<4; ++$a) {
		$r->[$b][$g] += $dam[$b][$a]*$mut[$a][$g];
	    }
	}
    }
#    print "Product matrix\n";
#    print_mat($r);
    return $r;
}



# Correction to qual score probability
#    A(b,g) = P(b,g)/\sum_g' P(b,g')
# where P(b,g)=P(b|g)P(g)
sub correction_matrix {
    my $r = $_[0];   # P(b|g)
    my $q = $_[1];   # P(g)

    my $A = [([(0,0,0,0)],[(0,0,0,0)],[(0,0,0,0)],[(0,0,0,0)])];

    for ($b=0; $b<4; ++$b) {
	$sum=0;
	for ($g=0; $g<4; ++$g) { $sum += $r->[$b][$g] * $q->[$g]; }
	for ($g=0; $g<4; ++$g) { $A->[$b][$g] = $r->[$b][$g] * $q->[$g]/$sum; }
    }
#    print "Correction matrix\n";
#    print_mat($A);
    return $A;
}


# The actual genomic probabilities for x - the called base and qual
# Returns the 4 base probs P(g|x)
sub genome_probs {
    my $A = $_[0];    # correction matrix A(b,g)
    my $p = $_[1];    # P(b|x)  4 x 1

    my $r = [(0,0,0,0)];

    for ($g=0; $g<4; ++$g) {
	for ($b=0; $b<4; ++$b) { $r->[$g] += $p->[$b] * $A->[$b][$g]; }
    }
    return $r;
}



# For a given qual score (number) and called base (number),
# return the scores of the four bases using the correction matrix $A
sub adjusted_scores {
    my $A = $_[0];    # correction matrix A(b,g)
    my $base = $_[1]; # called base
    my $qual = $_[2]; # quality score (interger between 0 and qmax)
    my $q = $_[3];    # P(g)

    my $e = 10.**(-0.1*$qual);
    $e=0.75 if ($e>0.75);
    my @p = ($e/3,$e/3,$e/3,$e/3);
    $p[$base] = 1.-$e;

    my $r = genome_probs($A,\@p);

    # Calculate log-scores
    for ($g=0; $g<4; ++$g) { $r->[$g] = $logConstant*log($r->[$g]/$q->[$g]); }

    return $r;
}







# Make a huge hash table of matrix scores for all possible
# combinations of base and qual score
# This needs to be done separately for each position in the
# 5' and 3' end

# Hash keys:
#   "CZ1" means base C quality Z, position 1 (5')
#   "Gu-2" means base G quality u, position -2 (3')
#   "A[" means base A quality [, at a position more than $Ndamage from both ends
# Position 0 is the first base of the read, ad position -0 (!) is the last


@MUT = mutation_matrix();
# First away from ends
# No damage, so R=MUT
$A = correction_matrix(\@MUT,\@Q);
for ($i=0; $i<=$qmax; ++$i) {
    for ($x=0; $x<4; ++$x) {
	$key = $alph[$x] . chr($qbase+$i);
	$scoreHash{$key} = adjusted_scores($A,$x,$i,\@Q);
#	for ($g=0; $g<4; ++$g) { print $scoreHash{$key}->[$g] ."  "; }
#	print "$key\n";
    }
}


# Do the same for 5' positions
for ($k=0; $k<$Ndamage; ++$k) {
    $R = damage_matrix($CT5[$k],$GA5[$k],\@MUT);
    $A = correction_matrix($R,\@Q);
    for ($i=0; $i<=$qmax; ++$i) {
	for ($x=0; $x<4; ++$x) {
	    $key = $alph[$x] . chr($qbase+$i) . $k;
	    $scoreHash{$key} = adjusted_scores($A,$x,$i,\@Q);
#	    for ($g=0; $g<4; ++$g) { print $scoreHash{$key}->[$g] ."  "; }
#	    print "$key\n";
	}
    }
}


# Do the same for 3' positions
for ($k=0; $k<$Ndamage; ++$k) {
    $R = damage_matrix($CT3[$k],$GA3[$k],\@MUT);
    $A = correction_matrix($R,\@Q);
    for ($i=0; $i<=$qmax; ++$i) {
	for ($x=0; $x<4; ++$x) {
	    $key = $alph[$x] . chr($qbase+$i) . "-" . $k;
	    $scoreHash{$key} = adjusted_scores($A,$x,$i,\@Q);
#	    for ($g=0; $g<4; ++$g) { print $scoreHash{$key}->[$g] ."  "; }
#	    print "$key\n";
	}
    }
}


# Now read in fastq file from stdin and convert to PSSMs
while (<>) {
    ++$k;
    $k=0 if (/^@/);
    print if ($k<2);   # Print ID and called sequence
    if ($k==1) {
	chop;
	@seq = split(//,uc($_));
    }
    if ($k==3) {
	print "&\n";
	chop;
	@qual = split(//,$_);
	$l = @qual;
	next if ($l<2*$Ndamage+2);
	@oline = ("","","","");
	for ($i=0;$i<$l;++$i) {
	    if ( $seq[$i] !~ /[ACGT]/) { $seq[$i]="A"; $qual[$i]=chr($qbase); }
	    $key = $seq[$i] . $qual[$i];
	    if ($i<$Ndamage) { $key .= $i; }
	    # Assume NO DAMAGE in full length reads
	    if ( $l<$maxlength ) {
		$j = $l-$i-1;
		if ($j<$Ndamage) { $key .= "-".$j; }
	    }
#	    print STDERR "$key $i $l\n";
	    $r = $scoreHash{$key};
	    for ($g=0; $g<4; ++$g) { $oline[$g] .= sprintf("%.2f ", $r->[$g]); }
	}
	for ($g=0; $g<4; ++$g) {
	    $oline[$g] =~s/ $/\n/;
	    print $oline[$g];
	}
    }
}






exit(0);

print "Mutation matrix P(a|g)\n";
print_mat(\@MUT);

$R = damage_matrix($CT5[0],$GA5[0],\@MUT);
print "product matrix P(b|g)\n";
print_mat($R);

$A = correction_matrix($R,\@Q);
print "Correction matrix A(b,g)\n";
print_mat($A);
