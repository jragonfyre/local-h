#! /usr/bin/env polymake --script
#
# Short description for local-h.pl
#
# Author jragonfyre <jragonfyre@jragonfyre>
# Version 0.1
# Copyright (C) 2017 jragonfyre <jragonfyre@jragonfyre>
# Modified On 2017-07-05 09:17
# Created  2017-07-05 09:17
#
use strict;
use warnings;
use application "polytope"

sub standard_simplex_vertices {
    my $n = shift;
    my (@vert_arr, @temp_arr);
    for my $i (0..$n) {
        $temp_arr[$i]=0;
    }
    for my $i (0..$n) {
        my @copy = @temp_arr;
        $copy[$i] = 1;
        push @vert_arr, \@copy;
    }
    return \@vert_arr;
}

sub pArr {
    my $aref = shift;
    print "[",join(", ",@{$aref}),"]\n";
}

sub pArrArr {
    my $aref = shift;
    print "[\n";
    for my $arrRef (@$aref) {
        print "[",join(", ",@{$arrRef}),"]\n";
    }
    print "]\n";
}

sub standard_simplex {
    my $n = shift;
    return new Polytope(POINTS => standard_simplex_vertices($n));
}

sub restriction {
    my $subdiv = shift;
    my $n = $subdiv->DIM;
    my @face = @{shift @_}; # face is an arrref of vertices on face.
                            # vertices are given by their index from 0 to $n 
    #pArr \@face;
    my @faceComplement = grep { my $v = $_; not(grep { $v == $_ } @face); } (0..$n);
    #pArr \@faceComplement;
    my $onface = sub {
        # takes an arrRef and returns true if the point with those coords is on the face given by $face.
        my $aref = shift;
        for my $i (@faceComplement) {
            if ($aref->[$i] != 0) {
                return 0;
            }
        }
        return 1;
    };
    my $indices = $subdiv->VERTEX_INDICES;
    #pArr($indices);
    my $coords = $subdiv->COORDINATES;
    #pArrArr($coords);
    my $taggedcoords = []; #[map { [$_,$coords->[$_]] } @{$indices}];
    for my $i (0..$#{ $indices }) {
        $taggedcoords->[$i]=[$i,$coords->[$indices->[$i]]];
    }
    # find the vertices on the face
    #pArrArr($taggedcoords);
    my $tCoordsOnFace = [ grep { &{$onface}(@{$_}[1]); } @{$taggedcoords} ];
    my $coordsOnFace = [ map { @{$_}[0]; } @{$tCoordsOnFace} ];
    #pArr($coordsOnFace);
    my $ans = Polymake::topaz::induced_subcomplex($subdiv, $coordsOnFace); # forgets the geometry
    return $ans;
}

sub subarrs {
    if (scalar(@_) == 0) {
        return ([]);
    }
    my $first = shift;
    my @rec = subarrs(@_);
    my @recp = map { my @temp = @{$_}; unshift @temp, $first; \@temp; } @rec;
    push(@rec, @recp);
    return @rec;
}

sub powMinusOne {
    my $n = shift;
    if($n % 2) {
        return -1;
    }
    return 1;
}

sub addToIndex {
    my $aref = shift;
    my $i = shift;
    my $val = shift;
    if (not defined $aref->[$i]) {
        $aref->[$i]=0;
    }
    $aref->[$i] = $aref->[$i] + $val;
}

sub getNRationalPoints {
    my $n = shift;
    my $dim = shift; # dimension of ambient space not simplex.
    if ($n < 0) {
        return [];
    }
    if ($dim==0) {
        return ($n==0)?[[]]:[];
    }
    if ($dim==1) {
        return [[$n]];
    }
    my @pts = ();
    for my $i (0..$n) {
        for my $pt (@{getNRationalPoints($n-$i, $dim-1)}) {
            push(@{$pt}, $i);
            push(@pts, $pt);
        }
    }
    return \@pts;
}

sub getBoundaryNRationalPointsRec {
    my $n = shift;
    my $dim = shift; # dimension of ambient space not simplex.
    if ($n < 0) {
        return [];
    }
    if ($dim==0) {
        return ($n==0)?[[]]:[];
    }
    if ($dim==1) {
        return ($n==0)?[[0]]:[];
    }
    my @pts = ();
    for my $pt (@{getNRationalPoints($n, $dim-1)}) {
        push(@{$pt},0);
        push(@pts, $pt);
    }
    for my $i (1..$n) {
        for my $pt (@{getBoundaryNRationalPointsRec($n-$i, $dim-1)}) {
            push(@{$pt}, $i);
            push(@pts, $pt);
        }
    }
    return \@pts;
}


sub getBoundaryNRationalPoints {
	my $n = shift;
    my $dim = shift;
    return new Matrix<Rational>(getBoundaryNRationalPointsRec($n,$dim));
}

sub getRandomVector {
    my $len = shift;
    my @ans = ();
    for my $i (1..$len) {
        push(@ans, rand());
    }
    return \@ans;
}

sub toBarycentricCoords {
    my @arr = @{shift @_};
    my $sum = 0;
    for my $val (@arr) {
        $sum += $val;
    }
    return [ map { $_ / $sum } @arr ]
}

sub regularSubdivisionOfStandardSimplex { 
    my @coords = @{shift @_};
    my $weightsRef = new Vector<Rational>(shift @_);
    my @baryCoords = map { toBarycentricCoords($_); } @coords;
    my $baryMat=new Matrix<Rational>(\@baryCoords);
    return new topaz::GeometricSimplicialComplex(COORDINATES=>$baryMat,
        INPUT_FACES=>regular_subdivision($baryMat, $weightsRef));
}

sub randomRegularSubdivOfStandardSimplex {
    my $coordsRef = shift;
    my $weightsRef = getRandomVector(scalar(@{$coordsRef}));
    return regularSubdivisionOfStandardSimplex($coordsRef, $weightsRef);
}

sub totallyRandomRegularSubdivOfStandardSimplex {
    my $dSimplex = shift;
    my $n = shift;
    my $coordsRef = getNRationalPoints($n, ($dSimplex + 1));
    my $weightsRef = getRandomVector(scalar(@{$coordsRef}));
    return regularSubdivisionOfStandardSimplex($coordsRef, $weightsRef);
}


# computes as alternating sum of H_VECTORs
sub local_h {
    my $subdiv = shift; # subdivision of a standard simplex of the appropriate dimension
    my $n = $subdiv->DIM;
    my $d = $n+1;
    my $ans = [];
    for my $face (subarrs(0..$n)) {
        if (scalar(@{$face}) == 0) {
            addToIndex($ans,0,powMinusOne($d));
        } else {
            my $rest = restriction($subdiv,$face);
            my @hvect = @{$rest->H_VECTOR};
            #pArr(\@hvect);
            for my $i (0..$#hvect) {
                addToIndex($ans, $i, powMinusOne($n-$rest->DIM)*$hvect[$i]);
            }
        }
    }
    return $ans;
}

sub generateABunch {
	my $dim = shift;
	my $n = shift;
	my $iter = shift;
	my @ans = ();
	for my $i (1..$iter) {
		my $rsub = randomRegularSubdivOfStandardSimplex(getBoundaryNRationalPoints($n,$dim+1));
		my $lh = local_h($rsub);
		my $sum=0;
		for my $j (2..$dim-1) {
			$sum+= $lh->[$j];
		}
		if($sum == 0) {
			push(@ans, $rsub);
		} 
		print "on iteration $i, found ", scalar(@ans), "\n";
	}
	return \@ans;
}

sub isCone {
	#print "isCone start";
    my $subdiv = shift;
    my $n = $subdiv->DIM;
    my $indices = $subdiv->VERTEX_INDICES;
    my $nverts = $subdiv->N_VERTICES;
    my $coords = $subdiv->COORDINATES;
    my $taggedcoords = []; #[map { [$_,$coords->[$_]] } @{$indices}];
    for my $i (0..$#{ $indices }) {
        $taggedcoords->[$i]=[$i,$coords->[$indices->[$i]]];
    }
    for my $i (0..$n){
	    my $onface = sub {
	        # takes an arrRef and returns true if the point with those coords is on the face given by $face.
	        my $aref = shift;
	        if ($aref->[$i] != 0) {
	            return 0;
	        }
	        return 1;
	    };
	    my $tCoordsOnFace = [ grep { &{$onface}(@{$_}[1]); } @{$taggedcoords} ];
	    if (scalar(@{$tCoordsOnFace})==$nverts-1){
	    	#print "isCone end";
	    	return 1;
	    }
    }
    #print "isCone end";
    return 0;
 }

sub generateGoodShit {
	my $dim = shift;
	my $n = shift;
	my $iter = shift;
	my @ans = ();
	for my $i (1..$iter) {
		my $rsub = randomRegularSubdivOfStandardSimplex(getBoundaryNRationalPoints($n,$dim+1));
		my $lh = local_h($rsub);
		my $sum=0;
		for my $j (2..$dim-1) {
			$sum+= $lh->[$j];
		}
		if($sum == 0 and not(isCone($rsub))) {
			push(@ans, $rsub);
		} 
		print "on iteration $i, found ", scalar(@ans), "\n";
		if (isCone($rsub)) {
			print " Look, sir -- CONES!\n";
		}
	}
	return \@ans;
}

sub printSubdiv {
	$subdiv=shift;
	print "F vector: ", $subdiv->F_VECTOR,"\n";
	print "H vector: ", $subdiv->H_VECTOR,"\n";
	print "local h vector: ";
	pArr(local_h($subdiv));
	print "\n";
	print "Facets: \n", $subdiv->FACETS,"\n";
	print "Vertex coords:\n";
	for my $i (0..$#{$subdiv->VERTEX_INDICES}) {
		print "$i: ";
		pArr($subdiv->COORDINATES->[$subdiv->VERTEX_INDICES->[$i]]);
	}
}

@letters = ("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z");

sub indexToPointName {
	use integer;
	my $n = shift;
	my $quot=$n/26;
	my $rem=$n % 26;
	return ($letters[$rem] . ("'" x $quot));
}

sub coordsToGeogebraString {
	my @coords = @{shift @_};
	my $scale = shift;
	my $ans = "(";
	shift @coords;
	my $lcrd = pop @coords;
	for my $crd (@coords) {
		$ans = $ans . ($crd * $scale) . ",";
	}
	return ($ans . ($lcrd * $scale) . ")");
}

sub printToGeogebra {
	my $subdiv = shift;
	my $scale =shift;
	my @pts = map { "\"" . indexToPointName($_) . "=" . coordsToGeogebraString($subdiv->COORDINATES->[$subdiv->VERTEX_INDICES->[$_]],$scale) . "\"" } (0..$subdiv->N_VERTICES-1);
	my @facets = map { "\"Pyramid[" . join(",", map { indexToPointName($_); } @{$_}) . "]\"" } @{$subdiv->FACETS};
	print "Execute[{", join(",",@pts),",",join(",",@facets),"}]";
}

sub getEdges {
	my @facet = @{shift @_};
	my @ans = ();
	for my $i (0..scalar(@facet)-1){
		for my $j ($i+1..scalar(@facet)-1){
			push(@ans,[$facet[$i],$facet[$j]]);
		}	
	}
	return \@ans;
}

sub 

sub getGraph {
	my $subdiv = shift;
	my $facets = $subdiv->FACETS;
	my @edgelists = map {
		getEdges($_);
	} @{$facets};
	my @edges = ();
	for my $ref (@edgelists) {
		for my $edge (@{$ref}) {
			push (@edges,$edge);
		}
	}
	@edges
	$subdiv->COORDINATES->[$subdiv->VERTEX_INDICES->[$_]]
}

