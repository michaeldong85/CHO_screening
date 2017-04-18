#################################################################
#cho_screening.pl
#This script calculates whether there are stacking TRP residues #
#in a given pdb file. If there exists such kind of TRP residues,# 
#this protein is considered as carbohydrate binding protein     #
#################################################################

our $pdb_name=$ARGV[0];   #input the pdb file

##########
#pdb path#
##########

open FILE, "$pdb_name.pdb" or die $!;
my @lines=<FILE>;

#step1: simplify the pdb file by getting rid of ligand atoms and hydrogen atoms
our @pro; my $id=0;
foreach my $string (@lines){
	if (substr($string,0,4) eq "ATOM")
	{
		$pro[$id]=$string;
		$id=$id+1;
	}
}

#step2: get the aromatic residue list; sometimes the aromatic residues are mislabelled, and the verification step makes sure all the aromatic residues are correct
my @resi_list=&resi_list();
my @aromatic_resi_list=&aromatic_resi_list(@resi_list);
my @rings=&aromatic_residue_verification(@aromatic_resi_list);

#step3: defining the cylinder parameters for TRP/TYR/PHE
#TRP: height cutoff 6.7A; raidus cutoff 2.8A
#TYR: height cutoff 6.6A; raidus cutoff 2.8A
#PHE: height cutoff 6.7A; raidus cutoff 2.9A
my %cutoff_h; my %cutoff_r;
$cutoff_h{"TRP"}=6.7;       
$cutoff_h{"TYR"}=6.6;
$cutoff_h{"PHE"}=6.7;
$cutoff_r{"TRP"}=2.8;       
$cutoff_r{"TYR"}=2.8;
$cutoff_r{"PHE"}=2.9;

#step4: calculating how many atoms (there are two sides for each ring, and the one with smaller number of atoms is chosen) reside in the cylinder region of each aromatic residue
my @ring_position_info; my @ring_atom_info;
for (my $m=0; $m<scalar(@rings); ++$m)
{
	my ($ring_name, $x1, $y1, $z1, $x2, $y2, $z2, $x3, $y3, $z3, $xo, $yo, $zo)=&ring_information($rings[$m]);
	my ($a, $b, $c, $d)=&plane($x1, $y1, $z1, $x2, $y2, $z2, $xo, $yo, $zo);
	$ring_position_info[$m]="$xo	$yo	$zo	$a	$b	$c	$d";
	my @all_atom_position=&all_atom_position($x, $y, $z, $xo, $yo, $zo, $a, $b, $c, $d, $rings[$m]);
	our $height_cutoff=$cutoff_h{$ring_name}; our $radius_cutoff=$cutoff_r{$ring_name};
	my ($min, $judge)=&weak_side(@all_atom_position);
	$ring_atom_info[$m]="$min	$judge";
}


#step5: getting the list of aromatic residues with 0 atoms in their cylinder space
my @stack_resi; my $id1; 
my @stack_trp_resi; my $id2;
for (my $m=0; $m<scalar(@rings); ++$m)
{
	my @elements=split("	", $ring_atom_info[$m]);
	if ($elements[0]==0)
	{
		$stack_resi[$id1]="$rings[$m]	$ring_position_info[$m]	$ring_atom_info[$m]";
		$id1=$id1+1;
	}
	if (($elements[0]==0) and (substr($rings[$m],0,3) eq "TRP"))
	{
		$stack_trp_resi[$id2]="$rings[$m]	$ring_position_info[$m]	$ring_atom_info[$m]";
		$id2=$id2+1;
	}
}

#step6: in this step, the TRP residues satisfying two conditions are selected: 
#(1) there are 0 atoms in the cylinder space of TRP residues
#(2) there are at least one other aromatic residues (TRP or TYR or PHE with 0 atoms in their specific cylinder space) within 13.9A distance
for (my $m=0; $m<scalar(@stack_trp_resi); ++$m)
{
	my ($name1, $xo1, $yo1, $zo1, $a1, $b1, $c1, $d1, $min1, $judge1)=split("	", $stack_trp_resi[$m]);
	my @hits; my $id=0;
	for (my $n=0; $n<scalar(@stack_resi); ++$n)
	{
		my ($name2, $xo2, $yo2, $zo2, $a2, $b2, $c2, $d2, $min2, $judge2)=split("	", $stack_resi[$n]);
		if ($name1 ne $name2)
		{
			my $distance=sqrt(($xo1-$xo2)**2+($yo1-$yo2)**2+($zo1-$zo2)**2);
			if ($distance<13.9)
			{
				$hits[$id]=$stack_resi[$id]; $id=$id+1;
			}
		}
	}
	if (scalar(@hits)>0)
	{
		print "$name1";
		print "\n";
	}
}







#for each aromatic ring, there are two sides, and the weak side is difined as the side with less atoms in the cylinder space
sub weak_side{
	my @data=@_;
	my $positive=0; my $negative=0; my $judge; my $min;
	for (my $m=0; $m<scalar(@data); ++$m)
	{
		my @elements=split("	", $data[$m]);
		if (($elements[0]>0) and ($elements[0]<$radius_cutoff) and ($elements[1]>0) and ($elements[1]<$height_cutoff))
		{
			$positive=$positive+1;
		}
		if (($elements[0]>0) and ($elements[0]<$radius_cutoff) and ($elements[1]>-1*$height_cutoff) and ($elements[1]<0))	
		{
			$negative=$negative+1;
		}	
	}
	if ($positive>$negative)
	{
		$min=$negative; $judge=-1;
	}
	else
	{
		$min=$positive; $judge=1;
	}
	return ($min, $judge);
}

#this sub calculates the relative radius and height of all the atoms in a pdb file to a specific aromatic ring
sub all_atom_position{
	my ($x, $y, $z, $xo, $yo, $zo, $a, $b, $c, $d, $string)=@_;
	my @all_atom_position; my $id=0;
	for (my $n=0; $n<scalar(@pro); ++$n)
	{
		if ((substr($pro[$n],0,4) eq "ATOM") and (substr($pro[$n],17,9) ne $string))
		{
			my $stringxxx=substr($pro[$n],6,20);
			my $x=substr($pro[$n],30,8)+0; my $y=substr($pro[$n],38,8)+0; my $z=substr($pro[$n],46,8)+0;
			my ($radius, $height)=&radius_height($x, $y, $z, $xo, $yo, $zo, $a, $b, $c, $d);
			$all_atom_position[$id]="$radius	$height";
			$id=$id+1;
		}
	}	
	return @all_atom_position;
}

#this sub calculates for a given atom (x,y,z), the relative radius and height to a specific aromatic ring
sub radius_height{
	my ($x, $y, $z, $xo, $yo, $zo, $a, $b, $c, $d)=@_;
	my $height=($a*$x+$b*$y+$c*$z+$d)/sqrt($a**2+$b**2+$c**2);
	my $distance=sqrt(($x-$xo)**2+($y-$yo)**2+($z-$zo)**2);
	my $value=abs($distance**2-(abs($height))**2);
	my $radius=sqrt($value);
	return ($radius, $height);
}

#this sub is to get the ring plane information (ax+by+cz+d=0)
sub plane{
	my ($x1, $y1, $z1, $x2, $y2, $z2, $x3, $y3, $z3)=@_;
	my $x0=($x1+$x2+$x3)/3; my $y0=($y1+$y2+$y3)/3; my $z0=($z1+$z2+$z3)/3;
	my $a=($y2-$y1)*($z3-$z1)-($y3-$y1)*($z2-$z1);
	my $b=($z2-$z1)*($x3-$x1)-($z3-$z1)*($x2-$x1);
	my $c=($x2-$x1)*($y3-$y1)-($x3-$x1)*($y2-$y1);
	my $d=-$a*$x0-$b*$y0-$c*$z0;
	return ($a, $b, $c, $d);
}
#this sub is to get the ring information including ring name , 3 atoms defining the ring plane (x1,y1,z1); (x2,y2,z2); (x3,y3,z3); and the centor of the ring (x0,y0,z0)
sub ring_information{
	my @array=@_; my $name=substr($array[0],0,3);
	my ($ring_name, $x1, $y1, $z1, $x2, $y2, $z2, $x3, $y3, $z3, $xo, $yo, $zo);
	$ring_name=substr($array[0],0,3);
	if ((substr($array[0],0,3) eq "HIS") or (substr($array[0],0,3) eq "PHE") or (substr($array[0],0,3) eq "TYR"))
	{
		open FILE, "pai_ring/1-$name.txt" or die $!;
		my @lines=<FILE>; chomp(@lines);
		($xo, $yo, $zo); my @x; my @y; my @z;
		for (my $m=0; $m<scalar(@lines); ++$m)
		{
			for (my $n=0; $n<scalar(@pro); ++$n)
			{
				if ((substr($pro[$n],13,3) eq substr($lines[$m],0,3)) and (substr($pro[$n],17,9) eq $array[0]))
				{
					$x[$m]=substr($pro[$n],30,8); $y[$m]=substr($pro[$n],38,8); $z[$m]=substr($pro[$n],46,8);
				}
			}
		}
		my $sum_x; my $sum_y; my $sum_z;
		for (my $p=0; $p<scalar(@x); ++$p)
		{
			$sum_x=$x[$p]+$sum_x; $sum_y=$sum_y+$y[$p]; $sum_z=$sum_z+$z[$p];
		}
		$xo=$sum_x/scalar(@x); $yo=$sum_y/scalar(@y); $zo=$sum_z/scalar(@z);
		($x1, $y1, $z1, $x2, $y2, $z2, $x3, $y3, $z3)=($x[0], $y[0], $z[0], $x[2], $y[2], $z[2], $x[4], $y[4], $z[4]);
	}
	else
	{
		open FILE, "pai_ring/1-TRP.txt" or die $!;
		my @lines1=<FILE>; chomp(@lines1);
		open FILE, "pai_ring/2-TRP.txt" or die $!;
		my @lines2=<FILE>; chomp(@lines2);
		my ($xo1, $yo1, $zo1); my @x1; my @y1; my @z1;
		for (my $m=0; $m<scalar(@lines1); ++$m)
		{
			for (my $n=0; $n<scalar(@pro); ++$n)
			{
				if ((substr($pro[$n],13,3) eq substr($lines1[$m],0,3)) and (substr($pro[$n],17,9) eq $array[0]))
				{
					$x1[$m]=substr($pro[$n],30,8); $y1[$m]=substr($pro[$n],38,8); $z1[$m]=substr($pro[$n],46,8);
				}
			}
		}		
                my $sum_x1; my $sum_y1; my $sum_z1;
                for (my $p=0; $p<scalar(@x1); ++$p)
                {
                        $sum_x1=$x1[$p]+$sum_x1; $sum_y1=$sum_y1+$y1[$p]; $sum_z1=$sum_z1+$z1[$p];
                }
                $xo1=$sum_x1/scalar(@x1); $yo1=$sum_y1/scalar(@y1); $zo1=$sum_z1/scalar(@z1);

		my ($xo2, $yo2, $zo2); my @x2; my @y2; my @z2;
		for (my $m=0; $m<scalar(@lines2); ++$m)
		{
			for (my $n=0; $n<scalar(@pro); ++$n)
			{
				if ((substr($pro[$n],13,3) eq substr($lines2[$m],0,3)) and (substr($pro[$n],17,9) eq $array[0]))
				{
					$x2[$m]=substr($pro[$n],30,8); $y2[$m]=substr($pro[$n],38,8); $z2[$m]=substr($pro[$n],46,8);
				}
			}
		}		
                my $sum_x2; my $sum_y2; my $sum_z2;
                for (my $p=0; $p<scalar(@x2); ++$p)
                {
                        $sum_x2=$x2[$p]+$sum_x2; $sum_y2=$sum_y2+$y2[$p]; $sum_z2=$sum_z2+$z2[$p];
                }
                $xo2=$sum_x2/scalar(@x2); $yo2=$sum_y2/scalar(@y2); $zo2=$sum_z2/scalar(@z2);
		($x1, $y1, $z1, $x2, $y2, $z2, $x3, $y3, $z3)=($x2[0], $y2[0], $z2[0], $x2[2], $y2[2], $z2[2], $x2[4], $y2[4], $z2[4]);		
		$xo=($xo1+$xo2)/2; $yo=($yo1+$yo2)/2; $zo=($zo1+$zo2)/2;
	}
	return ($ring_name, $x1, $y1, $z1, $x2, $y2, $z2, $x3, $y3, $z3, $xo, $yo, $zo);
}

#this sub is to verify the aromatic residues in pdb files in case of mislabeling
sub aromatic_residue_verification{
	my @list=@_;
	my @rings; my $id=0;
	for (my $m=0; $m<scalar(@list); ++$m)
	{
		my $ring_name=substr($list[$m],0,3);
		open FILE, "pai_ring/$ring_name.txt" or die $!;
		my @lines=<FILE>;
		my $judge=0;
		for (my $n=0; $n<scalar(@lines); ++$n)
		{
			for (my $p=0; $p<scalar(@pro); ++$p)
			{
				if ((substr($pro[$p],13,3) eq substr($lines[$n],0,3)) and (substr($pro[$p],17,9) eq $list[$m]))
				{
					$judge=$judge+1;
					last;
				}
			}
		}
		if ($judge==scalar(@lines))
		{
			$rings[$id]=$list[$m];
			$id=$id+1;
		}
	}
	return @rings;
}
#this sub is to get the list of the aromatic residues in pdb file
sub aromatic_resi_list{
	my @resi_list=@_;
	my @aromatic_resi_list; my $id=0;
	for (my $m=0; $m<scalar(@resi_list); ++$m)
	{
		if ((substr($resi_list[$m],0,3) eq "PHE") or (substr($resi_list[$m],0,3) eq "TRP") or (substr($resi_list[$m],0,3) eq "TYR"))
		{
			$aromatic_resi_list[$id]=$resi_list[$m];
			$id=$id+1;
		}
	}
	return (@aromatic_resi_list);
}

# this sub is to get the list of residues in pdb file
sub resi_list{
	my @resi_name; my $id=0;
	for (my $m=0; $m<scalar(@pro); ++$m)
	{
		if (substr($pro[$m],0,4) eq "ATOM")
		{
			$resi_name[$id]=substr($pro[$m],17,9);
			$id=$id+1;
		}
	}
	my @uniq=&uniq(@resi_name);
	return @uniq;
}


sub uniq {
    my %seen;
    grep !$seen{$_}++, @_;
}

