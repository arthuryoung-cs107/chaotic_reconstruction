#!/usr/bin/perl

# Starting frame
$s=0;

# Total beads
$n=100;$j=0;$m=0;
@t=();

$minx=1e6;$maxx=0;
$miny=1e6;$maxy=0;

open F,"${n}_ramp.txt" or die "Can't open input file\n";
<F>;
while(<F>) {
    ($tt,undef,$id,$xx,$yy)=split;
    
    break if int($id)!=$j;

    $t[$m++]=$tt,$j=0 if ++$j==$n;

    $minx=$xx if $xx<$minx;
    $maxx=$xx if $xx>$maxx;
    $miny=$yy if $yy<$miny;
    $maxy=$yy if $yy>$maxy;
    push @x,$xx;
    push @y,$yy;
};

$r=$m;
while(<F>) {
    $r++,$j=0 if ++$j==$n;
}
close F;

print "# $m records out of $r processed\n";
die "# No records to save\n" if $m<=$s;
$m=49900;

printf "# $minx $maxx $miny $maxy %g %g\n",0.5*($minx+$maxx),0.5*($miny+$maxy);
open G,">${n}_ramp.dat";
binmode G;
print G pack("ii",$n,$m-$s);
print G pack("f",$t[$_]) foreach ($s..($m-1));
print G pack("ff",$x[$_],$y[$_]) foreach ($s*$n..($m*$n-1));

open F,"${n}_ramp.the" or die "Can't open angle file\n";
$k=0;
<F> foreach 0..4;
while($k<$m) {
    ($ti,$ang)=split " ",<F>;

    die "$k $ti $t[$k]\n" unless $ti==$t[$k];
    print G pack("f",$ang) if $k>=$s;
    $k++;
}
close F;
close G;

exit;

$h[$_]=0 foreach 0..200;
foreach $f ($s..($m-1)) {
    foreach $i (0..($n-1)) {
        $xx=$x[$f*$n+$i];
        $yy=$y[$f*$n+$i];
        foreach $j (0..($n-1)) {
            $dx=$xx-$x[$f*$n+$j];
            $dy=$yy-$y[$f*$n+$j];
            $de=sqrt $dx*$dx+$dy*$dy;
            $h[int(4*$de+0.5)]++ if $de<50;
        }
    }
}

printf "%g $h[$_]\n",$_*0.25 foreach 0..200;
