#!/usr/bin/perl
use Getopt::Std;
use Math::Trig;
use Sys::Hostname;
$rc=180/pi;

getopts("de:hn:q:rs:vw");

# Print help message if -h option is specified
if($opt_h) {
    print "Usage: pov-movie.pl <switches> <snapshot-directory> <type> [<header-number>]\n\n";
    print "Switches:\n";
    print "-d             (Don't duplicate existing files)\n";
    print "-e <num>       (Only render every <num> frame)\n";
    print "-h             (Print this information)\n";
    print "-n <num>       (Render up to <num> threads)\n";
    print "-q <quality>   (Quality of rendering, 1=good, 3=extreme)\n";
    print "-r             (Render remotely in parallel)\n";
    print "-s <frame>     (Render a single frame)\n";
    print "-v             (Verbose output)\n";
    print "-w             (Disable making a movie)\n";
    exit 0;
}

die "Two or three arguments required" unless @ARGV==2 || @ARGV==3;

# Set variables used for remote processing
if($opt_r) {
    $rdir="esim/mesh/fiber";
    open A,"../config/rhosts" or die "Can't open remote hosts\n";
    @nlist=();@nthr=();
    while(<A>) {
        next if /^#/;
        @c=split;
        if($#c>=1) {
            push @nlist,$c[0];
            push @nthr,$c[1];
        }
    }
    close A;
    $nodes=$#nlist+1;
    $queue=$nodes==1?1:0;$h=0;
}

# Set miscellaneous variables, and those used to control which frames are
# rendered
$dr=$ARGV[0];
$verb=$opt_v?"":">/dev/null 2>/dev/null";
$every=$opt_e?$opt_e:1;
$a=defined $opt_s?$opt_s:0;
$opt_n=$opt_s if defined $opt_s;

# Read the first line of the POV-Ray header file to get the rendering
# dimensions. Assemble the POV-Ray flags.
$opt_q=1 unless defined $opt_q;
die "POV quality out of range\n" if $opt_q<0 || $opt_q>3;
open B,"pov_headers/h".(@ARGV==3?$ARGV[2]:"1").".pov"
    or die "Can't open POV-Ray header file\n";
$_=<B>;
m/^\/\/ W=(\d*) H=(\d*)/ or die "Can't read rendering size\n";
$pov_opts="+W$1 +H$2 ".
          ("","+R3 +A0.01 -J","+R6 +A0.001 -J","+R9 +A0.0001 -J")[$opt_q];

# Read, process, and store the rest of the POV-Ray header file
$gpn=0;
while(<B>) {
    $gp[$gpn++]=$_;
}
close B;$gpn--;

# Loop over the available frames
while(-e "$dr/$ARGV[1].$a") {

    # Assemble output filename and check for skipping/termination conditions
    $fn=sprintf "fr_%04d.png",$a;
    last if defined $opt_n && $a>$opt_n;
    $a++,next if defined $opt_d && -e "$dr/$fn";

    # Assemble the POV file
    $tf="rtemp$h.pov";
    open A,$opt_r?"| bzip2 -9 -c >$dr/$tf.bz2":">$dr/$tf" or die "Can't open temporary POV file\n";
    open B,"$dr/$ARGV[1].$a" or die "Can't open data file\n";
    ($t,$cx,$cy,$wall_sca)=split " ",<B>;
    foreach $i (0..$gpn) {
        $_=$gp[$i];
        s/DISHX/$cx/g;
        s/DISHY/$cy/g;
        s/WALL_SCA/$wall_sca/g;

        if(/^#include "sph\.pov"/) {
            while(<B>) {
                ($x,$y,$z,$q0,$q1,$q2,$q3)=split;
                print A "sphere{<$x,$y,$z>,0.5}\n" unless $x eq "nan";
            }
        } elsif(/^#include "ref.pov"/) {
            ($ef=$ARGV[1])=~s/^p/e/;
            if(open C,"$dr/$ef.$a") {
                <C>;
                while(<C>) {
                    ($x,$y,$z,$q0,$q1,$q2,$q3)=split;$z+=100;
                    print A "sphere{<$x,$y,$z>,0.05}\n" unless $x eq "nan";
                }
                close C;
            }
        } elsif(/^#include "rsph\.pov"/) {
            $q=$qq=0;
            while(<B>) {
                $q=1 if ++$q==45;
                $qq+=4;
                ($x,$y,$z,$q0,$q1,$q2,$q3)=split;
                $rx=$rc*atan2(2*($q0*$q1+$q2*$q3),1-2*($q1*$q1+$q2*$q2));
                $ry=$rc*asin(2*($q0*$q2-$q3*$q1));
                $rz=$rc*atan2(2*($q0*$q3+$q1*$q2),1-2*($q2*$q2+$q3*$q3));
                printf A "sphere{<0,0,0>,0.5 texture{T_Stone$q translate <$qq,0,0>} matrix<%.5g,%.5g,%.5g,%.5g,%.5g,%.5g,%.5g,%.5g,%.5g,$x,$y,$z>}\n",
                    1-2*($q2*$q2+$q3*$q3),2*($q1*$q2+$q3*$q0),2*($q1*$q3-$q2*$q0),
                    2*($q1*$q2-$q3*$q0),1-2*($q1*$q1+$q3*$q3),2*($q2*$q3+$q1*$q0),
                    2*($q1*$q3+$q2*$q0),2*($q2*$q3-$q1*$q0),1-2*($q1*$q1+$q2*$q2) ## conversion of quaternion to rotation matrix
                    unless $x eq "nan";
            }
        } else {
            print A;
        }
    }
    close A;
    close B;

    # Render the frame, forking jobs to remote processors if the "-r" option is
    # given
    $pov_cmd="nice -n 19 povray $tf -D +O$fn $pov_opts";
    if($opt_r) {

        # Send the POV-Ray file to a node for processing
        $hst=$nlist[$h];
        print "Frame $a to $hst\n";
        exec "rsync -q $dr/$tf.bz2 $hst:$rdir; ".
             "ssh $hst \"cd $rdir; bunzip2 -f $tf.bz2 ; $pov_cmd +WT$nthr[$h] \" $verb ; ".
             "rsync -q $hst:$rdir/$fn $dr ; ssh $hst \"rm $rdir/$fn $rdir/$tf\" " if ($pid[$h]=fork)==0;

        # Wait for one of the forked jobs to finish
        if ($queue) {
            $piddone=wait;$h=0;
            $h++ while $piddone!=$pid[$h]&&$h<$nodes;
            die "PID return error!\n" if $h>=$nodes;
        } else {
            $h++;$queue=1 if $h>=$nodes-1;
        }
    } else {

        # Run POV-Ray locally
        print "Frame $a\n";
        die if system "cd $dr; $pov_cmd $verb";
    }
    $a+=$every;
}

if($opt_r) {wait foreach 0..($queue?$nodes-1:$h);}

# Additional code to automatically make a movie
unless ($opt_w) {
    ($mf=$dr)=~s/\.odr//;
    unlink "${mf}_$ARGV[1].mov";
    system "ffmpeg -r 30 -y -i $dr/fr_%4d.png -preset fast -c:v libx265 -crf 17 -pix_fmt yuv420p -tag:v hvc1 -movflags faststart ${mf}_$ARGV[1].mov";
}
