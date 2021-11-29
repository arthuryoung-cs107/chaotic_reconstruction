#!/usr/bin/perl
@ARGV==1 or die "Syntax: ./arg.pl <output_dir>\n";

# Loop over the available parameter state files, taking into account that they
# are only saved intermittently
$n=0;$k=0;
while ($k<101) {
    if(-e "$ARGV[0]/st.$n") {

        # Open the state file in binary mode
        open A,"$ARGV[0]/st.$n" or die "Can't open stats file\n";
        binmode A;

        # Compute accumulators for each parameters
        $l=0;@sc=();@scc=();
        while( 52==read A,$b,52) {
            @c=unpack "fffffffffffff",$b;
            foreach (0..11) {
                $sc[$_]+=$c[$_+1];
                $scc[$_]+=$c[$_+1]*$c[$_+1];
            }
            $l++;
        }
        close A;

        # Compute the mean and standard deviation for each parameter, and print
        # the result
        @mc=();
        @mcc=();
        foreach (0..11) {
            $mc[$_]=$sc[$_]/$l;
            $mcc[$_]=sqrt($scc[$_]/$l-$mc[$_]*$mc[$_]);
        }
        print "$n @mc @mcc $l\n";
        $k=0;
    } else {
        $k++;
    }
    $n++;
}
