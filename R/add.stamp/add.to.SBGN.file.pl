use strict;
use List::Util qw[min max];
sub addStamp{
        my ($folder) = @_;
        # `mkdir $folder.without.stamp`;
        my $sbgnStampSolder = "C:/xiaoxi/work/papers/SBGNview/package/SBGN-ML.files/data/SBGN.with.stamp";
        my $sbgnFolder = "C:/xiaoxi/work/papers/SBGNview/package/SBGN-ML.files/data/SBGN";
        
        `mkdir $sbgnStampSolder/$folder`;
        print ($folder);
        my @files = `ls $sbgnFolder/$folder`;
        foreach my $file (@files){
                print($file);
                my $content = `cat $sbgnFolder/$folder/$file | grep 'bbox'`;
                #my @all_xs = $content =~ /x[ ]?=[ ]?(\d+\.\d+)[^\d]/g;
                my @all_xs = $content =~ /bbox.+x *= *"*(\d+\.?\d*)[^\d]/g;
                my $min_x = min(@all_xs) + 10;
                
                #my @all_ys = $content =~ /y[ ]?=(.+) /g;
                my @all_ys = $content =~ /bbox.+y *= *"*(\d+\.?\d*)[^\d]/g;
                #print(@all_xs[1..4]);
                my $min_y = min(@all_ys);
                my $max_y = max(@all_ys) + 70;
                print("\n\n>>>",$min_x,"<<<\n\n");
                print("\n\n>>>",$min_y,"<<<\n\n");
                print("\n\n>>>",$max_y,"<<<\n\n");
                
                
                my $stamp = '<glyph class="macromolecule" id="stamp" > <label text="SBGNhub Pathway Collection"/> <bbox w="180" h="20" x="20" y="'.$max_y.'"/> </glyph>';
                print "$stamp\n\n";
                
                $content =~s/\<\/map\>/\n $stamp \n\<\/map\>/g;
                # `cp $folder/$file $folder.without.stamp/$file`;
                open OUT ,">$sbgnStampSolder/$folder/$file";
                print OUT $content;
                close;
               # print "\n\n\n>>>,\n\n", $content," \n\n\n";
               # die;
        }
}
addStamp("pathwayCommons");
addStamp("MetaCyc");