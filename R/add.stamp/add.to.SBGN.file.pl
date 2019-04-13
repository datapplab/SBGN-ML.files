use strict;
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
                my $content = `cat $sbgnFolder/$folder/$file`;
                $content =~s/\<\/sbgn\>/\n \<SBGNview.collection\>\<\/SBGNview.collection\> \n\<\/sbgn\>/g;
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
addStamp("MetaCrop");