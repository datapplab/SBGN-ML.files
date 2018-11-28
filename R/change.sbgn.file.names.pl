use strict;
sub changeMetacycFileName{
        my $MetacycFolder = "C:/xiaoxi/work/papers/SBGNview/package/SBGN-ML.files/data/SBGN/MetaCyc/";
        my @MetaCycSbgnFiles = `ls $MetacycFolder `;
        foreach my $metacycFile (@MetaCycSbgnFiles){
                chomp $metacycFile;
                my $oldFile = $metacycFile;
                $metacycFile =~s/\.biopax\.sbgn//g;
                print "\n\n $MetacycFolder/$metacycFile\n\n";
                print "\n\n $MetacycFolder/$oldFile\n\n";
                my $cpCommand = "mv -T $MetacycFolder/$oldFile $MetacycFolder/$metacycFile";
                print "\n\n$cpCommand\n\n";
                `$cpCommand`;
        }
}

#changeMetacycFileName();

sub changePathwayCommonsFileName{
        #my $PathwayCommonsFolder = "C:/xiaoxi/work/papers/SBGNview/package/SBGN-ML.files/data/SBGN/pathwayCommons/";
        my $PathwayCommonsFolder = "C:/xiaoxi/work/output/layout.sbgn.svg.pathway.commons.pathways-sbgn.new.sbgn.with.svg";
        my @MetaCycSbgnFiles = `ls $PathwayCommonsFolder `;
        foreach my $metacycFile (@MetaCycSbgnFiles){
                chomp $metacycFile;
                my $oldFile = $metacycFile;
                $metacycFile =~s/\.xml\.sbgn//g;
                print "\n\n $PathwayCommonsFolder/$metacycFile\n\n";
                print "\n\n $PathwayCommonsFolder/$oldFile\n\n";
                my $cpCommand = "mv -T $PathwayCommonsFolder/$oldFile $PathwayCommonsFolder/$metacycFile";
                print "\n\n$cpCommand\n\n";
                `$cpCommand`;
        }
}

changePathwayCommonsFileName();