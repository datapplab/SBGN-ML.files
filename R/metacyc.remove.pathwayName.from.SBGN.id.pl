use strict;
my $SbgnFileFolder = "C:/xiaoxi/work/papers/SBGNview/package/SBGN-ML.files/data/SBGN/MetaCyc";
my @sbgnFiles = `ls $SbgnFileFolder |grep -v 'svg'`;
my $outputFolder = "C:/xiaoxi/work/papers/SBGNview/package/SBGN-ML.files/data/SBGN/MetaCyc.correct.id";
`mkdir $outputFolder`;
foreach my $sbgnFile (@sbgnFiles){
    print $sbgnFile;
    chomp $sbgnFile;
    my $pathwayName = $sbgnFile;
    $pathwayName =~s/\.sbgn//g;
    $pathwayName =~s/.biopax.sbgn.sbgn//g;
    my $sbgnContent = `cat $SbgnFileFolder/$sbgnFile`;
    $sbgnContent =~s/http___http___BioCyc_org__META_pathway-biopax_type_3_38object_${pathwayName}//g;
    $sbgnContent =~s/(id|source|target)="(SmallMolecule|Protein)([^"]+)"/$1="$2$3:\@:$pathwayName"/g;
    #$sbgnContent =~s /_[^"]+"/"/g;
    print $1,"\n";
    open OUT, ">$outputFolder/$pathwayName.sbgn";
    print OUT $sbgnContent  ;
    close OUT;
}
