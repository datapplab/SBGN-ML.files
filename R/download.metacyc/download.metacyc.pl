#use strict;
my @pathways = `cat /home/xiaoxi/0work/input/metacyc/metacyc.ids.txt`;
my $outputFolder = "/home/xiaoxi/0work/input/metacyc/biopax";
my $biopaxLevel = 2;
my $nPathways = scalar(@pathways);
my $i=0;
foreach my $pathway(@pathways){
	$i +=1;
	print "\n\n $i in $nPathways ",$pathway;
	$pathway =~s/\s+//g;
	next if (-e "$outputFolder/$pathway.biopax");
	my $jsTemplate = `cat /home/xiaoxi/0work/input/metacyc/download.metacyc.biopax.level3.js`;
	$jsTemplate =~s/pathwayCurrent/$pathway/g;
	$jsTemplate =~s/outFolder/$outputFolder/g;
	$jsTemplate =~s/biopax.level/$biopaxLevel/g;
	
	my $newJsFile = "$outputFolder.$pathway.js";
	open OUT,">$newJsFile";
	print OUT $jsTemplate;
	close OUT;
	print "downloading\n";
	`casperjs $newJsFile`;
	`rm $newJsFile `;
}
