var casper = require('casper').create();
var pathway = 'pathwayCurrent';
var outputFolder = 'outFolder';
var biopaxLevel = 'biopax.level';
casper.start('https://metacyc.org/META/NEW-IMAGE?type=NIL&object='+pathway+'&redirect=T', function() {
      var url = 'https://metacyc.org/META/pathway-biopax?type='+biopaxLevel+'&object='+pathway;
      this.download(url, outputFolder+'/'+pathway+'.biopax');
  });
  casper.run();