#!/usr/local/bin/php
<?php
ini_set('memory_limit', '3G');

$prg  = array_shift($argv);
$fst  = array_shift($argv);
$blat = array_shift($argv);

$self_blat = "$fst.self.psl";

system("$blat -maxIntron=0 -noHead -maxGap=0 $fst $fst $self_blat");

$seqs = readFasta($fst);

error(count($seqs)." seqs read from $fst",'i');

$xy=array();

$cp=0;

$buf = fopen($self_blat,"r");
$l = fgets($buf);
while (!feof($buf)) {
  $l=trim($l);
  if($l!='') {
    $a = preg_split("/\t/",$l); //print_r($a);
    if($a[0]==$a[10] && $a[0]==$a[14] && $a[8]=='-' && ($a[9]!=$a[13]) ) { 
      if(!isset($xy[$a[9]])) $xy[$a[9]]=array();
      $xy[$a[9]][$a[13]]=true;
      if(isset($keep[$a[13]]) && isset($rm[$a[9]])) {
		$cp++;
		fprintf(STDERR,"\r $cp complementary pair already set "); 
		}
      else {
		$keep[$a[9]]=true;
		if(isset($keep[$a[13]])) error($a[13]." set to keep",'e');
		else $rm[$a[13]]=true;
		}
      }
    }
  $l = fgets($buf);
  }
fclose($buf);

fprintf(STDERR,"\n"); 
error(count($rm)." to eliminate cos revcomp redundant",'i');

$coutctgs = 0;
$fon = "$fst.red1";
$fo = fopen($fon,'w');
foreach($seqs as $k => $s) {
  if(!isset($rm[$k])) {
    $coutctgs++;
    fprintf($fo,">$k\n".chunk_split($s,100,"\n"));
    }
  }
fclose($fo);
error("$coutctgs / ".count($seqs)." seqs in $fon",'i');

die;

////////////////////////////////////////////////////////////
function error($l,$d){
  switch($d) {
    case 'W':
    case 'w':
      fprintf(STDERR,"WARNING: $l\n");
      break;
    case 'E':
    case 'e':
      fprintf(STDERR,"ERROR: $l\n");
      die;
      break;
    case 'I':
    case 'i':
      fprintf(STDERR,"# $l\n");
      break;
    default:
      break;
    }
  }

////////////////////////////////////////////////////////////
function readFasta($fst){
  $iSeq = array();
  $buf = fopen($fst,"r");
  if(!$buf) error("openning fasta file ($fst)",'e');
  $ID=-1;
  $l = fgets($buf);
  while (!feof($buf)) {
    $l = trim($l);
    if($l != "") {
      if($l[0] != ">") {
		$iSeq[$ID] .= strtoupper($l); 
		}
      else { 
		$ID=substr($l,1);
		$iSeq[$ID] ="";
		}
      }
    $l = fgets($buf);
    }
  fclose($buf);
  return($iSeq);
  }

