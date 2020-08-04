#!/usr/local/bin/php
<?php
ini_set('memory_limit', '3G');

$prg      = array_shift($argv);
$minmatch = array_shift($argv); /// min overlap NTs  
$mnfovl   = array_shift($argv); /// fraction of overlap of the smaller seq
$fst      = array_shift($argv);
$tmppsl   = array_shift($argv);

$seqs = readFasta($fst);

error(count($seqs)." seqs read from $fst",'i');

$xy=array();
$cp=0;

$discard = array();
$cselft = 0;
$cmn = 0;
$c = 0; 

$buf = fopen($tmppsl,"r");
if(!$buf) error("openning psl file ($tmppsl)",'e');

error("reading matches from $tmppsl ...",'i');

$l = fgets($buf);
while (!feof($buf)) {
  $l=trim($l);
  if($l!='') {
    $a = preg_split("/\t/",$l);
    
    $nq = $a[9];
    $nt = $a[13];
    if($nq == $nt) {
      $cselft++;
      }
    if($nq != $nt && !isset($discard[$nq]) ) {
      $lq = $a[10];
      $lt = $a[14];
      
      $m = $a[0]; /// match 
      $str = $a[8];
      
      if($m > $minmatch) {
		$v = "K"; /// keep
		if($lq>=$lt) {
	  		$fovl = $m/$lt;
	  		if($fovl>=$mnfovl) {
	    		$discard[$nt]=true; /// set target to discard
	    		$v="D"; /// Discard
	    		}
	  		}
		$cmn++;
		}
      $c++;
      }
    }
  $l = fgets($buf);
  }

fclose($buf);
error(" $cmn / $c with match >= $minmatch nt loaded (same seq = $cselft). ".count($discard)." seqs to discard",'i');

foreach($seqs as $k => $s) {
  if(!isset($discard[$k])) {
    print ">$k\n".chunk_split($s,100,"\n");
    }
  }

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

