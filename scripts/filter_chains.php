#!/usr/local/bin/php
<?php
ini_set('memory_limit', '3G');

$prg   = array_shift($argv);
$mnlen = array_shift($argv);
$file  = array_shift($argv);

$uctgs = array();
$chns  = array();

$fi = file($file);
foreach($fi as $l) {
  $l=trim($l);
  if($l!='') {
    $a = preg_split("/\t/",$l); 
    if($a[0][0] == '#') {
      $ctg = $a[1];
      $ext = substr($ctg,-2); /// _R if reversed contig in chain
      if($ext=='_R') {	
		$ctg=substr($ctg,0,-2);
		}
      $uctgs[$ctg] = true;
      }
    elseif($a[0][0] == '>') {
      $chns[substr($a[0],1)] = array();
      }
    }
  }
  
$fst_chns = array_shift($argv);
$schns = readFasta($fst_chns);  

$fst_ctgs = array_shift($argv);
$sctgs = readFasta($fst_ctgs);  

foreach($schns as $k => $s) $schnslen[$k] = strlen($s); 
foreach($sctgs as $k => $s) $sctgslen[$k] = strlen($s); 

arsort($schnslen); 
arsort($sctgslen); 

$mnlenchns=0;
$mnlenctgs=0;

foreach($schnslen as $k => $l) { if($l>=$mnlen) $mnlenchns++; }
foreach($sctgslen as $k => $l) { if($l>=$mnlen) $mnlenctgs++; }

error(count($chns)." (=".count($schns).") chains, ".count($sctgs)." contigs",'i');
error(count($uctgs)." ctgs in chains out of ".count($sctgs),'i');
error("$mnlenctgs  contigs, $mnlenchns chains (>= $mnlen)",'i');

$coutchns = 0;
$fon = "$fst_chns.filter";
$fo = fopen($fon,'w');
foreach($schns as $k => $s) {
  if($schnslen[$k]>=$mnlen) {
    $coutchns++;
    fprintf($fo,">$k\n".chunk_split($s,100,"\n"));
    }
  }
fclose($fo);
error("$coutchns / ".count($schns)." chains seqs in $fon",'i');

$coutctgs = 0;
$fon = "$fst_ctgs.filter";
$fo = fopen($fon,'w');
foreach($sctgs as $k => $s) {
  if(!isset($uctgs[$k]) && $sctgslen[$k]>=$mnlen) {
    $coutctgs++;
    fprintf($fo,">$k\n".chunk_split($s,100,"\n"));
    }
  }
fclose($fo);
error("$coutctgs / ".count($sctgs)." contigs seqs in $fon",'i');


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

