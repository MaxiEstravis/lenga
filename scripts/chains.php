#!/usr/local/bin/php
<?php
ini_set('memory_limit', '8G');

$prg  = array_shift($argv);
$fst  = array_shift($argv);
$exonerate = array_shift($argv);
$blat = array_shift($argv);

/// prepare fasta 

$wfst = "$fst.fwd-rev";
if(!is_file($wfst)) {
  system("$exonerate/bin/fastarevcomp $fst > $fst.rev"); /// reverse complement 
  system("sed -i 's|\ \[revcomp\]|_R|g' $fst.rev");                          /// ID for blat (_R = reverse)
  system("$exonerate/bin/fastareformat $fst > $wfst");   /// each line 70 chars
  system("cat $fst.rev >> $wfst");
  }
error("FWD+REV txs: $wfst",'i');

/// self blat

$psl = "$wfst.self.psl";
if(!is_file($psl)) {
  system("$blat -maxIntron=0 -noHead -maxGap=0 $wfst $wfst $psl"); /// restrictive blat
  }
error("PSL aln: $psl",'i');

$ms=array();
$me=array();

$fi = file($psl);
foreach($fi as $l) {
  $l=trim($l);
  if($l!='' && $l[0]!='#') {
    $a = preg_split("/\t/",$l); //print_r($a);
    if(($a[9]!=$a[13]) && $a[8]=='+') {
      //         -----1---->
      // ----2---->
      if($a[11]==0 && $a[14]==$a[16]) {  
        $ms[$a[9]][$a[13]]=$a[0];      
	$me[$a[13]][$a[9]]=$a[0];      
	}
      // -----1---->
      //         ----2---->
      elseif($a[15]==0 && $a[10]==$a[12]) { 
	$me[$a[9]][$a[13]]=$a[0];      
	$ms[$a[13]][$a[9]]=$a[0];      
	}
      $tblat[$a[9]]=true;
      $tblat[$a[13]]=true;
      }
    }
  }

/// locate chain ends

foreach($ms as $x => $a) {
  if(!isset($me[$x])) $v=0;
  else $v=count($me[$x]); 
  $c=count($a);
  //print "$x\t$c\t$v\n";
  if($c>0 && $v==0) $seed[$x]=true;
  }
  
error("Total txs w/ aln + = ".count($tblat),'i');  
error("End seeds = ".count($seed),'i');

/// get score (cvrg) of seeds and tx alns

foreach($seed as $x => $v) {
  if(preg_match("/\_cov\_(.*)\_g\d+\_i\d+/",$x,$m)) $cov=$m[1];
  else $cov=0;
  $seed[$x]=$cov;
  }
arsort($seed);

foreach($tblat as $x => $v) {
  if(preg_match("/\_cov\_(.*)\_g\d+\_i\d+/",$x,$m)) $cov=$m[1];
  else $cov=0;
  $tblat[$x]=$cov;
  }

/// get a greedy path with max cvrg

error("making chains from seeds...",'i');
$chain = array();
$c=0;
$ct = count($seed);
foreach($seed as $x => $v) {
  $chain[$c] = array();
  $n = $x;
  //print "seed = $n\n";
  $chain[$c][$n]=$tblat[$n];
  while($n!='') { 
    $cvrg = array(); 
    if(!isset($ms[$n])) $n='';
    else {
      foreach($ms[$n] as $y => $v) {
		$cvrg[$y]=$tblat[$y];
		}
      arsort($cvrg);
      $keys = array_keys($cvrg);
      $n = array_shift($keys);
      if(isset($chain[$c][$n])) $n=''; /// WARNING cos recursive could be found
      else $chain[$c][$n]=$tblat[$n];
      }
    }
  $c++;
  }

/// just get chains with decreasing cvrg

error("filtering chains...",'i');

$goodchains = array();
foreach($chain as $k => $a) {
  $fail=0;
  $cvrg = array_shift($a);
  foreach($a as $c) {
    if($c>$cvrg) $fail=1;  
    $cvrg = $c;
    }
    $goodchains[]=$k;
  }
$cok = count($goodchains);  

/// load sequences from fasta file

$seqs = readFasta($wfst);

/// get chained sequences

$ochains = "$fst.chains_info";
$fo = fopen($ochains,'w');
$ofst = "$fst.chains";
$fof = fopen($ofst,'w');
$c=0;
foreach($goodchains as $k) {
  $c++;
  fprintf($fo,">chain_$c"."_$k\n"); /// ntxs, length
  $a = $chain[$k];
  $ct=0;
  foreach($a as $n => $cvrg) {
    $ct++;
    fprintf($fo,"# $ct\t$n\t$cvrg\n");
    }
  $keys = array_keys($a);
  $x = array_shift($keys);
  $aseq = $seqs[$x];
  foreach($keys as $y) {
    $ovl = $ms[$x][$y];
    fprintf($fo,"$x\t$y\t$ovl\n");
    $aseq = substr($seqs[$y],0,-$ovl).$aseq;
    $x=$y;
    }
  $len = strlen($aseq);  
  fprintf($fof,">chain_$c"."_len_$len"."_txs_$ct\n".chunk_split($aseq,100,"\n"));
  }
fclose($fo);
fclose($fof);
error("chain info file: $ochains",'i');
error("chain fasta file: $ofst",'i');
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
	//$ID++;
	$ID=substr($l,1);
	$iSeq[$ID] ="";
	}
      }
    $l = fgets($buf);
    }
  fclose($buf);
  return($iSeq);
  }
  
