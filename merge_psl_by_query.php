#!/usr/local/bin/php
<?php
ini_set('memory_limit', '3G');

$prg  = array_shift($argv);
$file  = array_shift($argv);

$qlen = array();
$tlen = array();
$qcvrg = array();

$fi = file($file);
foreach($fi as $l) {
  $l=trim($l);
  if($l!='' && $l[0]!='#') {
    $a = preg_split("/\t/",$l); // print_r($a);
    $q = $a[9];
    $ql = $a[10];
    $qlen[$q] = $ql;
    $t = $a[13];
    $tlen[$t] = $a[14];
    if(!isset($qcvrg[$q]))   $qcvrg[$q] =array();
    if(!isset($qcvrg[$q][$t])) $qcvrg[$q][$t] =array();
    $qcvrg[$q][$t][] =array($a[11],$a[12]);
    }
  }
 
// print_r($qcvrg);  
foreach($qcvrg as $q =>$b) {
foreach($b as $t => $a) {
  $len = $qlen[$q];
  $qc = str_repeat(0,$len);
  foreach($a as $b) {
    for($j=$b[0];$j<=$b[1];$j++) $qc[$j]='1';
    }
  //print "$q\t$t\t[$qc]\n";
  $cc = count_chars($qc, 1);
  //print " -- ".$cc[49]."  -- \n";
  $cvrg = $cc[49];
  $pcvrg = $cc[49]/$len;
  print "$q\t$t\t$len\t".$tlen[$t]."\t".$cc[49]."\t$pcvrg\n";
  }}

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
