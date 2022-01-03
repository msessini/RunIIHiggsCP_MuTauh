#! /usr/bin/perl
use Cwd;
use POSIX;
use POSIX qw(strftime);

#############################################
$numArgs = $#ARGV +1;
$ARGV[$argnum];

$UserID= POSIX::cuserid();
$UserIDCern=$UserID;
$time= strftime("%h_%d_%Y",localtime);
$SkimName = "SkimmedNtuples_$time";
if($ARGV[0] eq "--help" || $ARGV[0] eq ""){
    printf("\n ===========>>>>  Hello  $UserName ! <<<<=============\n\n");
    printf("\n Options:                                                                                                       ");
    printf("\n                                --TransferSkimed");
    printf("\n                                Will transfer Skiimed ntuples to Strasbourg Tier3 ");
    printf("\n                                This script is incomplete, contact Vladimir Cherepanov if you want to use it  \n\n\n");
    exit(0);  
} 
if($ARGV[$l] eq  "--DirName"){
    $l++;
    $SkimName=$ARGV[$l];
} 

my $dir = getcwd;

$temp= $set . $time;
$set=$temp;
 system(sprintf("  rm  junk0 "));
 system(sprintf("  rm  junk1 "));
 system(sprintf(" touch  junk0 "));
 system(sprintf(" touch  junk1 "));
$Tier3path = "/dpm/in2p3.fr/home/cms/phedex/store/user/cherepan/";
 if( $ARGV[0] eq "--TransferSkimed" ){
    
     system(sprintf(" rfmkdir $Tier3path$SkimName"));
     
     system(sprintf(" ls -l  $dir  >>  junk0 "));
     system(sprintf("cat junk0 | awk '{print \$9}' >& junk1"));
     open(DAT, "junk1");
     while ($item = <DAT>) {
	 chomp($item);
	 if( ($item =~ m/Set/g)){
	     push(@sets,$item);
	 }
     }
     system(sprintf("  rm  junk0 "));
     system(sprintf("  rm  junk1 "));
     $nsets = @sets;
     for($l=1;$l<=$nsets; $l++){
	 opendir(DIR, "$dir/Set_$l");
	 @skfiles = grep(/SKIMMED_NTUP/,readdir(DIR));
	 closedir(DIR);
	 if(@skfiles){  print "Skimmed File is available in $dir/Set_$l:   \n"; push(@Skimmed,$l);
	 } else { print "There is no SKIMMED ntuple in $dir/Set_$l \n"; push(@isNotSkimmed,$l);}
     }
     print "\n  Skimmed :  @Skimmed  \n\n\n   ";
     print " Not  Skimmed :  @isNotSkimmed  \n\n\n   ";
     foreach $s (@Skimmed){
	 open(INPUT,"$dir/Set_$s/Input.txt")  || die "can't open Input  file $dir/Set_$l/Input.txt" ;
	 while (<INPUT>) {
	     ($i1,$i2,$i3,$i4)=split(/ /,$_);
	     if($i1 eq "InputNtuples:" ){# print "Set_$s/  :  $i2  ";
		 ($p1,$p2,$p3,$p4,$p5,$p6,$p7,$p8)=split(/\//,$i2);
	     }
	 }
	 close(INPUT);
	 print "Transfer:  Set_$s/SKIMMED_NTUP.root   to   $SkimName/$p5/$p6/$p7/0000/SKIMMED_NTUP_$s.root \n\n\n";
	 system(sprintf(" rfmkdir $Tier3path$SkimName/$p5"));
	 system(sprintf(" rfmkdir $Tier3path$SkimName/$p5/$p6"));
	 system(sprintf(" rfmkdir $Tier3path$SkimName/$p5/$p6/$p7/"));
	 system(sprintf(" rfmkdir $Tier3path$SkimName/$p5/$p6/$p7/0000"));
	 system(sprintf(" rfcp  $dir/Set_$s/SKIMMED_NTUP.root  $Tier3path$SkimName/$p5/$p6/$p7/0000/SKIMMED_NTUP_$s.root"));
     }
}
