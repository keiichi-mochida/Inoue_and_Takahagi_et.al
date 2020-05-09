use strict;

my $Bd_sam=$ARGV[0];  # sam file (mapped to the B. distachyon genomem, Phytozome v3,1)
my $Bs_sam=$ARGV[1];  # sam file (mapped to the virtual B. stacei genome, Takahagi et al. GigaScience 2018)
my $out_dir=$ARGV[2]; # out dir
my $name=$ARGV[3];    # uniq name

my $BhBd="$out_dir/$name.BhBd.sam";
my $BhBs="$out_dir/$name.BhBs.sam";
my $Others="$out_dir/$name.Others.sam";

open(BHBD,">>$BhBd");
open(BHBS,">>$BhBs");
open(OTHERS,">>$Others");

my %Bd_reads=();
my %Bs_reads=();

my $all_reads=0;
my $BhBd_reads=0;
my $BhBs_reads=0;
my $Other_reads=0;
my $defective_reads=0;

open(A,"$Bd_sam");
while(my $line=<A>){
	$line=~s/\r//;
	chomp($line);
	
	if($line=~/^\@/){
		#print"$line\n";
		
		print BHBD "$line\n";
		print OTHERS "$line\n";
		
	}else{
		#print"$line\n";
		#print"$tmp[0]\n" if $line=~/RG:Z:NOID\.1/;
		#print"$tmp[0]\n" if $line=~/RG:Z:NOID\.2/;
		
		my @tmp=split("\t",$line);
		#print"$tmp[0]\n";
		
		my $key=$tmp[0];
		#print"$key\n";
		
		$line=~/MD:Z:.+/;
		my $nmi=&get_nmi($&);
		#print"$nmi\n";
		
		$Bd_reads{$key}[0]=$line if $line=~/RG:Z:NOID\.1/;
		$Bd_reads{$key}[1]=$line if $line=~/RG:Z:NOID\.2/;
		$Bd_reads{$key}[2]=$tmp[2] if $line=~/RG:Z:NOID\.1/;
		$Bd_reads{$key}[3]=$tmp[2] if $line=~/RG:Z:NOID\.2/;
		$Bd_reads{$key}[4]=$tmp[6] if $line=~/RG:Z:NOID\.1/;
		$Bd_reads{$key}[5]=$tmp[6] if $line=~/RG:Z:NOID\.2/;
		$Bd_reads{$key}[6]=$nmi if $line=~/RG:Z:NOID\.1/;
		$Bd_reads{$key}[7]=$nmi if $line=~/RG:Z:NOID\.2/;
		
	}
	
}
close(A);

open(B,"$Bs_sam");
while(my $line=<B>){
	$line=~s/\r//;
	chomp($line);
	
	if($line=~/^\@/){
		#print"$line\n";
		
		print BHBS "$line\n";
		
	}else{
		#print"$line\n";
		
		$all_reads++;
		
		my @tmp=split("\t",$line);
		#print"$tmp[0]\n";
		
		my $key=$tmp[0];
		#print"$key\n";
		
		$line=~/MD:Z:.+/;
		my $nmi=&get_nmi($&);
		#print"$nmi\n";
		
		if(exists($Bs_reads{$key})){
			
			my $tmp2_Bd_1=$Bd_reads{$key}[2];
			my $tmp2_Bd_2=$Bd_reads{$key}[3];
			my $tmp6_Bd_1=$Bd_reads{$key}[4];
			my $tmp6_Bd_2=$Bd_reads{$key}[5];
			
			my $tmp2_Bs_1="";
			my $tmp2_Bs_2="";
			my $tmp6_Bs_1="";
			my $tmp6_Bs_2="";
			if($line=~/RG:Z:NOID\.1/){
				
				$tmp2_Bs_1=$tmp[2];
				$tmp2_Bs_2=$Bs_reads{$key}[1];
				$tmp6_Bs_1=$tmp[6];
				$tmp6_Bs_2=$Bs_reads{$key}[2];
				
			}else{
				
				$tmp2_Bs_1=$Bs_reads{$key}[1];
				$tmp2_Bs_2=$tmp[2];
				$tmp6_Bs_1=$Bs_reads{$key}[2];
				$tmp6_Bs_2=$tmp[6];
				
			}
			
			my $check_Bd="";
			$check_Bd="ok" if $tmp2_Bd_1=~/Bd\d/ && $tmp2_Bd_2=~/Bd\d/ && $tmp6_Bd_1 eq "=" && $tmp6_Bd_2 eq "=";
			my $check_Bs="";
			$check_Bs="ok" if $tmp2_Bs_1=~/Bd\d/ && $tmp2_Bs_2=~/Bd\d/ && $tmp6_Bs_1 eq "=" && $tmp6_Bs_2 eq "=";
			
			if($check_Bd eq "ok" && $check_Bs eq "ok"){
				
				my $nmi_Bd_1=$Bd_reads{$key}[6];
				my $nmi_Bd_2=$Bd_reads{$key}[7];
				my $nmi_Bs_1="";
				my $nmi_Bs_2="";
				if($line=~/RG:Z:NOID\.1/){
					
					$nmi_Bs_1=$nmi;
					$nmi_Bs_2=$Bs_reads{$key}[3];
					
				}else{
					
					$nmi_Bs_1=$Bs_reads{$key}[3];
					$nmi_Bs_2=$nmi;
					
				}
				#print"$nmi_Bd_1\t$nmi_Bd_2\t$nmi_Bs_1\t$nmi_Bs_2\n";
				
				if($nmi_Bd_1 < $nmi_Bs_1 && $nmi_Bd_2 < $nmi_Bs_2){
					
					$BhBd_reads+=2;
					print BHBD "$Bd_reads{$key}[0]\n$Bd_reads{$key}[1]\n";
					
				}elsif($nmi_Bd_1 eq $nmi_Bs_1 && $nmi_Bd_2 < $nmi_Bs_2){
					
					$BhBd_reads+=2;
					print BHBD "$Bd_reads{$key}[0]\n$Bd_reads{$key}[1]\n";
					
				}elsif($nmi_Bd_1 < $nmi_Bs_1 && $nmi_Bd_2 eq $nmi_Bs_2){
					
					$BhBd_reads+=2;
					print BHBD "$Bd_reads{$key}[0]\n$Bd_reads{$key}[1]\n";
					
				}elsif($nmi_Bd_1 > $nmi_Bs_1 && $nmi_Bd_2 > $nmi_Bs_2){
					
					$BhBs_reads+=2;
					print BHBS "$Bs_reads{$key}[0]\n$line\n";
					
				}elsif($nmi_Bd_1 eq $nmi_Bs_1 && $nmi_Bd_2 > $nmi_Bs_2){
					
					$BhBs_reads+=2;
					print BHBS "$Bs_reads{$key}[0]\n$line\n";
					
				}elsif($nmi_Bd_1 > $nmi_Bs_1 && $nmi_Bd_2 eq $nmi_Bs_2){
					
					$BhBs_reads+=2;
					print BHBS "$Bs_reads{$key}[0]\n$line\n";
					
				}else{
					
					$Other_reads+=2;
					print OTHERS "$Bd_reads{$key}[0]\n$Bd_reads{$key}[1]\n";
					
				}
				
			}elsif($check_Bd eq "ok" && $check_Bs eq ""){
				
				$BhBd_reads+=2;
				print BHBD "$Bd_reads{$key}[0]\n$Bd_reads{$key}[1]\n";
				
			}elsif($check_Bd eq "" && $check_Bs eq "ok"){
				
				$BhBs_reads+=2;
				print BHBS "$Bs_reads{$key}[0]\n$line\n";
				
			}else{
				
				$defective_reads+=2;
				
			}
			
			delete($Bd_reads{$key});
			delete($Bs_reads{$key});
			
		}else{
			
			$Bs_reads{$key}[0]=$line;
			$Bs_reads{$key}[1]=$tmp[2];
			$Bs_reads{$key}[2]=$tmp[6];
			$Bs_reads{$key}[3]=$nmi;
			
		}
		
		
	}
	
}
close(B);

close(BHBD);
close(BHBS);
close(OTHERS);

#print"$name\t$all_reads\t$BhBd_reads\t$BhBs_reads\t$Other_reads\t$defective_reads\n";

sub get_nmi{
	
	my $nmi=$_[0];
	#print"$nmi\n";
	$nmi=~s/MD:Z://;
	#print"$nmi\n";
	$nmi=~s/\t.+//;
	#print"$nmi\n";
	$nmi=~s/\^\D+//g;
	#print"$nmi\n";
	$nmi=(()=$nmi=~/\D/g);
	#print"$nmi\n";
	
	return $nmi;
	
}
