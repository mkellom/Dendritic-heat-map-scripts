#!usr/bin/perl

$generations=$ARGV[0];
@bases=("A","T","G","C");
@codons=("GCA","GCC","GCG","GCT","TGC","TGT","GAC","GAT","GAA","GAG","TTC","TTT","GGA","GGC","GGG","GGT","CAC","CAT","ATA","ATC","ATT","AAA","AAG","CTA","CTC","CTG","CTT","TTA","TTG","ATG","AAC","AAT","CCA","CCC","CCG","CCT","CAA","CAG","AGA","AGG","CGA","CGC","CGG","CGT","AGC","AGT","TCA","TCC","TCG","TCT","ACA","ACC","ACG","ACT","GTA","GTC","GTG","GTT","TGG","TAC","TAT","TAA","TAG","TGA");
$switch=0;
while ($switch<=0){
	undef $seq;
	undef @seqs;
	undef @fakeA;
	undef @fakeB;
	undef @fakeC;
	undef @tempA;
	undef @tempB;
	undef @tempC;
	$seq.=$bases[rand @bases] for 1..100;
	push(@seqs,$seq);
	for($i=1;$i<=$generations;$i++){
		foreach $seq (@seqs){
			$pos=int(rand(100));
			@letters=split(//,$seq);
			splice @letters, $pos, 1;
			$mute=$bases[rand @bases];
			splice @letters, $pos, 0, "$mute";
			$seq1=join('',@letters);
			$As=()=$seq=~/A/g;
			$Ts=()=$seq=~/T/g;
			$Gs=()=$seq=~/G/g;
			$Cs=()=$seq=~/C/g;
			$gc=(($Gs+$Cs)/($As+$Ts+$Gs+$Cs))*100;
			if ($gc<=50){
				push(@fakeC,$seq);
				push(@fakeA,$seq);
			}else{
				push(@fakeC,$seq);
				push(@fakeB,$seq);
			}
			$As=()=$seq1=~/A/g;
			$Ts=()=$seq1=~/T/g;
			$Gs=()=$seq1=~/G/g;
			$Cs=()=$seq1=~/C/g;
			$gc=(($Gs+$Cs)/($As+$Ts+$Gs+$Cs))*100;
			if ($gc<=50){
				push(@fakeC,$seq1);
				push(@fakeA,$seq1);
			}else{
				push(@fakeC,$seq1);
				push(@fakeB,$seq1);
			}
			push(@newseqs,$seq);
			push(@newseqs,$seq1);
		}
		$sizeA=scalar @fakeA;
		$sizeB=scalar @fakeB;
		$sizeC=scalar @fakeC;
		$diff=abs($sizeA-$sizeB);
		if($diff<$sizeC/2){
			@seqs=@newseqs;
			undef @newseqs;
		}else{	
			print "Uneven Growth at Generation $i; Restarting\n";
			undef @fakeA;
			undef @fakeB;
			undef @fakeC;
			undef @newseqs;
			if ($i<=1){
				undef $seq;
				undef @seqs;
				$seq.=$bases[rand @bases] for 1..100;
				push(@seqs,$seq);
			}
			redo;
		}
	}
	$balance=$sizeC*0.05;
	$diff=abs($sizeA-$sizeB);
	print "$balance:$diff;";
	if ($diff<=$balance){
		print "\n";
		$fake1=$fakeC[0];
		$fake2=$fakeC[1];
		$As=()=$fake1=~/A/g;
		$Ts=()=$fake1=~/T/g;
		$Gs=()=$fake1=~/G/g;
		$Cs=()=$fake1=~/C/g;
		$gc1=(($Gs+$Cs)/($As+$Ts+$Gs+$Cs))*100;
		$As=()=$fake2=~/A/g;
		$Ts=()=$fake2=~/T/g;
		$Gs=()=$fake2=~/G/g;
		$Cs=()=$fake2=~/C/g;
		$gc2=(($Gs+$Cs)/($As+$Ts+$Gs+$Cs))*100;
		if (($gc1<=50 && $gc2>50) || ($gc1>50 && $gc2<=50)){
			print "Seq1 ${gc1}\%\nSeq2 ${gc2}\%\n";
			for($i=0;$i<=$generations;$i++){
				system("mkdir -p -m 777 $i");
				$exponent=2**$i;
				for($e=0;$e<=$exponent-1;$e++){
					$fake=$fakeC[$e];
					$As=()=$fake=~/A/g;
					$Ts=()=$fake=~/T/g;
					$Gs=()=$fake=~/G/g;
					$Cs=()=$fake=~/C/g;
					$gc=(($Gs+$Cs)/($As+$Ts+$Gs+$Cs))*100;
					if ($gc<=50){
						open(OUTFILE,">>${i}/fake$i.fasta");
						print OUTFILE ">fake_${e}_${gc};\n$fake\n";
						close OUTFILE;
						open(OUTFILE,">>${i}/fakeA$i.fasta");
						print OUTFILE ">fake_${e}_${gc};\n$fake\n";
						close OUTFILE;
						open(OUTFILE,">>${i}/fakeB$i.fasta");
						close OUTFILE;
					}else{
						open(OUTFILE,">>${i}/fake$i.fasta");
						print OUTFILE ">fake_${e}_${gc};\n$fake\n";
						close OUTFILE;
						open(OUTFILE,">>${i}/fakeB$i.fasta");
						print OUTFILE ">fake_${e}_${gc};\n$fake\n";
						close OUTFILE;
						open(OUTFILE,">>${i}/fakeA$i.fasta");
						close OUTFILE;
					}
				}
				print "Generation ${i} Sequences Written\n";
			}
			$switch=1;
		}else{
			print "Seq1 ${gc1}\%\nSeq2 ${gc2}\%\nNOT POLAR; Restarting\n";
		}
	}else{
		print " Difference exceeds limit; Restarting;\n";
	}
}
print "\n";
open(FILE,"<${generations}/fake${generations}.fasta");
@fastas=<FILE>;
close FASTA;
$fasta=join('',@fastas);
$total_population=()=$fasta=~/>/g;
$shifter=int($total_population*0.01);
sub log10{
	$n=shift;
	return log($n)/log(10);
}
open(FILE,"<faux.conf");
@lines=<FILE>;
close FILE;
$faux_draft=join('',@lines);
$min=1;
$max=-1;
for($g=0;$g<=$generations;$g++){
	system("mkdir -p -m 777 $g/75");
	undef @og_orders;
	undef @go_orders;
	print "\n";
	if ($g>=1){
		system("/mnt/datA/Applications/clustalo-1.2.0 -i ${g}/fake${g}.fasta -t DNA --infmt=fasta --full -o ${g}/fake${g}.aligned.fasta --outfmt=fasta --output-order=tree-order --threads 60 --verbose --force");
		open(FILE,"<${g}/fake${g}.aligned.fasta");
		@alignments=<FILE>;
		close FILE;
	}else{
		open(FILE,"<${g}/fake${g}.fasta");
		@alignments=<FILE>;
		close FILE;
	}
	foreach  $alignment (@alignments){
		next if ($alignment!~/>/);
		chomp $alignment;
		$alignment=~s/>//g;
		push(@og_orders,"$alignment");
	}
	@temp_orders=@og_orders;
	while (@temp_orders){
		$left=shift @temp_orders;
		push(@go_orders,$left);
		$right=pop @temp_orders;
		if ($right=~/fake/){
			push(@go_orders,$right);
		}
	}
	open(FILE,"<${g}/fake${g}.fasta");
	@lines=<FILE>;
	close FILE;
	$single=join('',@lines);
	@originals=split(/>/,$single);
	shift @originals;
	foreach $go (@go_orders){
		foreach $ori (@originals){
			if ($ori=~/$go/){
				open(OUTFILE,">>${g}/fake${g}.corrected.fasta");
				print OUTFILE ">$ori";
				close OUTFILE;
			}
		}
	}
	system("/mnt/datA/Applications/usearch7.0.1090_i86linux64 -derep_fulllength ${g}/fake${g}.corrected.fasta -output ${g}/derep${g}.fasta -sizeout -threads 60");
	system("/mnt/datA/Applications/usearch7.0.1090_i86linux64 -cluster_fast ${g}/fake${g}.corrected.fasta -id 0.75 -msaout $g/75/fake_0_75.msa -threads 60");
	open(FILE,"<$g/75/fake_0_75.msa");
	@lines=<FILE>;
	close FILE;
	$single=join('',@lines);
	$single=~s/-//g;
	@otus=split(/>\*/,$single);
	shift @otus;
	$q=0;
	foreach $otu (@otus){
		($seqs,$cons)=split(/\n>consensus/,$otu,2);
		open(OUTFILE,">>${g}/75/fake_${q}_75.fasta");
		print OUTFILE ">$seqs";
		close OUTFILE;
		open(OUTFILE,">>${g}/75/fake_75.clusters");
		print OUTFILE ">*$seqs\n";
		close OUTFILE;
		$q++;
	}
	for($i=76;$i<=99;$i++){
		system("mkdir -p -m 777 ${g}/$i");
		$f=$i-1;
		@files=<$g/$f/*>;
		$p=0;
		$q=0;
		foreach $file (@files){
			if ($file=~/\.fasta/){
				($name,$ext)=split(/\./,$file,2);
				open(FILE,"<$file");
				@lines=<FILE>;
				close FILE;
				$single=join('',@lines);
				@cluster_seqs=split(/>/,$single);
				shift @cluster_seqs;
				undef @go_orders;
				foreach $og (@og_orders){
					foreach $cluster_seq (@cluster_seqs){
						if ($cluster_seq=~/$og/){
							chomp $cluster_seq;
							push(@go_orders,$cluster_seq);
						}
					}
				}
				while (@go_orders){
					$left=shift @go_orders;
					open(OUTFILE,">>$name.corrected.fasta");
					print OUTFILE ">$left\n";
					close OUTFILE;
					$right=pop @go_orders;
					if ($right=~/fake/){
						open(OUTFILE,">>$name.corrected.fasta");
						print OUTFILE ">$right\n";
						close OUTFILE;
					}
				}
				system("/mnt/datA/Applications/usearch7.0.1090_i86linux64 -cluster_fast $name.corrected.fasta -id 0.$i -msaout ${g}/${i}/fake_${p}_${i}.msa -threads 60");
				open(FILE,"<${g}/${i}/fake_${p}_${i}.msa");
				@lines=<FILE>;
				close FILE;
				$single=join('',@lines);
				$single=~s/-//g;
				@otus=split(/>\*/,$single);
				shift @otus;
				foreach $otu (@otus){
					($seqs,$cons)=split(/\n>consensus/,$otu,2);
					open(OUTFILE,">>${g}/${i}/fake_${q}_${i}.fasta");
					print OUTFILE ">$seqs";
					close OUTFILE;
					open(OUTFILE,">>${g}/${i}/fake_${i}.clusters");
					print OUTFILE ">*$seqs\n";
					close OUTFILE;
					$q++;
				}
				$p++;
			}
		}
	}
	@files=<${g}/99/*>;
	$p=0;
	$q=0;
	system("mkdir -p -m 777 ${g}/100");
	foreach $file (@files){
		if ($file=~/\.fasta/){
			($name,$ext)=split(/\./,$file,2);
			open(FILE,"<$file");
			@lines=<FILE>;
			close FILE;
			$single=join('',@lines);
			@cluster_seqs=split(/>/,$single);
			shift @cluster_seqs;
			undef @go_orders;
			foreach $og (@og_orders){
				foreach $cluster_seq (@cluster_seqs){
					if ($cluster_seq=~/$og/){
						chomp $cluster_seq;
						push(@go_orders,$cluster_seq);
					}
				}
			}
			while (@go_orders){
				$left=shift @go_orders;
				open(OUTFILE,">>$name.corrected.fasta");
				print OUTFILE ">$left\n";
				close OUTFILE;
				$right=pop @go_orders;
				if ($right=~/fake/){
					open(OUTFILE,">>$name.corrected.fasta");
					print OUTFILE ">$right\n";
					close OUTFILE;
				}
			}
			system("/mnt/datA/Applications/usearch7.0.1090_i86linux64 -cluster_fast $name.corrected.fasta -id 1.0 -msaout ${g}/100/fake_${p}_100.msa -threads 60");
			open(FILE,"<${g}/100/fake_${p}_100.msa");
			@lines=<FILE>;
			close FILE;
			$single=join('',@lines);
			$single=~s/-//g;
			@otus=split(/>\*/,$single);
			shift @otus;
			foreach $otu (@otus){
				($seqs,$cons)=split(/\n>consensus/,$otu,2);
				open(OUTFILE,">>${g}/100/fake_${q}_100.fasta");
				print OUTFILE ">$seqs";
				close OUTFILE;
				open(OUTFILE,">>${g}/100/fake_100.clusters");
				print OUTFILE ">*$seqs\n";
				close OUTFILE;
				$q++;
			}
			$p++;
		}
	}
	@orders=@og_orders;
	system("mkdir -p -m 777 circos_inputs");
	system("mkdir -p -m 777 circos_inputs/$g");
	open(FILE,"<${g}/fake${g}.fasta");
	@fastas=<FILE>;
	close FILE;
	$fasta=join('',@fastas);
	$population=()=$fasta=~/>/g;
	for($i=75;$i<=100;$i++){
		@files=<${g}/$i/*>;
		$start=$shifter*2;
		$written="NULL";
		undef @neworders;
		open(FILE,"${g}/${i}/fake_${i}.clusters");
		@lines=<FILE>;
		close FILE;
		$single=join('',@lines);
		@clusters=split(/>\*/,$single);
		shift @clusters;
		foreach $order (@orders){
			next if ($written=~/$order/);
			$order=~s/;//g;
			foreach $cluster (@clusters){
				if ($cluster=~/$order/){
					chomp $cluster;
					$Acount=0;
					$Bcount=0;
					$total=0;
					@seqs=split(/>/,$cluster);
					foreach $order2 (@orders){
						foreach $seq (@seqs){
							if ($seq=~/$order2/){
								if($g>=1){
									$written.="$order2;";
									push(@neworders,$order2);
									($fake_tag,$code)=split(/\n/,$seq,2);
									$code=~s/\n//g;
									open(FILE,"<${g}/fakeA${g}.fasta");
									@checkAs=<FILE>;
									close FILE;
									$singleA=join('',@checkAs);
									if ($singleA=~/$code/){
										$Acount++;
									}
									open(FILE,"<${g}/fakeB${g}.fasta");
									@checkBs=<FILE>;
									close FILE;
									$singleB=join('',@checkBs);
									if ($singleB=~/$code/){
										$Bcount++;
									}
									open(FILE,"<${g}/derep${g}.fasta");
									@dereps=<FILE>;
									close FILE;
									$single=join('',@dereps);
									@dereps=split(/>/,$single);
									shift @dereps;
									foreach $derep (@dereps){
										$derep=~s/\n//g;
										if ($derep=~/$code/){
											($tag,$rep,$nucs)=split(/;/,$derep,3);
											$rep=~s/size=//g;
											$total=$total+$rep;
										}
									}
								}else{
									$written.="$order2;";
									push(@neworders,$order2);
									$total=1;
								}
							}
						}
					}
					$skew=log10(($Acount+1)/($Bcount+1));
					$skew=sprintf("%.3f",$skew);
					$abs_skew=abs($skew);
					$neg_skew=(0-$abs_skew);
					$pos_skew=$abs_skew;
					if ($neg_skew<=$min){
						$min=$neg_skew;
					}
					if ($pos_skew>=$max){
						$max=$pos_skew;
					}
					$total=$total*$total_population/$population;
					$stop=$start+$total;
					open(OUTFILE,">>circos_inputs/${g}/faux_$i.txt");
					print OUTFILE "faux $start $stop $skew\n";
					close OUTFILE;
					$start=$stop;
				}
			}
		}
		@orders=@neworders;
		print "${g}/$i Complete\n";
	}
}
print "\n";
system("mkdir -p -m 777 confs");
system("mkdir -p -m 777 circos_images");
$start=$start+$shifter;
$stop=$start+$shifter;
for($g=0;$g<=$generations;$g++){
	for($i=75;$i<=100;$i++){
		open(OUTFILE,">>circos_inputs/${g}/faux_$i.txt");
		print OUTFILE "faux 0 $shifter $min\n";
		close OUTFILE;
		open(OUTFILE,">>circos_inputs/${g}/faux_$i.txt");
		print OUTFILE "faux $start $stop $max\n";
		close OUTFILE;
	}
	open(OUTFILE,">>circos_inputs/${g}/faux_karyotype.txt");
	print OUTFILE "chr - faux faux 0 $stop black";
	close OUTFILE;
	print "${g} Frame Complete\n";
	$faux=$faux_draft;
	$faux=~s/generation/${g}/g;
	open(OUTFILE,">>confs/faux${g}.conf");
	print OUTFILE "$faux";
	close OUTFILE;
	print "${g} Configuration Complete\n\n";
	system("perl /mnt/datA/matt/dendritic_heat_maps/circos-0.64/bin/circos -conf /mnt/datA/matt/software/confs/faux${g}.conf");
	print "\n";
}
