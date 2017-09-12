#!usr/bin/perl

$leaps=$ARGV[0];
@bases=("A","T","G","C");
@linkages=("min","max","avg");
$switch="off";

# Creates FASTA files for mutation series {

while ($switch eq "off"){
	$seq.=$bases[rand @bases] for 1..100;
	$As=()=$seq=~/A/g;
	$Ts=()=$seq=~/T/g;
	$Gs=()=$seq=~/G/g;
	$Cs=()=$seq=~/C/g;
	$gc=(($Gs+$Cs)/($As+$Ts+$Gs+$Cs))*100;
	if ($gc==50){
		print "Initial Sequence GC Content: ${gc}\%\n$seq\n";
		$switch="on";
	}else{
		print "Initial Sequence GC Content: ${gc}\% ; Restarting\n";
		undef $seq;
	}
}
for($i=0;$i<=$leaps;$i++){
	for($q=0;$q<=99;$q++){
		if ($i>=1){
			$seq=$pop{$q};
			$pos=int(rand(100));
			@letters=split(//,$seq);
			splice @letters, $pos, 1;
			$mute=$bases[rand @bases];
			splice @letters, $pos, 0, "$mute";
			$seq=join('',@letters);
		}
		push(@news,"$q,$seq");
		$As=()=$seq=~/A/g;
		$Ts=()=$seq=~/T/g;
		$Gs=()=$seq=~/G/g;
		$Cs=()=$seq=~/C/g;
		$gc=(($Gs+$Cs)/($As+$Ts+$Gs+$Cs))*100;
		if ($gc<=50){
			push(@Agroup,$seq);
		}else{
			push(@Bgroup,$seq);
		}
	}
	if ($i>=1){
		$Asize=scalar @Agroup;
		$Bsize=scalar @Bgroup;
		$diff=abs($Asize-$Bsize);
		if($diff<=20){
			print "Mutation ${i} Done ; Difference=$diff\n";
			undef @Agroup;
			undef @Bgroup;
			foreach $new (@news){
				($k,$v)=split(/,/,$new,2);
				$pop{$k}=$v;
				$As=()=$v=~/A/g;
				$Ts=()=$v=~/T/g;
				$Gs=()=$v=~/G/g;
				$Cs=()=$v=~/C/g;
				$gc=(($Gs+$Cs)/($As+$Ts+$Gs+$Cs))*100;
				system("mkdir -p -m 777 ${i}");
				open(OUTFILE,">>${i}/fake${i}.fasta");
				print OUTFILE ">fake_${k}_${i}_${gc};\n$v\n";
				close OUTFILE;
			}
			undef @news;
		}else{
			print "Uneven Groups at Mutation $i ; Difference=$diff\n";
			undef @Agroup;
			undef @Bgroup;
			undef @news;
			redo;
		}
	}else{
		print "Mutation ${i} Done\n";
		undef @Agroup;
		undef @Bgroup;
		foreach $new (@news){
			($k,$v)=split(/,/,$new,2);
			$pop{$k}=$v;
			$As=()=$v=~/A/g;
			$Ts=()=$v=~/T/g;
			$Gs=()=$v=~/G/g;
			$Cs=()=$v=~/C/g;
			$gc=(($Gs+$Cs)/($As+$Ts+$Gs+$Cs))*100;
			system("mkdir -p -m 777 ${i}");
			open(OUTFILE,">>${i}/fake${i}.fasta");
			print OUTFILE ">fake_${k}_${i}_${gc};\n$v\n";
			close OUTFILE;
		}
		undef @news;
	}
}

# }
# log base 10 subroutine {

sub log10{
	$n=shift;
	return log($n)/log(10);
}

# }
# Sets min/max and multiple alignment for mutation series {

$min=1;
$max=-1;
for($i=0;$i<=$leaps;$i++){
	system("/mnt/datA/Applications/clustalo-1.2.0 -i ${i}/fake${i}.fasta -t DNA --infmt=fasta --full -o ${i}/fake${i}.aligned.fasta --outfmt=fasta --output-order=tree-order --threads 60 --verbose --force");
}
print "\nAlignments Complete\n";

# }
# Bottom-up clustering and Circos files creation {

open(FILE,"<${leaps}/fake${leaps}.fasta");
@fastas=<FILE>;
close FASTA;
$fasta=join('',@fastas);
$total_population=()=$fasta=~/>/g;
$shifter=sprintf("%.0f",$total_population*0.01);
system("mkdir -p -m 777 bottom_up");
foreach $linkage (@linkages){
	system("mkdir -p -m 777 bottom_up/${linkage}");
	system("mkdir -p -m 777 bottom_up/${linkage}/circos_inputs");
	for($i=0;$i<=$leaps;$i++){
		system("mkdir -p -m 777 bottom_up/${linkage}/${i}");
		undef @orders;
		print "\n";
		open(FILE,"<${i}/fake${i}.aligned.fasta");
		@lines=<FILE>;
		close FILE;
		foreach $order (@lines){
			chomp $order;
			if ($order=~/>/){
				$order=~s/>//g;
				push(@orders,$order);
			}
		}
		for($q=75;$q<=99;$q++){
			system("/mnt/datA/Applications/usearch8.0.1517_i86linux64 -cluster_agg ${i}/fake${i}.fasta -clusterout bottom_up/${linkage}/${i}/clusters${q}.txt -id 0.${q} -linkage $linkage -threads 60");
		}
		system("/mnt/datA/Applications/usearch8.0.1517_i86linux64 -cluster_agg ${i}/fake${i}.fasta -clusterout bottom_up/${linkage}/${i}/clusters100.txt -id 1.0 -linkage $linkage -threads 60");
		print "\n";
		system("mkdir -p -m 777 bottom_up/${linkage}/circos_inputs/${i}");
		open(FILE,"<${i}/fake${i}.fasta");
		@fastas=<FILE>;
		close FILE;
		$fasta=join('',@fastas);
		$population=()=$fasta=~/>/g;
		for($q=75;$q<=100;$q++){
			system("mkdir -p -m 777 bottom_up/${linkage}/${i}/${q}");
			$start=$shifter*2;
			open(FILE,"<bottom_up/${linkage}/${i}/clusters${q}.txt");
			@members=<FILE>;
			close FILE;
			foreach $member (@members){
				chomp $member;
				($cluster,$tag)=split(/\t/,$member,2);
				open(OUTFILE,">>bottom_up/${linkage}/${i}/${q}/cluster$cluster.txt");
				print OUTFILE "$tag\n";
				close OUTFILE;
			}
			@clusters=<bottom_up/${linkage}/${i}/${q}/*>;
			foreach $order (@orders){
				next if ($written=~/$order/);
				foreach $cluster (@clusters){
					open(FILE,"<$cluster");
					@seqs=<FILE>;
					close FILE;
					$single=join('',@seqs);
					if ($single=~/$order/){
						$Acount=0;
						$Bcount=0;
						$total=0;
						foreach $order2 (@orders){
							foreach $seq (@seqs){
								chomp $seq;
								if ($order2 eq $seq){
									$written.="$seq";
									push(@neworders,$seq);
									if ($i>=1){
										@parts=split(/_/,$seq);
										$gc=$parts[3];
										if($gc<=50){
											$Acount++;
										}else{
											$Bcount++;
										}
									}
								}
							}
						}
						if ($i>=1){
							$total=$Acount+$Bcount;
						}else{
							$total=scalar @seqs;
						}
						$skew=log10(($Acount+1)/($Bcount+1));
						undef $Acount;
						undef $Bcount;
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
						$stop=$start+$total;
						$stop=sprintf("%.0f",$stop);
						open(OUTFILE,">>bottom_up/${linkage}/circos_inputs/${i}/faux_$q.txt");
						print OUTFILE "faux $start $stop $skew\n";
						close OUTFILE;
						$start=$stop;
					}
				}
			}
			@orders=@neworders;
			undef @neworders;
			undef $written;
			print "bottom_up/${linkage}/${i}/${q} Complete\n";
		}
	}
	print "\n";
	open(FILE,"<faux.conf");
	@lines=<FILE>;
	close FILE;
	$faux_draft=join('',@lines);
	system("mkdir -p -m 777 bottom_up/${linkage}/confs");
	$start=$start+1;
	$stop=$start+1;
	for($g=0;$g<=$leaps;$g++){
		for($i=75;$i<=100;$i++){
			open(OUTFILE,">>bottom_up/${linkage}/circos_inputs/${g}/faux_$i.txt");
			print OUTFILE "faux 0 1 $min\n";
			close OUTFILE;
			open(OUTFILE,">>bottom_up/${linkage}/circos_inputs/${g}/faux_$i.txt");
			print OUTFILE "faux $start $stop $max\n";
			close OUTFILE;
		}
		open(OUTFILE,">>bottom_up/${linkage}/circos_inputs/${g}/faux_karyotype.txt");
		print OUTFILE "chr - faux faux 0 $stop black";
		close OUTFILE;
		print "bottom_up/${linkage}/${g} Frame Complete\n";
		$faux=$faux_draft;
		$faux=~s/generation/${g}/g;
		$faux=~s/linkage/${linkage}/g;
		$faux=~s/direction/bottom_up/g;
		open(OUTFILE,">>bottom_up/${linkage}/confs/faux${g}.conf");
		print OUTFILE "$faux";
		close OUTFILE;
		print "bottom_up/${linkage}/${g} Configuration Complete\n\n";
	}
}
print "\nBottom-Up Complete\n";

# }
# Top-down fake series clustering and Circos files creation {

open(FILE,"<faux.conf");
@lines=<FILE>;
close FILE;
$faux_draft=join('',@lines);
system("mkdir -p -m 777 top_down");
system("mkdir -p -m 777 top_down/circos_inputs");
for($g=0;$g<=$leaps;$g++){
	print "\n";
	undef @go_orders;
	system("mkdir -p -m 777 top_down/${g}");
	system("mkdir -p -m 777 top_down/${g}/75");
	open(FILE,"<${g}/fake${g}.aligned.fasta");
	@lines=<FILE>;
	close FILE;
	foreach $order (@lines){
		chomp $order;
		if ($order=~/>/){
			$order=~s/>//g;
			push(@og_orders,$order);
		}
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
	system("/mnt/datA/Applications/usearch7.0.1090_i86linux64 -derep_fulllength ${g}/fake${g}.corrected.fasta -output top_down/${g}/derep${g}.fasta -sizeout -threads 60");
	open(FILE,"<top_down/${g}/derep${g}.fasta");
	@lines=<FILE>;
	close FILE;
	$single=join('',@lines);
	$single=~s/A\n/A/g;
	$single=~s/T\n/T/g;
	$single=~s/G\n/G/g;
	$single=~s/C\n/C/g;
	$single=~s/>/\n>/g;
	@dereps=split(/>/,$single);
	shift @dereps;
	system("/mnt/datA/Applications/usearch7.0.1090_i86linux64 -cluster_fast ${g}/fake${g}.corrected.fasta -id 0.75 -msaout top_down/$g/75/fake_0_75.msa -threads 60");
	open(FILE,"<top_down/$g/75/fake_0_75.msa");
	@lines=<FILE>;
	close FILE;
	$single=join('',@lines);
	$single=~s/-//g;
	@otus=split(/>\*/,$single);
	shift @otus;
	$q=0;
	foreach $otu (@otus){
		($seqs,$cons)=split(/\n>consensus/,$otu,2);
		$seqs=~s/A\n/A/g;
		$seqs=~s/T\n/T/g;
		$seqs=~s/G\n/G/g;
		$seqs=~s/C\n/C/g;
		$seqs=~s/>/\n>/g;
		open(OUTFILE,">>top_down/${g}/75/fake_${q}_75.fasta");
		print OUTFILE ">$seqs";
		close OUTFILE;
		open(OUTFILE,">>top_down/${g}/75/fake_75.clusters");
		print OUTFILE ">*$seqs\n";
		close OUTFILE;
		$q++;
	}
	for($i=76;$i<=99;$i++){
		system("mkdir -p -m 777 top_down/${g}/$i");
		$f=$i-1;
		@files=<top_down/$g/$f/*>;
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
				system("/mnt/datA/Applications/usearch7.0.1090_i86linux64 -cluster_fast $name.corrected.fasta -id 0.$i -msaout top_down/${g}/${i}/fake_${p}_${i}.msa -threads 60");
				open(FILE,"<top_down/${g}/${i}/fake_${p}_${i}.msa");
				@lines=<FILE>;
				close FILE;
				$single=join('',@lines);
				$single=~s/-//g;
				@otus=split(/>\*/,$single);
				shift @otus;
				foreach $otu (@otus){
					($seqs,$cons)=split(/\n>consensus/,$otu,2);
					$seqs=~s/A\n/A/g;
					$seqs=~s/T\n/T/g;
					$seqs=~s/G\n/G/g;
					$seqs=~s/C\n/C/g;
					$seqs=~s/>/\n>/g;
					open(OUTFILE,">>top_down/${g}/${i}/fake_${q}_${i}.fasta");
					print OUTFILE ">$seqs";
					close OUTFILE;
					open(OUTFILE,">>top_down/${g}/${i}/fake_${i}.clusters");
					print OUTFILE ">*$seqs\n";
					close OUTFILE;
					$q++;
				}
				$p++;
			}
		}
	}
	@files=<top_down/${g}/99/*>;
	$p=0;
	$q=0;
	system("mkdir -p -m 777 top_down/${g}/100");
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
			system("/mnt/datA/Applications/usearch7.0.1090_i86linux64 -cluster_fast $name.corrected.fasta -id 1.0 -msaout top_down/${g}/100/fake_${p}_100.msa -threads 60");
			open(FILE,"<top_down/${g}/100/fake_${p}_100.msa");
			@lines=<FILE>;
			close FILE;
			$single=join('',@lines);
			$single=~s/-//g;
			@otus=split(/>\*/,$single);
			shift @otus;
			foreach $otu (@otus){
				($seqs,$cons)=split(/\n>consensus/,$otu,2);
				$seqs=~s/A\n/A/g;
				$seqs=~s/T\n/T/g;
				$seqs=~s/G\n/G/g;
				$seqs=~s/C\n/C/g;
				$seqs=~s/>/\n>/g;
				open(OUTFILE,">>top_down/${g}/100/fake_${q}_100.fasta");
				print OUTFILE ">$seqs";
				close OUTFILE;
				open(OUTFILE,">>top_down/${g}/100/fake_100.clusters");
				print OUTFILE ">*$seqs\n";
				close OUTFILE;
				$q++;
			}
			$p++;
		}
	}
	@orders=@og_orders;
	system("mkdir -p -m 777 top_down/circos_inputs");
	system("mkdir -p -m 777 top_down/circos_inputs/$g");
	open(FILE,"<${g}/grow${g}.fasta");
	@fastas=<FILE>;
	close FILE;
	$fasta=join('',@fastas);
	$population=()=$fasta=~/>/g;
	for($i=75;$i<=100;$i++){
		@files=<top_down/${g}/${i}/*>;
		$start=$shifter*2;
		$written="NULL";
		foreach $file (@files){
			if ($file=~/\.clusters/){
				open(FILE,"<$file");
				@lines=<FILE>;
				close FILE;
				$single=join('',@lines);
				@clusters=split(/>\*/,$single);
				shift @clusters;
				foreach $order (@orders){
					if ($written!~/$order/){
						foreach $cluster (@clusters){
							if ($cluster=~/$order/){
								$Acount=0;
								$Bcount=0;
								$total=0;
								@seqs=split(/>/,$cluster);
								foreach $order2 (@orders){
									foreach $seq (@seqs){
										chomp $seq;
										if ($seq=~/$order2/){
											$written.="$order2";
											push(@neworders,$order2);
											($tag,$code)=split(/\n/,$seq,2);
											foreach $derep (@dereps){
												if ($derep=~/$code/){
													@parts=split(/;/,$derep);
													$size=$parts[1];
													$size=~s/size=//g;
													$total+=$size;
													if ($g>=1){
														@parts=split(/_/,$tag);
														$gc=$parts[3];
														if($gc<=50){
															$Acount+=$size;
														}else{
															$Bcount+=$size;
														}
												
													}
												}
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
								$stop=$start+$total;
								$stop=sprintf("%.0f",$stop);
								open(OUTFILE,">>top_down/circos_inputs/${g}/faux_$i.txt");
								print OUTFILE "faux $start $stop $skew\n";
								close OUTFILE;
								$start=$stop;
							}
						}
					}
				}
			}
		}
		@orders=@neworders;
		undef @neworders;
		print "top_down/${g}/${i} Complete\n";
	}
	print "\n";
	system("mkdir -p -m 777 top_down/confs");
	$start=$start+1;
	$stop=$start+1;
	for($i=75;$i<=100;$i++){
		open(OUTFILE,">>top_down/circos_inputs/${g}/faux_$i.txt");
		print OUTFILE "faux 0 1 $min\n";
		close OUTFILE;
		open(OUTFILE,">>top_down/circos_inputs/${g}/faux_$i.txt");
		print OUTFILE "faux $start $stop $max\n";
		close OUTFILE;
	}
	open(OUTFILE,">>top_down/circos_inputs/${g}/faux_karyotype.txt");
	print OUTFILE "chr - faux faux 0 $stop black";
	close OUTFILE;
	print "top_down/${g} Frame Complete\n";
	$faux=$faux_draft;
	$faux=~s/generation/${g}/g;
	$faux=~s/\/linkage//g;
	$faux=~s/direction/top_down/g;
	open(OUTFILE,">>top_down/confs/faux${g}.conf");
	print OUTFILE "$faux";
	close OUTFILE;
	print "top_down/${g} Configuration Complete\n\n";
}
print "\nTop-Down Complete\n";

# }
# Fake series image creation {

system("mkdir -p -m 777 circos_images");
system("mkdir -p -m 777 circos_images/top_down");
system("mkdir -p -m 777 circos_images/bottom_up");
for($i=0;$i<=$leaps;$i++){
	print "\nCreating Images for Mutation ${i}\n\n";
	foreach $linkage (@linkages){
		system("mkdir -p -m 777 circos_images/bottom_up/${linkage}");
		system("perl /mnt/datA/matt/temp/circos-0.64/bin/circos -conf /mnt/datA/matt/temp/bottom_up/${linkage}/confs/faux${i}.conf");
		print "\nbottom_up/${i}/${linkage} Done\n\n";
	}
	system("perl /mnt/datA/matt/temp/circos-0.64/bin/circos -conf /mnt/datA/matt/temp/top_down/confs/faux${i}.conf");
	print "\ntop_down/${i} Done\n";
}

# }
