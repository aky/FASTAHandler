#!/usr/bin/perl
use strict;
use warnings;

#########################################################
# AMIT KUMAR YADAV (amit.yadav@thsti.res.in)
# Updated 28 April, 2023 Amit
#########################################################
# Moved timestamp and decoy subroutines here
# from Configreader & pepmass modules. Pepmass deleted.

####################################################################
#Date n timestamp for an output file
####################################################################
sub datetimestamp
{
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
	my $datestamp = sprintf"%4d-%02d-%02d_%02d.%02d.%02d",($year+1900),($mon+1),$mday,$hour,$min,$sec;
	return $datestamp;
}

#################################################################
# subroutine to get primary accession number of protein
#################################################################
sub get_id_n_accession {
	my ( $line ) = @_ ;
	chomp $line ;
	# IF-ELSIF LOOPS to capture proteins ids (gi/swissprot/uniprot ids)
	my ( @arr , $id , $desc);
	if (($line=~/^gi/)||($line=~/^decoy_gi/))#this captures GI-id of a protein
        {
                @arr = split (/\|/,$line);
                $id=$arr[0]."|".$arr[1];
		shift @arr;shift @arr;
		$desc=join"|",@arr;
        }
	elsif (($line=~/^sp/)||($line=~/^decoy_sp/))#this captures Swissprot ID
        {
                @arr = split (/\|/,$line);
                $id=$arr[0]."|".$arr[1];
		shift @arr;shift @arr;
		$desc=join"|",@arr;
        }
	# this is for parsing out IPI protein ids
	elsif (($line=~/^IPI/)||($line=~/^decoy_IPI/))
        {
                @arr = split (/\|/,$line);
                $id=$arr[0]."|".$arr[1];
		shift @arr;shift @arr;
		$desc=join "|",@arr;
        }
	elsif($line=~/\|/)
	{
		@arr = split (/\|/,$line);
                $id=$arr[0]."|".$arr[1];
		shift @arr;shift @arr;
		$desc=join "|",@arr;		
	}
	elsif($line=~/Contaminant/)
	{
		if ($line=~/\|/)
		{
			@arr = split (/\|/,$line);
			#check if 2nd el has some value
			if (!defined$arr[1])
			{
				$id=$arr[0];
				shift @arr;
				$desc=join " ",@arr;	
			}
			else
			{
			
				#print scalar(@arr)," el\n";<>;
				#print join"\n",@arr;<>;
				$id=$arr[0]."|".$arr[1];
				shift @arr;shift @arr;
				$desc=join "|",@arr;
			}
		}
		else
		{
			@arr = split (/\s/,$line);
			#check if 2nd el has some value
			if (!defined$arr[1])
			{
				$id=$arr[0];
				shift @arr;
				$desc=join " ",@arr;	
			}
			else
			{
			
				#print scalar(@arr)," el\n";<>;
				#print join"\n",@arr;<>;
				$id=$arr[0]." ".$arr[1];
				shift @arr;shift @arr;
				$desc=join " ",@arr;
			}
		}
	}
	else
        {
                @arr = split (/\||\s+/,$line);
		if (defined $arr[1])
		{
			$id=$arr[0];#."|".$arr[1];
			shift @arr;#shift @arr;
			$desc=join" ",@arr;
		}
		else
		{
			$id=$arr[0];
			shift @arr;
			$desc=join" ",@arr;
		}
	}
	
	#When no description, make id as description too
	if ( (!defined $desc) || ($desc eq '') || ($desc eq ' ') )
	{
		$desc=$id;
	}
	
	#print "$line:--\n...\t$id\n...\t$desc\n\n";<>;
	return ($id,$desc);
}


####################################################################
#Read FASTA file as a hash
####################################################################
sub read_fasta_hash{
	my($fasta)=@_;
	my %db;
	open FASTA,$fasta or die;
	my @fcont=<FASTA>;
	close FASTA;
	my $con=join"", @fcont;
	my@cont=split />/,$con;
	shift @cont; #remove first empty element
	foreach my $seqarr (@cont) #seqarr has header +seq multilines
        {
                my@temp=split "\n",$seqarr;
                my $header=shift @temp;
                my $seq=join"",@temp;
                $seq=~s/\s\r\n//g;
                $db{$header}=$seq;
        }
	#print scalar keys %db," proteins present in DB hash\n";
	return %db;
}

####################################################################
#Index FASTA file (For MassWiz SQLite version)
####################################################################
sub read_fasta_hash_index{
	my($fasta)=@_;
	my %db; #index->seq
	my %index;#index->header
	my $indexedfasta=$fasta;
	$indexedfasta=~s/\.fasta//;
	$indexedfasta.='_index.fasta';
	open FASTA,$fasta or die;
	my @fcont=<FASTA>;
	close FASTA;
	my $con=join"", @fcont;
	my@cont=split />/,$con;
	shift @cont; #remove first empty element
	my $i=1;
	open INDEX,">$indexedfasta" or die;
	foreach my $seqarr (@cont) #seqarr has header +seq multilines
	{
                my@temp=split "\n",$seqarr;
                my $header=shift @temp;
                my $seq=join"",@temp;
                $seq=~s/\s\r\n//g;
                $index{$i}=$header;
                $db{$i}=$seq;
                print INDEX ">$i\n$seq\n";
                $i++;
	}
	close INDEX;
	#print scalar keys %db," proteins present in DB hash\n";
	return (\%db,\%index,$indexedfasta);
}

####################################################################
#Read FASTA file as an Array, Each element is one header+multiline Sequence
####################################################################
sub read_fasta_array{
	my($fasta)=@_;
	my @db;
	open FASTA,$fasta or die;
	@db=<FASTA>;
	close FASTA;
	my $con=join"", @db;
	@db=split />/,$con;
	shift @db;#remove 1st empty element
	print scalar@db," sequences present in FASTA array\n";
	return @db;
}

####################################################################
# Takes in one element of FASTA array, returns Seq and Header
####################################################################
sub FetchSeqID{
	my($seqitem)=@_;
	my @lines=split "\n",$seqitem;
	my $prot_id= get_protein_id($lines[0]);
	shift @lines;
	my $protein= join"",@lines;
	return ($prot_id,$protein);
}
####################################################################
# Get Protein AccID & Description from Accession header line
####################################################################
sub get_accession_description{
	my($header)=@_;
	print "\$header=$header\n";
	if ($header=~m/\|/)
	{
                $header=~m/^(.+\|.+)\|(\w+)$/;
		my $id=$1;
		print "\t\$1=$1\n\t\$2=$2\n\n";<>;
                return($1,$2);
		
	}
	else
	{
                my @arr=split(/\s+/,$header);
                my $acc=shift @arr;
                my $desc=join " ",@arr;
		print "\t\$acc=$acc\n\t\$desc=$desc\n\n";<>;
                return($acc,$desc);
	}
}

sub ppm2da{
        my($val,$tol)=@_;#takes value and tol(in ppm)and returns max and min peak value to search
        my $minval=$val-($tol*$val/1000000);
        my $maxval=$val+($tol*$val/1000000);
        return ($minval,$maxval);
}

################################################################
# Takes a FASTA database
# creates reversed decoy database and returns a 
##################################################################
sub create_decoy_fasta
{
	my($infile,$outfile,$decoytag)=@_;
	if(!defined $decoytag)
	{
		$decoytag='decoy';#default. Edit in COnfigReader
	}
	#$outfile=~s/\.(\w+)$/_decoy.fasta/;
	#print"\nCreating decoy...\n";
	open OUT,">$outfile" or die "Error in reading DECOYDB file name from configuration file.You entered <",($outfile=~ s/\\/DATABASE\//),">\n";
	my %seq = read_fasta_hash($infile);	#reading fasta file
	foreach ( keys%seq )
	{
                chomp $_;
                my $protein= reverse($seq{$_});
                print OUT ">$decoytag"."_".$_,"\n";
                print OUT $protein,"\n";
	}
	close OUT;
	#print"DECOY CREATED\n";
	return $outfile;
}


1;

