#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Statistics::Descriptive;

my $splitgene;

GetOptions("help|h",
	"splitgene|splitgene=s" => \$splitgene,
	);

#die "No vcf file was provided!\n" if !$splitgene;
my $gene_anno_all = 18336;
my $goterm_num_file="GO.term.num";

 
# GO:0006366      457
# GO:0009986      455
# GO:0007155      448


my $gene2go_file="gene2GO";


# A1BG    GO:0003674
# A1BG    GO:0005576
# A1BG    GO:0005615



my $gene_list="protein-coding.gene";
# A1BG
 # A1CF
 # A2M

###########################
# read go term num file
###############################

our %goterm_num;
open IN, "<$goterm_num_file" or die "Open referene genome failed!\n";

my %hash;
while(my $line=<IN>)
{
	chomp($line);
	next if $line !~ /\S+/;	
	my @f = split(/\t/, $line);	
	$goterm_num{$f[0]}=$f[1]; #build hash: [go term]:number
}	
close IN;



#########################
# read gene2go file
###########################

my %gene2go;
open IN2, "<$gene2go_file" or die "Open referene genome failed!\n";
while(my $line=<IN2>)
{
	chomp($line);
	next if $line !~ /\S+/;	
	my @f = split(/\t/, $line);

	if (exists $gene2go{$f[0]})  #exist this gene
	{
		$gene2go{$f[0]}="$gene2go{$f[0]};$f[1]"; #append this new go term to the gene
	}
	else # haven't seen this gene 
	{
		$gene2go{$f[0]}=$f[1]; #make a new entry, [gene]:go
	}
}

close IN2;



########################
# loop through each gene pair
############################

open IN3, "<$gene_list" or die "Open referene genome failed!\n";
open IN4, "<$gene_list" or die "Open referene genome failed!\n";
while(my $line=<IN3>) # gene list
{
	chomp($line);
	next if $line !~ /\S+/;

	print "$line";

	seek IN4,0,0;  
	while(my $line2=<IN4>) # gene list
	{
		chomp($line2);
		next if $line2 !~ /\S+/;

		my $overlap=";0"; 
		my $e_go;
		
		
        #[gene]:go
		if(exists $gene2go{$line} && exists $gene2go{$line2}) # there are GO terms associated with gene 1 or gene 2
		{
        # call function overlap_goterm
			$overlap=&overlap_goterm($gene2go{$line},$gene2go{$line2});
		}

		if($overlap eq ";0")
		{
			$e_go=0;
		}
		else
		{
			$e_go=&e_go_analysis($overlap);
		} 
		print "\t$e_go";
	}

	print "\n";


}

close IN3;
close IN4;

# functions: sub

sub overlap_goterm()
{
	my ($go1, $go2) = @_; # list

	my $overlap_go="";

	my @go1_split=split(/\;/,$go1);
	my @go2_split=split(/\;/,$go2);

	for (my $i = 0; $i < @go1_split; $i++) 
	{
		for (my $j = 0; $j < @go2_split; $j++) 
		{
			if ($go1_split[$i] eq $go2_split[$j]) 
			{
				$overlap_go=$overlap_go.";".$go1_split[$i] #append overlap terms
			}
		}
	}

	if (length($overlap_go)==0) 
	{
		return ";0";
	}
	else
	{
		return $overlap_go;
	}
}



sub e_go_analysis()
{
	(my $overlap) = @_;
	my @go_term =split(/\;/,$overlap);

	{
		my $sum=0;
		for (my $i = 1; $i < @go_term; $i++) 
		{
			
			#print $goterm_num{$go_term[$i]}.":";
			$sum += -2*log($goterm_num{$go_term[$i]}/$gene_anno_all)/log(10);
		}
		return $sum;
	}
		
}
=a

foreach my $line1 (keys %hash)	  
{	
	
	#if($k>=$splitgene && $k<$splitgene+2000)
	{
		my @x=split("\t",$hash{$line1});
		foreach my $line2 (keys %hash)	  
		{	
			my @y=split("\t",$hash{$line2});
			my $pearson=&pearson_mutlib(\@x,\@y);
			print TO1 $pearson."\t"; 
		}
		print TO1 "\n";
		print TO2 $line1."\n";
	}
	$k++;
}

close TO1;
close TO2;


sub method()
{
	my @temp=@_;
	my $stat = Statistics::Descriptive::Full->new();
	$stat->add_data(\@temp);
	my $mean = $stat->mean();
	#my $standard_deviation=$stat->standard_deviation();
	my $variance=$stat->variance();
	$variance=$variance*(@temp-1)/@temp; ##
	return ($mean,$variance);
}

sub sum_xy()
{
	my ($ref_a, $ref_b) = @_;
	my @x = @{$ref_a};
	my @y = @{$ref_b};
	my $sum_xy = 0;

	for(my $i=0;$i<@x;$i++)
	{
		$sum_xy += $x[$i]*$y[$i];
	}

	return($sum_xy);
		
}

sub pearson()
{
	my ($ref_a, $ref_b) = @_;
	my @x = @{$ref_a};
	my @y = @{$ref_b};

	(my $x_mean, my $x_var) = &method(@x);
	(my $y_mean, my $y_var) = &method(@y);

	my $sum_xy = &sum_xy(\@x,\@y);

	my $cor=($sum_xy/@x-$x_mean*$y_mean);
	my $var_x_y=$x_var*$y_var;
	my $pearson_cor=$cor/sqrt($var_x_y);

	return($pearson_cor);
}

sub pearson_mutlib()
{
	my ($ref_a, $ref_b) = @_;
	my @x = @{$ref_a};
	my @y = @{$ref_b};


	(my $x_mean, my $x_var) = @x[0..1];
	(my $y_mean, my $y_var) = @y[0..1];

	my @x_other=@x[2..@x-1];
	my @y_other=@y[2..@y-1];

	my $sum_xy = &sum_xy(\@x_other,\@y_other);

	my $cor=($sum_xy/@x_other-$x_mean*$y_mean);
	my $var_x_y=$x_var*$y_var;
	my $pearson_cor=$cor/sqrt($var_x_y);

	return($pearson_cor);

}



