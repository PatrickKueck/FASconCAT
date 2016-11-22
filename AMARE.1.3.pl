#!/usr/bin/perl

#####################################################################################################################################################################

# written by Patrick Kück at August		2009 ;
# updated by Patrick Kück at August		2010 ;
# updated by Patrick Kück at September	2010 ; changed 18.10.10 ; changed 27.10.10 ; ; changed 19.12.10
# updated by Patrick Kück at May		2011 ;
# updated by Patrick Kück at June		2011 ;
# updated by Patrick Kück at February	2011 ;
# updated by Patrick Kück at september	2012 ; inclusion of the '-r' option to remain complete invariable 0|0 bins
# debugged by Patrick Kück at september 2012 ; sub readin_infile -> sub &file_check : last 8 single taxa not considered in the analyses if infile was in Beckman format
# debugged by Patrick Kück at september 2012 ; bug in abi-concatenation and print terminal ;


#####################################################################################################################################################################
# AMARE-Header and Global Variables
# ---------------------------------

use strict         ;
use File::Copy     ;
use Tie::File      ;
use Getopt::Std    ;

print  "\n\t---------------------------------------------" ;
print  "\n\t          Welcome to AMARE !     "             ;
print  "\n\t   A Software for AFLPmatrix reduction"        ;
print  "\n\twritten by Patrick Kueck (ZFMK Bonn, 2010/11)" ;
print  "\n\t---------------------------------------------" ;

my @infiles					= () ;                               # List of all given inputfile names (defined in subroutine 'ARGV')
my @input_lines				= () ;                               # Collects all input lines of given data (defined in section MAIN ANALYSES)
my @double_samples			= () ;                               # Collects all double sample BIN strings (defined in section MAIN ANALYSES)
my %input_of_file			= () ;                               # Key: Filename; Value: List of corresponding input lines (defined in subroutine 'readin_file')
my %missing_taxa_of_file	= () ;                               # Key: Filename; Value: Taxonnames which are found in other input files, but which are not included in the corresponding input file (defined in subroutine 'readin_file')
my $N_files_in			= () ;                               # Number of given inputfiles (defined in subroutine 'readin_infile')
my $datatype				= () ;                               # Datatype of inputfiles (ABI or BECKMAN) (defined in subroutine 'ARGV')
my $level_eq_neq			= () ;                               # BIN Threshold of repeatability of individual loci (defined in subroutine 'ARGV')
my $level_11_01			= 0  ;                               # BIN Threshold of repeatability of equal 11 bins in ratio to unequal bins (predefined as zero; not used yet)
my $level_distance		= () ;                               # Threshold for minimum allowed BIN distances (defined in subroutine 'ARGV')
my $error_type			= () ;                               # Type of chosen Error Rate (Euclidean or Jaccard) (defined in subroutine 'ARGV')
my $handling_bins			= 1  ;
my $infilename_main_log	= 'AMARE_supermatrix.txt' ;          # Name of the AMARE main log output file (predefined)



#####################################################################################################################################################################
# PRE-ANALYSES: ARGV, Input file read IN and check, concatenation of multiple infiles, and converting of ABI files to BECKMAN format
#####################################################################################################################################################################

#####################################
## ARGV read-in
&argv( \@ARGV , \@infiles , \$level_eq_neq , \$level_distance , \$datatype , \$error_type, \$handling_bins ) ;
#####################################

#####################################
## Read IN and check of single input files
&readin_infile( \@infiles , \%input_of_file , \$N_files_in , \$datatype , \%missing_taxa_of_file ) ;
#####################################

#####################################
## Generate AMARE subfolder
mkdir "AMARE_Matrices" , 0755 ;
mkdir "AMARE_logfiles" , 0755 ;
mkdir "AMARE_svgfiles" , 0755 ;
#####################################

#####################################
## Data converting of ABI to BECKMAN format and concatenation process for more than one input file
if ($N_files_in >= 2 ){
	
	if ( $datatype =~ /b/i ){										@input_lines = &concatenation( \@infiles , \%input_of_file , \%missing_taxa_of_file ) }
	if ( $datatype =~ /a/i ){	&convert_abi ( \%input_of_file ) ;	@input_lines = &concatenation( \@infiles , \%input_of_file , \%missing_taxa_of_file ) }
}
else { $infilename_main_log = $infiles[0] ;
	
	if ( $datatype =~ /a/i ){ &convert_abi ( \%input_of_file ) }
	@input_lines = exists ( $input_of_file{$infiles[0]} ) ? @{$input_of_file{$infiles[0]}} : ()
}
#####################################



#####################################################################################################################################################################
# MAIN ANALYSES: 
#####################################################################################################################################################################

my %results_of_taxcut_level	= () ; # Key: Infostring; Value: all deleted BIN positions and Replicate samples (defined in subroutine 'error_rate')
my %graph_line_of_tax		= () ; # Key: Replicate Sample Names; Value: consensus character string of each Replicate Sample for graphical output (defined in subroutine 'graphical_bin_states')

#####################################
## Check of equal fragment and sample length values
## BIN columns with unequal values are collected in @delete_columns
my @delete_columns = &fragment_sample_check( \$input_lines[2] , \$input_lines[3] ) ;
#####################################

#####################################
## Collect all double samples in @double_samples
for my $line ( @input_lines ){ if ( $line =~ /^a\w+,.*$/ ){ push @double_samples , $line } }
@double_samples = sort @double_samples ;
#####################################

############################################################################################################### Main BIN check loop
## For each BIN threshold perform data analyses with different replicate thresholds
for ( my $level_eq_ne = $level_eq_neq ; $level_eq_ne <= 95 ; $level_eq_ne++ ){ print "\n\n\tLevel Eq <=> Uneq: $level_eq_ne %" ;
	
	#TAXCUT:
	for ( my $level_tax_eq = 0 ; $level_tax_eq < 100 ; $level_tax_eq += 10 ){ 
		
		my	%deleted_tax		= () ;
			$level_eq_ne		= sprintf "%.1f", $level_eq_ne  ;
			$level_tax_eq		= sprintf "%.1f", $level_tax_eq ; 
		my	$N_del_samp_frag	= @delete_columns ;
		my	$outfile			= "AMARE_".$level_eq_ne."_".$level_tax_eq."_log.txt"   ;
		
		#####################################
		## Open OUT logfile
		open	LOG,	">AMARE_logfiles/$outfile" or warn "\n\t\nERROR: Can not print out /AMARE_logfiles/$outfile!" ;
		print	LOG		"AMARE LOG-FILE\n\nParameter-Setup\n----------------\n\tLevel 1|1+0|0/0|1 bin-ratio:\t$level_eq_ne\n\tLevel 1|1+0|0/0|1 taxon-ratio:\t$level_tax_eq\n\tLevel Bin-Distance:\t\t$level_distance\n\n" ;
		print	LOG		"Sample/fragment-Check\n---------------------\n\tNumber of deleted bins: $N_del_samp_frag\n\tDeleted bins: @delete_columns\n\n" ;
		print			"\n\tbin_taxon_check: level taxon restriction: $level_tax_eq" ;
		#####################################
		
		#####################################
		## Main Loop
		&loop( \@delete_columns , \@double_samples , \%deleted_tax , \$level_11_01 , \$level_eq_ne , \$level_tax_eq , \$level_distance , \$input_lines[4] , \$input_lines[5] , \%results_of_taxcut_level , \$error_type, \$handling_bins ) ;
		#####################################
		
		############################################################################################################### BIN-Check for graphical output after analyses
		## Compare BIN states between single Sample Replicates for graphical visualization after analysing under
		## define BIN and Replicate Threshold
		&graphical_bin_states(\@double_samples, \%graph_line_of_tax, \$level_eq_ne) ;
		#####################################
		
		close LOG ;
	}
}
#####################################

# Analyse Bin-Check results
&taxcut_analyse( \%results_of_taxcut_level , \@input_lines , \$infilename_main_log , \%graph_line_of_tax , \$error_type ) ;


##

#####################################################################################################################################################################
# END
#####################################################################################################################################################################


#####################################################################################################################################################################
# SUBROUTINES PRE-ANALYSES
#####################################################################################################################################################################

sub argv{
	
	#####################################################
	# Subroutine ARGV handles and checks all given command line options
	# If command line options are given incorrectly, AMARE finished with a synopsis description
	# If -help is typed via command line, AMARE finished with a more detailed help description
	#####################################################
	
	my $aref_argv				= $_[0] ; # User given command line options
	my $aref_infiles			= $_[1] ; # List of all given inputfile names (defined in subroutine 'ARGV')
	my $sref_level_t			= $_[2] ; # Threshold of repeatability of individual loci (defined in subroutine 'ARGV')
	my $sref_level_d			= $_[3] ; # Threshold for minimum allowed BIN distances (defined in subroutine 'ARGV')
	my $sref_datatype			= $_[4] ; # Datatype of inputfiles (ABI or BECKMAN) (defined in subroutine 'ARGV')
	my $sref_errortype		= $_[5] ; # Type of chosen Error Rate (Euclidean or Jaccard) (defined in subroutine 'ARGV')
	my $sref_handling_bins	= $_[6] ; # Additional option which allows to remain complete 00 bin states
	
	#####################################
	## This allows alternative input commands: not only '-command1', but also '- command 1', '- command 1-  command2', and so on...
	my ( $commandline ) = join "", @$aref_argv ;
	#####################################
	
	############################################################################################################### ARGV check and asignment
	## Assignment and check of single given commands
	if ( $commandline ne undef ){
		
		my @commands = split "-", $commandline ; shift @commands ;
		my %commands = () ;
		
		REPEAT_ARGV:
		for my $single_command ( sort @commands ){
			
			my @com_parts = split "" , $single_command ;
			my $first_com_part = shift @com_parts ; my $parameter = join "" , @com_parts ; $commands{$first_com_part}++ ;
			
			$first_com_part =~ /^h$/ and do { if ( $parameter =~ /elp/      ){ &help ; &synopsis ; exit										} else { print "\n\n\tCOMMAND-ERROR: Unknown command $parameter!\n"								; &synopsis ; exit } } ;
			$first_com_part =~ /^i$/ and do {	 push @$aref_infiles,											  $parameter ; next REPEAT_ARGV	} ;
			$first_com_part =~ /^d$/ and do { if ( $parameter =~ /^\d+\.\d+$/ ){ $$sref_level_d			= $parameter ; next REPEAT_ARGV	} else { print "\n\n\tCOMMAND-ERROR: Wrong Minimum BIN Distance Threshold: $parameter!\n"			; &synopsis ; exit } } ;
			$first_com_part =~ /^l$/ and do { if ( $parameter =~ /^\d+\.\d+$/ ){ $$sref_level_t			= $parameter ; next REPEAT_ARGV	} else { print "\n\n\tCOMMAND-ERROR: Wrong Minimum BIN Reliability Threshold: $parameter!\n"		; &synopsis ; exit } } ;
			$first_com_part =~ /^t$/ and do { if ( $parameter =~ /^b$|^a$/    ){ $$sref_datatype			= $parameter ; next REPEAT_ARGV	} else { print "\n\n\tCOMMAND-ERROR: Unknown type of specified table format: $parameter!\n"		; &synopsis ; exit } } ;
			$first_com_part =~ /^e$/ and do { if ( $parameter =~ /^e$|^j$/    ){ $$sref_errortype			= $parameter ; next REPEAT_ARGV	} else { print "\n\n\tCOMMAND-ERROR: Unknown type of specified ERROR RATE: $parameter!\n"			; &synopsis ; exit } } ;
			$first_com_part =~ /^r$/ and do { if ( $parameter =~ /^0$|^1$/    ){ $$sref_handling_bins	= $parameter ; next REPEAT_ARGV	} else { print "\n\n\tCOMMAND-ERROR: Unknown type of specified BIN state handling: $parameter!\n"	; &synopsis ; exit } } ;
			
			#####################################
			## If given command sign can not be found, finish AMARE with an error prompt and the synopsis
			print "\n\n\tCOMMAND-ERROR: Unknown option: $single_command!" ; &synopsis ; exit
			#####################################
		}
		
		#####################################
		## Command options  d, l, t, e have to be given once. If commands are not found or found multiple times, finish AMARE with an error prompt and the synopsis
		my @options = (qw/d l t e/) ;
		for ( @options ){ 
			
			if ( $commands{$_} > 1 ){ print "\n\n\tCOMMAND-ERROR: Multiple usage of option: -$_  !\n\n"    ; &synopsis ; exit } 
			if ( $commands{$_} < 1 ){ print "\n\n\tCOMMAND-ERROR: Missing parameter for option: -$_ !\n\n" ; &synopsis ; exit }
		}
		#####################################
	}
	
	else{ &synopsis ; exit }
	############################################################################################################### Synopsis and Help Menu
	
	sub synopsis{
		
		print  "\n\n\tSYNOPSIS\n\t--------------------------------------------" ;
		print  "\n\tStart (UNIX):\n\tperl AMARE_1.215.pl -i <I1> -l <T1> -d <T2> -t <DT> -e <ER> -r <0>\n" ;
		print  "\n\tStart (WINDOWS):\n\tAMARE_1.3beta.pl -i <I1> -l <T1> -d <T2> -t <DT> -e <ER> -r <0>\n" ;
		print  "\n\t<I1>: Inputfile 1 in .txt format" ;
		print  "\n\t<T1>: % Threshold of Minimum Bin Reliability <float>" ;
		print  "\n\t<T2>: % Threshold oof Minimum BIN-distance <float>" ;
		print  "\n\t<DT>: Type of data:\n\t\t<a> for GeneMapper table format\n\t\t<b> for BECKMAN table format\n\t\tA combination of both file formats is not possible\n" ;
		print  "\n\t<ER>: Type of error-rate calculation\n\t\t<e> for BONIN error-rate calculation\n\t\t<j> for JACCARD distance calculation\n\t\tA combination of both calculations is not possible\n" ; 
		print  "\n\t<0> : To  disable the default exclusion of\n\t\tcomplete (0,0) BIN states type -r 0\n" ;
		print  "\n\n\tExample 1 (1 GeneMapper infile under WINDOWS)\n\t\tAMARE_1.3beta.pl -i i1.txt -l 90.0 -d 0.25 -t a -e e" ;
		print  "\n\n\tExample 2 (2 BECKMAN infiles under UNIX)\n\t\tperl AMARE_1.3beta.pl -i i1.txt -i i2.txt -l 80.0 -d 0.15 -t b -e j\n\n" ;
		print  "\n\tHelp: (perl) AMARE_1.3beta.pl -help\n\t" ;
		print  "--------------------------------------------\n\n" ;
	}
	
	
	sub help{

	print
<<help

        ------------AMARE HELP---------------------
	AMARE is designed to optimize the signal-to-noise
	ratio in AFLP data sets (single or multiple primer
	combinations). It identifies potentially erroneous
	AFLP genotyping errors of GeneMapper (Applied Bio-
	systems) or CEQTM System Fragment Analysis v.9.0.25
	(Beckman Coulter) as the number of unreproducibly
	scored markers between replicated AFLP profiles of
	single individuals. AMARE concatenates multiple
	AFLP matrices and then analyzes the concatenated
	supermatrix in one process.
	
	For each primer combination /AFLP matrix type:
		'-i inputfile_I.txt' '-i inputfile_II.txt'
		
	To specify the input data type:
		'-t a' for GeneMapper table format
		'-t b' for BECKMAN table format
		
	Do not mix GeneMapper table format with BECKMAN 
	table format! AMARE cannot analyze both table
	formarts in one single run!
	
	AMARE considers 3 different kinds of thresholds. 
	Threshold 1 and 3 are specified by the user:
	
	1) The Minimum BIN Reliability Threshold (-l) 
	defines the relative number of reproducible BINs
	over all replicates which can range between 0 and
	1. The user specified BIN Reliability Threshold
	sets the acceptance value of the minimum number
	of reproducible (0,0) and (1,1) BIN states. 
	To define the Minimum BIN Reliability Threshold type:
		'-l float' (e.g. '-l 88.0)'
	2) The Replicate Reliability Threshold (automatically
	incorporated) is defined as the relative number of
	reproducible markers between replicates y of a single
	individual over all n bins
	3) The Minimum Bin distance threshold (-d) determines
	the minimum acceptable BIN distance between differ-
	ently sized BINs. To specify the Minimum Bin distance
	threshold type: 
		'-d float' (e.g. '-d 0.25')
	
	AMARE calculates 2 error-rates:
	1) BONIN distance (-e e):
		(Incongruent_BIN_states_01 + Incongruent_BIN_states_10) 
		/ (Incongruent_BIN_states_01 + Incongruent_BIN_states_10
		+ Congruent_BIN_states_11 + Congruent_BIN_states_00 )
	2) JACCARD distance (-e j):
		(Incongruent_BIN_states_01 + Incongruent_BIN_states_10) 
		/ (Incongruent_BIN_states_01 + Incongruent_BIN_states_10
		+ Congruent_BIN_states_11 )
	
	AMARE excludes all BINs which consist exclusively of 0|0 states
	per default. To turn this behaviour off type:
		'-r 0'
	
	For each individual threshold set, AMARE generates
	i) a single log file reporting the masking of bins
	and replicates and ii) a character matrix in text (.txt)
	and nexus (.nex) format, if error rate conditions
	and minimum number of remaining bins are met. A
	summary of all threshold sets and corresponding
	error rates are stored in the main log file.
	Additionally, AMARE plots a graphical overview of
	the original and each masked replicate matrix. 
	
	Further information about the algorithm, usage, 
	and options of AMARE are found in the AMARE manual.
	
	For further questions or bug reports, please write to
	patrick_kueck\@web.de
	
help
;

	}
}


sub readin_infile{
	
	#####################################################
	# Subroutine 'readin_file' reads IN all input files, Untie line feeds, and verifies correct file formats
	# Each inputfile is checked for:
	# - Correct line format, equal number of columns, double sample names, and not allowed equal sample names -> subroutine 'file_check'
	# - Defined and paired number of replicate samples -> subroutine 'readin_file'
	# - Missing taxa in single input files for concatenation process -> subroutine 'readin_file'
	#####################################################
	
	print "\n\tREAD-IN and File-Check..." ;
	
	my $aref_infiles					= $_[0]	; # List of all given inputfile names (defined in subroutine 'ARGV')
	my $href_in_of_file				= $_[1]	; # Key: Filename; Value: List of corresponding input lines (defined in subroutine 'readin_file')
	my $sref_N_files_in				= $_[2]	; # Number of inputfiles, needed to find possible missing taxa in single input files (defined in subroutine 'readin_infile')
	my $sref_type_data				= $_[3]	; # Datatype of inputfiles (ABI or BECKMAN) (defined in subroutine 'ARGV')
	my $href_missing_taxa_of_file	= $_[4]	; # Key: Filename; Value: Taxonnames which are found in other input files, but which are not included in the corresponding input file (defined in subroutine 'readin_file')
	my %N_names_all_files				= ()	; # Key: Sample Name; Value: List of Filenames which include corresponding sample name (defined in subroutine 'file_check')
	my %replicate_names					= ()	; # Key: Double-Sample Name; Value: Number of occurence in inputfiles (defined in subroutine 'file_check')
	
	############################################################################################################### Input Files
	## Untie line feeds, Read In of inputfiles, Check right format of each inputfile, 
	## store inputlines of correct inputfiles as hashreference, and count correct input files
	for my $file ( @$aref_infiles ){
		
		#####################################
		## Untie line feeds
		&tie( \$file ) ;
		#####################################
		
		#####################################
		## Read IN input lines of corresponding input file
		open IN, "<$file" or die "\n\n\tFILE-ERROR: Can not open $file !\n" ;
		chomp ( my @input_data = <IN> ) ;
		#####################################
		
		#####################################
		## Check right format of each inputfile
		&file_check( \@input_data , \$$sref_type_data, \$file , \%N_names_all_files, \%replicate_names ) ;
		#####################################
		
		#####################################
		## Store inputlines of correct inputfiles as hashreference
		push @{$href_in_of_file->{$file}} , @input_data ;
		#####################################
		
		#####################################
		## Count correct input files
		$$sref_N_files_in++
		#####################################
	}
	#####################################
	
	############################################################################################################### Replicate Samples
	## Check for defined and paired number of replicates
	my $N_taxanames =      keys %replicate_names ;
	my @replicates  = sort keys %replicate_names ; # Sorting needed for assignment of corresponding replicate names
	
	if ( $N_taxanames     == 0 ){ die "\n\n\tFILE-ERROR: Cannot find predefined replicates !\n\n"     }
	if ( $N_taxanames % 2 == 1 ){ die "\n\n\tFILE-ERROR: Impair number of predefined replicates!\n\n" }
	#####################################
	
	#####################################
	## Check for identical replicate names
	my @wrong_re_names = () ;
	for ( my $k=0; $k<=$N_taxanames; $k+=2 ){ my $l = $k+1 ;
		
		#####################################
		## Reject optional quotation marks in sample names
		$replicates[$k] =~ s/"//g ;
		$replicates[$l] =~ s/"//g ;
		#####################################
		
		#####################################
		## Remove the last sign (difference between equal sample replicates)
		my @rep_1 = split "", $replicates[$k] ; pop @rep_1 ; my $rep_name_1 = join "", @rep_1 ;
		my @rep_2 = split "", $replicates[$l] ; pop @rep_2 ; my $rep_name_2 = join "", @rep_2 ;
		#####################################
		
		#####################################
		## Store non-identic replicate names in list '@wrong_re_names'
		unless ( $rep_name_1 =~ /^$rep_name_2$/i ){ push @wrong_re_names, $replicates[$k]." and ".$replicates[$l] }
		#####################################
	}
	
	#####################################
	## If non-identic replicate names are found print OUT list of incorrect replicate names and finish AMARE with an error prompt
	if ( @wrong_re_names ){ for ( @wrong_re_names ){ print "\n\n\tFILE-ERROR: Unpaired replicates:\n\t$_\n" } print "\n\n\tExcept of the last sign, paired replicates have to be named identic!\n\te.g. aHomo_sapiens_1 and aHomo_sapiens_2\n\n"; exit }
	#####################################
	
	############################################################################################################### Missing Taxa
	## If multiple infiles are given, missing taxa in single input files have to be identified and stored as hashlist '$href_missing_taxa_of_file'
	for my $taxon ( keys %N_names_all_files ){ 
		
		#####################################
		## Store all filenames in which the corresponding taxon is included in a list
		my $aref_files_involved = exists ( $N_names_all_files{$taxon} ) ? \@{$N_names_all_files{$taxon}} : () ;
		#####################################
		
		#####################################
		## Compare number of inputfiles in which corresponding taxon occure with total number of input files
		## If number of inputfiles is not equal to the total number of input files, store corresponding taxon in hashlists of missing input files 
		my %seen_file = () ;
		unless ( @$aref_files_involved == $$sref_N_files_in ){ 
			
			for my $file ( @$aref_files_involved	){				$seen_file{$file}++	}
			for my $file ( keys %$href_in_of_file	){ unless (	$seen_file{$file}		){ push @{$href_missing_taxa_of_file->{$file}} , $taxon; print "\n\tTaxon $taxon not included in $file" } print "\n"; }
		}
		#####################################
	}
	#####################################
	
	############################################################################################################### Subroutines within 'sub readin_infile'
	sub tie{
		
		#####################################################
		# Subroutine 'tie' adapts all linefeeds to unix-format (\n)
		
		my $tie_file = $_[0] ;
		
		TIE:
		(tie ( my @data, 'Tie::File', $$tie_file )) ;
		
		die  "\n\t!FILE-ERROR!: $infiles[0] is empty!\n" if 0 == @data ;
		
		map { s/\r\n/\n/g } @data ;
		map { s/\r/\n/g   } @data ;
		
		untie @data ;
		#####################################################
	}
	
	sub file_check{
		
		#####################################################
		# Subroutine 'file_check' checks each each line ($aref_inlines) of given input files ($name_file) according to given datatype ($sref_d_type) for 
		# Correct line format, equal number of columns, double sample names, and equal sample names
		# Hashreferenz $href_Nname_all stores all filenames in which given sample names occure (important for concatenation process)
		# Hashreferenz $href_replicate stores double sample names
		#####################################################
		
		my $aref_inlines		= $_[0] ; # List of inputlines
		my $sref_d_type		= $_[1] ; # Given datatype of inputfile (BECKMAN or ABI)
		my $name_file			= $_[2] ; # Name of the inputfile
		my $href_Nname_all	= $_[3] ; # Key: Sample Name; Value: List of Filenames which include corresponding sample name (defined in subroutine 'file_check')
		my $href_replicate	= $_[4] ; # Key: Double-Sample Name; Value: Number of occurence in inputfiles (defined in subroutine 'file_check')
		
		############################################################################################################### ABI
		## If datatype = 'ABI'
		if ( $$sref_d_type =~ /a/ ){
			
			#####################################
			## Check line 1,2 and 3 for correct line format
			unless ( $aref_inlines->[0] =~ /^\s+(\d+\s)+(\d+)?$/ig                  ){ die "\n\n\tFILE-ERROR: Line 1 of table $$name_file is not in GeneMapper table format!\n\n" }
			unless ( $aref_inlines->[1] =~ /^Sample.Name\s("?\w+"?\s)+("?\w+"?)?$/i ){ die "\n\n\tFILE-ERROR: Line 2 of table $$name_file is not in GeneMapper table format!\n\n" }
			unless ( $aref_inlines->[2] =~ /Dye/i                                   ){ die "\n\n\tFILE-ERROR: Line 3 of table $$name_file is not in GeneMapper table format!\n\n" }
			#####################################
			
			
			for ( @$aref_inlines ){ $_ =~ s/"|'+//g }
			
			#####################################
			## split line 0 (BIN numbers) and line 1 (Sample names)
			my @N_cols_line_0  = split "\t" , $aref_inlines->[0] ;
			my @name_line      = split "\t" , $aref_inlines->[1] ; 
			#####################################
			
			#####################################
			## Check Bin lines for correct line format and column numbers
			for ( my $line=3; $line<=@$aref_inlines-1 ; $line++ ){ my $real_line = $line+1 ;
				
				#####################################
				## Check of column numbers
				( my $line_x = $aref_inlines->[$line] ) =~ s/\t/\tpk/g ;
				my @N_cols  = split "\t" , $line_x ; my $Ncols = @N_cols ;
				unless ( @N_cols == @N_cols_line_0 ){ die "\n\tFILE-ERROR: Unequal number of columns in line $real_line\n\tof GeneMapper table $$name_file\n\n" }
				#####################################
				
				#####################################
				## Check of line format
				unless ( $aref_inlines->[$line] =~ /^Allele \d+(\s+\d+(\(\d\))?|\s+\?)+(\s+)?$|^Size \d+(\s+\d+(.\d+)?)+(\s+)?$/ig ){ die "\n\tFILE_ERROR: Forbidden sign in line $real_line of GeneMapper table $$name_file!\n\n" }
				#####################################
			}
			#####################################
			
			#####################################
			## Check sample names and assign double sample names
			my %tax_names = () ;
			for my $name ( @name_line ){
				
				#####################################
				## If sample name starts with 'a' or '"a"' define it as double sample name
				unless ( $tax_names{$name} ){
					
					if    ( $name =~ /^a\w+$/ ){ $href_replicate->{$name}++ }
					push @{$href_Nname_all->{$name}} , $$name_file ; $tax_names{$name}++
				}
				else   { die "\n\n\tFILE-ERROR: Sample $name appears more than once in table $$name_file !\n\n" }
				#####################################
			}
			#####################################
		}
		#####################################
		
		############################################################################################################### BECKMAN
		## If datatype = 'BECKMAN'
		if ( $$sref_d_type =~ /b/ ){
			
			#####################################
			## Check infolines 0-8 for right format
			unless ( $aref_inlines->[0] =~ /^"?Bin"?(,"?\d+"?)+$/ig          ){ die "\n\tFILE_ERROR: Wrong format of Bin-line in Beckman Coulter table $$name_file!\n\n"       }
			unless ( $aref_inlines->[1] =~ /^"?Dye"?(,"?\w+"?)+$/ig          ){ die "\n\tFILE_ERROR: Wrong format of Dye-line in Beckman Coulter table $$name_file!\n\n"       }
			unless ( $aref_inlines->[2] =~ /^"?Samples"?(,"?\d+"?)+$/ig      ){ die "\n\tFILE_ERROR: Wrong format of Samples-line in Beckman Coulter table $$name_file!\n\n"   }
			unless ( $aref_inlines->[3] =~ /^"?Fragments"?(,"?\d+"?)+$/ig    ){ die "\n\tFILE_ERROR: Wrong format of Fragments-line in Beckman Coulter table $$name_file!\n\n" }
			unless ( $aref_inlines->[4] =~ /^"?XMin"?(,"?\d+(.\d+)?"?)+$/ig  ){ die "\n\tFILE_ERROR: Wrong format of XMin-line in Beckman Coulter table $$name_file!\n\n"      }
			unless ( $aref_inlines->[5] =~ /^"?XMax"?(,"?\d+(.\d+)?"?)+$/ig  ){ die "\n\tFILE_ERROR: Wrong format of Xmax-line in Beckman Coulter table $$name_file!\n\n"      }
			unless ( $aref_inlines->[6] =~ /^"?XMean"?(,"?\d+(.\d+)?"?)+$/ig ){ die "\n\tFILE_ERROR: Wrong format of XMean-line in Beckman Coulter table $$name_file!\n\n"     }
			unless ( $aref_inlines->[7] =~ /^"?XVar"?(,"?\d+(.\d+)?"?)+$/ig  ){ die "\n\tFILE_ERROR: Wrong format of XVar-line in Beckman Coulter table $$name_file!\n\n"      }
			unless ( $aref_inlines->[8] =~ /^"?YMean"?(,"?\d+(.\d+)?"?)+$/ig ){ die "\n\tFILE_ERROR: Wrong format of YMean-line in Beckman Coulter table $$name_file!\n\n"     }
			#####################################
			
			#####################################
			## split line 0 (BIN numbers)
			my @N_cols_line_0  = split "," , $aref_inlines->[0] ;
			#####################################
			
			#####################################
			## Check infolines 0-8 for equal number of columns
			for my $linenumber ( 0 .. 8 ){ my $real_line = $linenumber+1 ; 
				
				my @N_cols = split ",", $aref_inlines->[$linenumber] ;
				unless ( (@N_cols) == (@N_cols_line_0) ){ die "\n\tFILE_ERROR: Unequal number of columns in line $real_line\n\tof Beckman Coulter table $$name_file!\n\n" }
			}
			#####################################
			
			#####################################
			## Check sample lines for correct number of columns, format, sample names, and assign double sample names
			my %tax_names = () ;
			for ( my $line=9; $line<=@$aref_inlines-1 ; $line++ ){ my $real_line = $line+1 ;
				
				#####################################
				## Check sample lines for right format
				if ( $aref_inlines->[$line] =~ /^"?\w+"?(,\d|,\?)+$/ig ){
					
					#####################################
					## Check of equal sample names
					$aref_inlines->[$line] =~ s/"|'+//g ;
					my @parts = split "," , $aref_inlines->[$line] ;
					
					unless ( $tax_names{$parts[0]} ){
						
						#####################################
						## If sample name starts with 'a' or '"a"' define it as double sample name
						if   ( $parts[0] =~ /^a\w+$/ ){ $href_replicate->{$parts[0]}++ }
						#####################################
						
						$tax_names{$parts[0]}++ ; push @{$href_Nname_all->{$parts[0]}} , $$name_file ;
						
						#####################################
						## Check of correct column number
						my @N_cols = split ",", $aref_inlines->[$line] ;
						unless ( (@N_cols) == (@N_cols_line_0) ){ die "\n\tFILE_ERROR: Unequal number of columns in line $real_line\n\tin $$name_file!\n\n" }
						#for ( 1 .. @N_cols-1 ){ unless ( $N_cols[$_] =~ /^1$|^0$|^\?$/ ){ die "\n\tFILE_ERROR: Unallowed character state in line $real_line\n\tin $$name_file!\n\n" } }
						#####################################
					}
					else { die "\n\tFILE-ERROR: Sample $parts[0] appears more than once in table $$name_file!\n\n" }
					#####################################
				}
				else{ die "\n\tFILE-ERROR: Wrong format of of Bin-line $real_line in Beckman Coulter table $$name_file!\n\n" }
				#####################################
			}
			#####################################
		}
	}
}


sub convert_abi{
	
	#####################################################
	# Subroutine 'convert_abi' converts the format of ABI input files into BECKMAN format
	#####################################################
	
	my $href_in_of_file = $_[0] ; # Key: Filename; Value: List of corresponding input lines (defined in subroutine 'readin_file')
	
	for my $file ( keys %$href_in_of_file ){
		
		print "\n\tConvert GeneMapper-file $file..." ;
		
		#####################################
		## Store all input lines of corresponding input file in a array list
		my $aref_in_lines = exists ( $href_in_of_file->{$file} ) ? \@{$href_in_of_file->{$file}} : () ;
		#####################################
		
		#####################################
		## Extract and separate allele and size lines from the array list of input lines
		my @alleles = grep ( /^Allele/ , @$aref_in_lines ) ;
		my @sizes   = grep ( /^Size/   , @$aref_in_lines ) ;
		#####################################
		
		my $N_bins = @alleles ;
		
		#####################################
		## Reject possible tabstopsigns at each line end
		$aref_in_lines->[1] =~ s/\t+$//g ;
		#####################################
		
		############################################################################################################### Connection and Translation of single BIN states
		## Connnect each taxonname from the sampleline with corresponding and translated allele states to
		## Taxon\tBIN1\tBIN2...\tBINn\n
		my %bin_of_tax  = () ; # Key: Taxonname; Value: List of corresponding BINs
		my @taxon_names = split "\t" , $aref_in_lines->[1] ; shift @taxon_names;
		
		for my $allele ( @alleles ){
			
			#####################################
			## Sample all BINs from a BIN line as seperate list values and reject the Linedeclaration (First column)
			my @single_bins = split "\t" , $allele ; shift @single_bins ;
			#####################################
			
			#####################################
			## Sample all single BINs of each corresponding taxon as hashlist
			## If BINvalues are found in the ABI matrix, replace them with 1. if no entry is found, replace it with 0.
			## Questionmarks are taken over.
			for my $taxon_number ( 0 .. @taxon_names-1 ){
				
				my $single_bin = () ;
				
				if 		( $single_bins[$taxon_number]	=~ /\d+/ )	{ $single_bin =	1	}
				elsif 	( $single_bins[$taxon_number]	=~ /\?/  )	{ $single_bin =	'?'	}
				else 													{ $single_bin =	0	}
				
				#####################################
				## Store list of single bins as hashlist '$bin_of_tax': Key: Taxonname; Value: List of BIN states
				push @{$bin_of_tax{$taxon_names[$taxon_number]}} , $single_bin ;
				#####################################
			}
			#####################################
		}
		#####################################
		
		############################################################################################################### Identification and Convertion of Xmin and Xmax size values
		## convert size-lines to Xmin and Xmax lines
		my $max_value_line = "\"XMax\"" ;
		my $min_value_line = "\"XMin\"" ;
		
		#####################################
		## Remove declaration 'Size' from each line
		## If size line includes one or multiple size numbers, sort sizes
		## and join minimum and maximum size to the corresponding size line (min or max)
		for my $size_line ( @sizes ){
			
			#####################################
			## Split size line and remove line name
			my @single_sizes = split "\t" , $size_line ; shift @single_sizes ;
			#####################################
			
			if ( $size_line =~ /\d+/ ){
				
				#####################################
				## sort size values
				my @sort_sizes = sort {$a<=>$b} ( @single_sizes ) ;
				#####################################
				
				#####################################
				## Take first size value (minimum) and last size value (maximum)
				## It is important to push back the maximum value
				## because it could happen that only one size (max equal min) is given in a line
				my ( $min_value , $max_value ) = () ;
				until ( $max_value =~ /\d+/ ){ $max_value = pop   @sort_sizes } push @sort_sizes, $max_value ; 
				until ( $min_value =~ /\d+/ ){ $min_value = shift @sort_sizes }
				#####################################
				
				#####################################
				## Concatenation of inferred max and min BIN sizes
				$max_value_line = $max_value_line.",".$max_value ;
				$min_value_line = $min_value_line.",".$min_value ;
				#####################################
			}
			else{ die "\n\tFILE-ERROR: No size value found in $file!\n\n" }
		}
		#####################################
		
		############################################################################################################### Generate BECKMAN format
		## generate info header in BECKMAN format for amare-analyses
		## and supermatrix output
		my $bin_line		= "\"BIN\""			; # Dummy for AMARE and Supermatrix Output
		my $dye_line		= "\"Dye\""			; # Dummy for AMARE and Supermatrix Output
		my $sample_line	= "\"Samples\""		; # Dummy for AMARE and Supermatrix Output
		my $fragment_line	= "\"Fragments\""	; # Dummy for AMARE and Supermatrix Output
		my $xmean_line	= "\"XMean\""		; # Dummy for AMARE and Supermatrix Output
		my $xvar_line		= "\"XVar\""		; # Dummy for AMARE and Supermatrix Output
		my $ymean_line	= "\"YMean\""		; # Dummy for AMARE and Supermatrix Output
		
		for ( 1 .. $N_bins ){
			
			$bin_line			= $bin_line.",\"".$_."\""		;
			$dye_line			= $dye_line.",D4"				;
			$sample_line		= $sample_line.",".$_			;
			$fragment_line	= $fragment_line.",".$_		;
			$xmean_line		= $xmean_line.",\"".$_."\""	;
			$xvar_line			= $xvar_line.",\"".$_."\""	;
			$ymean_line		= $ymean_line.",\"".$_."\""
		}
		
		my @infile = ( $bin_line, $dye_line, $sample_line, $fragment_line, $min_value_line, $max_value_line, $xmean_line, $xvar_line, $ymean_line ) ;
		#####################################
		
		#####################################
		## Generate sample block in BECKMAN format
		## and store the whole converted input file
		## in a hashlistreference '$href_in_of_file': Key: Filename; Value: input lines
		for my $tax ( sort keys %bin_of_tax ){
			
			my @line_parts = exists ($bin_of_tax{$tax}) ? @{$bin_of_tax{$tax}} : () ; 
			my $tax_line   = join "," , @line_parts ; 
			
			push @infile , $tax.",".$tax_line ;
		} 
		
		@{$href_in_of_file->{$file}} = @infile
		####################################
	}
}


sub concatenation{
	
	#####################################################
	# Subroutine 'concatenation'
	# Concatenation of corresponding sequence strings from different input files
	# Missing sequence BIN strings of double and single samples are filled with '?'
	#####################################################
	
	my $href_data_of_file	= $_[1] ; # Key: Filename; Value: List of corresponding input lines (defined in subroutine 'readin_file')
	my $href_mis_tax_of_file	= $_[2] ; # Key: Filename; Value: Taxonnames which are found in other input files, but which are not included in the corresponding input file (defined in subroutine 'readin_file')
	
	my @head_bin		= "Bin"			;
	my @head_dye		= "Dye"			;
	my @head_sam		= "Samples"		;
	my @head_fra		= "Fragments"	;
	my @head_xmi		= "XMin"		;
	my @head_xma		= "XMax"		;
	my @head_xme		= "XMean"		;
	my @head_xva		= "XVar"		; 
	my @head_yme		= "YMean"		;
	my %conc_rep		= ()			; # Key: Taxon names of replicate samples; Value: concatenated sequence
	my %conc_res		= ()			; # Key: Taxon names of single samples; Value: concatenated sequence
	my %label_of_line	= ()			;
	my @infile			= ()			;
	
	############################################################################################################### String sampling
	## Sampling of corresponding sequence strings from different input files
	for my $file ( sort keys %$href_data_of_file ){
		
		print "\n\tConcatenation of file $file" ;
		
		#####################################
		## generate an array list of existing and missing sample lines
		my $aref_indata = exists( $href_data_of_file   ->{$file} ) ? \@{$href_data_of_file   ->{$file}} : () ;
		my $aref_mistax = exists( $href_mis_tax_of_file->{$file} ) ? \@{$href_mis_tax_of_file->{$file}} : () ;
		#####################################
		
		#####################################
		## Generate a list of characters to identify the total number of characters
		my @N_columns = split ","  , $aref_indata->[0] ;
		#####################################
		
		#####################################
		## Generate a sequence of missing data (?) with original character length
		my @N_parts   = () ; 
		for ( my $k=1 ; $k<= @N_columns-1 ; $k++ ){ push @N_parts , "?" }
		#####################################
		
		#####################################
		## identification of missing taxa and adding of missing data string to corresponding taxon sequence
		for my $taxa ( @$aref_mistax ){ $taxa =~ s/"|'+//g ;
			
			if 		( $taxa =~ /^"?a\w+"?$/ )	{ push @{$conc_rep{$taxa}} , @N_parts }  
			else 								{ push @{$conc_res{$taxa}} , @N_parts }
		}
		#####################################
		
		#####################################
		## Sampling of existing input lines in accordance to associated content
		## in seprated arrays or hashlists
		for my $line ( @$aref_indata ){
			
			my @line_parts = split "," , $line ;
			my $label      = shift @line_parts ;
			
			if 		( $label =~ /^"?Bin"?$/i		)	{ push @head_bin				, @line_parts }
			elsif 	( $label =~ /^"?Dye"?$/i		)	{ push @head_dye				, @line_parts }
			elsif 	( $label =~ /^"?Samples"?$/i	)	{ push @head_sam				, @line_parts }
			elsif 	( $label =~ /^"?Fragments"?$/i	)	{ push @head_fra				, @line_parts }
			elsif 	( $label =~ /^"?XMin"?$/i		)	{ push @head_xmi				, @line_parts }
			elsif 	( $label =~ /^"?XMax"?$/i		)	{ push @head_xma				, @line_parts }
			elsif 	( $label =~ /^"?XMean"?$/i		)	{ push @head_xme				, @line_parts }
			elsif 	( $label =~ /^"?XVar"?$/i		)	{ push @head_xva				, @line_parts }
			elsif 	( $label =~ /^"?YMean"?$/i		)	{ push @head_yme				, @line_parts }  
			elsif 	( $label =~ /^"?a\w+"?$/		)	{ push @{$conc_rep{$label}}	, @line_parts }  
			else 										{ push @{$conc_res{$label}}	, @line_parts }
		}
		#####################################
	}
	#####################################
	
	print "\n\tConcatenated Supermatrix: AMARE_supermatrix.txt" ;
	
	#####################################
	## Storing of all sampled head info strings in @infile
	my $line_bin = join ",", @head_bin ; push @infile, $line_bin ;
	my $line_dye = join ",", @head_dye ; push @infile, $line_dye ;
	my $line_sam = join ",", @head_sam ; push @infile, $line_sam ;
	my $line_fra = join ",", @head_fra ; push @infile, $line_fra ;
	my $line_xmi = join ",", @head_xmi ; push @infile, $line_xmi ;
	my $line_xma = join ",", @head_xma ; push @infile, $line_xma ;
	my $line_xme = join ",", @head_xme ; push @infile, $line_xme ;
	my $line_xva = join ",", @head_xva ; push @infile, $line_xva ;
	my $line_yme = join ",", @head_yme ; push @infile, $line_yme ;
	#####################################
	
	#####################################
	## Storing of all double sample sequences in @infile
	## and insertion of corresponding sample name
	for my $col_one ( sort keys %conc_rep ){
		
		my	$aref_rep = exists ($conc_rep{$col_one} ) ? \@{$conc_rep{$col_one}} : () ; unshift @$aref_rep , $col_one ;
		my	$line_rep = join "," , @$aref_rep ; push @infile , $line_rep ; 
	}
	#####################################
	
	#####################################
	## Storing of all single sample sequences in @infile
	## and insertion of corresponding sample name
	for my $col_one ( sort keys %conc_res ){
		
		my	$aref_res = exists ($conc_res{$col_one} ) ? \@{$conc_res{$col_one}} : () ; unshift @$aref_res , $col_one ;
		my	$line_res = join "," , @$aref_res ; push @infile , $line_res ; 
	}
	#####################################
	
	#####################################
	## Final concatenation and print OUT of the supermatrix
	my $outline = join "\n" , @infile ;

	open OUTconc, ">AMARE_Matrices/AMARE_supermatrix.txt" ; print OUTconc $outline ; close OUTconc ;
	#####################################
	
	#####################################
	## return concatenated string sequences, each string as separate element
	return @infile
	#####################################
}

#####################################################################################################################################################################
# SUBROUTINES MAIN-ANALYSES
#####################################################################################################################################################################

sub fragment_sample_check{
	
	#####################################################
	# Subroutine 'fragment_sample_check' compares single fragment and sample values
	# and put all BIN columns with unequal values to the remove_bin array list
	my $sref_samp		= $_[0] ; # Sample String of given input file
	my $sref_frag		= $_[1] ; # Fragment String of given input file
	
	my @samples			= split ",", $$sref_samp ;
	my @fragments		= split ",", $$sref_frag ;
	my @removed_bins	= () ;

	for ( my $i = 1 ; $i <= ( @samples-1 ) ; $i++ ){ unless ( $samples[$i] =~ $fragments[$i] ){ push @removed_bins , $i } }
	return @removed_bins ;
	#####################################
}


sub loop{
	
	#####################################################
	# Subroutine '
	#####################################################
	
	my $aref_pre_cut_pos		= $_[0]  ; # List of deletable BIN marker (defined in MAIN ANALYSES)
	my $aref_d_samples		= $_[1]  ; # Collects all double sample BIN strings (defined in section MAIN ANALYSES)
	my $href_deleted_tax		= $_[2]  ; # Key: Deleted Taxon; Value: Counter
	my $sref_lv_11_01			= $_[3]  ; # Threshold of BIN repeatability of equal 11 bins in ratio to unequal bins (predefined as zero; not used yet)
	my $sref_lv_eq_neq		= $_[4]  ; # Threshold of repeatability of individual loci (defined in subroutine 'ARGV')
	my $sref_lv_tax_cut		= $_[5]  ; # Threshold of Sample repeatability of equal 11 bins in ratio to unequal bins (predefined in 'MAIN ANALYSES')
	my $sref_lv_dist			= $_[6]  ; # Threshold for minimum allowed BIN distances (defined in subroutine 'ARGV')
	my $sref_li_minima		= $_[7]  ; # String of minimum BIN distances (defined in MAIN ANALYSES)
	my $sref_li_maxima		= $_[8]  ; # String of maximum BIN distances (defined in MAIN ANALYSES)
	my $href_res_tax_cut		= $_[9]  ; # Key: Infostring; Value: all deleted BIN positions and Replicate samples
	my $error_rate_type		= $_[10] ; # Type of chosen error rate (defined in ARGV)
	my $sref_handling_bins	= $_[11] ; # Type of remaining bin states, if 0 -> remain bins which only contain 0|0 bin states (defined in ARGV with -r option)
	
	my %bins_remain			= () ; # Key: Positions of remaining BINs; Value: Counter (defined in check_bin_states)
	my %deleted_bin			= () ; # Key: deleted BIN column; Value: Counter
	my %bins_equal_1		= () ; # Key: BIN position; Value: Number of equal 11 states in all double samples (defined in check_bin_states)
	my %bins_equal_0		= () ; # Key: BIN position; Value: Number of equal 00 states in all double samples (defined in check_bin_states)
	my %bins_unequal		= () ; # Key: BIN position; Value: Number of unequal states in all double samples (defined in check_bin_states)
	my %bins_missing		= () ; # Key: BIN position; Value: Number of equal 11 states in all double samples (defined in check_bin_states)
	my %tax_remain			= () ; # Key: Name of remaining taxa; Value: Counter (defined in check_bin_states)
	my %tax_equal_1			= () ; # Key: Taxonname of one double sample; Value: Number of equal 11 states (defined in check_bin_states)
	my %tax_equal_0			= () ; # Key: Taxonname of one double sample; Value: Number of equal 00 states (defined in check_bin_states)
	my %tax_unequal			= () ; # Key: Taxonname of one double sample; Value: Number of unequal states (defined in check_bin_states)
	my %tax_missing			= () ; # Key: Taxonname of one double sample; Value: Number of equal ?? states (defined in check_bin_states)
	my $N_unequal_01		= () ; # Total number of remaining unequal states
	my $N_equal_11		= () ; # Total number of remaining equal 11 states
	my $N_missing_data	= () ; # Total number of remaining missing states
	
	print LOG "LOOP-Start\n----------\n\n" ;
	
	#####################################
	## Assign each deleted BIN number of fragment and sample check to a Hashkey of '%deleted_bin'
	for ( @$aref_pre_cut_pos ){ $deleted_bin{$_}++  }
	#####################################
	
	############################################################################################################### Delete Taxa (horizontally)
	## Count all equal and unequal BIN states for each double sample
	&check_bin_states( \@$aref_d_samples, \%$href_deleted_tax, \%deleted_bin, \%tax_remain, \%bins_remain, \%tax_equal_1 , \%tax_equal_0 , \%tax_missing, \%tax_unequal, \'tax' ) ;
	#####################################
	
	#####################################
	## Obtain the number of Remaining double samples
	my $N_taxa_befor_delete = keys %tax_remain ; my @tax_names = sort {$a<=>$b} keys %tax_remain ;
	print LOG "Number Replicate-Samples:\t$N_taxa_befor_delete\n\tReplicates: @tax_names\n\n" ;
	#####################################
	
	#####################################
	## Obtain the number of Remaining BIN positions which are 
	## numerical sorted and stored in @bin_numbers
	my $N_bins_befor_delete = keys %bins_remain ; my @bin_numbers = sort {$a<=>$b} keys %bins_remain ;
	#####################################
	
	#####################################
	## Delete Taxa below given Replicate threshold
	&delete_taxa( \%tax_remain, \%tax_equal_1, \%tax_equal_0, \%tax_missing , \%$href_deleted_tax , \$N_bins_befor_delete , \$$sref_lv_tax_cut ) ; 
	my $N_taxa_after_delete = keys %tax_remain ; my @taxa = keys %tax_remain ;
	my $N_taxa_deleted      = $N_taxa_befor_delete - $N_taxa_after_delete ;
	print LOG "Number of deleted Replicates: $N_taxa_deleted\n\tNumber of remaining Replicates: $N_taxa_after_delete\n\tRemaining Replicates: @taxa\n\n" ;
	#####################################
	
	#####################################
	## If no taxa left after replicate check abort Loop and goto analyse loop
	## At least one replicate sample should be left after analyses with AMARE
	if   ( $N_taxa_after_delete == 0 ){ print LOG "LOOP-Abort\n--------\n\tNo Taxa left !!\n\n" ; &error_rate( \%$href_res_tax_cut , \'0' , \$$sref_lv_tax_cut , \'0', \'0' , \'0' , \'0' , \$$sref_lv_eq_neq , \'0' , \'0' , \$$error_rate_type ); goto LOOPEND }	
	#####################################
	
	############################################################################################################### Delete BINs (vertically)
	## Count all equal and unequal BIN states for each remaining BIN column
	&check_bin_states( \@$aref_d_samples, \%$href_deleted_tax, \%deleted_bin, \%tax_remain, \%bins_remain, \%bins_equal_1 , \%bins_equal_0 , \%bins_missing, \%bins_unequal, \'bin' ) ;
	#####################################
	
	#####################################
	## Obtain the number of Remaining BIN positions which are 
	## numerical sorted and stored in @bin_numbers
	my $N_bins_befor_delete = keys %bins_remain ; my @bin_numbers = sort {$a<=>$b} keys %bins_remain ;
	print LOG "Number of Bins:\t$N_bins_befor_delete\n\tBins: @bin_numbers\n" ;
	#####################################
	
	#####################################
	## Delete BINs below given BIN threshold or without equal 11 states
	&delete_bins( \%deleted_bin , \%bins_remain, \%bins_equal_1 , \%bins_equal_0 , \%bins_unequal , \%bins_missing, \$N_taxa_after_delete , \$N_unequal_01 , \$N_equal_11, \$$sref_lv_11_01 , \$$sref_lv_eq_neq , \$N_missing_data, \$$sref_handling_bins ) ; 
	my $N_bins_after_delete = keys %bins_remain ; @bin_numbers = sort {$a<=>$b} keys %bins_remain ;
	print LOG "Number of remaining Bins: $N_bins_after_delete\n\tRemaining Bins: @bin_numbers\n\n" ;
	#####################################
	
	#####################################
	## Sort and store all deleted BINs in @cut_bin_loop
	my @cut_bin_loop         = sort {$a<=>$b} keys %deleted_bin ;
	#####################################
	
	############################################################################################################### BIN-Distance check
	## Check the distances between single BIN states and remove all BINs below the given BIN-Distance Threshold
	my @cut_bin_all			= &bin_distance_check( \@cut_bin_loop, \%bins_remain, \$$sref_lv_dist, \$$sref_li_minima, \$$sref_li_maxima ) ;
	my $N_bin_remain_dist	= keys %bins_remain ; @bin_numbers = sort {$a<=>$b} keys %bins_remain ;
	print LOG "Number of remaining Bins: $N_bin_remain_dist\n\tRemaining Bins: @bin_numbers\n\n" ;
	
	#####################################
	## Inferre number of deleted BINs before and after Distance Check
	## for Run-Check in LOG file
	my $N_cut_bin_loop		= @cut_bin_loop	; 
	my $N_cut_bin_all			= @cut_bin_all	;
	print LOG "\n\tRun-Check\n\t---------\n\tNumber of deleted bins before distance-check: $N_cut_bin_loop\n\tNumber of deleted bins after distance-check: $N_cut_bin_all\n\t" ;
	#####################################
	
	#####################################
	## unless 5 BINs remain after distance check abort loop
	## 5 BIN marker should be remain after AMARE analyses
	unless ( ( $N_bin_remain_dist ) >= 3			){ print LOG "\nLOOP-Abort\n--------\n\tBin-Number below 5 !!\n\n"						; &error_rate( \%$href_res_tax_cut , \@cut_bin_all , \$$sref_lv_tax_cut , \$N_taxa_after_delete , \'0' , \'0' ,  \'0' , \$$sref_lv_eq_neq , \'0' , \'0' , \$$error_rate_type ) }
	#####################################
	
	#####################################
	## If no BINs are removed in Distance Check 
	## else repeat Replicate and BIN check 
	if 		( ( @cut_bin_loop ) == ( @cut_bin_all )	){ print LOG "\nLOOP-End\n--------\n\tNo further bins deleted in distance-check\n\n" 	; &error_rate( \%$href_res_tax_cut , \@cut_bin_all , \$$sref_lv_tax_cut , \$N_taxa_after_delete , \$N_bin_remain_dist , \$N_unequal_01 ,  \$N_equal_11 , \$$sref_lv_eq_neq , \%$href_deleted_tax , \$N_missing_data , \$$error_rate_type ) }
	else{ print LOG "\nLOOP-continue\n------------\n\tBins deleted in distance-check\n\n"													; &loop( \@cut_bin_all , \@$aref_d_samples , \%$href_deleted_tax , \$$sref_lv_11_01 , \$$sref_lv_eq_neq ,  \$$sref_lv_tax_cut , \$$sref_lv_dist , \$$sref_li_minima , \$$sref_li_maxima , \%$href_res_tax_cut , \$$error_rate_type, \$$sref_handling_bins ) } 
	#####################################
	
	LOOPEND:;
	############################################################################################################### Subroutines within 'Loop'
	
	sub check_bin_states{
		
		my $aref_double_taxa		= $_[0] ; # Collect all double sample BIN strings (defined in section MAIN ANALYSES)
		my $href_seen_tax_delete	= $_[1] ; # Key: Taxonname which will be excluded; Value: counter (defined in delete_taxa)
		my $href_seen_bin_delete	= $_[2] ; # Key: deleted BIN column; Value: Counter
		my $href_remain_tax		= $_[3] ; # Key: Name of remaining taxa; Value: Counter (defined in check_bin_states)
		my $href_remain_bin		= $_[4] ; # Key: Positions of remaining BINs; Value: Counter (defined in check_bin_states)
		my $href_equal_1			= $_[5] ; # Key: Taxonname or BIN position; Value: Number of equal 11 states (defined in check_bin_states)
		my $href_equal_0			= $_[6] ; # Key: Taxonname or BIN position; Value: Number of equal 00 states (defined in check_bin_states)
		my $href_missing			= $_[7] ; # Key: Taxonname or BIN position; Value: Number of equal ?? states (defined in check_bin_states)
		my $href_unequal			= $_[8] ; # Key: Taxonname or BIN position; Value: Number of unequal  states (defined in check_bin_states)
		my $sref_code				= $_[9] ; # Important to know which kind of states should be compared
		
		#####################################
		## Compare repeatability of paired double samples
		for ( my $i = 0 ; $i <= @$aref_double_taxa-1 ; $i += 2 ){ my $j = $i+1 ; 
			
			my @sample_1_states = split ",", $aref_double_taxa->[$i] ; 
			my @sample_2_states = split ",", $aref_double_taxa->[$j] ;
			
			############################################################################################################### Check repeatability of single BINstates between double samples
			## unless double samples are not already excluded or BIN states rejected
			unless ( $href_seen_tax_delete->{ $sample_1_states[0] } ){
				
				#####################################
				## For all BIN states...
				for ( my $k = 1 ; $k <= @sample_1_states-1 ; $k++ ){ 
					
					#####################################
					## unless BIN state is not already excluded
					unless ( $href_seen_bin_delete->{$k} ){
						
						#####################################
						## Count BIN state and associated Taxon
						$href_remain_bin->{$k}++ ; $href_remain_tax->{$sample_1_states[0]}++ ;
						#####################################
						
						if ( $$sref_code =~ /tax/ ){
							
							#####################################
							## Count total repeatability of taxon states
							if 		( ( $sample_1_states[$k] =~	1	) && ( $sample_2_states[$k] =~	1		) )	{ $href_equal_1->{ $sample_1_states[0]}++ }
							elsif	( ( $sample_1_states[$k] =~	0	) && ( $sample_2_states[$k] =~	0		) )	{ $href_equal_0->{ $sample_1_states[0]}++ }
							elsif	( ( $sample_1_states[$k] =~ /\?/	) || ( $sample_2_states[$k] =~	/\?/	) )	{ $href_missing->{ $sample_1_states[0]}++ }
							else 																						{ $href_unequal->{ $sample_1_states[0]}++ }
							#####################################
							
						}
						
						else{
							
							#####################################
							## Count total repeatability of BIN states
							if 		( ( $sample_1_states[$k] =~	1	) && ( $sample_2_states[$k] =~	1		) )	{ $href_equal_1->{ $k }++ }
							elsif 	( ( $sample_1_states[$k] =~	0	) && ( $sample_2_states[$k] =~	0		) )	{ $href_equal_0->{ $k }++ }
							elsif 	( ( $sample_1_states[$k] =~ /\?/	) || ( $sample_2_states[$k] =~	/\?/	) )	{ $href_missing->{ $k }++ }
							else 																						{ $href_unequal->{ $k }++ }
							#####################################
						}
					}
					#####################################
				}
				#####################################
			}
			#####################################
		}
		#####################################
		
		if ( $$sref_code =~ /tax/ ){ print LOG "\tReplicate-Check\n\t---------------\n\t" } else { print LOG "\n\tBIN-Check\n\t---------\n\t" }
	}
	
	sub delete_taxa{
		
		my $href_taxa_remain		= $_[0] ; # Key: Name of remaining taxa; Value: Counter (defined in check_bin_states)
		my $href_taxon_equal_1	= $_[1] ; # Key Taxonname of one double sample; Value: Number of equal 11 states (defined in check_bin_states)
		my $href_taxon_equal_0	= $_[2] ; # Key Taxonname of one double sample; Value: Number of equal 00 states (defined in check_bin_states)
		my $href_taxon_missing	= $_[3] ; # Key Taxonname of one double sample; Value: Number of equal ?? states (defined in check_bin_states)
		my $href_tax_deleted		= $_[4] ; # Key: Taxonname which will be excluded; Value: counter (defined in delete_taxa)
		my $sref_Nbin_delete		= $_[5] ; # Number of remaining BINs
		my $sref_lev_tax_cut		= $_[6] ; # Threshold of Sample repeatability of equal 11 bins in ratio to unequal bins (predefined in 'MAIN ANALYSES')
		
		#####################################
		## Delete taxa under given thresholds as haskey in '$href_taxa_remain'
		## if repeatability of BIN states between double samples is below the replicate threshold
		## Calculation: Number of equal states between double sample / Number of BINs - Number of missing data
		for my $tax ( keys %$href_taxa_remain ){
			
			my $percent_eq_tax = sprintf "%.1f", ( ( ( $href_taxon_equal_0->{$tax} + $href_taxon_equal_1->{$tax} ) / ( $$sref_Nbin_delete - $href_taxon_missing->{$tax} ) ) * 100 ) ;
			
			if ( ( $$sref_lev_tax_cut ) >= ( $percent_eq_tax ) ){ 
				
				$href_tax_deleted->{ $tax }++ ; delete $href_taxa_remain->{ $tax } }
		}
		#####################################
		
		print LOG "\n\tDeleted Replicates\n\t------------------\n\t" ;
	}
	
	sub delete_bins{
		
		my $href_del_bin_check	= $_[0]  ; # Key: Deleted BINs; Value Counter (defined in delete_bins)
		my $href_bins_remain		= $_[1]  ; # Key: Positions of remaining BINs; Value: Counter (defined in check_bin_states)
		my $href_bins_equal_1	= $_[2]  ; # Key: BIN position; Value: Number of equal 11 states in all double samples (defined in check_bin_states)
		my $href_bins_equal_0	= $_[3]  ; # Key: BIN position; Value: Number of equal 00 states in all double samples (defined in check_bin_states)
		my $href_bins_unequal	= $_[4]  ; # Key: BIN position; Value: Number of unequal states in all double samples (defined in check_bin_states)
		my $href_bins_missing	= $_[5]  ; # Key: BIN position; Value: Number of equal 11 states in all double samples (defined in check_bin_states)
		my $sref_N_remain_tax	= $_[6]  ; # Number of remaining taxa
		my $sref_N_unequal_01	= $_[7]  ; # Total number of remaining unequal states
		my $sref_N_equal_11		= $_[8]  ; # Total number of remaining equal 11 states
		my $sref_lev_11_01		= $_[9]  ; # Threshold of BIN repeatability of equal 11 bins in ratio to unequal bins (predefined as zero; not used yet)
		my $sref_lev_eq_neq		= $_[10] ; # Threshold of repeatability of individual loci (defined in subroutine 'ARGV')
		my $sref_missing_data	= $_[11] ; # Total number of remaining missing states
		my $sref_bin_handling	= $_[12] ; # Remain complete 0|0 bins (=0) or reject complete 0|0 bins (=1)
		
		
		#####################################
		## For LOG file output
		my $N_del_tota_bins		= 0 ;
		my $N_del_cons_00_bins	= 0 ;
		my $N_del_cons_01_bins	= 0 ;
		#####################################
		
		#####################################
		## Remove remaining bins if...
		COLUMN:
		for my $col ( sort {$a<=>$b} keys %$href_bins_remain ){
			
			############################################################################################################### Remove BIN all columns without 11 states
			
			## ...all BIN states include vertically only 00 states and missing data unless -r option = 0
			if 	( ( $href_bins_equal_0->{$col} == $$sref_N_remain_tax - $href_bins_missing->{$col} ) && ( $$sref_bin_handling == 1 ) ){ 
				
					$href_del_bin_check->{$col}++ ; delete $href_bins_remain->{$col} ; $N_del_tota_bins++ ;
				(	$href_bins_equal_0->{$col} == $$sref_N_remain_tax - $href_bins_missing->{$col} ) ? $N_del_cons_00_bins++ : $N_del_cons_01_bins++;
				next COLUMN;
			}
			elsif ( ( $href_bins_equal_0->{$col} == $$sref_N_remain_tax - $href_bins_missing->{$col} ) && ( $$sref_bin_handling == 0 ) ){ next COLUMN }
			#####################################
			
			############################################################################################################### Remove BIN all columns without 11 states and only 00 states
			
			## ...all BIN states include vertically only 00 and 01 states and missing data
			elsif 	( $href_bins_equal_0->{$col} + $href_bins_unequal->{$col} == $$sref_N_remain_tax - $href_bins_missing->{$col} ){ 
				
					$href_del_bin_check->{$col}++ ; delete $href_bins_remain->{$col} ; $N_del_tota_bins++ ;
				(	$href_bins_equal_0->{$col} == $$sref_N_remain_tax - $href_bins_missing->{$col} ) ? $N_del_cons_00_bins++ : $N_del_cons_01_bins++;
				next COLUMN
			}
			#####################################
			
			############################################################################################################### Remove BIN columns below BIN threshold
			## else calculate percentages of 11 and 00 states without incluing missing data
			## if percentage of 11 and 00 states is lower/equal as the given BIN threshold delete BIN column
			else{
				
				#####################################
				## Calculate single percentages of 11 and 00 states
				my $percent_equal_1 = sprintf "%.1f", ( ( $href_bins_equal_1->{$col} / ( $$sref_N_remain_tax - $href_bins_missing->{$col} ) ) * 100 ) ;
				my $percent_equal_0 = sprintf "%.1f", ( ( $href_bins_equal_0->{$col} / ( $$sref_N_remain_tax - $href_bins_missing->{$col} ) ) * 100 ) ;
				#####################################
				
				#####################################
				## Compare complete percentage value of equal 11 and 00 states with BIN threshold
				if  ( $$sref_lev_eq_neq >= $percent_equal_1 + $percent_equal_0 ){
					
					$href_del_bin_check->{$col}++ ; $N_del_tota_bins++ ; delete $href_bins_remain->{$col};
				}
				#####################################
				
				#####################################
				## Add Number of unequal, equal 11, and missing states to total number of identified unequal, equal 11, and missing states
				else{ $$sref_N_unequal_01 += $href_bins_unequal->{$col}; $$sref_N_equal_11 += $href_bins_equal_1->{$col}; $$sref_missing_data += $href_bins_missing->{$col} }
				#####################################
			}
			#####################################
		}
		#####################################
		
		print LOG "\n\n\tDeleted Bins\n\t------------\n\tNumber of deleted BIN without any congruent (1,1) BIN state\n\t(0,0): $N_del_cons_00_bins\n\t(0,1): $N_del_cons_01_bins\n\tTotal number of deleted bins: $N_del_tota_bins\n\t" ;
	}
	
	sub bin_distance_check{
		
		my $aref_del_bins		= $_[0] ; # list of deleted BIN positions after Replicate and BIN check
		my $href_remain_bins	= $_[1] ; # Key: Positions of remaining BINs; Value: Counter (defined in check_bin_states)
		my $sref_dist_lv		= $_[2] ; # Threshold for minimum allowed BIN distances (defined in subroutine 'ARGV')
		my $sref_Xmin_values	= $_[3] ; # String of minimum BIN distances (defined in MAIN ANALYSES)
		my $sref_Xmax_values	= $_[4] ; # String of maximum BIN distances (defined in MAIN ANALYSES)
		
		my %seen_cut			= () ; # Key: Deleted BINs; Value: Counter
		my @remain_pos			= () ; # List of remaining Xmin values
		my $del_bins_befor	= () ; # Number of previously deleted BINs
		
		#####################################
		## Store all Xmin and Xmax lengths of single BINs in corersponding arrays
		my @single_Xmin = split "," , $$sref_Xmin_values ;
		my @single_Xmax = split "," , $$sref_Xmax_values ;
		#####################################
		
		#####################################
		## Store deleted BIN positions as haskey to filter remaining Xmin Values
		for my $deleted_pos ( @$aref_del_bins ){ $seen_cut{$deleted_pos}++ } $del_bins_befor = keys %seen_cut ;
		#####################################
		
		#####################################
		## Store Xmin values of remaining BINs in @remaining positions
		for ( my $i = 1 ; $i <= @single_Xmin-1 ; $i++ ){ unless ( $seen_cut{$i} ){ push @remain_pos , $i }}
		#####################################
		
		#####################################
		## Calculate Distances between Xmin and Xmax values
		for ( my $i = 0 ; $i < @remain_pos-1   ; $i++ ){ my $j = $i+1 ;
			
			$single_Xmin[$remain_pos[$j]] =~ s/"//g ;
			$single_Xmax[$remain_pos[$i]] =~ s/"//g ;
			
			my   $distance = sprintf "%.3f" ,  abs ( $single_Xmax[$remain_pos[$i]] - $single_Xmin[$remain_pos[$j]] ) ; 
			
			#####################################
			## Delete BIN if inferred distance is below Distance-Threshold
			if ( $$sref_dist_lv >= $distance ){ $seen_cut{$remain_pos[$i]}++ ; $seen_cut{$remain_pos[$j]}++ , delete $href_remain_bins->{$i} , delete $href_remain_bins->{$j} }
			#####################################
		}
		#####################################
		
		#####################################
		## Store all deleted BINs in @rejected_pos_all
		my @rejected_pos_all = sort {$a<=>$b} keys %seen_cut ;
		#####################################
		
		#####################################
		## Calculate Number of deleted BINs caused by Distance error
		my $deleted_bins     = @rejected_pos_all - $del_bins_befor ;
		#####################################
		
		print LOG "\n\tDistance-Check\n\t--------------\n\tNumber of deleted Bins: $deleted_bins\n\t" ;
		
		return @rejected_pos_all
	}
	
	sub error_rate{
		
		print LOG "LOOP-Analysis\n-------------\n\t" ;
		
		my $href_result_taxcut	= $_[0]  ; # Key: Infostring; Value: all deleted BIN positions and Replicate samples
		my $aref_cuts_final		= $_[1]  ; # List of all deleted BINs
		my $sref_taxon_cut_level	= $_[2]  ; # Replicate-Threshold
		my $sref_N_remaining_tax	= $_[3]  ; # Number of remaining replicate taxa
		my $sref_N_remaining_bin	= $_[4]  ; # Number of remaining BINs
		my $sref_N_uneq_01		= $_[5]  ; # Number of unequal 01 BINs
		my $sref_N_eq_11			= $_[6]  ; # Number of equal 11 BINs
		my $sref_eq_neq_lv		= $_[7]  ; # BIN-Threshold
		my $href_taxa_deleted	= $_[8]  ; # Key: Deleted Taxon; Value: Counter
		my $sref_missing_data	= $_[9]  ; # Number of missing data
		my $sref_error_type		= $_[10] ; # Type of Error Rate
		
		my $info_string			= () ;
		my $error_rate			= () ;
		my $info_error			= () ;
		
		#####################################
		## Generate info string if less than 5 BINs or Replicate sample remain after analyses
		if ( $$sref_N_remaining_bin == 0 || $$sref_N_remaining_tax == 0 ){
			
			$info_string = $$sref_eq_neq_lv."\t".$$sref_taxon_cut_level."\t".$$sref_N_remaining_bin."\t0\t0\t0\t-1" ; 
			$href_result_taxcut->{$info_string} = -1 ; 
			
			print LOG "LOOP-Abort\n\t$info_string\n\n---END---" ; close LOG ;
		}
		#####################################
		
		#####################################
		## Calculate Error rate and generate infostring for output
		else{
			
			#####################################
			## Calculate absolute difference between remaining BINs and Sample Replicates
			## Calculate absolute number of remaining BIN states
			my $Nbin_Ntax_difference	= abs	( $$sref_N_remaining_tax - $$sref_N_remaining_bin ) ;
			my $N_character				=		( $$sref_N_remaining_tax * $$sref_N_remaining_bin ) ;
			#####################################
			
			############################################################################################################### Error Rate
			## Calculation of Error-Rate
			if 		( $$sref_error_type =~ /^e$/i ){ $error_rate = sprintf "%.5f" , ( $$sref_N_uneq_01 / ( $N_character - $$sref_missing_data ) ) ; $info_error = 'euclidean distance' } 
			elsif 	( $$sref_error_type =~ /^j$/i ){ $error_rate = sprintf "%.5f" , ( $$sref_N_uneq_01 / ( $$sref_N_uneq_01 + $$sref_N_eq_11  ) ) ; $info_error = 'jaccard distance'   } 
			#####################################
			
			#####################################
			## Generate infostring if BINs remain
			## Info string contains a list of all deleted BIN positions and Replicate samples
			$info_string		= $$sref_eq_neq_lv."\t".$$sref_taxon_cut_level."\t".$$sref_N_remaining_tax."\t".$$sref_N_remaining_bin."\t".$Nbin_Ntax_difference."\t".$N_character."\t".$error_rate ;
			my $cuts_string	= join "::" , @$aref_cuts_final ;
			my $taxa_string	= join "::" , keys %$href_taxa_deleted ; 
			
			$href_result_taxcut->{$info_string} = $cuts_string.";;".$taxa_string ;
			#####################################
			
			print LOG "eq_ne_level\ttax_level\tN_tax\tN_bin\tdiff\tN_char\t$info_error\n\t$info_string\n\n---END---" ; close LOG ;
		}
		#####################################
	}
}


sub graphical_bin_states{
	
	my $aref_double_taxa		= $_[0] ; # Collect all double sample BIN strings (defined in section MAIN ANALYSES)
	my $href_graphsymb_line	= $_[1] ; # Key: Repeatability-Threshold::DoubleSampleName Value: List of Results obtained from pairwise comparisons (defined in graphical_bin_states)
	my $sref_eq_ne_level		= $_[2] ; # Threshold of repeatability of individual loci (defined in subroutine 'ARGV')
	
	#####################################
	## Compare repeatability of paired double samples
	for ( my $i = 0 ; $i <= @$aref_double_taxa-1 ; $i += 2 ){ my $j = $i+1 ; 
		
		my @sample_1_states = split ",", $aref_double_taxa->[$i] ;
		my @sample_2_states = split ",", $aref_double_taxa->[$j] ;
		
		#####################################
		## for graphical output
		unless ( $href_graphsymb_line->{ $$sref_eq_ne_level.";;".$sample_1_states[0] } ){
			
			#####################################
			## Compare each BIN state between associated double samples
			## Store threshold of BIN repeatability together with the double sample name as haskey of '$href_graphsymb_line'
			## and the results of each comparison as value list
			## Equal 11 states -> 1
			## Equal 00 states -> 0
			## Equal ?? states -> ?
			## Unequal states  -> X
			for ( my $k = 1 ; $k <= @sample_1_states-1 ; $k++ ){ 
				
				if 		( ( $sample_1_states[$k] =~	1	) && ( $sample_2_states[$k] =~	1	) )	{ push @{ $href_graphsymb_line->{$$sref_eq_ne_level.";;".$sample_1_states[0]} } , '1' }
				elsif 	( ( $sample_1_states[$k] =~	0	) && ( $sample_2_states[$k] =~	0	) )	{ push @{ $href_graphsymb_line->{$$sref_eq_ne_level.";;".$sample_1_states[0]} } , '0' }
				elsif 	( ( $sample_1_states[$k] =~ /\?/	) || ( $sample_2_states[$k] =~ /\?/	) )	{ push @{ $href_graphsymb_line->{$$sref_eq_ne_level.";;".$sample_1_states[0]} } , '?' }
				else 																					{ push @{ $href_graphsymb_line->{$$sref_eq_ne_level.";;".$sample_1_states[0]} } , 'X' }
			}
			#####################################
		}
		#####################################
	}
}


sub taxcut_analyse{
	
	#####################################################
	## Subroutine 'taxcut_analyse' prints out all results
	## inferred from the AMARE analyses in a main log file
	## and resulting matrices as svg vector graphics
	#####################################################
	
	my $href_taxcut_result	= $_[0] ; # Key: Infostring; Value: all deleted BIN positions and Replicate samples (defined in subroutine 'error_rate')
	my $aref_inputfile		= $_[1] ; # Collects all input lines of given data (defined in section MAIN ANALYSES)
	my $sref_filename			= $_[2] ; # Name of the AMARE main log output file (predefined)
	my $href_linegraph_tax	= $_[3] ; # Key: Repeatability-Threshold::DoubleSampleName Value: List of Results obtained from pairwise comparisons (defined in graphical_bin_states)
	my $sref_error_typ		= $_[4] ; # Type of chosen Error Rate (Euclidean or Jaccard) (defined in subroutine 'ARGV')
	
	my $error_name			= () ;
	my %seen_eqlevel			= () ;
	my %highest_error_bin_ratio	= () ;
	my $loguser_filename		= "AMARE_main_log.txt" ;
	
	#####################################
	## Main log file - Header
	if ( $$sref_error_typ =~ /^e$/i ){ $error_name = 'euclidean distance' } else { $error_name = 'jaccard distance' }
	
	open  LOGuser, ">AMARE_logfiles/$loguser_filename" or warn "\nERROR: Can not print out AMARE_Matrices/$loguser_filename\n" ;
	print          "\n\n\tl-EQ\tl-T\tN_tax\tN_bin\tdiff\tN_char\t$error_name" ;
	print LOGuser  "AMARE_1.3beta.pl MAIN-LOG-File\n\nInputfile: $$sref_filename\n\n\tl-EQ\tl-T\tN_tax\tN_bin\tdiff\tN_char\t$error_name\tMatrix-Outfile:" ;
	#####################################
	
	#####################################
	## Print all info results in the main log file
	for my $info ( sort keys %$href_taxcut_result ){
		
		my @info_results = split "\t" , $info ;
		
		#####################################
		## Include newline signs between different AMARE BIN threshold results
		unless ( $seen_eqlevel{$info_results[0]} ){ print "\n" ; print LOGuser "\n\n" } ;
		$seen_eqlevel{$info_results[0]}++;
		#####################################
		
		#####################################
		## if error rate is below 0.1 and...
		if ( $info_results[6] <= 0.1 ){ 
			
			#####################################
			## unless -1 (which indicate that not
			## enough BINs or Replicate Samples are left
			unless ( $info_results[6] == -1 ){ 
				
				#####################################
				## print complete results in main log and resulting AFLP matrices as svg graphic plot
				print LOGuser "\t$info" ; print "\n\t$info" ; 
				
				&print_out ( \$href_taxcut_result->{$info}  , \$info_results[0]     , \@$aref_inputfile , \$info_results[1] ) ;
				&print_svg ( \$href_taxcut_result->{$info}  , \%$href_linegraph_tax , \$info_results[0] , \$info_results[1] ) ;
				#####################################
			}
			#####################################
			
			#####################################
			## print error prompt
			else   { print LOGuser "\t$info_results[0]\t$info_results[1]\tno remaining taxa or bins" ; print "\n\t$info_results[0]\t$info_results[1]\tno remaining taxa or bins" }
			#####################################
		} 
		
		else{ 
			
			#####################################
			## print error prompt
			print         "\t$info_results[0]\t$info_results[1]\t$error_name > 0.1\n" ;
			print LOGuser "\t$info_results[0]\t$info_results[1]\t$error_name > 0.1\n" 
			#####################################
		}
		#####################################
	}
	
	close LOGuser ;
	#####################################
	
	
	
	&end ;
	
	sub print_svg{
		
		#####################################################
		## Subroutine 'print_svg' prints out all resulting matrices
		## generated from the AMARE analyses as svg vector graphics
		#####################################################
		
		my $sref_taxcut_results	= $_[0]	; # Key: Infostring; Value: all deleted BIN positions and Replicate samples (defined in subroutine 'error_rate')
		my $href_graph_of_lv_tax	= $_[1]	; # Key: Repeatability-Threshold::DoubleSampleName Value: List of Results obtained from pairwise comparisons (defined in graphical_bin_states)
		my $sref_lv_eq_ne			= $_[2]	; # BIN-Threshold of repeatability of individual loci (defined in subroutine 'ARGV')
		my $sref_tax_cut_lv		= $_[3]	; # Replicate-Threshold of repeatability of individual loci (defined in subroutine 'ARGV')
		
		my %graph_values_of_tax		= ()	;
		
		my @cut_b_t	= split ";;" , $$sref_taxcut_results ;
		my @cut_bin	= split "::" , $cut_b_t[0] ;
		my @cut_tax	= split "::" , $cut_b_t[1] ;
		
		my %seen_del_bin			= ()	;
		my %seen_del_tax			= ()	;
		
		for my $bin ( @cut_bin ){ $seen_del_bin{$bin}++ }
		for my $tax ( @cut_tax ){ $seen_del_tax{$tax}++ }
		
		my @print_1					= ()	;
		my @print_1_head			= "\t"	;
		my @print_2					= ()	;
		my $line_print_2h			= ()	;
		my $N_char					= ()	;
		my %svg_tax_states_red		= ()	;
		my %svg_tax_states_unred	= ()	;
		
		#####################################################
		## Prepare all variables for svg.plot of reduced
		## and unreduced matrices
		for my $lveq_tax ( sort keys %$href_linegraph_tax ){ 
			
			my @parts_lv_tax = split ";;" , $lveq_tax ;
			if ( $$sref_lv_eq_ne =~ /^$parts_lv_tax[0]$/ ){ 
				
				my $aref_graph_values = exists ( $href_linegraph_tax->{$lveq_tax} ) ? \@{$href_linegraph_tax->{$lveq_tax}} : () ; $N_char = @$aref_graph_values ;
				
				#####################################################
				## SVG-unreduced
				push @{$svg_tax_states_unred{$parts_lv_tax[1]}} , @$aref_graph_values ;
				#####################################################
				
				#####################################################
				## SVG-reduced
				unless ( $seen_del_tax{$parts_lv_tax[1]} ){ 
					
					push my @print   , $parts_lv_tax[1] ;
					push my @print_h , " " ;
					
					for ( my $bin=0; $bin<=$N_char-1 ; $bin++ ){ my $bin_org = $bin ; $bin_org++;
						
						unless ( $seen_del_bin{$bin_org} ){
						
							push @print_h , $bin_org ;
							push @print   , $aref_graph_values->[$bin]
						}
					}
					
					shift @print ; push @{$svg_tax_states_red{$parts_lv_tax[1]}} , @print ;
					
					$line_print_2h = join "\t", @print_h ; 
				}
				#####################################################
			}
		}
		#####################################################
		
		for ( my $bin=1; $bin<=$N_char ; $bin++ ){ push @print_1_head , "$bin\t" } my @head_reduced = split "\t", $line_print_2h ;
		
		&svg_plot ( \%svg_tax_states_unred , \'unreduced' , \$$sref_lv_eq_ne , \@print_1_head                         ) ;
		&svg_plot ( \%svg_tax_states_red   , \'reduced'   , \$$sref_lv_eq_ne , \@head_reduced , \$$sref_tax_cut_lv ,  ) ;
		
		
		sub svg_plot{
			
			my $href_tax_states	= $_[0] ;
			my $sref_name_part	= $_[1] ;
			my $sref_eq_ne		= $_[2] ;
			my $aref_bin_numbers	= $_[3] ;
			my $sref_taxcut		= $_[4] ;
			
			my $filename = () ;
			
			if ( defined $$sref_taxcut ){ $filename = $$sref_eq_ne."_".$$sref_taxcut."_".$$sref_name_part.".svg" } else { $filename = $$sref_eq_ne."_".$$sref_name_part.".svg" }
			
			my $nrows		= keys	%$href_tax_states ;	$nrows		*= 10 ; 
			my $ncolumns	= @$aref_bin_numbers-1 ;	$ncolumns	*= 10 ;
			
			my $init_line	= '<?xml version="1.0" encoding="UTF-8" standalone="no"?>'	;
			my $gen_line 	= '<!-- created by AMARE.1.3beta.pl -->'					;
			
			open my $fh_matrix_out , ">AMARE_svgfiles/$filename"	; 

			my $width		= $ncolumns	+ 2 ;
			my $height		= $nrows		+ 2 ;
			my @taxa		= () ;

			print $fh_matrix_out <<FRAME;
$init_line
$gen_line

<svg transform="translate(10,70)"
   xmlns:svg="http://www.w3.org/2000/svg"
   xmlns="http://www.w3.org/2000/svg"
   version="1.3beta"
   width="$width"
   height="$height"
   id="svg2">

  <defs
     id="defs4" />

<rect
     width="$width"
     height="$height"
     x="0"
     y="0"
     style="opacity:1;fill:white;fill-opacity:0;stroke:black;stroke-width:1;stroke-miterlimit:10;stroke-dasharray:none;stroke-opacity:1"
     id="rect001" />


FRAME
			
			
			my $y = 1 ;
			
			for my $taxon ( sort keys %$href_tax_states ) {
				
				push @taxa , $taxon ;
				
				my @characs = exists $href_tax_states->{$taxon} ? @{$href_tax_states->{$taxon}} : () ;
				
				my $x = 1 ;
				
				for my $char ( @characs ) {
					
					my $id = int rand 100000 ;
					
					if ( $char =~ /^1$/ ) {
					
print $fh_matrix_out <<RECT;

<rect
     width="10"
     height="10"
     x="$x"
     y="$y"
     style="fill:navy;fill-opacity:1;stroke:black;stroke-width:0.1;stroke-miterlimit:10;stroke-dasharray:none;stroke-opacity:1"
     id="rect$id" />
  
RECT
					$x +=10
					}
					
					elsif ( $char =~ /^0$/ ) {
					
print $fh_matrix_out <<RECT;

<rect
     width="10"
     height="10"
     x="$x"
     y="$y"
     style="fill:royalblue;fill-opacity:1;stroke:black;stroke-width:0.1;stroke-miterlimit:10;stroke-dasharray:none;stroke-opacity:1"
     id="rect$id" />
  
RECT
					$x +=10
					}
					
					elsif ( $char =~ /^X$/ ) {
					
print $fh_matrix_out <<RECT;

<rect
     width="10"
     height="10"
     x="$x"
     y="$y"
     style="fill:red;fill-opacity:1;stroke:black;stroke-width:0.1;stroke-miterlimit:10;stroke-dasharray:none;stroke-opacity:1"
     id="rect$id" />
  
RECT
					$x +=10
					}
					
					
					else{
print $fh_matrix_out <<RECT;
<rect
     width="10"
     height="10"
     x="$x"
     y="$y"
     style="fill:white;fill-opacity:1;stroke:black;stroke-width:0.1;stroke-miterlimit:10;stroke-dasharray:none;stroke-opacity:1"
     id="rect$id" />
  
RECT
					$x +=10
					}
				}
				$y += 10
			}
			
			
			$y = 10 ;

#write taxa to matrix

print $fh_matrix_out <<TAXA;	
 <text
     x="-5"
     y="10"
     style="font-size:8px;font-style:normal;font-weight:normal;text-align:end;text-anchor:end;fill:black;fill-opacity:1;stroke:none;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;font-family:Bitstream Vera Sans"
    id="text1892"
    
TAXA

	print $fh_matrix_out    "xml:space=\"preserve\">" ;


for my $taxon ( @taxa ) {
	
	$taxon =~ s/_/ /g ;
	
	my $id = int rand 100000 ;
	
	print $fh_matrix_out "<tspan\n       x=\"-5\"\n       y=\"$y\"\n       style=\"text-align:end;text-anchor:end\"\n       id=\"tspan$id\">$taxon</tspan>" ;
	
	$y += 10 ;
	
}
print $fh_matrix_out <<TAXA;	
</text>

TAXA

my $x = 0     ;

for my $bin ( @$aref_bin_numbers ) {  
	
	my $id = int rand 100000 ;
	
print $fh_matrix_out <<BINS;
<text
     x="$x"
     y="-5"
     transform="rotate(-90,$x,-5)"     	
     style="font-size:8px;font-style:normal;font-weight:normal;fill:black;fill-opacity:1;stroke:none;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;font-family:Bitstream Vera Sans"
     id="text8747"
     xml:space="preserve"><tspan
       x="$x"
       y="-5"
       id="tspan$id">$bin</tspan></text>
BINS
	
	$x += 10 ;
	
}




print $fh_matrix_out <<FINISH;

</svg>


FINISH
		}
	}
	
	
	sub ratio_error_bin{
		
		my $sref_info_name		= $_[0] ;
		my $sref_l_eq_ne			= $_[1] ;
		my $sref_N_bins			= $_[2] ;
		my $sref_error_number	= $_[3] ;
		my $href_high_ratio		= $_[4] ;
		 
		my $ratio					= $$sref_N_bins * $$sref_error_number ;
		
		if ( $href_high_ratio->{$$sref_l_eq_ne} ){ 
			
			my   @values = split ";;" 		, $href_high_ratio->{$$sref_l_eq_ne} ;
			if ( $ratio  > $values[1] )	{ $href_high_ratio->{$$sref_l_eq_ne} = $$sref_info_name.";;".$ratio } 
		}
		
		else 								{ $href_high_ratio->{$$sref_l_eq_ne} = $$sref_info_name.";;".$ratio }
	}
	
	
	sub print_out{
		
		my $sref_cuts		= $_[0] ;
		my $sref_lv_eq	= $_[1] ;
		my $aref_infile	= $_[2] ;
		my $sref_lv_tcut	= $_[3] ;
		
		my %seen_cut_bin	= () ;
		my %seen_cut_tax	= () ;
		my @remain			= () ;
		
		my $filename = "AMARE_matrix_".$$sref_lv_eq."_".$$sref_lv_tcut ;
		
		my @chars	= split ","  , $aref_infile->[0] ; my $N_chars  = @chars ;
		my @cut_b_t	= split ";;" , $$sref_cuts ;
		my @cut_bin	= split "::" , $cut_b_t[0] ;
		my @cut_tax	= split "::" , $cut_b_t[1] ;
		
		print LOGuser "\t$filename\n" ;
		
		for ( @cut_bin						){				$seen_cut_bin{$_}++	}
		for ( @cut_tax						){				$seen_cut_tax{$_}++	}
		for ( 1..($N_chars-1), my $i=0	){ unless (	$seen_cut_bin{$i}		){ push @remain, $i } $i++ }
		
		my $char_left		= @remain-1 ; 
		my $Ndeleted_dtax	= ( ( @cut_tax * 2 ) ) ;
		my $taxa_left		= @$aref_infile - $Ndeleted_dtax - 9 ;
		
		open  OUTmat, ">AMARE_Matrices/${filename}.txt" or warn "\nERROR: Can not print out /AMARE_Matrices/${filename}.txt\n" ;
		open  OUTnex, ">AMARE_Matrices/${filename}.nex" or warn "\nERROR: Can not print out /AMARE_Matrices/${filename}.nex\n" ;
		print OUTnex  "#NEXUS\nbegin data;\ndimensions ntax=$taxa_left nchar=$char_left ;\nformat missing=?;\nmatrix\n"        ;
		
		my $tax_pair = 0 ; my $line_count = 0 ;
		
		LINE:
		for my $line ( @$aref_infile ){ 
			
			my @bases = split ",", $line           ;
			
			# delete double samples which are rejected in main analysis
			if ( $seen_cut_tax{$bases[0]} ){ $tax_pair++ ; next LINE }
			if ( $tax_pair =~ 1           ){ $tax_pair-- ; next LINE }
			
			# for remaining taxa delete rejected bins
			my @final =  map { $bases[$_] } @remain ; 
			
			# Matrix output
			my $final_mat = ( join ",", @final )."\n"  ;
			print OUTmat "$final_mat" ;
			
			# NEXUS output
			my @blanks    =  () ; 
			my $length_b  =  20 - length $final[0] ;
			
			for ( 0 .. $length_b ){ push @blanks, ' ' } $length_b = join "" , @blanks ;
			
			$final[0] =~ s/$final[0]/$final[0]$length_b/ ;
			$final[0] =~ s/$final[0]/$final[0]\t/ ;
			$final[0] =~ s/"//g ;
			
			my $final_nex = ( join "" , @final )."\n"  ;
			
			if ( $line_count >= 9 ){ print OUTnex "$final_nex" }
			
			$line_count++
		}
		
		print OUTnex ";\nend;\n" ;
		close OUTmat ;
		close OUTnex ;
	}
	
	
	sub end{
		
		TIMER:
		# set timer
		my ( $user, $system, $cuser, $csystem ) = times ;
		
		print "\n\n\tAMARE 1.3beta.pl Finished !\n\n\tReduced AFLP-Matrices printed out as NEXUS- + TEXT-Files\n\tin ./AMARE_Matrices/\n\n\tSingle- + Main Log-Files printed out as TEXT-Files\n\tin ./AMARE_Logs/\n\n\tSVG-Plots of all matrices printed\n\tin ./AMARE_SVGs/\n\n\tThanks for using AMARE.\n\tGOOD BYE !\n\t--------------------------------------------\n\n";
		
		print <<TIME;
		
			***  time used: $user sec  ***
		
TIME
		

		
		
		exit ;
	}
}


#####################################################################################################################################################################
#####################################################################################################################################################################




