#!/usr/bin/perl -w

use Getopt::Long;# conda install bioconda::perl-getopt-long
use strict;
use POSIX qw(system);
use XML::Simple;#conda install bioconda::perl-xml-simple
use File::Basename;#conda install bioconda::perl-file-path

$| = 1;

my $configureXmlFile = "";
my $bedFile          = "";
my $output           = "";
my $prefix           = "motifFinding";
my $help             = "";
my $org              = "hg38";
my $unique           = 1;
my $isClass          = 0;
our $extendLen    = 25;
our %geneTypeHash = ();
our $motifType    = "geneType";

if ( scalar(@ARGV) < 2 ) {
    usage();
}

my $optLong = GetOptions(
    "conf=s"   => \$configureXmlFile,    # configure file
    "bed=s"    => \$bedFile,             # mod file
    "pre=s"    => \$prefix,              # pre
    "org=s"    => \$org,                 # org
    "class=s"  => \$isClass,             # isClass
    "type=s"   => \$motifType,           # motifType
    "output=s" => \$output,              # output file
    "unique|u" => \$unique,              # keep unique ncRNAs
    "help|h"   => \$help
) or usage();

&usage() if $help;

sub usage
### usage function
{
    print
"Usage: perl $0 [options] --conf <configure file> --bed <bed from bedAnnotator> --output <output dir> --pre <prefix name>\n";
    print "--conf     conf file\n";
    print "--bed      bed from bedAnnotator\n";
    print "--pre      prefix name\n";
    print "--class    classify by the genomic regions\n";
    print
"--type     the type for motif enrichment[default=geneType], e.g. tsiType, shape\n";
    print "--output   output directory\n";
    print "--unique   get unique seqence with cdhit software\n";
    print "--org      organism, e.g. hg38\n";
    print "--help     help information\n";
    exit;
}

if ( $bedFile eq "" ) {
    print "Please use --bed option\n";
    usage();
}

if ( $output eq "" ) {
    print "Please use --output option\n";
    usage();
}

if ( $prefix eq "" ) {
    print "Please use --pre option\n";
    usage();
}

if ( $org eq "" ) {
    print "Please use --org option\n";
    usage();
}

if ( !( -e $output ) ) {
    &system("mkdir -p $output");
}

# $configureXmlFile = "motifConfigure_hg38.xml";
my $xmlConfHash = XMLin($configureXmlFile);

if ($unique) {
    $bedFile = &getUniqNcrna( $xmlConfHash, $bedFile, $output, $prefix, $org );
}


&motifFinder( $xmlConfHash, $bedFile, $output, $prefix, $org, $isClass,
    $motifType );

sub getUniqNcrna {
    my ( $confHash, $ncrnaFile, $outputDir, $prefix, $orgName ) = @_;
    my %rmHash     = ();
    my $ncrnaFasta = $outputDir . "/" . $prefix . "_gencode.ncRNA.fa";
    my $ncrnaFasta_cleaned = $outputDir . "/" . $prefix . "_gencode.ncRNA.fa.cleaned";
    my $cdhitFile  = $outputDir . "/" . $prefix . "_gencode.ncRNA.cdhit";
    my $genome     = $confHash->{'genomeFile'};

    my $extendFile = $outputDir . "/" . $prefix . ".ex10nt.ncRNA.bed";

    open( NCRNA,  "<$ncrnaFile" )  || die "can't open the $ncrnaFile\n";
    open( EXTEND, ">$extendFile" ) || die "can't open the $extendFile\n";
    while ( my $line = <NCRNA> ) {
        $line =~ s/\s+$//;
        next if ( $line =~ /\#/ );
        my @items      = split /\t/, $line;
        my $name       = $items[3];
        my $chromStart = $items[1];
        my $chromEnd   = $items[2];
        if ( $chromEnd - $chromStart < 20 ) {
            $chromStart = $chromStart - 10;
            $chromEnd   = $chromEnd + 10;
            $chromStart = 0 if ( $chromStart < 0 );
        }
        $items[1] = $chromStart;
        $items[2] = $chromEnd;
        print EXTEND join( "\t", @items ), "\n";
    }
    close(NCRNA);
    close(EXTEND);

    my $commandLine =
        $confHash->{'bedtools'} . " "
      . " getfasta -s -name -fi "
      . $confHash->{'genomeFile'}
      . " -bed "
      . $extendFile
      . " -fo "
      . $ncrnaFasta;
    print "$commandLine\n";
    &system($commandLine);

    open(my $FASTA, '<', $ncrnaFasta) or die "Cannot open FASTA file: $!";
    open(my $OUTPUT, '>', "$ncrnaFasta.cleaned") or die "Cannot open output file: $!";

    while (<$FASTA>) {
      if (/>/) {
        # header line
        my $header = $_;
        $header =~ s/::.*//;
        print $OUTPUT $header; 
      } else {
        # sequence line 
        print $OUTPUT $_;
      }
    }

    close($FASTA);
    close($OUTPUT);

    $commandLine =
        $confHash->{'cdHit'} . " "
      . $confHash->{'cdHitParameter'} . " -i "
      . $ncrnaFasta_cleaned . " -o "
      . $cdhitFile;
    print "$commandLine\n";
    &system($commandLine);

    my $clstrFile = $cdhitFile . ".clstr";
    open( CLUSTER, "<$clstrFile" ) || die "can't open the $clstrFile\n";

    while ( my $line = <CLUSTER> ) {
        $line =~ s/\s+$//;
        if ( $line =~ /at\s+\d+/ ) {
            my @items = split /\s+/, $line;
            my $id    = $items[2];
            $id =~ s/\>//;
            $id =~ s/\.+$//;
            $rmHash{$id} = $items[2];
        }
    }
    close(CLUSTER);

    my $rmCdHitFile =
      $outputDir . "/" . $prefix . ".identification.allTypeNcrna.cdhit.bed";
    open( NCRNA, "<$ncrnaFile" )   || die "can't open the $ncrnaFile\n";
    open( CDHIT, ">$rmCdHitFile" ) || die "can't open the $cdhitFile\n";
    while ( my $line = <NCRNA> ) {
        $line =~ s/\s+$//;
        my @items = split /\t/, $line;
        my $name  = $items[3];
        if ( defined( $rmHash{$name} ) ) {
            next;
        }
        print CDHIT $line, "\n";
    }
    close(CDHIT);
    close(NCRNA);
    return $rmCdHitFile;
}

sub motifFinder {
    my ( $confHash, $contigFile, $createDir, $prefix, $org, $isClass,
        $motifType )
      = @_;
    my %targetHash    = ();
    my %upRnaHash     = ();
    my %downRnaHash   = ();
    my %seqLenHash    = ();
    my %ncrnaHash     = ();
    my %ncrnaUpHash   = ();
    my %ncrnaDownHash = ();
    my $totalGeneNum  = 0;
    my $minGeneNum    = 30;

    #my $extendLen     = 20;
    my $extendSiteLen = 10;
    my $geneTypeIdx   = -1;
    my $shapeIdx      = -1;
    my $tsiIdx        = -1;
    my $siteNum       = 0;
    open( CONTIGRNA, "<$contigFile" ) || die "can't open the $contigFile\n";
    while ( my $line = <CONTIGRNA> ) {
        $line =~ s/\s+$//;
        my @items = split /\s+/, $line;
        if ( $line =~ /\#/ ) {
            for ( my $i = 0 ; $i < scalar(@items) ; $i++ ) {
                if ( $items[$i] eq "shape" && $motifType eq "shape" ) {
                    $shapeIdx = $i;
                }
                if ( $items[$i] eq "tsiType" && $motifType eq "tsiType" ) {
                    $tsiIdx = $i;
                }
                if ( $items[$i] eq "geneType" && $motifType eq "geneType" ) {
                    $geneTypeIdx = $i;
                }
            }
            next;
        }
        if ( scalar(@items) < 2 ) {
            warn "short line: $line\n";
            next;
        }
        my $chrom      = $items[0];
        my $chromStart = $items[1];
        my $chromEnd   = $items[2];
        my $name       = $items[3];
        my $score      = $items[4];

        if ( $score < 2 ) {

            # next;
        }

        my $strand    = $items[5];
        my $seq       = $items[24];
        my $geneType  = $items[$geneTypeIdx];
        my $shapeType = $items[$shapeIdx];
        my $tsiType   = $items[$tsiIdx];
        if ( defined($geneType) ) {
            $geneTypeHash{$name} = $geneType;
        }
        else {
            $geneTypeHash{$name} = "unknown";
        }

        my $seqLen = $chromEnd - $chromStart;

        my $upStart = ( $chromStart - $extendLen );
        my $upEnd   = ( $chromStart + $extendLen + 1 );
        $upStart = 0 if ( $upStart < 0 );
        my $upRNA =
            $chrom . "\t"
          . $upStart . "\t"
          . $upEnd . "\t"
          . $name . "\t"
          . $score . "\t"
          . $strand;
        my $downStart = ( $chromEnd - $extendLen - 1 );
        my $downEnd   = ( $chromEnd + $extendLen );
        $downStart = 0 if ( $downStart < 0 );
        my $downRNA =
            $chrom . "\t"
          . $downStart . "\t"
          . $downEnd . "\t"
          . $name . "\t"
          . $score . "\t"
          . $strand;

        if ( $strand eq "-" ) {
            $downRNA =
                $chrom . "\t"
              . $upStart . "\t"
              . $upEnd . "\t"
              . $name . "\t"
              . $score . "\t"
              . $strand;
            $upRNA =
                $chrom . "\t"
              . $downStart . "\t"
              . $downEnd . "\t"
              . $name . "\t"
              . $score . "\t"
              . $strand;
        }
        if ( $seqLen < 10 ) {
            $chromStart = $chromStart - $extendSiteLen;
            $chromEnd   = $chromEnd + $extendSiteLen;
            $chromStart = 0 if ( $chromStart < 0 );
        }
        my $ncRNA =
            $chrom . "\t"
          . $chromStart . "\t"
          . $chromEnd . "\t"
          . $name . "\t"
          . $score . "\t"
          . $strand;

        $totalGeneNum++;

        $ncrnaHash{$ncRNA} = $name;

        if ( $upStart != 0 && $downStart != 0 ) {
            $ncrnaUpHash{$upRNA}     = $name;
            $ncrnaDownHash{$downRNA} = $name;
        }

        $seqLenHash{$name} = $seqLen;

        my @classTypes = ();
        if ( $geneTypeIdx > 0 ) {
            push @classTypes, $geneType;
        }
        if ( $tsiIdx > 0 ) {
            push @classTypes, $tsiType;
        }
        if ( $shapeIdx > 0 ) {
            push @classTypes, $shapeType;
        }

        #@classTypes = ($geneType);    # discard classification of shape
        foreach my $classType (@classTypes) {
            if ( $classType =~ /\[/ ) {
                if ( length($classType) < 3 ) {
                    next;
                }
            }
            if ( $classType eq "." ) {
                next;
            }
            $classType =~ s/\[/A/g;
            $classType =~ s/\]/B/g;
            if ( defined( $targetHash{$classType} ) ) {
                push @{ $targetHash{$classType} }, $ncRNA;
                if ( $upStart != 0 && $downStart != 0 ) {
                    push @{ $upRnaHash{$classType} },   $upRNA;
                    push @{ $downRnaHash{$classType} }, $downRNA;
                }
            }
            else {
                $targetHash{$classType}  = [];
                $upRnaHash{$classType}   = [];
                $downRnaHash{$classType} = [];
                push @{ $targetHash{$classType} }, $ncRNA;

                if ( $upStart != 0 && $downStart != 0 ) {
                    push @{ $upRnaHash{$classType} },   $upRNA;
                    push @{ $downRnaHash{$classType} }, $downRNA;
                }
            }
        }
    }    ### while
    close(CONTIGRNA);

    my $allMotifSiteFile = $createDir . "/" . "allType.motif.sites.txt";
    my $outfp;
    open( $outfp, ">$allMotifSiteFile" )
      || die "can't open the $allMotifSiteFile\n";
    print $outfp
"#typeName\tmotifName\tmotifSeq\tdistType\tdist\tnegLg(pval)\tPercentage\tfoldChange\tdistTargetMotifNum\ttargetMotifNum\tdistBackgroundMotifNum\tbackgroundMotifNum\tmotifGeneType\n";
    my $motifDir     = "";
    my $sigMotifFile = "";

    #my $typeNum = scalar(keys %targetHash);
    if ($isClass) {    ### classify by genomic regions
        foreach my $type ( sort keys %targetHash ) {
            my @ncRNAs = @{ $targetHash{$type} };
            my $seqNum = scalar(@ncRNAs);
            if ( $seqNum >= $minGeneNum ) {
                warn "The number of geneType-" . $type, " is ", $seqNum, "\n";
                my $fullType = $type . "_gene";
                &generateHomerTypeMotifs(
                    $confHash, $createDir, $fullType,
                    $prefix,   $org,       \@ncRNAs
                );
                $motifDir =
                  $createDir . "/" . $fullType . "_HomerGenomeResults";
                &outputMotifToOneFile( $outfp, $motifDir, $prefix, $fullType,
                    "homer_" . $fullType );
            }
        }
        foreach my $type ( sort keys %upRnaHash ) {
            my @ncRNAs = @{ $upRnaHash{$type} };
            my $seqNum = scalar(@ncRNAs);
            if ( $seqNum >= $minGeneNum ) {
                warn "The number of upstream-" . $extendLen . "-" . $type,
                  " is ", $seqNum, "\n";
                my $fullType = $type . "_upstream" . $extendLen . "nt";
                &generateHomerTypeMotifs(
                    $confHash, $createDir, $fullType,
                    $prefix,   $org,       \@ncRNAs
                );
                $motifDir =
                  $createDir . "/" . $fullType . "_HomerGenomeResults";
                &outputMotifToOneFile( $outfp, $motifDir, $prefix, $fullType,
                    "homer_" . $fullType );
            }
        }
        foreach my $type ( sort keys %downRnaHash ) {
            my @ncRNAs = @{ $downRnaHash{$type} };
            my $seqNum = scalar(@ncRNAs);
            if ( $seqNum >= $minGeneNum ) {
                warn "The number of downstream-" . $extendLen . "-" . $type,
                  " is ", $seqNum, "\n";
                my $fullType = $type . "_downstream" . $extendLen . "nt";
                &generateHomerTypeMotifs(
                    $confHash, $createDir, $fullType,
                    $prefix,   $org,       \@ncRNAs
                );
                $motifDir =
                  $createDir . "/" . $fullType . "_HomerGenomeResults";
                &outputMotifToOneFile( $outfp, $motifDir, $prefix, $fullType,
                    "homer_" . $fullType );
            }
        }
    }
    else {
        my $type   = 'full';
        my @ncRNAs = sort keys %ncrnaHash;
        my $seqNum = scalar(@ncRNAs);
        if ( $seqNum >= $minGeneNum ) {
            warn "The number of geneType-" . $type, " is ", $seqNum, "\n";
            my $fullType = $type . "_gene";
            &generateHomerTypeMotifs( $confHash, $createDir, $fullType,
                $prefix, $org, \@ncRNAs );
            $motifDir = $createDir . "/" . $fullType . "_HomerGenomeResults";
            &outputMotifToOneFile( $outfp, $motifDir, $prefix, $fullType,
                "homer_" . $fullType );
        }
        @ncRNAs = sort keys %ncrnaUpHash;
        $seqNum = scalar(@ncRNAs);
        if ( $seqNum >= $minGeneNum ) {
            warn "The number of upstream-" . $extendLen . "-" . $type,
              " is ", $seqNum, "\n";
            my $fullType = $type . "_upstream" . $extendLen . "nt";
            &generateHomerTypeMotifs( $confHash, $createDir, $fullType,
                $prefix, $org, \@ncRNAs );
            $motifDir = $createDir . "/" . $fullType . "_HomerGenomeResults";
            &outputMotifToOneFile( $outfp, $motifDir, $prefix, $fullType,
                "homer_" . $fullType );
        }
        @ncRNAs = sort keys %ncrnaDownHash;
        $seqNum = scalar(@ncRNAs);
        if ( $seqNum >= $minGeneNum ) {
            warn "The number of downstream-" . $extendLen . "-" . $type,
              " is ", $seqNum, "\n";
            my $fullType = $type . "_downstream" . $extendLen . "nt";
            &generateHomerTypeMotifs( $confHash, $createDir, $fullType,
                $prefix, $org, \@ncRNAs );
            $motifDir = $createDir . "/" . $fullType . "_HomerGenomeResults";
            &outputMotifToOneFile( $outfp, $motifDir, $prefix, $fullType,
                "homer_" . $fullType );
        }
    }    # else end
    close($outfp);
}

sub outputMotifToOneFile {
    my ( $outfp, $motifDir, $prefix, $fullType, $typeName ) = @_;
    my $motifFile =
        $motifDir . "/"
      . $prefix . "."
      . $fullType
      . ".significant.positionSpecifc.motifs.txt";

    my $motifSeqFile =
        $motifDir . "/"
      . $prefix . "."
      . $fullType
      . ".significant.positionSpecifc.motifs.sequences.txt";

    my %seqTypeHash = ();
    if ( -e $motifSeqFile ) {
        open( MOTIFSEQ, "<$motifSeqFile" )
          || die "can't open the $motifSeqFile\n";
        while ( my $line = <MOTIFSEQ> ) {
            $line =~ s/\s+$//;
            next if ( $line =~ /\#/ );
            my @items     = split /\t/, $line;
            my $motifName = $items[0] . ":" . $items[2] . ":" . $items[3];
            my $geneType  = $items[-1];
            $seqTypeHash{$motifName}{$geneType}++;
        }
        close(MOTIFSEQ);
    }

    my %motifGeneTypeHash = ();
    foreach my $motifName ( keys %seqTypeHash ) {
        my %geneTypeHash = %{ $seqTypeHash{$motifName} };
        my $sum          = 0;
        foreach my $geneType ( keys %geneTypeHash ) {
            $sum += $geneTypeHash{$geneType};
        }
        my $geneTypeStr = "";
        foreach my $geneType (
            sort { $geneTypeHash{$b} <=> $geneTypeHash{$a} }
            keys %geneTypeHash
          )
        {
            my $num = $geneTypeHash{$geneType};
            my $per = sprintf( "%.2f", $num / $sum * 100 );
            if ( $geneTypeStr eq "" ) {
                $geneTypeStr = $geneType . ":" . $per . ":" . $num;
            }
            else {
                $geneTypeStr =
                  $geneTypeStr . "#" . $geneType . ":" . $per . ":" . $num;
            }
        }
        $motifGeneTypeHash{$motifName} = $geneTypeStr;
    }

    if ( -e $motifFile ) {
        open( MOTIF, "<$motifFile" ) || die "can't open the $motifFile\n";
        while ( my $line = <MOTIF> ) {
            $line =~ s/\s+$//;
            next if ( $line =~ /\#/ );
            my @items     = split /\t/, $line;
            my $motifName = $items[0] . ":" . $items[2] . ":" . $items[3];
            print $outfp $typeName, "\t", $line, "\t",
              $motifGeneTypeHash{$motifName}, "\n";
        }
        close(MOTIF);
    }
}

sub generateHomerTypeMotifs {
    my ( $confHash, $createDir, $type, $prefix, $org, $dataSet ) = @_;

    my $motifDir = $createDir . "/" . $type . "_HomerGenomeResults";
    if ( !( -e $motifDir ) ) {
        &system("mkdir -p $motifDir");
    }
    my $targetFile = $motifDir . "/" . $type . "_contigs.bed6";
    open( TARGET, ">$targetFile" )
      || die "can't open the $targetFile\n";
    foreach my $bedLine ( @{$dataSet} ) {
        print TARGET $bedLine, "\n";
    }
    close(TARGET);
    &findHomerMotifs( $confHash, $motifDir, $type, $prefix, $org, $targetFile );
}

sub findHomerMotifs {
    my ( $confHash, $motifDir, $type, $prefix, $org, $targetFile ) = @_;
    my $logFile = $motifDir . "/" . $prefix . "." . $type . ".log.txt";
    my $commandLine = $confHash->{'homer'}." "
      . $targetFile . " "
      . $confHash->{'genomeFile'} . " "
      . $motifDir . " "
      . $confHash->{'homerParameter'}
      . " 2>$logFile";
      # . " -norevopp -noknown -rna -len 4,5,6,7,8,9 -p 20 -size given -dumpFasta 2>$logFile";
    print "$commandLine\n";
    &system($commandLine);
    &findHomerMotifsites( $confHash, $motifDir, $prefix, $type );
    &findAllPositionMotifSites( $confHash, $motifDir, $prefix, $type );
}


sub findHomerMotifsites {
    my ( $confHash, $motifDir, $prefix, $type ) = @_;
    my $tgFile     = $motifDir . "/" . "target.fa";
    my $bgFile     = $motifDir . "/" . "background.fa";
    my $motifNum   = 35;
    my $similarNum = 20;

    my $sigMotifFile =
        $motifDir . "/"
      . $prefix . "."
      . $type
      . ".significant.positionSpecifc.motifs.txt";
    my $sigMotifSeqFile =
        $motifDir . "/"
      . $prefix . "."
      . $type
      . ".significant.positionSpecifc.motifs.sequences.txt";
    my $logFile =
        $motifDir . "/"
      . $prefix . "."
      . $type
      . ".significant.positionSpecifc.motifs.logs.txt";

    open( SIGMOTIF, ">$sigMotifFile" )
      || die "can't open the $sigMotifFile\n";
    open( SEQMOTIF, ">$sigMotifSeqFile" )
      || die "can't open the $sigMotifSeqFile\n";
    print SIGMOTIF
"#motifName\tmotifSeq\tdistType\tdist\tnegLg(pval)\tPercentage\tfoldChange\tdistTargetMotifNum\ttargetMotifNum\tdistBackgroundMotifNum\tbackgroundMotifNum\n";
    print SEQMOTIF
"#motifName\tmotifSeq\tdistType\tdist\tmotifInSeq\tscore\tsequence\tseqName\tseqType\n";

    for ( my $i = 1 ; $i < $motifNum ; $i++ ) {
        my $motifName = "motif" . $i . ".motif";
        my $motifFile = $motifDir . "/homerResults/" . $motifName;
        my $outfile =
          $motifDir . "/" . $type . "." . $motifName . ".significant.txt";
        my $seqfile =
          $motifDir . "/" . $type . "." . $motifName . ".significant.seq.txt";
        if ( -e $motifFile ) {
            my $commandLine = $confHash->{'findMotif'}.
            " -p -10 -f 2 -s 0.85 -r 0.03 -n 10 -t $tgFile -b $bgFile -m $motifFile -o $outfile -e $seqfile 2>$logFile";
            print "$commandLine\n";
            &system($commandLine);
            open( MOTIF, "<$outfile" ) || die "can't open the $outfile\n";
            while ( my $line = <MOTIF> ) {
                $line =~ s/\s+$//;
                if ( $line =~ /\#/ ) {
                    next;
                }
                else {
                    print SIGMOTIF $motifName, "\t", $line, "\n";
                }
            }
            close(MOTIF);
            open( SEQ, "<$seqfile" ) || die "can't open the $seqfile\n";
            while ( my $line = <SEQ> ) {
                $line =~ s/\s+$//;
                my @items   = split /\t/, $line;
                my $seqName = $items[$#items];
                my $seqType = "unknown";
                if ( defined( $geneTypeHash{$seqName} ) ) {
                    $seqType = $geneTypeHash{$seqName};
                }
                if ( $line =~ /\#/ ) {
                    next;
                }
                else {
                    print SEQMOTIF $motifName, "\t", $line, "\t",
                      $seqType, "\n";
                }
            }
            close(SEQ);
            $commandLine = "rm -f $outfile $seqfile";
            &system($commandLine);
        }
        for ( my $j = 1 ; $j < $similarNum ; $j++ ) {
            $motifName = "motif" . $i . ".similar" . $j . ".motif";
            $motifFile = $motifDir . "/homerResults/" . $motifName;
            $outfile =
              $motifDir . "/" . $type . "." . $motifName . ".significant.txt";
            $seqfile =
                $motifDir . "/"
              . $type . "."
              . $motifName
              . ".significant.seq.txt";
            if ( -e $motifFile ) {
                my $commandLine = $confHash->{'findMotif'}.
                " -p -10 -f 2 -s 0.85 -r 0.03 -n 10 -t $tgFile -b $bgFile -m $motifFile -o $outfile -e $seqfile 2>$logFile";
                print "$commandLine\n";
                &system($commandLine);
                open( MOTIF, "<$outfile" ) || die "can't open the $outfile\n";
                while ( my $line = <MOTIF> ) {
                    $line =~ s/\s+$//;
                    if ( $line =~ /\#/ ) {
                        next;
                    }
                    else {
                        print SIGMOTIF $motifName, "\t", $line, "\n";
                    }
                }
                close(MOTIF);
                open( SEQ, "<$seqfile" ) || die "can't open the $seqfile\n";
                while ( my $line = <SEQ> ) {
                    $line =~ s/\s+$//;
                    my @items   = split /\t/, $line;
                    my $seqName = $items[$#items];
                    my $seqType = "unknown";
                    if ( defined( $geneTypeHash{$seqName} ) ) {
                        $seqType = $geneTypeHash{$seqName};
                    }
                    if ( $line =~ /\#/ ) {
                        next;
                    }
                    else {
                        print SEQMOTIF $motifName, "\t", $line, "\t",
                          $seqType, "\n";
                    }
                }
                close(SEQ);
                $commandLine = "rm -f $outfile $seqfile";
                &system($commandLine);
            }
        }
    }
    close(SIGMOTIF);
    close(SEQMOTIF);
}

sub findAllPositionMotifSites {
    my ( $confHash, $motifDir, $prefix, $type ) = @_;
    my $tgFile     = $motifDir . "/" . "target.fa";
    my $bgFile     = $motifDir . "/" . "background.fa";
    my $motifNum   = 35;
    my $similarNum = 20;

    my $sigMotifFile =
        $motifDir . "/"
      . $prefix . "."
      . $type
      . ".significant.Allpositions.motifs.txt";
    my $sigMotifSeqFile =
        $motifDir . "/"
      . $prefix . "."
      . $type
      . ".significant.Allpositions.motifs.sequences.txt";
    my $logFile =
        $motifDir . "/"
      . $prefix . "."
      . $type
      . ".significant.Allpositions.motifs.logs.txt";

    open( SIGMOTIF, ">$sigMotifFile" )
      || die "can't open the $sigMotifFile\n";
    open( SEQMOTIF, ">$sigMotifSeqFile" )
      || die "can't open the $sigMotifSeqFile\n";
    print SIGMOTIF
"#motifName\tmotifSeq\tdistType\tdist\tnegLg(pval)\tPercentage\tfoldChange\tdistTargetMotifNum\ttargetMotifNum\tdistBackgroundMotifNum\tbackgroundMotifNum\n";
    print SEQMOTIF
"#motifName\tmotifSeq\tdistType\tdist\tmotifInSeq\tscore\tsequence\tseqName\tseqType\n";

    for ( my $i = 1 ; $i < $motifNum ; $i++ ) {
        my $motifName = "motif" . $i . ".motif";
        my $motifFile = $motifDir . "/homerResults/" . $motifName;
        my $outfile =
          $motifDir . "/" . $type . "." . $motifName . ".all.significant.txt";
        my $seqfile =
            $motifDir . "/"
          . $type . "."
          . $motifName
          . ".all.significant.seq.txt";
        if ( -e $motifFile ) {
            my $commandLine =$confHash->{'findMotif'}.
            " --all -p -1.0 -f 1.0 -s 0.85 -r 0.005 -n 10 -t $tgFile -b $bgFile -m $motifFile -o $outfile -e $seqfile 2>$logFile";
            print "$commandLine\n";
            &system($commandLine);
            open( MOTIF, "<$outfile" ) || die "can't open the $outfile\n";
            while ( my $line = <MOTIF> ) {
                $line =~ s/\s+$//;
                if ( $line =~ /\#/ ) {
                    next;
                }
                else {
                    print SIGMOTIF $motifName, "\t", $line, "\n";
                }
            }
            close(MOTIF);
            open( SEQ, "<$seqfile" ) || die "can't open the $seqfile\n";
            while ( my $line = <SEQ> ) {
                $line =~ s/\s+$//;
                my @items   = split /\t/, $line;
                my $seqName = $items[$#items];
                my $seqType = "unknown";
                if ( defined( $geneTypeHash{$seqName} ) ) {
                    $seqType = $geneTypeHash{$seqName};
                }
                if ( $line =~ /\#/ ) {
                    next;
                }
                else {
                    print SEQMOTIF $motifName, "\t", $line, "\t",
                      $seqType, "\n";
                }
            }
            close(SEQ);
            $commandLine = "rm -f $outfile $seqfile";
            &system($commandLine);
        }
        for ( my $j = 1 ; $j < $similarNum ; $j++ ) {
            $motifName = "motif" . $i . ".similar" . $j . ".motif";
            $motifFile = $motifDir . "/homerResults/" . $motifName;
            $outfile =
                $motifDir . "/"
              . $type . "."
              . $motifName
              . ".all.significant.txt";
            $seqfile =
                $motifDir . "/"
              . $type . "."
              . $motifName
              . ".all.significant.seq.txt";
            if ( -e $motifFile ) {
                my $commandLine =$confHash->{'findMotif'}.
                " --all -p -1.0 -f 1.0 -s 0.85 -r 0.005 -n 10 -t $tgFile -b $bgFile -m $motifFile -o $outfile -e $seqfile 2>$logFile";
                print "$commandLine\n";
                &system($commandLine);
                open( MOTIF, "<$outfile" ) || die "can't open the $outfile\n";
                while ( my $line = <MOTIF> ) {
                    $line =~ s/\s+$//;
                    if ( $line =~ /\#/ ) {
                        next;
                    }
                    else {
                        print SIGMOTIF $motifName, "\t", $line, "\n";
                    }
                }
                close(MOTIF);
                open( SEQ, "<$seqfile" ) || die "can't open the $seqfile\n";
                while ( my $line = <SEQ> ) {
                    $line =~ s/\s+$//;
                    my @items   = split /\t/, $line;
                    my $seqName = $items[$#items];
                    my $seqType = "unknown";
                    if ( defined( $geneTypeHash{$seqName} ) ) {
                        $seqType = $geneTypeHash{$seqName};
                    }
                    if ( $line =~ /\#/ ) {
                        next;
                    }
                    else {
                        print SEQMOTIF $motifName, "\t", $line, "\t",
                          $seqType, "\n";
                    }
                }
                close(SEQ);
                $commandLine = "rm -f $outfile $seqfile";
                &system($commandLine);
            }
        }
    }
    close(SIGMOTIF);
    close(SEQMOTIF);
}

