#!/usr/bin/env perl
#
#  FASTA Splitter  -  a script for partitioning a FASTA file into pieces
#
#  Version 0.2.6 (August 1, 2017)
#
#  Copyright (c) 2012-2017 Kirill Kryukov
#
#  This software is provided 'as-is', without any express or implied
#  warranty. In no event will the authors be held liable for any damages
#  arising from the use of this software.
#
#  Permission is granted to anyone to use this software for any purpose,
#  including commercial applications, and to alter it and redistribute it
#  freely, subject to the following restrictions:
#
#  1. The origin of this software must not be misrepresented; you must not
#     claim that you wrote the original software. If you use this software
#     in a product, an acknowledgment in the product documentation would be
#     appreciated but is not required.
#  2. Altered source versions must be plainly marked as such, and must not be
#     misrepresented as being the original software.
#  3. This notice may not be removed or altered from any source distribution.
#

use strict;
use File::Basename qw(basename);
use File::Path qw(make_path);
use Getopt::Long qw(:config pass_through);

$| = 1;

my ($script_name,$script_version,$script_date,$script_years) = ('fasta-splitter','0.2.6','2017-08-01','2012-2017');
my $start_time = time;

my @files = ();
my ($opt_n_parts,$opt_part_size,$opt_part_num_prefix,$opt_measure,$opt_line_len,$opt_eol,$out_dir,$nopad,$ver,$help);
GetOptions('n-parts=i'         => \$opt_n_parts,
           'part-size=i'       => \$opt_part_size,
           'part-num-prefix=s' => \$opt_part_num_prefix,
           'measure=s'         => \$opt_measure,
           'line-length=i'     => \$opt_line_len,
           'eol=s'             => \$opt_eol,
           'out-dir=s'         => \$out_dir,
           'nopad'             => \$nopad,
           'version'           => \$ver,
           'help'              => \$help);
for (my $i=0; $i<scalar(@ARGV); $i++)
{
    if (substr($ARGV[$i],0,1) eq '-' and $i < scalar(@ARGV)-1)
    {
        if ($ARGV[$i] eq '-n-parts-total'     ) { $opt_n_parts   = int($ARGV[++$i]); $opt_measure = 'all'; }
        if ($ARGV[$i] eq '-n-parts-sequence'  ) { $opt_n_parts   = int($ARGV[++$i]); $opt_measure = 'seq'; }
        if ($ARGV[$i] eq '-part-total-size'   ) { $opt_part_size = int($ARGV[++$i]); $opt_measure = 'all'; }
        if ($ARGV[$i] eq '-part-sequence-size') { $opt_part_size = int($ARGV[++$i]); $opt_measure = 'seq'; }
        if ($ARGV[$i] eq '-line-length') { $opt_line_len = int($ARGV[++$i]); }
        if ($ARGV[$i] eq '-eol'        ) { $opt_eol      = int($ARGV[++$i]); }
    }
    else { push @files, $ARGV[$i]; }
}

my $ver_str = "$script_name, version $script_version, $script_date\nCopyright (c) $script_years Kirill Kryukov\n";
my $help_str = qq{Usage: ${script_name} [options] <file>...
Options:
    --n-parts <N>        - Divide into <N> parts
    --part-size <N>      - Divide into parts of size <N>
    --measure (all|seq|count) - Specify whether all data, sequence length, or
                           number of sequences is used for determining part
                           sizes ('all' by default).
    --line-length        - Set output sequence line length, 0 for single line
                           (default: 60).
    --eol (dos|mac|unix) - Choose end-of-line character ('unix' by default).
    --part-num-prefix T  - Put T before part number in file names (def.: .part-)
    --out-dir            - Specify output directory.
    --nopad              - Don't pad part numbers with 0.
    --version            - Show version.
    --help               - Show help.
};

print (($ver ? $ver_str : ''), ($help ? $help_str : ''));
if (!defined($opt_n_parts) and !defined($opt_part_size) and !defined($opt_measure) and !defined($opt_line_len) and !defined($opt_eol))
{
    if (!$help and !$ver) { print $ver_str, $help_str; } exit;
}
if (!defined($opt_n_parts) and !defined($opt_part_size)) { die "Splitting method is not specified\nUse -h for help\n"; }
if (!@files) { die "File for splitting is not specified\n"; }

if (defined($opt_n_parts) and $opt_n_parts <= 0) { die "Non-positive number of parts\n"; }
if (defined($opt_part_size) and $opt_part_size <= 0) { die "Non-positive part size\n"; }
if (defined($opt_measure) and $opt_measure ne 'all' and $opt_measure ne 'seq' and $opt_measure ne 'count') { die "Unknown value of --measure option\n"; }
if (defined($opt_eol) and $opt_eol ne 'dos' and $opt_eol ne 'mac' and $opt_eol ne 'unix') { die "Unknown value of --eol option\n"; }
if (defined($out_dir))
{
    $out_dir =~ s/[\/\\]+$//;
    if (!-e $out_dir) { make_path($out_dir); }
    if (!-e $out_dir || !-d $out_dir) { die "Can't create output directory \"$out_dir\"\n"; }
    $out_dir .= '/';
}

my $n_parts = defined($opt_n_parts) ? $opt_n_parts : 0;
my $part_size = defined($opt_part_size) ? $opt_part_size : 0;
my $line_len = (defined($opt_line_len) and $opt_line_len >= 0) ? $opt_line_len : 60;
my $eol = defined($opt_eol) ? (($opt_eol eq 'dos') ? "\x0D\x0A" : ($opt_eol eq 'mac') ? "\x0D" : "\x0A") : "\x0A";
my $eol_len = length($eol);
my $measure = defined($opt_measure) ? (($opt_measure eq 'count') ? 0 : ($opt_measure eq 'seq') ? 1 : 2) : 2;
my $part_num_prefix = defined($opt_part_num_prefix) ? $opt_part_num_prefix : '.part-';
my @part_start = ();
my ($base,$ext,$num_len,$total_size);
my ($OUT,$name,$data,$written_total,$written_this_part,$part_end,$part);

foreach my $infile (@files) { split_file($infile); }

my $elapsed_time = time - $start_time;
print "All done, $elapsed_time second", (($elapsed_time==1)?'':'s'), " elapsed\n";

sub split_file
{
    my ($infile) = @_;
    if (!-e $infile or !-f $infile) { print "Can't find file \"$infile\"\n"; return; }
    print $infile;

    ($base,$ext) = (basename($infile),'');
    if ($base =~ /^(.+?)(\.[^\.]+)$/) { ($base,$ext) = ($1,$2); }

    @part_start = ();
    my ($n_seq,$total_seq_len,$n_parts_found) = (0,0,0);

    if ($part_size)
    {
        ($n_seq,$total_seq_len,$total_size,$n_parts_found) = get_file_size_and_part_boundaries($infile);
        if (!$n_parts) { print ": $n_seq sequences, $total_seq_len bp"; }
        print ' => ', ($n_parts ? 'extracting' : 'dividing into'), ' ', $n_parts_found, ' part', ($n_parts_found > 1 ? 's' : ''),
              " of <= $part_size ", ($measure ? (($measure > 1) ? 'bytes' : 'bp') : 'sequences'), "\n";
        open(my $IN,'<',$infile) or die "Error: Can't open file \"$infile\"\n";
        binmode $IN;
        $num_len = length($n_parts_found);
        $OUT = undef;
        my ($out_file,$part,$si,$buffer) = (undef,0,-1,'');
        while (<$IN>)
        {
            $_ =~ s/[\x0D\x0A]+$//;
            if (substr($_,0,1) eq '>')
            {
                if ($OUT)
                {
                    if ($line_len == 0) { if ($si >= 0) { print $OUT $eol; } }
                    elsif ($buffer ne '') { print $OUT $buffer, $eol; $buffer = ''; }
                }
                $si++;
                if ($si >= $part_start[$part+1])
                {
                    if ($OUT) { close $OUT; }
                    $part++;
                    if ($part > $n_parts_found) { last; }
                    $out_file = $out_dir . $base . $part_num_prefix . ($nopad ? $part : sprintf('%0*d',$num_len,$part)) . $ext;
                    open($OUT,'>',$out_file) or die "Can't create output file \"$out_file\"\n";
                    binmode $OUT;
                }
                print $OUT $_, $eol;
                next;
            }
            if ($line_len)
            {
                $buffer .= $_;
                while (length($buffer) >= $line_len) { print $OUT substr($buffer,0,$line_len,''), $eol; }
            }
            else { print $OUT $_; }
        }
        close $IN;
        if ($OUT)
        {
            if (!$line_len) { if ($si >= 0) { print $OUT $eol; } }
            elsif ($buffer ne '') { print $OUT $buffer, $eol; $buffer = ''; }
            close $OUT;
        }
    }
    else
    {
        ($n_seq,$total_seq_len,$total_size) = get_file_size($infile);
        print ": $n_seq sequences, $total_seq_len bp => dividing into $n_parts part", ($n_parts > 1 ? 's' : ''), " ";
        open(my $IN,'<',$infile) or die "Error: Can't open file \"$infile\"\n";
        binmode $IN;
        $num_len = length($n_parts);
        ($OUT,$name,$data,$written_total,$written_this_part,$part_end,$part) = (undef,undef,'',0,0,int($total_size / $n_parts),1);
        while(<$IN>)
        {
            $_ =~ s/[\x0D\x0A]+$//;
            if (substr($_,0,1) eq '>')
            {
                if (defined $name) { dump_seq(); }
                $name = $_; $data = ''; next;
            }
            $data .= $_;
        }
        if (defined $name) { dump_seq(); }
        close $IN;
        if ($OUT) { close $OUT; }
        print " OK\n";
    }
}

sub dump_seq
{
    my $slen = length($data);
    my $seq_size = seq_size(length($name),$slen);
    my $new_written_total = $written_total + $seq_size;
    if ( !$OUT or
         ($written_this_part and ($new_written_total > $part_end) and ($new_written_total - $part_end > $part_end - $written_total)) )
    {
        if ($OUT) { close $OUT; $part++; $part_end = int($total_size / $n_parts * $part) + 1; }

        my $part_file = $out_dir . $base . $part_num_prefix . ($nopad ? $part : sprintf('%0*d',$num_len,$part)) . $ext;

        open($OUT,'>',$part_file) or die "Error: Can't create file \"$part_file\"\n";
        binmode $OUT;
        $written_this_part = 0;
        print ".";
    }
    print $OUT $name, $eol;
    if ($line_len) { for (my $s=0; $s<$slen; $s+=$line_len) { print $OUT substr($data,$s,$line_len), $eol; } }
    else { print $OUT $data, $eol; }
    $written_this_part += $seq_size;
    $written_total += $seq_size;
}

sub get_file_size_and_part_boundaries
{
    my ($file) = @_;
    open(my $IN,'<',$file) or die "Error: Can't open file \"$file\"\n";
    binmode $IN;
    my ($nseq,$total_seq_length,$total_size,$n_parts_found,$this_part_size,$nlen,$slen,$stop) = (0,0,0,1,0,0,0,0);
    $part_start[1] = 0;
    while (<$IN>)
    {
        $_ =~ s/[\x0D\x0A]+$//;
        my $len = length($_);
        if (substr($_,0,1) eq '>')
        {
            if ($nlen)
            {
                my $seq_size = seq_size($nlen,$slen);
                if ($part_size and $this_part_size and ($this_part_size + $seq_size > $part_size))
                {
                    if ($n_parts and $n_parts_found == $n_parts) { $stop = 1; last; }
                    else { $this_part_size = $seq_size; $n_parts_found++; $part_start[$n_parts_found] = $nseq; }
                }
                else { $this_part_size += $seq_size; }
                $nseq++; $total_seq_length += $slen; $total_size += $seq_size; 
            }
            ($nlen,$slen) = ($len,0); next;
        }
        if ($nlen) { $slen += $len; }
    }
    if ($nlen and !$stop)
    {
        my $seq_size = seq_size($nlen,$slen);
        if ($part_size and $this_part_size and ($this_part_size + $seq_size > $part_size))
        {
            if ($n_parts and $n_parts_found == $n_parts) { $stop = 1; }
            else { $this_part_size = $seq_size; $n_parts_found++; $part_start[$n_parts_found] = $nseq; }
        }
        if (!$stop) { $nseq++; $total_seq_length += $slen; $total_size += $seq_size; }
    }
    close $IN;
    $part_start[$n_parts_found+1] = $nseq;
    return ($nseq,$total_seq_length,$total_size,$n_parts_found);
}

sub get_file_size
{
    my ($file) = @_;
    open(my $IN,'<',$file) or die "Error: Can't open file \"$file\"\n";
    binmode $IN;
    my ($nseq,$total_seq_length,$total_size,$nlen,$slen) = (0,0,0,0,0);
    while (<$IN>)
    {
        $_ =~ s/[\x0D\x0A]+$//;
        my $len = length($_);
        if (substr($_,0,1) eq '>')
        { 
            if ($nlen) { $nseq++; $total_seq_length += $slen; $total_size += seq_size($nlen,$slen); }
            ($nlen,$slen) = ($len,0); next;
        }
        if ($nlen) { $slen += $len; }
    }
    if ($nlen) { $nseq++; $total_seq_length += $slen; $total_size += seq_size($nlen,$slen); }
    close $IN;
    return ($nseq,$total_seq_length,$total_size);
}

sub seq_size
{
    my ($nlen,$slen) = @_;
    return ($measure == 0) ? 1 :
           ($measure == 1) ? $slen :
           $slen + $nlen + $eol_len*(1 + ($line_len ? int(($slen+$line_len-1)/$line_len) : 1));
}
