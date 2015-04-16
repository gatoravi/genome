package Genome::Model::ClinSeq::Command::IdentifyLoh;

use strict;
use warnings;
use File::Spec;
use Genome;

class Genome::Model::ClinSeq::Command::IdentifyLoh {
  is => ['Command::V2',
        'Genome::Model::ClinSeq::Util',
  ],
  has_input => [
    clinseq_build => {
        is => 'Genome::Model::Build::ClinSeq',
        doc => 'ClinSeq build to identify LOH regions in.
            [Either this or a somatic variation build is required.]',
        is_optional => 1,
    },
    outdir => {
        is => 'FilesystemPath',
        doc => 'Directory where output files will be written.',
    },
    somvar_build => {
        is => 'Genome::Model::Build::SomaticVariation',
        doc => 'SomVar build to identify LOH regions in.
        [Either this or a clinseq build is required.]',
        is_optional => 1,
    },
    minprobes=> {
        is => 'Number',
        is_optional => 1,
        doc => 'Minimum number of probes to call a loh segment.',
        default => 10,
    },
    segment_pc_cutoff => {
        is => 'Number',
        is_optional => 1,
        doc => 'Minimum segment mean.',
        default => 0.95,
    },
    germline => {
        is => 'FilesystemPath',
        is_optional => 1,
        doc => 'optional custom germline file in the somatic variation ' .
            'snvs.Germline.hc format',
    },
    loh => {
        is => 'FilesystemPath',
        is_optional => 1,
        doc => 'optional custom loh file in the somatic variation ' .
            'snvs.LOH.hc format',
    },
    prefix => {
        is => 'Text',
        is_optional => 1,
        doc => 'Optional output-file prefix. These files will be in the outdir.',
    },
    test => {
        is => 'Boolean',
        doc => 'set for test-cases',
        is_optional => 1,
        default => 0,
    },
    ],
    doc => 'Identify regions of LOH using clinseq tumor/normal pairs.',
};

sub help_synopsis {
    return <<EOS
genome model clin-seq identify-loh\\
    --outdir=/tmp/test/ \\
    --clinseq-build='a4abcd1313eb4376b59e68a9dd9d5ad2'
genome model clin-seq identify-loh\\
    --outdir=/tmp/test/ \\
    --somvar-build='9dc0385a1b634c9bb85eb2017b3a0c73'
EOS
}

sub help_detail {
    return <<EOS
Use the results from Varscan in the somatic-variation builds and identify
regions of LOH. Specifically, look at heterozygous single nucleotide variants
in the normal sample and compare with stretches of homozygosity in the tumor.

If a clinseq build is supplied as input, the tool attempts to use the results
from the WGS-somatic-variation build input, if these are unavailable the tool
looks for WEx-somatic-variation build.

Alternatively, instead of a clin-seq build you can supply a
somatic-variation build as an input.

EOS
}

sub __errors__ {
    my $self = shift;
    my @errors = $self->SUPER::__errors__(@_);
    unless (-e $self->outdir && -d $self->outdir) {
        push @errors, UR::Object::Tag->create(
            type => 'error',
            properties => ['outdir'],
            desc => "Outdir: " . $self->outdir .
                " not found or not a directory",
        );
    }
    return @errors;
}

#Get a somatic-variation build from input or from clinseq-build
sub resolve_somvar {
    my $self = shift;
    my $somvar_build = $self->somvar_build;
    unless($somvar_build) {
        my $clinseq_build = $self->clinseq_build;
        unless($clinseq_build) {
            die $self->error_message("Please pass a somatic-variation
                or clinseq build as input.");
        }
        $somvar_build = $self->get_best_somvar_build($clinseq_build);
    }
    return $somvar_build;
}

#Get Varscan SNVs from somatic-variation
sub get_germline_loh {
    my $self = shift;
    my $somvar_build = shift;
    my $variants = shift;
    my $dd =
        glob(File::Spec->join($somvar_build->data_directory,
                              "variants/snv/varscan-somatic-*/varscan-high-confidence-v1-*"));
    unless(-e $dd) {
        die $self->error_message("Unable to find varscan SNV directory.");
    }
    $self->status_message("varscan snvs directory found: $dd");
    $variants->{germline} = File::Spec->join($dd,
                                           "snvs.Germline.hc");
    $variants->{loh} = File::Spec->join($dd,
                                      "snvs.LOH.hc");
    if($self->germline) {
        $variants->{germline} = $self->germline;
    }
    if($self->loh) {
        $variants->{loh} = $self->loh;
    }
    unless(-e $variants->{germline} or -e $variants->{germline}) {
        die $self->error_message("Unable to find germline or LOH in somatic-variation");
    }
}

#Combine germline, LOH calls and sort them
sub combine_sort_snvs {
    my $self = shift;
    my $snv_prefix = shift;
    my $variants = shift;
    my $germline = $variants->{germline};
    my $loh = $variants->{loh};
    my $snv_combined = $snv_prefix . ".combined";
    my $snv_combined_sorted = $snv_prefix . ".combined.sorted";
    Genome::Sys->cat(
        input_files => [$germline, $loh],
        output_file => $snv_combined);
    my $sort = Genome::Model::Tools::Capture::SortByChrPos->create(
        input_file => $snv_combined,
        output_file => $snv_combined_sorted);
    $sort->execute;
    unlink $snv_combined;
    return $snv_combined_sorted;
}

#Segment the combined calls
sub segment_loh {
    my $self = shift;
    my $snvs = shift;
    my $loh_basename = shift;
    my $segmenter = Genome::Model::Tools::Varscan::LohSegments->
    create(
        variant_file => $snvs,
        output_basename => $loh_basename,
    );
    $segmenter->execute();
}

#Filter out LOH segments without sufficient support
sub filter_loh {
    my $self = shift;
    my $loh_basename = shift;
    my $segment_pc_cutoff = $self->segment_pc_cutoff;
    my $minprobes = $self->minprobes;
    my $loh_segments =
        $loh_basename . ".segments.cbs";
    my $loh_segments_filtered = $loh_segments . ".filtered";
    unless(-e $loh_segments) {
        die $self->error_message("Unable to find loh file $loh_segments");
    }
    my $awk_filter =
        "awk '\$5 > $segment_pc_cutoff && \$4 >= $minprobes' $loh_segments " .
            "> $loh_segments_filtered";
    Genome::Sys->shellcmd(cmd => $awk_filter);
}

#Get, Filter, Combine, Segment, Filter, Cleanup
sub execute {
    my $self = shift;
    my %variants;
    my $snv_prefix = $self->prefix || File::Spec->join(
        $self->outdir,
        "varscan.snps");
    my $somvar_build = $self->resolve_somvar;
    $self->get_germline_loh($somvar_build, \%variants);
    $self->combine_sort_snvs($snv_prefix, \%variants);
    my $loh_basename = File::Spec->join(
        $self->outdir,
        "varscan.output.loh");
    my $combined_sorted = $snv_prefix . ".combined.sorted";
    my $loh_segments = $self->segment_loh($combined_sorted, $loh_basename);
    $self->filter_loh($loh_basename);
    return 1;
}

1;

