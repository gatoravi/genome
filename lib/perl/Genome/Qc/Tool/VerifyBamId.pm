package Genome::Qc::Tool::VerifyBamId;

use strict;
use warnings;
use Genome;
use List::MoreUtils qw(any);

class Genome::Qc::Tool::VerifyBamId {
    is => 'Genome::Qc::Tool::WithVariationListVcf',
};


sub cmd_line {
    my $self = shift;

    my $cmd = $self->gmt_class->create($self->gmt_params);
    return $cmd->_get_cmd_list;
}

sub get_metrics {
    my $self = shift;

    my $file = $self->qc_metrics_file . '.selfSM';
    my $reader = Genome::Utility::IO::SeparatedValueReader->create(
        input => $file,
        separator => '\t',
        is_regex => 1,
    );

    my %metrics;
    while ( my $line = $reader->next ) {
        if ($line->{'#SEQ_ID'} eq $self->sample_name) {
            for my $column (@{$reader->headers}) {
                next if any { $_ eq $column } $self->_non_metric_columns;
                my $metric_name = join('-', $line->{RG}, $column);
                $metrics{$metric_name} = $line->{$column};
            }
        }
    }
    return %metrics;
}

sub gmt_class {
    return 'Genome::Model::Tools::VerifyBamId';
}

sub qc_metrics_file_accessor {
    return 'out_prefix';
}

sub _non_metric_columns {
    return qw(#SEQ_ID RG);
}

1;
