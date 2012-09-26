package Genome::File::Vep::Reader;

use Genome::File::Vep::Header;
use Genome::File::Vep::Entry;
use Data::Dumper;
use Genome;
use Carp qw/confess/;
use strict;
use warnings;

sub new {
    my ($class, $path) = @_;
    return $class->fhopen(Genome::Sys->open_file_for_reading($path), $path);
}

sub fhopen {
    my ($class, $fh, $name) = @_;
    $name |= "unknown file path";
    my $header_txt = $fh->getline;

    my $self = {
        name => $name,
        filehandle => $fh,
        header => new Genome::File::Vep::Header($header_txt),
        line_num => 1,
    };

    return bless $self, $class;
}

sub next {
    my $self = shift;
    while (my $line = $self->{filehandle}->getline) {
        ++$self->{line_num};
        chomp $line;
        # There are blank lines in vep files sometimes ._.
        next unless $line;
        return Genome::File::Vep::Entry->new($line);
    }
}

1;
