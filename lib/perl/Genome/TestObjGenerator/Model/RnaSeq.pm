package Genome::TestObjGenerator::Model::RnaSeq;
use Genome::TestObjGenerator::Model;
@ISA = (Genome::TestObjGenerator::Model);

use strict;
use warnings;
use Genome;
use Genome::TestObjGenerator::ProcessingProfile::RnaSeq;

my @required_params = ("subject_name", "subject_type");

sub create_processing_profile_id {
    my $p = Genome::TestObjGenerator::ProcessingProfile::RnaSeq->setup_object();
    return $p->id;
}

sub generate_obj {
    my $self = shift;
    my $m = Genome::Model::RnaSeq->create(@_);
    return $m;
}

sub get_required_params {
    my $class = shift;
    my $super = $class->SUPER::get_required_params;
    my @all = (@$super, @required_params);
    return \@all;
}
sub create_subject_name {
    return Genome::TestObjGenerator::Util::generate_name("test_subject");
}

sub create_subject_type {
    return "sample_name";
}
1;

