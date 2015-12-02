#!/usr/bin/env genome-perl

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
    $ENV{UR_USE_DUMMY_AUTOGENERATED_IDS} = 1;
};

use strict;
use warnings;

use above "Genome";

use Test::More;

use_ok('Genome::Sample') or die;

my $taxon = Genome::Taxon->__define__(name => '__TEST_TAXON__');
ok($taxon, 'define taxon');

my $source = Genome::Individual->__define__(name => '__TEST__IND__', taxon => $taxon);
ok($source, 'define source');

my %sample_params = (
    name             => 'full_name.test',
    common_name      => 'common',
    extraction_label => 'TCGA-1234-232-12',
    extraction_type  => 'genomic dna',
    extraction_desc  => 'This is a test',
    source           => $source,
    tissue_desc      => 'normal',
    age              => 99,
    body_mass_index  => 22.4,
);

my $sample = Genome::Sample->create(%sample_params);

$sample_params{name} = 'H_KA-673778G-S.15614';
$sample_params{extraction_label} = 'ABCD1234';
my $sample_special = Genome::Sample->create(%sample_params);

my @tcga_names = qw(TCGA-12-1234-03A-01D-0742-05 TCGA-12-1234-07A-02D-0332-08);

for my $tcga_name (@tcga_names) {
    my $attr = Genome::SubjectAttribute->create(
        subject_id      => $sample->id,
        attribute_label => 'external_name',
        attribute_value => $tcga_name,
        nomenclature    => 'TCGA',
    );
}


isa_ok($sample, 'Genome::Sample');

is($sample->name_in_vcf, "TCGA-1234-232-12", 'sample TCGA name');
is($sample_special->name_in_vcf, 'TCGA-AB-3012-03A-01D-0741-05', 'special sample TCGA name');
is($sample->resolve_tcga_patient_id, 'TCGA-12-1234', 'resolve_tcga_patient_id');

my @sample_tcga_names = $sample->get_tcga_names;
is_deeply(\@tcga_names, \@sample_tcga_names, 'get_tcga_names');

is($sample->subject_type, 'sample_name', 'subject type is organism sample');
is_deeply($sample->taxon, $source->taxon, 'taxon');
is($sample->age, 99, 'age');
is($sample->body_mass_index, 22.4, 'body_mass_index');

ok(!$sample->is_rna, 'sample is not rna');
$sample->extraction_type('cdna');
ok($sample->is_rna, 'sample is now rna');

done_testing();
