#!/usr/bin/env genome-perl

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
    $ENV{UR_USE_DUMMY_AUTOGENERATED_IDS} = 1;
    $ENV{NO_LSF} = 1;
}

use strict;
use warnings;

use above "Genome";
use Test::More;
use Test::Deep;
use File::Basename qw(basename);
use File::Spec;
use Genome::Annotation::TestHelpers qw(get_test_somatic_variation_build get_test_dir);
use Genome::Utility::Test qw(compare_ok);

my $pkg = 'Genome::Annotation::MasterCommand';
use_ok($pkg);

my $version = 1;
my $build_version = 1;
my $test_dir = get_test_dir($pkg, $version);

my $build = get_test_somatic_variation_build(version => $build_version,
    snvs_plan_file => File::Spec->join($test_dir, 'snvs_plan.yaml'));

my $output_dir = Genome::Sys->create_temp_directory;
my $cmd = $pkg->execute(build_id => $build->id, variant_type => 'snvs',
    output_directory => $output_dir);

my $expected_dir = File::Spec->join($test_dir, 'expected');
compare_dir_ok($output_dir, $expected_dir, 'All reports are as expected');

done_testing;

sub compare_dir_ok {
    my ($got_dir, $expected_dir, $message) = @_;

    my @got_files = map {basename($_)} glob(File::Spec->join($got_dir, '*'));
    my @expected_files = map {basename($_)} glob(File::Spec->join($expected_dir, '*'));

    cmp_bag(\@got_files, \@expected_files, 'Got all expected files');

    for my $filename (@got_files) {
        my $got = File::Spec->join($got_dir, $filename);
        my $expected = File::Spec->join($expected_dir, $filename);
        compare_ok($got, $expected, "File ($filename) is as expected");
    }
}