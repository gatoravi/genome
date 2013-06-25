#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
    $ENV{UR_USE_DUMMY_AUTOGENERATED_IDS} = 1;
};

use File::Path;
use File::Temp;
use File::Compare;
use Test::More;
use above 'Genome';
use Genome::SoftwareResult;

my $archos = `uname -a`;
if ($archos !~ /64/) {
    plan skip_all => "Must run from a 64-bit machine";
}

use_ok('Genome::Model::Tools::DetectVariants2::GatkSomaticIndel');

my $refbuild_id = 101947881;
my $ref_seq_build = Genome::Model::Build::ImportedReferenceSequence->get($refbuild_id);
ok($ref_seq_build, 'human36 reference sequence build') or die;

my $test_data = $ENV{GENOME_TEST_INPUTS} . "/Genome-Model-Tools-DetectVariants2-GatkSomaticIndel/inputs/2013-06-12";
my $expected_data = $ENV{GENOME_TEST_INPUTS} . "/Genome-Model-Tools-DetectVariants2-GatkSomaticIndel/expected_10";
my $tumor =  $test_data.'/true_positive_tumor_validation.bam';
my $normal = $test_data.'/true_positive_normal_validation.bam';

my $tmpbase = File::Temp::tempdir('GatkSomaticIndelXXXXX', CLEANUP => 1, TMPDIR => 1);
my $tmpdir = "$tmpbase/output";

my $gatk_somatic_indel = Genome::Model::Tools::DetectVariants2::GatkSomaticIndel->create(
        aligned_reads_input=>$tumor, 
        control_aligned_reads_input=>$normal,
        reference_build_id => $refbuild_id,
        output_directory => $tmpdir, 
        mb_of_ram => 3500,
        version => 5336,
        aligned_reads_sample => 'TEST_tumor',
        control_aligned_reads_sample => 'TEST_normal',
);

ok($gatk_somatic_indel, 'gatk_somatic_indel command created');
$gatk_somatic_indel->dump_status_messages(1);
my $rv = $gatk_somatic_indel->execute;
is($rv, 1, 'Testing for successful execution.  Expecting 1.  Got: '.$rv);

my @files = qw|     gatk_output_file
                    indels.hq
                    indels.hq.bed
                    indels.hq.v1.bed
                    indels.hq.v2.bed |;

for my $file (@files){
    my $expected_file = "$expected_data/$file";
    my $actual_file = "$tmpdir/$file";
    is(compare($actual_file,$expected_file),0,"Actual file is the same as the expected file: $file");
}

my $expected_vcf = "$expected_data/indels.vcf.gz";
my $output_vcf   = "$tmpdir/indels.vcf.gz";
my $expected = `zcat $expected_vcf | grep -v fileDate`;
my $output   = `zcat $output_vcf | grep -v fileDate`;
my $diff = Genome::Sys->diff_text_vs_text($output, $expected);
ok(!$diff, 'Actual file is the same as the expected file: indels.vcf.gz')
    or diag("diff results:\n" . $diff);
done_testing();
