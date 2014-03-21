#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
    $ENV{UR_USE_DUMMY_AUTOGENERATED_IDS} = 1;
    $ENV{NO_LSF} = 1;
};

use above "Genome";
use Test::More;

use_ok('Genome::Model::ClinSeq::Command::MicroarrayCnv') or die;
 
#Define the test where expected results are stored
my $expected_output_dir = $ENV{"GENOME_TEST_INPUTS"} . "Genome-Model-ClinSeq-Command-MicroarrayCnv/2013-12-30/";
ok(-e $expected_output_dir, "Found test dir: $expected_output_dir") or die;

my $temp_dir = Genome::Sys->create_temp_directory();
ok($temp_dir, "created temp directory: $temp_dir") or die;

#Run MicroarrayCNV on the 'apipe-test-clinseq-wer' model
my $clinseq_model = Genome::Model->get(name => 'apipe-test-clinseq-wer');
my $run_microarray_cnv = Genome::Model::ClinSeq::Command::MicroarrayCnv->create(outdir=>$temp_dir, clinseq_model=>$clinseq_model, test=>1, min_cnv_diff=>0.1);
$run_microarray_cnv->queue_status_messages(1);
$run_microarray_cnv->execute();

#Dump the output to a log file
my @output1 = $run_microarray_cnv->status_messages();
my $log_file = $temp_dir . "/RunMicroarrayCnv.log.txt";
my $log = IO::File->new(">$log_file");
$log->print(join("\n", @output1));
$log->close();
ok(-e $log_file, "Wrote message file from microarray-cnv to a log file: $log_file");

#Perform a diff between the stored results and those generated by this test
my @diff = `diff -r -x '*.log.txt' -x '*.pdf' -x '*.stderr' -x '*.stdout' -x '*.jpeg' $expected_output_dir $temp_dir`;
ok(@diff == 0, "Found only expected number of differences between expected results and test results")
or do {
  diag("expected: $expected_output_dir\nactual: $temp_dir\n");
  diag("differences are:");
  diag(@diff);
  my $diff_line_count = scalar(@diff);
  print "\n\nFound $diff_line_count differing lines\n\n";
  Genome::Sys->shellcmd(cmd => "rm -fr /tmp/last-run-microarray-cnview/");
  Genome::Sys->shellcmd(cmd => "mv $temp_dir /tmp/last-run-microarray-cnview");
};

done_testing();
