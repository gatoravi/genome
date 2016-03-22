#!/usr/bin/env genome-perl

#Written by Malachi Griffith

use strict;
use warnings;
use File::Basename;
use Cwd 'abs_path';

BEGIN {
    $ENV{UR_DBI_NO_COMMIT}               = 1;
    $ENV{UR_USE_DUMMY_AUTOGENERATED_IDS} = 1;
}

use above "Genome";
use Test::More tests => 7;  #One per 'ok', 'is', etc. statement below
use Genome::Model::ClinSeq::Command::RunCnView;
use Data::Dumper;

use_ok('Genome::Model::ClinSeq::Command::RunCnView') or die;

#Define the test where expected results are stored
my $expected_output_dir = Genome::Config::get('test_inputs') . "Genome-Model-ClinSeq-Command-RunCnView/2013-02-01/";
ok(-e $expected_output_dir, "Found test dir: $expected_output_dir") or die;

#Use an existing cnvhmm file for testing purposes, kept with the expected out
my $input_file = Genome::Config::get('test_inputs') . "Genome-Model-ClinSeq-Command-RunCnView/cnaseq.cnvhmm.input";
ok(-e $input_file, "Found input file: $input_file") or die;

#Create a temp dir for results
my $temp_dir = Genome::Sys->create_temp_directory();
ok($temp_dir, "created temp directory: $temp_dir");

#Get a somatic-variation build
my $build_id = 129399487;
my $build    = Genome::Model::Build->get($build_id);

#Create run-cn-view command and execute
#genome model clin-seq run-cn-view --outdir=/tmp/cnv/ --cnv-hmm-file=/gscmnt/gc1401/info/model_data/2889933976/build132760359/AML109/clonality/cnaseq.cnvhmm --test 129399487

my $cancer_annotation_db = Genome::Db->get("tgi/cancer-annotation/human/build37-20130401.1");
my $run_cn_view_cmd      = Genome::Model::ClinSeq::Command::RunCnView->create(
    outdir               => $temp_dir,
    cnv_hmm_file         => $input_file,
    test                 => 1,
    build                => $build,
    cancer_annotation_db => $cancer_annotation_db
);
$run_cn_view_cmd->queue_status_messages(1);
my $r1 = $run_cn_view_cmd->execute();
is($r1, 1, 'Testing for successful execution.  Expecting 1.  Got: ' . $r1);

#Dump the output to a log file
my @output1  = $run_cn_view_cmd->status_messages();
my $log_file = $temp_dir . "/RunCnView.log.txt";
my $log      = IO::File->new(">$log_file");
$log->print(join("\n", @output1));
ok(-e $log_file, "Wrote message file from run-cn-view to a log file: $log_file");

#The first time we run this we will need to save our initial result to diff against
#Genome::Sys->shellcmd(cmd => "cp -r -L $temp_dir/* $expected_output_dir");

#Perform a diff between the stored results and those generated by this test
my @diff = `diff -r -x '*.log.txt' -x '*.pdf' -x '*.stderr' -x '*.stdout' -x '*.jpeg' $expected_output_dir $temp_dir`;
ok(@diff == 0, "Found only expected number of differences between expected results and test results")
    or do {
    diag("expected: $expected_output_dir\nactual: $temp_dir\n");
    diag("differences are:");
    diag(@diff);
    my $diff_line_count = scalar(@diff);
    print "\n\nFound $diff_line_count differing lines\n\n";
    Genome::Sys->shellcmd(cmd => "rm -fr /tmp/last-run-cn-view/");
    Genome::Sys->shellcmd(cmd => "mv $temp_dir /tmp/last-run-cn-view");
    };

