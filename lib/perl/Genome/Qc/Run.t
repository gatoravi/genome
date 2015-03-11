#!/usr/bin/env genome-perl

BEGIN { 
    $ENV{UR_DBI_NO_COMMIT} = 1;
    $ENV{UR_USE_DUMMY_AUTOGENERATED_IDS} = 1;
}

use strict;
use warnings;

use above "Genome";
use Test::More;
use Genome::Test::Factory::InstrumentData::MergedAlignmentResult;
use Genome::Test::Factory::SoftwareResult::User;
use Sub::Install qw (reinstall_sub);

my $pkg = "Genome::Qc::Run";
use_ok($pkg);

{
    package TestTool1;

    use Genome;

    class TestTool1 {
        is => [$pkg],
        has => {param1 => {}},
    };

    sub cmd_line {
        my $self = shift;
        return ("echo", $self->param1);
    }

    sub supports_streaming {
        return 0;
    }

    sub get_metrics {
        return {
            metric1 => 1,
            metric2 => 2,
        };
    }
}

{
    package TestTool2;

    use Genome;

    class TestTool2 {
        is => [$pkg],
        has => {param2 => {}},
    };

    sub cmd_line {
        my $self = shift;
        return ("echo", $self->param2);
    }

    sub supports_streaming {
        return 0;
    }

    sub get_metrics {
        return {
            metricA => 1,
            metricB => 2,
        };
    }
}

use Genome::Qc::Config;
reinstall_sub({
    into => 'Genome::Qc::Config',
    as => 'get_commands_for_alignment_result',
    code => sub {
        return {test1 => {class => "TestTool1", params => {param1 => 1}},
            test2 => {class => "TestTool2", params => {param2 => 2}}};
    },
});

my $alignment_result = Genome::Test::Factory::InstrumentData::MergedAlignmentResult->setup_object();
my $command = $pkg->create(
    alignment_result => $alignment_result,
    %{Genome::Test::Factory::SoftwareResult::User->setup_user_hash},
);

ok($command->execute, "Command executes ok");
ok($command->output_result, "Result is defined after command executes");

my $expected_metrics = {
    metric1 => 1,
    metric2 => 2,
    metricA => 1,
    metricB => 2,
};
is_deeply({$command->output_result->get_metrics}, $expected_metrics, "Metrics correctly set on result");

done_testing;

