#!/gsc/bin/perl

BEGIN { 
    $ENV{UR_DBI_NO_COMMIT} = 1;
    $ENV{UR_USE_DUMMY_AUTOGENERATED_IDS} = 1;
}

use strict;
use warnings;

use above "Genome";
use Test::More;

use_ok("Genome::TestObjGenerator::ProcessingProfile::ReferenceAlignment");

my $p = Genome::TestObjGenerator::ProcessingProfile::ReferenceAlignment->setup_object();
ok($p->isa("Genome::ProcessingProfile::ReferenceAlignment"), "Generated a refalign pp");

done_testing;

