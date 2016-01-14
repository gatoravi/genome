#!/usr/bin/env genome-perl
use above "Genome";
use Test::More tests => 10;

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
    $ENV{UR_USE_DUMMY_AUTOGENERATED_IDS} = 1;
}

class Genome::Foo { roles => 'Genome::Role::Notable' };
ok(Genome::Foo->does('Genome::Role::Notable'), 'made a "Genome::Role::Notable" test class Genome::Foo');

my $o1;
subtest 'created Genome::Foo(100)' => sub {
    plan tests => 2;
    $o1 = Genome::Foo->create(100);
    ok($o1, "created a test notable object");
    is(count($o1->notes), 0, 'no notes at start');
};

subtest 'added note to Genome::Foo(100)' => sub {
    plan tests => 3;
    my $n1 = $o1->add_note(
        header_text => "head1",
        body_text => "body1",
    );
    ok($n1, "added note");
    is($n1->header_text, 'head1', 'header is okay');
    is($n1->body_text, 'body1', 'body is okay');
};

my $n2;
subtest 'added second note to Genome::Foo(100)' => sub {
    plan tests => 1;
    $n2 = $o1->add_note(
        header_text => 'head2',
        body_text => 'body2',
    );
    ok($n2, "added a 2nd note");
};

my $o2;
subtest 'created Genome::Foo(200)' => sub {
    plan tests => 2;
    $o2 = Genome::Foo->create(200);
    ok($o2, "made a 2nd notable object");
    is(count($o2->notes), 0, 'no notes at start');
};

subtest 'added note to Genome::Foo(200)' => sub {
    plan tests => 1;
    my $n3 = $o2->add_note(
        header_text => 'head3',
        body_text => 'body3',
    );
    ok($n3, "added a note to the 2nd object");
};

is(count($o1->notes), 2, 'got expected note count for object 1');
is(count($o2->notes), 1, 'got expected note count for object 2');

subtest 'query note' => sub {
    plan tests => 2;
    my $a1 = $o1->notes(header_text => 'head2');
    ok($a1,"got expected note");
    is($a1->body_text, 'body2', 'got correct header');
};

subtest 'remove note' => sub {
    plan tests => 2;
    ok($o1->remove_note($n2), "removed the 2nd note from object 1");
    is(count($o1->notes), 1, "got expected note count for object 1");
};

UR::Context->_sync_databases() or die;

sub count {
    my @list = @_;
    return scalar(@list);
}