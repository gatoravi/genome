#!/gsc/bin/perl

BEGIN { 
    $ENV{UR_DBI_NO_COMMIT} = 1;
    $ENV{UR_USE_DUMMY_AUTOGENERATED_IDS} = 1;
}

use strict;
use warnings;

use above "Genome";
use Test::More;
use Genome::Utility::Test qw(compare_ok);
use Genome::Test::Factory::Model::ReferenceAlignment;
use Genome::Test::Factory::Model::SomaticVariation;
use Genome::Test::Factory::Build;

my $class = "Genome::Model::Tools::Tcga::CreateSubmissionArchive";
use_ok($class);

my $base_dir = Genome::Utility::Test->data_dir_ok($class, "v1");

my @headers = $class->get_sdrf_headers;

my @expected_headers = (
        'Material Extract Name',
        'Material Comment [TCGA Barcode]',
        'Material Comment [is tumor]',
        'Material Material Type',
        'Material Annotation REF',
        'Material Comment [TCGA Genome Reference]',
        'Library Protocol REF',
        'Library Parameter Value [Vendor]',
        'Library Parameter Value [Catalog Name]',
        'Library Parameter Value [Catalog Number]',
        'Library Parameter Value [Annotation URL]',
        'Library Parameter Value [Product URL]',
        'Library Parameter Value [Target File URL]',
        'Library Parameter Value [Target File Format]',
        'Library Parameter Value [Target File Format Version]',
        'Library Parameter Value [Probe File URL]',
        'Library Parameter Value [Probe File Format]',
        'Library Parameter Value [Probe File Format Version]',
        'Library Parameter Value [Target Reference Accession]',
        'Sequencing Protocol REF',
        'Mapping Protocol REF',
        'Mapping Comment [Derived Data File REF]',
        'Mapping Comment [TCGA CGHub ID]',
        'Mapping Comment [TCGA CGHub metadata URL]',
        'Mapping Comment [TCGA Include for Analysis]',
        'Mapping2 Derived Data File',
        'Mapping2 Comment [TCGA Include for Analysis]',
        'Mapping2 Comment [TCGA Data Type]',
        'Mapping2 Comment [TCGA Data Level]',
        'Mapping2 Comment [TCGA Archive Name]',
        'Mapping2 Parameter Value [Protocol Min Base Quality]',
        'Mapping2 Parameter Value [Protocol Min Map Quality]',
        'Mapping2 Parameter Value [Protocol Min Tumor Coverage]',
        'Mapping2 Parameter Value [Protocol Min Normal Coverage]',
        'Variants Protocol REF',
        'Variants Derived Data File',
        'Variants Comment [TCGA Spec Version]',
        'Variants Comment [TCGA Include for Analysis]',
        'Variants Comment [TCGA Data Type]',
        'Variants Comment [TCGA Data Level]',
        'Variants Comment [TCGA Archive Name]',
        'Maf Protocol REF',
        'Maf Derived Data File',
        'Maf Comment [TCGA Spec Version]',
        'Maf Comment [TCGA Include for Analysis]',
        'Maf Comment [TCGA Data Type]',
        'Maf Comment [TCGA Data Level]',
        'Maf Comment [TCGA Archive Name]',
        'Validation Protocol REF',
        'Validation Derived Data File',
        'Validation Comment [TCGA Spec Version]',
        'Validation Comment [TCGA Include for Analysis]',
        'Validation Comment [TCGA Data Type]',
        'Validation Comment [TCGA Data Level]',
        'Validation Comment [TCGA Archive Name]'
        );

is_deeply(\@headers, \@expected_headers, "Headers generated correctly");

my %empty_row;
my $null_row = $class->fill_in_nulls(\%empty_row);
my @empty_keys = sort keys %$null_row;
my @sorted_headers = sort @expected_headers;
is_deeply(\@empty_keys, \@sorted_headers, "Empty row got filled in");

my $test_output = Genome::Sys->create_temp_file_path;
ok($class->print_sdrf($test_output, $null_row), "print_sdrf ran ok with a row of nulls");
compare_ok($test_output, $base_dir."/expected_null.sdrf", "null sdrf printed correctly");

my $test_somatic_build = Genome::Test::Factory::Model::SomaticVariation->setup_somatic_variation_build();
$test_somatic_build->normal_build->subject->common_name("normal");
$test_somatic_build->normal_build->subject->extraction_label("TCGA-1");
$test_somatic_build->normal_build->subject->source->upn("TCGA-UPN-A");
$test_somatic_build->normal_build->model->target_region_set_name("11111001 capture chip set");
$test_somatic_build->normal_build->data_directory($base_dir."/refalign_dir");

$test_somatic_build->tumor_build->subject->common_name("tumor");
$test_somatic_build->tumor_build->subject->extraction_label("TCGA-2");
$test_somatic_build->tumor_build->model->target_region_set_name("SeqCap EZ Human Exome v2.0");
$test_somatic_build->tumor_build->data_directory($base_dir."/refalign_dir2");

$test_somatic_build->data_directory("$base_dir/somvar_dir");

my $cghub_ids = Genome::Sys->create_temp_file_path;
`echo "CGHub_ID\tTCGA_Name\tBAM_path\ncghub1\tTCGA-1\t/dev/null\ncghub2\tTCGA-2\t/dev/null" > $cghub_ids`;

is_deeply($class->load_cghub_info($cghub_ids, "TCGA_Name"), {"TCGA-1" => "cghub1", "TCGA-2" => "cghub2"}, "CGHub info loaded correctly");
is($class->resolve_cghub_id($test_somatic_build->normal_build, $cghub_ids), "cghub1", "CGHub called correctly");

my $sample_1 = {
   ID => {content => "TCGA_1"},
   SampleUUID => {content => "3958t6"},
   SampleTCGABarcode => {content => "TCGA_1"},
};

my $row1 = $class->create_vcf_row($test_somatic_build->normal_build, "test_archive", undef, $cghub_ids, "snvs.vcf", $sample_1);
my $row2 = $class->create_vcf_row($test_somatic_build->normal_build, "test_archive", undef, $cghub_ids, "indels.vcf", $sample_1);
my $row3 = $class->create_maf_row($test_somatic_build->normal_build, "test_archive", "/test/maf/path", undef, $cghub_ids, $sample_1);
my $row4 = $class->create_maf_row($test_somatic_build->tumor_build, "test_archive", "/test/maf/path", undef, $cghub_ids, $sample_1);

my $output_sdrf = Genome::Sys->create_temp_file_path;
ok($class->print_sdrf($output_sdrf, ($row1, $row2, $row3, $row4)), "sdrf printed");

my %protocol_db = (
    "library preparation" => [
        {name => "libraryprep1", description => "First library prep protocol"}
    ],
    "nucleic acid sequencing" => [
        {name => "sequencing1", description => "First sequencing protocol"},
    ],
    "sequence alignment" => [
        {name => "alignment1", description => "First mapping protocol"},
    ],
    "variant calling" => [
        {name => "variants1", description => "First variant detection protocol"},
    ],
    "mutation filtering and annotation" => [
        {name => "maf1", description => "First filtering protocol"},
    ],
);
my $output_idf = Genome::Sys->create_temp_file_path;
ok($class->print_idf($output_idf, \%protocol_db), "Print idf called successfully");
compare_ok($output_idf, $base_dir."/expected.idf", "idf printed as expected");

my $archive_output_dir = Genome::Sys->create_temp_directory;
my $cmd = Genome::Model::Tools::Tcga::CreateSubmissionArchive->create(
    models => [$test_somatic_build->model],
    output_dir => $archive_output_dir,
    archive_name => "test_archive",
    archive_version => "1.0.0",
    cghub_id_file => $cghub_ids,
    create_archive => 1,
);
ok($cmd, "Command created");
ok($cmd->execute, "Command executed");

is($class->resolve_maf_protocol, "genome.wustl.edu:maf_creation:data_consolidation:01", "Maf protocol resolved correctly");
is($class->resolve_mapping_protocol($test_somatic_build->normal_build), "genome.wustl.edu:alignment:".$test_somatic_build->normal_build->processing_profile->id.":01", "Mapping protocol resolved correctly");
is($class->resolve_library_protocol, "genome.wustl.edu:DNA_extraction:Illumina_DNASeq:01", "Library protocol resolved correctly");
is($class->resolve_variants_protocol($test_somatic_build->normal_build), "genome.wustl.edu:variant_calling:".$test_somatic_build->normal_build->processing_profile->id.":01", "Variants protocol defined correctly");

is_deeply([$class->resolve_capture_reagent($test_somatic_build->normal_build)], ["Nimblegen", "Nimblegen EZ Exome v3.0", "06465692001"], "Capture reagent resolved correctly");

done_testing;
