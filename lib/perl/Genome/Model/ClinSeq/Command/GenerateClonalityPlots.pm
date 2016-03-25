package Genome::Model::ClinSeq::Command::GenerateClonalityPlots;
#Written by Malachi Griffith and Nate Dees
#Updated into tool form by Scott Smith

use strict;
use warnings;
use Genome;

class Genome::Model::ClinSeq::Command::GenerateClonalityPlots {
    is        => ['Command::V2', 'Genome::Model::ClinSeq::Util'],
    has_input => [
        somatic_var_build => {
            is          => 'Genome::Model::Build',
            id_by       => 'somatic_var_build_id',
            doc         => 'Build ID for a somatic variation model',
            is_optional => 1
        },

        misc_annotation_db => {
            is  => 'Genome::Db::Tgi::MiscAnnotation',
            doc => 'misc annotation for the reference sequence, supplying centromere ideogram and gaps'
        },

        common_name => {
            is  => 'Text',
            doc => 'Human readable name for the patient / sample comparison'
        },

        output_dir => {
            is  => 'Text',
            doc => 'Directory to place temp files and results'
        },
    ],
    has_param => [
        read_counts => {
            is          => 'Text',
            is_optional => 1,
            doc         => 'Instead of generating read counts, use this pre-prepared readcount file'
        },

        verbose => {
            is            => 'Boolean',
            is_optional   => 1,
            default_value => 0,
            doc           => 'To display more output, set this flag'
        },

        limit => {
            is          => 'Number',
            is_optional => 1,
            doc         => 'Limit the number of SNVs to the first N (mostly for testing).'
        },

        chromosome => {
            is          => 'Text',
            is_optional => 1,
            doc         => 'Limit analysis to variants on a specified chromosome'
        },

        min_tumor_cov => {
            is            => 'Number',
            is_optional   => 1,
            default_value => 20,
            doc => 'Specify a minimum tumor coverage level for variants (you will still get another plot with all data)'
        },

        max_tumor_cov => {
            is            => 'Number',
            is_optional   => 1,
            default_value => 1000000,
            doc => 'Specify a maximum tumor coverage level for variants (you will still get another plot with all data)'
        },

        max_normal_vaf => {
            is            => 'Number',
            is_optional   => 1,
            default_value => 3,
            doc           => 'Specify a maximum allowed normal VAF (you will still get another plot with all data)'
        },

        bam_readcount_version => {
            is  => 'Text',
            doc => 'The version of bam-readcounts to use',
        },
    ],
    has_output => [
        cnv_hmm_file => {
            is          => 'FilesystemPath',
            is_optional => 1,
        },
        cnv_hq_file => {
            is          => 'FilesystemPath',
            is_optional => 1,
        },
    ],
    doc => "creates clonality plots"
};

sub help_synopsis {
    return <<INFO
  genome model clin-seq generate-clonality-plots --somatic-var-build=129396826  --output-dir=/tmp/generate_clonality_plots/ --common-name='AML54'  --verbose
INFO
}

sub execute {
    my $self = shift;

    #This script running a series of commands obtained from Nate Dees that results in the creation of a clonality plot (.pdf)
    my $somatic_var_build = $self->somatic_var_build;
    my $output_dir        = $self->output_dir;
    my $common_name       = $self->common_name;
    my $verbose           = $self->verbose;
    my $limit             = $self->limit;
    my $chromosome        = $self->chromosome;

    if (not defined $limit) {
        $self->debug_message("limit is not defined");
    }
    else {
        $self->debug_message("limit is $limit");
    }

    if (not defined $chromosome) {
        $self->debug_message("chromosome filter is not defined");
    }
    else {
        $self->debug_message("chromosome filter is $chromosome");
    }

    $self->debug_message("Creating clonality plot for $common_name");

    #TODO: Replace all of this with a new process that gets variants from a unified clin-seq BAM read counts result
    #TODO: This step should just run the clonality tool in different ways on different input files.  All of this hacky file manipulation should be removed

    #Get somatic variation effects dir, tumor bam and normal bam from a somatic variation model ID
    my %data_paths;
    my $is_copycat = $self->_is_copycat_somvar($somatic_var_build);
    $self->get_data_paths(\%data_paths, $somatic_var_build, $is_copycat, $output_dir);
    my $somatic_effects_dir = $data_paths{effects_dir};

    #Make sure the specified parameters are correct
    $somatic_effects_dir = $self->checkDir(
        '-dir'   => $somatic_effects_dir,
        '-clear' => "no"
    );
    $output_dir = $self->checkDir(
        '-dir'   => $output_dir,
        '-clear' => "no"
    );

    #Step 1 - gather the tier 1-3 snv files from the build:
    my $tier1_snv_file = $somatic_effects_dir . "snvs.hq.novel.tier1.v2.bed";
    my $tier2_snv_file = $somatic_effects_dir . "snvs.hq.novel.tier2.v2.bed";
    my $tier3_snv_file = $somatic_effects_dir . "snvs.hq.novel.tier3.v2.bed";
    my $cp_cmd         = "cp $tier1_snv_file $tier2_snv_file $tier3_snv_file $output_dir";
    if ($verbose) {$self->debug_message("$cp_cmd")}
    Genome::Sys->shellcmd(cmd => $cp_cmd);

    #Step 2 - put them together in one file:
    my $snv_file = $output_dir . "allsnvs.hq.novel.tier123.v2.bed";
    my $cat_cmd  = "cat $output_dir" . "snvs* > $snv_file";

    if ($verbose) {$self->debug_message("$cat_cmd");}
    Genome::Sys->shellcmd(cmd => $cat_cmd);

    #Apply the chromosome filter if specified
    if (defined $chromosome or defined $limit) {
        $self->warning_message("limiting SNVs to only those on chromosome: $chromosome") if defined $chromosome;
        $self->warning_message("limiting SNVs to a max of $limit")                       if defined $limit;

        my $file        = $snv_file . ".tmp";
        my $snv_fh      = Genome::Sys->open_file_for_reading($snv_file);
        my $filtered_fh = Genome::Sys->open_file_for_writing($file);
        my $s           = 0;
        while (<$snv_fh>) {
            if (!defined($chromosome) or $_ =~ /^$chromosome\s+/) {
                print $filtered_fh $_;
                $s++;
            }

            last if defined($limit) and $s >= $limit;
        }
        $snv_fh->close;
        $filtered_fh->close;
        $self->debug_message("Filtered down to $s variants on chromosome: $chromosome");
        unlink($snv_file);
        Genome::Sys->move_file($file, $snv_file);
    }

    #Step 3 - take it out of bed format to be fed into bam-readcounts:
    my $adapted_file = "$output_dir" . "allsnvs.hq.novel.tier123.v2.bed.adapted";
    my $awk_cmd =
        "awk \'{OFS=\"\\t\";FS=\"\\t\";}{print \$1,\$3,\$3,\$4}\' $snv_file | sed \'s/\\//\\t/g\' > $adapted_file";
    if ($verbose) {$self->debug_message("$awk_cmd");}
    Genome::Sys->shellcmd(cmd => $awk_cmd);

    #Define the BAM files.
    #The 'old' method supplied both Tumor & Normal coverage and both would be used to assess a minimum coverage cutoff for plotting
    #The 'new' method uses only the Tumor coverage
    my $tumor_bam  = $data_paths{tumor_bam};
    my $normal_bam = $data_paths{normal_bam};

    #Step 4 - run bam readcounts and assess the particular reads for the reference and variant and print out details about the numbers of reads and the percentages for multiple bam files:
    my $readcounts_outfile;
    if ($self->read_counts) {
        $readcounts_outfile = $self->read_counts;
    }
    else {
        $readcounts_outfile = "$adapted_file" . ".readcounts";
        my $read_counts_cmd = Genome::Model::ClinSeq::Command::GenerateClonalityPlots::Readcounts->create(
            sites_file            => $adapted_file,
            bam_files             => ["Tumor:$tumor_bam", "Normal:$normal_bam"],
            reference_build       => $somatic_var_build->reference_sequence_build,
            output_file           => $readcounts_outfile,
            bam_readcount_version => $self->bam_readcount_version,
        );
        $read_counts_cmd->execute
            or die $self->error_message('Failed to generate readcounts');
    }

    #Step 5 - create a varscan-format file from these outputs:
    my $readcounts_varscan_file = "$readcounts_outfile" . ".varscan";
    $self->create_pseudo_varscan_file($adapted_file, $readcounts_outfile, $readcounts_varscan_file);

    #TODO: Replace steps 4-5 above by using the following script:
    #TODO: Once you know this is working, use $readcounts_clonality_outfile instead of $readcounts_varscan_file in the clonality commands below.  Then comment out steps 2-5 above
    #gmt validation prepare-wgs-for-clonality-plot --help
    #USAGE
    # gmt validation prepare-wgs-for-clonality-plot --output-file=? --snv-file=? [--bam-file=?]
    #    [--genome-build=?] [--min-mapping-quality=?] [--output-readcounts-file=?] [--readcounts-file=?]
    #Use the optional --bam-file input so that readcounts are generated for you.

    #my $readcounts_clonality_outfile = "$output_dir"."readcounts.clonality";
    #my $readcounts_formatted_outfile = "$output_dir"."readcounts.formatted";
    #my $prepare_cmd = Genome::Model::Tools::Validation::PrepareWgsForClonalityPlot->create(output_file=>$readcounts_clonality_outfile, snv_file=>$adapted_file, bam_file=>$tumor_bam, genome_build=>$data_paths{reference_fasta}, output_readcounts_file=>$readcounts_formatted_outfile);
    #$prepare_cmd->execute();

    #Step 6 - Take the cnvs.hq file from the somatic-variation build, and run the cna-seg tool to create known regions of copy-number
    #Specify config file paths for hg19/build37
    #gmt copy-number cna-seg --copy-number-file=/gscmnt/ams1184/info/model_data/2875816457/build111674790/variants/cnvs.hq  --min-markers=4  --detect-somatic  --centromere-file=/gscmnt/sata186/info/medseq/kchen/work/SolexaCNV/scripts/centromere.hg19.csv  --gap-file=/gscmnt/sata186/info/medseq/kchen/work/SolexaCNV/scripts/hg19gaps.csv  --output-file=hg1.cnvhmm

    #Make a copy of the cnvs.hq file
    my $cnvs_output_path = "${output_dir}cnvs.hq";
    Genome::Sys->copy_file($data_paths{cnvs_hq}, $cnvs_output_path);
    chmod 0660, $cnvs_output_path;

    my $misc_annotation_db = $self->misc_annotation_db;
    my $centromere_file    = $misc_annotation_db->data_set_path("centromere.csv");
    my $gap_file           = $misc_annotation_db->data_set_path("gaps.csv");

    #if ($data_paths{display_name} =~ /NCBI\-human\-build36/){
    #    $centromere_file = $clinseq_annotations_dir . "hg18/ideogram/centromere.hg18.csv";
    #    $gap_file = $clinseq_annotations_dir . "hg18/ideogram/hg18gaps.csv";
    #}elsif($data_paths{display_name} =~ /GRCh37\-lite\-build37/){
    #    $centromere_file = $clinseq_annotations_dir . "hg19/ideogram/centromere.hg19.csv";
    #    $gap_file = $clinseq_annotations_dir . "hg19/ideogram/hg19gaps.csv";
    #}else{
    #    $self->error_message("Unrecognized build - unable to identify centromere and gapfiles, you will need to generate these and place in the appropriate location");
    #    return;
    #}

    my $cnvhmm_file = "$output_dir" . "cnaseq.cnvhmm";

    my $cnaseg_cmd;
    if (defined($chromosome)) {
        $cnaseg_cmd = Genome::Model::Tools::CopyNumber::CnaSeg->create(
            copy_number_file => $data_paths{cnvs_hq},
            min_markers      => 4,
            detect_somatic   => 1,
            centromere_file  => $centromere_file,
            gap_file         => $gap_file,
            output_file      => $cnvhmm_file,
            chromosome       => $chromosome
        );
    }
    else {
        $cnaseg_cmd = Genome::Model::Tools::CopyNumber::CnaSeg->create(
            copy_number_file => $data_paths{cnvs_hq},
            min_markers      => 4,
            detect_somatic   => 1,
            centromere_file  => $centromere_file,
            gap_file         => $gap_file,
            output_file      => $cnvhmm_file
        );
    }
    $cnaseg_cmd->execute();

    my $varscan_file = $readcounts_varscan_file;
    #my $varscan_file = $readcounts_clonality_outfile;

    #Step 7 - then, put the cna-seg and varscan-format snv file together in this clonality tool:
    #gmt validation clonality-plot     --cnvhmm-file     /gscuser/ndees/103/wgs/SV_somatic/CNV/aml103.cnvhmm     --output-image     aml103.clonality.pdf     --r-script-output-file     clonality.R     --varscan-file     allsnvs.hq.novel.tier123.v2.bed.adapted.readcounts.varscan     --analysis-type     wgs     --sample-id     'AML103'     --positions-highlight     IL2RA_NF1_positions

    #gmt validation clonality-plot  --cnvhmm-file='/gscmnt/sata132/techd/mgriffit/hg1/clonality/hg1.cnvhmm'  --output-image hg1.clonality.pdf  --r-script-output-file clonality.R  --varscan-file allsnvs.hq.novel.tier123.v2.bed.adapted.readcounts.varscan  --analysis-type wgs  --sample-id 'HG1'

    #Step 7-A. Without clusters
    my $output_image_file1a = "$output_dir" . "$common_name" . ".clonality.pdf";
    my $r_script_file       = "$output_dir" . "clonality.R";
    my $uc_common_name      = uc($common_name);
    if (-s $varscan_file) {
        my $clonality_cmd1a = Genome::Model::Tools::Validation::ClonalityPlot->create(
            cnvhmm_file          => $cnvhmm_file,
            output_image         => $output_image_file1a,
            r_script_output_file => $r_script_file,
            varscan_file         => $varscan_file,
            analysis_type        => 'wgs',
            sample_id            => $uc_common_name
        );
        $clonality_cmd1a->execute();

        my $output_image_file1b = "$output_dir" . "$common_name" . ".clonality.cn2.pdf";
        my $clonality_cmd1b     = Genome::Model::Tools::Validation::ClonalityPlot->create(
            cnvhmm_file          => $cnvhmm_file,
            output_image         => $output_image_file1b,
            r_script_output_file => $r_script_file,
            varscan_file         => $varscan_file,
            analysis_type        => 'wgs',
            sample_id            => $uc_common_name,
            plot_only_cn2        => 1
        );
        $clonality_cmd1b->execute();
    }
    else {
        $self->warning_message("Empty snv file: $varscan_file");
    }
    #Step 7-B. With clusters - not working - deactivate until the new version of clonality is available
    #my $clustered_data_output_file = $output_dir . $common_name . ".clustered.data.tsv";
    #my $output_image_file2 = "$output_dir"."$common_name".".clonality.clusters.pdf";
    #my $clonality_cmd2 = "gmt validation clonality-plot  --cnvhmm-file=$cnvhmm_file  --output-image=$output_image_file2  --r-script-output-file=$r_script_file  --varscan-file=$varscan_file  --analysis-type=wgs  --sample-id='$uc_common_name'  --plot-clusters  --clustered-data-output-file=$clustered_data_output_file";
    #if ($verbose){$self->debug_message("$clonality_cmd2");}
    #TODO: until the --plot-clusters functionality is more stable, allow a failed exit code
    #Genome::Sys->shellcmd(cmd => $clonality_cmd2, allow_failed_exit_code => 1);

    #Step 7-C.  Without clusters but after filtering the input file to remove those with low coverage in tumor or high VAF in normal
    #normal coverage = col5+col6; normal vaf = col7; tumor coverage = col9+10; tumor vaf = col11
    my $filtered_file = $varscan_file . ".filt";

    my $varscan_fh  = Genome::Sys->open_file_for_reading($varscan_file);
    my $filtered_fh = Genome::Sys->open_file_for_writing($filtered_file);
    while (<$varscan_fh>) {
        my @line       = split("\t", $_);
        my $normal_cov = $line[4] + $line[5];
        my $normal_vaf = $line[6];
        my $tumor_cov  = $line[8] + $line[9];
        my $tumor_vaf  = $line[10];
        next
            if ($tumor_cov < $self->min_tumor_cov
            || $tumor_cov > $self->max_tumor_cov
            || $normal_vaf > $self->max_normal_vaf);
        print $filtered_fh $_;
    }
    $varscan_fh->close;
    $filtered_fh->close;

    if (-s $filtered_file) {
        my $output_image_file3a = "$output_dir" . "$common_name" . ".clonality.filtered_snvs.pdf";
        my $clonality_cmd3a     = Genome::Model::Tools::Validation::ClonalityPlot->create(
            cnvhmm_file          => $cnvhmm_file,
            output_image         => $output_image_file3a,
            r_script_output_file => $r_script_file,
            varscan_file         => $filtered_file,
            analysis_type        => 'wgs',
            sample_id            => $uc_common_name
        );
        $clonality_cmd3a->execute();

        my $output_image_file3b = "$output_dir" . "$common_name" . ".clonality.filtered_snvs.cn2.pdf";
        my $clonality_cmd3b     = Genome::Model::Tools::Validation::ClonalityPlot->create(
            cnvhmm_file          => $cnvhmm_file,
            output_image         => $output_image_file3b,
            r_script_output_file => $r_script_file,
            varscan_file         => $filtered_file,
            analysis_type        => 'wgs',
            sample_id            => $uc_common_name,
            plot_only_cn2        => 1
        );
        $clonality_cmd3b->execute();
    }
    else {
        $self->warning_message("Empty filtered snv file: $filtered_file");
    }

    #Keep the files that were needed to run the cna-seg and clonality plot steps so that someone can rerun with different parameters

    #Define the output parameter for the cnvhmm so that it can be fed into downstream steps in a workflow
    if (not $is_copycat) {
        $self->cnv_hmm_file($cnvhmm_file);
    }
    else {
        my $copycat_cnvhmm_file = $output_dir . "cnaseq.cnvhmm";
        $self->create_copycat_cnvhmm_file($somatic_var_build, $copycat_cnvhmm_file);
        $self->cnv_hmm_file($copycat_cnvhmm_file);
    }
    $self->cnv_hq_file($cnvs_output_path);

    return 1;
}

sub get_data_paths {
    my $self              = shift;
    my $data_paths        = shift;
    my $somatic_var_build = shift;
    my $is_copycat        = shift;
    my $output_dir        = shift;
    $data_paths->{root_dir}    = $somatic_var_build->data_directory . "/";
    $data_paths->{effects_dir} = $data_paths->{root_dir} . "effects/";
    if (not $is_copycat) {
        $data_paths->{cnvs_hq} = $data_paths->{root_dir} . "variants/cnvs.hq";
    }
    else {
        $data_paths->{cnvs_hq} = $self->create_copycat_cnvhq_file($somatic_var_build, $output_dir);
    }
    $data_paths->{normal_bam} = $somatic_var_build->normal_bam;
    $data_paths->{tumor_bam}  = $somatic_var_build->tumor_bam;
    my $reference_build = $somatic_var_build->reference_sequence_build;
    $data_paths->{reference_fasta} = $reference_build->full_consensus_path('fa');
    $data_paths->{display_name}    = $reference_build->__display_name__;
}

sub create_pseudo_varscan_file {
    my $self            = shift;
    my $snvs_file       = shift;
    my $readcounts_file = shift;
    my $output_file     = shift;

    my $output_fh = Genome::Sys->open_file_for_writing($output_file);

    my %tumHash;
    my $rc_fh = Genome::Sys->open_file_for_reading($readcounts_file);
    while (my $line = $rc_fh->getline) {
        chomp($line);
        my @fields = split("\t", $line);
        $tumHash{$fields[0] . "|" . $fields[1]} = $line;
    }
    $rc_fh->close;

    my $snvs_fh = Genome::Sys->open_file_for_reading($snvs_file);
    while (my $line = $snvs_fh->getline) {
        chomp($line);
        my @fields = split("\t", $line);
        if (exists($tumHash{$fields[0] . "|" . $fields[1]})) {
            my @tum = split("\t", $tumHash{$fields[0] . "|" . $fields[1]});

            print $output_fh join("\t",
                @fields[0 .. 1, 3 .. 4],
                @tum[10 .. 12],
                "NULL", @tum[5 .. 7],
                "NULL", "Somatic", (("NULL") x 6)),
                "\n";

        }
    }
    $snvs_fh->close;

    return 1;
}

1;

