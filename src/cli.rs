// src/cli.rs
// Argument parser and CLI (Command Line Interface)
// for BreakFast

use clap::{App, Arg, SubCommand, AppSettings};

pub fn build_cli() -> App<'static, 'static> {
    App::new("Breakfast")
        .version("0.1")
        .about("\nBreakfast is a toolkit for detecting chromosomal rearrangements based on whole genome sequencing data.")
        .setting(AppSettings::ArgRequiredElseHelp)
        .setting(AppSettings::DisableHelpSubcommand)
        .setting(AppSettings::UnifiedHelpMessage)
        .setting(AppSettings::VersionlessSubcommands)
        .setting(AppSettings::SubcommandRequiredElseHelp)

        .subcommand(SubCommand::with_name("detect")
            .arg(Arg::with_name("bam_file")
                .required(true).takes_value(true)
            	.help("BAM file that will be analyzed for rearrangements"))
            .arg(Arg::with_name("genome")
                .required(true).takes_value(true)
            	.help("Path to a Bowtie index"))
            .arg(Arg::with_name("anchor-len")
                .short("a").long("anchor-len")
                .takes_value(true).value_name("N").default_value("30")
                .help("Anchor length for split read analysis"))
            .arg(Arg::with_name("max-frag-len")
                .short("f").long("max-frag-len")
                .takes_value(true).value_name("N").default_value("5000")
                .help("Maximum fragment length")))

        .subcommand(SubCommand::with_name("filter")
            .arg(Arg::with_name("sv_file")
                .required(true).takes_value(true))
            .arg(Arg::with_name("blacklist")
            	.long("blacklist").value_name("path").default_value("")
                .takes_value(true))
            .arg(Arg::with_name("min-reads")
                .short("r").long("min-reads")
                .required(true)
                .takes_value(true).value_name("N")
                .help("Minimum number of breakpoint-spanning reads required to accept a rearrangement")))

        .subcommand(SubCommand::with_name("annotate")
            .arg(Arg::with_name("sv_file")
                .required(true)
                .takes_value(true))
            .arg(Arg::with_name("bed_file")
                .required(true)
                .takes_value(true)))

        .subcommand(SubCommand::with_name("blacklist")
            .arg(Arg::with_name("sv_files")
                .required(true)
                .takes_value(true)
                .multiple(true))
            .arg(Arg::with_name("freq-above")
                .long("freq-above")
                .takes_value(true).value_name("FREQ")
                .help("Minimum frequency at which a variant must be present among the control samples to be considered a false positive [default: 0].")))

}
