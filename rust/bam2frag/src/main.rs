use anyhow::Result;
use bam2frag::{
    cli::{Cli, Command},
    config, output, pipeline,
};
use clap::Parser;

fn main() -> Result<()> {
    tracing_subscriber::fmt()
        .with_env_filter(
            tracing_subscriber::EnvFilter::try_from_default_env()
                .unwrap_or_else(|_| tracing_subscriber::EnvFilter::new("info")),
        )
        .with_ansi(false)
        .init();
    let cli = Cli::parse();
    match cli.command {
        Command::Extract(args) => pipeline::run_extract(args),
        Command::ValidateConfig(args) => config::run_validate(args),
        Command::Inspect(args) => output::inspect_output(args),
        Command::Version => {
            println!("{}", env!("CARGO_PKG_VERSION"));
            Ok(())
        }
    }
}
