use anyhow::{anyhow, Result};
use std::path::PathBuf;
use structopt::StructOpt;

#[derive(StructOpt, Debug)]
#[structopt(
    about = "A tool for identifying obiquitous variants in bcf files.",
    name = "ubiquitous-variant-identifier"
)]
pub struct UbiquitousVariantIdentifier {
    /// bcf files to be inspected.
    #[structopt(long, short = "b", parse(from_os_str))]
    pub(crate) bcf_paths: Vec<PathBuf>,

    /// Threshold for minimum number of variant occurences to be considered ubiquitous.
    #[structopt(long, short = "t", default_value = "0.7")]
    pub(crate) threshold: f32,
}

pub(crate) trait ValidateArguments {
    fn validate(&self) -> Result<()>;
}

impl ValidateArguments for UbiquitousVariantIdentifier {
    fn validate(&self) -> Result<()> {
        if self.bcf_paths.len() < 2 {
            return Err(anyhow!(
                "You need to specify at least two bcf files to analyse."
            ));
        }
        if !(0.0..1.0).contains(&self.threshold) {
            return Err(anyhow!("Threshold must be a float between 0 and 1."));
        }
        Ok(())
    }
}
