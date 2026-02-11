use clap::Parser;

/// Run two-tier chi-squared homogeneity test on genotype data.
#[derive(Parser, Debug)]
#[command(name = "chi2mc", about)]
pub struct Args {
    /// Path to the input TSV file containing genotype data
    pub input_file: String,

    /// Path to the output TSV file for results
    pub output_file: String,

    /// Minimum expected count threshold for using standard test
    #[arg(long, default_value_t = 5)]
    pub min_expected_count: usize,

    /// Number of Monte Carlo simulations for small samples
    #[arg(long, default_value_t = 10000)]
    pub n_simulations: usize,

    /// Random seed for reproducibility
    #[arg(long, default_value_t = 42)]
    pub random_seed: u64,

    /// Number of parallel workers to use (default: all CPU cores)
    #[arg(long)]
    pub n_workers: Option<usize>,

    /// Show detailed progress with chi-squared and p-values for each probeset
    #[arg(long)]
    pub verbose_progress: bool,

    /// Number of probesets per batch (kept for compatibility)
    #[arg(long, default_value_t = 1000)]
    pub batch_size: usize,
}
