mod cli;
mod io;
mod monte_carlo;
mod probeset;
mod report;
mod stats;

use std::path::Path;
use std::time::Instant;

use clap::Parser;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;

use cli::Args;
use probeset::ProbesetResult;

fn main() {
    let args = Args::parse();

    println!("Two-Tier Chi-Squared Homogeneity Test (OPTIMIZED)");
    println!("{}", "=".repeat(80));
    println!("Input file: {}", args.input_file);
    println!("Output file: {}", args.output_file);

    let n_cpus = num_cpus::get();
    println!("Available CPU cores: {n_cpus}");
    let n_workers = args.n_workers.unwrap_or(n_cpus);
    println!("Processing probesets in parallel using {n_workers} worker(s)...");
    println!("Rust Monte Carlo engine: ENABLED (native binary)");

    if args.show_progress || args.verbose_progress {
        println!("Progress tracking: ENABLED (using indicatif progress bar)");
    }

    println!("Note: Using {n_workers} parallel workers for processing");
    println!("Note: Batch size = {} probesets per worker", args.batch_size);
    println!();

    // Configure rayon thread pool
    rayon::ThreadPoolBuilder::new()
        .num_threads(n_workers)
        .build_global()
        .ok();

    // Load and group input data
    let groups = match io::read_input(Path::new(&args.input_file)) {
        Ok(g) => g,
        Err(e) => {
            eprintln!("Error reading input file: {e}");
            std::process::exit(1);
        }
    };

    let n_total = groups.len();
    println!("Loaded {n_total} probesets from input file.");
    println!();

    // Collect into a Vec for parallel iteration
    let probeset_entries: Vec<_> = groups.into_iter().collect();

    let start = Instant::now();

    // Optional progress bar
    let pb = if args.show_progress || args.verbose_progress {
        let bar = ProgressBar::new(n_total as u64);
        bar.set_style(
            ProgressStyle::with_template(
                "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta})",
            )
            .unwrap()
            .progress_chars("#>-"),
        );
        Some(bar)
    } else {
        None
    };

    // Process probesets in parallel with rayon
    let results: Vec<ProbesetResult> = probeset_entries
        .par_iter()
        .map(|(probeset_id, records)| {
            let result = probeset::process_probeset(
                probeset_id,
                records,
                args.min_expected_count,
                args.n_simulations,
                args.random_seed,
            );
            if let Some(ref bar) = pb {
                bar.inc(1);
                if args.verbose_progress {
                    let msg = format_progress_message(&result);
                    bar.println(msg);
                }
            }
            result
        })
        .collect();

    if let Some(ref bar) = pb {
        bar.finish_with_message("done");
    }

    let elapsed = start.elapsed().as_secs_f64();

    // Sort results by probeset_id
    let mut results = results;
    results.sort_by(|a, b| a.probeset_id.cmp(&b.probeset_id));

    // Display results summary
    println!("Results:");
    println!("{}", "=".repeat(80));
    println!(
        "shape: ({}, 10)\n",
        results.len()
    );

    // Write output TSV
    let output_path = Path::new(&args.output_file);
    if let Err(e) = io::write_output(output_path, &results) {
        eprintln!("Error writing output file: {e}");
        std::process::exit(1);
    }
    println!("Results saved to: {}", args.output_file);

    // Console summary
    println!("\n{}", "=".repeat(80));
    println!("Summary:");
    println!("Total probesets tested: {}", results.len());

    let n_standard = results.iter().filter(|r| r.method == "standard").count();
    println!("Standard chi-squared tests: {n_standard}");

    let n_mc = results.iter().filter(|r| r.method == "monte_carlo").count();
    println!("Monte Carlo tests: {n_mc}");

    // Excluded genotypes
    let excluded: Vec<_> = results
        .iter()
        .filter(|r| r.excluded_genotypes.is_some())
        .collect();
    if !excluded.is_empty() {
        println!("Probesets with excluded genotypes: {}", excluded.len());
        for (i, r) in excluded.iter().enumerate() {
            if i < 10 {
                println!(
                    "  - {}: excluded {} (all zeros)",
                    r.probeset_id,
                    r.excluded_genotypes.as_deref().unwrap()
                );
            } else if i == 10 {
                println!(
                    "  ... and {} more (see summary report)",
                    excluded.len() - 10
                );
                break;
            }
        }
    }

    // Skipped probesets
    let skipped: Vec<_> = results.iter().filter(|r| r.method == "skipped").collect();
    if !skipped.is_empty() {
        println!("Skipped probesets (errors): {}", skipped.len());
        for r in &skipped {
            println!(
                "  - {}: {}",
                r.probeset_id,
                r.error.as_deref().unwrap_or("Unknown error")
            );
        }
    }

    // Significant results
    let valid: Vec<_> = results.iter().filter(|r| r.method != "skipped").collect();
    if !valid.is_empty() {
        let n_significant = valid
            .iter()
            .filter(|r| r.pvalue.map_or(false, |p| p < 0.05))
            .count();
        println!("Significant results (p < 0.05): {n_significant}");
    }

    println!("Processing time: {elapsed:.2} seconds");

    // Generate summary report
    let summary_file = args.output_file.replace(".tsv", "_summary.txt");
    if let Err(e) =
        report::write_summary_report(Path::new(&summary_file), &args.input_file, &args.output_file, elapsed, &results)
    {
        eprintln!("Error writing summary report: {e}");
    } else {
        println!("\nDetailed summary report saved to: {summary_file}");
    }
}

fn format_progress_message(result: &ProbesetResult) -> String {
    if result.method == "skipped" {
        return format!(
            "{}: SKIPPED - {}",
            result.probeset_id,
            result.error.as_deref().unwrap_or("Unknown error")
        );
    }

    let chi2 = result.chi2_statistic.unwrap_or(f64::NAN);
    let pval = result.pvalue.unwrap_or(f64::NAN);

    let pval_str = if pval < 0.0001 {
        format!("{pval:.2e}")
    } else {
        format!("{pval:.6}")
    };

    let sig_marker = if pval < 0.001 {
        " ***"
    } else if pval < 0.01 {
        " **"
    } else if pval < 0.05 {
        " *"
    } else {
        ""
    };

    format!(
        "{}: \u{03c7}\u{00b2}={chi2:.4}, p={pval_str}{sig_marker} ({})",
        result.probeset_id, result.method,
    )
}
