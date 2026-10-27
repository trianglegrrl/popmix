//! DATES for Python: admixture dates from the decay of weighted LD.

use crate::compute;
use crate::genotypes::PyGenotypeMatrix;
use popmix_core::python::input_error;
use pyo3::prelude::*;
use pyo3::types::PyDict;

/// Result of a single LD bin from DATES.
#[pyclass(skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyLdBin {
    #[pyo3(get)]
    distance_cm: f64,
    #[pyo3(get)]
    weighted_ld: f64,
    #[pyo3(get)]
    z_score: f64,
    #[pyo3(get)]
    correlation: f64,
    #[pyo3(get)]
    regression_coef: f64,
    #[pyo3(get)]
    n_pairs: usize,
    #[pyo3(get)]
    std_err: f64,
}

impl From<&popmix_core::LdBin> for PyLdBin {
    fn from(b: &popmix_core::LdBin) -> Self {
        PyLdBin {
            distance_cm: b.distance_cm,
            weighted_ld: b.weighted_ld,
            z_score: b.z_score,
            correlation: b.correlation,
            regression_coef: b.regression_coef,
            n_pairs: b.n_pairs,
            std_err: b.std_err,
        }
    }
}

/// Result of exponential decay fitting.
#[pyclass(skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyDecayFit {
    #[pyo3(get)]
    amplitude: f64,
    #[pyo3(get)]
    decay_rate: f64,
    #[pyo3(get)]
    intercept: f64,
    #[pyo3(get)]
    r_squared: f64,
    #[pyo3(get)]
    decay_rate_se: f64,
    #[pyo3(get)]
    rss: f64,
}

impl From<&popmix_core::DecayFit> for PyDecayFit {
    fn from(f: &popmix_core::DecayFit) -> Self {
        PyDecayFit {
            amplitude: f.amplitude,
            decay_rate: f.decay_rate,
            intercept: f.intercept,
            r_squared: f.r_squared,
            decay_rate_se: f.decay_rate_se,
            rss: f.rss,
        }
    }
}

/// Result of DATES admixture dating analysis.
#[pyclass]
pub(crate) struct PyDatesResult {
    #[pyo3(get)]
    test_pop: String,
    #[pyo3(get)]
    source1: String,
    #[pyo3(get)]
    source2: String,
    #[pyo3(get)]
    generations: f64,
    #[pyo3(get)]
    generations_se: f64,
    #[pyo3(get)]
    years_bp: f64,
    #[pyo3(get)]
    years_bp_se: f64,
    #[pyo3(get)]
    n_snps: usize,
    #[pyo3(get)]
    n_chromosomes: usize,
    ld_bins: Vec<popmix_core::LdBin>,
    fit: popmix_core::DecayFit,
    jackknife: Vec<popmix_core::dates::JackknifeBlock>,
}

#[pymethods]
impl PyDatesResult {
    /// Per-chromosome jackknife blocks as (chromosome, weight, generations),
    /// matching the DATES `.jin` file. Empty if the jackknife did not run.
    #[getter]
    fn jackknife_by_chromosome(&self) -> Vec<(u8, f64, f64)> {
        self.jackknife
            .iter()
            .map(|b| (b.chromosome, b.weight, b.generations))
            .collect()
    }

    /// Get the LD decay bins.
    #[getter]
    fn ld_bins(&self) -> Vec<PyLdBin> {
        self.ld_bins.iter().map(PyLdBin::from).collect()
    }

    /// Get the exponential fit result.
    #[getter]
    fn fit(&self) -> PyDecayFit {
        PyDecayFit::from(&self.fit)
    }

    /// Get LD decay data as a list of dictionaries for easy pandas conversion.
    fn to_ld_dataframe<'py>(&self, py: Python<'py>) -> PyResult<Vec<Bound<'py, PyDict>>> {
        self.ld_bins
            .iter()
            .map(|b| {
                let d = PyDict::new(py);
                d.set_item("distance_cm", b.distance_cm)?;
                d.set_item("weighted_ld", b.weighted_ld)?;
                d.set_item("z_score", b.z_score)?;
                d.set_item("correlation", b.correlation)?;
                d.set_item("regression_coef", b.regression_coef)?;
                d.set_item("n_pairs", b.n_pairs as f64)?;
                d.set_item("std_err", b.std_err)?;
                Ok(d)
            })
            .collect()
    }

    fn __repr__(&self) -> String {
        format!(
            "DatesResult(test_pop='{}', generations={:.1}±{:.1}, years_bp={:.0}±{:.0}, n_snps={})",
            self.test_pop,
            self.generations,
            self.generations_se,
            self.years_bp,
            self.years_bp_se,
            self.n_snps
        )
    }
}

/// Estimate admixture time using the DATES algorithm.
///
/// DATES estimates when admixture occurred by measuring the decay of
/// linkage disequilibrium with genetic distance. LD from recent
/// admixture decays exponentially, with the decay rate proportional
/// to the number of generations since admixture.
///
/// Args:
///     geno: Genotype matrix with all relevant populations
///     test_pop: Name of the test (admixed) population
///     source1: Name of first source population
///     source2: Name of second source population
///     generation_time: Generation time in years (default: 29)
///     min_distance_cm: Minimum distance for fitting, DATES lovalfit (default: 0.45)
///     max_distance_cm: Maximum distance, DATES maxdis x 100 (default: 100.0)
///     bin_width_cm: Width of distance bins, DATES binsize x 100 (default: 0.1)
///     run_jackknife: Run jackknife for SE estimation (default: True)
///     qbin: DATES qbin. 0 visits every SNP pair; > 0 uses an FFT mesh of
///         qbin cells per bin, which is faster and slightly approximate
///         (default: 0, the DATES default)
///     compat: Reproduce a reference implementation's known quirks instead
///         of correcting them. None (default) corrects them all;
///         "dates-v753" gives output identical to DATES v753, for checking
///         against it. See DEVIATIONS.md for what each quirk changes.
///
/// Returns:
///     DatesResult with estimated admixture timing. ``generations`` and
///     ``generations_se`` are the weighted jackknife estimate and SE that
///     DATES reports; ``fit.decay_rate`` is the full-data fit.
///
/// Example:
///     >>> geno = pm.read_packedancestrymap("/path/to/data")
///     >>> result = pm.dates(geno, "French", "Yamnaya", "Anatolia_N")
///     >>> print(f"Admixture time: {result.generations:.0f} generations")
///     >>> print(f"Years BP: {result.years_bp:.0f}")
#[pyfunction]
#[pyo3(signature = (geno, test_pop, source1, source2, generation_time=None, min_distance_cm=None, max_distance_cm=None, bin_width_cm=None, run_jackknife=None, qbin=None, compat=None))]
#[allow(clippy::too_many_arguments)]
pub(crate) fn dates(
    py: Python<'_>,
    geno: &PyGenotypeMatrix,
    test_pop: &str,
    source1: &str,
    source2: &str,
    generation_time: Option<f64>,
    min_distance_cm: Option<f64>,
    max_distance_cm: Option<f64>,
    bin_width_cm: Option<f64>,
    run_jackknife: Option<bool>,
    qbin: Option<usize>,
    compat: Option<&str>,
) -> PyResult<PyDatesResult> {
    compute(py, move || {
        let compat = match compat {
            None => popmix_core::dates::DatesCompat::CORRECTED,
            Some("dates-v753") => popmix_core::dates::DatesCompat::V753,
            Some(other) => {
                return Err(input_error(format!(
                    "unknown compat {other:?} for dates; expected None or \"dates-v753\""
                )))
            }
        };
        let mut ld_options = popmix_core::LdOptions::default();
        if let Some(q) = qbin {
            ld_options.qbin = q;
        }
        if let Some(bw) = bin_width_cm {
            ld_options.bin_width_cm = bw;
        }
        if let Some(max_d) = max_distance_cm {
            ld_options.max_distance_cm = max_d;
        }

        let mut fit_options = popmix_core::dates::FitOptions::default();
        if let Some(bw) = bin_width_cm {
            fit_options.bin_size_cm = bw;
        }
        if let Some(min_d) = min_distance_cm {
            fit_options.min_distance_cm = min_d;
        }
        if let Some(max_d) = max_distance_cm {
            fit_options.max_distance_cm = max_d;
        }

        let opts = {
            let mut o = popmix_core::DatesOptions::default();
            o.ld_options = ld_options;
            o.fit_options = fit_options;
            o.generation_time = generation_time.unwrap_or(29.0);
            o.run_jackknife = run_jackknife.unwrap_or(true);
            o.compat = compat;
            o
        };

        let result = popmix_core::dates(&geno.inner, None, test_pop, source1, source2, &opts)?;

        Ok(PyDatesResult {
            test_pop: result.test_pop.to_string(),
            source1: result.source1.to_string(),
            source2: result.source2.to_string(),
            generations: result.generations,
            generations_se: result.generations_se,
            years_bp: result.years_bp,
            years_bp_se: result.years_bp_se,
            n_snps: result.n_snps,
            n_chromosomes: result.n_chromosomes,
            ld_bins: result.ld_bins,
            fit: result.fit,
            jackknife: result.jackknife_by_chromosome,
        })
    })
}
