//! DATES for Python: admixture dates from the decay of weighted LD.

use crate::compute;
use crate::genotypes::read_populations;
use popmix_core::python::input_error;
use pyo3::prelude::*;
use pyo3::types::PyDict;

/// Result of a single LD bin from DATES.
#[pyclass(module = "popmix._core", skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyLdBin {
    #[pyo3(get)]
    distance_cm: f64,
    #[pyo3(get)]
    weighted_ld: f64,
    #[pyo3(get, name = "z")]
    z_score: f64,
    #[pyo3(get)]
    correlation: f64,
    #[pyo3(get)]
    regression_coef: f64,
    #[pyo3(get)]
    n_pairs: usize,
    #[pyo3(get, name = "se")]
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
#[pyclass(module = "popmix._core", skip_from_py_object)]
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
#[pyclass(module = "popmix._core")]
pub(crate) struct PyDatesResult {
    #[pyo3(get, name = "admixed")]
    test_pop: String,
    source1: String,
    source2: String,
    #[pyo3(get)]
    generations: f64,
    #[pyo3(get)]
    generations_se: f64,
    #[pyo3(get)]
    years_bp: f64,
    #[pyo3(get)]
    years_bp_se: f64,
    /// The assumptions behind the years: years per generation and the
    /// samples' age (years BP), with their SEs.
    #[pyo3(get)]
    generation_time: f64,
    #[pyo3(get)]
    generation_time_se: f64,
    #[pyo3(get)]
    sample_age: f64,
    #[pyo3(get)]
    sample_age_se: f64,
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
    /// The two source populations.
    #[getter]
    fn sources(&self) -> Vec<String> {
        vec![self.source1.clone(), self.source2.clone()]
    }

    /// The LD bins as a pandas DataFrame, one row per bin.
    fn to_ld_dataframe<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        let rows = self
            .ld_bins
            .iter()
            .map(|b| {
                let d = PyDict::new(py);
                d.set_item("distance_cm", b.distance_cm)?;
                d.set_item("n_pairs", b.n_pairs)?;
                d.set_item("weighted_ld", b.weighted_ld)?;
                d.set_item("se", b.std_err)?;
                d.set_item("z", b.z_score)?;
                d.set_item("correlation", b.correlation)?;
                d.set_item("regression_coef", b.regression_coef)?;
                Ok(d)
            })
            .collect::<PyResult<Vec<_>>>()?;
        py.import("pandas")?.getattr("DataFrame")?.call1((rows,))
    }

    fn __repr__(&self) -> String {
        format!(
            "DatesResult(admixed='{}', generations={:.1}±{:.1}, years_bp={:.0}±{:.0}, n_snps={})",
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
///     prefix: The genotype dataset (only the three populations' samples
///         are used; their ploidy is detected as for every analysis)
///     admixed: The admixed population (DATES's test population)
///     sources: Its two source populations
///     generation_time: Generation time in years (default: 29)
///     generation_time_se: Its standard error in years (default: 0)
///     sample_age: Age of the admixed samples in years BP (default: 0),
///         so that ``years_bp`` is before present, not before they lived
///     sample_age_se: Its standard error in years (default: 0)
///     min_distance_cm: Minimum distance for fitting, DATES lovalfit (default: 0.45)
///     max_distance_cm: Maximum distance, DATES maxdis x 100 (default: 100.0)
///     bin_width_cm: Width of distance bins, DATES binsize x 100 (default: 0.1)
///     jackknife: Leave-one-chromosome-out jackknife for the SE (default: True)
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
///     >>> result = pm.dates("/path/to/data", "French", ["Yamnaya", "Anatolia_N"])
///     >>> print(f"Admixture time: {result.generations:.0f} generations")
#[pyfunction]
#[pyo3(signature = (prefix, admixed, sources, *, generation_time=29.0, generation_time_se=0.0, sample_age=0.0, sample_age_se=0.0, min_distance_cm=0.45, max_distance_cm=100.0, bin_width_cm=0.1, jackknife=true, qbin=0, compat=None))]
#[allow(clippy::too_many_arguments)]
pub(crate) fn dates(
    py: Python<'_>,
    prefix: &str,
    admixed: &str,
    sources: Vec<String>,
    generation_time: f64,
    generation_time_se: f64,
    sample_age: f64,
    sample_age_se: f64,
    min_distance_cm: f64,
    max_distance_cm: f64,
    bin_width_cm: f64,
    jackknife: bool,
    qbin: usize,
    compat: Option<&str>,
) -> PyResult<PyDatesResult> {
    let [source1, source2] = sources.as_slice() else {
        return Err(input_error(format!(
            "DATES takes two sources, got {}",
            sources.len()
        )));
    };
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
        let mut pops = vec![admixed.to_string()];
        for p in [source1, source2] {
            if !pops.contains(p) {
                pops.push(p.clone());
            }
        }
        let geno = read_populations(prefix, &pops)?;

        let mut opts = popmix_core::DatesOptions::default();
        opts.ld_options.qbin = qbin;
        opts.ld_options.bin_width_cm = bin_width_cm;
        opts.ld_options.max_distance_cm = max_distance_cm;
        opts.fit_options.bin_size_cm = bin_width_cm;
        opts.fit_options.min_distance_cm = min_distance_cm;
        opts.fit_options.max_distance_cm = max_distance_cm;
        opts.generation_time = generation_time;
        opts.generation_time_se = generation_time_se;
        opts.sample_age = sample_age;
        opts.sample_age_se = sample_age_se;
        opts.run_jackknife = jackknife;
        opts.compat = compat;

        let result = popmix_core::dates(&geno, None, admixed, source1, source2, &opts)?;

        Ok(PyDatesResult {
            test_pop: result.test_pop.to_string(),
            source1: result.source1.to_string(),
            source2: result.source2.to_string(),
            generations: result.generations,
            generations_se: result.generations_se,
            years_bp: result.years_bp,
            years_bp_se: result.years_bp_se,
            generation_time,
            generation_time_se,
            sample_age,
            sample_age_se,
            n_snps: result.n_snps,
            n_chromosomes: result.n_chromosomes,
            ld_bins: result.ld_bins,
            fit: result.fit,
            jackknife: result.jackknife_by_chromosome,
        })
    })
}
