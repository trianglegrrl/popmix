//! Thread pool configuration for parallel computation.
//!
//! Popmix computes on its own pool of threads, 4 by default: entry points
//! run their work through [`install`]. Work run outside it uses rayon's
//! global pool.

use std::sync::{Arc, RwLock};

use crate::error::{PopmixError, Result};

/// Default number of threads for parallel operations.
pub const DEFAULT_THREADS: usize = 4;

/// Most threads a pool may have. rayon allocates for every thread it is
/// asked for, so a mistyped count would exhaust memory and abort.
pub const MAX_THREADS: usize = 1024;

/// The pool popmix computes on, and its number of threads.
static POOL: RwLock<Option<(usize, Arc<rayon::ThreadPool>)>> = RwLock::new(None);

/// Use `threads` threads (default [`DEFAULT_THREADS`]) for the computations
/// run through [`install`], replacing the pool used so far; it can be called
/// at any time.
///
/// # Returns
///
/// `Ok(true)`; an error if the count is 0 or above [`MAX_THREADS`], or the
/// threads cannot be started.
///
/// # Example
///
/// ```
/// use popmix_core::init_thread_pool;
///
/// init_thread_pool(Some(8)).unwrap();
/// assert_eq!(popmix_core::threading::install(rayon::current_num_threads), 8);
/// ```
pub fn init_thread_pool(threads: Option<usize>) -> Result<bool> {
    let n = threads.unwrap_or(DEFAULT_THREADS);
    if n == 0 || n > MAX_THREADS {
        return Err(PopmixError::InvalidThreadCount(n));
    }
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(n)
        .build()
        .map_err(|e| PopmixError::ThreadPool(format!("cannot start {n} threads: {e}")))?;

    *POOL
        .write()
        .unwrap_or_else(std::sync::PoisonError::into_inner) = Some((n, Arc::new(pool)));
    Ok(true)
}

/// The number of threads popmix computes with ([`DEFAULT_THREADS`] until
/// [`init_thread_pool`] sets it); asking starts no thread.
pub fn current_threads() -> usize {
    POOL.read()
        .unwrap_or_else(std::sync::PoisonError::into_inner)
        .as_ref()
        .map_or(DEFAULT_THREADS, |(n, _)| *n)
}

/// Run `f` on popmix's pool, so its parallel work uses the configured number
/// of threads (the pool is started with the default on first use).
pub fn install<R: Send>(f: impl FnOnce() -> R + Send) -> R {
    let configured = POOL
        .read()
        .unwrap_or_else(std::sync::PoisonError::into_inner)
        .as_ref()
        .map(|(_, p)| Arc::clone(p));
    let Some(pool) = configured else {
        // a failure to start threads surfaces on the next explicit init
        if init_thread_pool(None).is_err() {
            return f();
        }
        return install(f);
    };

    pool.install(f)
}

/// `zero() + part(0) + part(1) + ... + part(n_groups - 1)`, in that order,
/// computing `threads` parts at a time.
pub fn ordered_sum<T: Send>(
    n_groups: usize,
    zero: impl Fn() -> T,
    part: impl Fn(usize) -> T + Sync,
    add: impl Fn(&mut T, T),
) -> T {
    use rayon::prelude::*;
    let mut acc = zero();
    let batch = rayon::current_num_threads().max(1);
    for start in (0..n_groups).step_by(batch) {
        let parts: Vec<T> = (start..(start + batch).min(n_groups))
            .into_par_iter()
            .map(&part)
            .collect();
        for p in parts {
            add(&mut acc, p);
        }
    }
    acc
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default_threads() {
        assert_eq!(DEFAULT_THREADS, 4);
    }

    #[test]
    fn test_invalid_thread_count() {
        let result = init_thread_pool(Some(0));
        assert!(matches!(result, Err(PopmixError::InvalidThreadCount(0))));
    }

    #[test]
    fn too_many_threads_is_an_error_before_any_is_started() {
        for n in [MAX_THREADS + 1, usize::MAX] {
            assert!(matches!(
                init_thread_pool(Some(n)),
                Err(PopmixError::InvalidThreadCount(_))
            ));
        }
    }

    #[test]
    fn computations_run_on_the_configured_threads_and_it_can_change() {
        // the only test that configures the pool, as it is global
        assert!(init_thread_pool(Some(3)).unwrap());
        assert_eq!(current_threads(), 3);
        assert_eq!(install(rayon::current_num_threads), 3);
        assert!(init_thread_pool(Some(5)).unwrap());
        assert_eq!(install(rayon::current_num_threads), 5);
        // nested use stays on the pool
        assert_eq!(install(|| install(rayon::current_num_threads)), 5);
    }

    #[test]
    fn test_current_threads() {
        // Just verify it returns a reasonable value
        let n = current_threads();
        assert!(n >= 1);
    }
}
