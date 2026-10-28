//! Numbers as C's `printf("%.*g")` and Python's `f"{x:.12g}"` write them,
//! for tables written in Rust that must be the same text as those Python
//! writes (TASKS.md F14).

/// `x` with `precision` significant digits, as `%.{precision}g`: fixed
/// notation when the exponent is from -4 to precision - 1, else
/// `d.ddde±XX`; trailing zeros and a trailing point removed; `inf`, `-inf`,
/// `nan`.
pub fn format_g(x: f64, precision: usize) -> String {
    if x.is_nan() {
        return "nan".to_string();
    }
    if x.is_infinite() {
        return if x > 0.0 { "inf" } else { "-inf" }.to_string();
    }
    let p = precision.max(1);
    if x == 0.0 {
        return if x.is_sign_negative() { "-0" } else { "0" }.to_string();
    }
    // the exponent after rounding to p digits
    let sci = format!("{:.*e}", p - 1, x);
    let (mantissa, exp) = sci.split_once('e').expect("exponent");
    let exp: i32 = exp.parse().expect("integer exponent");
    if exp < -4 || exp >= p as i32 {
        let m = strip(mantissa);
        let sign = if exp < 0 { '-' } else { '+' };
        format!("{m}e{sign}{:02}", exp.abs())
    } else {
        let decimals = (p as i32 - 1 - exp).max(0) as usize;
        strip(&format!("{x:.decimals$}")).to_string()
    }
}

fn strip(s: &str) -> &str {
    if s.contains('.') {
        s.trim_end_matches('0').trim_end_matches('.')
    } else {
        s
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn as_printf_writes_them() {
        for (x, want) in [
            (0.0, "0"),
            (1.0, "1"),
            (-0.25, "-0.25"),
            (0.1 + 0.2, "0.3"),
            (123456789012.0, "123456789012"),
            (1234567890123.0, "1.23456789012e+12"),
            (0.0001, "0.0001"),
            (0.00001234, "1.234e-05"),
            (-1.5e-300, "-1.5e-300"),
            (99999999999.99, "100000000000"),
            (999999999999.5, "1e+12"),
            (f64::NEG_INFINITY, "-inf"),
        ] {
            assert_eq!(format_g(x, 12), want, "{x}");
        }
    }
}
