use super::*;
use approx::assert_relative_eq;

fn el(pairs: &[(&str, &str)]) -> Vec<(String, String)> {
    pairs
        .iter()
        .map(|(a, b)| (a.to_string(), b.to_string()))
        .collect()
}

#[test]
fn test_vertex_leaf_and_admixture_edge_order_follow_igraph() {
    // igraph numbers vertices by first appearance; incident_edges lists
    // an admixture node's incoming edges by source vertex, not edge order.
    let g = EdgeListGraph::parse(&el(&[
        ("R", "x"),
        ("R", "a"),
        ("x", "adm"),
        ("a", "adm"),
        ("x", "L1"),
        ("a", "L2"),
        ("adm", "L3"),
        ("R", "y"),
        ("y", "adm2"),
        ("x", "adm2"),
        ("adm2", "L4"),
    ]))
    .unwrap();
    let names: Vec<&str> = g.leaves.iter().map(|&v| g.vertices[v].as_str()).collect();
    assert_eq!(names, ["L1", "L2", "L3", "L4"]);
    // adm: x->adm (edge 2) first; adm2: x->adm2 (edge 9) before y->adm2 (edge 8)
    assert_eq!(g.admix, vec![[2, 3], [9, 8]]);
    assert_eq!(g.normal, vec![0, 1, 4, 5, 6, 7, 10]);
}

#[test]
fn test_invalid_graphs_are_rejected() {
    assert!(EdgeListGraph::parse(&el(&[("R", "A"), ("S", "B")])).is_err());
    assert!(EdgeListGraph::parse(&el(&[("R", "a"), ("a", "b"), ("b", "a"), ("b", "L")])).is_err());
    assert!(EdgeListGraph::parse(&el(&[
        ("R", "a"),
        ("R", "b"),
        ("R", "c"),
        ("a", "m"),
        ("b", "m"),
        ("c", "m"),
        ("m", "L")
    ]))
    .is_err());
}

#[test]
fn test_pwts_split_ancestry_at_admixture() {
    // R -> a -> A; R -> b -> B; a -> m, b -> m (m = w a + (1 - w) b); m -> M
    let g = EdgeListGraph::parse(&el(&[
        ("R", "a"),
        ("R", "b"),
        ("a", "A"),
        ("b", "B"),
        ("a", "m"),
        ("b", "m"),
        ("m", "M"),
    ]))
    .unwrap();
    let p = Paths::new(&g).pwts(&[0.3]);
    // drift edges: R-a, R-b, a-A, b-B, m-M; leaves A, B, M
    let expect = [
        [1.0, 0.0, 0.3],
        [0.0, 1.0, 0.7],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ];
    for (e, row) in expect.iter().enumerate() {
        for (l, &v) in row.iter().enumerate() {
            assert_relative_eq!(p[(e, l)], v, max_relative = 1e-15);
        }
    }
}

#[test]
fn test_box_qp_meets_kkt_conditions() {
    // Oracle: the KKT conditions characterise the unique optimum.
    let n = 6;
    let m = DMatrix::from_fn(n, n, |i, j| ((i * 7 + j * 3) as f64 * 0.9).sin());
    let d = &m * m.transpose() + DMatrix::identity(n, n) * 0.1;
    let q = DVector::from_fn(n, |i, _| ((i as f64) * 1.3).cos() * 2.0);
    let lo = DVector::zeros(n);
    let hi = DVector::from_element(n, 0.8);
    let x = box_qp(&d, &q, &lo, &hi).unwrap();
    let g = &d * &x - &q;
    let mut at_bound = 0;
    for i in 0..n {
        assert!(x[i] >= lo[i] && x[i] <= hi[i]);
        if x[i] == lo[i] {
            at_bound += 1;
            assert!(g[i] >= -1e-12, "at lower bound with g {}", g[i]);
        } else if x[i] == hi[i] {
            at_bound += 1;
            assert!(g[i] <= 1e-12, "at upper bound with g {}", g[i]);
        } else {
            assert!(g[i].abs() < 1e-12, "free with g {}", g[i]);
        }
    }
    assert!(at_bound > 0, "test should exercise the bounds");
    // inside the box it is the unconstrained solution
    let wide = DVector::from_element(n, 1e9);
    let x = box_qp(&d, &q, &(-&wide), &wide).unwrap();
    let exact = d.clone().lu().solve(&q).unwrap();
    assert!((x - exact).amax() < 1e-12);
}

#[test]
fn test_f3_stats_are_the_jackknife_of_block_f3() {
    // Four populations, five blocks of unequal length, f2 block means
    // from an additive tree-like metric plus noise.
    let names = ["P0", "P1", "P2", "P3"];
    let lengths = vec![30, 12, 25, 40, 18];
    let nb = lengths.len();
    let means: Vec<Vec<Vec<f64>>> = (0..4)
        .map(|i| {
            (0..4)
                .map(|j| {
                    (0..nb)
                        .map(|b| {
                            if i == j {
                                0.0
                            } else {
                                let (a, c) = (i.min(j), i.max(j));
                                0.01 * (a + c) as f64
                                    + 0.002 * ((a * 5 + c * 3 + b * 7) as f64).sin()
                            }
                        })
                        .collect()
                })
                .collect()
        })
        .collect();
    let counts: Vec<Vec<Vec<f64>>> = (0..4)
        .map(|_| {
            (0..4)
                .map(|_| lengths.iter().map(|&l| l as f64).collect())
                .collect()
        })
        .collect();
    let pops: Vec<crate::genotypes::PopulationId> = names
        .iter()
        .map(|&n| crate::genotypes::PopulationId::new(n))
        .collect();
    let f2 = F2Blocks::from_block_means(pops, lengths.clone(), &means, &counts).unwrap();
    let st = precompute_f3(&f2, &names, 1e-5).unwrap();

    // f3(P0; Pi, Pj) = (f2(0,i) + f2(0,j) - f2(i,j)) / 2, pairs i <= j
    let pairs = f3_pairs(4);
    assert_eq!(pairs, vec![(1, 1), (1, 2), (1, 3), (2, 2), (2, 3), (3, 3)]);
    let lf: Vec<f64> = lengths.iter().map(|&l| l as f64).collect();
    for (r, &(i, j)) in pairs.iter().enumerate() {
        let blocks: Vec<f64> = (0..nb)
            .map(|b| (means[0][i][b] + means[0][j][b] - means[i][j][b]) / 2.0)
            .collect();
        let (est, se) = crate::fstats::block_jackknife(&blocks, &lf).unwrap();
        assert_relative_eq!(st.est[r], est, max_relative = 1e-12);
        assert_relative_eq!(st.var[(r, r)], se * se, max_relative = 1e-12);
    }
    // f3(P0; Pi, Pi) is f2(P0, Pi)
    let (f2_01, _, _) = f2.f2("P0", "P1").unwrap();
    assert_relative_eq!(st.est[0], f2_01, max_relative = 1e-12);

    // ppinv inverts the covariance with diag_f3 * trace added
    let k = st.var.nrows();
    let mut reg = st.var.clone();
    let add = 1e-5 * reg.trace();
    for i in 0..k {
        reg[(i, i)] += add;
    }
    let id = &st.ppinv * reg;
    assert!((id - DMatrix::identity(k, k)).amax() < 1e-8);
}
