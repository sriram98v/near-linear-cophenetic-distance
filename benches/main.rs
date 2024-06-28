use nlcd::nlcd::near_linear_cophenetic_distance::NearLinearCopheneticDistance;
use phylo::prelude::*;
use phylo::tree::SimpleRootedTree;

const NUM_TAXA: usize = 4000;
const NORM: u32 = 2;

fn main() {
    // Run registered benchmarks.
    divan::main();
}

#[divan::bench]
fn benchmark_nlcd(bencher: divan::Bencher) {
    bencher
        .with_inputs(|| {
            fn depth(
                tree: &SimpleRootedTree,
                node_id: <SimpleRootedTree as RootedTree>::NodeID,
            ) -> f32 {
                EulerWalk::get_node_depth(tree, node_id) as f32
            }

            let mut t1 = SimpleRootedTree::yule(NUM_TAXA).unwrap();
            let mut t2 = SimpleRootedTree::yule(NUM_TAXA).unwrap();
            t1.precompute_constant_time_lca();
            t2.precompute_constant_time_lca();
            t1.set_zeta(depth);
            t2.set_zeta(depth);
            (t1, t2)
        })
        .bench_refs(|(t1, t2)| {
            t1.cophen_dist(t2, NORM);
        });
}

#[divan::bench]
fn benchmark_naive(bencher: divan::Bencher) {
    bencher
        .with_inputs(|| {
            fn depth(
                tree: &SimpleRootedTree,
                node_id: <SimpleRootedTree as RootedTree>::NodeID,
            ) -> f32 {
                EulerWalk::get_node_depth(tree, node_id) as f32
            }

            let mut t1 = SimpleRootedTree::yule(NUM_TAXA).unwrap();
            let mut t2 = SimpleRootedTree::yule(NUM_TAXA).unwrap();
            t1.precompute_constant_time_lca();
            t2.precompute_constant_time_lca();
            t1.set_zeta(depth);
            t2.set_zeta(depth);
            (t1, t2)
        })
        .bench_refs(|(t1, t2)| {
            t1.cophen_dist_naive(t2, NORM);
        });
}
