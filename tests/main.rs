use nlcd::nlcd::near_linear_cophenetic_distance::NearLinearCopheneticDistance;
use phylo::iter::node_iter::{Ancestors, EulerWalk};
use phylo::prelude::*;
use phylo::tree::distances::{PathFunction, CopheneticDistance};
use phylo::tree::{simple_rtree::RootedTree, SimpleRootedTree};

#[test]
fn cophenetic_dist() {
    let norm = 2;

    fn depth(tree: &SimpleRootedTree, node_id: <<SimpleRootedTree as RootedTree>::Node as RootedTreeNode>::NodeID)->f32
    {
        tree.depth(node_id) as f32
    }

    let mut t1 = SimpleRootedTree::yule(100).unwrap();
    let mut t2 = SimpleRootedTree::yule(100).unwrap();

    t1.precompute_constant_time_lca();
    t2.precompute_constant_time_lca();

    t1.set_zeta(depth);
    t2.set_zeta(depth);

    assert_eq!(t1.cophen_dist_naive(&t2, norm), t1.cophen_dist(&t2, norm));
}
