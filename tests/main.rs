use nlcd::nlcd::near_linear_cophenetic_distance::NearLinearCopheneticDistance;
use phylo::iter::node_iter::{Ancestors, EulerWalk};
use phylo::tree::distances::PathFunction;
use phylo::tree::io::Newick;
use phylo::tree::ops::CopheneticDistance;
use phylo::tree::{simple_rtree::RootedTree, SimpleRootedTree};

#[test]
fn cophenetic_dist() {
    let norm = 2;

    fn depth(tree: &SimpleRootedTree, node_id: <SimpleRootedTree as RootedTree>::NodeID) -> f32 {
        tree.depth(node_id) as f32
    }
    // let t1_input_str: String = String::from("(B,(C,(A,D)));");
    // let t2_input_str: String = String::from("(A,(C,(B,D)));");
    // let mut t1 = SimpleRootedTree::from_newick(t1_input_str.as_bytes());
    // let mut t2 = SimpleRootedTree::from_newick(t2_input_str.as_bytes());

    let mut t1 = SimpleRootedTree::yule(100).unwrap();
    let mut t2 = SimpleRootedTree::yule(100).unwrap();
    dbg!(format!("{}", t1.to_newick()));
    dbg!(format!("{}", t2.to_newick()));

    t1.precompute_constant_time_lca();
    t2.precompute_constant_time_lca();

    t1.set_zeta(depth);
    t2.set_zeta(depth);

    assert_eq!(t1.cophen_dist_naive(&t2, norm), t1.cophen_dist(&t2, norm));
}
