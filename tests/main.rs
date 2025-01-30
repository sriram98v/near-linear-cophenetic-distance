use itertools::Itertools;
use nlcd::nlcd::near_linear_cophenetic_distance::{LcaMap, NearLinearCopheneticDistance};
use phylo::iter::node_iter::{Ancestors, EulerWalk};
use phylo::prelude::*;
use phylo::tree::distances::{PathFunction, CopheneticDistance};
use phylo::tree::{simple_rtree::RootedTree, DemoTree};
use rand::{distributions::Uniform as Unif, Rng}; // 0.6.5

#[test]
fn cophenetic_dist() {

    fn depth(tree: &DemoTree, node_id: <<DemoTree as RootedTree>::Node as RootedTreeNode>::NodeID)->f32
    {
        tree.depth(node_id) as f32
    }

    let norms = (1..10).collect_vec();


    let mut t1 = DemoTree::yule(100);
    let mut t2 = DemoTree::yule(100);

    t1.precompute_constant_time_lca();
    t2.precompute_constant_time_lca();

    t1.set_zeta(depth);
    t2.set_zeta(depth);

    let t1_lca = LcaMap::from_tree(&t1);
    let t2_lca = LcaMap::from_tree(&t2);    

    for norm in norms.iter(){
        // let pt = DemoTree::pascal_triangle(*norm);
        dbg!((t1.naive_cophen_dist(&t2, &t1_lca, &t2_lca, *norm), t1.nl_cophen_dist(&t2, &t1_lca, &t2_lca, *norm)));
        assert!((t1.cophen_dist(&t2, *norm)-t1.nl_cophen_dist(&t2, &t1_lca, &t2_lca, *norm)).abs()<0.01);
    }
}

#[test]
fn n_choose_k(){
    let n = 100;
    for k in 0..n{
        assert_eq!(DemoTree::n_choose_k(n, k), DemoTree::n_choose_k(n, n-k));
    }
}

#[test]
fn lca_map() {
    let mut t1 = DemoTree::unif(1000);
    t1.precompute_constant_time_lca();
    let lca_map = LcaMap::from_tree(&t1);

    for l1 in t1.get_leaf_ids(){
        for l2 in t1.get_leaf_ids(){
            assert_eq!(*lca_map.get_lca(l1, l2), t1.get_lca_id(vec![l1,l2].as_slice()));
        }
    }

}

#[test]
fn pascal() {
    let norm = 10_u32;

    let pt = DemoTree::pascal_triangle(norm);

    dbg!(&pt);

    for k in 0..norm+1{
        dbg!(k);
        assert_eq!(pt[k as usize], DemoTree::n_choose_k(norm, k) as u32);
    }

}


#[test]
fn seq_product(){
    let mut rng = rand::thread_rng();
    let range = Unif::new(0, 20);
    let vec_len = 10_usize;

    let mut alpha: Vec<u64> = (0..vec_len).map(|_| rng.sample(&range)).collect();
    let mut beta: Vec<u64> = (0..vec_len).map(|_| rng.sample(&range)).collect();
    
    alpha.sort();
    beta.sort();
    
    let alpha: Vec<f32> = alpha.into_iter().map(|x| x as f32).collect();
    let beta: Vec<f32> = beta.into_iter().map(|x| x as f32).collect();

    let mut out: f32 = 0.0;
    for i in 0..vec_len{
        for j in 0..vec_len{
            out += (alpha[i]-beta[j]).abs();
        }
    }

    let pt = DemoTree::pascal_triangle(1);

    let seq_product_out = DemoTree::seq_product(alpha, beta, &pt, 1);

    assert_eq!(seq_product_out, out);
}
