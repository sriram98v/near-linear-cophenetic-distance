use std::{fmt::{Debug, Display}, hash::Hash};
use fxhash::FxHashSet as HashSet;
use fxhash::FxHashMap as HashMap;

use num::{Float, NumCast, Signed, Zero};
use phylo::prelude::*;
use itertools::Itertools;

// pub struct TaxaIntersections<T>
// where
//     T: RootedMetaTree,
// {
//     a_int_a_hat: Vec<<T as NearLinearCopheneticDistance>::Meta>,
//     a_int_b_hat: Vec<<T as NearLinearCopheneticDistance>::Meta>,
//     b_int_a_hat: Vec<<T as NearLinearCopheneticDistance>::Meta>,
//     b_int_b_hat: Vec<<T as NearLinearCopheneticDistance>::Meta>,
// }


pub trait NearLinearCopheneticDistance: CopheneticDistance
where
    <Self as RootedTree>::Node: RootedMetaNode<Meta=<Self as CopheneticDistance>::Meta> + RootedZetaNode,
    <<Self as RootedTree>::Node as RootedZetaNode>::Zeta: Signed + Clone + NumCast + std::iter::Sum + Debug + Display + Float + PartialOrd + Copy + Send,
{
    type Meta: Display + Debug + Eq + PartialEq + Clone + Ord + Hash + Send + Sync;

    // helper functions

    /// Returns taxa present in upper tree.
    fn upper_tree_taxa(&self)->impl Iterator<Item=<Self as CopheneticDistance>::Meta>
    {
        let lower_tree_taxa = self.lower_tree_taxa().collect::<HashSet<_>>();
        self.get_leaves().map(|x| x.get_taxa().unwrap()).filter(move |x| !lower_tree_taxa.contains(x))
    }

    /// Returns taxa present in lower tree.
    fn lower_tree_taxa(&self)->impl Iterator<Item=<Self as CopheneticDistance>::Meta>
    {
        let median_node = self.get_median_node_id();
        self.get_cluster(median_node.clone()).into_iter().filter(|x| x.get_taxa().is_some()).map(|x| x.get_taxa().unwrap())
    }

    /// Returns taxa present in upper tree.
    fn upper_tree_leaves(&self)->impl Iterator<Item=<Self as RootedTree>::Node>
    {
        let lower_tree_leaf_ids = self.lower_tree_leaves().map(|x| x.get_id()).collect::<HashSet<_>>();
        self.get_leaves().filter(move |x| !lower_tree_leaf_ids.contains(&x.get_id())).map(|x| x.clone())
    }

    /// Returns leaves present in lower tree.
    fn lower_tree_leaves(&self)->impl Iterator<Item=<Self as RootedTree>::Node>
    {
        let median_node = self.get_median_node_id();
        self.get_cluster(median_node.clone())
    }
    
    /// Returns lower tree.
    fn lower_tree(&self)->Self
    {
        let lower_tree_taxa = HashSet::from_iter(self.lower_tree_taxa());
        self.contract_tree(&lower_tree_taxa.iter().map(|x| self.get_taxa_node_id(x).unwrap()).collect_vec())

    }

    /// Returns upper tree.
    fn upper_tree(&self)->Self
    {
        let upper_tree_taxa = HashSet::from_iter(self.upper_tree_taxa());
        self.contract_tree(&upper_tree_taxa.iter().map(|x| self.get_taxa_node_id(x).unwrap()).collect_vec())
    }

    fn get_node_sibling_id(&self, node_id: &<Self as RootedTree>::NodeID)-> <Self as RootedTree>::NodeID
    {
        let node_parent_id = self.get_node_parent_id(node_id.clone()).expect("Node has no siblings");
        let siblings = self.get_node_children_ids(node_parent_id).collect_vec();
        if &siblings[0]==node_id{
            return siblings[1]
        }
        else{
            return siblings[0];
        }
    }

    /// Returns value of n choose k.
    fn n_choose_k(n: u32, k: u32)-> i32
    {
        fn multiplicative_form(n: u32, k: u32)->i32
        {
            (1..k+1).map(|x| (n+1-x) as f32/(x as f32)).product::<f32>() as i32
        }
    
        match k <= (n as f32/2 as f32) as u32{
            true =>{multiplicative_form(n, k)},
            false =>{multiplicative_form(n-k, n-k)},
        }
    }

    /// Returns the Cophenetic distance between two trees in O(pnlog^2n) time. 
    fn cophen_dist(&self, tree: &Self, norm: u32)-><<Self as RootedTree>::Node as RootedZetaNode>::Zeta
    {
        let mut ops: Vec<<<Self as RootedTree>::Node as RootedZetaNode>::Zeta> = vec![];

        let binding1 = self.get_taxa_space().into_iter().collect::<HashSet<<Self as CopheneticDistance>::Meta>>();
        let binding2 = tree.get_taxa_space().into_iter().collect::<HashSet<<Self as CopheneticDistance>::Meta>>();
        let taxa_set = binding1.intersection(&binding2).map(|x| x.clone()).collect_vec();

        self.populate_op_vec(tree, norm, taxa_set.clone(), &mut ops);

        let distance = ops.into_iter().sum::<<<Self as RootedTree>::Node as RootedZetaNode>::Zeta>();
    
        let taxa_distance = taxa_set.iter()
                            .map(|x| {
                                let zeta_1 = self.get_zeta_taxa(x);
                                let zeta_2 = tree.get_zeta_taxa(x);
                                return (zeta_1-zeta_2).abs().powi(norm as i32)
                            }).sum();
        return (distance + taxa_distance).powf(<<<Self as RootedTree>::Node as RootedZetaNode>::Zeta as NumCast>::from(norm).unwrap().powi(-1));

    }

    /// Populates vector with divide and conquer distances.
    fn populate_op_vec(&self, tree: &Self, norm: u32, taxa_set: Vec<<Self as RootedMetaTree>::Meta>, op_vec: &mut Vec<<<Self as RootedTree>::Node as RootedZetaNode>::Zeta>)
    {
        let t = self.get_median_node_id();
        let t_hat = tree.get_median_node_id();
        
        let b: HashSet<<Self as CopheneticDistance>::Meta> = HashSet::from_iter(self.get_cluster(t).filter(|x| x.get_taxa().is_some()).map(|x| x.get_taxa().unwrap()).filter(|x| taxa_set.contains(x)));
        let b_hat: HashSet<<Self as CopheneticDistance>::Meta> = HashSet::from_iter(tree.get_cluster(t_hat).filter(|x| x.get_taxa().is_some()).map(|x| x.get_taxa().unwrap()).filter(|x| taxa_set.contains(x)));

        let a: HashSet<<Self as CopheneticDistance>::Meta> = HashSet::from_iter(self.get_taxa_space()).difference(&b).filter(|x| taxa_set.contains(x)).map(|x| x.clone()).collect();
        let a_hat: HashSet<<Self as CopheneticDistance>::Meta> = HashSet::from_iter(tree.get_taxa_space()).difference(&b_hat).filter(|x| taxa_set.contains(x)).map(|x| x.clone()).collect();

        let a_int_a_hat = a.intersection(&a_hat).map(|x| x.clone()).collect_vec();
        let a_int_b_hat = a.intersection(&b_hat).map(|x| x.clone()).collect_vec();
        let b_int_a_hat = b.intersection(&a_hat).map(|x| x.clone()).collect_vec();
        let b_int_b_hat = b.intersection(&b_hat).map(|x| x.clone()).collect_vec();

        let double_mix_distance = self.distance_double_mix_type(tree, norm, &t, &t_hat, &a_int_a_hat, &a_int_b_hat, &b_int_a_hat, &b_int_b_hat);
        let single_mix_distance = self.distance_single_mix_type(tree, norm, &t, &t_hat, &a_int_a_hat, &a_int_b_hat, &b_int_a_hat, &b_int_b_hat);

        op_vec.push(double_mix_distance);
        op_vec.push(single_mix_distance);

        if taxa_set.len()>2{        
            if a_int_a_hat.len()>1{
                let self_tree = self.contract_tree(&a_int_a_hat.iter().map(|x| self.get_taxa_node_id(x).unwrap()).collect_vec());
                let new_tree = tree.contract_tree(&a_int_a_hat.iter().map(|x| tree.get_taxa_node_id(x).unwrap()).collect_vec());
                self_tree.populate_op_vec(&new_tree, norm, a_int_a_hat, op_vec);
            }
    
            if a_int_b_hat.len()>1{
                let self_tree = self.contract_tree(&a_int_b_hat.iter().map(|x| self.get_taxa_node_id(x).unwrap()).collect_vec());
                let new_tree = tree.contract_tree(&a_int_b_hat.iter().map(|x| tree.get_taxa_node_id(x).unwrap()).collect_vec());
                self_tree.populate_op_vec(&new_tree, norm, a_int_b_hat, op_vec);
            }
    
            if b_int_b_hat.len()>1{
                let self_tree = self.contract_tree(&b_int_b_hat.iter().map(|x| self.get_taxa_node_id(x).unwrap()).collect_vec());
                let new_tree = tree.contract_tree(&b_int_b_hat.iter().map(|x| tree.get_taxa_node_id(x).unwrap()).collect_vec());
                self_tree.populate_op_vec(&new_tree, norm, b_int_b_hat, op_vec);
            }
    
            if b_int_a_hat.len()>1{
                let self_tree = self.contract_tree(&b_int_a_hat.iter().map(|x| self.get_taxa_node_id(x).unwrap()).collect_vec());
                let new_tree = tree.contract_tree(&b_int_a_hat.iter().map(|x| tree.get_taxa_node_id(x).unwrap()).collect_vec());
                self_tree.populate_op_vec(&new_tree, norm, b_int_a_hat, op_vec);
            }
        }
    }

    /// Returns ordered iterator used in double mix type cases
    fn get_cntr(&self, leaf_set: HashSet<Self::NodeID>)->Vec<<<Self as RootedTree>::Node as RootedZetaNode>::Zeta>
    {
        // line 5 in algo 1
        let mut gamma: Vec<<<Self as RootedTree>::Node as RootedZetaNode>::Zeta> = Vec::new();
        let median_node_id = self.get_median_node_id();
        // line 3 in algo 1
        let mut median_path = self.root_to_node(median_node_id.clone()).into_iter().map(|x| (x.get_id(), 0)).collect::<HashMap<_,_>>();
        for node_id in leaf_set{
            // line 4 in algo 1
            median_path.entry(self.get_lca_id(&vec![node_id, median_node_id.clone()])).and_modify(|x| *x+=1);
        }
        for node in self.root_to_node(median_node_id.clone()).into_iter(){
            let c = median_path.get(&node.get_id()).cloned().unwrap();
            for _ in 0..c{
                gamma.push(self.get_zeta(node.get_id()).unwrap())
            }
        }
        gamma
    }

    /// Returns seqPrd^p(\alpha,\beta) used for double mix type cases.
    fn seq_product(mut alpha: Vec<<<Self as RootedTree>::Node as RootedZetaNode>::Zeta>, mut beta:Vec<<<Self as RootedTree>::Node as RootedZetaNode>::Zeta>, norm: u32)-><<Self as RootedTree>::Node as RootedZetaNode>::Zeta
    {
        if alpha.is_empty() || beta.is_empty(){
            return <<<Self as RootedTree>::Node as RootedZetaNode>::Zeta as Zero>::zero();
        }
        if alpha.last().unwrap()<&alpha[0]{
            alpha.reverse();
        }
        if beta.last().unwrap()<&beta[0]{
            beta.reverse();
        }
        if alpha.last().unwrap()>beta.last().unwrap(){
            std::mem::swap(&mut alpha, &mut beta);
        }
        let sigma = Self::precompute_alpha_sums(&alpha, &beta, norm);
        let final_out = (0..beta.len()).map(|j| {
            let term_1 = (0..norm+1).map(|l| {
                let t1 = <<<Self as RootedTree>::Node as RootedZetaNode>::Zeta as NumCast>::from(Self::n_choose_k(norm, l)).unwrap();
                let t2 = beta[j].powi(l as i32);
                let t3 = <<<Self as RootedTree>::Node as RootedZetaNode>::Zeta as NumCast>::from((-1 as i32).pow(norm-l)).unwrap();
                let t4 = sigma[j][(norm-l) as usize];
                return t1*t2*t3*t4;
            }).sum::<<<Self as RootedTree>::Node as RootedZetaNode>::Zeta>();
            
            let term_2 = (0..norm+1).map(|l| {
                let t1 = <<<Self as RootedTree>::Node as RootedZetaNode>::Zeta as NumCast>::from(Self::n_choose_k(norm, l)).unwrap();
                let t2 = beta[j].powi((norm-l) as i32);
                let t3 = <<<Self as RootedTree>::Node as RootedZetaNode>::Zeta as NumCast>::from((-1 as i32).pow(norm-l)).unwrap();
                let t4 = alpha.iter().map(|a| a.powi(l as i32)).sum::<<<Self as RootedTree>::Node as RootedZetaNode>::Zeta>()-sigma[j][l as usize];
                return t1*t2*t3*t4;
            }).sum::<<<Self as RootedTree>::Node as RootedZetaNode>::Zeta>();

            return term_1+term_2;
        }).sum::<<<Self as RootedTree>::Node as RootedZetaNode>::Zeta>();

        final_out
    }

    /// Algorithm 3 for precomputation of alpha sums
    fn precompute_alpha_sums(alpha: &Vec<<<Self as RootedTree>::Node as RootedZetaNode>::Zeta>, beta: &Vec<<<Self as RootedTree>::Node as RootedZetaNode>::Zeta>, norm: u32) -> Vec<Vec<<<Self as RootedTree>::Node as RootedZetaNode>::Zeta>>
    {
        let mut j = 0;
        let mut i = 0;
        let k = alpha.len();
        let m = beta.len();

        let mut sigma = vec![vec![<<Self as RootedTree>::Node as RootedZetaNode>::Zeta::zero(); (norm+1) as usize]; m];

        while j < m {
            match i < k && alpha[i]<=beta[j]{
                true => {
                    (0..norm+1).for_each(|l| sigma[j][l as usize] = sigma[j][l as usize] + alpha[i].powi(l as i32));
                    i += 1;
                },
                false => {
                    j += 1;
                    if j < m{
                        (0..norm+1).for_each(|l| sigma[j][l as usize] = sigma[j-1][l as usize]);
                    }
                },
            }
        }

        sigma
    }

    /// This method generates the distance contributed by all taxa pairs 
    /// that are present in different subtrees in both trees(raised to the p^{th} power).
    /// 
    /// This includes the following assignments: AB|A'B', AB|B'A'
    fn distance_double_mix_type(&self, 
        tree: &Self, 
        norm: u32, 
        t: &<Self as RootedTree>::NodeID, 
        t_hat: &<Self as RootedTree>::NodeID, 
        a_int_a_hat: &Vec<<Self as CopheneticDistance>::Meta>,
        a_int_b_hat: &Vec<<Self as CopheneticDistance>::Meta>,
        b_int_a_hat: &Vec<<Self as CopheneticDistance>::Meta>,
        b_int_b_hat: &Vec<<Self as CopheneticDistance>::Meta>)-><<Self as RootedTree>::Node as RootedZetaNode>::Zeta
    {
        let alpha = self.get_cntr(a_int_b_hat.iter().map(|x| self.get_taxa_node_id(&x).unwrap()).collect::<HashSet<Self::NodeID>>());
        let beta = tree.get_cntr(b_int_a_hat.iter().map(|x| tree.get_taxa_node_id(&x).unwrap()).collect::<HashSet<Self::NodeID>>());

        // AB|A'B'
        let b_int_b_hat_len = b_int_b_hat.len();
        let dd2 = a_int_a_hat.iter().map(|x| x.clone())
            .map(|x| {
                let t_lca_id = self.get_lca_id(&vec![self.get_taxa_node_id(&x).unwrap(), t.clone()]);
                let t_hat_lca_id = tree.get_lca_id(&vec![tree.get_taxa_node_id(&x).unwrap(), t_hat.clone()]);
                let zeta_1 = self.get_zeta(t_lca_id).unwrap();
                let zeta_2 = tree.get_zeta(t_hat_lca_id).unwrap();
                return (zeta_1-zeta_2).abs().powi(norm as i32)
            })
            .sum::<<<Self as RootedTree>::Node as RootedZetaNode>::Zeta>() * <<<Self as RootedTree>::Node as RootedZetaNode>::Zeta as NumCast>::from(b_int_b_hat_len).unwrap();

        return Self::seq_product(alpha, beta, norm) + dd2;
    }

    /// This method generates the distance contributed by all taxa pairs 
    /// that are present in the same subtree in exactly one of the two trees(raised to the p^{th} power).
    /// 
    /// This includes the following assignments: AA|A'B', AA|B'A', BB|A'B', BB|B'A', BA|B'B', BA|A'A', AB|B'B', AB|A'A'.
    fn distance_single_mix_type(&self,
        tree: &Self,
        norm: u32,
        t: &<Self as RootedTree>::NodeID, 
        t_hat: &<Self as RootedTree>::NodeID, 
        a_int_a_hat: &Vec<<Self as CopheneticDistance>::Meta>,
        a_int_b_hat: &Vec<<Self as CopheneticDistance>::Meta>,
        b_int_a_hat: &Vec<<Self as CopheneticDistance>::Meta>,
        b_int_b_hat: &Vec<<Self as CopheneticDistance>::Meta>)-><<Self as RootedTree>::Node as RootedZetaNode>::Zeta
    {
        if self.num_taxa()<=2{
            return <<<Self as RootedTree>::Node as RootedZetaNode>::Zeta>::zero();
        }
        // AA|A'B'
        let d1 = Self::single_mix_xxxy(self, tree, t, t_hat, norm, a_int_a_hat, true);
        // AB|A'A'
        let d2 = Self::single_mix_xxxy(tree, self, t_hat, t, norm, a_int_a_hat, true);
        // BB|A'B'
        let d3 = Self::single_mix_xxxy(self, tree, t, t_hat, norm, b_int_a_hat, false);
        // AB|B'B'
        let d4 = Self::single_mix_xxxy(tree, self, t_hat, t, norm, a_int_b_hat, false);

        return d1+d2+d3+d4;
    }

    fn preprocess_single_mix_even(t1: &Self,
        t2: &Self, 
        t2_median: &<Self as RootedTree>::NodeID,
        norm: u32,
        taxa_set: &Vec<<Self as CopheneticDistance>::Meta>,
        sigma: &mut HashMap<<Self as RootedTree>::NodeID,Vec<<<Self as RootedTree>::Node as RootedZetaNode>::Zeta>>,
        counterpart_count: &mut HashMap<<Self as RootedTree>::NodeID,u32>,
        upper_mixed: bool)
        {
            let leaf_iter = match upper_mixed{
                true => t1.upper_tree_taxa().collect_vec(),
                false => t1.lower_tree_taxa().collect_vec(),
            };
            for leaf in leaf_iter{
                    match taxa_set.contains(&leaf){
                        true => {
                            let t1_node_id = t1.get_taxa_node_id(&leaf).unwrap();
                            let t2_node_id = t2.get_taxa_node_id(&leaf).unwrap();
                            let lca_x_t_hat = t2.get_lca_id(&vec![t2_node_id.clone(), t2_median.clone()]);
                            let beta = t2.get_zeta(lca_x_t_hat).unwrap();
                            sigma.insert(t1_node_id, (0..norm+1).map(|l| beta.powi(l as i32)).collect_vec());
                        },
                        false => {
                            let node_id = t1.get_taxa_node_id(&leaf).unwrap();
                            counterpart_count.entry(node_id).and_modify(|c| *c += 1);
                        },    
                }
            }
        }


    fn preprocess_single_mix_odd(t1: &Self,
        t2: &Self, 
        t1_median: &<Self as RootedTree>::NodeID,
        t2_median: &<Self as RootedTree>::NodeID,
        norm: u32,
        taxa_set: &Vec<<Self as CopheneticDistance>::Meta>,
        sigma_pos: &mut HashMap<<Self as RootedTree>::NodeID,Vec<<<Self as RootedTree>::Node as RootedZetaNode>::Zeta>>,
        sigma_neg: &mut HashMap<<Self as RootedTree>::NodeID,Vec<<<Self as RootedTree>::Node as RootedZetaNode>::Zeta>>,
        delta: &mut HashMap<<Self as RootedTree>::NodeID,Vec<<<Self as RootedTree>::Node as RootedZetaNode>::Zeta>>,
        counterpart_count: &mut HashMap<<Self as RootedTree>::NodeID,u32>,
        upper_mixed: bool)
        {
            let leaf_iter = match upper_mixed{
                true => t1.upper_tree_taxa().collect_vec(),
                false => t1.lower_tree_taxa().collect_vec(),
            };
            for leaf in leaf_iter{
                    match taxa_set.contains(&leaf){
                        true => {
                            let t2_node_id = t2.get_taxa_node_id(&leaf).unwrap();
                            let lca_x_t_hat = t2.get_lca_id(&vec![t2_node_id.clone(), t2_median.clone()]);
                            let beta = t2.get_zeta(lca_x_t_hat).unwrap();
                            let t1_node_id = t1.get_taxa_node_id(&leaf).unwrap();

                            let t1_node_parent_id = t1.get_node_parent_id(t1_node_id).unwrap();

                            if beta <= t1.get_zeta(t1_node_id).unwrap(){
                                //find omega_x
                                let mut v_ancestors = vec![];
                                for node in t1.node_to_root(t1_node_id).into_iter(){
                                    v_ancestors.push(node.get_id());
                                    if &node.get_id()==t1_median{
                                        break;
                                    }
                                }
                                    
                                let omega_x = v_ancestors.into_iter().min_by(|w,y| {
                                    t1.get_zeta(w.clone()).unwrap().partial_cmp(&t1.get_zeta(y.clone()).unwrap()).unwrap()
                                }).unwrap();
                                // set omega_x.delta
                                delta.entry(omega_x)
                                    .and_modify(|e| {
                                        for l in 0..norm+1{
                                            e[l as usize] = e[l as usize]+beta.powi(l as i32);
                                        }
                                    });
                            }
                            if beta <= t1.get_zeta(t1_node_parent_id).unwrap(){
                                sigma_pos.insert(t1_node_id, (0..norm+1).map(|l| beta.powi(l as i32)).collect_vec());
                            }
                            else{
                                sigma_neg.insert(t1_node_id, (0..norm+1).map(|l| beta.powi(l as i32)).collect_vec());
                            }
                        },
                        false => {
                            let node_id = t1.get_taxa_node_id(&leaf).unwrap();
                            counterpart_count.entry(node_id).and_modify(|c| *c += 1);
                        },    
                }
            }
        }


    // this method solves AA|A'B'
    fn single_mix_xxxy(t1: &Self, 
        t2: &Self,
        t1_median: &<Self as RootedTree>::NodeID,
        t2_median: &<Self as RootedTree>::NodeID, 
        norm: u32, 
        taxa_set: &Vec<<Self as CopheneticDistance>::Meta>,
        upper_mixed: bool)-><<Self as RootedTree>::Node as RootedZetaNode>::Zeta
    {
        let subtree_nodes = match upper_mixed{
            true => {        
                let lower_tree_nodes = t1.postord(t1_median.clone()).map(|x| x.get_id()).collect_vec();
                t1.postord(t1.get_root_id()).map(|x| x.get_id()).filter(|x| !lower_tree_nodes.contains(x)).collect_vec()
            },
            false => {
                t1.postord(t1_median.clone()).map(|x| x.get_id()).collect_vec()
            },
        };

        let mut kappa: HashMap<_,_> = t1.get_node_ids().map(|x| (x, <<Self as RootedTree>::Node as RootedZetaNode>::Zeta::zero())).collect();
        let mut counterpart_count: HashMap<_,_> = t1.get_node_ids().map(|x| (x, 0_u32)).collect();

        match norm%2{
            0 => {
                // setting sigma to zeros; kappa already set.
                let mut sigma: HashMap<_,_> = t1.get_node_ids().map(|x| (x, vec![<<Self as RootedTree>::Node as RootedZetaNode>::Zeta::zero();norm as usize+1])).collect();        
                // Preprocessing loop
                Self::preprocess_single_mix_even(t1, t2, &t2_median, norm, taxa_set, &mut sigma, &mut counterpart_count, upper_mixed);

                for v_id in subtree_nodes.iter(){
                    if v_id!=&t1.get_root_id() && !t1.is_leaf(&v_id){
                        // calculate v_sigma
                        let v_value = t1.get_node_children_ids(v_id.clone())
                            .map(|x| sigma.get(&x).cloned().unwrap())
                            .fold(vec![<<Self as RootedTree>::Node as RootedZetaNode>::Zeta::zero();norm as usize +1], |acc, x| {
                                acc.iter().zip(x).map(|(a,b)| *a+b).collect_vec()
                            });
                        sigma.insert(v_id.clone(), v_value);
                        let v_counterpart_count = t1.get_node_children_ids(v_id.clone()).map(|chid| counterpart_count.get(&chid).unwrap()).sum();
                        counterpart_count.insert(v_id.clone(), v_counterpart_count);
                    }
                }


                for v_id in subtree_nodes.iter(){
                    if v_id!=&t1.get_root_id(){
                        kappa.entry(v_id.clone()).and_modify(|v_kappa| {
                            let v_sibling_id = t1.get_node_sibling_id(&v_id);
                            let v_parent_id = t1.get_node_parent_id(v_id.clone()).unwrap();
                            let v_counterpart_count = <<<Self as RootedTree>::Node as RootedZetaNode>::Zeta as NumCast>::from(counterpart_count.get(&v_sibling_id).cloned().unwrap()).unwrap();
                            let summation_term = (0..norm+1).map(|l| {
                                let term1 = <<<Self as RootedTree>::Node as RootedZetaNode>::Zeta as NumCast>::from(Self::n_choose_k(norm, l)).unwrap();
                                let term2 = <<<Self as RootedTree>::Node as RootedZetaNode>::Zeta as NumCast>::from((-1_i32).pow(norm-l)).unwrap();
                                let term3 = (t1.get_zeta(v_parent_id).unwrap()).powi((norm-l) as i32);
                                let term4 = sigma.get(&v_id).unwrap()[l as usize];
                                return term1*term2*term3*term4;
                            }).sum::<<<Self as RootedTree>::Node as RootedZetaNode>::Zeta>();
                            *v_kappa = v_counterpart_count*summation_term;
                        });
                    }
                }
            }
            _ => {
                let mut sigma_pos: HashMap<_,_> = t1.get_node_ids().map(|x| (x, vec![<<Self as RootedTree>::Node as RootedZetaNode>::Zeta::zero();norm as usize+1])).collect();
                let mut sigma_neg: HashMap<_,_> = t1.get_node_ids().map(|x| (x, vec![<<Self as RootedTree>::Node as RootedZetaNode>::Zeta::zero();norm as usize+1])).collect();        
                let mut delta: HashMap<_,_> = t1.get_node_ids().map(|x| (x, vec![<<Self as RootedTree>::Node as RootedZetaNode>::Zeta::zero(); norm as usize+1])).collect();
                
                Self::preprocess_single_mix_odd(t1, t2, &t1_median, &t2_median, norm, taxa_set, &mut sigma_pos, &mut sigma_neg, &mut delta, &mut counterpart_count, upper_mixed);

                for v_id in subtree_nodes.iter(){
                    if v_id!=&t1.get_root_id() && !t1.is_leaf(&v_id){
                        let v_children = t1.get_node_children_ids(v_id.clone()).collect_vec();
                        let v_left_child_id = v_children[0];
                        let v_right_child_id = v_children[1];
                        let v_delta = delta.get(&v_id).cloned().unwrap();
                        // calculate v_sigma_pos
                        let v_left_sigma_pos = sigma_pos.get(&v_left_child_id).cloned().unwrap();
                        let v_right_sigma_pos = sigma_pos.get(&v_right_child_id).cloned().unwrap();
                        sigma_pos.entry(v_id.clone())
                            .and_modify(|e| {
                                for l in 0..norm+1{
                                    e[l as usize] = v_left_sigma_pos[l as usize]+v_right_sigma_pos[l as usize]-v_delta[l as usize];
                                }
                            });
                        // calculate v_sigma_neg
                        let v_left_sigma_neg = sigma_neg.get(&v_left_child_id).cloned().unwrap();
                        let v_right_sigma_neg = sigma_neg.get(&v_right_child_id).cloned().unwrap();
                        sigma_neg.entry(v_id.clone())
                            .and_modify(|e| {
                                for l in 0..norm+1{
                                    e[l as usize] = v_left_sigma_neg[l as usize]+v_right_sigma_neg[l as usize]+v_delta[l as usize];
                                }
                            }
                        );
                        let v_counterpart_count = t1.get_node_children_ids(v_id.clone()).map(|chid| counterpart_count.get(&chid).unwrap()).sum();
                        counterpart_count.insert(v_id.clone(), v_counterpart_count);    
                    }
                }

                for v_id in subtree_nodes{
                    if v_id!=t1.get_root_id(){
                        kappa.entry(v_id).and_modify(|v_kappa| {
                            let v_sibling_id = t1.get_node_sibling_id(&v_id);
                            let v_parent_id = t1.get_node_parent_id(v_id).unwrap();
                            let v_parent_zeta = t1.get_zeta(v_parent_id).unwrap();
                            let v_counterpart_count = <<<Self as RootedTree>::Node as RootedZetaNode>::Zeta as NumCast>::from(counterpart_count.get(&v_sibling_id).cloned().unwrap()).unwrap();
                            let summation_term = (0..norm+1).map(|l| {
                                let term1 = <<<Self as RootedTree>::Node as RootedZetaNode>::Zeta as NumCast>::from(Self::n_choose_k(norm, l)).unwrap();
                                let term2 = <<<Self as RootedTree>::Node as RootedZetaNode>::Zeta as NumCast>::from((-1_i32).pow(norm-l)).unwrap();
                                
                                let term3_1 = v_parent_zeta.powi(l as i32)*sigma_pos.get(&v_id).unwrap()[(norm-l) as usize];
                                let term3_2 = v_parent_zeta.powi((norm-l) as i32)*sigma_neg.get(&v_id).unwrap()[l as usize];
                                
                                let term3 = term3_1+term3_2;

                                return term1*term2*term3;
                            }).sum::<<<Self as RootedTree>::Node as RootedZetaNode>::Zeta>();
                            *v_kappa = v_counterpart_count*summation_term;
                        });
                    }
                }
            }
        }
        return kappa.into_values().sum();
    }
}