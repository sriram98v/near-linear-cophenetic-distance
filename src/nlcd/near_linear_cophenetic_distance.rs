use fxhash::FxHashMap as HashMap;
use fxhash::FxHashSet as HashSet;

use std::{
    fmt::{Debug, Display},
    ops::{AddAssign, IndexMut, Index},
};

use vers_vecs::BitVec;
use itertools::Itertools;
use num::{Float, NumCast, Signed, Zero};
use phylo::prelude::*;

#[derive(Debug, Clone)]
pub struct LcaMap<T: RootedTree + EulerWalk>{
    data: HashMap<(TreeNodeID<T>, TreeNodeID<T>), TreeNodeID<T>>,
}

impl<T: RootedTree + EulerWalk> LcaMap<T>{
    pub fn from_tree(tree: &T)->Self{
        // let leaf_max_id = tree.get_leaf_ids().max().unwrap()+1;
        let mut data = vec![].into_iter().collect::<HashMap<(TreeNodeID<T>, TreeNodeID<T>),TreeNodeID<T>>>();
        for l1 in tree.get_node_ids(){
            for l2 in tree.get_node_ids(){
                let lca = tree.get_lca_id(vec![l1,l2].as_slice());
                data.insert((l1, l2), lca);
                data.insert((l2, l1), lca);
            }
        }
        LcaMap { data: data }
    }

    pub fn get_lca(&self, l1: TreeNodeID<T>, l2: TreeNodeID<T>)->&TreeNodeID<T>{
        // dbg!(&self.data, &(l1, l2), self.data.get(&(l1, l2)));
        self.data.get(&(l1, l2)).unwrap()
    }
}

#[derive(Debug, Clone)]
pub enum NlcdAttribute<T: Float + NumCast + Signed + Zero> {
    Norms(Vec<T>),
    Count(u32),
    Kappa(T),
}

impl<T: Float + NumCast + Signed + Zero> NlcdAttribute<T> {
    fn reset(&mut self) {
        match self {
            Self::Norms(c) => {
                for i in c {
                    *i = T::zero();
                }
            }
            Self::Count(c) => {
                *c = 0;
            }
            Self::Kappa(c) => {
                *c = T::zero();
            }
        }
    }
}

impl<T: Float + NumCast + Signed + Zero> AddAssign for NlcdAttribute<T> {
    fn add_assign(&mut self, other: Self) {
        match other {
            NlcdAttribute::Norms(o) => match self {
                NlcdAttribute::Norms(c) => {
                    if c.len() != o.len() {
                        panic!("Vectors of different lenghts cannot be element-wise added!")
                    }
                    for i in 0..c.len() {
                        c[i] = c[i] + o[i];
                    }
                }
                _ => {
                    panic!("Attempted to add scalar to vec!")
                }
            },
            NlcdAttribute::Count(o) => match self {
                NlcdAttribute::Count(c) => {
                    *c += o;
                }
                _ => {
                    panic!("Attempted to add vec to scalar!")
                }
            },
            NlcdAttribute::Kappa(o) => match self {
                NlcdAttribute::Kappa(c) => {
                    *c = *c + o;
                }
                _ => {
                    panic!("Attempted to add scalar to vec!")
                }
            },
        }
    }
}

pub enum NlcdAttributeType {
    Sigma,
    SigmaPos,
    SigmaNeg,
    Delta,
    CounterpartCount,
    Kappa,
}

pub trait NlcdNodeAttributes<T: Float + NumCast + Signed + Zero>:
    IndexMut<NlcdAttributeType, Output = NlcdAttribute<T>>
{
    fn reset(&mut self) {
        for i in vec![
            NlcdAttributeType::Sigma,
            NlcdAttributeType::SigmaPos,
            NlcdAttributeType::SigmaNeg,
            NlcdAttributeType::Delta,
            NlcdAttributeType::CounterpartCount,
            NlcdAttributeType::Kappa,
        ]
        .into_iter()
        {
            self[i].reset();
        }
    }

    fn get_sigma(&self) -> &[T] {
        match &self[NlcdAttributeType::Sigma] {
            NlcdAttribute::Norms(c) => c.as_slice(),
            _ => panic!("Oops"),
        }
    }
    fn get_sigma_pos(&self) -> &[T] {
        match &self[NlcdAttributeType::SigmaPos] {
            NlcdAttribute::Norms(c) => c.as_slice(),
            _ => panic!("Oops"),
        }
    }
    fn get_sigma_neg(&self) -> &[T] {
        match &self[NlcdAttributeType::SigmaNeg] {
            NlcdAttribute::Norms(c) => c.as_slice(),
            _ => panic!("Oops"),
        }
    }
    fn get_delta(&self) -> &[T] {
        match &self[NlcdAttributeType::Delta] {
            NlcdAttribute::Norms(c) => c.as_slice(),
            _ => panic!("Oops"),
        }
    }
    fn get_counterpart_count(&self) -> u32 {
        match &self[NlcdAttributeType::CounterpartCount] {
            NlcdAttribute::Count(c) => *c,
            _ => panic!("Oops"),
        }
    }
    fn get_kappa(&self) -> &T {
        match &self[NlcdAttributeType::Kappa] {
            NlcdAttribute::Kappa(c) => c,
            _ => panic!("Oops"),
        }
    }

    fn set_sigma(&mut self, value: Vec<T>) {
        self[NlcdAttributeType::Sigma] = NlcdAttribute::Norms(value);
    }
    fn set_sigma_pos(&mut self, value: Vec<T>) {
        self[NlcdAttributeType::SigmaPos] = NlcdAttribute::Norms(value);
    }
    fn set_sigma_neg(&mut self, value: Vec<T>) {
        self[NlcdAttributeType::SigmaNeg] = NlcdAttribute::Norms(value);
    }
    fn set_delta(&mut self, value: Vec<T>) {
        self[NlcdAttributeType::Delta] = NlcdAttribute::Norms(value);
    }
    fn set_counterpart_count(&mut self, value: u32) {
        self[NlcdAttributeType::CounterpartCount] = NlcdAttribute::Count(value);
    }
    fn set_kappa(&mut self, value: T) {
        self[NlcdAttributeType::Kappa] = NlcdAttribute::Kappa(value);
    }
}

pub trait NlcdTreeAttributes<U, T: Float + NumCast + Signed + Zero>:
    IndexMut<U, Output: NlcdNodeAttributes<T>>
{
    fn get_sigma<'a>(&'a self, node_id: U) -> &'a [T] where <Self as Index<U>>::Output: 'a {
        match &self[node_id][NlcdAttributeType::Sigma] {
            NlcdAttribute::Norms(c) => c.as_slice(),
            _ => panic!("Oops"),
        }
    }

    fn get_sigma_pos<'a>(&'a self, node_id: U) -> &'a [T] where <Self as Index<U>>::Output: 'a {
        match &self[node_id][NlcdAttributeType::SigmaPos] {
            NlcdAttribute::Norms(c) => c.as_slice(),
            _ => panic!("Oops"),
        }
    }

    fn get_sigma_neg<'a>(&'a self, node_id: U) -> &'a [T] where <Self as Index<U>>::Output: 'a {
        match &self[node_id][NlcdAttributeType::SigmaNeg] {
            NlcdAttribute::Norms(c) => c.as_slice(),
            _ => panic!("Oops"),
        }
    }

    fn get_delta<'a>(&'a self, node_id: U) -> &'a [T] where <Self as Index<U>>::Output: 'a {
        match &self[node_id][NlcdAttributeType::Delta] {
            NlcdAttribute::Norms(c) => c.as_slice(),
            _ => panic!("Oops"),
        }
    }

    fn get_counterpart_count(&self, node_id: U) -> u32 {
        match &self[node_id][NlcdAttributeType::CounterpartCount] {
            NlcdAttribute::Count(c) => *c,
            _ => panic!("Oops"),
        }
    }
    fn get_kappa(&self, node_id: U) -> T {
        match &self[node_id][NlcdAttributeType::Kappa] {
            NlcdAttribute::Kappa(c) => *c,
            _ => panic!("Oops"),
        }
    }

    fn increment_count(&mut self, node_id: U) {
        self[node_id][NlcdAttributeType::CounterpartCount] += NlcdAttribute::Count(1);
    }

    fn increment_sigma(&mut self, node_id: U, value: Vec<T>) {
        self[node_id][NlcdAttributeType::Sigma] += NlcdAttribute::Norms(value);
    }

    fn increment_sigma_pos(&mut self, node_id: U, value: Vec<T>) {
        self[node_id][NlcdAttributeType::SigmaPos] += NlcdAttribute::Norms(value);
    }

    fn increment_sigma_neg(&mut self, node_id: U, value: Vec<T>) {
        self[node_id][NlcdAttributeType::SigmaNeg] += NlcdAttribute::Norms(value);
    }

    fn increment_delta(&mut self, node_id: U, value: Vec<T>) {
        self[node_id][NlcdAttributeType::Delta] += NlcdAttribute::Norms(value);
    }

    fn set_sigma(&mut self, node_id: U, value: Vec<T>) {
        self[node_id][NlcdAttributeType::Sigma] = NlcdAttribute::Norms(value);
    }

    fn set_sigma_pos(&mut self, node_id: U, value: Vec<T>) {
        self[node_id][NlcdAttributeType::SigmaPos] = NlcdAttribute::Norms(value);
    }

    fn set_sigma_neg(&mut self, node_id: U, value: Vec<T>) {
        self[node_id][NlcdAttributeType::SigmaNeg] = NlcdAttribute::Norms(value);
    }

    fn set_delta(&mut self, node_id: U, value: Vec<T>) {
        self[node_id][NlcdAttributeType::Delta] = NlcdAttribute::Norms(value);
    }

    fn set_counterpart_count(&mut self, node_id: U, value: u32) {
        self[node_id][NlcdAttributeType::CounterpartCount] = NlcdAttribute::Count(value);
    }
    fn set_kappa(&mut self, node_id: U, value: T) {
        self[node_id][NlcdAttributeType::Kappa] = NlcdAttribute::Kappa(value);
    }

    fn reset_node(&mut self, node_id: U) {
        self[node_id].reset()
    }
}

pub trait NearLinearCopheneticDistance: CopheneticDistance
where
    <Self as RootedTree>::Node: RootedMetaNode + RootedZetaNode,
    <<Self as RootedTree>::Node as RootedZetaNode>::Zeta: Signed
        + Clone
        + NumCast
        + std::iter::Sum
        + Debug
        + Display
        + Float
        + PartialOrd
        + Copy
        + Send,
{
    // helper functions
    fn tree_attributes<'a>(
        &'a self,
        norm: u32,
    ) -> impl NlcdTreeAttributes<
        TreeNodeID<Self>,
        TreeNodeZeta<Self>,
    >;

    /// Returns taxa present in upper tree.
    fn upper_tree_taxa(&self, median_node_id: TreeNodeID<Self>) -> impl Iterator<Item = TreeNodeMeta<Self>> {
        let lower_tree_taxa = self.lower_tree_taxa(median_node_id).collect::<HashSet<_>>();
        self.get_leaf_ids()
            .map(|x| self.get_node_taxa_cloned(x).unwrap())
            .filter(move |x| !lower_tree_taxa.contains(x))
    }

    /// Returns taxa present in lower tree.
    fn lower_tree_taxa(&self, median_node_id: TreeNodeID<Self>) -> impl Iterator<Item = TreeNodeMeta<Self>> {
        self.get_cluster_ids(median_node_id)
            .filter_map(|x| self.get_node_taxa_cloned(x))
    }

    /// Returns taxa present in upper tree.
    fn upper_tree_leaves(&self, median_node_id: TreeNodeID<Self>) -> impl Iterator<Item = <Self as RootedTree>::Node> {
        let lower_tree_leaf_ids = self
            .lower_tree_leaves(median_node_id)
            .map(|x| x.get_id())
            .collect::<HashSet<_>>();
        self.get_leaves()
            .filter(move |x| !lower_tree_leaf_ids.contains(&x.get_id()))
            .cloned()
    }

    /// Returns leaves present in lower tree.
    fn lower_tree_leaves(&self, median_node_id: TreeNodeID<Self>) -> impl Iterator<Item = &<Self as RootedTree>::Node> {
        self.get_cluster(median_node_id)
    }

    /// Returns lower tree.
    fn lower_tree(&self, median_node_id: TreeNodeID<Self>) -> Self;


    /// Returns upper tree.
    fn upper_tree(&self, median_node_id: TreeNodeID<Self>) -> Self;

    
    fn get_node_sibling_id(
        &self,
        node_id: &TreeNodeID<Self>,
    ) -> TreeNodeID<Self> {
        let node_parent_id = self
            .get_node_parent_id(*node_id)
            .expect("Node has no siblings");
        let siblings = self.get_node_children_ids(node_parent_id).collect_vec();
        if &siblings[0] == node_id {
            siblings[1]
        } else {
            siblings[0]
        }
    }

    /// Returns value of n choose k.
    fn n_choose_k(n: u32, k: u32) -> i32 {
        fn multiplicative_form(n: u32, k: u32) -> i32 {
            (1..k + 1)
                .map(|x| (n + 1 - x) as f32 / (x as f32))
                .product::<f32>() as i32
        }

        match k <= (n as f32 / 2_f32) as u32 {
            true => multiplicative_form(n, k),
            false => multiplicative_form(n, n - k),
        }
    }

    /// Returns pascal triangle for precomputed binomial terms
    fn pascal_triangle(norm: u32)->Vec<Vec<u32>>{
        let mut pt = vec![vec![0;norm as usize+1]; norm as usize+1];
        pt[0][0]=1;
        for n in 1..norm as usize+1{
            pt[n][0] = 1;
            pt[n][n] = 1;
            for k in 1..n{
                let n_choose_k = pt[n-1][k-1]+pt[n-1][k];
                pt[n][k] = n_choose_k;
            }
        }
        pt
    }

    //  Naive algorithm
    fn naive_cophen_dist(
        &self,
        tree: &Self,
        self_lca_map: &LcaMap<Self>,
        tree_lca_map: &LcaMap<Self>,
        norm: u32,
    ) -> TreeNodeZeta<Self> {
        let mut distances: TreeNodeZeta<Self> = <TreeNodeZeta<Self>>::zero();

        let taxa_set = self
            .get_taxa_space()
            // .map(|x| x.clone())
            .collect_vec();

        let num_taxa = taxa_set.len();

        let mut seen_pairs: HashSet<(usize, usize)> = vec![].into_iter().collect();

        for i in 0..num_taxa{
            for j in 0..num_taxa{
                match seen_pairs.contains(&(i,j)){
                    true => {},
                    false => {
                        seen_pairs.insert((i,j));
                        seen_pairs.insert((j,i));

                        let lca_1 = self_lca_map.get_lca(self.get_taxa_node_id(taxa_set[i]).unwrap(), self.get_taxa_node_id(taxa_set[j]).unwrap());
                        let lca_2 = tree_lca_map.get_lca(tree.get_taxa_node_id(taxa_set[i]).unwrap(), tree.get_taxa_node_id(taxa_set[j]).unwrap());
                        let zeta_1 = self.get_zeta(*lca_1).unwrap();
                        let zeta_2 = tree.get_zeta(*lca_2).unwrap();
                        let zeta = (zeta_1 - zeta_2).abs();
                        if norm==0{
                            distances = match zeta >= distances{
                                true => {zeta}
                                false => {distances}
                            };
                        }
                        else if norm == 1 {
                            distances = distances + zeta;
                        }
                        else {
                            distances = distances + zeta.powi(norm as i32);
                        }
        
                    },
                }
            }
        }

        distances.powf(
            <TreeNodeZeta<Self> as NumCast>::from(norm)
                .unwrap()
                .powi(-1))
    }

    /// Returns the Cophenetic distance between two trees in O(nlog^2n + pnlogn) time.
    fn nl_cophen_dist(
        &self,
        tree: &Self,
        self_lca_map: &LcaMap<Self>,
        tree_lca_map: &LcaMap<Self>,
        // pascal_triangle: &Vec<Vec<u32>>,
        norm: u32,
    ) -> TreeNodeZeta<Self> {
        let mut distances: TreeNodeZeta<Self> = <TreeNodeZeta<Self>>::zero();

        let taxa_set = self
            .get_taxa_space()
            .map(|x| x.clone())
            .collect::<HashSet<TreeNodeMeta<Self>>>();


        let mut taxa_map: HashMap<&TreeNodeMeta<Self>, usize> = vec![].into_iter().collect();
        let mut taxa_map_rev: HashMap<usize,TreeNodeMeta<Self>> = vec![].into_iter().collect();


        for (idx, x) in taxa_set.iter().enumerate(){
            taxa_map.insert(x, idx);
            taxa_map_rev.insert(idx,x.clone());
        }

        let (self_node_cluster_bitsets, self_postord_node_ids, t) = Self::cluster_bitsets_and_median(self, &taxa_map);
        let (tree_node_cluster_bitsets, tree_postord_node_ids, t_hat) = Self::cluster_bitsets_and_median(tree, &taxa_map);

        let self_lower_set: HashSet<_> = self.postord_ids(t).collect();
        let self_lower_node_ids_postord = self.postord_ids(t).collect_vec();
        let self_upper_node_ids_postord = self_postord_node_ids.iter().filter(|x| !self_lower_set.contains(x))
        .map(|x| x.clone())
        .collect_vec();

        let tree_lower_set: HashSet<_> = tree.postord_ids(t_hat).collect();
        let tree_lower_node_ids_postord = tree.postord_ids(t_hat).collect_vec();
        let tree_upper_node_ids_postord = tree_postord_node_ids.iter().filter(|x| !tree_lower_set.contains(x))
        .map(|x| x.clone())
        .collect_vec();

        let mut self_node_attributes = self.tree_attributes(norm);
        let mut tree_node_attributes = tree.tree_attributes(norm);

        // Precompute pascal triangle
        let pt = Self::pascal_triangle(norm);

        self.populate_op_vec(
            tree,
            self_lca_map,
            tree_lca_map,
            t,
            t_hat,
            self_node_cluster_bitsets,
            tree_node_cluster_bitsets,
            self_postord_node_ids,
            self_lower_node_ids_postord, 
            self_upper_node_ids_postord,
            tree_postord_node_ids,
            tree_lower_node_ids_postord, 
            tree_upper_node_ids_postord,
            taxa_map_rev,
            &pt,
            norm,
            &taxa_set,
            &mut distances,
            &mut self_node_attributes,
            &mut tree_node_attributes,
        );

        distances.powf(
            <TreeNodeZeta<Self> as NumCast>::from(norm)
                .unwrap()
                .powi(-1),
        )
    }

    /// returns bitvec representation of all vertices and vertex ids in postorder.
    fn cluster_bitsets_and_median(tree: &Self, taxa_map: &HashMap<&TreeNodeMeta<Self>, usize>)-> (HashMap<TreeNodeID<Self>, BitVec>, Vec<TreeNodeID<Self>>, TreeNodeID<Self>){
        let num_leaves = taxa_map.len();

        let mut tree_postord_node_ids = Vec::with_capacity(2*num_leaves);
        let mut tree_node_cluster_bitsets: HashMap<TreeNodeID<Self>, BitVec> = vec![].into_iter().collect();

        for n_id in tree.postord_ids(tree.get_root_id()){
            tree_postord_node_ids.push(n_id);
            let mut bitset = BitVec::from_zeros(num_leaves);
            match tree.is_leaf(n_id){
                true => {
                    let node_taxa = tree.get_node_taxa(n_id).unwrap();
                    let taxa_id = taxa_map.get(node_taxa).unwrap();
                    bitset.flip_bit(*taxa_id);
                    tree_node_cluster_bitsets.insert(n_id, bitset);
                },
                false => {
                    tree.get_node_children_ids(n_id)
                        .map(|c_id| tree_node_cluster_bitsets.get(&c_id).unwrap())
                        .for_each(|c_bitset| {let _ = bitset.apply_mask_or(c_bitset);});
                    tree_node_cluster_bitsets.insert(n_id, bitset);
                },
            };
        }

        let mut median_node_id: TreeNodeID<Self> = tree.get_root_id();
        loop {
            median_node_id = tree.get_node_children_ids(median_node_id)
                .max_by(|x, y| {
                    let x_cluster_size = tree_node_cluster_bitsets.get(x).unwrap().count_ones();
                    let y_cluster_size = tree_node_cluster_bitsets.get(y).unwrap().count_ones();
                    x_cluster_size.cmp(&y_cluster_size)
                })
                .unwrap();
            if tree_node_cluster_bitsets.get(&median_node_id).unwrap().count_ones() as usize <= (num_leaves / 2) {
                break;
            }
        }
        (tree_node_cluster_bitsets, tree_postord_node_ids, median_node_id)
    }

    fn get_upper_lower_sets<'a>(t: &TreeNodeID<Self>, t_hat: &TreeNodeID<Self>, self_node_cluster_bitsets: &'a HashMap<TreeNodeID<Self>, BitVec>, tree_node_cluster_bitsets: &'a HashMap<TreeNodeID<Self>, BitVec>) -> (&'a BitVec, &'a BitVec, BitVec, BitVec, BitVec, BitVec, BitVec, BitVec){

        let b_bitvec = self_node_cluster_bitsets.get(&t).unwrap();
        let b_hat_bitvec = tree_node_cluster_bitsets.get(&t_hat).unwrap();
        let a_bitvec = b_bitvec.mask_xor(&BitVec::from_ones(b_bitvec.len())).unwrap().to_bit_vec();
        let a_hat_bitvec = b_hat_bitvec.mask_xor(&BitVec::from_ones(b_bitvec.len())).unwrap().to_bit_vec();

        let a_int_a_hat_bitvec = a_bitvec.mask_and(&a_hat_bitvec).unwrap().to_bit_vec();
        let a_int_b_hat_bitvec = a_bitvec.mask_and(&b_hat_bitvec).unwrap().to_bit_vec();
        let b_int_a_hat_bitvec = b_bitvec.mask_and(&a_hat_bitvec).unwrap().to_bit_vec();
        let b_int_b_hat_bitvec = b_bitvec.mask_and(&b_hat_bitvec).unwrap().to_bit_vec();

        (b_bitvec, b_hat_bitvec, a_bitvec, a_hat_bitvec, a_int_a_hat_bitvec, a_int_b_hat_bitvec, b_int_a_hat_bitvec, b_int_b_hat_bitvec)
    }

    fn bitvec_to_set(bitvec: &BitVec, taxa_map_rev: &HashMap<usize, TreeNodeMeta<Self>>)->HashSet<TreeNodeMeta<Self>>{
        (0..taxa_map_rev.len()).filter(|idx| bitvec.is_bit_set_unchecked(*idx)).map(|x| taxa_map_rev.get(&x).cloned().unwrap()).collect()
    }

    fn create_taxa_maps(taxa_set: &HashSet<TreeNodeMeta<Self>>)->(HashMap<&TreeNodeMeta<Self>, usize>, HashMap<usize,TreeNodeMeta<Self>>){
        let mut taxa_map: HashMap<&TreeNodeMeta<Self>, usize> = vec![].into_iter().collect();
        let mut taxa_map_rev: HashMap<usize,TreeNodeMeta<Self>> = vec![].into_iter().collect();

        for (idx, x) in taxa_set.iter().enumerate(){
            taxa_map.insert(x, idx);
            taxa_map_rev.insert(idx,x.clone());
        }

        (taxa_map, taxa_map_rev)
    }

    /// Combining tree contraction, median vertex, and cluster bit vectors for each vertex of tree.
    fn contract_and_median(tree: &Self, taxa_map: &HashMap<&TreeNodeMeta<Self>, usize>, tree_postord_node_ids: &[TreeNodeID<Self>])->(Self, TreeNodeID<Self>, HashMap<TreeNodeID<Self>, BitVec>, Vec<TreeNodeID<Self>>, Vec<TreeNodeID<Self>>, Vec<TreeNodeID<Self>>);

    /// Populates vector with divide and conquer distances.
    fn populate_op_vec(
        &self,
        tree: &Self,
        self_lca_map: &LcaMap<Self>,
        tree_lca_map: &LcaMap<Self>,
        t: TreeNodeID<Self>,
        t_hat: TreeNodeID<Self>,
        self_node_cluster_bitsets: HashMap<TreeNodeID<Self>, BitVec>,
        tree_node_cluster_bitsets: HashMap<TreeNodeID<Self>, BitVec>,
        self_postord_node_ids: Vec<TreeNodeID<Self>>,
        self_lower_postord_node_ids: Vec<TreeNodeID<Self>>,
        self_upper_postord_node_ids: Vec<TreeNodeID<Self>>,
        tree_postord_node_ids: Vec<TreeNodeID<Self>>,
        tree_lower_postord_node_ids: Vec<TreeNodeID<Self>>,
        tree_upper_postord_node_ids: Vec<TreeNodeID<Self>>,
        taxa_map_rev: HashMap<usize,TreeNodeMeta<Self>>,
        pascal_triangle: &Vec<Vec<u32>>,
        norm: u32,
        taxa_set: &HashSet<TreeNodeMeta<Self>>,
        distances: &mut TreeNodeZeta<Self>,
        self_node_attributes: &mut impl NlcdTreeAttributes<
            TreeNodeID<Self>,
            TreeNodeZeta<Self>,
        >,
        tree_node_attributes: &mut impl NlcdTreeAttributes<
            TreeNodeID<Self>,
            TreeNodeZeta<Self>,
        >,
    ) {
        let (b_bitvec, b_hat_bitvec, a_bitvec, a_hat_bitvec, a_int_a_hat_bitvec, a_int_b_hat_bitvec, b_int_a_hat_bitvec, b_int_b_hat_bitvec) = Self::get_upper_lower_sets(&t, &t_hat, &self_node_cluster_bitsets, &tree_node_cluster_bitsets);

        let b: HashSet<TreeNodeMeta<Self>> = Self::bitvec_to_set(&b_bitvec, &taxa_map_rev);
        let b_hat: HashSet<TreeNodeMeta<Self>> = Self::bitvec_to_set(&b_hat_bitvec, &taxa_map_rev);
        let a: HashSet<TreeNodeMeta<Self>> = Self::bitvec_to_set(&a_bitvec, &taxa_map_rev);
        let a_hat: HashSet<TreeNodeMeta<Self>> = Self::bitvec_to_set(&a_hat_bitvec, &taxa_map_rev);

        let a_int_a_hat: HashSet<TreeNodeMeta<Self>> = Self::bitvec_to_set(&a_int_a_hat_bitvec, &taxa_map_rev);
        let a_int_b_hat: HashSet<TreeNodeMeta<Self>> = Self::bitvec_to_set(&a_int_b_hat_bitvec, &taxa_map_rev);
        let b_int_a_hat: HashSet<TreeNodeMeta<Self>> = Self::bitvec_to_set(&b_int_a_hat_bitvec, &taxa_map_rev);
        let b_int_b_hat: HashSet<TreeNodeMeta<Self>> = Self::bitvec_to_set(&b_int_b_hat_bitvec, &taxa_map_rev);
        
        let double_mix_distance = self.distance_double_mix_type(
            tree,
            self_lca_map,
            tree_lca_map,
            pascal_triangle,
            norm,
            &t,
            &t_hat,
            &a_int_a_hat,
            &a_int_b_hat,
            &b_int_a_hat,
            &b_int_b_hat,
        );

        let single_mix_distance = self.distance_single_mix_type(
            tree,
            self_lca_map,
            tree_lca_map,
            pascal_triangle,
            norm,
            &t,
            &t_hat,
            &a,
            &b,
            &a_hat,
            &b_hat,
            self_lower_postord_node_ids.as_slice(),
            self_upper_postord_node_ids.as_slice(),
            tree_lower_postord_node_ids.as_slice(),
            tree_upper_postord_node_ids.as_slice(),
            self_node_attributes,
            tree_node_attributes,
            &a_int_a_hat,
            &a_int_b_hat,
            &b_int_a_hat,
        );

        *distances = *distances + double_mix_distance+single_mix_distance;

        if taxa_set.len() > 2 {
            for taxa_subset in [&a_int_a_hat, &a_int_b_hat, &b_int_a_hat, &b_int_b_hat]{
                if taxa_subset.len() > 1 {
                    let (taxa_map, taxa_map_rev) = Self::create_taxa_maps(taxa_subset);
                    let (self_tree, t, self_node_cluster_bitsets, self_postord_node_ids, self_lower_postord, self_upper_postord) = Self::contract_and_median(&self, &taxa_map, self_postord_node_ids.as_slice());
                    let (new_tree, t_hat, tree_node_cluster_bitsets, tree_postord_node_ids, tree_lower_postord, tree_upper_postord) = Self::contract_and_median(&tree, &taxa_map, tree_postord_node_ids.as_slice());
                    self_tree.populate_op_vec(
                        &new_tree,
                        self_lca_map,
                        tree_lca_map,
                        t,
                        t_hat,
                        self_node_cluster_bitsets,
                        tree_node_cluster_bitsets,
                        self_postord_node_ids,
                        self_lower_postord, 
                        self_upper_postord,
                        tree_postord_node_ids,
                        tree_lower_postord, 
                        tree_upper_postord,
                        taxa_map_rev,
                        pascal_triangle,
                        norm,
                        &taxa_set,
                        distances,
                        self_node_attributes,
                        tree_node_attributes,
                    );
                }
                if taxa_subset.len()==1{
                    let x = taxa_subset.iter().next().unwrap();
                    let zeta_1 = self.get_zeta_taxa(&x);
                    let zeta_2 = tree.get_zeta_taxa(&x);
                    *distances = *distances + (zeta_1 - zeta_2).abs().powi(norm as i32);
                }
            }
        }
    }

    /// Returns ordered iterator used in double mix type cases
    fn get_cntr(
        &self,
        self_lca_map: &LcaMap<Self>,
        median_node_id: TreeNodeID<Self>,
        leaf_set: HashSet<TreeNodeID<Self>>,
    ) -> Vec<TreeNodeZeta<Self>> {
        // line 5 in algo 1
        let mut gamma: Vec<TreeNodeZeta<Self>> = Vec::with_capacity(2*leaf_set.len());
        // line 3 in algo 1
        let med_path = self.root_to_node_ids(median_node_id).collect_vec();
        let mut median_path = med_path
            .iter()
            .map(|x| (x.clone(), 0))
            .collect::<HashMap<_, _>>();
        for node_id in leaf_set {
            // line 4 in algo 1
            median_path
                .entry(*self_lca_map.get_lca(node_id, median_node_id))
                .and_modify(|x| *x += 1);
        }
        for node_id in med_path {
            let c = median_path.get(&node_id).cloned().unwrap();
            for _ in 0..c {
                gamma.push(self.get_zeta(node_id).unwrap())
            }
        }
        gamma.into()
    }

    /// Returns seqPrd^p(\alpha,\beta) used for double mix type cases.
    fn seq_product(
        mut alpha: Vec<TreeNodeZeta<Self>>,
        mut beta: Vec<TreeNodeZeta<Self>>,
        pascal_triangle: &Vec<Vec<u32>>,
        norm: u32,
    ) -> TreeNodeZeta<Self> {
        if alpha.is_empty() || beta.is_empty() {
            return <TreeNodeZeta<Self> as Zero>::zero();
        }
        if alpha.last().unwrap() < &alpha[0] {
            alpha.reverse();
        }
        if beta.last().unwrap() < &beta[0] {
            beta.reverse();
        }
        if alpha.last().unwrap() > beta.last().unwrap() {
            std::mem::swap(&mut alpha, &mut beta);
        }

        let sigma = Self::precompute_alpha_sums(&alpha, &beta, norm);

        let final_out = (0..beta.len()).map(|j| {
            let term_1 = (0..norm+1).map(|l| {
                let t1 = <TreeNodeZeta<Self> as NumCast>::from(pascal_triangle[norm as usize][l as usize]).unwrap();
                let t2 = beta[j].powi(l as i32);
                let t3 = <TreeNodeZeta<Self> as NumCast>::from((-1_i32).pow(norm-l)).unwrap();
                let t4 = sigma[j][(norm-l) as usize];
                t1*t2*t3*t4
            }).sum::<TreeNodeZeta<Self>>();

            let term_2 = (0..norm+1).map(|l| {
                let t1 = <TreeNodeZeta<Self> as NumCast>::from(pascal_triangle[norm as usize][l as usize]).unwrap();
                let t2 = beta[j].powi((norm-l) as i32);
                let t3 = <TreeNodeZeta<Self> as NumCast>::from((-1_i32).pow(norm-l)).unwrap();
                let t4 = sigma[beta.len()-1][l as usize]-sigma[j][l as usize];
                t1*t2*t3*t4
            }).sum::<TreeNodeZeta<Self>>();

            term_1+term_2
        }).sum::<TreeNodeZeta<Self>>();

        final_out
    }

    /// Algorithm 2 for precomputation of alpha sums
    fn precompute_alpha_sums(
        alpha: &[TreeNodeZeta<Self>],
        beta: &[TreeNodeZeta<Self>],
        norm: u32,
    ) -> Vec<Vec<TreeNodeZeta<Self>>> {
        let mut j = 0;
        let mut i = 0;
        let k = alpha.len();
        let m = beta.len();

        let mut sigma = vec![vec![<<<Self as RootedTree>::Node as RootedZetaNode>::Zeta>::zero();(norm + 1) as usize];m];

        while j < m {
            match i < k && alpha[i] <= beta[j] {
                true => {
                    (0..norm + 1).for_each(|l| {
                        sigma[j][l as usize] = sigma[j][l as usize] + alpha[i].powi(l as i32)
                    });
                    i += 1;
                }
                false => {
                    j += 1;
                    if j < m {
                        (0..norm + 1).for_each(|l| sigma[j][l as usize] = sigma[j - 1][l as usize]);
                    }
                }
            }
        }

        sigma
    }

    /// This method generates the distance contributed by all taxa pairs
    /// that are present in different subtrees in both trees(raised to the p^{th} power).
    ///
    /// This includes the following assignments: AB|A'B', AB|B'A'
    #[allow(clippy::too_many_arguments)]
    fn distance_double_mix_type(
        &self,
        tree: &Self,
        self_lca_map: &LcaMap<Self>,
        tree_lca_map: &LcaMap<Self>,
        pascal_triangle: &Vec<Vec<u32>>,
        norm: u32,
        t: &TreeNodeID<Self>,
        t_hat: &TreeNodeID<Self>,
        a_int_a_hat: &HashSet<TreeNodeMeta<Self>>,
        a_int_b_hat: &HashSet<TreeNodeMeta<Self>>,
        b_int_a_hat: &HashSet<TreeNodeMeta<Self>>,
        b_int_b_hat: &HashSet<TreeNodeMeta<Self>>,
    ) -> TreeNodeZeta<Self> {
        let alpha = self.get_cntr(
                self_lca_map,
            *t,
            a_int_b_hat
                .iter()
                .map(|x| self.get_taxa_node_id(x).unwrap())
                .collect::<HashSet<TreeNodeID<Self>>>(),
        );
        let beta = tree.get_cntr(
                tree_lca_map,
                *t_hat,
            b_int_a_hat
                .iter()
                .map(|x| tree.get_taxa_node_id(x).unwrap())
                .collect::<HashSet<TreeNodeID<Self>>>(),
        );

        // AB|A'B'
        let b_int_b_hat_len = b_int_b_hat.len();
        let dd2 = a_int_a_hat
            .iter()
            .map(|x| {
                let t_lca_id = self_lca_map.get_lca(self.get_taxa_node_id(x).unwrap(), *t);
                let t_hat_lca_id =
                    tree_lca_map.get_lca(tree.get_taxa_node_id(x).unwrap(), *t_hat);
                let zeta_1 = self.get_zeta(*t_lca_id).unwrap();
                let zeta_2 = tree.get_zeta(*t_hat_lca_id).unwrap();
                (zeta_1 - zeta_2).abs().powi(norm as i32)
            })
            .sum::<TreeNodeZeta<Self>>()
            * <TreeNodeZeta<Self> as NumCast>::from(
                b_int_b_hat_len,
            )
            .unwrap();

        Self::seq_product(alpha, beta, pascal_triangle, norm) + dd2
    }

    /// This method generates the distance contributed by all taxa pairs
    /// that are present in the same subtree in exactly one of the two trees(raised to the p^{th} power).
    ///
    /// This includes the following assignments: AA|A'B', AA|B'A', BB|A'B', BB|B'A', BA|B'B', BA|A'A', AB|B'B', AB|A'A'.
    #[allow(clippy::too_many_arguments)]
    fn distance_single_mix_type(
        &self,
        tree: &Self,
        self_lca_map: &LcaMap<Self>,
        tree_lca_map: &LcaMap<Self>,
        pascal_triangle: &Vec<Vec<u32>>,
        norm: u32,
        t: &TreeNodeID<Self>,
        t_hat: &TreeNodeID<Self>,
        t1_upper_taxa: &HashSet<TreeNodeMeta<Self>>,
        t1_lower_taxa: &HashSet<TreeNodeMeta<Self>>,
        t2_upper_taxa: &HashSet<TreeNodeMeta<Self>>,
        t2_lower_taxa: &HashSet<TreeNodeMeta<Self>>,
        t1_lower_postord_node_ids: &[TreeNodeID<Self>],
        t1_upper_postord_node_ids: &[TreeNodeID<Self>],
        t2_lower_postord_node_ids: &[TreeNodeID<Self>],
        t2_upper_postord_node_ids: &[TreeNodeID<Self>],
        self_node_attributes: &mut impl NlcdTreeAttributes<
            TreeNodeID<Self>,
            TreeNodeZeta<Self>,
        >,
        tree_node_attributes: &mut impl NlcdTreeAttributes<
            TreeNodeID<Self>,
            TreeNodeZeta<Self>,
        >,
        a_int_a_hat: &HashSet<TreeNodeMeta<Self>>,
        a_int_b_hat: &HashSet<TreeNodeMeta<Self>>,
        b_int_a_hat: &HashSet<TreeNodeMeta<Self>>,
    ) -> TreeNodeZeta<Self> {
        // AA|A'B'
        let d1 = Self::single_mix_xxxy(
            self,
            tree,
            tree_lca_map,
            t,
            t_hat,
            self_node_attributes,
            t1_upper_taxa,
            t1_lower_taxa,
            t1_lower_postord_node_ids,
            t1_upper_postord_node_ids,
            pascal_triangle,
            norm,
            a_int_a_hat,
            true,
        );
        // AB|A'A'
        let d2 = Self::single_mix_xxxy(
            tree,
            self,
            self_lca_map,
            t_hat,
            t,
            tree_node_attributes,
            t2_upper_taxa,
            t2_lower_taxa,
            t2_lower_postord_node_ids,
            t2_upper_postord_node_ids,
            pascal_triangle,
            norm,
            a_int_a_hat,
            true,
        );
        // BB|A'B'
        let d3 = Self::single_mix_xxxy(
            self,
            tree,
            tree_lca_map,
            t,
            t_hat,
            self_node_attributes,
            t1_upper_taxa,
            t1_lower_taxa,
            t1_lower_postord_node_ids,
            t1_upper_postord_node_ids,
            pascal_triangle,
            norm,
            b_int_a_hat,
            false,
        );
        // AB|B'B'
        let d4 = Self::single_mix_xxxy(
            tree,
            self,
            self_lca_map,
            t_hat,
            t,
            tree_node_attributes,
            t2_upper_taxa,
            t2_lower_taxa,
            t2_lower_postord_node_ids,
            t2_upper_postord_node_ids,
            pascal_triangle,
            norm,
            a_int_b_hat,
            false,
        );
        d1 + d2 + d3 + d4
    }

    #[allow(clippy::too_many_arguments)]
    fn preprocess_single_mix_even(
        t1: &Self,
        t2: &Self,
        t2_lca_map: &LcaMap<Self>,
        t2_median: &TreeNodeID<Self>,
        t1_node_attributes: &mut impl NlcdTreeAttributes<
            TreeNodeID<Self>,
            TreeNodeZeta<Self>,
        >,
        t1_upper_taxa: &HashSet<TreeNodeMeta<Self>>,
        t1_lower_taxa: &HashSet<TreeNodeMeta<Self>>,
        norm: u32,
        taxa_set: &HashSet<TreeNodeMeta<Self>>,
        upper_mixed: bool,
    ) {
        let leaf_iter = match upper_mixed {
            true => t1_upper_taxa,
            false => t1_lower_taxa,
        };
        for leaf in leaf_iter {
            match taxa_set.contains(&leaf) {
                true => {
                    let t1_node_id = t1.get_taxa_node_id(&leaf).unwrap();
                    let t2_node_id = t2.get_taxa_node_id(&leaf).unwrap();
                    let lca_x_t_hat = t2_lca_map.get_lca(t2_node_id, *t2_median);
                    let beta = t2.get_zeta(*lca_x_t_hat).unwrap();
                    t1_node_attributes.set_sigma(
                        t1_node_id,
                        (0..norm + 1).map(|l| beta.powi(l as i32)).collect_vec(),
                    );
                }
                false => {
                    let node_id = t1.get_taxa_node_id(&leaf).unwrap();
                    t1_node_attributes.increment_count(node_id);
                }
            }
        }
    }

    fn binary_search_omega(
        t1: &Self,
        node_ids: Vec<TreeNodeID<Self>>, 
        t2_zeta: TreeNodeZeta<Self>
    )->TreeNodeID<Self>
    {
        let mut high = node_ids.len()-1;
        let mut low = 0;
        let default = node_ids[high];

        while low <= high{
            let mid = low + (high - low)/2;
    
            if mid+1==node_ids.len(){
                return node_ids[mid]
            }

            if t1.get_zeta(node_ids[mid]).unwrap() >= t2_zeta && t1.get_zeta(node_ids[mid+1]).unwrap() < t2_zeta{
                return node_ids[mid];
            }
        
            if t1.get_zeta(node_ids[mid]).unwrap() > t2_zeta{
                low = mid + 1;
            }
        
            else{
                high = mid - 1;
            }
        }

        return default;

    }

    #[allow(clippy::too_many_arguments)]
    fn preprocess_single_mix_odd(
        t1: &Self,
        t2: &Self,
        t2_lca_map: &LcaMap<Self>,
        t1_median: &TreeNodeID<Self>,
        t2_median: &TreeNodeID<Self>,
        t1_node_attributes: &mut impl NlcdTreeAttributes<
            TreeNodeID<Self>,
            TreeNodeZeta<Self>,
        >,
        t1_upper_taxa: &HashSet<TreeNodeMeta<Self>>,
        t1_lower_taxa: &HashSet<TreeNodeMeta<Self>>,
        norm: u32,
        taxa_set: &HashSet<TreeNodeMeta<Self>>,
        upper_mixed: bool,
    ) {
        let leaf_iter = match upper_mixed {
            true => t1_upper_taxa,
            false => t1_lower_taxa,
        };
        for leaf in leaf_iter.iter() {
            match taxa_set.contains(&leaf) {
                true => {
                    let t2_node_id = t2.get_taxa_node_id(leaf).unwrap();
                    let lca_x_t_hat = t2_lca_map.get_lca(t2_node_id,*t2_median);
                    let beta = t2.get_zeta(*lca_x_t_hat).unwrap();
                    let t1_node_id = t1.get_taxa_node_id(leaf).unwrap();

                    let t1_node_parent_id = t1.get_node_parent_id(t1_node_id).unwrap();

                    if beta <= t1.get_zeta(t1_node_id).unwrap() {
                        //find omega_x
                        let mut v_ancestors = vec![];
                        for node_id in t1.node_to_root_ids(t1_node_id).into_iter() {
                            v_ancestors.push(node_id);
                            if node_id == *t1_median {
                                break;
                            }
                        }

                        let omega_x = Self::binary_search_omega(t1, v_ancestors, beta);

                        t1_node_attributes.increment_delta(
                            omega_x,
                            (0..norm + 1).map(|l| beta.powi(l as i32)).collect_vec(),
                        );
                    }
                    if beta <= t1.get_zeta(t1_node_parent_id).unwrap() {
                        t1_node_attributes.increment_sigma_pos(
                            t1_node_id,
                            (0..norm + 1).map(|l| beta.powi(l as i32)).collect_vec(),
                        );
                    } else {
                        t1_node_attributes.increment_sigma_neg(
                            t1_node_id,
                            (0..norm + 1).map(|l| beta.powi(l as i32)).collect_vec(),
                        );
                    }
                }
                false => {
                    let node_id = t1.get_taxa_node_id(leaf).unwrap();
                    t1_node_attributes.increment_count(node_id);
                }
            }
        }
    }

    // this method solves AA|A'B'
    #[allow(clippy::too_many_arguments)]
    fn single_mix_xxxy(
        t1: &Self,
        t2: &Self,
        t2_lca_map: &LcaMap<Self>,
        t1_median: &TreeNodeID<Self>,
        t2_median: &TreeNodeID<Self>,
        t1_node_attributes: &mut impl NlcdTreeAttributes<
            TreeNodeID<Self>,
            TreeNodeZeta<Self>,
        >,
        t1_upper_taxa: &HashSet<TreeNodeMeta<Self>>,
        t1_lower_taxa: &HashSet<TreeNodeMeta<Self>>,
        t1_lower_postord_node_ids: &[TreeNodeID<Self>],
        t1_upper_postord_node_ids: &[TreeNodeID<Self>],
        pascal_triangle: &Vec<Vec<u32>>,
        norm: u32,
        taxa_set: &HashSet<TreeNodeMeta<Self>>,
        upper_mixed: bool,
    ) -> TreeNodeZeta<Self> {
        let subtree_nodes = match upper_mixed {
            true => {
                t1_upper_postord_node_ids
            }
            false => t1_lower_postord_node_ids,
        };

        for n_id in t1_upper_postord_node_ids.iter().chain(t1_lower_postord_node_ids){
            t1_node_attributes.reset_node(*n_id);
        }

        match norm % 2 {
            0 => {
                // Preprocessing loop
                Self::preprocess_single_mix_even(
                    t1,
                    t2,
                    t2_lca_map,
                    t2_median,
                    t1_node_attributes,
                    t1_upper_taxa,
                    t1_lower_taxa,
                    norm,
                    taxa_set,
                    upper_mixed,
                );

                for v_id in subtree_nodes.iter() {
                    if v_id != &t1.get_root_id() && !t1.is_leaf(*v_id) {
                        // calculate v_sigma
                        let v_value = t1
                            .get_node_children_ids(*v_id)
                            .map(|x| t1_node_attributes[x].get_sigma())
                            .fold(
                                vec![<<<Self as RootedTree>::Node as RootedZetaNode>::Zeta>::zero();norm as usize + 1],
                                |acc, x| acc.iter().zip(x).map(|(a, b)| *a + *b).collect_vec(),
                            );
                        t1_node_attributes.set_sigma(*v_id, v_value);
                        let v_counterpart_count = t1
                            .get_node_children_ids(*v_id)
                            .map(|chid| t1_node_attributes[chid].get_counterpart_count())
                            .sum();
                        t1_node_attributes[*v_id].set_counterpart_count(v_counterpart_count);
                    }
                }

                for v_id in subtree_nodes.iter() {
                    if v_id != &t1.get_root_id() {
                        let v_sibling_id = t1.get_node_sibling_id(v_id);
                        let v_parent_id = t1.get_node_parent_id(*v_id).unwrap();
                        let v_counterpart_count = <TreeNodeZeta<Self> as NumCast>::from(t1_node_attributes[v_sibling_id].get_counterpart_count()).unwrap();
                        let summation_term = (0..norm+1).map(|l| {
                            let term1 = <TreeNodeZeta<Self> as NumCast>::from(pascal_triangle[norm as usize][l as usize]).unwrap();
                            let term2 = <TreeNodeZeta<Self> as NumCast>::from((-1_i32).pow(norm-l)).unwrap();
                            let term3 = (t1.get_zeta(v_parent_id).unwrap()).powi((norm-l) as i32);
                            let term4 = t1_node_attributes[*v_id].get_sigma()[l as usize];
                            term1*term2*term3*term4
                        }).sum::<TreeNodeZeta<Self>>();
                        t1_node_attributes[*v_id].set_kappa(v_counterpart_count * summation_term)
                    }
                }
            }
            _ => {
                Self::preprocess_single_mix_odd(
                    t1,
                    t2,
                    t2_lca_map,
                    t1_median,
                    t2_median,
                    t1_node_attributes,
                    t1_upper_taxa,
                    t1_lower_taxa,
                    norm,
                    taxa_set,
                    upper_mixed,
                );

                for v_id in subtree_nodes.iter() {
                    if v_id != &t1.get_root_id() && !t1.is_leaf(*v_id) {
                        let v_children = t1.get_node_children_ids(*v_id).collect_vec();
                        let v_left_child_id = v_children[0];
                        let v_right_child_id = v_children[1];
                        let v_delta = t1_node_attributes.get_delta(*v_id).to_vec();
                        // calculate v_sigma_pos
                        let v_left_sigma_pos = t1_node_attributes.get_sigma_pos(v_left_child_id);
                        let v_right_sigma_pos = t1_node_attributes.get_sigma_pos(v_right_child_id);
                        let v_sigma_pos = (0..norm + 1)
                            .map(|l| {
                                v_left_sigma_pos[l as usize] + v_right_sigma_pos[l as usize]
                                    - v_delta[l as usize]
                            })
                            .collect_vec();
                        t1_node_attributes.set_sigma_pos(*v_id, v_sigma_pos);
                        let v_left_sigma_neg = t1_node_attributes.get_sigma_neg(v_left_child_id);
                        let v_right_sigma_neg = t1_node_attributes.get_sigma_neg(v_right_child_id);
                        let v_sigma_neg = (0..norm + 1)
                            .map(|l| {
                                v_left_sigma_neg[l as usize]
                                    + v_right_sigma_neg[l as usize]
                                    + v_delta[l as usize]
                            })
                            .collect_vec();
                        t1_node_attributes.set_sigma_neg(*v_id, v_sigma_neg);
                        let v_counterpart_count = t1
                            .get_node_children_ids(*v_id)
                            .map(|chid| t1_node_attributes[chid].get_counterpart_count())
                            .sum();
                        t1_node_attributes.set_counterpart_count(*v_id, v_counterpart_count);
                    }
                }

                for v_id in subtree_nodes {
                    if v_id != &t1.get_root_id() {
                        let v_sibling_id = t1.get_node_sibling_id(v_id);
                        let v_parent_id = t1.get_node_parent_id(*v_id).unwrap();
                        let v_parent_zeta = t1.get_zeta(v_parent_id).unwrap();
                        let v_counterpart_count = <TreeNodeZeta<Self> as NumCast>::from(t1_node_attributes.get_counterpart_count(v_sibling_id)).unwrap();
                        let summation_term = (0..norm+1).map(|l| {
                            let term1 = <TreeNodeZeta<Self> as NumCast>::from(pascal_triangle[norm as usize][l as usize]).unwrap();
                            let term2 = <TreeNodeZeta<Self> as NumCast>::from((-1_i32).pow(norm-l)).unwrap();

                            let term3_1 = v_parent_zeta.powi(l as i32)*t1_node_attributes.get_sigma_pos(*v_id)[(norm-l) as usize];
                            let term3_2 = v_parent_zeta.powi((norm-l) as i32)*t1_node_attributes.get_sigma_neg(*v_id)[l as usize];

                            let term3 = term3_1+term3_2;

                            term1*term2*term3
                        }).sum::<TreeNodeZeta<Self>>();

                        t1_node_attributes[*v_id].set_kappa(v_counterpart_count * summation_term)
                    }
                }
            }
        }
        return t1
            .get_node_ids()
            .map(|x| *t1_node_attributes[x].get_kappa())
            .sum();
    }
}