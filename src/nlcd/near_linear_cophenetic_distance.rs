use fxhash::FxHashMap as HashMap;
use fxhash::FxHashSet as HashSet;
use std::{
    fmt::{Debug, Display},
    hash::Hash,
    ops::{AddAssign, IndexMut, Index},
};

use itertools::Itertools;
use num::{Float, NumCast, Signed, Zero};
use phylo::prelude::*;

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
    <Self as RootedTree>::Node:
        RootedMetaNode<Meta = <Self as CopheneticDistance>::Meta> + RootedZetaNode,
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
    type Meta: Display + Debug + Eq + PartialEq + Clone + Ord + Hash + Send + Sync;

    // helper functions
    fn tree_attributes(
        &self,
        norm: u32,
    ) -> impl NlcdTreeAttributes<
        <Self as RootedTree>::NodeID,
        <<Self as RootedTree>::Node as RootedZetaNode>::Zeta,
    >;

    /// Returns taxa present in upper tree.
    fn upper_tree_taxa(&self, median_node_id: <Self as RootedTree>::NodeID) -> impl Iterator<Item = <Self as CopheneticDistance>::Meta> {
        let lower_tree_taxa = self.lower_tree_taxa(median_node_id).collect::<HashSet<_>>();
        self.get_leaves()
            .map(|x| x.get_taxa().unwrap())
            .filter(move |x| !lower_tree_taxa.contains(x))
    }

    /// Returns taxa present in lower tree.
    fn lower_tree_taxa(&self, median_node_id: <Self as RootedTree>::NodeID) -> impl Iterator<Item = <Self as CopheneticDistance>::Meta> {
        self.get_cluster(median_node_id)
            .filter_map(|x| x.get_taxa())
    }

    /// Returns taxa present in upper tree.
    fn upper_tree_leaves(&self, median_node_id: <Self as RootedTree>::NodeID) -> impl Iterator<Item = <Self as RootedTree>::Node> {
        let lower_tree_leaf_ids = self
            .lower_tree_leaves(median_node_id)
            .map(|x| x.get_id())
            .collect::<HashSet<_>>();
        self.get_leaves()
            .filter(move |x| !lower_tree_leaf_ids.contains(&x.get_id()))
            .cloned()
    }

    /// Returns leaves present in lower tree.
    fn lower_tree_leaves(&self, median_node_id: <Self as RootedTree>::NodeID) -> impl Iterator<Item = <Self as RootedTree>::Node> {
        self.get_cluster(median_node_id)
    }

    /// Returns lower tree.
    fn lower_tree(&self, median_node_id: <Self as RootedTree>::NodeID) -> Self {
        let lower_tree_taxa = HashSet::from_iter(self.lower_tree_taxa(median_node_id));
        self.contract_tree(
            &lower_tree_taxa
                .iter()
                .map(|x| self.get_taxa_node_id(x).unwrap())
                .collect_vec(),
        )
    }

    /// Returns upper tree.
    fn upper_tree(&self, median_node_id: <Self as RootedTree>::NodeID) -> Self {
        let upper_tree_taxa = HashSet::from_iter(self.upper_tree_taxa(median_node_id));
        self.contract_tree(
            &upper_tree_taxa
                .iter()
                .map(|x| self.get_taxa_node_id(x).unwrap())
                .collect_vec(),
        )
    }

    fn get_node_sibling_id(
        &self,
        node_id: &<Self as RootedTree>::NodeID,
    ) -> <Self as RootedTree>::NodeID {
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
            false => multiplicative_form(n - k, n - k),
        }
    }

    /// Returns the Cophenetic distance between two trees in O(pnlog^2n) time.
    fn cophen_dist(
        &self,
        tree: &Self,
        norm: u32,
    ) -> <<Self as RootedTree>::Node as RootedZetaNode>::Zeta {
        // let mut ops: Vec<<<Self as RootedTree>::Node as RootedZetaNode>::Zeta> = vec![];

        let mut distances: <<Self as RootedTree>::Node as RootedZetaNode>::Zeta = <<<Self as RootedTree>::Node as RootedZetaNode>::Zeta>::zero();

        let binding1 = self
            .get_taxa_space()
            .collect::<HashSet<<Self as CopheneticDistance>::Meta>>();
        let binding2 = tree
            .get_taxa_space()
            .collect::<HashSet<<Self as CopheneticDistance>::Meta>>();
        let taxa_set: HashSet<_> = binding1.intersection(&binding2).collect();

        let mut self_node_attributes = self.tree_attributes(norm);
        let mut tree_node_attributes = tree.tree_attributes(norm);

        self.populate_op_vec(
            tree,
            norm,
            &taxa_set,
            &mut distances,
            &mut self_node_attributes,
            &mut tree_node_attributes,
        );

        let taxa_distance = taxa_set
            .iter()
            .map(|x| {
                let zeta_1 = self.get_zeta_taxa(x);
                let zeta_2 = tree.get_zeta_taxa(x);
                (zeta_1 - zeta_2).abs().powi(norm as i32)
            })
            .sum();
        (distances + taxa_distance).powf(
            <<<Self as RootedTree>::Node as RootedZetaNode>::Zeta as NumCast>::from(norm)
                .unwrap()
                .powi(-1),
        )
    }

    /// Populates vector with divide and conquer distances.
    fn populate_op_vec(
        &self,
        tree: &Self,
        norm: u32,
        taxa_set: &HashSet<&<Self as RootedMetaTree>::Meta>,
        distances: &mut <<Self as RootedTree>::Node as RootedZetaNode>::Zeta,
        self_node_attributes: &mut impl NlcdTreeAttributes<
            <Self as RootedTree>::NodeID,
            <<Self as RootedTree>::Node as RootedZetaNode>::Zeta,
        >,
        tree_node_attributes: &mut impl NlcdTreeAttributes<
            <Self as RootedTree>::NodeID,
            <<Self as RootedTree>::Node as RootedZetaNode>::Zeta,
        >,
    ) {
        let t = self.get_median_node_id();
        let t_hat = tree.get_median_node_id();

        let b: HashSet<<Self as CopheneticDistance>::Meta> = HashSet::from_iter(
            self.get_cluster(t)
                .filter_map(|x| x.get_taxa())
                .filter(|x| taxa_set.contains(x))
                // .map(|x| &x)
        );
        let b_hat: HashSet<<Self as CopheneticDistance>::Meta> = HashSet::from_iter(
            tree.get_cluster(t_hat)
                .filter_map(|x| x.get_taxa())
                .filter(|x| taxa_set.contains(x))
                // .map(|x| &x)
        );

        let a: HashSet<<Self as CopheneticDistance>::Meta> =
            HashSet::from_iter(self.get_taxa_space())
                .difference(&b)
                .filter(|x| taxa_set.contains(x))
                .cloned()
                .collect();
        let a_hat: HashSet<<Self as CopheneticDistance>::Meta> =
            HashSet::from_iter(tree.get_taxa_space())
                .difference(&b_hat)
                .filter(|x| taxa_set.contains(x))
                .cloned()
                .collect();

        let a_int_a_hat: HashSet<&<Self as CopheneticDistance>::Meta> = a.intersection(&a_hat).collect();
        let a_int_b_hat: HashSet<&<Self as CopheneticDistance>::Meta> = a.intersection(&b_hat).collect();
        let b_int_a_hat: HashSet<&<Self as CopheneticDistance>::Meta> = b.intersection(&a_hat).collect();
        let b_int_b_hat: HashSet<&<Self as CopheneticDistance>::Meta> = b.intersection(&b_hat).collect();

        let double_mix_distance = self.distance_double_mix_type(
            tree,
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
            norm,
            &t,
            &t_hat,
            self_node_attributes,
            tree_node_attributes,
            &a_int_a_hat,
            &a_int_b_hat,
            &b_int_a_hat,
        );

        *distances = *distances + double_mix_distance+single_mix_distance;

        if taxa_set.len() > 2 {
            if a_int_a_hat.len() > 1 {
                let self_tree = self.contract_tree(
                    &a_int_a_hat
                        .iter()
                        .map(|x| self.get_taxa_node_id(x).unwrap())
                        .collect_vec(),
                );
                let new_tree = tree.contract_tree(
                    &a_int_a_hat
                        .iter()
                        .map(|x| tree.get_taxa_node_id(x).unwrap())
                        .collect_vec(),
                );
                self_tree.populate_op_vec(
                    &new_tree,
                    norm,
                    &a_int_a_hat,
                    distances,
                    self_node_attributes,
                    tree_node_attributes,
                );
            }

            if a_int_b_hat.len() > 1 {
                let self_tree = self.contract_tree(
                    &a_int_b_hat
                        .iter()
                        .map(|x| self.get_taxa_node_id(x).unwrap())
                        .collect_vec(),
                );
                let new_tree = tree.contract_tree(
                    &a_int_b_hat
                        .iter()
                        .map(|x| tree.get_taxa_node_id(x).unwrap())
                        .collect_vec(),
                );
                self_tree.populate_op_vec(
                    &new_tree,
                    norm,
                    &a_int_b_hat,
                    distances,
                    self_node_attributes,
                    tree_node_attributes,
                );
            }

            if b_int_b_hat.len() > 1 {
                let self_tree = self.contract_tree(
                    &b_int_b_hat
                        .iter()
                        .map(|x| self.get_taxa_node_id(x).unwrap())
                        .collect_vec(),
                );
                let new_tree = tree.contract_tree(
                    &b_int_b_hat
                        .iter()
                        .map(|x| tree.get_taxa_node_id(x).unwrap())
                        .collect_vec(),
                );
                self_tree.populate_op_vec(
                    &new_tree,
                    norm,
                    &b_int_b_hat,
                    distances,
                    self_node_attributes,
                    tree_node_attributes,
                );
            }

            if b_int_a_hat.len() > 1 {
                let self_tree = self.contract_tree(
                    &b_int_a_hat
                        .iter()
                        .map(|x| self.get_taxa_node_id(x).unwrap())
                        .collect_vec(),
                );
                let new_tree = tree.contract_tree(
                    &b_int_a_hat
                        .iter()
                        .map(|x| tree.get_taxa_node_id(x).unwrap())
                        .collect_vec(),
                );
                self_tree.populate_op_vec(
                    &new_tree,
                    norm,
                    &b_int_a_hat,
                    distances,
                    self_node_attributes,
                    tree_node_attributes,
                );
            }
        }
    }

    /// Returns ordered iterator used in double mix type cases
    fn get_cntr(
        &self,
        median_node_id: <Self as RootedTree>::NodeID,
        leaf_set: HashSet<Self::NodeID>,
    ) -> Vec<<<Self as RootedTree>::Node as RootedZetaNode>::Zeta> {
        // line 5 in algo 1
        let mut gamma: Vec<<<Self as RootedTree>::Node as RootedZetaNode>::Zeta> = Vec::new();
        // line 3 in algo 1
        let mut median_path = self
            .root_to_node(median_node_id)
            .into_iter()
            .map(|x| (x.get_id(), 0))
            .collect::<HashMap<_, _>>();
        for node_id in leaf_set {
            // line 4 in algo 1
            median_path
                .entry(self.get_lca_id(&vec![node_id, median_node_id]))
                .and_modify(|x| *x += 1);
        }
        for node in self.root_to_node(median_node_id).into_iter() {
            let c = median_path.get(&node.get_id()).cloned().unwrap();
            for _ in 0..c {
                gamma.push(self.get_zeta(node.get_id()).unwrap())
            }
        }
        gamma
    }

    /// Returns seqPrd^p(\alpha,\beta) used for double mix type cases.
    fn seq_product(
        mut alpha: Vec<<<Self as RootedTree>::Node as RootedZetaNode>::Zeta>,
        mut beta: Vec<<<Self as RootedTree>::Node as RootedZetaNode>::Zeta>,
        norm: u32,
    ) -> <<Self as RootedTree>::Node as RootedZetaNode>::Zeta {
        if alpha.is_empty() || beta.is_empty() {
            return <<<Self as RootedTree>::Node as RootedZetaNode>::Zeta as Zero>::zero();
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
                let t1 = <<<Self as RootedTree>::Node as RootedZetaNode>::Zeta as NumCast>::from(Self::n_choose_k(norm, l)).unwrap();
                let t2 = beta[j].powi(l as i32);
                let t3 = <<<Self as RootedTree>::Node as RootedZetaNode>::Zeta as NumCast>::from((-1_i32).pow(norm-l)).unwrap();
                let t4 = sigma[j][(norm-l) as usize];
                t1*t2*t3*t4
            }).sum::<<<Self as RootedTree>::Node as RootedZetaNode>::Zeta>();

            let term_2 = (0..norm+1).map(|l| {
                let t1 = <<<Self as RootedTree>::Node as RootedZetaNode>::Zeta as NumCast>::from(Self::n_choose_k(norm, l)).unwrap();
                let t2 = beta[j].powi((norm-l) as i32);
                let t3 = <<<Self as RootedTree>::Node as RootedZetaNode>::Zeta as NumCast>::from((-1_i32).pow(norm-l)).unwrap();
                let t4 = alpha.iter().map(|a| a.powi(l as i32)).sum::<<<Self as RootedTree>::Node as RootedZetaNode>::Zeta>()-sigma[j][l as usize];
                t1*t2*t3*t4
            }).sum::<<<Self as RootedTree>::Node as RootedZetaNode>::Zeta>();

            term_1+term_2
        }).sum::<<<Self as RootedTree>::Node as RootedZetaNode>::Zeta>();

        final_out
    }

    /// Algorithm 3 for precomputation of alpha sums
    fn precompute_alpha_sums(
        alpha: &[<<Self as RootedTree>::Node as RootedZetaNode>::Zeta],
        beta: &[<<Self as RootedTree>::Node as RootedZetaNode>::Zeta],
        norm: u32,
    ) -> Vec<Vec<<<Self as RootedTree>::Node as RootedZetaNode>::Zeta>> {
        let mut j = 0;
        let mut i = 0;
        let k = alpha.len();
        let m = beta.len();

        let mut sigma = vec![
            vec![
                <<Self as RootedTree>::Node as RootedZetaNode>::Zeta::zero();
                (norm + 1) as usize
            ];
            m
        ];

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
        norm: u32,
        t: &<Self as RootedTree>::NodeID,
        t_hat: &<Self as RootedTree>::NodeID,
        a_int_a_hat: &HashSet<&<Self as CopheneticDistance>::Meta>,
        a_int_b_hat: &HashSet<&<Self as CopheneticDistance>::Meta>,
        b_int_a_hat: &HashSet<&<Self as CopheneticDistance>::Meta>,
        b_int_b_hat: &HashSet<&<Self as CopheneticDistance>::Meta>,
    ) -> <<Self as RootedTree>::Node as RootedZetaNode>::Zeta {
        let alpha = self.get_cntr(
            *t,
            a_int_b_hat
                .iter()
                .map(|x| self.get_taxa_node_id(x).unwrap())
                .collect::<HashSet<Self::NodeID>>(),
        );
        let beta = tree.get_cntr(
                *t_hat,
            b_int_a_hat
                .iter()
                .map(|x| tree.get_taxa_node_id(x).unwrap())
                .collect::<HashSet<Self::NodeID>>(),
        );

        // AB|A'B'
        let b_int_b_hat_len = b_int_b_hat.len();
        let dd2 = a_int_a_hat
            .iter()
            .map(|x| {
                let t_lca_id = self.get_lca_id(&vec![self.get_taxa_node_id(x).unwrap(), *t]);
                let t_hat_lca_id =
                    tree.get_lca_id(&vec![tree.get_taxa_node_id(x).unwrap(), *t_hat]);
                let zeta_1 = self.get_zeta(t_lca_id).unwrap();
                let zeta_2 = tree.get_zeta(t_hat_lca_id).unwrap();
                (zeta_1 - zeta_2).abs().powi(norm as i32)
            })
            .sum::<<<Self as RootedTree>::Node as RootedZetaNode>::Zeta>()
            * <<<Self as RootedTree>::Node as RootedZetaNode>::Zeta as NumCast>::from(
                b_int_b_hat_len,
            )
            .unwrap();

        Self::seq_product(alpha, beta, norm) + dd2
    }

    /// This method generates the distance contributed by all taxa pairs
    /// that are present in the same subtree in exactly one of the two trees(raised to the p^{th} power).
    ///
    /// This includes the following assignments: AA|A'B', AA|B'A', BB|A'B', BB|B'A', BA|B'B', BA|A'A', AB|B'B', AB|A'A'.
    #[allow(clippy::too_many_arguments)]
    fn distance_single_mix_type(
        &self,
        tree: &Self,
        norm: u32,
        t: &<Self as RootedTree>::NodeID,
        t_hat: &<Self as RootedTree>::NodeID,
        self_node_attributes: &mut impl NlcdTreeAttributes<
            <Self as RootedTree>::NodeID,
            <<Self as RootedTree>::Node as RootedZetaNode>::Zeta,
        >,
        tree_node_attributes: &mut impl NlcdTreeAttributes<
            <Self as RootedTree>::NodeID,
            <<Self as RootedTree>::Node as RootedZetaNode>::Zeta,
        >,
        a_int_a_hat: &HashSet<&<Self as CopheneticDistance>::Meta>,
        a_int_b_hat: &HashSet<&<Self as CopheneticDistance>::Meta>,
        b_int_a_hat: &HashSet<&<Self as CopheneticDistance>::Meta>,
    ) -> <<Self as RootedTree>::Node as RootedZetaNode>::Zeta {
        if self.num_taxa() <= 2 {
            return <<<Self as RootedTree>::Node as RootedZetaNode>::Zeta>::zero();
        }
        // AA|A'B'
        let d1 = Self::single_mix_xxxy(
            self,
            tree,
            t,
            t_hat,
            self_node_attributes,
            norm,
            a_int_a_hat,
            true,
        );
        // AB|A'A'
        let d2 = Self::single_mix_xxxy(
            tree,
            self,
            t_hat,
            t,
            tree_node_attributes,
            norm,
            a_int_a_hat,
            true,
        );
        // BB|A'B'
        let d3 = Self::single_mix_xxxy(
            self,
            tree,
            t,
            t_hat,
            self_node_attributes,
            norm,
            b_int_a_hat,
            false,
        );
        // AB|B'B'
        let d4 = Self::single_mix_xxxy(
            tree,
            self,
            t_hat,
            t,
            tree_node_attributes,
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
        t1_median: &<Self as RootedTree>::NodeID,
        t2_median: &<Self as RootedTree>::NodeID,
        t1_node_attributes: &mut impl NlcdTreeAttributes<
            <Self as RootedTree>::NodeID,
            <<Self as RootedTree>::Node as RootedZetaNode>::Zeta,
        >,
        norm: u32,
        taxa_set: &HashSet<&<Self as CopheneticDistance>::Meta>,
        upper_mixed: bool,
    ) {
        let leaf_iter = match upper_mixed {
            true => t1.upper_tree_taxa(*t1_median).collect_vec(),
            false => t1.lower_tree_taxa(*t1_median).collect_vec(),
        };
        for leaf in leaf_iter {
            match taxa_set.contains(&leaf) {
                true => {
                    let t1_node_id = t1.get_taxa_node_id(&leaf).unwrap();
                    let t2_node_id = t2.get_taxa_node_id(&leaf).unwrap();
                    let lca_x_t_hat = t2.get_lca_id(&vec![t2_node_id, *t2_median]);
                    let beta = t2.get_zeta(lca_x_t_hat).unwrap();
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

    #[allow(clippy::too_many_arguments)]
    fn preprocess_single_mix_odd(
        t1: &Self,
        t2: &Self,
        t1_median: &<Self as RootedTree>::NodeID,
        t2_median: &<Self as RootedTree>::NodeID,
        t1_node_attributes: &mut impl NlcdTreeAttributes<
            <Self as RootedTree>::NodeID,
            <<Self as RootedTree>::Node as RootedZetaNode>::Zeta,
        >,
        norm: u32,
        taxa_set: &HashSet<&<Self as CopheneticDistance>::Meta>,
        upper_mixed: bool,
    ) {
        let leaf_iter = match upper_mixed {
            true => t1.upper_tree_taxa(*t1_median).collect_vec(),
            false => t1.lower_tree_taxa(*t1_median).collect_vec(),
        };
        for leaf in leaf_iter.iter() {
            match taxa_set.contains(&leaf) {
                true => {
                    let t2_node_id = t2.get_taxa_node_id(leaf).unwrap();
                    let lca_x_t_hat = t2.get_lca_id(&vec![t2_node_id, *t2_median]);
                    let beta = t2.get_zeta(lca_x_t_hat).unwrap();
                    let t1_node_id = t1.get_taxa_node_id(leaf).unwrap();

                    let t1_node_parent_id = t1.get_node_parent_id(t1_node_id).unwrap();

                    if beta <= t1.get_zeta(t1_node_id).unwrap() {
                        //find omega_x
                        let mut v_ancestors = vec![];
                        for node in t1.node_to_root(t1_node_id).into_iter() {
                            v_ancestors.push(node.get_id());
                            if &node.get_id() == t1_median {
                                break;
                            }
                        }

                        let omega_x = v_ancestors
                            .into_iter()
                            .min_by(|w, y| {
                                t1.get_zeta(*w)
                                    .unwrap()
                                    .partial_cmp(&t1.get_zeta(*y).unwrap())
                                    .unwrap()
                            })
                            .unwrap();
                        // set omega_x.delta
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
        t1_median: &<Self as RootedTree>::NodeID,
        t2_median: &<Self as RootedTree>::NodeID,
        t1_node_attributes: &mut impl NlcdTreeAttributes<
            <Self as RootedTree>::NodeID,
            <<Self as RootedTree>::Node as RootedZetaNode>::Zeta,
        >,
        norm: u32,
        taxa_set: &HashSet<&<Self as CopheneticDistance>::Meta>,
        upper_mixed: bool,
    ) -> <<Self as RootedTree>::Node as RootedZetaNode>::Zeta {
        let subtree_nodes = match upper_mixed {
            true => {
                let lower_tree_nodes = t1.postord(*t1_median).map(|x| x.get_id()).collect_vec();
                t1.postord(t1.get_root_id())
                    .map(|x| x.get_id())
                    .filter(|x| !lower_tree_nodes.contains(x))
                    .collect_vec()
            }
            false => t1.postord(*t1_median).map(|x| x.get_id()).collect_vec(),
        };

        // resetting all node attributes
        t1.get_node_ids()
            .for_each(|node_id| t1_node_attributes.reset_node(node_id));

        match norm % 2 {
            0 => {
                // Preprocessing loop
                Self::preprocess_single_mix_even(
                    t1,
                    t2,
                    t1_median,
                    t2_median,
                    t1_node_attributes,
                    norm,
                    taxa_set,
                    upper_mixed,
                );

                for v_id in subtree_nodes.iter() {
                    if v_id != &t1.get_root_id() && !t1.is_leaf(v_id) {
                        // calculate v_sigma
                        let v_value = t1
                            .get_node_children_ids(*v_id)
                            .map(|x| t1_node_attributes[x].get_sigma())
                            .fold(
                                vec![
                                    <<Self as RootedTree>::Node as RootedZetaNode>::Zeta::zero();
                                    norm as usize + 1
                                ],
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
                        let v_counterpart_count = <<<Self as RootedTree>::Node as RootedZetaNode>::Zeta as NumCast>::from(t1_node_attributes[v_sibling_id].get_counterpart_count()).unwrap();
                        let summation_term = (0..norm+1).map(|l| {
                            let term1 = <<<Self as RootedTree>::Node as RootedZetaNode>::Zeta as NumCast>::from(Self::n_choose_k(norm, l)).unwrap();
                            let term2 = <<<Self as RootedTree>::Node as RootedZetaNode>::Zeta as NumCast>::from((-1_i32).pow(norm-l)).unwrap();
                            let term3 = (t1.get_zeta(v_parent_id).unwrap()).powi((norm-l) as i32);
                            let term4 = t1_node_attributes[*v_id].get_sigma()[l as usize];
                            term1*term2*term3*term4
                        }).sum::<<<Self as RootedTree>::Node as RootedZetaNode>::Zeta>();
                        t1_node_attributes[*v_id].set_kappa(v_counterpart_count * summation_term)
                    }
                }
            }
            _ => {
                Self::preprocess_single_mix_odd(
                    t1,
                    t2,
                    t1_median,
                    t2_median,
                    t1_node_attributes,
                    norm,
                    taxa_set,
                    upper_mixed,
                );

                for v_id in subtree_nodes.iter() {
                    if v_id != &t1.get_root_id() && !t1.is_leaf(v_id) {
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
                    if v_id != t1.get_root_id() {
                        let v_sibling_id = t1.get_node_sibling_id(&v_id);
                        let v_parent_id = t1.get_node_parent_id(v_id).unwrap();
                        let v_parent_zeta = t1.get_zeta(v_parent_id).unwrap();
                        let v_counterpart_count = <<<Self as RootedTree>::Node as RootedZetaNode>::Zeta as NumCast>::from(t1_node_attributes.get_counterpart_count(v_sibling_id)).unwrap();
                        let summation_term = (0..norm+1).map(|l| {
                            let term1 = <<<Self as RootedTree>::Node as RootedZetaNode>::Zeta as NumCast>::from(Self::n_choose_k(norm, l)).unwrap();
                            let term2 = <<<Self as RootedTree>::Node as RootedZetaNode>::Zeta as NumCast>::from((-1_i32).pow(norm-l)).unwrap();

                            let term3_1 = v_parent_zeta.powi(l as i32)*t1_node_attributes.get_sigma_pos(v_id)[(norm-l) as usize];
                            let term3_2 = v_parent_zeta.powi((norm-l) as i32)*t1_node_attributes.get_sigma_neg(v_id)[l as usize];

                            let term3 = term3_1+term3_2;

                            term1*term2*term3
                        }).sum::<<<Self as RootedTree>::Node as RootedZetaNode>::Zeta>();
                        t1_node_attributes[v_id].set_kappa(v_counterpart_count * summation_term)
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