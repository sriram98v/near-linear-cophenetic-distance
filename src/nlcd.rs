pub mod near_linear_cophenetic_distance;

use std::ops::{Index, IndexMut};

use fxhash::FxHashMap as HashMap;
use fxhash::FxHashSet as HashSet;
use itertools::Itertools;
use near_linear_cophenetic_distance::{
    NearLinearCopheneticDistance, NlcdAttribute, NlcdAttributeType, NlcdNodeAttributes,
    NlcdTreeAttributes,
};
use phylo::prelude::*;
use phylo::{tree::SimpleRootedTree, node::NodeID};

#[derive(Debug, Clone)]
pub struct NodeAttributes {
    sigma: NlcdAttribute<f32>,
    sigma_pos: NlcdAttribute<f32>,
    sigma_neg: NlcdAttribute<f32>,
    delta: NlcdAttribute<f32>,
    counterpart_count: NlcdAttribute<f32>,
    kappa: NlcdAttribute<f32>,
}

impl NodeAttributes {
    fn new(norm: u32) -> NodeAttributes {
        NodeAttributes {
            sigma: NlcdAttribute::Norms(vec![0_f32; (norm as usize) + 1]),
            sigma_pos: NlcdAttribute::Norms(vec![0_f32; (norm as usize) + 1]),
            sigma_neg: NlcdAttribute::Norms(vec![0_f32; (norm as usize) + 1]),
            delta: NlcdAttribute::Norms(vec![0_f32; (norm as usize) + 1]),
            counterpart_count: NlcdAttribute::Count(0_u32),
            kappa: NlcdAttribute::Kappa(0_f32),
        }
    }
}

impl Index<NlcdAttributeType> for NodeAttributes {
    type Output = NlcdAttribute<f32>;

    fn index(&self, index: NlcdAttributeType) -> &Self::Output {
        match index {
            NlcdAttributeType::Sigma => &self.sigma,
            NlcdAttributeType::SigmaPos => &self.sigma_pos,
            NlcdAttributeType::SigmaNeg => &self.sigma_neg,
            NlcdAttributeType::Delta => &self.delta,
            NlcdAttributeType::CounterpartCount => &self.counterpart_count,
            NlcdAttributeType::Kappa => &self.kappa,
        }
    }
}

impl IndexMut<NlcdAttributeType> for NodeAttributes {
    fn index_mut(&mut self, index: NlcdAttributeType) -> &mut Self::Output {
        match index {
            NlcdAttributeType::Sigma => &mut self.sigma,
            NlcdAttributeType::SigmaPos => &mut self.sigma_pos,
            NlcdAttributeType::SigmaNeg => &mut self.sigma_neg,
            NlcdAttributeType::Delta => &mut self.delta,
            NlcdAttributeType::CounterpartCount => &mut self.counterpart_count,
            NlcdAttributeType::Kappa => &mut self.kappa,
        }
    }
}

impl NlcdNodeAttributes<f32> for NodeAttributes {}

#[derive(Debug)]
pub struct TreeNodeAttributes {
    attr: HashMap<NodeID, NodeAttributes>,
}

impl TreeNodeAttributes {
    fn new(tree: &SimpleRootedTree, norm: u32) -> TreeNodeAttributes {
        TreeNodeAttributes {
            attr: tree.get_nodes().map(|x| (x.get_id(), NodeAttributes::new(norm))).collect(),
        }
        // let max_id = tree.get_node_ids().max().unwrap();
        // TreeNodeAttributes {
        //     attr: vec![NodeAttributes::new(norm); max_id + 1],
        // }
    }
}

impl Index<NodeID> for TreeNodeAttributes {
    type Output = NodeAttributes;

    fn index(&self, index: NodeID) -> &Self::Output {
        self.attr.get(&index).unwrap()
        // &self.attr[index]
    }
}

impl IndexMut<NodeID> for TreeNodeAttributes {
    fn index_mut(&mut self, index: NodeID) -> &mut Self::Output {
        self.attr.get_mut(&index).unwrap()
        // &mut self.attr[index]
    }
}

impl<'a> NlcdTreeAttributes<NodeID, f32> for TreeNodeAttributes {
}

impl<'a> NearLinearCopheneticDistance<'a> for SimpleRootedTree {

    fn tree_attributes(
        &self,
        norm: u32,
    ) -> impl NlcdTreeAttributes<usize,f32>
    {
        TreeNodeAttributes::new(self, norm)
    }

    /// Returns lower tree.
    fn lower_tree(&'a self, median_node_id: TreeNodeID<'a, Self>) -> Self {
        let lower_tree_taxa = HashSet::from_iter(self.lower_tree_taxa(median_node_id));
        self.contract_tree(
            &lower_tree_taxa
                .iter()
                .map(|x| self.get_taxa_node_id(x).unwrap())
                .collect_vec(),
        )
    }

    /// Returns upper tree.
    fn upper_tree(&'a self, median_node_id: TreeNodeID<'a, Self>) -> Self {
        let upper_tree_taxa = HashSet::from_iter(self.upper_tree_taxa(median_node_id));
        self.contract_tree(
            &upper_tree_taxa
                .iter()
                .map(|x| self.get_taxa_node_id(x).unwrap())
                .collect_vec(),
        )
    }

    // fn populate_op_vec(
    //     &'a self,
    //     tree: &'a Self,
    //     norm: u32,
    //     taxa_set: &HashSet<TreeNodeMeta<'a, Self>>,
    //     distances: &mut TreeNodeZeta<'a, Self>,
    //     self_node_attributes: &mut impl NlcdTreeAttributes<
    //         TreeNodeID<'a, Self>,
    //         TreeNodeZeta<'a, Self>,
    //     >,
    //     tree_node_attributes: &mut impl NlcdTreeAttributes<
    //         TreeNodeID<'a, Self>,
    //         TreeNodeZeta<'a, Self>,
    //     >,
    // ) {
    //     let t = self.get_median_node_id();
    //     let t_hat = tree.get_median_node_id();

    //     let self_postord_node_ids = self.postord_ids(self.get_root_id()).collect_vec();
    //     let tree_postord_node_ids = tree.postord_ids(tree.get_root_id()).collect_vec();

    //     let b: HashSet<TreeNodeMeta<'a, Self>> = HashSet::from_iter(
    //         self.get_cluster_ids(t)
    //             .filter_map(|x| self.get_node_taxa(x).cloned())
    //             .filter(|x| taxa_set.contains(x))
    //             // .map(|x| &x)
    //     );
    //     let b_hat: HashSet<TreeNodeMeta<'a, Self>> = HashSet::from_iter(
    //         tree.get_cluster_ids(t_hat)
    //             .filter_map(|x| tree.get_node_taxa(x).cloned())
    //             .filter(|x| taxa_set.contains(x))
    //             // .map(|x| &x)
    //     );

    //     let a: HashSet<TreeNodeMeta<'a, Self>> =
    //         HashSet::from_iter(self.get_taxa_space().map(|x| x.clone()))
    //             .difference(&b)
    //             .filter(|x| taxa_set.contains(x.clone()))
    //             .cloned()
    //             .collect();
    //     let a_hat: HashSet<TreeNodeMeta<'a, Self>> =
    //         HashSet::from_iter(tree.get_taxa_space().map(|x| x.clone()))
    //             .difference(&b_hat)
    //             .filter(|x| taxa_set.contains(x.clone()))
    //             .cloned()
    //             .collect();

    //     let a_int_a_hat: HashSet<TreeNodeMeta<'a, Self>> = a.intersection(&a_hat).map(|x| x.clone()).collect();
    //     let a_int_b_hat: HashSet<TreeNodeMeta<'a, Self>> = a.intersection(&b_hat).map(|x| x.clone()).collect();
    //     let b_int_a_hat: HashSet<TreeNodeMeta<'a, Self>> = b.intersection(&a_hat).map(|x| x.clone()).collect();
    //     let b_int_b_hat: HashSet<TreeNodeMeta<'a, Self>> = b.intersection(&b_hat).map(|x| x.clone()).collect();
        
    //     if taxa_set.len() > 2 {
    //         if a_int_a_hat.len() > 1 {
    //             let self_tree = self.contract_tree(
    //                 &a_int_a_hat
    //                     .iter()
    //                     .map(|x| self.get_taxa_node_id(x).unwrap())
    //                     .collect_vec(),
    //             ).clone();
    //             let new_tree = tree.contract_tree(
    //                 &a_int_a_hat
    //                     .iter()
    //                     .map(|x| tree.get_taxa_node_id(x).unwrap())
    //                     .collect_vec(),
    //             );
    //             self_tree.populate_op_vec(
    //                 &new_tree,
    //                 norm,
    //                 &a_int_a_hat,
    //                 distances,
    //                 self_node_attributes,
    //                 tree_node_attributes,
    //             );
    //         }

    //         if a_int_b_hat.len() > 1 {
    //             let self_tree = self.contract_tree(
    //                 &a_int_b_hat
    //                     .iter()
    //                     .map(|x| self.get_taxa_node_id(x).unwrap())
    //                     .collect_vec(),
    //             );
    //             let new_tree = tree.contract_tree(
    //                 &a_int_b_hat
    //                     .iter()
    //                     .map(|x| tree.get_taxa_node_id(x).unwrap())
    //                     .collect_vec(),
    //             );
    //             self_tree.populate_op_vec(
    //                 &new_tree,
    //                 norm,
    //                 &a_int_b_hat,
    //                 distances,
    //                 self_node_attributes,
    //                 tree_node_attributes,
    //             );
    //         }

    //         if b_int_b_hat.len() > 1 {
    //             let self_tree = self.contract_tree(
    //                 &b_int_b_hat
    //                     .iter()
    //                     .map(|x| self.get_taxa_node_id(x).unwrap())
    //                     .collect_vec(),
    //             );
    //             let new_tree = tree.contract_tree(
    //                 &b_int_b_hat
    //                     .iter()
    //                     .map(|x| tree.get_taxa_node_id(x).unwrap())
    //                     .collect_vec(),
    //             );
    //             self_tree.populate_op_vec(
    //                 &new_tree,
    //                 norm,
    //                 &b_int_b_hat,
    //                 distances,
    //                 self_node_attributes,
    //                 tree_node_attributes,
    //             );
    //         }

    //         if b_int_a_hat.len() > 1 {
    //             let self_tree = self.contract_tree(
    //                 &b_int_a_hat
    //                     .iter()
    //                     .map(|x| self.get_taxa_node_id(x).unwrap())
    //                     .collect_vec(),
    //             );
    //             let new_tree = tree.contract_tree(
    //                 &b_int_a_hat
    //                     .iter()
    //                     .map(|x| tree.get_taxa_node_id(x).unwrap())
    //                     .collect_vec(),
    //             );
    //             self_tree.populate_op_vec(
    //                 &new_tree,
    //                 norm,
    //                 &b_int_a_hat,
    //                 distances,
    //                 self_node_attributes,
    //                 tree_node_attributes,
    //             );
    //         }
    //     }
    //     let double_mix_distance = self.distance_double_mix_type(
    //         tree,
    //         norm,
    //         &t,
    //         &t_hat,
    //         &a_int_a_hat,
    //         &a_int_b_hat,
    //         &b_int_a_hat,
    //         &b_int_b_hat,
    //     );
    //     let single_mix_distance = self.distance_single_mix_type(
    //         tree,
    //         norm,
    //         &t,
    //         &t_hat,
    //         self_postord_node_ids.as_slice(),
    //         tree_postord_node_ids.as_slice(),
    //         self_node_attributes,
    //         tree_node_attributes,
    //         &a_int_a_hat,
    //         &a_int_b_hat,
    //         &b_int_a_hat,
    //     );

    //     *distances = *distances + double_mix_distance+single_mix_distance;
    // }
    
}
