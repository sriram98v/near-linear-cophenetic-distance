pub mod near_linear_cophenetic_distance;

use std::ops::{Index, IndexMut};

use fxhash::FxHashMap as HashMap;
use near_linear_cophenetic_distance::{
    NearLinearCopheneticDistance, NlcdAttribute, NlcdAttributeType, NlcdNodeAttributes,
    NlcdTreeAttributes,
};
use phylo::prelude::*;
use phylo::tree::SimpleRootedTree;

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
    attr: HashMap<usize, NodeAttributes>,
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

impl Index<<SimpleRootedTree as RootedTree>::NodeID> for TreeNodeAttributes {
    type Output = NodeAttributes;

    fn index(&self, index: <SimpleRootedTree as RootedTree>::NodeID) -> &Self::Output {
        self.attr.get(&index).unwrap()
        // &self.attr[index]
    }
}

impl IndexMut<<SimpleRootedTree as RootedTree>::NodeID> for TreeNodeAttributes {
    fn index_mut(&mut self, index: <SimpleRootedTree as RootedTree>::NodeID) -> &mut Self::Output {
        self.attr.get_mut(&index).unwrap()
        // &mut self.attr[index]
    }
}

impl NlcdTreeAttributes<<SimpleRootedTree as RootedTree>::NodeID, f32> for TreeNodeAttributes {}

impl NearLinearCopheneticDistance for SimpleRootedTree {
    type Meta = String;

    fn tree_attributes(
        &self,
        norm: u32,
    ) -> impl NlcdTreeAttributes<
        <Self as RootedTree>::NodeID,
        <<Self as RootedTree>::Node as phylo::prelude::RootedZetaNode>::Zeta,
    > {
        TreeNodeAttributes::new(self, norm)
    }
}
