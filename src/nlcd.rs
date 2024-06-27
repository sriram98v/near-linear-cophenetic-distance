pub mod near_linear_cophenetic_distance;

use std::ops::{Index, IndexMut};

use near_linear_cophenetic_distance::{NearLinearCopheneticDistance, NlcdNodeAttributes, NlcdTreeAttributes, NlcdAttributeType};
use phylo::prelude::*;
use phylo::tree::SimpleRootedTree;

#[derive(Debug, Clone)]
pub struct NodeAttributes{
    sigma: Vec<f32>,
    sigma_pos: Vec<f32>,
    sigma_neg: Vec<f32>,
    delta: Vec<f32>,
}

impl NodeAttributes{
    fn new(norm: u32)->NodeAttributes
    {
        NodeAttributes{
            sigma: vec![0_f32; (norm as usize)+1],
            sigma_pos: vec![0_f32; (norm as usize)+1],
            sigma_neg: vec![0_f32; (norm as usize)+1],
            delta: vec![0_f32; (norm as usize)+1],
        }
    }
}

impl Index<NlcdAttributeType> for NodeAttributes{
    type Output = Vec<f32>;

    fn index(&self, index: NlcdAttributeType) -> &Self::Output {
        match index{
            NlcdAttributeType::Sigma => {&self.sigma}
            NlcdAttributeType::SigmaPos => {&self.sigma_pos}
            NlcdAttributeType::SigmaNeg => {&self.sigma_neg}
            NlcdAttributeType::Delta => {&self.delta}
        }
    }
}

impl IndexMut<NlcdAttributeType> for NodeAttributes{
    fn index_mut(&mut self, index: NlcdAttributeType) -> &mut Self::Output {
        match index{
            NlcdAttributeType::Sigma => {&mut self.sigma}
            NlcdAttributeType::SigmaPos => {&mut self.sigma_pos}
            NlcdAttributeType::SigmaNeg => {&mut self.sigma_neg}
            NlcdAttributeType::Delta => {&mut self.delta}
        }
    }
}

impl NlcdNodeAttributes<f32> for NodeAttributes{
    fn reset(&mut self) {
        for i in 0..self.sigma.len(){
            self.sigma[i] = 0_f32;
        }
        for i in 0..self.sigma_pos.len(){
            self.sigma_pos[i] = 0_f32;
        }
        for i in 0..self.sigma_neg.len(){
            self.sigma_neg[i] = 0_f32;
        }
        for i in 0..self.delta.len(){
            self.delta[i] = 0_f32;
        }
    }
}

pub struct TreeNodeAttributes{
    attr: Vec<NodeAttributes>,
}

impl TreeNodeAttributes{
    fn new(tree: &SimpleRootedTree, norm: u32)-> TreeNodeAttributes
    {
        let max_id = tree.get_node_ids().max().unwrap();
        TreeNodeAttributes{
            attr: vec![NodeAttributes::new(norm);max_id+1]
        }
    }
}

impl Index<<SimpleRootedTree as RootedTree>::NodeID> for TreeNodeAttributes{
    type Output = NodeAttributes;

    fn index(&self, index: <SimpleRootedTree as RootedTree>::NodeID) -> &Self::Output {
        &self.attr[index]
    }
}

impl IndexMut<<SimpleRootedTree as RootedTree>::NodeID> for TreeNodeAttributes{
    fn index_mut(&mut self, index: <SimpleRootedTree as RootedTree>::NodeID) -> &mut Self::Output {
        &mut self.attr[index]
    }
}

impl NlcdTreeAttributes<<SimpleRootedTree as RootedTree>::NodeID, f32> for TreeNodeAttributes{}



impl NearLinearCopheneticDistance for SimpleRootedTree{
    type Meta = String;

    fn tree_attributes(&self, norm: u32)->impl NlcdTreeAttributes<<Self as RootedTree>::NodeID, <<Self as RootedTree>::Node as phylo::prelude::RootedZetaNode>::Zeta> {
        TreeNodeAttributes::new(self, norm)
    }

}