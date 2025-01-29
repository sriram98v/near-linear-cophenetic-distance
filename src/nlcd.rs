pub mod near_linear_cophenetic_distance;

use std::ops::{Index, IndexMut};
use std::usize;

use fxhash::FxHashMap as HashMap;
use fxhash::FxHashSet as HashSet;
use itertools::Itertools;
use near_linear_cophenetic_distance::{
    NearLinearCopheneticDistance, NlcdAttribute, NlcdAttributeType, NlcdNodeAttributes,
    NlcdTreeAttributes,
};
use phylo::prelude::*;
use vers_vecs::BitVec;

#[derive(Debug, Clone)]
pub struct NodeAttributes<Z: NodeWeight> {
    sigma: NlcdAttribute<Z>,
    sigma_pos: NlcdAttribute<Z>,
    sigma_neg: NlcdAttribute<Z>,
    delta: NlcdAttribute<Z>,
    counterpart_count: NlcdAttribute<Z>,
    kappa: NlcdAttribute<Z>,
}

impl<Z: NodeWeight> NodeAttributes<Z> {
    fn new(norm: u32) -> NodeAttributes<Z> {
        NodeAttributes {
            sigma: NlcdAttribute::Norms(vec![Z::zero(); (norm as usize) + 1]),
            sigma_pos: NlcdAttribute::Norms(vec![Z::zero(); (norm as usize) + 1]),
            sigma_neg: NlcdAttribute::Norms(vec![Z::zero(); (norm as usize) + 1]),
            delta: NlcdAttribute::Norms(vec![Z::zero(); (norm as usize) + 1]),
            counterpart_count: NlcdAttribute::Count(0_u32),
            kappa: NlcdAttribute::Kappa(Z::zero()),
        }
    }
}

impl<Z: NodeWeight> Index<NlcdAttributeType> for NodeAttributes<Z> {
    type Output = NlcdAttribute<Z>;

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

impl<Z: NodeWeight> IndexMut<NlcdAttributeType> for NodeAttributes<Z> {
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

impl<Z: NodeWeight> NlcdNodeAttributes<Z> for NodeAttributes<Z> {}

#[derive(Debug)]
pub struct TreeNodeAttributes<T: RootedTree, Z: NodeWeight> {
    attr: HashMap<TreeNodeID<T>, NodeAttributes<Z>>,
}

impl<T, Z> TreeNodeAttributes<T, Z>
where 
    T: RootedTree,
    Z: NodeWeight
{
    fn new<'a>(tree: &'a T, norm: u32) -> TreeNodeAttributes<T, Z> {
        TreeNodeAttributes {
            attr: tree.get_node_ids().map(|x| (x, NodeAttributes::new(norm))).collect(),
        }
        // let max_id = tree.get_node_ids().max().unwrap();
        // TreeNodeAttributes {
        //     attr: vec![NodeAttributes::new(norm); max_id + 1],
        // }
    }
}

impl<T: RootedTree, Z:NodeWeight> Index<TreeNodeID<T>> for TreeNodeAttributes<T, Z> {
    type Output = NodeAttributes<Z>;

    fn index(&self, index: TreeNodeID<T>) -> &Self::Output {
        self.attr.get(&index).unwrap()
        // &self.attr[index]
    }
}

impl<T: RootedTree, Z: NodeWeight> IndexMut<TreeNodeID<T>> for TreeNodeAttributes<T, Z> {
    fn index_mut(&mut self, index: TreeNodeID<T>) -> &mut Self::Output {
        self.attr.get_mut(&index).unwrap()
        // &mut self.attr[index]
    }
}

impl<'a, T: RootedTree, Z: NodeWeight> NlcdTreeAttributes<TreeNodeID<T>, Z> for TreeNodeAttributes<T, Z> {
}

impl<T: NodeTaxa, W: EdgeWeight, Z: NodeWeight> NearLinearCopheneticDistance for SimpleRootedTree<T,W,Z> {

    fn tree_attributes(
        &self,
        norm: u32,
    ) -> impl NlcdTreeAttributes<usize,Z>
    {
        TreeNodeAttributes::new(self, norm)
    }

    /// Returns lower tree.
    fn lower_tree(&self, median_node_id: TreeNodeID<Self>) -> Self {
        let lower_tree_taxa = HashSet::from_iter(self.lower_tree_taxa(median_node_id));
        self.contract_tree(
            &lower_tree_taxa
                .iter()
                .map(|x| self.get_taxa_node_id(x).unwrap())
                .collect_vec(),
        ).unwrap()
    }

    /// Returns upper tree.
    fn upper_tree(&self, median_node_id: TreeNodeID<Self>) -> Self {
        let upper_tree_taxa = HashSet::from_iter(self.upper_tree_taxa(median_node_id));
        self.contract_tree(
            &upper_tree_taxa
                .iter()
                .map(|x| self.get_taxa_node_id(x).unwrap())
                .collect_vec(),
        ).unwrap()
    }

    fn contract_and_median(tree: &Self, taxa_map: &HashMap<&TreeNodeMeta<Self>, usize>, tree_postord_node_ids: &[TreeNodeID<Self>])->(Self, TreeNodeID<Self>, HashMap<TreeNodeID<Self>, BitVec>, Vec<TreeNodeID<Self>>, Vec<TreeNodeID<Self>>, Vec<TreeNodeID<Self>>) {
        let mut leaf_ids = Vec::with_capacity(taxa_map.len());

        for x in taxa_map.keys(){
            leaf_ids.push(tree.get_taxa_node_id(*x).unwrap());
        }

        let new_tree_root_id = tree.get_lca_id(leaf_ids.as_slice());
        let new_tree_root_parent_id = tree.get_node_parent_id(new_tree_root_id);
        let mut node_map: Vec<Option<Self::Node>> = vec![None; tree.get_capacity()];
        node_map[new_tree_root_id] = Some(tree.get_lca(leaf_ids.as_slice()).clone());

        let mut leaf_id_set = vec![false; tree.get_capacity()];
        for id in leaf_ids.iter() {
            leaf_id_set[*id] = true;
        }

        let mut vertex_clusters: HashMap<TreeNodeID<Self>, BitVec> = vec![].into_iter().collect();

        let mut remove_list = vec![false; tree.get_capacity()];
        for orig_node_id in tree_postord_node_ids{
            if Some(*orig_node_id)==new_tree_root_parent_id{
                break;
            }
            let mut node = tree.get_node(*orig_node_id).cloned().unwrap();
            let mut bitset = BitVec::from_zeros(taxa_map.len());
            match node.is_leaf() {
                true => {
                    if leaf_id_set[node.get_id()] {
                        node_map[node.get_id()] = Some(node.clone());
                        let node_taxa = node.get_taxa().unwrap();
                        let taxa_id = taxa_map.get(node_taxa).unwrap();
                        bitset.flip_bit_unchecked(*taxa_id);
                        vertex_clusters.insert(node.get_id(), bitset);    
                    }
                }
                false => {
                    let node_children_ids = node.get_children().collect_vec();
                    for child_id in node_children_ids.iter() {
                        match node_map[*child_id].is_some() {
                            true => {}
                            false => node.remove_child(child_id),
                        }
                    }
                    let node_children_ids = node.get_children().collect_vec();
                    match node_children_ids.len() {
                        0 => {}
                        1 => {
                            // the node is a unifurcation
                            // node should be added to both node_map and remove_list
                            // if child of node is already in remove list, attach node children to node first
                            let child_node_id = node_children_ids[0];

                            if remove_list[child_node_id] {
                                node.remove_child(&child_node_id);
                                let grandchildren_ids = node_map[child_node_id]
                                    .as_mut()
                                    .unwrap()
                                    .get_children()
                                    .collect_vec();
                                for grandchild_id in grandchildren_ids {
                                    node_map[grandchild_id]
                                        .as_mut()
                                        .unwrap()
                                        .set_parent(Some(node.get_id()));
                                    node.add_child(grandchild_id);
                                }
                            }
                            let n_id = node.get_id();
                            remove_list[n_id] = true;
                            node_map[n_id] = Some(node.clone());
                            let c_bitset = vertex_clusters.get(&child_node_id).unwrap();
                            let _ = bitset.apply_mask_or(c_bitset);
                            vertex_clusters.insert(n_id, bitset);
                        }
                        _ => {
                            // node has multiple children
                            // for each child, suppress child if child is in remove list
                            node_children_ids.into_iter().for_each(|chid| {
                                if remove_list[chid] {
                                    // suppress chid
                                    // remove chid from node children
                                    // children of chid are node grandchildren
                                    // add grandchildren to node children
                                    // set grandchildren parent to node
                                    node.remove_child(&chid);
                                    let node_grandchildren = node_map[chid]
                                        .as_mut()
                                        .unwrap()
                                        .get_children()
                                        .collect_vec();
                                    for grandchild in node_grandchildren {
                                        node.add_child(grandchild);
                                        node_map[grandchild]
                                            .as_mut()
                                            .unwrap()
                                            .set_parent(Some(node.get_id()))
                                    }
                                }
                            });
                            if node.get_id() == new_tree_root_id {
                                node.set_parent(None);
                            }
                            node_map[node.get_id()] = Some(node.clone());
                            tree.get_node_children_ids(*orig_node_id)
                                .map(|c_id| vertex_clusters.get(&c_id).unwrap())
                                .for_each(|c_bitset| {let _ = bitset.apply_mask_or(c_bitset);});
                            vertex_clusters.insert(*orig_node_id, bitset);
                        }
                    };
                }
            }
        };
        remove_list.into_iter().enumerate().for_each(|(n_id, x)| {
            if x {
                node_map[n_id] = None;
                vertex_clusters.remove(&n_id);
            }
        });

        let new_nodes = node_map.into_iter().flatten().collect_vec();

        let mut new_tree = tree.clone();
        new_tree.set_root(new_tree_root_id);
        new_tree.clear();
        new_tree.set_nodes(new_nodes.into_iter());

        let mut median_node_id: TreeNodeID<Self> = new_tree.get_root_id();
        loop {
            median_node_id = new_tree.get_node_children_ids(median_node_id)
                .max_by(|x, y| {
                    let x_cluster_size = vertex_clusters.get(x).unwrap().count_ones();
                    let y_cluster_size = vertex_clusters.get(y).unwrap().count_ones();
                    x_cluster_size.cmp(&y_cluster_size)
                })
                .unwrap();
            if vertex_clusters.get(&median_node_id).unwrap().count_ones() as usize <= (leaf_ids.len() / 2) {
                break;
            }
        }

        let nodes_postord = new_tree.postord_ids(new_tree.get_root_id()).collect_vec();

        let lower_set: HashSet<_> = new_tree.postord_ids(median_node_id).collect();
        let lower_node_ids_postord = new_tree.postord_ids(median_node_id).collect_vec();
        let upper_node_ids_postord = nodes_postord.iter().filter(|x| !lower_set.contains(x))
        .map(|x| x.clone())
        .collect_vec();


        (new_tree, median_node_id, vertex_clusters, nodes_postord, lower_node_ids_postord, upper_node_ids_postord)

        // (new_tree, median_node_id, vertex_clusters, vec![], vec![])
    }


}
