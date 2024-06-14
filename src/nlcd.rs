pub mod near_linear_cophenetic_distance;

use near_linear_cophenetic_distance::NearLinearCopheneticDistance;
use phylo::tree::SimpleRootedTree;

impl NearLinearCopheneticDistance for SimpleRootedTree{
    type Meta = String;
}