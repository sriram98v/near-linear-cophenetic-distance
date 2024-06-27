pub mod near_linear_cophenetic_distance;

use fxhash::FxHashSet as HashSet;

use near_linear_cophenetic_distance::NearLinearCopheneticDistance;
use phylo::tree::SimpleRootedTree;

// pub enum TaxaIntersectionType{
//     UpperUpper,
//     UpperLower,
//     LowerUpper,
//     LowerLower,
// }

// pub struct TaxaSets{
//     a: HashSet<<SimpleRootedTree as NearLinearCopheneticDistance>::Meta>, 
//     b: HashSet<<SimpleRootedTree as NearLinearCopheneticDistance>::Meta>,
//     a_hat: HashSet<<SimpleRootedTree as NearLinearCopheneticDistance>::Meta>, 
//     b_hat: HashSet<<SimpleRootedTree as NearLinearCopheneticDistance>::Meta>,
// }

// pub struct TaxaIntersections{
//     a_int_a_hat: Vec<<SimpleRootedTree as NearLinearCopheneticDistance>::Meta>,
//     a_int_b_hat: Vec<<SimpleRootedTree as NearLinearCopheneticDistance>::Meta>,
//     b_int_a_hat: Vec<<SimpleRootedTree as NearLinearCopheneticDistance>::Meta>,
//     b_int_b_hat: Vec<<SimpleRootedTree as NearLinearCopheneticDistance>::Meta>,
// }


impl NearLinearCopheneticDistance for SimpleRootedTree{
    type Meta = String;
}