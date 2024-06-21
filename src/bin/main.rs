use nlcd::nlcd::near_linear_cophenetic_distance::NearLinearCopheneticDistance;
use std::fs::File;
use std::io::Write;
use phylo::prelude::*;
use clap::{arg, Command};
use phylo::tree::SimpleRootedTree;
use itertools::Itertools;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;

fn main(){
    let matches = Command::new("Generalized suffix tree")
        .version("1.0")
        .author("Sriram Vijendran <vijendran.sriram@gmail.com>")
        .subcommand(Command::new("repr")
            .about("Reproduce results from article")
            .arg(arg!(-k --norm <NORM> "nth norm")
                .value_parser(clap::value_parser!(u32))
            )
            .arg(arg!(-n --num_trees <NUM_TREES> "number of trees")
                .required(true)
                .value_parser(clap::value_parser!(usize))
            )
            .arg(arg!(-x --num_taxa <NUM_TAXA> "number of taxa")
                .required(true)
                .value_parser(clap::value_parser!(usize))
            )
            .arg(arg!(-t --threads <THREADS> "number of threads")
                .value_parser(clap::value_parser!(usize))
            )
            .arg(arg!(-o --out_file <OUT_FILE> "output file")
                .required(true)
                .value_parser(clap::value_parser!(String))
            )
        )
        .about("CLI tool for quick tree operations")
        .get_matches();

        match matches.subcommand(){
            Some(("repr",  sub_m)) => {            
                // Returns node depth
                fn depth(tree: &SimpleRootedTree, node_id: <SimpleRootedTree as RootedTree>::NodeID)->f32
                {
                    EulerWalk::get_node_depth(tree, node_id) as f32
                }

                // Returns node cluster size
                fn size(tree: &SimpleRootedTree, node_id: <SimpleRootedTree as RootedTree>::NodeID)->f32
                {
                    tree.get_cluster_size(node_id) as f32
                }

                // Sets zeta to be node heights
                fn set_node_heights(tree: &mut SimpleRootedTree)
                {
                    let node_post_ord = tree.postord(tree.get_root_id()).map(|x| x.get_id()).collect_vec();
                    for node_id in node_post_ord{
                        match tree.is_leaf(&node_id){
                            true => {tree.get_node_mut(node_id).unwrap().set_zeta(Some(0_f32))},
                            false => {
                                let max_height = tree.get_node_children_ids(node_id)
                                    .map(|x| tree.get_zeta(x).unwrap() as u32)
                                    .max()
                                    .unwrap_or(0);
                                tree.get_node_mut(node_id).unwrap().set_zeta(Some((max_height+1) as f32));
                            },
                        }
                    }
                }

                // Generates trees as per parameters
                fn generate_trees(num_trees: &usize, num_taxa: &usize, dist_type: &str)->Vec<SimpleRootedTree>
                {
                    (0..*num_trees).map(|_| {
                        match dist_type{
                            "yule" => {
                                let mut tree = SimpleRootedTree::yule(*num_taxa).unwrap();
                                tree.precompute_constant_time_lca();
                                return tree;
                            },
                            _ => {
                                let mut tree = SimpleRootedTree::unif(*num_taxa).unwrap();
                                tree.precompute_constant_time_lca();
                                return tree;
                            },
                        };
                    }).collect_vec()
                }

                fn set_treeset_zeta(tree_set: &mut Vec<SimpleRootedTree>, zeta_type: &str)
                {
                    for tree in tree_set.iter_mut(){
                        match zeta_type {
                            "depth" => {
                                tree.set_zeta(depth);
                            },
                            "size" => {
                                tree.set_zeta(size);
                            },
                            _ => {
                                set_node_heights(tree);
                            },

                        }
                    }
                }

                // Returns distribution of pairwise cophenetic distances for a set of trees
                fn distance_distribution(tree_set: &Vec<SimpleRootedTree>, norm: u32)->Vec<f32>
                {
                    tree_set.iter().combinations(2)
                    .par_bridge()
                    .map(|x| {
                        let t1 = x[0];
                        let t2 = x[1];
                        let dist = t1.cophen_dist_naive(t2, norm);
                        return dist;
                    }).collect::<Vec<f32>>()
                }
                
                let dist_types = ["uniform", "yule"];
                let zeta_types = ["depth", "size", "height"];
                let norm = sub_m.get_one::<u32>("norm").expect("required");
                let num_trees = sub_m.get_one::<usize>("num_trees").expect("required");
                let num_taxa = sub_m.get_one::<usize>("num_taxa").expect("required");

                // todo: add functionality to set number of threads
                let num_threads: usize = *sub_m.get_one::<usize>("threads").unwrap_or(&1);
                let mut output_file = File::create(sub_m.get_one::<String>("out_file").expect("required")).unwrap();
                println!("Number of trees: {}", num_trees);
                println!("Number of taxa per tree: {}", num_taxa);
                println!("Norm: {}", norm);

                let mut all_dists = vec![];

                let num_steps = dist_types.len()*zeta_types.len()*(norm.clone() as usize);

                let pb = ProgressBar::new(num_steps as u64);
                pb.set_style(ProgressStyle::with_template("{spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] ({eta})")
                    .unwrap()
                    .progress_chars("#>-"));


                for dist_type in dist_types{
                    let mut tree_set = generate_trees(num_trees, num_taxa, dist_type);
                    for zeta_type in zeta_types{
                        set_treeset_zeta(&mut tree_set, zeta_type);
                        for p in 1..norm+1{
                            let dist = distance_distribution(&tree_set, p).iter().map(|x| x.to_string()).join(",");
                            let run_string = format!("{},{},{}:{}", dist_type, zeta_type, p, dist);
                            all_dists.push(run_string);
                            pb.inc(1);
                        }
                    }
                }
                
                output_file.write_all(all_dists.join("\n").as_bytes()).unwrap();            
            },
            _ => {
                println!("No option selected! Refer help page (-h flag)");
            }
        }
}