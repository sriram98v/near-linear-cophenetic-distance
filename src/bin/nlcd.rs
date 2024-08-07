use clap::{arg, Command};
use std::fs;
use indicatif::{ProgressBar, ProgressStyle};
use itertools::Itertools;
use nlcd::nlcd::near_linear_cophenetic_distance::NearLinearCopheneticDistance;
use phylo::prelude::*;
use phylo::tree::SimpleRootedTree;
use rayon::prelude::*;
use std::fs::File;
use std::io::Write;
use std::time::{Duration, Instant};

fn main() {
    let matches = Command::new("Generalized suffix tree")
        .version("1.0")
        .author("Sriram Vijendran <vijendran.sriram@gmail.com>")
        .subcommand(
            Command::new("repr-emp")
                .about("Reproduce empirical results from article")
                .arg(arg!(-p --norm <NORM> "nth norm").required(true).value_parser(clap::value_parser!(u32)))
                .arg(
                    arg!(-n --num_trees <NUM_TREES> "number of trees")
                        .required(true)
                        .value_parser(clap::value_parser!(usize)),
                )
                .arg(
                    arg!(-x --num_taxa <NUM_TAXA> "number of taxa")
                        .required(true)
                        .value_parser(clap::value_parser!(usize)),
                )
                .arg(
                    arg!(-o --out_file <OUT_FILE> "output file")
                        .required(true)
                        .value_parser(clap::value_parser!(String)),
                ),
        )
        .subcommand(
            Command::new("repr-sca")
                .about("Reproduce scalability results from article")
                .arg(arg!(-p --norm <NORM> "nth norm").required(true).value_parser(clap::value_parser!(u32)))
                .arg(
                    arg!(-s --num_start_taxa <NUM_START_TAXA> "number of starting taxa")
                        .required(true)
                        .value_parser(clap::value_parser!(usize)),
                )
                .arg(
                    arg!(-e --num_end_taxa <NUM_END_TAXA> "number of ending taxa")
                        .required(true)
                        .value_parser(clap::value_parser!(usize)),
                )
                .arg(
                    arg!(-x --step_size <STEP_SIZE> "Step size")
                        .required(true)
                        .value_parser(clap::value_parser!(usize)),
                )
                .arg(
                    arg!(-i --num_iter <NUM_ITER> "Step size")
                        .required(true)
                        .value_parser(clap::value_parser!(u32)),
                )
                .arg(
                    arg!(-o --out_file <OUT_FILE> "output file")
                        .required(true)
                        .value_parser(clap::value_parser!(String)),
                ),
        )
        .subcommand(
            Command::new("test")
                .arg(arg!(-p --norm <NORM> "nth norm").required(true).value_parser(clap::value_parser!(u32)))
                .arg(
                    arg!(-x --num_taxa <NUM_TAXA> "number of starting taxa")
                        .required(true)
                        .value_parser(clap::value_parser!(usize)),
                )
        )
        .subcommand(
            Command::new("dist")
                .arg(
                    arg!(-p --norm <NORM> "nth norm")
                        .required(true)
                        .value_parser(clap::value_parser!(u32))
                )
                .arg(
                    arg!(-i --input_file <FILE_PATH> "Input tree file in Newick format")
                        .required(true)
                        .value_parser(clap::value_parser!(String)),
                )
        )
        .about("CLI tool for quick tree operations")
        .get_matches();
    

    match matches.subcommand() {
        Some(("repr-emp", sub_m)) => {
            // Returns node depth
            fn depth(
                tree: &SimpleRootedTree,
                node_id:  <<SimpleRootedTree as RootedTree>::Node as RootedTreeNode>::NodeID,
            ) -> f32 {
                EulerWalk::get_node_depth(tree, node_id) as f32
            }

            // Returns node cluster size
            fn size(
                tree: &SimpleRootedTree,
                node_id:  <<SimpleRootedTree as RootedTree>::Node as RootedTreeNode>::NodeID,
            ) -> f32 {
                tree.get_cluster_size(node_id) as f32
            }

            // Sets zeta to be node heights
            fn set_node_heights(tree: &mut SimpleRootedTree) {
                let node_post_ord = tree
                    .postord_ids(tree.get_root_id())
                    // .map(|x| x.get_id())
                    .collect_vec();
                for node_id in node_post_ord {
                    match tree.is_leaf(node_id) {
                        true => tree.get_node_mut(node_id).unwrap().set_zeta(Some(0_f32)),
                        false => {
                            let max_height = tree
                                .get_node_children_ids(node_id)
                                .map(|x| tree.get_zeta(x).unwrap() as u32)
                                .max()
                                .unwrap_or(0);
                            tree.get_node_mut(node_id)
                                .unwrap()
                                .set_zeta(Some((max_height + 1) as f32));
                        }
                    }
                }
            }

            // Generates trees as per parameters
            fn generate_trees(
                num_trees: &usize,
                num_taxa: &usize,
                dist_type: &str,
            ) -> Vec<SimpleRootedTree> {
                (0..*num_trees)
                    .map(|_| match dist_type {
                        "yule" => {
                            let mut tree = SimpleRootedTree::yule(*num_taxa).unwrap();
                            tree.precompute_constant_time_lca();
                            tree
                        }
                        _ => {
                            let mut tree = SimpleRootedTree::unif(*num_taxa).unwrap();
                            tree.precompute_constant_time_lca();
                            tree
                        }
                    })
                    .collect_vec()
            }

            fn set_treeset_zeta(tree_set: &mut [SimpleRootedTree], zeta_type: &str) {
                for tree in tree_set.iter_mut() {
                    match zeta_type {
                        "depth" => {
                            tree.set_zeta(depth);
                        }
                        "size" => {
                            tree.set_zeta(size);
                        }
                        _ => {
                            set_node_heights(tree);
                        }
                    }
                }
            }

            // Returns distribution of pairwise cophenetic distances for a set of trees
            fn distance_distribution(tree_set: &[SimpleRootedTree], norm: u32) -> Vec<f32> {
                tree_set
                    .iter()
                    .combinations(2)
                    .par_bridge()
                    .map(|x| {
                        let t1 = x[0];
                        let t2 = x[1];

                        t1.cophen_dist_naive(t2, norm)
                    })
                    .collect::<Vec<f32>>()
            }

            // Empirical study results
            let dist_types = ["uniform", "yule"];
            let zeta_types = ["depth", "size", "height"];
            let norm = sub_m.get_one::<u32>("norm").expect("required");
            let num_trees = sub_m.get_one::<usize>("num_trees").expect("required");
            let num_taxa = sub_m.get_one::<usize>("num_taxa").expect("required");

            let mut output_file =
                File::create(sub_m.get_one::<String>("out_file").expect("required")).unwrap();
            println!("Number of trees: {}", num_trees);
            println!("Number of taxa per tree: {}", num_taxa);
            println!("Norm: {}", norm);

            let mut all_dists = vec![];

            let num_steps = dist_types.len() * zeta_types.len() * (*norm as usize);

            let pb = ProgressBar::new(num_steps as u64);
            pb.set_style(
                ProgressStyle::with_template(
                    "{spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] ({eta})",
                )
                .unwrap()
                .progress_chars("#>-"),
            );

            for dist_type in dist_types {
                let mut tree_set = generate_trees(num_trees, num_taxa, dist_type);
                for zeta_type in zeta_types {
                    set_treeset_zeta(&mut tree_set, zeta_type);
                    for p in 1..norm + 1 {
                        let dist = distance_distribution(&tree_set, p)
                            .iter()
                            .map(|x| x.to_string())
                            .join(",");
                        let run_string = format!("{},{},{}:{}", dist_type, zeta_type, p, dist);
                        all_dists.push(run_string);
                        pb.inc(1);
                    }
                }
            }

            output_file
                .write_all(all_dists.join("\n").as_bytes())
                .unwrap();
        }

        Some(("repr-sca", sub_m)) => {
            fn depth(
                tree: &SimpleRootedTree,
                node_id:  <<SimpleRootedTree as RootedTree>::Node as RootedTreeNode>::NodeID,
            ) -> f32 {
                EulerWalk::get_node_depth(tree, node_id) as f32
            }

            fn median(numbers: &mut [Duration]) -> Duration {
                numbers.sort();
                let mid = numbers.len() / 2;
                numbers[mid]
            }            

            fn mean_runtime_naive(
                trees: &[(SimpleRootedTree, SimpleRootedTree)],
                num_iter: &u32,
                norm: u32,
            ) -> String {
                trees
                    .iter()
                    .map(|(t1, t2)| {
                        let mut total_runtime = vec![];
                        for _ in 0..*num_iter {
                            let now = Instant::now();
                            let _ = t1.cophen_dist_naive(t2, norm);
                            total_runtime.push(now.elapsed());
                        }
                        median(total_runtime.as_mut())
                    })
                    .map(|x| format!("{}", x.as_millis()))
                    .collect_vec()
                    .join(",")
            }

            fn mean_runtime_nlcd(
                trees: &[(SimpleRootedTree, SimpleRootedTree)],
                num_iter: &u32,
                norm: u32,
            ) -> String {
                trees
                    .iter()
                    .map(|(t1, t2)| {
                        let mut total_runtime = vec![];
                        for _ in 0..*num_iter {
                            let now = Instant::now();
                            let _ = t1.cophen_dist(t2, norm);
                            total_runtime.push(now.elapsed());
                        }
                        median(total_runtime.as_mut())
                    })
                    .map(|x| format!("{}", x.as_millis()))
                    .collect_vec()
                    .join(",")
            }

            let norm = sub_m.get_one::<u32>("norm").expect("required");
            let num_start_taxa = sub_m.get_one::<usize>("num_start_taxa").expect("required");
            let step_size = sub_m.get_one::<usize>("step_size").expect("required");
            let num_end_taxa =
                sub_m.get_one::<usize>("num_end_taxa").expect("required") + step_size;
            let num_iter = sub_m.get_one::<u32>("num_iter").expect("required");

            let mut output_file =
                File::create(sub_m.get_one::<String>("out_file").expect("required")).unwrap();

            // Generating trees
            let trees: Vec<(SimpleRootedTree, SimpleRootedTree)> = (*num_start_taxa..num_end_taxa)
                .step_by(*step_size)
                .map(|x| {
                    dbg!(&x);
                    let mut t1 = SimpleRootedTree::yule(x).unwrap();
                    let mut t2 = SimpleRootedTree::yule(x).unwrap();
                    t1.precompute_constant_time_lca();
                    t2.precompute_constant_time_lca();
                    t1.set_zeta(depth);
                    t2.set_zeta(depth);
                    (t1, t2)
                })
                .collect_vec();

            for i in 1..norm+1{
                println!("Computing distance for norm={i}");
                let naive = format!("naive_{}:{}", i, mean_runtime_naive(&trees, num_iter, i));
                let nlcd = format!("nlcd_{}:{}\n", i, mean_runtime_nlcd(&trees, num_iter, i));
                output_file
                    .write_all(
                    [naive, nlcd]
                            .join("\n")
                            .as_bytes(),
                    )
                    .unwrap();
            }

        }

        Some(("test", sub_m)) => {
            fn depth(
                tree: &SimpleRootedTree,
                node_id:  <<SimpleRootedTree as RootedTree>::Node as RootedTreeNode>::NodeID,
            ) -> f32 {
                EulerWalk::get_node_depth(tree, node_id) as f32
            }

            let norm = sub_m.get_one::<u32>("norm").expect("required");
            let x = sub_m.get_one::<usize>("num_taxa").expect("required");

            let mut t1 = SimpleRootedTree::yule(*x).unwrap();
            let mut t2 = SimpleRootedTree::yule(*x).unwrap();
            t1.precompute_constant_time_lca();
            t2.precompute_constant_time_lca();
            t1.set_zeta(depth);
            t2.set_zeta(depth);

            let mut total_runtime = Duration::from_secs(0);
            for _ in 0..20 {
                let now = Instant::now();
                let _ = t1.cophen_dist_naive(&t2, *norm);
                total_runtime += now.elapsed();
            }

            println!("Naive: {:?}", (total_runtime / 20).as_secs());

            let mut total_runtime = Duration::from_secs(0);
            for _ in 0..20 {
                let now = Instant::now();
                let _ = t1.cophen_dist(&t2, *norm);
                total_runtime += now.elapsed();
            }

            println!("Nlcd: {:?}", (total_runtime / 20).as_secs());
        }
        Some(("dist", sub_m)) => {
            fn depth(
                tree: &SimpleRootedTree,
                node_id:  <<SimpleRootedTree as RootedTree>::Node as RootedTreeNode>::NodeID,
            ) -> f32 {
                EulerWalk::get_node_depth(tree, node_id) as f32
            }

            let norm = sub_m.get_one::<u32>("norm").expect("required");
            let input_file = matches.get_one::<String>("input_file").expect("tree-file argument required");

                let contents = fs::read_to_string(input_file)
                    .expect("Should have been able to read the file");

                let tree_strings = contents.split("\n").collect_vec();

                let mut t1 = SimpleRootedTree::from_newick(tree_strings[0].as_bytes());
                let mut t2 = SimpleRootedTree::from_newick(tree_strings[1].as_bytes());

                t1.precompute_constant_time_lca();
                t2.precompute_constant_time_lca();

                t1.set_zeta(depth);
                t2.set_zeta(depth);

                println!("Norm: {}, Cophenetic-distance: {}", norm, t1.cophen_dist_naive(&t2, *norm));
        }
        _ => {
            println!("No option selected! Refer help page (-h flag)");
        }
    }
}
