use clap::{arg, Command};
use std::fs;
use indicatif::{ProgressBar, ProgressStyle};
use itertools::Itertools;
use nlcd::nlcd::near_linear_cophenetic_distance::{LcaMap, NearLinearCopheneticDistance};
use phylo::prelude::*;
use phylo::tree::{DemoTree, PhyloTree};
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
                    arg!(-t --num_threads <NUM_THREADS> "number of threads. Set to 0 to use all threads. Default=1")
                        .required(false)
                        .default_value("1")
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
                    arg!(-n --num_taxa <NUM_START_TAXA> "number of taxa")
                        .required(true)
                        .value_parser(clap::value_parser!(usize)),
                )
                .arg(
                    arg!(-i --num_iter <NUM_ITER> "Step size")
                        .required(true)
                        .value_parser(clap::value_parser!(u32)),
                )
                .arg(
                    arg!(-t --num_threads <NUM_THREADS> "number of threads. Set to 0 to use all threads. Default=1")
                        .required(false)
                        .default_value("1")
                        .value_parser(clap::value_parser!(usize)),
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
                .arg(
                    arg!(-m --method <METHOD> "One of size, depth, or height")
                        .required(true)
                        .default_value("depth")
                        .value_parser(clap::value_parser!(String)),
                )
                .arg(
                    arg!(-w --weighted <WEIGHTED> "Use edge weights")
                        .required(true)
                        .default_value("false")
                        .value_parser(clap::value_parser!(bool)),
                )
        )
        .about("CLI tool for quick tree operations")
        .get_matches();
    

    match matches.subcommand() {
        Some(("repr-emp", sub_m)) => {
            // Returns node depth
            fn depth(
                tree: &DemoTree,
                node_id:  <<DemoTree as RootedTree>::Node as RootedTreeNode>::NodeID,
            ) -> f32 {
                EulerWalk::get_node_depth(tree, node_id) as f32
            }

            // Returns node cluster size
            fn size(
                tree: &DemoTree,
                node_id:  <<DemoTree as RootedTree>::Node as RootedTreeNode>::NodeID,
            ) -> f32 {
                tree.get_cluster_size(node_id) as f32
            }

            // Sets zeta to be node heights
            fn set_node_heights(tree: &mut DemoTree) {
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
            fn generate_tree(
                num_taxa: &usize,
                dist_type: &str,
                zeta_type: &str,
            ) -> DemoTree {
                let mut tree = match dist_type {
                    "yule" => {
                        DemoTree::yule(*num_taxa)
                        
                    }
                    _ => {
                        DemoTree::unif(*num_taxa)
                    }
                };
                tree.precompute_constant_time_lca();
                match zeta_type {
                    "depth" => {
                        tree.set_zeta(depth);
                    }
                    "size" => {
                        tree.set_zeta(size);
                    }
                    _ => {
                        set_node_heights(&mut tree);
                    }
                }    
                tree
            }

            // Empirical study results
            let dist_types = ["uniform", "yule"];
            let zeta_types = ["depth", "size", "height"];
            let norm = sub_m.get_one::<u32>("norm").expect("required");
            let num_trees = sub_m.get_one::<usize>("num_trees").expect("required");
            let num_taxa = sub_m.get_one::<usize>("num_taxa").expect("required");
            let num_threads = sub_m.get_one::<usize>("num_threads").unwrap();
            rayon::ThreadPoolBuilder::new().num_threads(*num_threads).build_global().unwrap();

            let mut output_file =
                File::create(sub_m.get_one::<String>("out_file").expect("required")).unwrap();
            println!("Number of trees: {}", num_trees);
            println!("Number of taxa per tree: {}", num_taxa);
            println!("Norm: {}", norm);

            
            let mut all_dists = vec![];

            let num_steps = dist_types.len() * zeta_types.len() * (*norm as usize + 4);

            let pb = ProgressBar::new(num_steps as u64);
            pb.set_style(
                ProgressStyle::with_template(
                    "{spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] ({eta})",
                )
                .unwrap()
                .progress_chars("#>-"),
            );

            for dist_type in dist_types {
                for zeta_type in zeta_types {
                    for p in 1..norm + 1 {
                        let dist = (0..*num_trees).par_bridge()
                            .map(|_| {
                                let t1 = generate_tree(num_taxa, dist_type, zeta_type);
                                let t2 = generate_tree(num_taxa, dist_type, zeta_type);
                                t1.cophen_dist(&t2, p)
                            })
                            .collect::<Vec<_>>();
                        let out_str = dist.into_iter().join(",");
                        let run_string = format!("{},{},{}:{}", dist_type, zeta_type, p, out_str);
                        all_dists.push(run_string);
                        pb.inc(1);
                    }
                    let dist = (0..*num_trees).par_bridge()
                            .map(|_| {
                                let t1 = generate_tree(num_taxa, dist_type, zeta_type);
                                let t2 = generate_tree(num_taxa, dist_type, zeta_type);
                                t1.cophen_dist(&t2, 0)
                            })
                            .collect::<Vec<_>>();
                        let out_str = dist.into_iter().join(",");
                        let run_string = format!("{},{},inf:{}", dist_type, zeta_type, out_str);
                        all_dists.push(run_string);
                        pb.inc(1);
                    
                    let dist = (0..*num_trees).par_bridge()
                        .map(|_| {
                            let t1 = generate_tree(num_taxa, dist_type, zeta_type);
                            let t2 = generate_tree(num_taxa, dist_type, zeta_type);
                            t1.cophen_dist(&t2, 20)
                        })
                        .collect::<Vec<_>>();
                    let out_str = dist.into_iter().join(",");
                    let run_string = format!("{},{},20:{}", dist_type, zeta_type, out_str);
                    all_dists.push(run_string);
                    pb.inc(1);

                    let dist = (0..*num_trees).par_bridge()
                        .map(|_| {
                            let t1 = generate_tree(num_taxa, dist_type, zeta_type);
                            let t2 = generate_tree(num_taxa, dist_type, zeta_type);
                            t1.cophen_dist(&t2, 50)
                        })
                        .collect::<Vec<_>>();
                    let out_str = dist.into_iter().join(",");
                    let run_string = format!("{},{},50:{}", dist_type, zeta_type, out_str);
                    all_dists.push(run_string);
                    pb.inc(1);

                    let dist = (0..*num_trees).par_bridge()
                        .map(|_| {
                            let t1 = generate_tree(num_taxa, dist_type, zeta_type);
                            let t2 = generate_tree(num_taxa, dist_type, zeta_type);
                            t1.cophen_dist(&t2, 100)
                        })
                        .collect::<Vec<_>>();
                    let out_str = dist.into_iter().join(",");
                    let run_string = format!("{},{},100:{}", dist_type, zeta_type, out_str);
                    all_dists.push(run_string);
                    pb.inc(1);
                }
            }

            output_file
                .write_all(all_dists.join("\n").as_bytes())
                .unwrap();
        }

        Some(("repr-sca", sub_m)) => {
            fn depth(
                tree: &DemoTree,
                node_id:  <<DemoTree as RootedTree>::Node as RootedTreeNode>::NodeID,
            ) -> f32 {
                EulerWalk::get_node_depth(tree, node_id) as f32
            }

            fn runtimes_naive(
                trees: &[(DemoTree, DemoTree)],
                norm: u32,
            ) -> Vec<f32> {
                trees
                    .par_iter()
                    .map(|(t1, t2)| {
                        let t1_lca = LcaMap::from_tree(t1);
                        let t2_lca = LcaMap::from_tree(t2);   
                        let now = Instant::now();
                        let _ = t1.naive_cophen_dist(t2, &t1_lca, &t2_lca, norm);
                        now.elapsed().as_secs_f32()
                    })
                    .collect::<Vec<f32>>()
            }

            fn runtimes_nlcd(
                trees: &[(DemoTree, DemoTree)],
                norm: u32,
            ) -> Vec<f32> {
                trees
                .par_iter()
                .map(|(t1, t2)| {
                        let t1_lca = LcaMap::from_tree(t1);
                        let t2_lca = LcaMap::from_tree(t2);   
                        let now = Instant::now();
                        let _ = t1.nl_cophen_dist(t2, &t1_lca, &t2_lca, norm);
                        now.elapsed().as_secs_f32()
                    })
                    .collect::<Vec<f32>>()
            }


            // let norm = sub_m.get_one::<u32>("norm").expect("required");
            // let num_taxa = sub_m.get_one::<usize>("num_taxa").expect("required");
            let num_iter = sub_m.get_one::<u32>("num_iter").expect("required");
            let num_threads = sub_m.get_one::<usize>("num_threads").unwrap();
            rayon::ThreadPoolBuilder::new().num_threads(*num_threads).build_global().unwrap();


            let mut output_file =
                File::create(sub_m.get_one::<String>("out_file").expect("required")).unwrap();

            let mut fname = sub_m.get_one::<String>("out_file").cloned().expect("required");
            fname.push('1');

            let mut output_file2 =
                File::create(fname).unwrap();

            // let mut lines: Vec<String> = vec![];


            // for taxa_size in (25..601).step_by(25){
            
            //     println!("Generating trees for n={}", taxa_size);
            //     // Generating trees
            //     let trees: Vec<(DemoTree, DemoTree)> = (0..*num_iter)
            //         .map(|_| {
            //             // dbg!(&x);
            //             let mut t1 = DemoTree::yule(taxa_size);
            //             let mut t2 = DemoTree::yule(taxa_size);
            //             t1.precompute_constant_time_lca();
            //             t2.precompute_constant_time_lca();
            //             t1.set_zeta(depth);
            //             t2.set_zeta(depth);
            //             (t1, t2)
            //         })
            //         .collect_vec();

            //     for p in [1,2,5,10,20,50,100]{
            //         println!("computing distances for p={}", p);

            //         let naive = format!("naive-{}-{}:{}", taxa_size, p, runtimes_naive(&trees, p).iter().map(|x| x.to_string()).join(","));
            //         let nlcd = format!("nlcd-{}-{}:{}\n", taxa_size, p, runtimes_nlcd(&trees, p).iter().map(|x| x.to_string()).join(","));

            //         // lines.push(naive);
            //         // lines.push(nlcd);

            //         output_file
            //             .write_all(
            //             [naive,nlcd]
            //                     .join("\n")
            //                     .as_bytes(),
            //             )
            //             .unwrap();
            //     }
            // }   
            for taxa_size in (1000..20001).step_by(1000){
            
                println!("Generating trees for n={}", taxa_size);
                // Generating trees
                let trees: Vec<(DemoTree, DemoTree)> = (0..*num_iter)
                    .map(|_| {
                        // dbg!(&x);
                        let mut t1 = DemoTree::yule(taxa_size);
                        let mut t2 = DemoTree::yule(taxa_size);
                        t1.precompute_constant_time_lca();
                        t2.precompute_constant_time_lca();
                        t1.set_zeta(depth);
                        t2.set_zeta(depth);
                        (t1, t2)
                    })
                    .collect_vec();

                // let out = vec![1,2,5,10,20,50,100].par_iter().map(|p| {
                //         let naive = format!("naive-{}-{}:{}", taxa_size, p, runtimes_naive(&trees, *p).iter().map(|x| x.to_string()).join(","));
                //         let nlcd = format!("nlcd-{}-{}:{}\n", taxa_size, p, runtimes_nlcd(&trees, *p).iter().map(|x| x.to_string()).join(","));
                //         vec![naive, nlcd].into_iter().join("\n")
                //     })
                //     .collect::<Vec<_>>();

                // output_file2
                //     .write_all(
                //     out.into_iter()
                //             .join("\n")
                //             .as_bytes(),
                //     )
                //     .unwrap();

                for p in [1,2,5,10,20,50,100]{
                    println!("computing distances for p={}", p);

                    let naive = format!("naive-{}-{}:{}", taxa_size, p, runtimes_naive(&trees, p).iter().map(|x| x.to_string()).join(","));
                    let nlcd = format!("nlcd-{}-{}:{}\n", taxa_size, p, runtimes_nlcd(&trees, p).iter().map(|x| x.to_string()).join(","));

                    // lines.push(naive);
                    // lines.push(nlcd);

                    output_file2
                        .write_all(
                        [naive,nlcd]
                                .join("\n")
                                .as_bytes(),
                        )
                        .unwrap();
                }
            }            

        }

        Some(("test", sub_m)) => {
            fn depth(
                tree: &DemoTree,
                node_id:  <<DemoTree as RootedTree>::Node as RootedTreeNode>::NodeID,
            ) -> f32 {
                EulerWalk::get_node_depth(tree, node_id) as f32
            }

            let norm = sub_m.get_one::<u32>("norm").expect("required");
            let x = sub_m.get_one::<usize>("num_taxa").expect("required");

            let mut t1 = DemoTree::yule(*x);
            let mut t2 = DemoTree::yule(*x);
            t1.precompute_constant_time_lca();
            t2.precompute_constant_time_lca();
            t1.set_zeta(depth).unwrap();
            t2.set_zeta(depth).unwrap();

            let mut total_runtime = Duration::from_secs(0);
            for _ in 0..20 {
                let now = Instant::now();
                let _ = t1.cophen_dist(&t2, *norm);
                total_runtime += now.elapsed();
            }

            println!("Naive: {:?}", (total_runtime / 20).as_secs());

            let mut total_runtime = Duration::from_secs(0);
            for _ in 0..20 {
                let now = Instant::now();
                // let _ = t1.nl_cophen_dist(&t2, *norm);
                total_runtime += now.elapsed();
            }

            println!("Nlcd: {:?}", (total_runtime / 20).as_secs());
        }
        Some(("dist", sub_m)) => {
            // Returns node depth
            fn depth(
                tree: &PhyloTree,
                node_id:  <<PhyloTree as RootedTree>::Node as RootedTreeNode>::NodeID,
            ) -> f32 {
                EulerWalk::get_node_depth(tree, node_id) as f32
            }

            // Returns node cluster size
            fn size(
                tree: &PhyloTree,
                node_id:  <<PhyloTree as RootedTree>::Node as RootedTreeNode>::NodeID,
            ) -> f32 {
                tree.get_cluster_size(node_id) as f32
            }

            // Sets zeta to be node heights
            fn set_node_heights(tree: &mut PhyloTree) {
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


            let norm = sub_m.get_one::<u32>("norm").expect("required");
            let input_file = matches.get_one::<String>("input_file").expect("tree-file argument required");
            let method = matches.get_one::<String>("method").expect("method required").as_str();
            let weighted = matches.get_one::<bool>("weighted").expect("weighted required");

            let contents = fs::read_to_string(input_file)
                .expect("Should have been able to read the file");

            let tree_strings = contents.split("\n").collect_vec();

            let mut t1 = PhyloTree::from_newick(tree_strings[0].as_bytes()).unwrap();
            let mut t2 = PhyloTree::from_newick(tree_strings[1].as_bytes()).unwrap();

            t1.precompute_constant_time_lca();
            t2.precompute_constant_time_lca();

            match method{
                "size" => {
                    t1.set_zeta(size);
                    t2.set_zeta(size);        
                }
                "height" => {
                    set_node_heights(&mut t1);
                    set_node_heights(&mut t2);        
                }
                _ => {
                    t1.set_zeta(depth);
                    t2.set_zeta(depth);        
                }
            };

            let t1_lca = LcaMap::from_tree(&t1);
            let t2_lca = LcaMap::from_tree(&t2);   

            println!("Norm: {}, Cophenetic-distance: {}", norm, t1.nl_cophen_dist(&t2, &t1_lca, &t2_lca, *norm));
        }
        _ => {
            println!("No option selected! Refer help page (-h flag)");
        }
    }
}
