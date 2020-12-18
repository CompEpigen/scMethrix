use clap::{App, Arg};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

pub fn args() -> Cmdargs {
  let matches = App::new("cpjoin")
    .version("1.0")
    .about("Outer join tabix bedgraph files from WGBS")
    .author("Anand Mayakonda")
    .arg(
      Arg::with_name("bdgs")
        .short("f")
        .value_name("FILE")
        .help("tabix bedgraph files (Required)")
        .takes_value(true)
        .required(true)
        .multiple(true)
        .display_order(1),
    )
	.arg(
	  Arg::with_name("header")
        .short("e")
        .help("add headers to output")
        .takes_value(false)
        .required(false)
        .multiple(false),
    )
	.arg(
      Arg::with_name("verbose")
        .short("v")
        .help("switch off verbose")
        .takes_value(false)
        .required(false)
        .multiple(false),
    )
    .arg(
      Arg::with_name("chrsize")
        .short("c")
        .help("tab seperatd genome file of format <chromName><TAB><chromSize> (Required)")
        .takes_value(true)
        .required(true)
        .multiple(false)
        .display_order(1),
    )
    .arg(
      Arg::with_name("chrs")
        .short("r")
        .help("Restrict to only these chromosome")
        .takes_value(true)
        .required(false)
        .multiple(true),
    )
    .arg(
      Arg::with_name("binsize")
        .short("s")
        .help("bin sizes to process. Default 0 (i.e, whole chromosme is processed)")
        .takes_value(true)
        .required(false)
        .multiple(false),
    )
    .arg(
      Arg::with_name("beta")
        .short("b")
        .help("column index containing beta values. Default 4")
        .takes_value(true)
        .required(false)
        .multiple(false)
        .display_order(3),
    )
    .arg(
      Arg::with_name("umeth")
        .short("u")
        .help("column index containing un-methylated reads. 0 (Not given)")
        .takes_value(true)
        .required(false)
        .multiple(false)
        .display_order(5),
    )
    .arg(
      Arg::with_name("meth")
        .short("m")
        .help("column index containing methylated reads. Default 0 (Not given)")
        .takes_value(true)
        .required(false)
        .multiple(false)
        .display_order(4),
    )
    .get_matches();

  let bdgs: Vec<_> = matches.values_of("bdgs").unwrap().collect();
  let chromsize: &str = matches.value_of("chrsize").unwrap();

  let mut chroms: Vec<String> = Vec::new();
  if matches.is_present("chrs") {
    if let Some(chrs) = matches.values_of("chrs") {
      for o in chrs {
        chroms.push(o.to_string());
      }
    }
    chroms.sort();
    chroms.dedup();
  }

  let mut umeth: &str = &"0".to_string(); //default value for unmethyltaed col index
  if matches.is_present("umeth") == true {
    umeth = matches.value_of("umeth").unwrap();
  }

  let mut verbose: bool = true;
  if matches.is_present("verbose") == true {
    verbose = false
  }
  
  let mut header: bool = false;
  if matches.is_present("header") == true {
    header = true
  }

  let mut meth: &str = &"0".to_string(); //default value for methyltaed col index
  if matches.is_present("meth") == true {
    meth = matches.value_of("meth").unwrap();
  }

  let mut beta: &str = &"4".to_string(); //default value for beta-value col index
  if matches.is_present("beta") == true {
    beta = matches.value_of("beta").unwrap();
  }

  let mut bin: &str = &"0".to_string(); //default value for bin size (0 = whole chromosome)
  if matches.is_present("binsize") == true {
    bin = matches.value_of("binsize").unwrap();
  }

  Cmdargs::new(bdgs, chromsize, bin, umeth, meth, beta, verbose, header, chroms)
}

#[derive(Debug)]
pub struct Cmdargs {
  pub files: Vec<String>,
  pub chrdb: HashMap<String, u64>,
  pub bin: u64,
  pub umethr: u8,
  pub methr: u8,
  pub beta: u8,
  pub verbose: bool,
  pub header: bool,
  pub chroms: Vec<String>,
}

impl Cmdargs {
  pub fn new(
    f: Vec<&str>,
    c: &str,
    bin: &str,
    um: &str,
    m: &str,
    b: &str,
    verbose: bool,
	header: bool,
    chroms: Vec<String>,
  ) -> Cmdargs {
    let mut files: Vec<String> = Vec::new();
    for b in f {
      files.push(b.to_string());
    }

    //static files: Vec<String> = f.to_vec();

    let bin: u64 = bin.parse().unwrap();

    let mut chrdb: HashMap<String, u64> = HashMap::new();
    let file = File::open(c).expect("Could not open");
    let fr = BufReader::new(file);
    let umethr: u8 = um
      .parse()
      .expect(&format!("Can not convert to {} to u8", um));
    let methr: u8 = m.parse().expect(&format!("Can not convert to {} to u8", m));
    let beta: u8 = b.parse().expect(&format!("Can not convert to {} to u8", b));

    for line in fr.lines() {
      let line = line.unwrap();
      let linespl: Vec<&str> = line.split('\t').collect();
      let chrsize: u64 = linespl[1].parse().unwrap();
      chrdb.insert(linespl[0].to_string(), chrsize);
    }

    Cmdargs {
      files,
      chrdb,
      bin,
      umethr,
      methr,
      beta,
      verbose,
	  header,
      chroms
    }
  }
}
