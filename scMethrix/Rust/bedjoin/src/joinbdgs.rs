use crate::collectargs;
use indicatif::ProgressBar;
use rust_htslib::tbx::{self, Read};
use std::collections::BTreeMap;
use std::collections::HashMap;
use std::iter::StepBy;
use std::str::from_utf8;
use regex::Regex;

pub fn inputs() {
    let args = collectargs::args();
    mergebdgs(args);
}

pub fn mergebdgs(args: collectargs::Cmdargs) {
    let mut chrlist: Vec<String> = Vec::new();
    let filelist: Vec<String> = args.files.to_vec();
	
	let re = Regex::new(r"/(\w*).").unwrap();
	let files: Vec<_> = filelist.iter().map(|x| re.captures(x).unwrap().get(1).unwrap().as_str()).collect::<Vec<_>>(); 
	
	if args.header == true {
		println!("chr\tstart\tend\t{}", files.join("\t"));
	}
	
    if args.chroms.is_empty() {
        chrlist = getchrlist(&args.files);
    } else {
        chrlist = args.chroms.to_vec();
    }

	let mut tripletmatrix = TripletMatrix::new(
		filelist,
		chrlist,
	);

    let mut chrtree = Chrtree::new(
        args.chrdb,
        args.methr,
        args.umethr,
        args.beta,
        args.verbose,
        args.bin,
    );

    let pb = ProgressBar::new(chrlist.len() as u64);
    //pb.set_style(sty);
    for chr in chrlist {
        pb.inc(1);
        pb.set_message(&chr);
        if chrtree.bin > 0 {
            chrtree.getrangedb(&filelist, &chr);
        } else {
            chrtree.processchr(&chr, &filelist);
        }
    }
    pb.finish_with_message("Done");
}

pub fn getchrlist(files: &Vec<String>) -> Vec<String> {
    let mut chrlist: Vec<String> = Vec::new();
    for f in files {
        let rdr = tbx::Reader::from_path(f).expect(&format!("Could not open {}", f));
        for s in rdr.seqnames() {
            chrlist.push(s);
        }
    }
    chrlist.sort();
    chrlist.dedup();
    chrlist
}

#[derive(Debug)]
struct Chrtree {
    chrlens: HashMap<String, u64>,
    m: u8,
    u: u8,
    b: u8,
    v: bool,
    bin: u64,
}



impl Chrtree {
    fn new(chrlens: HashMap<String, u64>, m: u8, u: u8, b: u8, v: bool, bin: u64) -> Chrtree {
        Chrtree {
            chrlens,
            m,
            u,
            b,
            v,
            bin,
        }
    }

    fn processchr(&mut self, chrname: &String, files: &Vec<String>) -> HashMap<String, String> {
        let mut chrdb = Bdgdb::new(&files.len());
        let mut fidx = 0;

        let mut chravgmeth: HashMap<String, String> = HashMap::new();

        //let pb = ProgressBar::new(files.len() as u64);

        for f in files {
            //pb.set_message(&f);
            //pb.inc(1);
            // if self.v == true {
            //     eprintln!("Processing {}", &f);
            // }
            let mut rdr = tbx::Reader::from_path(f).expect(&format!("Could not open {}", f));
            rdr.set_threads(8).unwrap();
            let x = rdr.tid(chrname);
            let tid = match x {
                Ok(tid) => tid,
                Err(_) => 1000,
            };
            let mut amstr: String = String::from("NA");
            let mut ncpgs: String = String::from("NA");

            if tid != 1000 {
                //let chrsize = self.getchrsize(chrname);

                let chrsize = match self.chrlens.get(chrname) {
                    Some(x) => x,
                    None => panic!("Could not find the desired chromosome in chrom sizes file"),
                };

                if rdr.fetch(tid, 0, *chrsize).is_ok() == true {
                    let mut avgmeth: Vec<f32> = Vec::new();
                    for rec in rdr.records() {
                        let rec = rec.unwrap();
                        let x: f32 = chrdb
                            .update2(rec, fidx, &self.m, &self.u, &self.b)
                            .parse()
                            .unwrap();
                        avgmeth.push(x);
                    }
                    ncpgs = avgmeth.len().to_string();
                    let am: f32 = avgmeth.iter().sum();
                    //am = am / chrdb.db.len() as f32;
                    amstr = (am / (chrdb.db.len() as f32)).to_string() + "|" + &ncpgs;
                }
            } else {
                amstr = amstr + "|" + &ncpgs;
            }
            chravgmeth.insert(f.to_string(), amstr.clone());
            fidx = fidx + 1;
        }
        chrdb.print(chrname);
        drop(chrdb);
        chravgmeth
    }
    pub fn getrangedb(&mut self, files: &Vec<String>, chrname: &String) {
        let chrsize = match self.chrlens.get(chrname) {
            Some(x) => x,
            None => panic!("Cound not find!"),
        };
        //chrdb.get(&chr).unwrap();
        let mut chrranges = makeranges(chrsize, &self.bin);
        let mut start = chrranges.next();
        loop {
            let end = chrranges.next();
            //
            let end2: u64 = match end {
                Some(x) => x,
                None => break,
            };
            let start2: u64 = start.unwrap() as u64;
            // if self.v == true {
            //     eprintln!("Processing: {}:{}-{}", &chrname, &start2, &end2);
            // }
            let mut chrdb = Bdgdb::new(&files.len());
            let mut fidx = 0;
            for f in files {
                // if self.v == true {
                //     eprintln!("Processing {}", &f);
                // }
                let mut rdr = tbx::Reader::from_path(f).expect(&format!("Could not open {}", f));
                rdr.set_threads(8).unwrap();
                let x = rdr.tid(&chrname);
                let tid = match x {
                    Ok(tid) => tid,
                    Err(_) => 1000,
                };
                if tid != 1000 {
                    if rdr.fetch(tid, start2, end2).is_ok() == true {
                        //let mut avgmeth: Vec<f32> = Vec::new();
                        for rec in rdr.records() {
                            let rec = rec.unwrap();
                            let _x: f32 = chrdb
                                .update2(rec, fidx, &self.m, &self.u, &self.b)
                                .parse()
                                .unwrap();
                            //avgmeth.push(x);
                        }
                    }
                }
                fidx = fidx + 1;
            }
            chrdb.print(&chrname);
            drop(chrdb);
            start = end;
        }
    }
}

struct Bdgdb {
    db: BTreeMap<u64, Vec<String>>,
    nfiles: usize,
}

impl Bdgdb {
    fn new(nfiles: &usize) -> Bdgdb {
        let init = vec![String::new(); *nfiles];
        let mut db: BTreeMap<u64, Vec<String>> = BTreeMap::new();
        let nfiles = *nfiles;
        db.insert(0, init);
        Bdgdb { db, nfiles }
    }

    pub fn update2(&mut self, rec: Vec<u8>, idx: usize, m: &u8, u: &u8, b: &u8) -> String {
        let rec = from_utf8(&rec).unwrap();
        let linespl: Vec<&str> = rec.split('\t').collect();
        let k: u64 = linespl[1].parse().unwrap();
        let mut v = String::new();
        if m > &0 && u > &0 {
            let m = (*m - 1) as usize;
            let u = (*u - 1) as usize;
            let mreads: f32 = linespl[m].parse().unwrap();
            let ureads: f32 = linespl[u].parse().unwrap();
            v = (mreads / (mreads + ureads)).to_string();
        } else {
            let b = (*b - 1) as usize;
            v = linespl[b].parse().unwrap();
        }

        if let Some(x) = self.db.get_mut(&k) {
            if let Some(elem) = x.get_mut(idx) {
                *elem = v.clone();
            }
        } else {
            let mut init = vec![String::from("NA"); self.nfiles];
            if let Some(elem) = init.get_mut(idx) {
                *elem = v.clone();
            }
            self.db.insert(k, init);
        }
        v
    }

    pub fn print(&mut self, chr: &str) {
        for (k, v) in &self.db {
            let cols_str: Vec<_> = v.iter().map(ToString::to_string).collect();
			
			

			
            if k != &0 {
                println!("{}\t{}\t{}\t{}", chr, k, k, cols_str.join("\t"));
            }
        }
    }
}

pub fn makeranges(size: &u64, bin: &u64) -> StepBy<std::ops::Range<u64>> {
    let chrange = (0..*size).step_by(*bin as usize);
    chrange
}

#[derive(Debug)]
struct TripletMatrix {
	dimchr: Vec<String>,         
	dimheaders: Vec<String>,   	//Start\tEnd\t[all cell filenames]
	i: Vec<u64>,					//Row (Chr)
	j: Vec<u8>,					//Col (Header)
	x: Vec<u8>,				//Value
}

impl TripletMatrix {
	fn new(dimchr: Vec<String>, dimheaders: Vec<String>) -> TripletMatrix {
		let dimchr = dimchr;
		let dimheaders = dimheaders;
		let mut i: Vec<u64> = Vec::new();
		let mut j: Vec<u8> = Vec::new();
		let mut x: Vec<u8> = Vec::new();
		TripletMatrix {dimchr, dimheaders, i, j, x}
		
	}
		
	fn add_triplet(&mut self, row: u64, col: u8, val: u8) {
		self.i.push(row);
		self.j.push(col);
		self.x.push(val);	
	}
}