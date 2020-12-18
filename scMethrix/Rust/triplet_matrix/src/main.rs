extern crate hdf5;

use hdf5::types::VarLenArray;

fn main() -> hdf5::Result<()> {
    
    let mut files: Vec<String> = Vec::new();
    let mut chrs: Vec<String> = Vec::new();

    files.push("cell1".to_string());
    files.push("cell2".to_string());
    files.push("cell3".to_string());

    chrs.push("chr1".to_string());
    chrs.push("chr2".to_string());
    chrs.push("chr3".to_string());

	let mut tripletmatrix = TripletMatrix::new(
		chrs,
		&files,
	);
	
	println!("Files: {:?}", files);
	println!("Headers: {:?}", tripletmatrix.dimheaders);
	println!("Chrs: {:?}", tripletmatrix.dimchrs);
	
	tripletmatrix.add_triplet(1,1,100);
	tripletmatrix.add_triplet(2,2,100);
	tripletmatrix.add_triplet(3,3,100);
	
	println!("Tuple 0: {:?}",tripletmatrix.get_triplet(0));
	println!("Tuple 1: {:?}",tripletmatrix.get_triplet(1));
	println!("Tuple 2: {:?}",tripletmatrix.get_triplet(2));
	
	let hdf5file = hdf5::File::create("triplet_matrix.h5")?;
	let group = hdf5file.create_group("triplet_matrix")?;
	
	//let dimchrs = group.new_dataset::<String>().create("dimchrs", tripletmatrix.dimchrs.len())?;
	//dimchrs.write(tripletmatrix.dimchrs)?;
	//let i = group.new_dataset::<VarLenArray<u64>>().create("i", tripletmatrix.i.len())?;
	//i.write(&tripletmatrix.i)?;
	
	let mut jdata = VarLenArray::from_slice(&tripletmatrix.j[..]);
	println!("Jdata: {:?}",jdata);	
		
	let j = group.new_dataset::<VarLenArray<u8>>().create("j", (jdata.len()))?;
	j.write(jdata)?;
	
	//let x = group.new_dataset::<FixedAscii<u8>>().create("x", tripletmatrix.x.len())?;
	//x.write(&tripletmatrix.x)?;
	
	Ok(())
}

struct TripletMatrix {
	dimchrs: Vec<String>,         
	dimheaders: Vec<String>,   	//Start\tEnd\t[all cell filenames]
	i: Vec<u64>,				//Row (Chr)
	j: Vec<u8>,					//Col (Header)
	x: Vec<u8>,				    //Value
}

impl TripletMatrix {
	fn new(chrs: Vec<String>, headers: &Vec<String>) -> TripletMatrix {
	
	    let mut ranges: Vec<String> = Vec::new();
	    ranges.push("start".to_string());
	    ranges.push("end".to_string());
	    ranges.extend_from_slice(&headers);
	
		let dimchrs = chrs;
		let dimheaders = ranges;
		let i: Vec<u64> = Vec::new();
		let j: Vec<u8> = Vec::new();
		let x: Vec<u8> = Vec::new();
		TripletMatrix {dimchrs, dimheaders, i, j, x}
		
	}
		
	fn add_triplet(&mut self, row: u64, col: u8, val: u8) {
		self.i.push(row);
		self.j.push(col);
		self.x.push(val);	
	}
	
	fn get_triplet(&mut self, index: usize ) -> (u64,u8,u8) {
	   let tuple = (self.i[index], self.j[index], self.x[index]);  
	   return tuple
	}
}
