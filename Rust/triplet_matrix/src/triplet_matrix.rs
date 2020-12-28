fn main() {

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
		let mut i: Vec<u64> = Vec::new();
		let mut j: Vec<u8> = Vec::new();
		let mut x: Vec<u8> = Vec::new();
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