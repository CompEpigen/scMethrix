mod collectargs;
mod joinbdgs;

//./target/debug/cptabix -f chr1_v2.bdg.gz -c hg38.chrsizes -b 10

fn main() {
    joinbdgs::inputs();
}
