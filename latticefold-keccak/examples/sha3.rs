use ark_bls12_381::Fr;
use ark_relations::r1cs::ConstraintSystemRef;
use ark_std::{rand::Rng, vec::Vec};
use arkworks_keccak::{
    constraints::{KeccakCircuit, KeccakMode},
    util::{bytes_to_bitvec, shake_256},
};
use latticefold_keccak::circuit::z_split;

fn main() {
    let mut rng = ark_std::rand::thread_rng();
    let preimage_length_bytes = rng.gen_range(1..=256);
    let preimage: Vec<u8> = (0..preimage_length_bytes).map(|_| rng.r#gen()).collect();
    let d: usize = rng.gen_range(100..=4032);
    println!("input length: {} bits", preimage.len());
    println!("d: {}", d);
    let expected = shake_256(&preimage, d / 8); // change
    let preimage = bytes_to_bitvec::<Fr>(&preimage);
    let circuit = KeccakCircuit::<Fr>::init_circuit(
        preimage.clone(),
        expected.to_vec(),
        KeccakMode::Shake256,
        d,
    );
    let cs: ConstraintSystemRef<Fr> = circuit.gen_cs();
    let is_satisfied = cs.is_satisfied().unwrap();
    assert!(is_satisfied);

    let (x_len, w_len, x_r1cs, w_r1cs) = z_split(cs);
    assert_eq!(x_r1cs.len(), preimage_length_bytes * 8 + expected.len() * 8);
}
