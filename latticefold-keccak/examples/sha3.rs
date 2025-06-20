use std::{fmt::Debug, time::Instant};

use ark_bls12_381::Fr;
use ark_relations::r1cs::ConstraintSystemRef;
use ark_serialize::{CanonicalSerialize, Compress};
use ark_std::{rand::Rng, vec::Vec, UniformRand};
use arkworks_keccak::{
    constraints::{KeccakCircuit, KeccakMode},
    util::{bytes_to_bitvec, shake_256},
};
use cyclotomic_rings::{
    challenge_set::LatticefoldChallengeSet,
    rings::{BabyBearChallengeSet, BabyBearRingNTT, SuitableRing},
};
use latticefold::{
    arith::{
        ccs::get_test_dummy_degree_three_ccs_non_scalar, r1cs::get_test_dummy_z_split_ntt, Arith,
        Witness, CCCS, CCS, LCCCS,
    },
    commitment::AjtaiCommitmentScheme,
    nifs::{
        linearization::{LFLinearizationProver, LinearizationProver},
        NIFSProver, NIFSVerifier,
    },
    transcript::poseidon::PoseidonTranscript,
};
use latticefold_keccak::circuit::{z_split, setup_environment};
include!(concat!(env!("CARGO_MANIFEST_DIR"), "/../target/debug/build/latticefold-735affedff5e3c44/out/examples_generated.rs"));

type RqNTT = BabyBearRingNTT;
type CS = BabyBearChallengeSet;
type T = PoseidonTranscript<RqNTT, CS>;




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

    println!("Setting up example environment...");

    println!("Decomposition parameters:");
    
        println!("\tB: {}", BabyBearExampleDP::B);
        println!("\tL: {}", BabyBearExampleDP::L);
        println!("\tB_SMALL: {}", BabyBearExampleDP::B_SMALL);
        println!("\tK: {}", BabyBearExampleDP::K);

    let (acc, wit_acc, cm_i, wit_i, ccs, scheme) =
        setup_environment::<C, RqNTT, BabyBearExampleDP, W_BABYBEAR, CS, Fr>(cs);
}
