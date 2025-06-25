use std::{fmt::Debug, time::Instant};

use ark_bls12_381::Fr;
use ark_relations::r1cs::ConstraintSystemRef;
use ark_serialize::{CanonicalSerialize, Compress};
use ark_std::{rand::Rng, vec::Vec, UniformRand};
use arkworks_keccak::{
    constraints::{KeccakCircuit, KeccakMode},
    util::{bytes_to_bitvec, sha3_256, shake_128, shake_256},
};
use cyclotomic_rings::{
    challenge_set::LatticefoldChallengeSet,
    rings::{GoldilocksChallengeSet, GoldilocksRingNTT, SuitableRing},
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
use latticefold_keccak::circuit::{setup_environment, z_split};
include!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../target/debug/build/latticefold-735affedff5e3c44/out/examples_generated.rs"
));

type RqNTT = GoldilocksRingNTT;
type CS = GoldilocksChallengeSet;
type T = PoseidonTranscript<RqNTT, CS>;

fn main() {
    let mut rng = ark_std::rand::thread_rng();
    let preimage_length_bytes = rng.gen_range(1..=10);
    let preimage: Vec<u8> = (0..preimage_length_bytes).map(|_| rng.r#gen()).collect();
    let d: usize = rng.gen_range(50..=100);
    println!("input length: {} bytes", preimage.len());
    println!("d: {} bits", d);
    let expected = shake_128(&preimage, d / 8); // change
    let preimage = bytes_to_bitvec::<Fr>(&preimage);
    let circuit = KeccakCircuit::<Fr>::init_circuit(
        preimage.clone(),
        expected.to_vec(),
        KeccakMode::Shake128,
        d,
    );
    let cs: ConstraintSystemRef<Fr> = circuit.gen_cs();
    let is_satisfied = cs.is_satisfied().unwrap();
    assert!(is_satisfied);

    let (x_len, w_len, x_r1cs, w_r1cs) = z_split(cs.clone());
    assert_eq!(x_r1cs.len(), preimage_length_bytes * 8 + expected.len() * 8);

    println!("Setting up example environment...");

    println!("Decomposition parameters:");

    println!("\tB: {}", GoldilocksExampleDP::B);
    println!("\tL: {}", GoldilocksExampleDP::L);
    println!("\tB_SMALL: {}", GoldilocksExampleDP::B_SMALL);
    println!("\tK: {}", GoldilocksExampleDP::K);

    let (acc, wit_acc, cm_i, wit_i, ccs, scheme) =
        setup_environment::<C, RqNTT, GoldilocksExampleDP, CS, Fr>(cs, w_len*GoldilocksExampleDP::L);

    let mut prover_transcript = PoseidonTranscript::<RqNTT, CS>::default();
    let mut verifier_transcript = PoseidonTranscript::<RqNTT, CS>::default();
    println!("Generating proof...");
    let start = Instant::now();
    let (_, _, proof) = NIFSProver::<C, RqNTT, GoldilocksExampleDP, T>::prove(
        &acc,
        &wit_acc,
        &cm_i,
        &wit_i,
        &mut prover_transcript,
        &ccs,
        &scheme,
        w_len*GoldilocksExampleDP::L
    )
    .unwrap();
    let duration = start.elapsed();
    println!("Proof generated in {:?}", duration);

    let mut serialized_proof = Vec::new();

    println!("Serializing proof (with compression)...");
    proof
        .serialize_with_mode(&mut serialized_proof, Compress::Yes)
        .unwrap();
    let compressed_size = serialized_proof.len();
    println!(
        "Proof size (with compression) size: {}",
        humansize::format_size(compressed_size, humansize::BINARY)
    );

    println!("Serializing proof (without compression)...");
    proof
        .serialize_with_mode(&mut serialized_proof, Compress::No)
        .unwrap();
    let uncompressed_size = serialized_proof.len();
    println!(
        "Proof (without compression) size: {}\n",
        humansize::format_size(uncompressed_size, humansize::BINARY)
    );

    println!("Verifying proof");
    let start = Instant::now();
    NIFSVerifier::<C, RqNTT, GoldilocksExampleDP, T>::verify(
        &acc,
        &cm_i,
        &proof,
        &mut verifier_transcript,
        &ccs,
    )
    .unwrap();
    let duration = start.elapsed();
    println!("Proof verified in {:?}", duration);
}
