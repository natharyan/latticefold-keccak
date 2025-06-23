#![allow(dead_code)]

use std::{fmt::Debug, time::Instant};

use ark_serialize::{CanonicalSerialize, Compress};
use ark_std::{vec::Vec, UniformRand};
use cyclotomic_rings::{
    challenge_set::LatticefoldChallengeSet,
    rings::{StarkChallengeSet, StarkRingNTT, SuitableRing},
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

include!(concat!(env!("OUT_DIR"), "/examples_generated.rs"));

#[allow(dead_code)]
pub fn wit_and_ccs_gen_degree_three_non_scalar<
    const X_LEN: usize,
    const C: usize, // rows
    const WIT_LEN: usize,
    const W: usize, // columns
    P: DecompositionParams,
    R: Clone + UniformRand + Debug + SuitableRing,
>(
    r1cs_rows: usize,
) -> (
    CCCS<C, R>,
    Witness<R>,
    CCS<R>,
    AjtaiCommitmentScheme<C, R>,
) {
    let mut rng = ark_std::test_rng();

    let new_r1cs_rows = if P::L == 1 && (WIT_LEN > 0 && (WIT_LEN & (WIT_LEN - 1)) == 0) {
        r1cs_rows - 2
    } else {
        r1cs_rows // This makes a square matrix but is too much memory
    };
    let (one, x_ccs, w_ccs) = get_test_dummy_z_split_ntt::<R, X_LEN, WIT_LEN>();

    let mut z = vec![one];
    z.extend(&x_ccs);
    z.extend(&w_ccs);
    let ccs: CCS<R> =
        get_test_dummy_degree_three_ccs_non_scalar::<R, X_LEN, WIT_LEN, W>(&z, P::L, new_r1cs_rows);
    ccs.check_relation(&z).expect("R1CS invalid!");

    let scheme: AjtaiCommitmentScheme<C, R> = AjtaiCommitmentScheme::rand(&mut rng, W);
    let wit: Witness<R> = Witness::from_w_ccs::<P>(w_ccs);

    let cm_i: CCCS<C, R> = CCCS {
        cm: wit.commit::<C, P>(&scheme, W).unwrap(),
        x_ccs,
    };

    (cm_i, wit, ccs, scheme)
}

#[allow(clippy::type_complexity)]
fn setup_example_environment<
    const X_LEN: usize,
    const C: usize,
    RqNTT: SuitableRing,
    DP: DecompositionParams,
    const W: usize,
    const WIT_LEN: usize,
    CS: LatticefoldChallengeSet<RqNTT>,
>() -> (
    LCCCS<C, RqNTT>,
    Witness<RqNTT>,
    CCCS<C, RqNTT>,
    Witness<RqNTT>,
    CCS<RqNTT>,
    AjtaiCommitmentScheme<C, RqNTT>,
) {
    let r1cs_rows = X_LEN + WIT_LEN + 1;

    let (cm_i, wit, ccs, scheme) =
        wit_and_ccs_gen_degree_three_non_scalar::<X_LEN, C, WIT_LEN, W, DP, RqNTT>(r1cs_rows);

    let rand_w_ccs: Vec<RqNTT> = (0..WIT_LEN).map(|i| RqNTT::from(i as u64)).collect();
    let wit_acc = Witness::from_w_ccs::<DP>(rand_w_ccs);

    let mut transcript = PoseidonTranscript::<RqNTT, CS>::default();

    let (acc, _) = LFLinearizationProver::<_, PoseidonTranscript<RqNTT, CS>>::prove(
        &cm_i,
        &wit_acc,
        &mut transcript,
        &ccs,
    )
    .expect("Failed to generate linearization proof");

    (acc, wit_acc, cm_i, wit, ccs, scheme)
}

type RqNTT = StarkRingNTT;
type CS = StarkChallengeSet;
type T = PoseidonTranscript<RqNTT, CS>;

fn main() {
    println!("Setting up example environment...");

    println!("Decomposition parameters:");
    println!("\tB: {}", StarkPrimeExampleDP::B);
    println!("\tL: {}", StarkPrimeExampleDP::L);
    println!("\tB_SMALL: {}", StarkPrimeExampleDP::B_SMALL);
    println!("\tK: {}", StarkPrimeExampleDP::K);

    let (acc, wit_acc, cm_i, wit_i, ccs, scheme) = setup_example_environment::<
        X_LEN_STARK,
        C_STARK,
        RqNTT,
        StarkPrimeExampleDP,
        W_STARK,
        WIT_LEN_STARK,
        CS,
    >();

    let mut prover_transcript = PoseidonTranscript::<RqNTT, CS>::default();
    let mut verifier_transcript = PoseidonTranscript::<RqNTT, CS>::default();
    println!("Generating proof...");
    let start = Instant::now();

    let (_, _, proof) = NIFSProver::<C_STARK, RqNTT, StarkPrimeExampleDP, T>::prove(
        &acc,
        &wit_acc,
        &cm_i,
        &wit_i,
        &mut prover_transcript,
        &ccs,
        &scheme,
        W_STARK
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
        "Proof (without compression) size: {}",
        humansize::format_size(uncompressed_size, humansize::BINARY)
    );

    println!("Verifying proof");
    let start = Instant::now();
    NIFSVerifier::<C_STARK, RqNTT, StarkPrimeExampleDP, T>::verify(
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
