use std::{fmt::Debug, time::Instant};

use ark_bls12_381::Fr;
use ark_relations::r1cs::{ConstraintSystemRef, Field};
use ark_serialize::{CanonicalSerialize, Compress};
use ark_std::{rand::Rng, vec::Vec, UniformRand};
use arkworks_keccak::{
    constraints::{KeccakCircuit, KeccakMode},
    util::{bytes_to_bitvec, keccak256, sha3_256, shake_128, shake_256},
};
use cyclotomic_rings::{
    challenge_set::LatticefoldChallengeSet,
    rings::{BabyBearChallengeSet, BabyBearRingNTT, SuitableRing},
};
use latticefold::{
    arith::{Arith, Witness, CCCS, CCS, LCCCS},
    commitment::AjtaiCommitmentScheme,
    nifs::{
        linearization::{LFLinearizationProver, LinearizationProver},
        NIFSProver, NIFSVerifier,
    },
    transcript::poseidon::PoseidonTranscript,
};
use latticefold::decomposition_parameters::DecompositionParams;

pub fn z_split_lengths<F: Field>(cs: ConstraintSystemRef<F>) -> (usize, usize) {
    let x_len = cs.num_instance_variables() - 1;
    let w_len = cs.num_witness_variables();
    println!("Public inputs: {}", x_len);
    println!("Witnesses: {}", w_len);
    (x_len,w_len)
}

// fn r1cs_to_ccs<
//     const X_LEN: usize,
//     const C: usize, // rows
//     const WIT_LEN: usize,
//     const W: usize, // columns
//     P: DecompositionParams,
//     R: Clone + UniformRand + Debug + SuitableRing,
//     >(
//     r1cs_rows: usize,
//     ) -> (
//     CCCS<C, R>,
//     Witness<R>,
//     CCS<R>,
//     AjtaiCommitmentScheme<C, W, R>,
// ) {

// }

// fn setup_environment<
//     const X_LEN: usize,
//     const C: usize,
//     RqNTT: SuitableRing,
//     DP: DecompositionParams,
//     const W: usize,
//     const WIT_LEN: usize,
//     CS: LatticefoldChallengeSet<RqNTT>,
//     F: Field,
//     >(cs: ConstraintSystemRef<F>) -> (
//     LCCCS<C, RqNTT>,
//     Witness<RqNTT>,
//     CCCS<C, RqNTT>,
//     Witness<RqNTT>,
//     CCS<RqNTT>,
//     AjtaiCommitmentScheme<C, W, RqNTT>,
// ) {
//     // let (X_LEN, WIT_LEN) = r1cs_dimensions(cs);
//     let r1cs_rows = X_LEN + WIT_LEN + 1;
//     let (cm_i, wit, ccs, scheme) = z_split_lengths<X_LEN,C,WIT_LEN, W, DP, RqNTT>(r1cs_rows);

// }